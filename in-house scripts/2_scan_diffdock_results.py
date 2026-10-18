"""
Score DiffDock's top pose for every docked drug.

WHAT THIS PROGRAM DOES
    Walks a DiffDock output tree and, for each ligand folder, measures that
    drug's rank-1 (highest-confidence) pose against the receptor:
      - how close it gets to a chosen target residue
      - how close it gets to an optional group of residues
      - whether it physically clashes with the receptor
      - how many distinct receptor residues it touches

    Only the rank-1 pose is loaded and scored. Lower-ranked poses are
    ignored because no downstream step uses them, which also keeps the
    run fast.

INPUT
    1. receptor_pdb   The same receptor PDB that was given to DiffDock.
    2. results_root   A DiffDock output tree laid out like this:

           <root>/
             drug_name/
               rank1_confidence-0.39.sdf
               rank1_reverseprocess.pdb      <- ignored
               rank2_confidence-0.77.sdf
               ...

    The ligand folder name becomes the drug name. The confidence score is
    read from the filename, so "confidence-0.39" becomes -0.39. DiffDock
    sometimes writes rank 1 twice, as a bare "rank1.sdf" and as
    "rank1_confidence-X.XX.sdf" holding the identical pose; the
    confidence-named file is preferred so the score is available.

OUTPUT
    One CSV row per drug:

        drug                 ligand folder name
        top_confidence       DiffDock confidence for the rank-1 pose
        top_clash            True if any ligand atom lies within 1.5 A of a
                             receptor atom, which is not a physically
                             plausible pose
        top_dist_target      distance in A from the pose to --residue
        top_n_contacts       distinct receptor residues within --cutoff
        top_min_dist_any     closest approach in A to any receptor atom
        n_heavy_atoms        heavy-atom count of the pose
        top_dist_<label>     distance in A to the --extra-residues group,
                             present only when that option is used

HOW TO RUN
    Requires:
        pip install rdkit numpy scipy pandas

    Basic run:
        python3 2_scan_diffdock_results.py receptor.pdb /path/to/results \
            --out report.csv

    Also measure distance to a group of residues, such as CDR2/Surface1:
        python3 2_scan_diffdock_results.py receptor.pdb ./data \
            --residue 47 --extra-residues 65 67 68 77 \
            --extra-label surface1_cdr2 --out report.csv

    Options:
        --residue         target residue number (default 47)
        --chain           restrict the target residue to a single chain
        --extra-residues  additional residues to measure as one group
        --extra-label     column label for that group (default ps_site)
        --cutoff          contact distance in A (default 4.5)
        --sort-by         top_confidence or top_dist_target
        --out             output CSV path (default diffdock_report.csv)

IMPORTANT
    The residue numbers given here must match the numbering inside
    receptor_pdb. Predicted structures are often numbered from 1 rather
    than in UniProt numbering, so check the receptor before relying on a
    residue number taken from the literature.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")   # silence RDKit sanitization warnings

# Matches "...confidence-0.39.sdf" and "...confidence1.42.sdf".
CONF_RE = re.compile(r"confidence(-?\d+(?:\.\d+)?)", re.IGNORECASE)
RANK_RE = re.compile(r"rank(\d+)", re.IGNORECASE)

# A ligand atom closer than this to any receptor atom means the pose is
# not physically plausible.
CLASH_DIST = 1.5  # angstroms

OUTPUT_COLUMNS = [
    "drug", "top_confidence", "top_clash", "top_dist_target",
    "top_n_contacts", "top_min_dist_any", "n_heavy_atoms",
]  # The extra-residue-group column (e.g. top_dist_surface1_cdr2) is added
   # after this list, once its label is known at run time.


# ----------------------------------------------------------------------
# Receptor
# ----------------------------------------------------------------------

def load_receptor(pdb_path):
    """Read receptor atoms from a PDB file.

    Returns (coords, info), where info holds one
    (chain, resnum, resname, atomname) tuple per atom.
    """
    coords, info = [], []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM"):
                coords.append([float(line[30:38]),
                               float(line[38:46]),
                               float(line[46:54])])
                info.append((line[21],
                             int(line[22:26]),
                             line[17:20].strip(),
                             line[12:16].strip()))
    if not coords:
        raise ValueError(f"No ATOM records read from {pdb_path}")
    return np.array(coords), info


def residue_atom_indices(info, resnum, chain=None):
    return [i for i, (c, r, _, _) in enumerate(info)
            if r == resnum and (chain is None or c == chain)]


# ----------------------------------------------------------------------
# Ligand poses
# ----------------------------------------------------------------------

def parse_pose_filename(path):
    """Extract (rank, confidence) from a DiffDock SDF filename."""
    name = Path(path).name
    rank_m = RANK_RE.search(name)
    conf_m = CONF_RE.search(name)
    rank = int(rank_m.group(1)) if rank_m else None
    conf = float(conf_m.group(1)) if conf_m else None
    return rank, conf


def load_ligand_coords(sdf_path):
    """Return heavy-atom coordinates of the first molecule in an SDF."""
    supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=True, sanitize=False)
    mol = next((m for m in supplier if m is not None), None)
    if mol is None:
        return None, None
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i))
                       for i in range(mol.GetNumAtoms())])
    return coords, mol


def find_ligand_dirs(root):
    """Find every directory that directly contains a rank*.sdf file."""
    root = Path(root)
    dirs = set()
    for sdf in root.rglob("*.sdf"):
        if RANK_RE.search(sdf.name):
            dirs.add(sdf.parent)
    return sorted(dirs)


def find_top_pose_path(lig_dir):
    """Pick the file holding DiffDock's rank-1 pose.

    Matches on the parsed rank number, not on a filename prefix. A glob
    such as 'rank1*.sdf' would also match 'rank10_confidence-X.XX.sdf',
    silently selecting one of the worst poses instead of the best.

    When DiffDock has written rank 1 twice (a bare 'rank1.sdf' plus a
    'rank1_confidence-X.XX.sdf' holding the same pose), the
    confidence-named file is preferred so the score is available.

    Only filenames are read here, so coordinates are loaded for one file
    per drug rather than for every sampled pose.

    Falls back to the highest-confidence pose available when no rank-1
    file exists, which happens for a small number of drugs.
    """
    sdf_paths = [p for p in lig_dir.glob("*.sdf") if RANK_RE.search(p.name)]
    if not sdf_paths:
        return None, None

    parsed = [(p,) + parse_pose_filename(p) for p in sdf_paths]  # (path, rank, conf)

    rank1 = [t for t in parsed if t[1] == 1]
    if rank1:
        with_conf = [t for t in rank1 if t[2] is not None]
        chosen = with_conf[0] if with_conf else rank1[0]
        return chosen[0], chosen[2]

    with_conf = [t for t in parsed if t[2] is not None]
    if not with_conf:
        return None, None
    chosen = max(with_conf, key=lambda t: t[2])
    return chosen[0], chosen[2]


# ----------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------

def analyse_ligand_dir(lig_dir, rec_coords, info, target_idx,
                       extra_idx_map, contact_cutoff):
    """Score the rank-1 pose for one drug.

    Returns a dict of results, or None when no usable pose was found.
    """
    top_path, top_confidence = find_top_pose_path(lig_dir)
    if top_path is None:
        return None

    coords, mol = load_ligand_coords(top_path)
    if coords is None or len(coords) == 0:
        return None

    d_all = cdist(coords, rec_coords)
    dist_target = float(d_all[:, target_idx].min()) if target_idx else float("nan")
    min_dist_any = float(d_all.min())

    extras = {}
    for label, idx in extra_idx_map.items():
        extras[label] = float(d_all[:, idx].min()) if idx else float("nan")

    close = np.where(d_all.min(axis=0) < contact_cutoff)[0]
    n_contacts = len({(info[i][0], info[i][1]) for i in close})  # distinct residues

    row = {
        "drug": lig_dir.name,
        "top_confidence": top_confidence,
        "top_clash": bool(min_dist_any < CLASH_DIST),
        "top_dist_target": round(dist_target, 2),
        "top_n_contacts": n_contacts,
        "top_min_dist_any": round(min_dist_any, 2),
        "n_heavy_atoms": mol.GetNumAtoms(),
    }
    for label in extra_idx_map:
        row[f"top_dist_{label}"] = (
            round(extras[label], 2) if not np.isnan(extras[label]) else np.nan)

    return row


# ----------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("receptor_pdb",
                   help="The SAME receptor PDB that was given to DiffDock")
    p.add_argument("results_root",
                   help="Root directory containing drug_name/ folders")
    p.add_argument("--residue", type=int, default=47,
                   help="Target residue number (default: 47)")
    p.add_argument("--chain", default=None,
                   help="Restrict the target residue to one chain "
                        "(5UD7 has six copies of residue 47)")
    p.add_argument("--extra-residues", nargs="*", type=int, default=[],
                   help="Additional residues to report distance to, as a group "
                        "(e.g. the CDR2/Surface1 group)")
    p.add_argument("--extra-label", default="ps_site",
                   help="Column label for the extra-residue group")
    p.add_argument("--cutoff", type=float, default=4.5,
                   help="Contact distance cutoff, angstroms (default: 4.5)")
    p.add_argument("--sort-by", default="top_confidence",
                   choices=["top_confidence", "top_dist_target"],
                   help="Column to sort the report by (default: "
                        "top_confidence, descending)")
    p.add_argument("--out", default="diffdock_report.csv")
    args = p.parse_args()

    rec_coords, info = load_receptor(args.receptor_pdb)
    chains = sorted({c for c, _, _, _ in info})
    print(f"Receptor: {args.receptor_pdb}")
    print(f"  {len(rec_coords)} atoms, chains: {chains}")

    target_idx = residue_atom_indices(info, args.residue, args.chain)
    if not target_idx:
        sys.exit(f"ERROR: residue {args.residue}"
                 f"{' chain ' + args.chain if args.chain else ''} "
                 f"not found in receptor.")

    copies = sorted({info[i][0] for i in target_idx})
    resname = info[target_idx[0]][2]
    print(f"  target: {resname}{args.residue}, "
          f"{len(target_idx)} atoms across chain(s) {copies}")
    if args.chain is None and len(copies) > 1:
        print("  NOTE: multiple copies present; distance is to the nearest copy. "
              "Use --chain to restrict to one.")

    extra_idx_map = {}
    if args.extra_residues:
        idx = [i for i, (c, r, _, _) in enumerate(info)
               if r in args.extra_residues]
        extra_idx_map[args.extra_label] = idx
        print(f"  extra group '{args.extra_label}': "
              f"{len(args.extra_residues)} residues, {len(idx)} atoms")

    lig_dirs = find_ligand_dirs(args.results_root)
    print(f"\nFound {len(lig_dirs)} ligand folders under {args.results_root}\n")
    if not lig_dirs:
        sys.exit("No rank*.sdf files found. Check the path.")

    rows, failed = [], []
    for n, d in enumerate(lig_dirs, 1):
        try:
            row = analyse_ligand_dir(d, rec_coords, info, target_idx,
                                     extra_idx_map, args.cutoff)
            if row:
                rows.append(row)
            else:
                failed.append(str(d))
        except Exception as e:          # log and continue; one bad folder
                                        # must not stop the whole batch
            failed.append(f"{d}  ({e})")
        if n % 100 == 0:
            print(f"  processed {n}/{len(lig_dirs)}")

    if not rows:
        sys.exit("No poses could be analysed.")

    df = pd.DataFrame(rows)

    # Fixed column order: drug, then confidence/clash/distance statistics,
    # with the extra-residue-group distance last when present.
    out_cols = list(OUTPUT_COLUMNS)
    for label in extra_idx_map:
        out_cols.append(f"top_dist_{label}")
    df = df[out_cols]

    # Sorting direction depends on the column: shorter distance is better,
    # higher confidence is better.
    ascending = "dist" in args.sort_by
    df = df.sort_values(args.sort_by, ascending=ascending).reset_index(drop=True)
    df.to_csv(args.out, index=False)

    print(f"\nWrote {len(df)} rows to {args.out}")
    if failed:
        print(f"{len(failed)} folder(s) could not be read:")
        for f in failed[:10]:
            print(f"   {f}")
        if len(failed) > 10:
            print(f"   ... and {len(failed) - 10} more")

    print(f"\nTop 15 by {args.sort_by}:")
    cols = ["drug", "top_confidence", "top_dist_target", "top_clash", "top_n_contacts"]
    for label in extra_idx_map:
        cols.append(f"top_dist_{label}")
    print(df[cols].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
