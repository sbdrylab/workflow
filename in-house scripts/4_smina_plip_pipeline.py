#!/usr/bin/env python3
"""
Refine each candidate pose and measure its contacts with the target pocket.

WHAT THIS PROGRAM DOES
    Runs after 3_filter_candidates.py. For every candidate drug it:

      1. Finds that drug's DiffDock rank-1 pose.
      2. Relaxes the pose with smina (--local_only --minimize). This is a
         small local relaxation, not a re-dock: it eases minor steric
         strain and returns a more reliable binding-affinity estimate
         along with the adjusted pose.
      3. Merges the relaxed ligand and the receptor into one complex PDB.
      4. Runs PLIP on that complex and counts the interactions that
         involve the target residue group: hydrogen bonds, hydrophobic
         contacts, salt bridges and pi-stacking.
      5. Writes one row per drug.

    This is the step that tests whether a pose actually interacts with
    the pocket, as opposed to merely sitting close to it.

INPUT
    --receptor        Protonated receptor PDB, used to build the complex
                      for PLIP. Must match whatever --receptor-pdbqt was
                      generated from.
    --receptor-pdbqt  Receptor in PDBQT form, used by smina.
    --results-root    DiffDock results tree, e.g. results/trem2_wt
    --candidates      CSV with a "drug" column, from 3_filter_candidates.py

OUTPUT
    A CSV with one row per drug and these nine columns:

        drug                               folder name of the drug
        status                             ok, no_pose_found, or
                                           smina_minimize_failed
        pose_in_box                        True if the pose stayed inside
                                           the fixed search box
        smina_minimized_affinity_kcal_mol  refined affinity, more negative
                                           means stronger predicted binding
        plip_hbonds_target                 hydrogen bonds to the group
        plip_hydrophobic_target            hydrophobic contacts to the group
        plip_saltbridge_target             salt bridges to the group
        plip_pistack_target                pi-stacking contacts to the group
        plip_target_residues_contacted     which target residues were hit,
                                           semicolon separated

    Intermediate files for each drug (smina log, relaxed pose, complex
    PDB, PLIP output) are kept under --workdir for inspection.

HOW TO RUN
    Use a separate conda environment from the DiffDock one, since smina,
    openbabel and PLIP bring their own dependency stack:

        conda create -n smina_plip -c conda-forge python=3.10 smina openbabel -y
        conda activate smina_plip
        pip install plip pdb2pqr lxml pandas

    Prepare the receptor once, before running this script. PDB2PQR adds
    hydrogens and PROPKA assigns titration states at pH 7.4:

        pdb2pqr30 --ff=AMBER --with-ph=7.4 --titration-state-method=propka \
            --pdb-output=receptor_protonated.pdb \
            receptor.pdb receptor.pqr
        obabel receptor_protonated.pdb -O receptor.pdbqt -xr

    Then run:

        python3 4_smina_plip_pipeline.py \
            --receptor receptor_protonated.pdb \
            --receptor-pdbqt receptor.pdbqt \
            --results-root result_directory \
            --candidates trem2_surface1_candidates.csv \
            --top-n 400 \
            --residues 47 65 67 68 77 \
            --center 49.3 32.1 -12.1 \
            --box-size 35 35 35 \
            --out smina_plip_report.csv

THREE SETTINGS THAT MUST BE CHECKED FOR EVERY RECEPTOR
    1. --top-n defaults to 20, which silently processes only the first 20
       candidates. Set it higher than the number of rows in --candidates
       so the whole pool is scored. When comparing two receptors, both
       runs must process all of their candidates, or the two funnels
       cannot be compared.

    2. --center must be recomputed for each receptor. The default
       (49.3, 32.1, -12.1) is the midpoint of residues 47/65/67/68/77 in
       one particular structure. Predicted models are each placed in
       their own coordinate frame, so a centre taken from a different
       structure will point at empty space.

    3. --residues must match the numbering used inside --receptor.
       Predicted models are often numbered from 1 rather than in UniProt
       numbering, in which case the literature residue numbers are wrong
       for that file and must be re-mapped by sequence alignment first.

    Nothing in this script can detect any of these three mistakes. Each
    one produces plausible-looking numbers that are measured against the
    wrong place.

"""

import argparse
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

AFFINITY_RE = re.compile(r"Affinity:\s*([+-]?\d+\.?\d*)")
# smina prints the affinity in two different shapes:
#   --score_only:              "Affinity: -5.73 (kcal/mol)"
#   --local_only --minimize:   "Affinity: -3.71561  0.00000 (kcal/mol)"
# The regex above takes the first number after "Affinity:" in either case,
# rather than requiring it to sit directly before "(kcal/mol)".
#
# The pattern below is a fallback for smina builds that print the
# "mode | affinity | dist from best mode" table instead.
TABLE_ROW_RE = re.compile(r"^\s*1\s+(-?\d+\.\d+)\s+\d+\.\d+\s+\d+\.\d+",
                           re.MULTILINE)


def parse_smina_affinity(text):
    """Try both smina output formats and return the first match."""
    m = AFFINITY_RE.search(text)
    if m:
        return float(m.group(1))
    m = TABLE_ROW_RE.search(text)
    if m:
        return float(m.group(1))
    return None


# Default search box, centred on residues 47/65/67/68/77 in one specific
# receptor. Recompute --center for any other structure.
#
# The box is deliberately generous rather than tight around the residue
# group. Candidates are selected on their closest atom reaching the
# pocket, but the rest of a drug-sized ligand can extend well beyond that
# atom in any direction. Because --local_only refines a pose where it
# already sits and cannot move it to another site, a wide box costs
# nothing: it only needs to avoid clipping part of the ligand, not to
# constrain the search.
DEFAULT_CENTER = (49.3, 32.1, -12.1)
DEFAULT_BOX_SIZE = (35.0, 35.0, 35.0)


def read_sdf_coords(sdf_path):
    """Read atom coordinates from the first molecule of a V2000 SDF.

    Deliberately minimal: just enough to check a pose against the search
    box, so this script does not need RDKit as a dependency.
    """
    coords = []
    with open(sdf_path) as f:
        lines = f.readlines()
    if len(lines) < 4:
        return coords
    try:
        n_atoms = int(lines[3][0:3])
    except ValueError:
        return coords
    for line in lines[4:4 + n_atoms]:
        try:
            x, y, z = float(line[0:10]), float(line[10:20]), float(line[20:30])
            coords.append((x, y, z))
        except ValueError:
            continue
    return coords


def check_pose_in_box(sdf_path, center, box_size):
    """Check whether every atom of a pose lies inside the search box.

    With a fixed box rather than --autobox_ligand, a pose that has drifted
    away from the target region will not be minimized or scored correctly,
    so the caller flags it instead of trusting the numbers.
    """
    coords = read_sdf_coords(sdf_path)
    if not coords:
        return True  # couldn't check; don't block on it
    lo = [c - s / 2 for c, s in zip(center, box_size)]
    hi = [c + s / 2 for c, s in zip(center, box_size)]
    for (x, y, z) in coords:
        if not (lo[0] <= x <= hi[0] and lo[1] <= y <= hi[1] and lo[2] <= z <= hi[2]):
            return False
    return True


RANK_NUM_RE = re.compile(r"^rank(\d+)", re.IGNORECASE)


def find_top_pose(results_root, drug):
    """Find this drug's DiffDock rank-1 pose SDF.

    Matches on the parsed rank number, not on a filename prefix. A glob
    such as 'rank1*.sdf' also matches 'rank10_confidence-X.XX.sdf', and
    that collision once made this pipeline use every drug's worst pose
    instead of its best.

    Prefers the confidence-named file over a bare rank1.sdf, matching the
    behaviour of 2_scan_diffdock_results.py.
    """
    d = Path(results_root) / drug
    cands = []
    for p in d.glob("*.sdf"):
        m = RANK_NUM_RE.match(p.name)
        if m and int(m.group(1)) == 1:
            cands.append(p)
    cands = sorted(cands)
    if not cands:
        return None
    with_conf = [c for c in cands if "confidence" in c.name.lower()]
    return with_conf[0] if with_conf else cands[0]


def run_smina(receptor_pdbqt, ligand_sdf, out_dir, tag, minimize=False,
              center=DEFAULT_CENTER, box_size=DEFAULT_BOX_SIZE):
    """Run smina in score_only or local_only+minimize mode.

    Uses a fixed search box rather than --autobox_ligand, so every
    candidate is scored against the same region instead of a box that
    re-centres on wherever each pose happens to sit.

    Returns (affinity in kcal/mol or None, path to the output SDF or None).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / f"{tag}.log"
    cmd = ["smina", "-r", str(receptor_pdbqt), "-l", str(ligand_sdf),
           "--log", str(log_path)]
    out_sdf = None
    if minimize:
        out_sdf = out_dir / f"{tag}_minimized.sdf"
        cmd += ["--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(box_size[0]),
                "--size_y", str(box_size[1]),
                "--size_z", str(box_size[2]),
                "--local_only", "--minimize", "--minimize_iters", "2000",
                "-o", str(out_sdf)]
    else:
        cmd += ["--score_only"]

    result = subprocess.run(cmd, capture_output=True, text=True)
    text = result.stdout + "\n" + result.stderr
    if log_path.exists():
        text += "\n" + log_path.read_text()

    affinity = parse_smina_affinity(text)
    if affinity is None:
        print(f"    WARNING: could not parse smina affinity for {tag} "
              f"(stderr tail: {result.stderr[-200:]!r})", file=sys.stderr)
    return affinity, (out_sdf if (out_sdf and out_sdf.exists()) else None)


def merge_complex(receptor_pdb, ligand_pdb, out_pdb, lig_chain="X",
                   lig_resname="LIG"):
    """Join receptor and ligand into one complex PDB for PLIP.

    Ligand atoms are written as HETATM records relabelled onto
    lig_chain and lig_resname.
    """
    rec_lines = []
    with open(receptor_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM", "TER")):
                rec_lines.append(line.rstrip("\n"))

    lig_lines = []
    with open(ligand_pdb) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                line = line.rstrip("\n").ljust(80)
                line = ("HETATM" + line[6:17] +
                        f"{lig_resname:>3}" + line[20:21] +
                        f"{lig_chain:>1}" + f"{1:>4}" + line[26:])
                lig_lines.append(line)

    with open(out_pdb, "w") as f:
        f.write("\n".join(rec_lines) + "\n")
        f.write("TER\n")
        f.write("\n".join(lig_lines) + "\n")
        f.write("TER\nEND\n")


def run_plip(complex_pdb, out_dir):
    """Run PLIP on one complex PDB.

    Returns the path to its XML report, or None when PLIP produced none,
    which happens when it detects no binding site at all, for example if
    the ligand ended up far from the protein.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["plip", "-f", str(complex_pdb), "-o", str(out_dir), "-x", "-q"]
    subprocess.run(cmd, capture_output=True, text=True)
    xml_candidates = list(out_dir.glob("*.xml"))
    return xml_candidates[0] if xml_candidates else None


OUTPUT_COLUMNS = [
    "drug", "status", "pose_in_box", "smina_minimized_affinity_kcal_mol",
    "plip_hbonds_target", "plip_hydrophobic_target",
    "plip_saltbridge_target", "plip_pistack_target",
    "plip_target_residues_contacted",
]

EMPTY_PLIP_FIELDS = {
    "plip_hbonds_target": None,
    "plip_hydrophobic_target": None,
    "plip_saltbridge_target": None,
    "plip_pistack_target": None,
    "plip_target_residues_contacted": "",
}


def parse_plip_xml(xml_path, target_residues):
    """Count the interactions involving target_residues in a PLIP report.

    PLIP reports interactions across the whole complex. Those are read
    here because they are needed to work out which target residues were
    contacted, but only the counts at the target residues are written to
    the output row; the later filtering steps use nothing else.
    """
    if xml_path is None or not Path(xml_path).exists():
        return dict(EMPTY_PLIP_FIELDS)

    tree = ET.parse(xml_path)
    root = tree.getroot()

    def collect(tag):
        rows = []
        for site in root.iter("bindingsite"):
            for group in site.iter(tag):
                for item in group:
                    resnr_el = item.find("resnr")
                    if resnr_el is not None and resnr_el.text:
                        rows.append(int(resnr_el.text))
        return rows

    hbonds = collect("hydrogen_bonds")
    hydrophobic = collect("hydrophobic_interactions")
    saltbridges = collect("salt_bridges")
    pistack = collect("pi_stacks")

    target_residues = set(target_residues)
    all_hits = hbonds + hydrophobic + saltbridges + pistack
    contacted = sorted(set(r for r in all_hits if r in target_residues))

    return {
        "plip_hbonds_target": sum(1 for r in hbonds if r in target_residues),
        "plip_hydrophobic_target": sum(1 for r in hydrophobic if r in target_residues),
        "plip_saltbridge_target": sum(1 for r in saltbridges if r in target_residues),
        "plip_pistack_target": sum(1 for r in pistack if r in target_residues),
        "plip_target_residues_contacted": ";".join(str(r) for r in contacted),
    }


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--receptor", required=True,
                   help="Protonated receptor PDB (e.g. "
                        "trem2_wt_protonated.pdb from the pdb2pqr step) -- "
                        "used to build the PLIP complex file, so it should "
                        "match what --receptor-pdbqt was generated from")
    p.add_argument("--receptor-pdbqt", required=True,
                   help="Pre-generated receptor PDBQT -- see prep step in "
                        "the script docstring")
    p.add_argument("--results-root", required=True,
                   help="DiffDock results root, e.g. results/trem2_wt")
    p.add_argument("--candidates", required=True,
                   help="CSV with a 'drug' column, e.g. "
                        "trem2_surface1_candidates.csv")
    p.add_argument("--top-n", type=int, default=20,
                   help="Only process the first N rows of --candidates "
                        "(it should already be sorted by confidence)")
    p.add_argument("--residues", nargs="+", type=int,
                   default=[47, 65, 67, 68, 77],
                   help="Target residue numbers for the PLIP contact "
                        "summary (default: the Surface1/CDR2 group)")
    p.add_argument("--center", nargs=3, type=float, default=list(DEFAULT_CENTER),
                   metavar=("X", "Y", "Z"),
                   help="Fixed smina search-box center. Default is the "
                        "bounding-box midpoint of residues 47/65/67/68/77 "
                        "in trem2_wt.pdb (49.3 32.1 -12.1). Recompute if "
                        "you change the residue group or the receptor.")
    p.add_argument("--box-size", nargs=3, type=float, default=list(DEFAULT_BOX_SIZE),
                   metavar=("X", "Y", "Z"),
                   help="Fixed smina search-box size, in Angstroms "
                        "(default: 35 35 35 -- generous enough that a "
                        "drug-sized ligand near the pocket isn't clipped; "
                        "safe with --local_only since the pose isn't "
                        "searched globally within this box).")
    p.add_argument("--workdir", default="smina_plip_work",
                   help="Scratch directory for intermediate files")
    p.add_argument("--out", default="smina_plip_report.csv")
    args = p.parse_args()

    for exe in ("smina", "plip", "obabel"):
        if shutil.which(exe) is None:
            sys.exit(f"ERROR: '{exe}' not found on PATH. Activate the "
                      f"smina_plip conda env first (see script docstring).")

    df = pd.read_csv(args.candidates)
    drugs = df["drug"].head(args.top_n).tolist()
    print(f"Running smina + PLIP on {len(drugs)} candidates...\n")

    workdir = Path(args.workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for n, drug in enumerate(drugs, 1):
        print(f"[{n}/{len(drugs)}] {drug}")
        pose = find_top_pose(args.results_root, drug)
        if pose is None:
            print("    no rank1 pose found, skipping")
            rows.append({"drug": drug, "status": "no_pose_found"})
            continue

        drug_dir = workdir / drug.replace(" ", "_").replace("/", "_")
        drug_dir.mkdir(exist_ok=True)

        in_box = check_pose_in_box(pose, args.center, args.box_size)
        if not in_box:
            print("    WARNING: pose extends outside the fixed docking "
                  "box -- smina minimize/score below may be unreliable "
                  "for this drug (consider --box-size if this happens a lot)")

        min_affinity, min_sdf = run_smina(args.receptor_pdbqt, pose, drug_dir,
                                           "minimized", minimize=True,
                                           center=args.center,
                                           box_size=args.box_size)

        plip_fields = dict(EMPTY_PLIP_FIELDS)
        if min_sdf is not None:
            lig_pdb = drug_dir / "minimized_ligand.pdb"
            subprocess.run(["obabel", str(min_sdf), "-O", str(lig_pdb)],
                            capture_output=True, text=True)
            if lig_pdb.exists():
                complex_pdb = drug_dir / "complex.pdb"
                merge_complex(args.receptor, lig_pdb, complex_pdb)
                xml_path = run_plip(complex_pdb, drug_dir / "plip_out")
                plip_fields = parse_plip_xml(xml_path, args.residues)

        rows.append({
            "drug": drug,
            "status": "ok" if min_sdf is not None else "smina_minimize_failed",
            "pose_in_box": in_box,
            "smina_minimized_affinity_kcal_mol": min_affinity,
            **plip_fields,
        })

    out_df = pd.DataFrame(rows).reindex(columns=OUTPUT_COLUMNS)
    out_df.to_csv(args.out, index=False)
    print(f"\nWrote {len(out_df)} rows to {args.out}")
    if "smina_minimized_affinity_kcal_mol" in out_df.columns:
        ok = out_df.dropna(subset=["smina_minimized_affinity_kcal_mol"])
        ok = ok.sort_values("smina_minimized_affinity_kcal_mol")
        print("\nTop by smina minimized affinity (more negative = stronger):")
        cols = ["drug", "smina_minimized_affinity_kcal_mol",
                "plip_hbonds_target", "plip_target_residues_contacted"]
        print(ok[cols].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
