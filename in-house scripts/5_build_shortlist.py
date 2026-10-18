#!/usr/bin/env python3
"""
Build the final ranked shortlist of candidate drugs.

WHAT THIS PROGRAM DOES
    Turns the smina/PLIP report into a shortlist, rather than simply
    sorting every row by affinity.

    The reason this step exists: sitting near the target pocket is not the
    same as interacting with it. Most candidates that pass the earlier
    proximity screen turn out to form no hydrogen bond or salt bridge to
    the target residues at all; they are merely nearby. Ranking by
    affinity alone, or by DiffDock confidence alone, puts those
    "near but not touching" poses straight at the top.

    So the filters run in this order:

      1. Drop unusable rows: status is not "ok", or the pose fell outside
         the search box, or the affinity is missing.
      2. Drop rows with a non-negative affinity. After local relaxation a
         value at or above zero means the pose never resolved into
         anything favourable, which is a failed pose rather than a weak
         binder. PLIP cannot tell the difference, because it only checks
         geometry.
      3. Apply the hard engagement filter: keep only rows with at least
         one hydrogen bond or salt bridge to the target residue group.
         This is the step that enforces engagement over proximity.
      4. Optionally apply the extra filters described below.
      5. Rank the survivors by affinity, or by a consensus of affinity and
         DiffDock confidence when --diffdock-candidates is supplied.

INPUT
    --smina-plip           Report CSV from 4_smina_plip_pipeline.py.
    --diffdock-candidates  Optional. Candidates CSV from
                           3_filter_candidates.py, used to add DiffDock
                           confidence and switch to consensus ranking.
    --smiles-csv           Required only with --require-cns-druglike. Maps
                           drug name to SMILES, so molecular properties
                           can be computed; the smina/PLIP report carries
                           no structures.

OUTPUT
    A CSV holding one row per surviving drug:

        shortlist_position, drug, smina_minimized_affinity_kcal_mol,
        plip_hbonds_target, plip_hydrophobic_target,
        plip_saltbridge_target, plip_pistack_target,
        plip_target_residues_contacted,
        n_interaction_types, n_target_residues_contacted,
        mw, logp, tpsa, hbd          (only with --require-cns-druglike)
        top_confidence               (only with --diffdock-candidates)

    A funnel summary is printed showing how many rows survived each step.

    Columns deliberately left out: status and pose_in_box, which are
    constant across survivors by construction; the protein-wide PLIP
    counts and the raw smina score, which no filter here uses; and the
    intermediate rank columns smina_rank, diffdock_rank and
    consensus_rank_sum.

    Those rank columns are still computed internally to decide the sort
    order, so shortlist_position reflects the same consensus ordering as
    before. They are simply not written out, because a homemade composite
    score shipped as a column invites the question of how it was
    validated, which is hard to answer for an ad hoc ranking heuristic. If
    the ranking needs to be shown in a paper, recompute it from this CSV's
    affinity and confidence columns rather than shipping it precomputed.

HOW TO RUN
    Requires:
        pip install pandas
        pip install rdkit        (only for --require-cns-druglike)

    Affinity ranking only:
        python3 5_build_shortlist.py \
            --smina-plip smina_plip_report.csv \
            --out trem2_final_shortlist.csv

    With DiffDock confidence folded in for consensus ranking:
        python3 5_build_shortlist.py \
            --smina-plip smina_plip_report.csv \
            --diffdock-candidates trem2_surface1_candidates.csv \
            --out trem2_final_shortlist.csv

    With the CNS drug-likeness filter, as used for the TREM2 screen:
        python3 5_build_shortlist.py \
            --smina-plip smina_plip_report.csv \
            --diffdock-candidates trem2_surface1_candidates.csv \
            --require-cns-druglike --smiles-csv fda_approved_drugs.csv \
            --out trem2_final_shortlist.csv

OPTIONAL FILTERS, AND HOW TO DESCRIBE THEM
    --require-cns-druglike
        Off by default, so it must be passed explicitly. Keeps only drugs
        within MW <= 450, logP 1 to 4, TPSA <= 90 A^2 and at most 3
        hydrogen-bond donors. These are medicinal-chemistry desirability
        ranges adapted from the CNS MPO framework, not measured
        permeability values, and should be described that way.

    --min-interaction-types / --min-target-residues
        Both default to 1, which means off. Raising either to 2 demands
        more than a single, possibly incidental contact. Both are
        reasoned heuristics rather than thresholds established in the
        literature; say so in Methods if either is used.

    --require-specific-residue RESNUM
        Off by default. Requires contact with one named residue. For
        TREM2, Arg47 is structurally load-bearing in the CDR2/Surface1
        network (Sudom et al. 2018, J Biol Chem 293(33):12634-12646), so
        this one is literature-motivated, unlike the two above.

        The number must match the numbering inside the receptor used
        earlier in the pipeline. Predicted models are often numbered from
        1, in which case the residue number from the literature is not
        the number in that file.
"""

import argparse

import pandas as pd

# Final output column order. Anything present in --smina-plip but absent
# from this list is dropped before writing --out, whether it came from an
# older untrimmed report or is an intermediate used only while filtering.
#
# Columns that depend on an optional flag (mw/logp/tpsa/hbd from
# --require-cns-druglike, top_confidence from --diffdock-candidates) are
# skipped when that flag was not used, since they will not exist.
#
# smina_rank, diffdock_rank and consensus_rank_sum are deliberately absent:
# they are computed internally to decide sort order but are not written
# out. See the module docstring for why.
FINAL_COLUMNS = [
    "shortlist_position", "drug", "smina_minimized_affinity_kcal_mol",
    "plip_hbonds_target", "plip_hydrophobic_target",
    "plip_saltbridge_target", "plip_pistack_target",
    "plip_target_residues_contacted",
    "n_interaction_types", "n_target_residues_contacted",
    "mw", "logp", "tpsa", "hbd",
    "top_confidence",
]


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--smina-plip", required=True,
                   help="smina_plip_report.csv from smina_plip_pipeline.py")
    p.add_argument("--diffdock-candidates", default=None,
                   help="Optional: trem2_surface1_candidates.csv from "
                        "scan_diffdock_results.py, to add DiffDock "
                        "top_confidence and a consensus rank")
    p.add_argument("--require-in-box", action="store_true", default=True,
                   help="Drop rows where the pose fell outside the fixed "
                        "smina docking box (default: on)")
    p.add_argument("--allow-out-of-box", dest="require_in_box",
                   action="store_false",
                   help="Keep rows even if pose_in_box was False "
                        "(not recommended -- their smina numbers may be "
                        "unreliable)")
    p.add_argument("--min-target-hbonds", type=int, default=0,
                   help="Minimum PLIP hydrogen bonds to the target "
                        "residue group (default: 0, i.e. not required "
                        "on its own -- see --require-hbond-or-saltbridge)")
    p.add_argument("--min-target-saltbridges", type=int, default=0,
                   help="Minimum PLIP salt bridges to the target "
                        "residue group (default: 0, same note as above)")
    p.add_argument("--require-hbond-or-saltbridge", action="store_true",
                   default=True,
                   help="HARD FILTER (default: on): require "
                        "plip_hbonds_target >= 1 OR "
                        "plip_saltbridge_target >= 1. This is what "
                        "enforces actual target-residue engagement, not "
                        "just proximity. Combines with --min-target-* "
                        "via OR if you raise those above 0/1.")
    p.add_argument("--no-require-hbond-or-saltbridge",
                   dest="require_hbond_or_saltbridge", action="store_false",
                   help="Turn off the hard contact filter (not "
                        "recommended -- see script docstring for why)")

    p.add_argument("--min-interaction-types", type=int, default=1,
                   help="Require at least this many DISTINCT PLIP "
                        "interaction types to the target residue group, "
                        "counted across {hydrogen bond, hydrophobic "
                        "contact, salt bridge, pi-stack} (default: 1, "
                        "i.e. off / same as before). Set to 2 to require "
                        "more than a single, potentially incidental "
                        "contact. NOTE: this is a conservative filtering "
                        "heuristic, not a threshold established in the "
                        "literature -- describe it as such in Methods.")
    p.add_argument("--min-target-residues", type=int, default=1,
                   help="Require PLIP-detected contact with at least "
                        "this many DISTINCT residues in the target group "
                        "(default: 1, i.e. off / same as before). Set to "
                        "2 to require engagement with more than one "
                        "residue of the CDR2/Surface1 network. Same "
                        "caveat as --min-interaction-types: a reasoned "
                        "heuristic, not a literature-established cutoff.")
    p.add_argument("--require-specific-residue", type=int, default=None,
                   metavar="RESNUM",
                   help="Require PLIP-detected contact with this exact "
                        "residue number (e.g. 47 for Arg47). Off by "
                        "default. Arg47 is structurally load-bearing for "
                        "the CDR2/Surface1 network -- see Sudom et al. "
                        "2018, J Biol Chem 293(33):12634-12646 -- so "
                        "'--require-specific-residue 47' is a "
                        "literature-motivated filter, unlike the two "
                        "above.")

    p.add_argument("--require-cns-druglike", action="store_true", default=False,
                   help="Require CNS-drug-likeness thresholds (MW, logP, "
                        "TPSA, H-bond donors). Needs --smiles-csv, since "
                        "smina_plip_report.csv has no structures to "
                        "compute these from.")
    p.add_argument("--smiles-csv", default=None,
                   help="CSV mapping drug name -> SMILES (e.g. your "
                        "original DiffDock input library file). Required "
                        "if --require-cns-druglike is set.")
    p.add_argument("--smiles-drug-col", default="drug",
                   help="Column name for the drug identifier in "
                        "--smiles-csv (default: 'drug'). Must match the "
                        "'drug' values in --smina-plip exactly.")
    p.add_argument("--smiles-col", default="smiles",
                   help="Column name for the SMILES string in "
                        "--smiles-csv (default: 'smiles')")
    p.add_argument("--max-mw", type=float, default=450.0,
                   help="CNS filter: maximum molecular weight (default: 450)")
    p.add_argument("--logp-min", type=float, default=1.0,
                   help="CNS filter: minimum calculated logP (default: 1.0)")
    p.add_argument("--logp-max", type=float, default=4.0,
                   help="CNS filter: maximum calculated logP (default: 4.0)")
    p.add_argument("--max-tpsa", type=float, default=90.0,
                   help="CNS filter: maximum topological polar surface "
                        "area, in A^2 (default: 90)")
    p.add_argument("--max-hbd", type=int, default=3,
                   help="CNS filter: maximum H-bond donor count (default: 3)")
    p.add_argument("--out", default="trem2_final_shortlist.csv")
    args = p.parse_args()

    df = pd.read_csv(args.smina_plip)
    n_start = len(df)

    df = df[df["status"] == "ok"].copy()
    n_after_status = len(df)

    if args.require_in_box and "pose_in_box" in df.columns:
        df = df[df["pose_in_box"] == True]  # noqa: E712
    n_after_box = len(df)

    df = df.dropna(subset=["smina_minimized_affinity_kcal_mol"])
    n_after_affinity = len(df)

    # A non-negative affinity means local relaxation never resolved the
    # pose into anything favourable. That is a failed pose rather than a
    # weak binder, and PLIP cannot tell the difference because it only
    # checks geometry. Exclude these explicitly instead of hoping the
    # contact filter below removes them, since a bad pose can still
    # satisfy a hydrogen-bond or salt-bridge geometry test.
    df = df[df["smina_minimized_affinity_kcal_mol"] < 0]
    n_after_affinity_sign = len(df)

    contact_mask = pd.Series(False, index=df.index)
    if args.require_hbond_or_saltbridge:
        contact_mask = contact_mask | (df["plip_hbonds_target"] >= 1)
        contact_mask = contact_mask | (df["plip_saltbridge_target"] >= 1)
    if args.min_target_hbonds > 0:
        contact_mask = contact_mask | (df["plip_hbonds_target"] >= args.min_target_hbonds)
    if args.min_target_saltbridges > 0:
        contact_mask = contact_mask | (df["plip_saltbridge_target"] >= args.min_target_saltbridges)

    if args.require_hbond_or_saltbridge or args.min_target_hbonds > 0 or args.min_target_saltbridges > 0:
        df = df[contact_mask]
    n_after_contact = len(df)

    # Optional filters: interaction-type diversity, residue diversity and a
    # specific-residue requirement. All are derived from columns already in
    # the smina/PLIP report, so no extra input is needed.
    interaction_cols = ["plip_hbonds_target", "plip_hydrophobic_target",
                         "plip_saltbridge_target", "plip_pistack_target"]
    interaction_cols = [c for c in interaction_cols if c in df.columns]
    df["n_interaction_types"] = sum((df[c] >= 1).astype(int) for c in interaction_cols)

    def _residue_set(s):
        if pd.isna(s) or str(s).strip() == "":
            return set()
        return {r.strip() for r in str(s).split(";") if r.strip()}

    df["n_target_residues_contacted"] = df["plip_target_residues_contacted"].apply(
        lambda s: len(_residue_set(s)))

    if args.min_interaction_types > 1:
        df = df[df["n_interaction_types"] >= args.min_interaction_types]
    n_after_interaction_types = len(df)

    if args.min_target_residues > 1:
        df = df[df["n_target_residues_contacted"] >= args.min_target_residues]
    n_after_residue_diversity = len(df)

    if args.require_specific_residue is not None:
        target_str = str(args.require_specific_residue)
        df = df[df["plip_target_residues_contacted"].apply(
            lambda s: target_str in _residue_set(s))]
    n_after_specific_residue = len(df)

    # CNS drug-likeness. Needs SMILES from an external file, because the
    # smina/PLIP report carries no structures.
    n_after_cns = None
    if args.require_cns_druglike:
        if not args.smiles_csv:
            raise SystemExit("ERROR: --require-cns-druglike needs "
                              "--smiles-csv (a file mapping drug name -> "
                              "SMILES).")
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors
        except ImportError:
            raise SystemExit("ERROR: --require-cns-druglike needs rdkit. "
                              "Install with: pip install rdkit")

        smi_df = pd.read_csv(args.smiles_csv)
        if args.smiles_drug_col not in smi_df.columns or args.smiles_col not in smi_df.columns:
            raise SystemExit(f"ERROR: --smiles-csv must have columns "
                              f"'{args.smiles_drug_col}' and "
                              f"'{args.smiles_col}'. Found: "
                              f"{list(smi_df.columns)}")
        smi_map = dict(zip(smi_df[args.smiles_drug_col], smi_df[args.smiles_col]))

        mw_list, logp_list, tpsa_list, hbd_list, missing = [], [], [], [], []
        for drug in df["drug"]:
            smi = smi_map.get(drug)
            mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
            if mol is None:
                missing.append(drug)
                mw_list.append(None); logp_list.append(None)
                tpsa_list.append(None); hbd_list.append(None)
                continue
            mw_list.append(Descriptors.MolWt(mol))
            logp_list.append(Descriptors.MolLogP(mol))
            tpsa_list.append(Descriptors.TPSA(mol))
            hbd_list.append(Descriptors.NumHDonors(mol))

        df["mw"] = mw_list
        df["logp"] = logp_list
        df["tpsa"] = tpsa_list
        df["hbd"] = hbd_list

        if missing:
            print(f"WARNING: no SMILES found for {len(missing)} drug(s) "
                  f"in --smiles-csv -- these will FAIL the CNS filter: "
                  f"{missing}", file=__import__("sys").stderr)

        cns_mask = (
            df["mw"].notna()
            & (df["mw"] <= args.max_mw)
            & (df["logp"] >= args.logp_min) & (df["logp"] <= args.logp_max)
            & (df["tpsa"] <= args.max_tpsa)
            & (df["hbd"] <= args.max_hbd)
        )
        df = df[cns_mask]
        n_after_cns = len(df)

    df = df.sort_values("smina_minimized_affinity_kcal_mol", ascending=True)

    if args.diffdock_candidates:
        dd = pd.read_csv(args.diffdock_candidates)
        keep_cols = ["drug", "top_confidence"]
        keep_cols = [c for c in keep_cols if c in dd.columns]
        df = df.merge(dd[keep_cols], on="drug", how="left")

        df["smina_rank"] = df["smina_minimized_affinity_kcal_mol"].rank(ascending=True)
        df["diffdock_rank"] = df["top_confidence"].rank(ascending=False)
        df["consensus_rank_sum"] = df["smina_rank"] + df["diffdock_rank"]
        df = df.sort_values("consensus_rank_sum", ascending=True)

    df.insert(0, "shortlist_position", range(1, len(df) + 1))

    keep = [c for c in FINAL_COLUMNS if c in df.columns]
    dropped = [c for c in df.columns if c not in keep]
    df = df[keep]

    df.to_csv(args.out, index=False)

    print("Funnel:")
    print(f"  smina/PLIP input rows:              {n_start}")
    print(f"  status == 'ok':                     {n_after_status}")
    if args.require_in_box:
        print(f"  pose_in_box == True:                {n_after_box}")
    print(f"  smina affinity parsed:              {n_after_affinity}")
    print(f"  smina affinity < 0 (valid pose):    {n_after_affinity_sign}")
    if args.require_hbond_or_saltbridge or args.min_target_hbonds > 0 or args.min_target_saltbridges > 0:
        print(f"  PLIP target H-bond/salt-bridge hit: {n_after_contact}")
    if args.min_interaction_types > 1:
        print(f"  >= {args.min_interaction_types} distinct interaction types: {n_after_interaction_types}")
    if args.min_target_residues > 1:
        print(f"  >= {args.min_target_residues} distinct target residues:   {n_after_residue_diversity}")
    if args.require_specific_residue is not None:
        print(f"  contacts residue {args.require_specific_residue}:                  {n_after_specific_residue}")
    if n_after_cns is not None:
        print(f"  CNS drug-like (MW/logP/TPSA/HBD):    {n_after_cns}")
    if dropped:
        print(f"\nDropped {len(dropped)} column(s) not in the final "
              f"schema: {dropped}")
    print(f"\nWrote {len(df)} rows to {args.out}  (columns: {keep})")

    if len(df):
        print("\nFinal shortlist:")
        cols = ["drug", "smina_minimized_affinity_kcal_mol",
                "plip_hbonds_target", "plip_saltbridge_target",
                "n_interaction_types", "n_target_residues_contacted",
                "plip_target_residues_contacted"]
        if "top_confidence" in df.columns:
            cols.insert(1, "top_confidence")
        if "mw" in df.columns:
            cols += ["mw", "logp", "tpsa", "hbd"]
        print(df[cols].to_string(index=False))
    else:
        print("\nNo candidates survived the filters. Consider "
              "--no-require-hbond-or-saltbridge or --allow-out-of-box "
              "to see what's being excluded and why.")


if __name__ == "__main__":
    main()
