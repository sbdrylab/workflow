#!/usr/bin/env python3
"""
Filter the DiffDock pose report down to a candidate pool.

WHAT THIS PROGRAM DOES
    Takes the full per-drug report written by 2_scan_diffdock_results.py
    and keeps only the drugs whose top pose is worth scoring further.

    A drug is kept when both are true:
      1. top_clash is False
         The pose has no ligand atom within 1.5 A of a receptor atom, so
         it is at least physically possible.
      2. its distance to the target residue group is <= --cutoff
         At least one ligand atom comes within that many angstroms of the
         residue group measured by 2_scan_diffdock_results.py.

    Survivors are sorted by DiffDock confidence, highest first, and
    numbered with a candidate_rank column.

    This is a proximity screen only. It shows that a pose sits near the
    pocket, not that it actually interacts with it; the interaction test
    happens in the next stage of the pipeline.

INPUT
    --report    The CSV written by 2_scan_diffdock_results.py. It must
                contain a top_clash column and the distance column named
                by --dist-column.

OUTPUT
    A CSV with candidate_rank followed by these columns:

        drug, top_confidence, top_clash, top_dist_target,
        top_n_contacts, top_min_dist_any, n_heavy_atoms,
        and the --dist-column distance

    Any other columns in the input report are dropped, so the output
    schema stays the same no matter which version of the report is fed in.
    Dropped column names are printed, along with a funnel summary showing
    how many drugs went in and how many survived.

HOW TO RUN
    Requires:
        pip install pandas

    python3 3_filter_candidates.py \
        --report trem2_surface1_report.csv \
        --cutoff 4.0 \
        --out trem2_surface1_candidates.csv

    Options:
        --report        input report CSV (required)
        --cutoff        proximity cutoff in angstroms (default 8.0)
        --dist-column   distance column to filter on
                        (default top_dist_surface1_cdr2, which matches the
                        --extra-label used when running
                        2_scan_diffdock_results.py)
        --out           output CSV path
"""

import argparse

import pandas as pd

# Columns kept in the output, in this order. candidate_rank is inserted
# separately at position 0, and the --dist-column distance is appended after
# these unless it is already one of them.
KEEP_COLUMNS = [
    "drug", "top_confidence", "top_clash", "top_dist_target",
    "top_n_contacts", "top_min_dist_any", "n_heavy_atoms",
]


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--report", required=True,
                   help="Full report CSV from scan_diffdock_results.py "
                        "(e.g. trem2_surface1_report.csv)")
    p.add_argument("--dist-column", default="top_dist_surface1_cdr2",
                   help="Column holding the top pose's distance to the "
                        "target residue group (default: "
                        "top_dist_surface1_cdr2 -- matches the "
                        "--extra-label used when scan_diffdock_results.py "
                        "was run; change this if you used a different "
                        "--extra-label)")
    p.add_argument("--cutoff", type=float, default=8.0,
                   help="Proximity cutoff in Angstroms (default: 8.0 -- "
                        "this is the value actually used to build "
                        "trem2_surface1_candidates.csv)")
    p.add_argument("--out", default="trem2_surface1_candidates.csv")
    args = p.parse_args()

    df = pd.read_csv(args.report)
    n_start = len(df)

    if "top_clash" not in df.columns:
        raise SystemExit(f"ERROR: 'top_clash' column not found in {args.report}. "
                          f"Is this the output of scan_diffdock_results.py?")
    if args.dist_column not in df.columns:
        raise SystemExit(f"ERROR: '{args.dist_column}' column not found in "
                          f"{args.report}. Available columns: "
                          f"{list(df.columns)}")

    cand = df[(df["top_clash"] == False) &  # noqa: E712
              (df[args.dist_column] <= args.cutoff)].copy()
    n_after_filter = len(cand)

    cand = cand.sort_values("top_confidence", ascending=False).reset_index(drop=True)

    # Keep only the current schema, dropping any extra columns an older
    # report may carry, so the output is identical whichever version of the
    # report was supplied.
    keep = [c for c in KEEP_COLUMNS if c in cand.columns]
    if args.dist_column not in keep:
        keep.append(args.dist_column)
    dropped = [c for c in cand.columns if c not in keep]
    cand = cand[keep]

    cand.insert(0, "candidate_rank", range(1, len(cand) + 1))
    cand.to_csv(args.out, index=False)

    print("Funnel:")
    print(f"  input drugs (scan_diffdock_results.py report): {n_start}")
    print(f"  top_clash == False AND {args.dist_column} <= {args.cutoff}: {n_after_filter}")
    if dropped:
        print(f"\nDropped {len(dropped)} column(s) not in the trimmed schema: {dropped}")
    print(f"\nWrote {len(cand)} candidates to {args.out}  "
          f"(columns: {['candidate_rank'] + keep})")

    if len(cand):
        print("\nTop 10 by top_confidence:")
        cols = ["candidate_rank", "drug", "top_confidence",
                "top_dist_target", args.dist_column]
        cols = [c for c in cols if c in cand.columns]
        print(cand[cols].head(10).to_string(index=False))


if __name__ == "__main__":
    main()
