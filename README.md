# workflow
computational workflow to indentify CDR2 engaging molecules


## Layout

in-house scripts. The five pipeline stages, plus all the result files

## workflow

Each stage is a standalone script with its own `--help`. Every stage writes a
file, so any step can be re-run without repeating the ones before it.

**1. Build the drug library** — reduces FDA Drugs@FDA products to 1,040
discrete small molecules, excluding biologics, peptides, radiotracers,
mixtures and inorganics.

```bash
python3 scripts/1_build_small_molecules.py \
    --input-dir /path/to/fda_files --output fda_approved_drugs.csv
```

**2. Dock, then score each top pose** — run DiffDock separately, then:

```bash
python3 scripts/2_scan_diffdock_results.py structures/trem2_wt_af3.pdb results/diffdock_wt \
    --residue 27 --extra-residues 45 47 48 57 \
    --extra-label surface1_cdr2 --out wt_surface1_report.csv
```

**3. Keep plausible, nearby poses**

```bash
python3 scripts/3_filter_candidates.py \
    --report wt_surface1_report.csv --cutoff 4.0 --out wt_surface1_candidates.csv
```

**4. Refine and profile contacts** — needs smina, PLIP and OpenBabel:

```bash
python3 scripts/4_smina_plip_pipeline.py \
    --receptor wt_protonated.pdb --receptor-pdbqt wt.pdbqt \
    --results-root results/diffdock_wt --candidates wt_surface1_candidates.csv \
    --top-n 400 --residues 27 45 47 48 57 \
    --center 9.1 5.2 -4.8 --box-size 35 35 35 --out wt_smina_plip_report.csv
```

**5. Build the shortlist**

```bash
python3 scripts/5_build_shortlist.py \
    --smina-plip wt_smina_plip_report.csv \
    --diffdock-candidates wt_surface1_candidates.csv \
    --require-cns-druglike --smiles-csv fda_approved_drugs.csv \
    --out wt_final_shortlist.csv
```

## Three settings that must be checked for every receptor

None of these can be detected at runtime. Each one produces plausible-looking
numbers measured against the wrong thing.

1. **Residue numbering.** AlphaFold numbers each model from 1, not in UniProt
   numbering. In these models the pocket is local residues **27, 45, 47, 48,
   57**, corresponding to UniProt 47, 65, 67, 68, 77. The mutated residue is
   local 27.
2. **Box centre.** Each predicted model sits in its own coordinate frame.
   Wild type: `9.1 5.2 -4.8`. R47H: `1.4 7.9 -6.3`. A centre borrowed from
   another structure points at empty space.
3. **`--top-n`.** It defaults to 20 and silently truncates. Set it above the
   candidate count, and make sure both arms process all of theirs or the
   funnels are not comparable.

## Requirements

Two environments, because the docking and scoring stacks conflict:

```bash
# stages 1, 2, 3, 5
pip install rdkit numpy scipy pandas

# stage 4
conda create -n smina_plip -c conda-forge python=3.10 smina openbabel -y
conda activate smina_plip
pip install plip pdb2pqr lxml pandas
```

DiffDock itself is installed separately; see its own documentation.

## Data provenance

- Drug library: [FDA Drugs@FDA](https://www.fda.gov/drugs/drug-approvals-and-databases/drugsfda-data-files) (public domain)
- Structures resolved via [openFDA GSRS](https://open.fda.gov/) and [PubChem](https://pubchem.ncbi.nlm.nih.gov/)
- Receptors predicted with [AlphaFold3](https://alphafoldserver.com/) from UniProt Q9NZC2 residues 21–130

Tools: [DiffDock](https://github.com/gcorso/DiffDock),
[smina](https://sourceforge.net/projects/smina/),
[PLIP](https://github.com/pharmai/plip),
[PDB2PQR/PROPKA](https://github.com/Electrostatics/pdb2pqr),
[Open Babel](https://openbabel.org/), [RDKit](https://www.rdkit.org/).

