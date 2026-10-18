#!/usr/bin/env python3
"""
Build a small-molecule structure table from FDA-approved drug products.

WHAT THIS PROGRAM DOES
    Reads the FDA Drugs@FDA product files, reduces them to unique active
    ingredients, and resolves each ingredient to one well-defined chemical
    structure suitable for molecular docking.

    The pipeline:
      1. Keep only Prescription and Over-the-counter products.
      2. Split combination products into their individual ingredients.
      3. Drop ingredients that cannot be represented as a single discrete
         small molecule: undefined mixtures, radiolabeled tracers,
         biologics, therapeutic peptides, inorganic/mineral non-ligands,
         excipients and formulation materials, and diagnostic agents.
      4. Normalize salt, hydrate and simple prodrug-ester names to the
         active moiety that should actually be docked.
      5. Look up each structure in openFDA GSRS, falling back to PubChem.
      6. Standardize salts to the largest organic fragment, apply a
         molecular-weight cutoff, and deduplicate by InChIKey.

INPUT
    A folder of tab-delimited FDA Drugs@FDA files:
        Products.txt
        MarketingStatus.txt
        MarketingStatus_Lookup.txt

    Marketing status codes (from MarketingStatus_Lookup.txt):
        1 Prescription              kept
        2 Over-the-counter          kept
        3 Discontinued              excluded
        4 None (Tentative Approval) excluded
        5 For Further Manufacturing excluded

OUTPUT
    fda_approved_drugs.csv
        One row per unique structure:
        ActiveIngredient, isomeric SMILES, InChIKey, Source

    <output_stem>_review.csv
        Audit trail: canonicalized aliases, why each excluded ingredient
        was dropped, and any structure or name collisions.

    structure_lookup_cache_*.json
        Cached API lookups, reused on later runs.

    Exit code 0 on success, 2 if unresolved hard-QC collisions are found
    (affected rows are withheld from the CSV either way).

HOW TO RUN
    Requires RDKit:
        mamba install -c conda-forge rdkit

    Run:
        python3 1_build_small_molecules.py \
            --input-dir /path/to/fda_files \
            --output fda_approved_drugs.csv

    Frequently used options:
        --max-molecular-weight  MW cutoff after standardization (default 500)
        --cache                 lookup cache path; start a fresh cache after
                                changing any canonicalization rule
        --openfda-api-key       optional, also read from OPENFDA_API_KEY
        --allow-qc-warnings     write output despite QC collisions, exit 0
        --no-openfda            skip the GSRS lookup
        --no-pubchem            skip the PubChem fallback

NOTE ON THE FILTERS
    The exclusion rules below are project-scope decisions for building a
    conventional small-molecule docking library. Excluding an ingredient
    says it is unsuitable for this particular library, not that it is not
    a legitimate FDA-approved product. Every filter runs before any
    network lookup, so excluded ingredients never consume an API call.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.Chem import inchi as rd_inchi
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError as exc:  # pragma: no cover - depends on user environment
    raise SystemExit(
        "RDKit is required. Install it with conda/mamba, e.g.:\n"
        "  mamba install -c conda-forge rdkit\n"
    ) from exc


DEFAULT_INPUT_DIR = Path(
    "/Path/dafdat20260811"
)
OPENFDA_SUBSTANCE_URL = "https://api.fda.gov/other/substance.json"
PUBCHEM_PUG_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# Marketing status scope: Prescription (1) and Over-the-counter (2) only.
# Per FDA's MarketingStatus_Lookup.txt:
#   1 Prescription | 2 Over-the-counter | 3 Discontinued |
#   4 None (Tentative Approval) | 5 For Further Manufacturing Use
PRESCRIPTION_OTC_STATUS_IDS = {"1", "2"}


@dataclass
class StructureResult:
    isomeric_smiles: str
    inchikey: str
    molecular_weight: float
    source: str
    source_name: str
    source_substance_class: str = ""


@dataclass
class CachedLookup:
    status: str
    result: dict[str, Any] | None = None
    reason: str = ""


def clean_text(value: str) -> str:
    value = value.replace("\xa0", " ")
    value = re.sub(r"\s+", " ", value)
    return value.strip()


def key_text(value: str) -> str:
    return clean_text(value).upper()


# ==============================================================================
# COMBINATION-PRODUCT SPLITTING
# ==============================================================================
# FDA data separates combination ingredients with ";", "||", ",", " AND ",
# "/", and occasionally nothing at all. Order matters: the exception cases
# below are checked before generic splitting, because a naive split would
# get them wrong.

# Comma is ambiguous: most comma-containing entries are genuine multi-drug
# combinations ("CARBIDOPA, LEVODOPA"), but some use a comma to attach a
# purely descriptive qualifier to a SINGLE substance -- crystal/particle
# form, isotope grade, or synthesis type. These must NOT be split.
COMMA_QUALIFIERS = {
    "CONJUGATED", "CONJUGATED SYNTHETIC A", "CONJUGATED SYNTHETIC B",
    "ESTERIFIED", "GLACIAL", "MACROCRYSTALLINE", "MICROCRYSTALLINE",
    "ULTRAMICROSIZE", "MICROSIZE", "ULTRAMICROCRYSTALLINE", "NF", "USP",
    "CHORIONIC",
}
ISOTOPE_COMMA_PATTERN = re.compile(r"^[A-Z]{1,3}-\d+[A-Z]?$")

# Multi-drug combinations written with no separator at all. A scan of long
# separator-free entries found this to be the only genuine case, so it is
# handled as an explicit map rather than by name-matching heuristics.
KNOWN_UNPUNCTUATED_COMBOS = {
    "DOLUTEGRAVIR LAMIVUDINE TENOFOVIR ALAFENAMIDE":
        ["DOLUTEGRAVIR", "LAMIVUDINE", "TENOFOVIR ALAFENAMIDE"],
}


def _comma_split_is_qualifier(active_ingredient: str) -> bool:
    """True if this is a single substance with a comma-attached qualifier,
    not a real multi-drug combination."""
    if "," not in active_ingredient:
        return False
    parts = [clean_text(p) for p in active_ingredient.split(",")]
    if len(parts) == 2:
        tail = parts[1].upper()
        if tail in COMMA_QUALIFIERS or ISOTOPE_COMMA_PATTERN.match(parts[1]):
            return True
    return False


def split_active_ingredients(active_ingredient: str) -> list[str]:
    text = clean_text(active_ingredient)

    if key_text(text) in KNOWN_UNPUNCTUATED_COMBOS:
        return list(KNOWN_UNPUNCTUATED_COMBOS[key_text(text)])

    if "/" in text:
        parts = [p for p in text.split("/")]
        return [
            piece
            for part in parts
            for piece in split_active_ingredients(part)
            if piece
        ]

    if ";" in text or "||" in text:
        parts = re.split(r"\s*;\s*|\s*\|\|\s*", text)
        return [
            piece
            for part in parts
            for piece in split_active_ingredients(part)
            if piece
        ]

    if re.search(r"\bAND\b", text, re.IGNORECASE):
        parts = re.split(r"\s+AND\s+", text, flags=re.IGNORECASE)
        return [
            piece
            for part in parts
            for piece in split_active_ingredients(part)
            if piece
        ]

    if "," in text:
        if _comma_split_is_qualifier(text):
            return [clean_text(text)] if clean_text(text) else []
        return [clean_text(p) for p in text.split(",") if clean_text(p)]

    return [clean_text(text)] if clean_text(text) else []


# ==============================================================================
# EXCLUSION FILTERS: undefined mixtures, radiolabeled tracers, biologics
# ==============================================================================
# These run before any network lookup. Each rule was derived from real FDA
# active-ingredient data; the comments note the cases that shaped it.

UNDEFINED_MIXTURES = {
    "ESTROGENS, CONJUGATED", "ESTROGENS, ESTERIFIED",
    "ESTROGENS, CONJUGATED SYNTHETIC A", "ESTROGENS, CONJUGATED SYNTHETIC B",
    "CONJUGATED ESTROGENS",
}


def is_undefined_mixture(name: str) -> bool:
    return key_text(name) in UNDEFINED_MIXTURES


# ==============================================================================
# TREM2 R47H PRE-DOCKING ELIGIBILITY FILTERS
# ==============================================================================
# Removes FDA active ingredients that do not belong in a conventional
# small-molecule docking library. Exact-name sets are preferred over broad
# substring rules wherever a name fragment could also appear inside a
# legitimate organic drug name.

# Chemically undefined or heterogeneous materials for which one reproducible
# molecular graph cannot faithfully represent the marketed active ingredient.
NON_DISCRETE_SUBSTANCES = {
    "AIR", "AIR POLYMER-TYPE A", "AMINO ACIDS", "BIRCH TRITERPENES",
    "FISH OIL", "FISH OIL TRIGLYCERIDES", "INVERT SUGAR", "MEDICAL AIR",
    "MEDIUM CHAIN TRIGLYCERIDES", "OLIVE OIL",
    "OMEGA-3-ACID ETHYL ESTERS", "OMEGA-3-ACID ETHYL ESTERS TYPE A",
    "OMEGA-3-CARBOXYLIC ACIDS", "OMEGA-3ACID ETHYL ESTERS",
    "PYRETHRINS", "RAUWOLFIA SERPENTINA ROOT", "SAFFLOWER OIL",
    "SINECATECHINS", "SOYBEAN OIL", "VERATRUM VIRIDE ROOT",
}

# Bulk polymers/resins/excipients or formulation-type materials that may be
# patient-facing active ingredients but are not meaningful TREM2 ligands.
EXCIPIENT_FORMULATION_SUBSTANCES = {
    "CELLULOSE SODIUM PHOSPHATE", "CHOLESTYRAMINE", "COLESEVELAM",
    "COLESTIPOL", "HYDROXYPROPYL CELLULOSE", "ICODEXTRIN",
    "PATIROMER SORBITEX", "POLYETHYLENE GLYCOL 3350",
    "POLYTETRAFLUOROETHYLENE", "SEVELAMER", "SEVELAMER CARBONATE",
    "SIMETHICONE", "SIMETHICONE-CELLULOSE",
    "SODIUM POLYSTYRENE SULFONATE", "SODIUM ZIRCONIUM CYCLOSILICATE",
    "TALC", "TYLOXAPOL", "WHITE PETROLATUM",
}

# Diagnostic / imaging agents that can have perfectly valid organic structures
# and therefore would otherwise survive a structure-only filter.
DIAGNOSTIC_IMAGING_SUBSTANCES = {
    "ACETRIZOATE", "BRILLIANT BLUE G", "CALCIUM METRIZOATE", "DIATRIZOATE",
    "ETHIODIZED OIL", "EVANS BLUE", "FERUMOXIDES", "FERUMOXSIL",
    "FERUMOXYTOL", "FLUORESCEIN", "GADOBENATE DIMEGLUMINE", "GADOBUTROL",
    "GADODIAMIDE", "GADOFOSVESET TRISODIUM", "GADOPENTETATE DIMEGLUMINE",
    "GADOPICLENOL", "GADOQUATRANE", "GADOTERATE", "GADOTERIDOL",
    "GADOVERSETAMIDE", "GADOXETATE", "INDIGOTINDISULFONATE",
    "INDOCYANINE GREEN", "IOCETAMIC ACID", "IODAMIDE", "IODIPAMIDE",
    "IODIXANOL", "IODOXAMATE", "IOHEXOL", "IOMEPROL", "IOPAMIDOL",
    "IOPANOIC ACID", "IOPHENDYLATE", "IOPROMIDE", "IOTHALAMATE",
    "IOTROLAN", "IOVERSOL", "IOXAGLATE", "IOXILAN", "IPODATE",
    "ISOSULFAN BLUE", "MANGAFODIPIR TRISODIUM", "MEGLUMINE METRIZOATE",
    "METRIZAMIDE", "METRIZOATE", "METRIZOIC ACID", "PAFOLACIANINE",
    "PERFLEXANE", "PERFLUBRON", "PERFLUOROHEXYLOCTANE",
    "PERFLUOROPOLYMETHYLISOPROPYL ETHER", "PERFLUTREN", "PROPYLIODONE",
    "SULFUR HEXAFLUORIDE LIPID-TYPE A MICROSPHERES", "TRYPAN BLUE",
    "TYROPANOATE",
}

# Elemental/mineral/simple inorganic products and bulk non-ligands.  Do NOT use
# substring matching for words such as SODIUM or CALCIUM: e.g. DICLOFENAC
# SODIUM is a valid organic drug whose counterion should be stripped instead.
INORGANIC_NON_LIGAND_SUBSTANCES = {
    "ALUMINUM ACETATE", "ALUMINUM HYDROXIDE", "AMMONIUM", "AMMONIUM CHLORIDE", "ARSENIC TRIOXIDE",
    "BARIUM", "CALCIUM", "CARBON DIOXIDE", "CHROMIC", "COPPER", "CUPRIC", "FERRIC",
    "FERRIC AMMONIUM", "FERRIC CARBOXYMALTOSE", "FERRIC DERISOMALTOSE",
    "FERRIC HEXACYANOFERRATE(II)", "FERRIC OXYHYDROXIDE",
    "FERRIC PYROPHOSPHATE", "FERROUS", "GALLIUM NITRATE", "HELIUM",
    "HYDROGEN PEROXIDE", "IRON DEXTRAN", "IRON SUCROSE", "MAGNESIUM",
    "MAGNESIUM ACETATE", "MAGNESIUM ACETATE TETRAHYDRATE",
    "MAGNESIUM CARBONATE", "MAGNESIUM HYDROXIDE", "MAGNESIUM SULFATE",
    "MAGNESIUM TRISILICATE", "MANGANESE", "MANGANESE CHLORIDE TETRAHYDRATE",
    "NITRIC OXIDE", "NITROGEN", "NITROUS OXIDE", "OXYGEN", "PHOSPHORIC ACID",
    "POTASSIUM", "POTASSIUM ACETATE", "POTASSIUM BITARTRATE",
    "POTASSIUM PERCHLORATE", "POTASSIUM PHOSPHATE", "SELENIOUS ACID",
    "SELENIUM SULFIDE", "SODIUM", "SODIUM ACETATE", "SODIUM CARBONATE",
    "SODIUM GLYCEROPHOSPHATE", "SODIUM MONOFLUOROPHOSPHATE", "SODIUM NITRITE",
    "SODIUM PHOSPHATE", "SODIUM SUCCINATE", "SODIUM SULFATE",
    "SODIUM THIOSULFATE", "SULFUR", "ZINC", "ZINC ACETATE",
    "ZINC SULFATE HEPTAHYDRATE",
    # simple mineral salts/non-ligands observed in the output audit
    "CALCIUM CARBONATE", "CALCIUM CITRATE", "CALCIUM GLUCONATE",
    "FERROUS GLUCONATE", "FERRIC CITRATE", "LANTHANUM CARBONATE",
    "LITHIUM", "LITHIUM CARBONATE", "LITHIUM CITRATE",
    "MAGNESIUM CITRATE", "MAGNESIUM GLUCONATE",
    "POTASSIUM CITRATE", "POTASSIUM GLUCONATE",
    "SODIUM BICARBONATE", "SODIUM CITRATE", "SODIUM GLUCONATE",
    "ZINC GLUCONATE",
}

# Labels that are not usable active-ingredient identities by themselves.
INVALID_GENERIC_LABELS = {
    "ANHYDROUS", "BETA", "DIBASIC", "DIMETHYL", "HEPTAHYDRATE", "METHYL",
    "MONOBASIC", "MONOHYDRATE", "MONOMETHYL",
}

# Therapeutic peptides/peptidomimetics too large or mechanistically unsuitable
# for the intended conventional small-molecule library.  The suffix rule is
# deliberately narrow and backed by exceptions rather than treating every
# peptide-like substring as biologic.
PEPTIDE_EXACT_NAMES = {
    "ANGIOTENSIN II", "ANGIOTENSIN II ACETATE",
    "BREMELANOTIDE", "EFTILAGIMOD ALPHA", "ELOBIXIBAT",  # harmless if absent
    "EPTIFIBATIDE", "EXENATIDE", "GLUCAGON", "ICATIBANT", "LANREOTIDE",
    "LEUPROLIDE", "LINACLOTIDE", "LIRAGLUTIDE", "LIXISENATIDE",
    "OCTREOTIDE", "PASIREOTIDE", "PLECANATIDE", "PRAMLINTIDE",
    "SEMAGLUTIDE", "SETMELANOTIDE", "TEDUGLUTIDE", "TIRZEPATIDE",
    "TRIPTORELIN", "VOSORITIDE", "ZICONOTIDE",
}
PEPTIDE_TIDE_FALSE_POSITIVE_EXCEPTIONS: set[str] = set()


# metal/coordination products for which the metal is part of the drug's
# molecular identity (or the marketed active is a metal-carbohydrate/metal-
# ligand complex).  These must NOT be converted to the "largest organic
# fragment", because doing so creates a different molecule.
INTEGRAL_METAL_NAME_EXCLUSIONS = {
    "AURANOFIN", "BISMUTH SUBCITRATE", "BISMUTH SUBSALICYLATE",
    "CARBOPLATIN", "CISPLATIN", "COPPER HISTIDINATE", "FERRIC MALTOL",
    "OXALIPLATIN", "SILVER SULFADIAZINE", "SODIUM FERRIC GLUCONATE COMPLEX",
    "SODIUM NITROPRUSSIDE", "SODIUM AUROTHIOMALATE",
}

# Metals that conventional organic-ligand docking preparation should not
# silently delete from a connected molecular graph.  Ordinary salts such as
# diclofenac sodium are normally represented as disconnected fragments and are
# therefore not rejected by the structural test below.
METAL_ATOMIC_NUMBERS = {
    3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
    37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
    55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
    72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
}


def is_integral_metal_name_exclusion(name: str) -> bool:
    return key_text(strip_nonidentity_qualifier(name)) in INTEGRAL_METAL_NAME_EXCLUSIONS


def is_non_discrete_substance(name: str) -> bool:
    upper = key_text(strip_nonidentity_qualifier(name))
    if upper in NON_DISCRETE_SUBSTANCES:
        return True
    # Conservative generic patterns for obviously heterogeneous oils/extracts.
    if upper.endswith(" OIL") or upper.endswith(" OILS"):
        return True
    if " TRIGLYCERIDES" in upper:
        return True
    return False


def is_excipient_or_formulation(name: str) -> bool:
    return key_text(strip_nonidentity_qualifier(name)) in EXCIPIENT_FORMULATION_SUBSTANCES


def is_diagnostic_or_imaging(name: str) -> bool:
    upper = key_text(strip_nonidentity_qualifier(name))
    if upper in DIAGNOSTIC_IMAGING_SUBSTANCES:
        return True
    # Gadolinium contrast INNs are highly distinctive and all are outside this
    # project's therapeutic docking scope.
    if upper.startswith("GADO"):
        return True
    return False


def is_inorganic_or_non_ligand(name: str) -> bool:
    return key_text(strip_nonidentity_qualifier(name)) in INORGANIC_NON_LIGAND_SUBSTANCES


def is_invalid_generic_label(name: str) -> bool:
    return key_text(strip_nonidentity_qualifier(name)) in INVALID_GENERIC_LABELS


def is_therapeutic_peptide(name: str) -> bool:
    upper = key_text(strip_nonidentity_qualifier(name))
    if upper in PEPTIDE_EXACT_NAMES:
        return True
    tokens = re.split(r"[\s,;]+", upper)
    for tok in tokens:
        tok_clean = tok.split("-")[0]
        if tok_clean.endswith("TIDE") and tok_clean not in PEPTIDE_TIDE_FALSE_POSITIVE_EXCEPTIONS:
            return True
        # Common peptide endings not captured by -tide.
        if tok_clean.endswith("RELIN") and tok_clean in {"GOSERELIN", "HISTRELIN", "NAFARELIN", "TRIPTORELIN"}:
            return True
    return False



def trem2_exclusion_reason(name: str) -> str | None:
    """Return a project-scope exclusion reason for a raw OR canonical name.

    v6 deliberately calls this both before structure lookup and again after
    active-moiety canonicalization.  This prevents a transformation from
    creating a forbidden/generic identity such as METHYL or DIMETHYL after the
    raw ingredient already passed the pre-filter.
    """
    if is_undefined_mixture(name):
        return "undefined mixture"
    if is_non_discrete_substance(name):
        return "non-discrete/heterogeneous substance"
    if is_multi_organic_component_exclusion(name):
        return "defined multi-organic-component ingredient"
    if is_integral_metal_name_exclusion(name):
        return "integral metal/coordination drug not suitable for conventional docking"
    if is_inorganic_or_non_ligand(name):
        return "inorganic/mineral/non-ligand"
    if is_diagnostic_or_imaging(name):
        return "diagnostic/imaging agent"
    if is_excipient_or_formulation(name):
        return "excipient/formulation material"
    if is_invalid_generic_label(name):
        return "invalid/generic active-ingredient label"
    # reporting priority: when categories overlap, report the most specific
    # scientifically informative reason first.
    if is_radiolabeled_tracer(name):
        return "radiolabeled diagnostic/radiopharmaceutical"
    if is_biologic(name):
        return "biologic"
    if is_therapeutic_peptide(name):
        return "therapeutic peptide/non-small-molecule"
    return None



# ==============================================================================
# ACTIVE-MOIETY CANONICALIZATION
# ==============================================================================
# The library holds molecular active moieties, not dosage-form salt names.
# This layer prevents two failure modes:
#   1) keeping a large organic counter-ion (notably PAMOATE) just because it
#      is bigger than the drug cation, and
#   2) merging distinct drugs because a name lookup returned one structure.
#
# Exact overrides handle complexes, prodrugs and legacy names where a blind
# suffix rule would be unsafe. Unambiguous salt and hydrate suffixes are then
# stripped iteratively, e.g. FORMOTEROL FUMARATE DIHYDRATE -> FORMOTEROL.

NONIDENTITY_COMMA_QUALIFIERS = {
    "USP", "NF", "MACROCRYSTALLINE", "MICROCRYSTALLINE", "MICROSIZE",
    "ULTRAMICROSIZE", "ULTRAMICROCRYSTALLINE", "GLACIAL",
}


def strip_nonidentity_qualifier(name: str) -> str:
    """Remove FDA descriptive comma qualifiers that do not change identity."""
    text = clean_text(name)
    if "," not in text:
        return text
    head, tail = text.rsplit(",", 1)
    if key_text(tail) in NONIDENTITY_COMMA_QUALIFIERS:
        return clean_text(head)
    return text


# Defined multi-organic-component actives for which choosing a single "largest
# fragment" would misrepresent the approved ingredient.  For a one-ligand-per-
# row TREM2 docking library these are safer to exclude than to silently choose
# one component.  Dimenhydrinate, for example, contains diphenhydramine AND
# 8-chlorotheophylline as a defined salt pair.
MULTI_ORGANIC_COMPONENT_EXCLUSIONS = {
    "DIMENHYDRINATE",
}


def is_multi_organic_component_exclusion(name: str) -> bool:
    return key_text(strip_nonidentity_qualifier(name)) in MULTI_ORGANIC_COMPONENT_EXCLUSIONS


# Exact redirects to the molecular active moiety used for docking.  These are
# intentionally explicit because the relationship is not safely inferable from
# a generic salt suffix alone.
ACTIVE_MOIETY_OVERRIDES = {
    "DEHYDRATED ALCOHOL": "ALCOHOL",
    "TOCOPHEROL ACETATE": "ALPHA-TOCOPHEROL ACETATE",
    "AMINOPHYLLINE": "THEOPHYLLINE",
    "CALCIUM OXYBATE": "OXYBATE",
    "MAGNESIUM OXYBATE": "OXYBATE",
    "POTASSIUM OXYBATE": "OXYBATE",
    "SODIUM OXYBATE": "OXYBATE",
    "DIVALPROEX SODIUM": "VALPROIC ACID",
    "VALPROATE SODIUM": "VALPROIC ACID",
    "GLYCOPYRROLATE": "GLYCOPYRRONIUM",
    "PENICILLIN G BENZATHINE": "PENICILLIN G",
    "PENICILLIN G PROCAINE": "PENICILLIN G",
    "PENICILLIN G POTASSIUM": "PENICILLIN G",
    "PENICILLIN G SODIUM": "PENICILLIN G",
    "MICONAZOLE NITRATE": "MICONAZOLE",
    # corrections exposed by full-output review
    "ETHACRYNATE": "ETHACRYNIC ACID",
    "DICLOFENAC EPOLAMINE": "DICLOFENAC",
    "DIAZOXIDE CHOLINE": "DIAZOXIDE",
    "ELTROMBOPAG CHOLINE": "ELTROMBOPAG",
    "PERINDOPRIL ERBUMINE": "PERINDOPRIL",
    "CHOLINE FENOFIBRATE": "FENOFIBRIC ACID",
    "SODIUM ASCORBATE": "ASCORBIC ACID",
    "SODIUM PHENYLACETATE": "PHENYLACETIC ACID",
    "SODIUM PHENYLBUTYRATE": "PHENYLBUTYRIC ACID",
    "CHOLINE SALICYLATE": "SALICYLIC ACID",
    "SODIUM SALICYLATE": "SALICYLIC ACID",
    "TROLAMINE SALICYLATE": "SALICYLIC ACID",
    # targeted corrections from the 100-drug audit
    "AMMONIUM LACTATE": "LACTIC ACID",
    "ALISKIREN HEMIFUMARATE": "ALISKIREN",
}

# Query names used when the desired canonical display name is not the most
# reliable compound-name query.  PubChem/GSRS recognize OXYBATE, but the parent
# chemical identity is 4-hydroxybutyric acid; querying that parent avoids
# keeping a metal counter-ion.
ACTIVE_MOIETY_QUERY_OVERRIDES = {
    "OXYBATE": "4-HYDROXYBUTYRIC ACID",
}

# safety valve: these FDA ingredient names end in a word that can also be a
# counterion name, but here that word is part of the drug molecule itself.
# They must never be shortened by suffix stripping.
COUNTERION_STRIP_PROTECTED_EXACT = {
    "DIMETHYL FUMARATE",
    "MONOMETHYL FUMARATE",
    "METHYL SALICYLATE",
}

# Stereochemistry-sensitive identities that are unsafe to resolve by an
# ambiguous name search.  PubChem CID lookup gives a deterministic structure.
# Keep these CIDs version-controlled and review them if PubChem changes the
# corresponding record semantics.
PUBCHEM_CID_OVERRIDES = {
    "ARFORMOTEROL": 3083544,
    "FORMOTEROL": 3410,
    "DEXMETHYLPHENIDATE": 154101,
    "METHYLPHENIDATE": 4158,
    "LEVONORGESTREL": 13109,
    "NORGESTREL": 16051930,
}

# Suffixes that are unambiguously formulation counter-ions/hydrates when they
# occur as the FINAL words of an ingredient name.  Ambiguous covalent ester
# suffixes such as ACETATE, PROPIONATE, DECANOATE, CYPIONATE and VALERATE are
# deliberately excluded and remain governed by explicit maps below.
UNAMBIGUOUS_COUNTERION_SUFFIXES = tuple(sorted({
    "LAURYL SULFATE",
    "TETRAHYDROCHLORIDE", "TRIHYDROCHLORIDE", "DIHYDROCHLORIDE", "MONOHYDROCHLORIDE",
    "HYDROCHLORIDE", "HYDROBROMIDE",
    "HEMIHYDRATE", "MONOHYDRATE", "DIHYDRATE", "TRIHYDRATE", "HYDRATE",
    "DIPOTASSIUM", "DISODIUM",
    "BENZOATE", "BESYLATE", "BESILATE", "DIBESYLATE", "MALEATE", "DIMALEATE", "MALATE", "FUMARATE", "DIFUMARATE",
    "TARTRATE", "BITARTRATE", "CITRATE", "MESYLATE", "DIMESYLATE", "MESILATE", "TOSYLATE", "TOSILATE",
    "SULFATE", "SULPHATE", "BROMIDE", "CHLORIDE", "IODIDE", "HYCLATE",
    "GLUCONATE", "LACTATE", "SALICYLATE", "CAMSYLATE", "CAMSILATE",
    "EDISYLATE", "ISETHIONATE", "NAPSYLATE", "PAMOATE",
    "POLYGALACTURONATE", "GLUCEPTATE", "HIPPURATE",
    "MEGLUMINE", "TROMETHAMINE", "DIOLAMINE", "OLAMINE", "ARGININE", "LYSINE",
    "TANNATE", "POLISTIREX", "XINAFOATE", "BISULFATE", "MONOSULFATE",
    "ASPARTATE", "SACCHARATE",
    "CALCIUM", "MAGNESIUM", "POTASSIUM", "SODIUM", "ZINC", "LITHIUM",
}, key=lambda x: (-len(x), x)))


def strip_unambiguous_counterion_suffixes(name: str) -> str:
    """Conservatively strip terminal counterion/hydrate words.

    v6 refuses a strip when the full FDA name is a protected chemical identity
    or when stripping would create an invalid/generic/inorganic name.  This
    specifically prevents failures such as DIMETHYL FUMARATE -> DIMETHYL and
    METHYL SALICYLATE -> METHYL.
    """
    text = clean_text(name)
    if key_text(text) in COUNTERION_STRIP_PROTECTED_EXACT:
        return text
    changed = True
    while changed and text:
        changed = False
        upper = key_text(text)
        for suffix in UNAMBIGUOUS_COUNTERION_SUFFIXES:
            token = " " + suffix
            if upper.endswith(token) and len(text) > len(token):
                candidate = clean_text(text[:-len(token)])
                if not candidate:
                    return text
                # Never let name-standardization manufacture a generic/mineral
                # identity that would not itself be a valid docking drug.
                if (
                    is_invalid_generic_label(candidate)
                    or is_inorganic_or_non_ligand(candidate)
                    or is_integral_metal_name_exclusion(candidate)
                ):
                    return text
                text = candidate
                if key_text(text) in COUNTERION_STRIP_PROTECTED_EXACT:
                    return text
                changed = True
                break
    return text


# Identities where stereochemistry or racemate status must stay distinct.
# Queried PubChem-first, since a name search collapsed these pairs onto one
# structure.
STEREOCHEMISTRY_SENSITIVE_CANONICAL_NAMES = set(PUBCHEM_CID_OVERRIDES)


def canonical_active_moiety_name(ingredient: str) -> str:
    """Return one clean molecular active-moiety name for docking/output."""
    name = strip_nonidentity_qualifier(ingredient)
    key = key_text(name)

    # Existing explicit prodrug/counterion maps are defined later in this file;
    # look them up dynamically at call time so this function can live with the
    # other pre-filter helpers.
    prodrug_map = globals().get("PRODRUG_ESTER_MAP", {})
    counterion_map = globals().get("SIMPLE_COUNTERION_MAP", {})
    if key in prodrug_map:
        return clean_text(prodrug_map[key])
    if key in counterion_map:
        return clean_text(counterion_map[key])
    if key in ACTIVE_MOIETY_OVERRIDES:
        return clean_text(ACTIVE_MOIETY_OVERRIDES[key])

    return strip_unambiguous_counterion_suffixes(name)


def canonical_lookup_name(ingredient: str) -> str:
    canonical = canonical_active_moiety_name(ingredient)
    return ACTIVE_MOIETY_QUERY_OVERRIDES.get(key_text(canonical), canonical)

# Radiolabeled diagnostic tracers. Covers the three isotope-notation styles
# that appear in real data; each one is needed to catch entries the others
# miss:
#   hyphenated:      "SODIUM IODIDE I-131", "LUTETIUM LU-177 ..."
#   space-separated: "LUTETIUM LU 177 DOTATATE"
#   no separator:    "TECHNETIUM TC99M ALBUMIN AGGREGATED"
ISOTOPE_PATTERN = re.compile(
    r"\b[A-Z]{1,2}-\d{2,3}[A-Z]?\b"      # hyphenated: I-131, LU-177
    r"|\b[A-Z]{1,2} \d{2,3}[A-Z]?\b"     # space-separated: LU 177
    r"|\bTC\d{2,3}M\b"                   # no-separator: TC99M
    r"|\bTECHNETIUM\b"
    r"|\bKIT$",
    re.IGNORECASE,
)


def is_radiolabeled_tracer(name: str) -> bool:
    return bool(ISOTOPE_PATTERN.search(name.upper()))


# USAN reserves the "-ase" suffix for enzymes, but a few small-molecule
# names end the same way. No false positives were found in the current
# data; this set remains as a safety valve for future entries.
ENZYME_SUFFIX_FALSE_POSITIVE_EXCEPTIONS: set[str] = set()


def is_biologic(name: str) -> bool:
    """
    Detect biologics by name. Plain substring/suffix containment is used
    deliberately for most stems here, NOT word-boundary (\\b) or start/
    end-anchored regexes -- repeated testing against real FDA data found
    that concatenated INN names defeat boundary/anchor checks:
        "ETANERCEPT"              -- no boundary before "CEPT"
        "ROPEGINTERFERON"         -- no boundary before "INTERFERON"
        "CHORIOGONADOTROPIN ALFA" -- no boundary before "GONADOTROPIN"
        "LONAPEGSOMATROPIN-TCGD"  -- no boundary before "SOMATROPIN", and
                                     "PEG" isn't at the string START either
        "EFBEMALENOGRASTIM"       -- not a FIL/SAR/PEGFIL-prefixed name
        "CIPAGLUCOSIDASE ALFA-ATGA" -- "-ASE" not at the very end of the
                                     whole string (a biosimilar suffix
                                     follows it)
    These stems are long/distinctive enough that plain substring
    containment does not introduce false positives in practice (verified
    empirically) -- with ONE confirmed exception: a blind "TROPIN"
    substring WOULD false-positive on real small molecules ATROPINE,
    BENZTROPINE, HOMATROPINE, ANISOTROPINE (unrelated tropane-alkaloid
    "-tropin/e" root) -- so TROPIN-family hormones are kept as an
    EXPLICIT list, not a general substring rule.
    """
    upper = key_text(name)

    # Antibody therapeutics: monoclonal antibodies, biosimilars and
    # antibody-drug conjugates, where a payload name follows the antibody
    # (BRENTUXIMAB VEDOTIN). Checking every whitespace-separated token rather
    # than just the string's end covers conjugates without listing payloads.
    tokens = re.split(r"[\s,;]+", upper)
    for tok in tokens:
        tok_clean = tok.rstrip(".")
        if re.search(r"MAB(-[A-Za-z]{3,6})?$", tok_clean):
            return True

    # Fusion proteins, cytokines, insulins, gonadotropins, and other
    # protein-hormone stems long/distinctive enough for safe substring use.
    substring_stems = [
        "CEPT", "INTERFERON", "INSULIN",
        "FOLLITROPIN", "UROFOLLITROPIN", "EPOETIN",
        "GONADOTROPIN", "SOMATROPIN", "THYROTROPIN", "CORTICOTROPIN",
        "COSYNTROPIN",
        "MENOTROPINS",
        "ALBUMIN HUMAN", "FIBRINOGEN", "THROMBIN", "PROTHROMBIN COMPLEX",
        "ANTIHEMOPHILIC FACTOR",
    ]
    if any(stem in upper for stem in substring_stems):
        return True

    # "-STIM" (thrombopoietin/GM-CSF/G-CSF receptor agonists) subsumes the
    # older narrower "GRASTIM"-only check and also catches ROMIPLOSTIM,
    # SARGRAMOSTIM (neither contains the literal substring "GRASTIM").
    for tok in tokens:
        tok_clean = tok.split("-")[0]
        if tok_clean.endswith("STIM"):
            return True

    # Interleukin-class biologics (USAN stems "-leukin", "-kinra").
    for tok in tokens:
        tok_clean = tok.split("-")[0]
        if tok_clean.endswith("LEUKIN") or tok_clean.endswith("KINRA"):
            return True

    # Recombinant/biologic enzymes: USAN reserves "-ase" for enzymes.
    # Checking any token (not just the string's end) catches biosimilar-
    # suffixed names like "CIPAGLUCOSIDASE ALFA-ATGA". Verified
    # empirically against real data: zero false positives.
    for tok in tokens:
        tok_clean = tok.split("-")[0]
        if (
            len(tok_clean) > 3
            and tok_clean.endswith("ASE")
            and tok_clean not in ENZYME_SUFFIX_FALSE_POSITIVE_EXCEPTIONS
        ):
            return True

    if re.search(r"\bFACTOR (VIII|IX|VIIA)\b", upper):
        return True

    # Botulinum toxin, in both the concatenated INN style
    # (ONABOTULINUMTOXINA) and the older spaced format (Botulinum Toxin Type
    # B). A check for only the concatenated form misses the latter.
    if "BOTULINUMTOXIN" in upper or "BOTULINUM TOXIN" in upper:
        return True

    # PEGylated biologics/peptides where "PEG" is a genuine prefix at the
    # START of the name (e.g. "PEGCETACOPLAN") -- catches PEGylated drugs
    # with no other recognizable protein stem. Mid-string "PEG" on a
    # fused name (e.g. "LONAPEGSOMATROPIN") is instead caught above via
    # the SOMATROPIN substring check.
    if re.match(r"^PEG[A-Z]", upper):
        return True

    # Antisense oligonucleotide and siRNA therapeutics: large synthetic
    # oligonucleotides, not docking candidates. Matched per token because real
    # entries append a salt suffix after the oligo stem (FOMIVIRSEN SODIUM),
    # which would otherwise hide the stem at the end of the string.
    for tok in tokens:
        tok_clean = tok.split("-")[0]
        if re.search(r"[A-Z]+SEN$", tok_clean) or tok_clean.endswith("SIRAN"):
            return True

    return False


# ==============================================================================
# SIMPLE ESTER/PRODRUG NAME REDIRECTION
# ==============================================================================
# Fragment standardization strips ionic salts, but cannot resolve a
# covalently bonded ester prodrug to its free parent: splitting separates
# disconnected ions and never breaks a covalent bond. This map redirects
# known simple prodrug esters to query the active parent instead. The
# original ingredient name is still written to the output CSV.
#
# Esters where the modified form is itself the dominant active entity
# (clobetasol propionate, medroxyprogesterone acetate, ulipristal acetate,
# triamcinolone acetonide) are deliberately absent from this map and are
# looked up under their own names.
PRODRUG_ESTER_MAP = {
    "HYDROCORTISONE SODIUM SUCCINATE": "HYDROCORTISONE",
    "HYDROCORTISONE SODIUM PHOSPHATE": "HYDROCORTISONE",
    "HYDROCORTISONE ACETATE": "HYDROCORTISONE",
    "HYDROCORTISONE CYPIONATE": "HYDROCORTISONE",
    "CORTISONE ACETATE": "CORTISONE",
    "PREDNISOLONE ACETATE": "PREDNISOLONE",
    "PREDNISOLONE SODIUM PHOSPHATE": "PREDNISOLONE",
    "METHYLPREDNISOLONE ACETATE": "METHYLPREDNISOLONE",
    "METHYLPREDNISOLONE SODIUM SUCCINATE": "METHYLPREDNISOLONE",
    "DEXAMETHASONE ACETATE": "DEXAMETHASONE",
    "DEXAMETHASONE SODIUM PHOSPHATE": "DEXAMETHASONE",
    "BETAMETHASONE ACETATE": "BETAMETHASONE",
    "BETAMETHASONE SODIUM PHOSPHATE": "BETAMETHASONE",
    "FLUDROCORTISONE ACETATE": "FLUDROCORTISONE",
    "PARAMETHASONE ACETATE": "PARAMETHASONE",
    "FLUOROMETHOLONE ACETATE": "FLUOROMETHOLONE",
    "CHLORAMPHENICOL SODIUM SUCCINATE": "CHLORAMPHENICOL",
    "TESTOSTERONE ENANTHATE": "TESTOSTERONE",
    "TESTOSTERONE CYPIONATE": "TESTOSTERONE",
    "TESTOSTERONE PROPIONATE": "TESTOSTERONE",
    "TESTOSTERONE UNDECANOATE": "TESTOSTERONE",
    "ESTRADIOL VALERATE": "ESTRADIOL",
    "ESTRADIOL CYPIONATE": "ESTRADIOL",
    "ESTRADIOL ACETATE": "ESTRADIOL",
    "NANDROLONE DECANOATE": "NANDROLONE",
    "NANDROLONE PHENPROPIONATE": "NANDROLONE",
    "DESOXYCORTICOSTERONE ACETATE": "DESOXYCORTICOSTERONE",
    "DESOXYCORTICOSTERONE PIVALATE": "DESOXYCORTICOSTERONE",
    "HALOPERIDOL DECANOATE": "HALOPERIDOL",
    "FLUPHENAZINE DECANOATE": "FLUPHENAZINE",
    "FLUPHENAZINE ENANTHATE": "FLUPHENAZINE",
    "NORETHINDRONE ACETATE": "NORETHINDRONE",  # rapidly deacetylated in vivo (unlike MPA)
    "ABIRATERONE ACETATE": "ABIRATERONE",       # inert prodrug mask; abiraterone is the active CYP17 inhibitor
    "ESLICARBAZEPINE ACETATE": "ESLICARBAZEPINE",
    "ERYTHROMYCIN ETHYLSUCCINATE": "ERYTHROMYCIN",
    "ETOPOSIDE PHOSPHATE": "ETOPOSIDE",
    "FLUDARABINE PHOSPHATE": "FLUDARABINE",
    "CLINDAMYCIN PHOSPHATE": "CLINDAMYCIN",
    "OMACETAXINE MEPESUCCINATE": "OMACETAXINE",
    # additional clearly defined prodrug -> active-parent mappings
    "ADEFOVIR DIPIVOXIL": "ADEFOVIR",
    "ARIPIPRAZOLE LAUROXIL": "ARIPIPRAZOLE",
}

# Simple counter-ion salts using PHOSPHATE/ACETATE/SUCCINATE/etc as a
# plain solubilizing salt (not a covalent ester prodrug) -- common for
# modern small-molecule and peptide drugs. Covalently identical to the
# free parent minus a genuinely ionic counter-ion; listed explicitly to
# fill gaps where the exact salt-form name might not be indexed under
# that name in GSRS/PubChem.
SIMPLE_COUNTERION_MAP = {
    "CODEINE PHOSPHATE": "CODEINE", "ANILERIDINE PHOSPHATE": "ANILERIDINE",
    "ANTAZOLINE PHOSPHATE": "ANTAZOLINE", "CHLOROQUINE PHOSPHATE": "CHLOROQUINE",
    "DISOPYRAMIDE PHOSPHATE": "DISOPYRAMIDE", "HISTAMINE PHOSPHATE": "HISTAMINE",
    "OSELTAMIVIR PHOSPHATE": "OSELTAMIVIR", "PRIMAQUINE PHOSPHATE": "PRIMAQUINE",
    "SITAGLIPTIN PHOSPHATE": "SITAGLIPTIN", "SONIDEGIB PHOSPHATE": "SONIDEGIB",
    "RUXOLITINIB PHOSPHATE": "RUXOLITINIB", "INFIGRATINIB PHOSPHATE": "INFIGRATINIB",
    "DEURUXOLITINIB PHOSPHATE": "DEURUXOLITINIB", "OSILODROSTAT PHOSPHATE": "OSILODROSTAT",
    "LENIOLISIB PHOSPHATE": "LENIOLISIB", "CARVEDILOL PHOSPHATE": "CARVEDILOL",
    "AMIFAMPRIDINE PHOSPHATE": "AMIFAMPRIDINE", "TEDIZOLID PHOSPHATE": "TEDIZOLID",
    "DESVENLAFAXINE SUCCINATE": "DESVENLAFAXINE", "DOXYLAMINE SUCCINATE": "DOXYLAMINE",
    "FROVATRIPTAN SUCCINATE": "FROVATRIPTAN", "LASMIDITAN SUCCINATE": "LASMIDITAN",
    "LOXAPINE SUCCINATE": "LOXAPINE", "METOPROLOL SUCCINATE": "METOPROLOL",
    "MOBOCERTINIB SUCCINATE": "MOBOCERTINIB", "PRUCALOPRIDE SUCCINATE": "PRUCALOPRIDE",
    "RIBOCICLIB SUCCINATE": "RIBOCICLIB", "SOLIFENACIN SUCCINATE": "SOLIFENACIN",
    "SUMATRIPTAN SUCCINATE": "SUMATRIPTAN", "TAFENOQUINE SUCCINATE": "TAFENOQUINE",
    "BAZEDOXIFENE ACETATE": "BAZEDOXIFENE", "CASPOFUNGIN ACETATE": "CASPOFUNGIN",
    "REZAFUNGIN ACETATE": "REZAFUNGIN", "LEFAMULIN ACETATE": "LEFAMULIN",
    "FLECAINIDE ACETATE": "FLECAINIDE", "GUANABENZ ACETATE": "GUANABENZ",
    "MAFENIDE ACETATE": "MAFENIDE",
}

# Ester forms deliberately kept separate from any "parent", because
# redirecting them would be chemically wrong. Not used by lookup logic;
# listed so these are not mistakenly added to the maps above.
PROTECTED_DISTINCT_ACTIVITY_DO_NOT_REDIRECT = {
    "CLOBETASOL PROPIONATE", "BETAMETHASONE DIPROPIONATE", "BETAMETHASONE VALERATE",
    "BECLOMETHASONE DIPROPIONATE", "ALCLOMETASONE DIPROPIONATE",
    "HALOBETASOL PROPIONATE", "DIFLORASONE DIACETATE",
    "FLUOCINOLONE ACETONIDE", "TRIAMCINOLONE ACETONIDE",
    "TRIAMCINOLONE DIACETATE", "TRIAMCINOLONE HEXACETONIDE",
    "FLUTICASONE PROPIONATE", "FLUTICASONE FUROATE", "MOMETASONE FUROATE",
    "HYDROCORTISONE BUTYRATE", "HYDROCORTISONE VALERATE",
    "FLUMETHASONE PIVALATE", "CLOCORTOLONE PIVALATE",
    "MEDROXYPROGESTERONE ACETATE", "MEGESTROL ACETATE",
    "SEGESTERONE ACETATE", "ULIPRISTAL ACETATE",
    "ISOSORBIDE DINITRATE", "ISOSORBIDE MONONITRATE",
}


def canonical_query_name(ingredient: str) -> str:
    """Return the molecular active-moiety name actually queried in v6."""
    return canonical_lookup_name(ingredient)


POLYMER_NAME_RE = re.compile(
    r"\b("
    r"POLYETHYLENE GLYCOL|POLYETHYLENE|POLYVINYL|POLYOXYL|POLYOXYETHYLENE|"
    r"POLYSORBATE|POLOXAMER|CARBOMER|POVIDONE|CROSPOVIDONE|"
    r"CELLULOSE|HYPROMELLOSE|DEXTRAN|STARCH|ALGINATE|GELATIN|PECTIN|"
    r"XANTHAN GUM|GUAR GUM|DIMETHICONE|SILICONE|SILOXANE|"
    r"HEPARIN|ENOXAPARIN|DALTEPARIN|TINZAPARIN|"
    r"PEPTIDE|POLYPEPTIDE|PROTEIN|ANTIBODY|IMMUNOGLOBULIN|"
    r"OLIGONUCLEOTIDE|NUCLEOTIDE|NUCLEIC ACID|RNA|DNA"
    r")\b",
    flags=re.IGNORECASE,
)


def looks_polymer_like(name: str) -> bool:
    return bool(POLYMER_NAME_RE.search(name))


def load_tsv(path: Path) -> list[dict[str, str]]:
    # The provided FDA files contain Latin-1 bytes in some rows.
    with path.open(newline="", encoding="latin-1") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def selected_status_ids(args: argparse.Namespace) -> set[str]:
    # Fixed to Prescription (1) + Over-the-counter (2) only for this
    # project's scope -- see PRESCRIPTION_OTC_STATUS_IDS and the module
    # docstring for the full MarketingStatusID -> description mapping.
    return set(PRESCRIPTION_OTC_STATUS_IDS)


def load_product_ingredients(
    input_dir: Path,
    args: argparse.Namespace,
) -> tuple[list[tuple[str, str]], list[dict[str, str]]]:
    """Load in-scope FDA ingredients and preserve a complete prelookup audit.

    Returns:
        product_rows:
            Unique (DrugName, FDA ingredient) pairs that survive the v6/v7
            prelookup chemical eligibility rules and therefore proceed to
            structure lookup.
        prelookup_audit:
            One terminal PRELOOKUP_EXCLUDED audit row per unique excluded FDA
            ingredient within the selected Prescription/OTC marketing scope.

    V7 does not change the v6 chemistry rules; it only makes all prelookup
    exclusions visible in the audit CSV.
    """
    products_path = input_dir / "Products.txt"
    marketing_path = input_dir / "MarketingStatus.txt"

    if not products_path.exists():
        raise FileNotFoundError(f"Missing required file: {products_path}")
    if not marketing_path.exists():
        raise FileNotFoundError(f"Missing required file: {marketing_path}")

    products = load_tsv(products_path)
    marketing_rows = load_tsv(marketing_path)
    keep_statuses = selected_status_ids(args)

    status_by_product: dict[tuple[str, str], str] = {
        (row["ApplNo"], row["ProductNo"]): row["MarketingStatusID"]
        for row in marketing_rows
    }

    seen: set[tuple[str, str]] = set()
    output: list[tuple[str, str]] = []
    prelookup_audit_by_key: dict[str, dict[str, str]] = {}

    skipped_undefined_mixture = 0
    skipped_non_discrete = 0
    skipped_multi_organic_component = 0
    skipped_inorganic_non_ligand = 0
    skipped_diagnostic_imaging = 0
    skipped_excipient_formulation = 0
    skipped_invalid_label = 0
    skipped_therapeutic_peptide = 0
    skipped_radiolabeled_tracer = 0
    skipped_biologic = 0

    def record_prelookup_exclusion(ingredient: str, reason: str) -> None:
        key = key_text(ingredient)
        prelookup_audit_by_key.setdefault(
            key,
            {
                "Disposition": "PRELOOKUP_EXCLUDED",
                "OriginalFDAIngredient": ingredient,
                "CanonicalActiveIngredient": "",
                "InChIKey": "",
                "Details": reason,
            },
        )

    for product in products:
        product_key = (product["ApplNo"], product["ProductNo"])
        status_id = status_by_product.get(product_key)
        if status_id not in keep_statuses:
            continue

        drug_name = clean_text(product.get("DrugName", ""))
        active_ingredient = clean_text(product.get("ActiveIngredient", ""))
        if not drug_name or not active_ingredient:
            continue

        ingredients = (
            [active_ingredient]
            if args.keep_combination_ingredients
            else split_active_ingredients(active_ingredient)
        )
        for ingredient in ingredients:
            # Filter order is fixed so chemistry behaviour stays stable; only the
            # terminal reason is recorded.
            reason = None
            if is_undefined_mixture(ingredient):
                skipped_undefined_mixture += 1
                reason = "undefined mixture"
            elif is_non_discrete_substance(ingredient):
                skipped_non_discrete += 1
                reason = "non-discrete/heterogeneous substance"
            elif is_multi_organic_component_exclusion(ingredient):
                skipped_multi_organic_component += 1
                reason = "defined multi-organic-component ingredient"
            elif is_integral_metal_name_exclusion(ingredient):
                skipped_inorganic_non_ligand += 1
                reason = "integral metal/coordination drug not suitable for conventional docking"
            elif is_inorganic_or_non_ligand(ingredient):
                skipped_inorganic_non_ligand += 1
                reason = "inorganic/mineral/non-ligand"
            elif is_diagnostic_or_imaging(ingredient):
                skipped_diagnostic_imaging += 1
                reason = "diagnostic/imaging agent"
            elif is_excipient_or_formulation(ingredient):
                skipped_excipient_formulation += 1
                reason = "excipient/formulation material"
            elif is_invalid_generic_label(ingredient):
                skipped_invalid_label += 1
                reason = "invalid/generic active-ingredient label"
            # reporting priority: radiopharmaceutical > biologic > peptide.
            elif is_radiolabeled_tracer(ingredient):
                skipped_radiolabeled_tracer += 1
                reason = "radiolabeled diagnostic/radiopharmaceutical"
            elif is_biologic(ingredient):
                skipped_biologic += 1
                reason = "biologic"
            elif is_therapeutic_peptide(ingredient):
                skipped_therapeutic_peptide += 1
                reason = "therapeutic peptide/non-small-molecule"

            if reason is not None:
                record_prelookup_exclusion(ingredient, reason)
                continue

            dedupe_key = (key_text(drug_name), key_text(ingredient))
            if dedupe_key in seen:
                continue
            seen.add(dedupe_key)
            output.append((drug_name, ingredient))

    print(
        "Pre-lookup filter: skipped "
        f"{skipped_undefined_mixture} undefined-mixture, "
        f"{skipped_non_discrete} non-discrete, "
        f"{skipped_multi_organic_component} multi-organic-component, "
        f"{skipped_inorganic_non_ligand} inorganic/non-ligand, "
        f"{skipped_diagnostic_imaging} diagnostic/imaging, "
        f"{skipped_excipient_formulation} excipient/formulation, "
        f"{skipped_invalid_label} invalid-label, "
        f"{skipped_therapeutic_peptide} therapeutic-peptide, "
        f"{skipped_radiolabeled_tracer} radiolabeled-tracer, "
        f"{skipped_biologic} biologic ingredient occurrences "
        "before any structure lookup.",
        file=sys.stderr,
    )

    prelookup_audit = sorted(
        prelookup_audit_by_key.values(),
        key=lambda row: key_text(row["OriginalFDAIngredient"]),
    )
    return (
        sorted(output, key=lambda row: (key_text(row[0]), key_text(row[1]))),
        prelookup_audit,
    )


def load_cache(path: Path) -> dict[str, CachedLookup]:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as handle:
        raw = json.load(handle)
    cache: dict[str, CachedLookup] = {}
    for key, value in raw.items():
        result = value.get("result")
        if value.get("status") == "hit" and isinstance(result, dict):
            add_missing_molecular_weight(result)
        cache[key] = CachedLookup(
            status=value.get("status", "miss"),
            result=result,
            reason=value.get("reason", ""),
        )
    return cache


def save_cache(path: Path, cache: dict[str, CachedLookup]) -> None:
    serializable = {key: asdict(value) for key, value in sorted(cache.items())}
    tmp_path = path.with_suffix(path.suffix + ".tmp")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(serializable, handle, indent=2, sort_keys=True)
        handle.write("\n")
    tmp_path.replace(path)


class JsonClient:
    def __init__(self, delay_seconds: float, timeout_seconds: float) -> None:
        self.delay_seconds = delay_seconds
        self.timeout_seconds = timeout_seconds
        self._last_request_at = 0.0

    def get_json(
        self,
        url: str,
        params: dict[str, str] | None = None,
        retries: int = 3,
    ) -> tuple[dict[str, Any] | None, str]:
        params = params or {}
        full_url = url
        if params:
            full_url = f"{url}?{urllib.parse.urlencode(params)}"

        for attempt in range(retries + 1):
            self._respect_rate_limit()
            request = urllib.request.Request(
                full_url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": "redrug-fda-structures/0.1",
                },
            )
            try:
                with urllib.request.urlopen(request, timeout=self.timeout_seconds) as response:
                    body = response.read().decode("utf-8")
                return json.loads(body), full_url
            except urllib.error.HTTPError as exc:
                if exc.code == 404:
                    return None, full_url
                if exc.code in {429, 500, 502, 503, 504} and attempt < retries:
                    time.sleep(min(60.0, 2.0 ** attempt + self.delay_seconds))
                    continue
                return None, full_url
            except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
                if attempt < retries:
                    time.sleep(min(60.0, 2.0 ** attempt + self.delay_seconds))
                    continue
                return None, full_url

        return None, full_url

    def _respect_rate_limit(self) -> None:
        elapsed = time.monotonic() - self._last_request_at
        if elapsed < self.delay_seconds:
            time.sleep(self.delay_seconds - elapsed)
        self._last_request_at = time.monotonic()


def get_nested(record: dict[str, Any], *names: str) -> Any:
    for name in names:
        if name in record:
            return record[name]
    return None


def normalize_substance_class(record: dict[str, Any]) -> str:
    value = get_nested(record, "substanceClass", "substance_class")
    return clean_text(str(value or "")).lower()


def extract_names(record: dict[str, Any]) -> list[str]:
    names: list[str] = []
    preferred = get_nested(record, "_name", "preferred_name", "approvalIDDisplay")
    if preferred:
        names.append(str(preferred))
    for name_record in record.get("names", []) or []:
        if isinstance(name_record, dict) and name_record.get("name"):
            names.append(str(name_record["name"]))
        elif isinstance(name_record, str):
            names.append(name_record)
    return [clean_text(name) for name in names if clean_text(name)]


def exact_weight_from_smiles(smiles: str) -> float | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return float(Descriptors.ExactMolWt(mol))


def add_missing_molecular_weight(result: dict[str, Any]) -> None:
    if result.get("molecular_weight") is not None:
        return
    smiles = result.get("isomeric_smiles")
    if not smiles:
        return
    molecular_weight = exact_weight_from_smiles(str(smiles))
    if molecular_weight is not None:
        result["molecular_weight"] = molecular_weight



def structure_has_integral_metal_issue(smiles: str) -> bool:
    """Detect metal chemistry that largest-fragment stripping would corrupt.

    Safe ordinary salts usually have one carbon-containing organic fragment and
    a disconnected metal counterion.  Reject when a metal is in the SAME
    carbon-containing fragment, or when a metal coexists with multiple organic
    fragments (a common representation for coordination/metal-ligand complexes).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False
    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    organic_fragments = 0
    any_metal = False
    for frag in fragments:
        atomic_numbers = {atom.GetAtomicNum() for atom in frag.GetAtoms()}
        has_carbon = 6 in atomic_numbers
        has_metal = bool(atomic_numbers & METAL_ATOMIC_NUMBERS)
        any_metal = any_metal or has_metal
        if has_carbon:
            organic_fragments += 1
        if has_carbon and has_metal:
            return True
    if any_metal and organic_fragments > 1:
        return True
    return False


def parent_identity_structure(smiles: str) -> tuple[str, str] | None:
    """Return a charge-normalized parent SMILES/InChIKey for identity dedup.

    This does NOT attempt physiological protonation.  It only makes common acid/
    base salt representations comparable (e.g. ascorbate vs ascorbic acid).
    Docking protonation/state preparation belongs downstream.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        mol = rdMolStandardize.Cleanup(mol)
        mol = rdMolStandardize.Uncharger().uncharge(mol)
        Chem.SanitizeMol(mol)
    except Exception:
        return None
    parent_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    parent_inchikey = rd_inchi.MolToInchiKey(mol)
    if not parent_smiles or not parent_inchikey:
        return None
    return parent_smiles, parent_inchikey


def standardize_structure(smiles: str, keep_whole: bool = False) -> tuple[str, str, float] | None:
    """
    Standardize a raw (possibly multi-fragment/salt) SMILES to a single
    reportable structure.

    keep_whole=False (default): keep only the LARGEST organic (carbon-
    containing) fragment, discarding disconnected counter-ions. Correct
    for true salts of an organic drug (e.g. "ERLOTINIB HYDROCHLORIDE" ->
    erlotinib; "AMPHETAMINE SULFATE" -> amphetamine), where the counter-
    ion is pharmacologically interchangeable.

    keep_whole=True: keep the FULL multi-fragment structure (all ions
    together) instead of discarding anything. Required for simple
    mineral salts where a carbon-containing anion (acetate, carbonate,
    gluconate, citrate, lactate) is paired with DIFFERENT metal cations
    that are each their OWN clinically distinct drug -- e.g. "CALCIUM
    CARBONATE" (antacid), "LANTHANUM CARBONATE" (dialysis phosphate
    binder), and "LITHIUM CARBONATE" (bipolar disorder mood stabilizer)
    would otherwise all collapse to bare "carbonate" under the default
    rule, since carbonate (CO3^2-) contains carbon and gets kept while
    the metal -- which is the part that actually matters here -- gets
    discarded as if it were just an inert counter-ion. See
    METAL_SALT_IDENTITY_CRITICAL below for the full list and rationale;
    this was found via direct verification against real output, not
    assumed in advance.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    # defense-in-depth: never silently turn an integral-metal drug into an
    # unrelated organic ligand by deleting the metal.
    if structure_has_integral_metal_issue(smiles):
        return None

    if keep_whole:
        isomeric_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        inchikey = rd_inchi.MolToInchiKey(mol)
        if not isomeric_smiles or not inchikey:
            return None
        molecular_weight = float(Descriptors.ExactMolWt(mol))
        return isomeric_smiles, inchikey, molecular_weight

    fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    candidates = []
    for fragment in fragments:
        atoms = list(fragment.GetAtoms())
        has_carbon = any(atom.GetAtomicNum() == 6 for atom in atoms)
        if not has_carbon:
            continue
        heavy_atom_count = sum(1 for atom in atoms if atom.GetAtomicNum() > 1)
        exact_weight = Descriptors.ExactMolWt(fragment)
        candidates.append((heavy_atom_count, exact_weight, fragment))

    if not candidates:
        return None

    fragment = max(candidates, key=lambda item: (item[0], item[1]))[2]
    isomeric_smiles = Chem.MolToSmiles(fragment, canonical=True, isomericSmiles=True)
    inchikey = rd_inchi.MolToInchiKey(fragment)
    if not isomeric_smiles or not inchikey:
        return None
    molecular_weight = float(Descriptors.ExactMolWt(fragment))
    return isomeric_smiles, inchikey, molecular_weight


# Backwards-compatible alias -- other call sites in this file use this name.
def standardize_to_largest_organic_fragment(smiles: str) -> tuple[str, str, float] | None:
    return standardize_structure(smiles, keep_whole=False)


# ==============================================================================
# METAL SALT IDENTITY OVERRIDE
# ==============================================================================
# Mineral salts where a carbon-containing anion would be kept as "the
# organic fragment" while the metal -- the part that defines the drug -- is
# discarded. Calcium carbonate (antacid), lanthanum carbonate (phosphate
# binder) and lithium carbonate (mood stabilizer) all collapsed to bare
# "carbonate" before this override existed.
#
# Deliberately not listed: salts where the anion is the active species and
# the metal is interchangeable, for which the default largest-fragment rule
# is already correct.
#   - AMPHETAMINE SULFATE, AMPICILLIN TRIHYDRATE: organic drug plus a simple
#     counter-ion.
#   - CALCIUM/MAGNESIUM/POTASSIUM/SODIUM OXYBATE: the oxybate anion is
#     active and the mixed metals exist to reduce sodium content, so
#     merging them is correct.
#
# This list covers the anion families seen so far (acetate, carbonate,
# gluconate, citrate, lactate) and is not exhaustive. If another metal salt
# of one of these families collapses unexpectedly, add it here.
METAL_SALT_IDENTITY_CRITICAL = {
    # Acetate salts -- different metals are different drugs
    "CALCIUM ACETATE", "MAGNESIUM ACETATE", "POTASSIUM ACETATE",
    "SODIUM ACETATE", "ZINC ACETATE",
    # Carbonate salts
    "CALCIUM CARBONATE", "LANTHANUM CARBONATE", "LITHIUM CARBONATE",
    "MAGNESIUM CARBONATE", "SODIUM CARBONATE",
    # Gluconate salts
    "CALCIUM GLUCONATE", "SODIUM GLUCONATE", "FERROUS GLUCONATE",
    "ZINC GLUCONATE", "POTASSIUM GLUCONATE", "MAGNESIUM GLUCONATE",
    # Citrate salts
    "CALCIUM CITRATE", "MAGNESIUM CITRATE", "POTASSIUM CITRATE",
    "SODIUM CITRATE", "FERRIC CITRATE", "LITHIUM CITRATE",
    # Lactate salts
    "AMMONIUM LACTATE", "SODIUM LACTATE", "CALCIUM LACTATE",
    "POTASSIUM LACTATE",
}


def metal_salt_identity_is_critical(ingredient: str) -> bool:
    return key_text(ingredient) in METAL_SALT_IDENTITY_CRITICAL


def openfda_query_variants(ingredient: str) -> list[str]:
    variants = [ingredient]
    without_parentheses = re.sub(r"\([^)]*\)", " ", ingredient)
    without_parentheses = clean_text(without_parentheses)
    if without_parentheses and key_text(without_parentheses) != key_text(ingredient):
        variants.append(without_parentheses)
    return unique_preserve_order(variants)


# Counter-ion words stripped to build an alternate query if the exact name
# is not indexed. Includes "RESIN COMPLEX", an ion-exchange formulation
# suffix equivalent to POLISTIREX: same molecule, different formulation.
#
# Deliberately excluded: ACETATE, PROPIONATE, PHOSPHATE, SUCCINATE,
# ENANTHATE, CYPIONATE, VALERATE, CAPROATE, PIVALATE, DECANOATE, BUTYRATE,
# FUROATE, ACETONIDE, NITRATE. Each is sometimes a plain salt and sometimes
# the drug's own distinct identity, so they are governed by the explicit
# maps above rather than a blind strip-and-retry that could silently query
# the wrong structure for something like MEDROXYPROGESTERONE ACETATE.
_PUBCHEM_SALT_WORDS = (
    "MONOHYDROCHLORIDE", "DIHYDROCHLORIDE", "HYDROCHLORIDE", "HYDROBROMIDE",
    "HEMIHYDRATE", "MONOHYDRATE", "DIHYDRATE", "TRIHYDRATE", "ANHYDROUS",
    "SODIUM", "POTASSIUM", "CALCIUM", "MAGNESIUM", "ZINC", "LITHIUM",
    "BESYLATE", "BESILATE", "MALEATE", "MALATE", "FUMARATE", "TARTRATE",
    "CITRATE", "MESYLATE", "MESILATE", "TOSYLATE", "TOSILATE",
    "SULFATE", "SULPHATE", "BROMIDE", "CHLORIDE", "IODIDE", "HYCLATE",
    "GLUCONATE", "LACTATE", "SALICYLATE", "CAMSYLATE", "CAMSILATE",
    "EDISYLATE", "ISETHIONATE", "NAPSYLATE", "PAMOATE",
    "POLYGALACTURONATE", "GLUCEPTATE", "HIPPURATE",
    "MEGLUMINE", "TROMETHAMINE", "DIOLAMINE", "OLAMINE", "ARGININE", "LYSINE",
    "TANNATE", "POLISTIREX", "XINAFOATE", "DISODIUM", "BISULFATE",
    "MONOSULFATE", "MEPESUCCINATE", "ASPARTATE", "SACCHARATE",
)
_PUBCHEM_SALT_PATTERN = re.compile(
    r"\b(" + "|".join(_PUBCHEM_SALT_WORDS) + r")\b", flags=re.IGNORECASE
)
_RESIN_COMPLEX_PATTERN = re.compile(r"\s+RESIN COMPLEX\s*$", flags=re.IGNORECASE)


def pubchem_query_variants(ingredient: str) -> list[str]:
    variants = openfda_query_variants(ingredient)

    no_resin_complex = clean_text(_RESIN_COMPLEX_PATTERN.sub("", ingredient))
    if no_resin_complex and key_text(no_resin_complex) != key_text(ingredient):
        variants.append(no_resin_complex)

    desalted = clean_text(_PUBCHEM_SALT_PATTERN.sub(" ", ingredient))
    if desalted and key_text(desalted) != key_text(ingredient):
        variants.append(desalted)

    return unique_preserve_order(variants)


def unique_preserve_order(values: list[str]) -> list[str]:
    seen: set[str] = set()
    output: list[str] = []
    for value in values:
        key = key_text(value)
        if not value or key in seen:
            continue
        seen.add(key)
        output.append(value)
    return output


def lookup_openfda(
    ingredient: str,
    client: JsonClient,
    api_key: str | None,
    keep_whole: bool = False,
) -> StructureResult | None:
    for query in openfda_query_variants(ingredient):
        params = {"search": f'names.name:"{query}"', "limit": "10"}
        if api_key:
            params["api_key"] = api_key

        payload, _ = client.get_json(OPENFDA_SUBSTANCE_URL, params)
        if not payload or not payload.get("results"):
            continue

        for record in sorted(payload["results"], key=lambda item: score_openfda_result(item, query), reverse=True):
            # Safeguard: openFDA search can return related names.
            # Never accept a record unless the queried canonical name is an
            # exact registered/synonym name on that record.
            record_name_keys = {key_text(n) for n in extract_names(record)}
            if key_text(query) not in record_name_keys:
                continue
            substance_class = normalize_substance_class(record)
            if substance_class == "polymer" or record.get("polymer"):
                return None

            structure = record.get("structure") or {}
            if not isinstance(structure, dict):
                continue

            raw_smiles = get_nested(structure, "smiles", "SMILES")
            if not raw_smiles:
                continue

            standardized = standardize_structure(str(raw_smiles), keep_whole=keep_whole)
            if standardized is None:
                continue

            uuid = clean_text(str(record.get("uuid", "")))
            if uuid:
                source = f"pFDA/GSRS: https://precision.fda.gov/ginas/app/ui/substances/{uuid}"
            else:
                source = f"FDA/openFDA GSRS: {OPENFDA_SUBSTANCE_URL}"

            source_names = extract_names(record)
            source_name = source_names[0] if source_names else query
            return StructureResult(
                isomeric_smiles=standardized[0],
                inchikey=standardized[1],
                molecular_weight=standardized[2],
                source=source,
                source_name=source_name,
                source_substance_class=substance_class,
            )
    return None


def score_openfda_result(record: dict[str, Any], query: str) -> tuple[int, int, int]:
    names = [key_text(name) for name in extract_names(record)]
    query_key = key_text(query)
    substance_class = normalize_substance_class(record)
    structure = record.get("structure") or {}
    has_structure = isinstance(structure, dict) and bool(get_nested(structure, "smiles", "SMILES"))
    return (
        1 if query_key in names else 0,
        1 if substance_class == "chemical" else 0,
        1 if has_structure else 0,
    )



def lookup_pubchem_cid(
    cid: int,
    canonical_name: str,
    client: JsonClient,
    keep_whole: bool = False,
) -> StructureResult | None:
    """Deterministic PubChem lookup for identities where name search is unsafe."""
    for property_name in ("IsomericSMILES", "SMILES"):
        url = (
            f"{PUBCHEM_PUG_URL}/compound/cid/{cid}/property/"
            f"{property_name},InChIKey/JSON"
        )
        payload, full_url = client.get_json(url)
        if not payload:
            continue
        properties = payload.get("PropertyTable", {}).get("Properties", [])
        if not properties:
            continue
        row = properties[0]
        raw_smiles = row.get(property_name) or row.get("SMILES") or row.get("IsomericSMILES")
        if not raw_smiles:
            continue
        standardized = standardize_structure(str(raw_smiles), keep_whole=keep_whole)
        if standardized is None:
            continue
        return StructureResult(
            isomeric_smiles=standardized[0],
            inchikey=standardized[1],
            molecular_weight=standardized[2],
            source=f"PubChem CID {cid}: https://pubchem.ncbi.nlm.nih.gov/compound/{cid}",
            source_name=canonical_name,
        )
    return None


def lookup_pubchem(ingredient: str, client: JsonClient, keep_whole: bool = False) -> StructureResult | None:
    for query in pubchem_query_variants(ingredient):
        encoded_name = urllib.parse.quote(query, safe="")
        for property_name in ("IsomericSMILES", "SMILES"):
            url = (
                f"{PUBCHEM_PUG_URL}/compound/name/{encoded_name}/property/"
                f"{property_name},InChIKey/JSON"
            )
            payload, full_url = client.get_json(url)
            if not payload:
                continue

            properties = (
                payload.get("PropertyTable", {}).get("Properties", [])
                if isinstance(payload, dict)
                else []
            )
            if not properties:
                continue

            row = properties[0]
            raw_smiles = row.get(property_name) or row.get("SMILES") or row.get("IsomericSMILES")
            if not raw_smiles:
                continue

            standardized = standardize_structure(str(raw_smiles), keep_whole=keep_whole)
            if standardized is None:
                continue

            cid = row.get("CID")
            source = (
                f"PubChem: https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
                if cid
                else f"PubChem PUG REST: {full_url}"
            )
            return StructureResult(
                isomeric_smiles=standardized[0],
                inchikey=standardized[1],
                molecular_weight=standardized[2],
                source=source,
                source_name=query,
            )
    return None


def lookup_structure(
    ingredient: str,
    client: JsonClient,
    api_key: str | None,
    args: argparse.Namespace,
) -> CachedLookup:
    if looks_polymer_like(ingredient):
        return CachedLookup(status="skip", reason="polymer-like ingredient name")

    # run project-scope filters on the RAW identity, then canonicalize and
    # run the SAME filters again.  Canonicalization is not allowed to create an
    # invalid or excluded active-moiety identity.
    raw_reason = trem2_exclusion_reason(ingredient)
    if raw_reason:
        return CachedLookup(status="skip", reason=f"pre-canonicalization exclusion: {raw_reason}")

    canonical_name = canonical_active_moiety_name(ingredient)
    canonical_reason = trem2_exclusion_reason(canonical_name)
    if canonical_reason:
        return CachedLookup(
            status="skip",
            reason=f"post-canonicalization exclusion ({canonical_name}): {canonical_reason}",
        )

    query_ingredient = canonical_lookup_name(ingredient)

    # This script no longer preserves whole mineral salts.  Mineral products are
    # excluded by project-scope rules; ordinary organic drug salts are reduced
    # to their parent organic fragment.
    keep_whole = False

    result = None
    canonical_key = key_text(canonical_name)
    cid_override = PUBCHEM_CID_OVERRIDES.get(canonical_key)

    if cid_override is not None:
        if args.no_pubchem:
            return CachedLookup(
                status="skip",
                reason=(
                    f"stereochemistry-sensitive identity {canonical_name} requires "
                    f"deterministic PubChem CID {cid_override}; --no-pubchem was supplied"
                ),
            )
        result = lookup_pubchem_cid(cid_override, canonical_name, client, keep_whole=False)
        # Do NOT silently fall back to an ambiguous name search for CID-pinned
        # stereochemistry-sensitive drugs.
        if result is None:
            return CachedLookup(
                status="miss",
                reason=f"deterministic PubChem CID {cid_override} lookup failed for {canonical_name}",
            )
    else:
        if not args.no_openfda:
            result = lookup_openfda(query_ingredient, client, api_key, keep_whole=keep_whole)
        if result is None and not args.no_pubchem:
            result = lookup_pubchem(query_ingredient, client, keep_whole=keep_whole)

    if result is None:
        return CachedLookup(status="miss", reason="no non-polymer organic structure found")
    if result.source_substance_class.lower() == "polymer":
        return CachedLookup(status="skip", reason="source substance class is polymer")

    # Normalize charge state only for identity/dedup representation.  This is a
    # dataset identity step, NOT a physiological protonation prediction.
    parent = parent_identity_structure(result.isomeric_smiles)
    if parent is None:
        return CachedLookup(status="skip", reason="could not build normalized parent identity")
    parent_smiles, parent_inchikey = parent
    parent_mol = Chem.MolFromSmiles(parent_smiles)
    if parent_mol is None:
        return CachedLookup(status="skip", reason="normalized parent SMILES is invalid")
    if any(atom.GetAtomicNum() in METAL_ATOMIC_NUMBERS for atom in parent_mol.GetAtoms()):
        return CachedLookup(status="skip", reason="normalized parent still contains an integral metal")
    if "." in parent_smiles:
        return CachedLookup(status="skip", reason="normalized parent remains multicomponent")

    parent_mw = float(Descriptors.ExactMolWt(parent_mol))
    if parent_mw > args.max_molecular_weight:
        return CachedLookup(
            status="skip",
            result=asdict(result),
            reason=(
                f"normalized parent molecular weight {parent_mw:.2f} Da exceeds "
                f"{args.max_molecular_weight:.2f} Da"
            ),
        )

    result.isomeric_smiles = parent_smiles
    result.inchikey = parent_inchikey
    result.molecular_weight = parent_mw
    return CachedLookup(status="hit", result=asdict(result))


def passes_output_filters(cached: CachedLookup, max_molecular_weight: float) -> bool:
    if cached.status != "hit" or not cached.result:
        return False
    add_missing_molecular_weight(cached.result)
    molecular_weight = cached.result.get("molecular_weight")
    if molecular_weight is None:
        return False
    return float(molecular_weight) <= max_molecular_weight


def _audit_disposition_for_cached(cached: CachedLookup | None) -> tuple[str, str]:
    """Map a cache outcome to one terminal v7 audit disposition."""
    if cached is None:
        return "LOOKUP_MISS", "no cache record was produced for this ingredient"
    if cached.status == "miss":
        return "LOOKUP_MISS", cached.reason or "structure lookup failed"
    if cached.status == "skip":
        reason = cached.reason or "structure/QC exclusion"
        lower = reason.lower()
        if "molecular weight" in lower or "exceeds" in lower:
            return "MW_EXCLUDED", reason
        if "polymer" in lower:
            return "POLYMER_EXCLUDED", reason
        if "integral metal" in lower or "coordination" in lower:
            return "INTEGRAL_METAL_EXCLUDED", reason
        if "pre-canonicalization exclusion" in lower or "post-canonicalization exclusion" in lower:
            return "PRELOOKUP_EXCLUDED", reason
        return "STRUCTURE_QC_EXCLUDED", reason
    if cached.status != "hit":
        return "STRUCTURE_QC_EXCLUDED", cached.reason or f"unexpected cache status: {cached.status}"
    return "", ""



def selected_cache_statistics(
    ingredients: list[str],
    cache: dict[str, CachedLookup],
    max_molecular_weight: float,
) -> tuple[int, int, int, int]:
    """Return hit/miss/skip/usable counts for ONLY the requested FDA ingredients.

    Canonical parent aliases may also exist as extra keys in the cache.  Those
    implementation-detail keys must never inflate run statistics.
    """
    selected_values = [
        cache.get(key_text(ingredient))
        for ingredient in ingredients
    ]
    hit_count = sum(
        1 for value in selected_values
        if value is not None and value.status == "hit"
    )
    miss_count = sum(
        1 for value in selected_values
        if value is not None and value.status == "miss"
    )
    skip_count = sum(
        1 for value in selected_values
        if value is not None and value.status == "skip"
    )
    usable_hit_count = sum(
        1 for value in selected_values
        if value is not None and passes_output_filters(value, max_molecular_weight)
    )
    return hit_count, miss_count, skip_count, usable_hit_count


def write_csv(
    output_path: Path,
    product_rows: list[tuple[str, str]],
    cache: dict[str, CachedLookup],
    max_molecular_weight: float,
    review_output_path: Path | None = None,
    prelookup_audit: list[dict[str, str]] | None = None,
) -> tuple[int, int, dict[str, int]]:
    """Write the v7 docking library and complete per-ingredient audit.

    The v6 chemistry/QC behavior is unchanged.  V7 adds:
      * corrected output column names: "isomeric SMILES" and "InChIKey"
      * exactly one terminal disposition row for every unique FDA ingredient
        that either entered the selected Prescription/OTC scope and was
        prefiltered, or survived prefiltering and proceeded to lookup.

    Rows involved in unresolved identity/structure collisions remain withheld
    from the final docking CSV.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    compounds: dict[tuple[str, str], dict[str, Any]] = {}
    identity_to_names: dict[str, set[str]] = {}
    canonical_to_identity_keys: dict[str, set[str]] = {}
    ingredient_to_group: dict[str, tuple[str, str]] = {}
    ingredient_canonical: dict[str, str] = {}
    ingredient_result_inchikey: dict[str, str] = {}
    qc_exclusions_by_ingredient: dict[str, dict[str, str]] = {}

    unique_ingredients = sorted({ingredient for _, ingredient in product_rows}, key=key_text)

    # First, classify cache-level failures/skips and build candidate groups for hits.
    for ingredient in unique_ingredients:
        cache_key = key_text(ingredient)
        cached = cache.get(cache_key)
        disposition, details = _audit_disposition_for_cached(cached)
        canonical_name = canonical_active_moiety_name(ingredient)
        ingredient_canonical[cache_key] = canonical_name

        if disposition:
            qc_exclusions_by_ingredient[cache_key] = {
                "Disposition": disposition,
                "OriginalFDAIngredient": ingredient,
                "CanonicalActiveIngredient": canonical_name if canonical_name else "",
                "InChIKey": (
                    str((cached.result or {}).get("inchikey", ""))
                    if cached is not None and cached.result
                    else ""
                ),
                "Details": details,
            }
            continue

        assert cached is not None and cached.result is not None

        # Preserve the post-canonicalization QC gate.
        canonical_reason = trem2_exclusion_reason(canonical_name)
        if canonical_reason:
            qc_exclusions_by_ingredient[cache_key] = {
                "Disposition": "PRELOOKUP_EXCLUDED",
                "OriginalFDAIngredient": ingredient,
                "CanonicalActiveIngredient": canonical_name,
                "InChIKey": "",
                "Details": f"post-canonicalization exclusion: {canonical_reason}",
            }
            continue

        result = cached.result or {}
        smiles = str(result.get("isomeric_smiles", ""))
        parent = parent_identity_structure(smiles)
        if parent is None:
            qc_exclusions_by_ingredient[cache_key] = {
                "Disposition": "STRUCTURE_QC_EXCLUDED",
                "OriginalFDAIngredient": ingredient,
                "CanonicalActiveIngredient": canonical_name,
                "InChIKey": str(result.get("inchikey", "")),
                "Details": "could not build normalized parent identity",
            }
            continue

        parent_smiles, parent_inchikey = parent
        ingredient_result_inchikey[cache_key] = parent_inchikey
        mol = Chem.MolFromSmiles(parent_smiles)
        if mol is None or "." in parent_smiles:
            qc_exclusions_by_ingredient[cache_key] = {
                "Disposition": "STRUCTURE_QC_EXCLUDED",
                "OriginalFDAIngredient": ingredient,
                "CanonicalActiveIngredient": canonical_name,
                "InChIKey": parent_inchikey,
                "Details": f"invalid or multicomponent normalized output structure: {parent_smiles}",
            }
            continue
        if any(atom.GetAtomicNum() in METAL_ATOMIC_NUMBERS for atom in mol.GetAtoms()):
            qc_exclusions_by_ingredient[cache_key] = {
                "Disposition": "INTEGRAL_METAL_EXCLUDED",
                "OriginalFDAIngredient": ingredient,
                "CanonicalActiveIngredient": canonical_name,
                "InChIKey": parent_inchikey,
                "Details": f"normalized parent still contains integral metal: {parent_smiles}",
            }
            continue

        canonical_key = key_text(canonical_name)
        group_key = (canonical_key, parent_inchikey)
        ingredient_to_group[cache_key] = group_key
        compound = compounds.setdefault(
            group_key,
            {
                "canonical_name": canonical_name,
                "aliases": set(),
                "isomeric_smiles": parent_smiles,
                "sources": set(),
            },
        )
        compound["aliases"].add(ingredient)
        compound["sources"].add(str(result.get("source", "")))
        identity_to_names.setdefault(parent_inchikey, set()).add(canonical_name)
        canonical_to_identity_keys.setdefault(canonical_key, set()).add(parent_inchikey)

    collision_identity_keys = {
        k for k, vals in identity_to_names.items() if len({key_text(v) for v in vals}) > 1
    }
    collision_canonical_keys = {
        k for k, vals in canonical_to_identity_keys.items() if len(vals) > 1
    }

    # Final docking CSV: only collision-free groups.
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["ActiveIngredient", "isomeric SMILES", "InChIKey", "Source"],
        )
        writer.writeheader()
        written = 0
        for (canonical_key, inchikey), compound in sorted(
            compounds.items(),
            key=lambda item: (key_text(item[1]["canonical_name"]), item[0][1]),
        ):
            if inchikey in collision_identity_keys or canonical_key in collision_canonical_keys:
                continue
            writer.writerow(
                {
                    "ActiveIngredient": compound["canonical_name"],
                    "isomeric SMILES": compound["isomeric_smiles"],
                    "InChIKey": inchikey,
                    "Source": "; ".join(sorted(x for x in compound["sources"] if x)),
                }
            )
            written += 1

    # Build one terminal disposition per unique FDA ingredient.
    audit_by_key: dict[str, dict[str, str]] = {}
    for row in prelookup_audit or []:
        audit_by_key[key_text(row["OriginalFDAIngredient"])] = dict(row)

    for ingredient in unique_ingredients:
        cache_key = key_text(ingredient)
        if cache_key in qc_exclusions_by_ingredient:
            audit_by_key[cache_key] = qc_exclusions_by_ingredient[cache_key]
            continue

        group_key = ingredient_to_group.get(cache_key)
        canonical_name = ingredient_canonical.get(cache_key, canonical_active_moiety_name(ingredient))
        parent_inchikey = ingredient_result_inchikey.get(cache_key, "")
        if group_key is None:
            # Defensive catch: a hit should have been assigned a group or an
            # explicit QC disposition above.
            audit_by_key[cache_key] = {
                "Disposition": "STRUCTURE_QC_EXCLUDED",
                "OriginalFDAIngredient": ingredient,
                "CanonicalActiveIngredient": canonical_name,
                "InChIKey": parent_inchikey,
                "Details": "internal audit fallback: hit was not assigned to an output group",
            }
            continue

        canonical_key, inchikey = group_key
        if inchikey in collision_identity_keys or canonical_key in collision_canonical_keys:
            reasons = []
            if inchikey in collision_identity_keys:
                reasons.append(
                    "same normalized parent InChIKey maps to multiple canonical drug names"
                )
            if canonical_key in collision_canonical_keys:
                reasons.append(
                    "same canonical drug name maps to multiple normalized parent InChIKeys"
                )
            audit_by_key[cache_key] = {
                "Disposition": "HARD_QC_WITHHELD",
                "OriginalFDAIngredient": ingredient,
                "CanonicalActiveIngredient": canonical_name,
                "InChIKey": inchikey,
                "Details": "WITHHELD FROM FINAL CSV: " + "; ".join(reasons),
            }
            continue

        if key_text(ingredient) == key_text(canonical_name):
            disposition = "KEPT"
            details = "FDA ingredient retained as the normalized docking active moiety"
        else:
            disposition = "CANONICALIZED"
            details = "FDA salt/prodrug/qualifier/form name mapped to normalized docking active moiety"

        audit_by_key[cache_key] = {
            "Disposition": disposition,
            "OriginalFDAIngredient": ingredient,
            "CanonicalActiveIngredient": canonical_name,
            "InChIKey": inchikey,
            "Details": details,
        }

    if review_output_path is None:
        review_output_path = output_path.with_name(output_path.stem + "_review.csv")
    review_output_path.parent.mkdir(parents=True, exist_ok=True)

    audit_fieldnames = [
        "Disposition",
        "OriginalFDAIngredient",
        "CanonicalActiveIngredient",
        "InChIKey",
        "Details",
    ]
    with review_output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_fieldnames)
        writer.writeheader()
        for row in sorted(
            audit_by_key.values(),
            key=lambda r: (key_text(r["OriginalFDAIngredient"]), r["Disposition"]),
        ):
            writer.writerow(row)

    disposition_counts: dict[str, int] = {}
    for row in audit_by_key.values():
        disposition = row["Disposition"]
        disposition_counts[disposition] = disposition_counts.get(disposition, 0) + 1

    hard_qc_issues = disposition_counts.get("HARD_QC_WITHHELD", 0)
    return written, hard_qc_issues, disposition_counts


def run_self_test() -> int:
    """Offline sanity checks for rules; makes no network requests."""
    canonical_cases = {
        "AMLODIPINE BESYLATE": "AMLODIPINE",
        "AMPICILLIN TRIHYDRATE": "AMPICILLIN",
        "BETAMETHASONE SODIUM PHOSPHATE": "BETAMETHASONE",
        "BUPROPION HYDROBROMIDE": "BUPROPION",
        "CALCIUM OXYBATE": "OXYBATE",
        "CHLORPHENIRAMINE POLISTIREX": "CHLORPHENIRAMINE",
        "DEXMETHYLPHENIDATE HYDROCHLORIDE": "DEXMETHYLPHENIDATE",
        "METHYLPHENIDATE HYDROCHLORIDE": "METHYLPHENIDATE",
        "DIVALPROEX SODIUM": "VALPROIC ACID",
        "FORMOTEROL FUMARATE DIHYDRATE": "FORMOTEROL",
        "HYDROXYZINE PAMOATE": "HYDROXYZINE",
        "IMIPRAMINE PAMOATE": "IMIPRAMINE",
        "OLANZAPINE PAMOATE": "OLANZAPINE",
        "NITROFURANTOIN, MACROCRYSTALLINE": "NITROFURANTOIN",
        "PENICILLIN G PROCAINE": "PENICILLIN G",
        "PITAVASTATIN CALCIUM": "PITAVASTATIN",
        "SITAGLIPTIN PHOSPHATE": "SITAGLIPTIN",
        "TESTOSTERONE CYPIONATE": "TESTOSTERONE",
        # output-audit corrections
        "ETHACRYNATE": "ETHACRYNIC ACID",
        "DICLOFENAC EPOLAMINE": "DICLOFENAC",
        "DIAZOXIDE CHOLINE": "DIAZOXIDE",
        "ELTROMBOPAG CHOLINE": "ELTROMBOPAG",
        "PERINDOPRIL ERBUMINE": "PERINDOPRIL",
        "CHOLINE FENOFIBRATE": "FENOFIBRIC ACID",
        "SODIUM ASCORBATE": "ASCORBIC ACID",
        # Names that MUST remain real drug identities, not be suffix-truncated.
        "METHYL SALICYLATE": "METHYL SALICYLATE",
        "DIMETHYL FUMARATE": "DIMETHYL FUMARATE",
        "MONOMETHYL FUMARATE": "MONOMETHYL FUMARATE",
        # targeted 100-drug audit corrections
        "ADEFOVIR DIPIVOXIL": "ADEFOVIR",
        "ARIPIPRAZOLE LAUROXIL": "ARIPIPRAZOLE",
        "AMMONIUM LACTATE": "LACTIC ACID",
        "ALISKIREN HEMIFUMARATE": "ALISKIREN",
    }
    failures: list[str] = []
    for original, expected in canonical_cases.items():
        observed = canonical_active_moiety_name(original)
        if key_text(observed) != key_text(expected):
            failures.append(f"canonical: {original!r} -> {observed!r}; expected {expected!r}")

    exclusion_cases = {
        "Carbon Dioxide, USP",
        "SODIUM FERRIC GLUCONATE COMPLEX",
        "CARBOPLATIN",
        "OXALIPLATIN",
        "AURANOFIN",
        "FERRIC MALTOL",
        "SILVER SULFADIAZINE",
        "CALCIUM CARBONATE",
        "LITHIUM CARBONATE",
        "SODIUM BICARBONATE",
        "AMMONIUM CHLORIDE",
        "ANGIOTENSIN II",
        "ANGIOTENSIN II ACETATE",
    }
    for name in exclusion_cases:
        if trem2_exclusion_reason(name) is None:
            failures.append(f"filter: {name!r} was not excluded")

    if not is_multi_organic_component_exclusion("DIMENHYDRINATE"):
        failures.append("filter: DIMENHYDRINATE was not excluded")

    # exact classification regression checks.
    if trem2_exclusion_reason("AMMONIUM CHLORIDE") != "inorganic/mineral/non-ligand":
        failures.append("classification: AMMONIUM CHLORIDE was not classified as inorganic/non-ligand")
    if trem2_exclusion_reason("ANGIOTENSIN II ACETATE") != "therapeutic peptide/non-small-molecule":
        failures.append("classification: ANGIOTENSIN II ACETATE was not classified as a therapeutic peptide")

    # Structural metal defense: disconnected sodium salt is okay; integral Pt is not.
    if structure_has_integral_metal_issue("CC(=O)[O-].[Na+]"):
        failures.append("metal: ordinary disconnected sodium salt was falsely rejected")
    if not structure_has_integral_metal_issue("C[Pt]C"):
        failures.append("metal: integral organoplatinum structure was not detected")

    # Charge-normalized identity should merge an acid with its conjugate-base form.
    acid = parent_identity_structure("CC(=O)O")
    acetate = parent_identity_structure("CC(=O)[O-]")
    if not acid or not acetate or acid[1] != acetate[1]:
        failures.append("identity: acetic acid/acetate did not normalize to the same parent InChIKey")

    # audit/reporting regression checks.
    if is_therapeutic_peptide("CANGRELOR"):
        failures.append("classification: CANGRELOR was still mislabeled as a peptide")
    if trem2_exclusion_reason("GALLIUM GA-68 EDOTREOTIDE") != "radiolabeled diagnostic/radiopharmaceutical":
        failures.append("classification: radiolabeled peptide tracer did not prioritize radiopharmaceutical")
    if trem2_exclusion_reason("CAPROMAB PENDETIDE") != "biologic":
        failures.append("classification: CAPROMAB PENDETIDE did not prioritize biologic over peptide")

    # CID-pinned sensitive identities must remain distinct.
    if PUBCHEM_CID_OVERRIDES["ARFORMOTEROL"] == PUBCHEM_CID_OVERRIDES["FORMOTEROL"]:
        failures.append("CID override: ARFORMOTEROL and FORMOTEROL are not distinct")
    if PUBCHEM_CID_OVERRIDES["LEVONORGESTREL"] == PUBCHEM_CID_OVERRIDES["NORGESTREL"]:
        failures.append("CID override: LEVONORGESTREL and NORGESTREL are not distinct")
    if PUBCHEM_CID_OVERRIDES["DEXMETHYLPHENIDATE"] == PUBCHEM_CID_OVERRIDES["METHYLPHENIDATE"]:
        failures.append("CID override: DEXMETHYLPHENIDATE and METHYLPHENIDATE are not distinct")

    if failures:
        print("SELF-TEST FAILED", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1
    print(f"SELF-TEST PASSED: {len(canonical_cases) + len(exclusion_cases) + 12} checks", file=sys.stderr)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create fda_approved_drugs.csv with active ingredients, isomeric SMILES, InChIKey, and source."
    )
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run offline canonicalization/filter/QC/audit sanity checks and exit.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=DEFAULT_INPUT_DIR,
        help=f"Folder containing FDA txt files. Default: {DEFAULT_INPUT_DIR}",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("fda_approved_drugs.csv"),
        help="Output CSV path. Default: fda_approved_drugs.csv",
    )
    parser.add_argument(
        "--review-output",
        type=Path,
        default=None,
        help=(
            "Optional audit/review CSV. Default: <output_stem>_review.csv. "
            "Lists canonicalized aliases and any structure/name collisions."
        ),
    )
    parser.add_argument(
        "--allow-qc-warnings",
        action="store_true",
        help=(
            "Write output even when hard QC collisions are detected and exit 0. "
            "affected rows are withheld from the final CSV either way."
        ),
    )
    parser.add_argument(
        "--cache",
        type=Path,
        default=Path("structure_lookup_cache.json"),
        help="JSON cache for API lookups. Default: structure_lookup_cache.json.",
    )
    parser.add_argument(
        "--openfda-api-key",
        default=os.getenv("OPENFDA_API_KEY"),
        help="Optional openFDA API key. Also read from OPENFDA_API_KEY.",
    )
    parser.add_argument(
        "--keep-combination-ingredients",
        action="store_true",
        help="Do not split multi-ingredient combination active ingredients.",
    )
    parser.add_argument(
        "--max-molecular-weight",
        type=float,
        default=500.0,
        help=(
            "Maximum exact molecular weight in Daltons after salt/"
            "multicomponent standardization. Default: 500.0."
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Process only the first N unique ingredients, useful for testing.",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.30,
        help="Minimum delay between API requests in seconds. Default: 0.30.",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=30.0,
        help="HTTP timeout in seconds. Default: 30.",
    )
    parser.add_argument("--no-openfda", action="store_true", help="Skip FDA/openFDA GSRS lookup.")
    parser.add_argument("--no-pubchem", action="store_true", help="Skip PubChem fallback lookup.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.self_test:
        return run_self_test()
    if not args.no_openfda and not args.openfda_api_key:
        print(
            "Warning: no openFDA API key detected. Keyless use is limited; "
            "set OPENFDA_API_KEY for a full multi-query run.",
            file=sys.stderr,
        )

    product_rows, prelookup_audit = load_product_ingredients(args.input_dir, args)
    unique_ingredients = sorted({ingredient for _, ingredient in product_rows}, key=key_text)
    if args.limit:
        unique_ingredients = unique_ingredients[: args.limit]
        selected_set = set(unique_ingredients)
        product_rows = [
            row for row in product_rows if row[1] in selected_set
        ]
        # Test-mode audit scope: report only the selected post-prefilter
        # ingredients.  Full runs (no --limit) still include the complete
        # prelookup-exclusion audit.
        prelookup_audit_for_output: list[dict[str, str]] = []
        print(
            f"Audit scope: --limit {args.limit} active; review CSV is restricted "
            f"to the {len(unique_ingredients)} selected post-prefilter ingredients. "
            "Global prelookup exclusions are summarized in the log only.",
            file=sys.stderr,
        )
    else:
        prelookup_audit_for_output = prelookup_audit

    cache = load_cache(args.cache)
    client = JsonClient(delay_seconds=args.delay, timeout_seconds=args.timeout)

    total = len(unique_ingredients)
    for index, ingredient in enumerate(unique_ingredients, start=1):
        cache_key = key_text(ingredient)

        if cache_key not in cache:
            # If this ingredient redirects to a prodrug ester's active parent
            # (canonical_query_name), and that PARENT has already been looked
            # up and cached under its own name earlier in this run, reuse that
            # result directly instead of spending a second API call.
            canonical = canonical_query_name(ingredient)
            canonical_key = key_text(canonical)
            if canonical_key != cache_key and canonical_key in cache:
                cache[cache_key] = cache[canonical_key]
            else:
                cache[cache_key] = lookup_structure(
                    ingredient, client, args.openfda_api_key, args
                )

                # Alias the result under the parent's own cache key too, so a
                # LATER ingredient redirecting to the same parent can reuse it.
                if canonical_key != cache_key:
                    cache.setdefault(canonical_key, cache[cache_key])

        if index == 1 or index % 25 == 0 or index == total:
            hit_count, miss_count, skip_count, usable_hit_count = selected_cache_statistics(
                unique_ingredients[:index], cache, args.max_molecular_weight
            )
            print(
                f"[{index}/{total}] selected ingredients: "
                f"{hit_count} hits, {miss_count} misses, {skip_count} skipped, "
                f"{usable_hit_count} usable",
                file=sys.stderr,
            )
            save_cache(args.cache, cache)

    save_cache(args.cache, cache)
    written, hard_qc_issues, disposition_counts = write_csv(
        args.output,
        product_rows,
        cache,
        args.max_molecular_weight,
        args.review_output,
        prelookup_audit_for_output,
    )

    hit_count, miss_count, skip_count, usable_hit_count = selected_cache_statistics(
        unique_ingredients, cache, args.max_molecular_weight
    )

    print(f"Input product/ingredient rows: {len(product_rows)}", file=sys.stderr)
    print(f"Unique active ingredients checked: {total}", file=sys.stderr)
    print(f"Structure hits: {hit_count}", file=sys.stderr)
    print(
        f"Usable hits at <= {args.max_molecular_weight:.2f} Da: {usable_hit_count}",
        file=sys.stderr,
    )
    print(f"Misses: {miss_count}", file=sys.stderr)
    print(f"Skipped: {skip_count}", file=sys.stderr)
    print(f"Unique compound rows written: {written}", file=sys.stderr)
    review_path = args.review_output or args.output.with_name(args.output.stem + "_review.csv")
    print(f"Output: {args.output}", file=sys.stderr)
    print(f"Review/audit: {review_path}", file=sys.stderr)
    if disposition_counts:
        print("Audit dispositions:", file=sys.stderr)
        for disposition, count in sorted(disposition_counts.items()):
            print(f"  {disposition}: {count}", file=sys.stderr)
    if hard_qc_issues:
        print(
            f"HARD QC: {hard_qc_issues} unresolved collision-affected ingredient(s) recorded in {review_path}. "
            "Affected rows were withheld from the final CSV.",
            file=sys.stderr,
        )
        if not args.allow_qc_warnings:
            return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
