# final csv

#!/usr/bin/env python3
import os
import re
import pandas as pd

# ----------------- USER PATHS -----------------
CHAIN_MAP_PATH  = "in/rcsb_chain_mapping.csv"          # chain mapping (excel paste ok)
COMPLEXES_PATH  = "in/classII_complexes.csv"           # class II complexes table
SALTBRIDGE_PATH = "tcr-pmhc/tcr_pmhc_saltbridges_all.tsv"               # your tcr-peptide salt bridge output
OUT_PATH        = "out/saltbridges_with_metadata_tcr_mphc.csv"  # output
# ----------------------------------------------

def read_table_autodetect(path: str) -> pd.DataFrame:
    """Read CSV/TSV by sniffing delimiter (works for Excel tab-pastes)."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"File not found: {path}")
    df = pd.read_csv(path, sep=None, engine="python", dtype=str)
    df.columns = [c.strip() for c in df.columns]
    for c in df.columns:
        df[c] = df[c].astype(str).str.strip()
    return df

def normalize_colname(c: str) -> str:
    c2 = c.strip().replace("<BR>", " ")
    c2 = re.sub(r"\s+", " ", c2)
    return c2

def normalize_pdb_id(x: str) -> str:
    return str(x).strip().upper()

def parse_resseq(resid_str: str):
    """
    Your resid looks like: "(' ', 30, ' ')" or similar.
    Extract the middle integer (resseq). Return None if not parseable.
    """
    s = str(resid_str)
    # grab the first integer we see
    m = re.search(r"(-?\d+)", s)
    return int(m.group(1)) if m else None

def find_pdb_col(df: pd.DataFrame) -> str:
    # Accept common variants: "PDB ID", "pdb_id", etc.
    for c in df.columns:
        key = c.lower().replace(" ", "").replace("_", "")
        if key in {"pdbid", "pdb"}:
            return c
    if "pdb_id" in df.columns:
        return "pdb_id"
    raise ValueError(f"Could not find a PDB ID column. Columns: {list(df.columns)}")

def main():
    os.makedirs(os.path.dirname(OUT_PATH) or ".", exist_ok=True)

    # ---------- Load chain mapping ----------
    chain_df = read_table_autodetect(CHAIN_MAP_PATH)
    chain_df.columns = [normalize_colname(c) for c in chain_df.columns]

    required_chain_cols = ["pdb_id", "chain_auth", "chain_asym", "label", "description"]
    missing = [c for c in required_chain_cols if c not in chain_df.columns]
    if missing:
        raise ValueError(
            f"{CHAIN_MAP_PATH} is missing columns: {missing}\n"
            f"Found columns: {list(chain_df.columns)}"
        )

    chain_df["pdb_id"] = chain_df["pdb_id"].map(normalize_pdb_id)
    chain_df["chain_auth"] = chain_df["chain_auth"].astype(str).str.strip()
    chain_df["chain_asym"] = chain_df["chain_asym"].astype(str).str.strip()
    chain_df["label"] = chain_df["label"].astype(str).str.strip().str.upper()
    chain_df["description"] = chain_df["description"].astype(str).str.strip()

    # We'll join on chain_auth by default (since your saltbridge file uses single-letter chain IDs)
    chain_key = "chain_auth"

    # ---------- Load complexes metadata ----------
    comp_df = read_table_autodetect(COMPLEXES_PATH)
    comp_df.columns = [normalize_colname(c) for c in comp_df.columns]

    pdb_col = find_pdb_col(comp_df)
    comp_df[pdb_col] = comp_df[pdb_col].map(normalize_pdb_id)
    comp_df = comp_df.rename(columns={pdb_col: "pdb_id"})

    # ---------- Load salt bridges ----------
    sb_df = read_table_autodetect(SALTBRIDGE_PATH)
    sb_df.columns = [normalize_colname(c) for c in sb_df.columns]

    required_sb_cols = ["pdb_id", "tcr_chain", "tcr_resid", "pmhc_chain", "pmhc_resid", "min_distance"]
    missing = [c for c in required_sb_cols if c not in sb_df.columns]
    if missing:
        raise ValueError(
            f"{SALTBRIDGE_PATH} is missing columns: {missing}\n"
            f"Found columns: {list(sb_df.columns)}"
        )

    sb_df["pdb_id"] = sb_df["pdb_id"].map(normalize_pdb_id)
    sb_df["tcr_chain"] = sb_df["tcr_chain"].astype(str).str.strip()
    sb_df["pmhc_chain"] = sb_df["pmhc_chain"].astype(str).str.strip()

    # Parse residue numbers for easier downstream work
    sb_df["tcr_resseq"] = sb_df["tcr_resid"].apply(parse_resseq)
    sb_df["peptide_resseq"] = sb_df["pmhc_resid"].apply(parse_resseq)

    # Ensure min_distance numeric if possible
    sb_df["min_distance"] = pd.to_numeric(sb_df["min_distance"], errors="coerce")

    # ---------- Join chain mapping for TCR chain ----------
    tcr_map = chain_df.rename(columns={
        chain_key: "tcr_chain",
        "chain_asym": "tcr_chain_asym",
        "label": "tcr_chain_label",
        "description": "tcr_chain_description",
    })[["pdb_id", "tcr_chain", "tcr_chain_asym", "tcr_chain_label", "tcr_chain_description"]]

    merged = sb_df.merge(tcr_map, on=["pdb_id", "tcr_chain"], how="left")

    # ---------- Join chain mapping for PEPTIDE chain ----------
    pep_map = chain_df.rename(columns={
        chain_key: "pmhc_chain",
        "chain_asym": "pmhc_chain_asym",
        "label": "pmhc_chain_label",
        "description": "pmhc_chain_description",
    })[["pdb_id", "pmhc_chain", "pmhc_chain_asym", "pmhc_chain_label", "pmhc_chain_description"]]

    merged = merged.merge(pep_map, on=["pdb_id", "pmhc_chain"], how="left")

    # ---------- Join complexes metadata (structure-level info) ----------
    merged = merged.merge(comp_df, on="pdb_id", how="left")

    # ---------- Nice column order ----------
    first = [
        "pdb_id",
        "tcr_chain", "tcr_chain_asym", "tcr_chain_label", "tcr_chain_description",
        "tcr_resid", "tcr_resseq",
        "pmhc_chain", "pmhc_chain_asym", "pmhc_chain_label", "pmhc_chain_description",
        "pmhc_resid", "peptide_resseq",
        "min_distance",
    ]
    cols = [c for c in first if c in merged.columns] + [c for c in merged.columns if c not in first]
    merged = merged[cols]

    # ---------- Basic QC warnings ----------
    # If the mapping didn't match some chains, labels/descriptions will be NaN/None-like strings.
    # Print a short report.
    missing_tcr = merged["tcr_chain_label"].isna().sum() if "tcr_chain_label" in merged.columns else 0
    missing_pep = merged["pmhc_chain_label"].isna().sum() if "pmhc_chain_label" in merged.columns else 0
    if missing_tcr or missing_pep:
        print("WARNING: Some salt-bridge chains did not match the chain-mapping file.")
        if missing_tcr:
            print(f"  Missing TCR chain mapping rows: {missing_tcr}")
        if missing_pep:
            print(f"  Missing peptide chain mapping rows: {missing_pep}")
        print("Check whether your mapping uses chain_auth vs chain_asym, or if chain IDs differ.\n")

    # ---------- Write ----------
    merged.to_csv(OUT_PATH, index=False)
    print(f"Wrote: {OUT_PATH}")
    print(f"Salt-bridge rows: {len(merged)}")

if __name__ == "__main__":
    main()
