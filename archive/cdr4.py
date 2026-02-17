import re
import csv
import requests

RCSB_GQL = "https://data.rcsb.org/graphql"

Q = """
query($id:String!){
  entry(entry_id:$id){
    polymer_entities{
      rcsb_polymer_entity{pdbx_description}
      polymer_entity_instances{
        rcsb_polymer_entity_instance_container_identifiers{
          asym_id auth_asym_id
        }
      }
    }
  }
}
"""

# --- regex rules (simple + ordered) ---
RE_TCR   = re.compile(r"\b(t[\-\s]?cell receptor|tcr)\b", re.I)
RE_ALPHA = re.compile(r"\balpha\b|\bα\b", re.I)
RE_BETA  = re.compile(r"\bbeta\b|\bβ\b", re.I)

RE_MHC = re.compile(r"\b(mhc|hla|histocompatibility)\b", re.I)

# Peptide/antigen is tricky because MHC descriptions can include “antigen”.
# So: only call PEPTIDE if NOT MHC and it mentions peptide/epitope.
RE_PEPTIDE = re.compile(r"\b(peptide|epitope)\b", re.I)
RE_ANTIGEN = re.compile(r"\bantigen(ic)?\b", re.I)

def classify(desc: str) -> str:
    d = desc or ""

    # TCR first
    if RE_TCR.search(d):
        if RE_ALPHA.search(d) and not RE_BETA.search(d):
            return "TCR_ALPHA"
        if RE_BETA.search(d) and not RE_ALPHA.search(d):
            return "TCR_BETA"
        return "TCR"

    # MHC next
    if RE_MHC.search(d):
        return "MHC"

    # Peptide/antigen last (guard against MHC+antigen strings)
    if RE_PEPTIDE.search(d) or (RE_ANTIGEN.search(d) and not RE_MHC.search(d)):
        return "PEPTIDE"

    return "OTHER"

def fetch_entry(pdb_id: str):
    r = requests.post(RCSB_GQL, json={"query": Q, "variables": {"id": pdb_id}}, timeout=30)
    r.raise_for_status()
    return r.json()

def pdb_chain_rows(pdb_id: str):
    j = fetch_entry(pdb_id)
    ents = (j.get("data", {}).get("entry") or {}).get("polymer_entities") or []
    rows = []

    for e in ents:
        desc = ((e.get("rcsb_polymer_entity") or {}).get("pdbx_description")) or ""
        label = classify(desc)

        for inst in (e.get("polymer_entity_instances") or []):
            ids = (inst.get("rcsb_polymer_entity_instance_container_identifiers") or {})
            rows.append({
                "pdb_id": pdb_id,
                "chain_auth": ids.get("auth_asym_id") or "",
                "chain_asym": ids.get("asym_id") or "",
                "label": label,
                "description": desc,
            })
    return rows

def write_mapping_csv(pdb_ids, out_csv="rcsb_chain_mapping.csv"):
    all_rows = []
    for pdb in pdb_ids:
        pdb = pdb.strip().upper()
        if not pdb:
            continue
        try:
            all_rows.extend(pdb_chain_rows(pdb))
        except Exception as e:
            all_rows.append({
                "pdb_id": pdb,
                "chain_auth": "",
                "chain_asym": "",
                "label": "ERROR",
                "description": str(e),
            })

    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pdb_id", "chain_auth", "chain_asym", "label", "description"])
        w.writeheader()
        w.writerows(all_rows)

    return out_csv

# ---- example ----
if __name__ == "__main__":
    pdb_ids = ["9EJH", "1FYT", "2WBJ"]
    out = write_mapping_csv(pdb_ids, "rcsb_chain_mapping.csv")
    print("Wrote:", out)
