from Bio.PDB import PDBList, PDBParser, NeighborSearch
import os
import csv

# SIFTS + UniProt
import gzip
import io
import xml.etree.ElementTree as ET
import requests
import re

# --------- USER SETTINGS ---------
pdb_csv = "in/classII_complexes.csv"
pdb_column = "PDB ID"

salt_bridge_cutoff_A = 4.0
include_histidine = False

# ---- directories ----
pdb_dir = "pdb"
out_dir = "out"
os.makedirs(pdb_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)

# ---- load pdb ids from CSV ----
pdb_ids = []
with open(pdb_csv, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        pdb_ids.append(row[pdb_column].strip())

# Charged sidechain atoms (heavy atoms)
acid_atoms = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
base_atoms = {"LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"]}
if include_histidine:
    base_atoms["HIS"] = ["ND1", "NE2"]

# ---- combined output files ----
positions_all_path = os.path.join(out_dir, "tcr_positions_all.csv")
bridges_all_path   = os.path.join(out_dir, "tcr_pmhc_saltbridges_all.tsv")
fasta_all_path     = os.path.join(out_dir, "tcr_sequences.fa")
chainmap_all_path  = os.path.join(out_dir, "chain_map.tsv")

write_pos_header = not os.path.exists(positions_all_path)
write_br_header  = not os.path.exists(bridges_all_path)
write_cm_header  = not os.path.exists(chainmap_all_path)

pos_all = open(positions_all_path, "a")
br_all  = open(bridges_all_path, "a")
fa_all  = open(fasta_all_path, "a")
cm_all  = open(chainmap_all_path, "a")

if write_pos_header:
    pos_all.write("pdb_id,pdb_chain,uniprot_acc,uniprot_pos,resname,has_salt_bridge\n")
if write_br_header:
    br_all.write("pdb_id\ttcr_chain\ttcr_resid\tpmhc_chain\tpmhc_resid\tmin_distance\n")
if write_cm_header:
    cm_all.write("pdb_id\tpdb_chain\tuniprot_acc\trole\tatom_len\n")


# ---------------- SIFTS / UniProt helpers ----------------

def _strip_ns(tag: str) -> str:
    return tag.split("}", 1)[-1]

def _sifts_cache_path(pdb_id: str) -> str:
    pdb_id = pdb_id.lower()
    subdir = pdb_id[1:3]
    return os.path.join(pdb_dir, f"sifts_{pdb_id}_{subdir}.xml.gz")

def download_sifts_xml_gz(pdb_id: str) -> bytes:
    pdb_id = pdb_id.lower()
    cache_path = _sifts_cache_path(pdb_id)
    if os.path.exists(cache_path):
        with open(cache_path, "rb") as f:
            return f.read()

    subdir = pdb_id[1:3]
    url = f"https://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/{subdir}/{pdb_id}.xml.gz"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    data = r.content

    with open(cache_path, "wb") as f:
        f.write(data)
    return data

def parse_sifts_chain_map(xml_gz_bytes: bytes, chain_id: str):
    """
    Returns:
      uniprot_acc (str|None)
      pdb_to_unp (dict): key=(resseq:int, icode:str|None) -> unp_pos:int
    """
    chain_id = chain_id.upper()
    with gzip.GzipFile(fileobj=io.BytesIO(xml_gz_bytes)) as gz:
        xml_bytes = gz.read()
    root = ET.fromstring(xml_bytes)

    uniprot_acc = None
    pdb_to_unp = {}

    for residue in root.iter():
        if _strip_ns(residue.tag) != "residue":
            continue

        pdb_ref = None
        unp_ref = None

        for xref in residue:
            if _strip_ns(xref.tag) != "crossRefDb":
                continue
            db = xref.attrib.get("dbSource")

            if db == "PDB":
                if xref.attrib.get("dbChainId", "").upper() == chain_id:
                    pdb_ref = xref
            elif db == "UniProt":
                unp_ref = xref

        if pdb_ref is None or unp_ref is None:
            continue

        if uniprot_acc is None:
            uniprot_acc = unp_ref.attrib.get("dbAccessionId")

        pdb_resnum = pdb_ref.attrib.get("dbResNum")
        if pdb_resnum is None:
            continue
        try:
            resseq = int(pdb_resnum)
        except ValueError:
            continue

        icode = pdb_ref.attrib.get("dbInsCode")
        if icode == "" or icode is None:
            icode = None

        unp_resnum = unp_ref.attrib.get("dbResNum")
        if unp_resnum is None:
            continue
        try:
            unp_pos = int(unp_resnum)
        except ValueError:
            continue

        pdb_to_unp[(resseq, icode)] = unp_pos

    return uniprot_acc, pdb_to_unp

def fetch_uniprot_sequence(uniprot_acc: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}.fasta"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    seq_lines = [ln.strip() for ln in r.text.splitlines() if ln.strip() and not ln.startswith(">")]
    return "".join(seq_lines)

# ---- NEW: disk-cached UniProt TXT fetch (avoids Ctrl-C stalls) ----
uniprot_cache_dir = os.path.join(pdb_dir, "uniprot_txt")
os.makedirs(uniprot_cache_dir, exist_ok=True)

def fetch_uniprot_txt(uniprot_acc: str) -> str:
    cache_path = os.path.join(uniprot_cache_dir, f"{uniprot_acc}.txt")
    if os.path.exists(cache_path):
        with open(cache_path, "r") as f:
            return f.read()

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_acc}.txt"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    txt = r.text

    with open(cache_path, "w") as f:
        f.write(txt)
    return txt

# ---- FIXED role classification: use GN Name=... only; do NOT use "immunoglobulin" substring ----
GN_NAME_RE = re.compile(r"^GN\s+Name=([^;]+);", re.M)

def classify_role_from_uniprot(uniprot_acc: str | None) -> str:
    """
    Returns one of: TCR, MHC, IG, OTHER, NA

    Avoids the 'immunoglobulin-like domain' trap by never classifying IG
    from the word 'immunoglobulin' alone.
    """
    if not uniprot_acc:
        return "NA"

    try:
        txt = fetch_uniprot_txt(uniprot_acc)
    except Exception:
        return "OTHER"

    m = GN_NAME_RE.search(txt)
    gene = m.group(1).upper() if m else ""

    # MHC-II (human + mouse)
    if gene.startswith("HLA-D") or gene.startswith("H2-"):
        return "MHC"

    # TCR
    if gene.startswith("TR"):
        return "TCR"

    # Antibody / Ig chains (only if gene indicates it)
    if gene.startswith(("IGH", "IGK", "IGL")):
        return "IG"

    # fallback: explicit antibody word (safe)
    up = txt.upper()
    if "ANTIBODY" in up:
        return "IG"

    return "OTHER"


def atom_chain_len(chain) -> int:
    n = 0
    for r in chain:
        if r.get_id()[0] == " ":
            n += 1
    return n


# ---------------- MAIN LOOP ----------------

try:
    for pdb_id in pdb_ids:
        try:
            pdb_path = PDBList().retrieve_pdb_file(pdb_id.lower(), pdir=pdb_dir, file_format="pdb")
            structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_path)
            model = next(structure.get_models())

            # Download SIFTS once per PDB
            try:
                sifts_xml_gz = download_sifts_xml_gz(pdb_id)
            except Exception as e:
                print(f"{pdb_id}: ERROR could not fetch SIFTS -> {type(e).__name__}: {e}")
                continue

            # build chain -> uniprot + role map
            available_chains = [c.id for c in model]
            chain_uniprot = {}
            chain_role = {}
            chain_len = {}

            for chain_id in available_chains:
                chain_len[chain_id] = atom_chain_len(model[chain_id])

                uniprot_acc, _ = parse_sifts_chain_map(sifts_xml_gz, chain_id)
                chain_uniprot[chain_id] = uniprot_acc
                chain_role[chain_id] = classify_role_from_uniprot(uniprot_acc)

                cm_all.write(f"{pdb_id}\t{chain_id}\t{uniprot_acc or 'NA'}\t{chain_role[chain_id]}\t{chain_len[chain_id]}\n")

            # Define TCR and pMHC sides from roles
            tcr_chains_present = [cid for cid in available_chains if chain_role.get(cid) == "TCR"]
            mhc_chains = [cid for cid in available_chains if chain_role.get(cid) == "MHC"]

            # Add peptide-like short chains to pMHC side (but exclude TCR)
            peptide_chains = [cid for cid in available_chains if 5 <= chain_len.get(cid, 0) <= 40 and chain_role.get(cid) != "TCR"]
            pmhc_chains_present = list(dict.fromkeys(mhc_chains + peptide_chains))

            if not tcr_chains_present or not mhc_chains:
                diag = ", ".join([f"{cid}:{chain_uniprot.get(cid) or 'NA'}:{chain_role.get(cid) or 'NA'}:{chain_len.get(cid, 0)}"
                                  for cid in available_chains])
                print(f"{pdb_id}: missing TCR/MHC after UniProt classification. chains={diag}")
                continue

            # get all charged atoms in TCR and pMHC
            tcr_charged_atoms = []
            pmhc_charged_atoms = []

            def add_charged_atoms(chain_id, target_list):
                chain = model[chain_id]
                for residue in chain:
                    if residue.get_id()[0] != " ":
                        continue
                    resname = residue.get_resname()

                    if resname in acid_atoms:
                        for atom_name in acid_atoms[resname]:
                            if atom_name in residue:
                                target_list.append({
                                    "chain": chain_id,
                                    "residue": residue,
                                    "atom": residue[atom_name],
                                    "charge": "acid",
                                })

                    if resname in base_atoms:
                        for atom_name in base_atoms[resname]:
                            if atom_name in residue:
                                target_list.append({
                                    "chain": chain_id,
                                    "residue": residue,
                                    "atom": residue[atom_name],
                                    "charge": "base",
                                })

            for chain_id in tcr_chains_present:
                add_charged_atoms(chain_id, tcr_charged_atoms)

            for chain_id in pmhc_chains_present:
                add_charged_atoms(chain_id, pmhc_charged_atoms)

            all_atoms = [x["atom"] for x in (tcr_charged_atoms + pmhc_charged_atoms)]
            if not all_atoms:
                print(f"{pdb_id}: no charged atoms found")
                continue

            ns = NeighborSearch(all_atoms)
            pmhc_atom_lookup = {id(x["atom"]): x for x in pmhc_charged_atoms}

            salt_bridge_pairs = {}
            tcr_residues_with_saltbridge = {cid: [] for cid in tcr_chains_present}

            def add_residue_once(chain_id, residue):
                if residue not in tcr_residues_with_saltbridge[chain_id]:
                    tcr_residues_with_saltbridge[chain_id].append(residue)

            # find salt bridges
            for tcr_item in tcr_charged_atoms:
                tcr_atom   = tcr_item["atom"]
                tcr_res    = tcr_item["residue"]
                tcr_chain  = tcr_item["chain"]
                tcr_charge = tcr_item["charge"]

                for nb_atom in ns.search(tcr_atom.coord, salt_bridge_cutoff_A, level="A"):
                    if nb_atom is tcr_atom:
                        continue

                    pmhc_item = pmhc_atom_lookup.get(id(nb_atom))
                    if pmhc_item is None:
                        continue

                    if pmhc_item["charge"] == tcr_charge:
                        continue

                    pmhc_res   = pmhc_item["residue"]
                    pmhc_chain = pmhc_item["chain"]

                    if pmhc_res is tcr_res:
                        continue

                    distance = float(tcr_atom - nb_atom)

                    key = (tcr_chain, tcr_res.get_id(), pmhc_chain, pmhc_res.get_id())
                    if (key not in salt_bridge_pairs) or (distance < salt_bridge_pairs[key]):
                        salt_bridge_pairs[key] = distance

                    add_residue_once(tcr_chain, tcr_res)

            # export full UniProt FASTA + UniProt-coordinate mask (ONLY for role=TCR)
            for chain_id in tcr_chains_present:
                uniprot_acc, pdb_to_unp = parse_sifts_chain_map(sifts_xml_gz, chain_id)
                if not uniprot_acc or not pdb_to_unp:
                    print(f"{pdb_id} chain {chain_id}: WARN no UniProt mapping found in SIFTS")
                    continue

                try:
                    unp_seq = fetch_uniprot_sequence(uniprot_acc)
                except Exception as e:
                    print(f"{pdb_id} chain {chain_id}: WARN UniProt fetch failed ({uniprot_acc}) -> {type(e).__name__}: {e}")
                    continue

                mask_chars = ["."] * len(unp_seq)

                for residue in tcr_residues_with_saltbridge[chain_id]:
                    het, resseq, icode = residue.get_id()
                    if het != " ":
                        continue
                    icode = None if (icode == " " or icode == "" or icode is None) else icode
                    unp_pos = pdb_to_unp.get((resseq, icode))
                    if unp_pos is None:
                        continue
                    if 1 <= unp_pos <= len(unp_seq):
                        mask_chars[unp_pos - 1] = "*"
                        pos_all.write(f"{pdb_id},{chain_id},{uniprot_acc},{unp_pos},{residue.get_resname()},1\n")

                fa_all.write(f">{pdb_id}|pdb_chain={chain_id}|uniprot={uniprot_acc}|role=TCR\n")
                fa_all.write(unp_seq + "\n")
                fa_all.write(f">{pdb_id}|pdb_chain={chain_id}|uniprot={uniprot_acc}|role=TCR|saltbridge_mask_uniprot\n")
                fa_all.write("".join(mask_chars) + "\n")

            # export bridges
            for (tc, tr, pc, pr), d in sorted(salt_bridge_pairs.items(), key=lambda x: x[1]):
                br_all.write(f"{pdb_id}\t{tc}\t{tr}\t{pc}\t{pr}\t{d:.3f}\n")

            print(f"{pdb_id}: TCR chains={tcr_chains_present} MHC chains={mhc_chains} peptide-ish={peptide_chains} bridges={len(salt_bridge_pairs)}")

        except Exception as e:
            print(f"{pdb_id}: ERROR -> {type(e).__name__}: {e}")

except KeyboardInterrupt:
    print("\nInterrupted by user (Ctrl-C). Closing files cleanly.")

finally:
    pos_all.close()
    br_all.close()
    fa_all.close()
    cm_all.close()

print("\nWrote combined files:")
# print(f"  {positions_all_path}")
# print(f"  {bridges_all_path}")
# print(f"  {fasta_all_path}")
# print(f"  {chainmap_all_path}")
