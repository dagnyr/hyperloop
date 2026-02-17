from Bio.PDB import PDBList, PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import PPBuilder
import os
import csv

# --------- USER SETTINGS ---------
pdb_csv = "in/classII_complexes.csv"
pdb_column = "PDB ID"

# RCSB chain mapping CSV (comma-separated) with columns:
# pdb_id,chain_auth,chain_asym,label,description
mapping_csv = "in/rcsb_chain_mapping.csv"
mapping_chain_column = "chain_auth"   # change to "chain_asym" if Biopython chain IDs match asym IDs instead
mapping_label_column = "label"

salt_bridge_cutoff_A = 4.0
include_histidine = False

# If a PDB is missing from mapping_csv, either skip it or fall back to defaults
skip_if_no_mapping = True
default_tcr_chain_ids = ["A", "B"]
default_peptide_chain_ids = ["E"]  # <-- NEW default (adjust to your common case)
# ---------------------------------

# ---- directories ----
pdb_dir = "pdb"
out_dir = "out"
os.makedirs(pdb_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)

# ---- load pdb ids from main CSV ----
pdb_ids = []
with open(pdb_csv, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        pdb_ids.append(row[pdb_column].strip().upper())

# ---- load chain mapping CSV into a dict: pdb_id -> {"tcr": [...], "peptide": [...]} ----
chain_map = {}
with open(mapping_csv, newline="") as f:
    reader = csv.DictReader(f)  # comma-separated by default

    # normalize headers once (handles accidental whitespace like "pdb_id ")
    if reader.fieldnames:
        reader.fieldnames = [h.strip() for h in reader.fieldnames]

    for row in reader:
        # normalize keys too, in case whitespace existed
        row = {k.strip(): v for k, v in row.items()}

        pid = row["pdb_id"].strip().upper()
        chain_id = row[mapping_chain_column].strip()
        label = row[mapping_label_column].strip().upper()

        if pid not in chain_map:
            chain_map[pid] = {"tcr": set(), "peptide": set()}

        if label.startswith("TCR_"):  # TCR_ALPHA, TCR_BETA, etc.
            chain_map[pid]["tcr"].add(chain_id)
        elif label in {"PEPTIDE", "ANTIGEN", "EPITOPE"}:
            chain_map[pid]["peptide"].add(chain_id)

# convert sets -> sorted lists
for pid in chain_map:
    chain_map[pid]["tcr"] = sorted(chain_map[pid]["tcr"])
    chain_map[pid]["peptide"] = sorted(chain_map[pid]["peptide"])

# Charged sidechain atoms (heavy atoms)
acid_atoms = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
base_atoms = {"LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"]}
if include_histidine:
    base_atoms["HIS"] = ["ND1", "NE2"]

# ---- combined output files ----
positions_all_path = os.path.join(out_dir, "tcr_positions_all.csv")
bridges_all_path   = os.path.join(out_dir, "tcr_peptide_saltbridges_all.tsv")  # <-- renamed
fasta_all_path     = os.path.join(out_dir, "tcr_sequences.fa")

write_pos_header = not os.path.exists(positions_all_path)
write_br_header  = not os.path.exists(bridges_all_path)

pos_all = open(positions_all_path, "a")
br_all  = open(bridges_all_path, "a")
fa_all  = open(fasta_all_path, "a")

if write_pos_header:
    pos_all.write("pdb_id,chain,seq_pos,resname,has_salt_bridge\n")
if write_br_header:
    br_all.write("pdb_id\ttcr_chain\ttcr_resid\tpeptide_chain\tpeptide_resid\tmin_distance\n")

# ---- loop over pdbs ----
for pdb_id in pdb_ids:

    try:
        # choose chains for this PDB from mapping (or fallback)
        if pdb_id in chain_map:
            tcr_chain_ids = chain_map[pdb_id]["tcr"]
            peptide_chain_ids = chain_map[pdb_id]["peptide"]
        else:
            if skip_if_no_mapping:
                print(f"{pdb_id}: skipping (no chain mapping found in {mapping_csv})")
                continue
            tcr_chain_ids = default_tcr_chain_ids
            peptide_chain_ids = default_peptide_chain_ids

        if not tcr_chain_ids or not peptide_chain_ids:
            print(f"{pdb_id}: skipping (empty chain mapping) tcr={tcr_chain_ids} peptide={peptide_chain_ids}")
            continue

        # load pdb
        pdb_path = PDBList().retrieve_pdb_file(pdb_id.lower(), pdir=pdb_dir, file_format="pdb")
        structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_path)
        model = next(structure.get_models())

        available_chains = {c.id for c in model}
        tcr_chains_present     = [c for c in tcr_chain_ids if c in available_chains]
        peptide_chains_present = [c for c in peptide_chain_ids if c in available_chains]

        if not tcr_chains_present:
            print(f"{pdb_id}: no mapped TCR chains present. mapped={tcr_chain_ids} available={sorted(available_chains)}")
            continue
        if not peptide_chains_present:
            print(f"{pdb_id}: no mapped PEPTIDE chains present. mapped={peptide_chain_ids} available={sorted(available_chains)}")
            continue

        # get all charged atoms in TCR and peptide
        tcr_charged_atoms = []
        peptide_charged_atoms = []

        def add_charged_atoms(chain_id, target_list):
            chain = model[chain_id]
            for residue in chain:
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

        for chain_id in peptide_chains_present:
            add_charged_atoms(chain_id, peptide_charged_atoms)

        # neighbor search (only TCR + peptide atoms)
        all_atoms = [x["atom"] for x in (tcr_charged_atoms + peptide_charged_atoms)]
        ns = NeighborSearch(all_atoms)
        peptide_atom_lookup = {id(x["atom"]): x for x in peptide_charged_atoms}

        salt_bridge_pairs = {}
        tcr_residues_with_saltbridge = {cid: [] for cid in tcr_chains_present}

        def add_residue_once(chain_id, residue):
            if residue not in tcr_residues_with_saltbridge[chain_id]:
                tcr_residues_with_saltbridge[chain_id].append(residue)

        # find salt bridges: TCR ↔ peptide only
        for tcr_item in tcr_charged_atoms:
            tcr_atom   = tcr_item["atom"]
            tcr_res    = tcr_item["residue"]
            tcr_chain  = tcr_item["chain"]
            tcr_charge = tcr_item["charge"]

            for nb_atom in ns.search(tcr_atom.coord, salt_bridge_cutoff_A, level="A"):
                if nb_atom is tcr_atom:
                    continue

                pep_item = peptide_atom_lookup.get(id(nb_atom))
                if pep_item is None:
                    continue

                if pep_item["charge"] == tcr_charge:
                    continue

                pep_res   = pep_item["residue"]
                pep_chain = pep_item["chain"]

                if pep_res is tcr_res:
                    continue

                distance = float(tcr_atom - nb_atom)

                key = (tcr_chain, tcr_res.get_id(), pep_chain, pep_res.get_id())
                if (key not in salt_bridge_pairs) or (distance < salt_bridge_pairs[key]):
                    salt_bridge_pairs[key] = distance

                add_residue_once(tcr_chain, tcr_res)

        # export positions + FASTA (still for TCR chains, but now mask = bridges to peptide)
        ppb = PPBuilder()
        for chain_id in tcr_chains_present:
            chain = model[chain_id]
            peptides = ppb.build_peptides(chain)
            if not peptides:
                continue

            sequence = "".join(str(p.get_sequence()) for p in peptides)
            residues_in_seq_order = []
            for p in peptides:
                residues_in_seq_order.extend(list(p))

            mask_chars = []
            for i, residue in enumerate(residues_in_seq_order, start=1):
                has_bridge = residue in tcr_residues_with_saltbridge[chain_id]
                mask_chars.append("*" if has_bridge else ".")

                pos_all.write(
                    f"{pdb_id},{chain_id},{i},{residue.get_resname()},{int(has_bridge)}\n"
                )

            fa_all.write(f">{pdb_id}|TCR_chain={chain_id}\n")
            fa_all.write(sequence + "\n")
            fa_all.write(f">{pdb_id}|TCR_chain={chain_id}|saltbridge_mask\n")
            fa_all.write("".join(mask_chars) + "\n")

        # export bridges (TCR ↔ peptide)
        for (tc, tr, pc, pr), d in sorted(salt_bridge_pairs.items(), key=lambda x: x[1]):
            br_all.write(f"{pdb_id}\t{tc}\t{tr}\t{pc}\t{pr}\t{d:.3f}\n")

        print(f"{pdb_id}: {len(salt_bridge_pairs)} salt-bridge residue pairs (TCR↔peptide)")

    except Exception as e:
        print(f"{pdb_id}: ERROR -> {type(e).__name__}: {e}")

pos_all.close()
br_all.close()
fa_all.close()

print(f"\nWrote combined files:")
print(f"  {positions_all_path}")
print(f"  {bridges_all_path}")
print(f"  {fasta_all_path}")
