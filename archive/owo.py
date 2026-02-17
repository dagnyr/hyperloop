from Bio.PDB import PDBList, PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import PPBuilder
import os
import csv

# --------- USER SETTINGS ---------
pdb_csv = "in/classII_complexes.csv"      # CSV containing PDB IDs
pdb_column = "PDB ID"                    # column name containing the IDs

tcr_chain_ids  = ["A", "B"]
pmhc_chain_ids = ["C", "D", "E"]
salt_bridge_cutoff_A = 4.0

include_histidine = False

# ---- directories ----
pdb_dir = "pdb"     # downloaded PDBs go here
out_dir = "out"     # output files go here
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
fasta_all_path     = os.path.join(out_dir, "tcr_sequences.fa")   # <-- NEW

write_pos_header = not os.path.exists(positions_all_path)
write_br_header  = not os.path.exists(bridges_all_path)

pos_all = open(positions_all_path, "a")
br_all  = open(bridges_all_path, "a")
fa_all  = open(fasta_all_path, "a")  # <-- NEW

if write_pos_header:
    pos_all.write("pdb_id,chain,seq_pos,resname,has_salt_bridge\n")
if write_br_header:
    br_all.write("pdb_id\ttcr_chain\ttcr_resid\tpmhc_chain\tpmhc_resid\tmin_distance\n")

# ---- loop over pdbs ----
for pdb_id in pdb_ids:

    try:
        # load pdb
        pdb_path = PDBList().retrieve_pdb_file(pdb_id.lower(), pdir=pdb_dir, file_format="pdb")
        structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_path)
        model = next(structure.get_models())

        available_chains = {c.id for c in model}
        tcr_chains_present  = [c for c in tcr_chain_ids if c in available_chains]
        pmhc_chains_present = [c for c in pmhc_chain_ids if c in available_chains]

        if not tcr_chains_present:
            print(f"{pdb_id}: no TCR chains found from {tcr_chain_ids}. available={sorted(available_chains)}")
            continue
        if not pmhc_chains_present:
            print(f"{pdb_id}: no pMHC chains found from {pmhc_chain_ids}. available={sorted(available_chains)}")
            continue

        # get all charged atoms in TCR and pMHC
        tcr_charged_atoms = []
        pmhc_charged_atoms = []

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

        for chain_id in pmhc_chains_present:
            add_charged_atoms(chain_id, pmhc_charged_atoms)

        # neighbor search
        all_atoms = [x["atom"] for x in (tcr_charged_atoms + pmhc_charged_atoms)]
        ns = NeighborSearch(all_atoms)
        pmhc_atom_lookup = {id(x["atom"]): x for x in pmhc_charged_atoms}

        salt_bridge_pairs = {}

        # NOTE: use tcr_chains_present here so lookups never KeyError
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

        # export positions + FASTA (append to combined files)
        ppb = PPBuilder()
        for chain_id in tcr_chains_present:
            chain = model[chain_id]
            peptides = ppb.build_peptides(chain)
            if not peptides:
                continue

            # build sequence + residue list in the same order
            sequence = "".join(str(p.get_sequence()) for p in peptides)
            residues_in_seq_order = []
            for p in peptides:
                residues_in_seq_order.extend(list(p))

            # mask + positions CSV
            mask_chars = []
            for i, residue in enumerate(residues_in_seq_order, start=1):
                has_bridge = residue in tcr_residues_with_saltbridge[chain_id]
                mask_chars.append("*" if has_bridge else ".")

                pos_all.write(
                    f"{pdb_id},{chain_id},{i},{residue.get_resname()},{int(has_bridge)}\n"
                )

            # NEW: write FASTA to file (sequence + mask)
            fa_all.write(f">{pdb_id}|TCR_chain={chain_id}\n")
            fa_all.write(sequence + "\n")
            fa_all.write(f">{pdb_id}|TCR_chain={chain_id}|saltbridge_mask\n")
            fa_all.write("".join(mask_chars) + "\n")

        # export bridges (append to one combined TSV)
        for (tc, tr, pc, pr), d in sorted(salt_bridge_pairs.items(), key=lambda x: x[1]):
            br_all.write(f"{pdb_id}\t{tc}\t{tr}\t{pc}\t{pr}\t{d:.3f}\n")

        print(f"{pdb_id}: {len(salt_bridge_pairs)} salt-bridge residue pairs")

    except Exception as e:
        print(f"{pdb_id}: ERROR -> {type(e).__name__}: {e}")

pos_all.close()
br_all.close()
fa_all.close()

print(f"\nWrote combined files:")
print(f"  {positions_all_path}")
print(f"  {bridges_all_path}")
print(f"  {fasta_all_path}")
