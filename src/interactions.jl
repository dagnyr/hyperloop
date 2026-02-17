using BioStructures
using DataFrames
using LinearAlgebra
using ProgressMeter

### -------------------------
### Constants

const hbond_acceptors = Set(["ALAO", "ARGO", "ASNO", "ASPO", "CYSO", "GLNO", "GLUO", "GLYO", "HISO", "ILEO", "LEUO", "LYSO", "METO", "PHEO", "PROO", "SERO", "THRO", "TRPO", "TYRO", "VALO", "ASNOD1", "ASNND2", "ASPOD1", "ASPOD2", "GLNOE1", "GLNNE2", "GLUOE1", "GLUOE2", "HISND1", "HISCE1", "HISNE2", "HISCD2", "METSD", "CYSSG", "SEROG", "THROG1", "TYROH", "ALAOXT", "ARGOXT", "ASNOXT", "ASPOXT", "CYSOXT", "GLNOXT", "GLUOXT", "GLYOXT", "HISOXT", "ILEOXT", "LEUOXT", "LYSOXT", "METOXT", "PHEOXT", "PROOXT", "SEROXT", "THROXT", "TRPOXT", "TYROXT", "VALOXT"])

const hbond_donors = Set(["ALAN", "ARGN", "ASNN", "ASPN", "CYSN", "GLNN", "GLUN", "GLYN", "HISN", "ILEN", "LEUN", "LYSN", "METN", "PHEN", "SERN", "THRN", "TRPN", "TYRN", "VALN", "ARGNE", "ARGNH1", "ARGNH2", "ASNND2", "ASNOD1", "CYSSG", "GLNNE2", "GLNOE1", "HISND1", "HISCE1", "HISNE2", "HISCD2", "LYSNZ", "SEROG", "THROG1", "TRPNE1", "TYROH"])

const pos_ionisable = Set(["ARGNE", "ARGCZ", "ARGNH1", "ARGNH2", "HISCG", "HISND1", "HISCE1", "HISNE2", "HISCD2", "LYSNZ"])

const neg_ionisable = Set(["ASPOD1", "ASPOD2", "GLUOE1", "GLUOE2"])

const hydrophobes = Set(["ALACB", "ARGCB", "ARGCG", "ASNCB", "ASPCB", "CYSCB", "GLNCB", "GLNCG", "GLUCB", "GLUCG", "HISCB", "ILECB", "ILECG1", "ILECD1", "ILECG2", "LEUCB", "LEUCG", "LEUCD1", "LEUCD2", "LYSCB", "LYSCG", "LYSCD", "METCB", "METCG", "METSD", "METCE", "PHECB", "PHECG", "PHECD1", "PHECD2", "PHECE1", "PHECE2", "PHECZ", "PROCB", "PROCG", "THRCG2", "TRPCB", "TRPCG", "TRPCD2", "TRPCE3", "TRPCZ3", "TRPCH2", "TRPCZ2", "TYRCG", "TYRCD1", "TYRCD2", "TYRCE1", "TYRCE2", "VALCB", "VALCG1", "VALCG2"])

const aromatic_residues = Set(["PHE", "TYR", "TRP", "HIS"])

### -------------------------
### Helpers

isheavy(at) = BioStructures.element(at) != "H"

function log_interaction!(db, pdb_id, int_type, dist, idA, classA, resA, idB, classB, resB)
    push!(db, (
        pdb_id, int_type, dist,
        idA, classA, resname(resA), resnumber(resA),
        idB, classB, resname(resB), resnumber(resB)
    ))
end

### -------------------------
### Geometry

function get_ring_atoms(res)
    valid_names = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "NE1", "CE3", "CZ2", "CZ3", "CH2", "ND1", "NE2"]
    return [at for at in res if atomname(at) in valid_names && isheavy(at)]
end

function get_centroid(ring_atoms)
    if length(ring_atoms) == 0
        return nothing
    end
    c_list = [BioStructures.coords(at) for at in ring_atoms]
    return sum(c_list) / length(c_list)
end

function get_normal_vector(ring_atoms)
    if length(ring_atoms) < 3
        return nothing
    end
    p1 = BioStructures.coords(ring_atoms[1])
    p2 = BioStructures.coords(ring_atoms[2])
    p3 = BioStructures.coords(ring_atoms[3])

    normal = cross(p2 - p1, p3 - p1)
    n_len = norm(normal)
    return n_len > 0 ? normal / n_len : nothing
end

function get_inter_ring_angle(normA, normB)
    cos_theta = clamp(dot(normA, normB), -1.0, 1.0)
    return acos(abs(cos_theta)) * (180 / π)
end

### -------------------------
### Interactions

@inline atom_id_str(res, atom) = string(resname(res), atomname(atom))

function check_gliph_contact!(db, pdb_id, dist, idA, classA, resA, idB, classB, resB)
    if dist > 5.0
        return false
    end
    log_interaction!(db, pdb_id, "GLIPH Contact", dist, idA, classA, resA, idB, classB, resB)
    return true
end

function check_salt_bridge!(db, pdb_id, dist, id_str_A, id_str_B, idA, classA, resA, idB, classB, resB)
    if dist > 4.0
        return false
    end
    if (id_str_A in pos_ionisable && id_str_B in neg_ionisable) || (id_str_B in pos_ionisable && id_str_A in neg_ionisable)
        log_interaction!(db, pdb_id, "Salt Bridge", dist, idA, classA, resA, idB, classB, resB)
        return true
    end
    return false
end

function check_hydrophobic!(db, pdb_id, dist, id_str_A, id_str_B, idA, classA, resA, idB, classB, resB)
    if dist > 4.0
        return false
    end
    if (id_str_A in hydrophobes) && (id_str_B in hydrophobes)
        log_interaction!(db, pdb_id, "Hydrophobic", dist, idA, classA, resA, idB, classB, resB)
        return true
    end
    return false
end

function check_polar_hbond!(db, pdb_id, dist, id_str_A, id_str_B, idA, classA, resA, idB, classB, resB)
    if dist > 3.5
        return false
    end
    if (id_str_A in hbond_donors && id_str_B in hbond_acceptors) || (id_str_B in hbond_donors && id_str_A in hbond_acceptors)
        log_interaction!(db, pdb_id, "Polar H-Bond", dist, idA, classA, resA, idB, classB, resB)
        return true
    end
    return false
end

function check_atom_chemistry!(db, pdb_id, atomA, atomB, resA, resB, idA, classA, idB, classB)
    dist = BioStructures.distance(atomA, atomB)
    if dist > 5.0
        return
    end

    check_gliph_contact!(db, pdb_id, dist, idA, classA, resA, idB, classB, resB)

    id_str_A = atom_id_str(resA, atomA)
    id_str_B = atom_id_str(resB, atomB)

    check_salt_bridge!(db, pdb_id, dist, id_str_A, id_str_B, idA, classA, resA, idB, classB, resB)
    check_hydrophobic!(db, pdb_id, dist, id_str_A, id_str_B, idA, classA, resA, idB, classB, resB)
    check_polar_hbond!(db, pdb_id, dist, id_str_A, id_str_B, idA, classA, resA, idB, classB, resB)
end

function check_pi_stacking!(db, pdb_id, resA, resB, idA, classA, idB, classB)
    if !(resname(resA) in aromatic_residues && resname(resB) in aromatic_residues)
        return
    end

    ringA = get_ring_atoms(resA)
    ringB = get_ring_atoms(resB)

    centA = get_centroid(ringA)
    centB = get_centroid(ringB)
    if isnothing(centA) || isnothing(centB)
        return
    end

    dist = norm(centA - centB)
    if dist > 6.0
        return
    end

    normA = get_normal_vector(ringA)
    normB = get_normal_vector(ringB)
    if isnothing(normA) || isnothing(normB)
        return
    end

    angle = get_inter_ring_angle(normA, normB)

    if angle <= 30.0
        log_interaction!(db, pdb_id, "pi-pi (Parallel)", dist, idA, classA, resA, idB, classB, resB)
    elseif 60.0 <= angle <= 90.0
        log_interaction!(db, pdb_id, "pi-pi (T-Shaped)", dist, idA, classA, resA, idB, classB, resB)
    end
end

### -------------------------
### Loading and Main Logic

function get_annotations(pdb_id::String)
    println("Downloading structure and gathering annotations for $pdb_id...")

    filepath = downloadpdb(pdb_id, format=MMCIFFormat) # load structure
    mmcif_dict = MMCIFDict(filepath)
    struc = MolecularStructure(mmcif_dict)

    entity_ids = get(mmcif_dict, "_entity.id", String[]) # make entity dictionary
    entity_descs = get(mmcif_dict, "_entity.pdbx_description", String[])
    entity_types = get(mmcif_dict, "_entity.type", String[])

    entity_map = Dict(id => (desc, type) for (id, desc, type) in zip(entity_ids, entity_descs, entity_types))

    atom_chains = get(mmcif_dict, "_atom_site.label_asym_id", String[]) # entity to chain dictionary
    atom_entities = get(mmcif_dict, "_atom_site.label_entity_id", String[])

    unique_chain_mapping = Dict(zip(atom_chains, atom_entities))

    chain_annotations = Dict{String,String}()

    for (chain_id, ent_id) in unique_chain_mapping
        !haskey(entity_map, ent_id) && continue

        raw_desc, ent_type = entity_map[ent_id]
        desc_low = lowercase(raw_desc)

        category = if ent_type == "water"
            "Solvent"
        elseif ent_type == "non-polymer"
            "Ligand/Ion"
        elseif occursin("microglobulin", desc_low)
            "beta-2-microglobulin"
        elseif occursin(r"alpha|tcr a|receptor a", desc_low)
            "TCR-alpha"
        elseif occursin(r"beta|tcr b|receptor b", desc_low)
            "TCR-beta"
        elseif occursin(r"mhc|hla|class i", desc_low)
            "MHC"
        else
            "Antigen/Protein"
        end

        chain_annotations[chain_id] = "[$category] $raw_desc"
    end

    return struc, chain_annotations

end

function get_interactions!(db, pdb_id, struc, annotations)
    chains = collect(struc[1])

    for (i, chainA) in pairs(chains)
        idA = chainid(chainA)
        roleA = get(annotations, idA, nothing)
        if roleA === nothing || occursin("Solvent", roleA)
            continue
        end

        for j in (i+1):length(chains)
            chainB = chains[j]
            idB = chainid(chainB)
            roleB = get(annotations, idB, nothing)
            if roleB === nothing || occursin("Solvent", roleB)
                continue
            end

            for resA in chainA, resB in chainB
                check_pi_stacking!(db, pdb_id, resA, resB, idA, roleA, idB, roleB)

                for atomA in resA
                    isheavy(atomA) || continue
                    for atomB in resB
                        isheavy(atomB) || continue
                        check_atom_chemistry!(db, pdb_id, atomA, atomB, resA, resB, idA, roleA, idB, roleB)
                    end
                end
            end
        end
    end
end

### -------------------------
### Executing

pdb_ids = ["1AO7", "2BAN"]
interaction_db = DataFrame(
    PDB_ID=String[],
    Interaction_Type=String[],
    Distance=Float64[],
    Chain_1=String[],
    Chain_1_Class=String[],
    Resname_1=String[],
    Resnum_1=Int[],
    Chain_2=String[],
    Chain_2_Class=String[],
    Resname_2=String[],
    Resnum_2=Int[]
)

@showprogress "Processing PDBs: " for pdb in pdb_ids

    struc, annotations = get_annotations(pdb)

    get_interactions!(interaction_db, pdb, struc, annotations)

end
