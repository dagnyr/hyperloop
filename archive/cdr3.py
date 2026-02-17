import requests

RCSB_GQL = "https://data.rcsb.org/graphql"

def rcsb_chain_map(pdb_ids):
    q = """
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
    }"""
    out = {}
    for pdb in pdb_ids:
        pdb = pdb.upper()
        j = requests.post(RCSB_GQL, json={"query": q, "variables": {"id": pdb}}, timeout=30).json()
        ents = (j.get("data", {}).get("entry", {}) or {}).get("polymer_entities", []) or []
        chains = []
        for e in ents:
            desc = (e.get("rcsb_polymer_entity", {}) or {}).get("pdbx_description", "")
            for inst in e.get("polymer_entity_instances", []) or []:
                ids = inst.get("rcsb_polymer_entity_instance_container_identifiers", {}) or {}
                chains.append({
                    "asym": ids.get("asym_id"),
                    "auth": ids.get("auth_asym_id"),
                    "desc": desc
                })
        out[pdb] = chains
    return out

# example
m = rcsb_chain_map(["9EJH", "1FYT"])
for pdb, chains in m.items():
    print(pdb)
    for c in chains:
        print(" ", c["auth"] or c["asym"], "->", c["desc"])
