from genbank_to_python import parse_genbank
from felix_schema import Role, BioPart

mock_gb_file = """LOCUS       pMOCK_01               300 bp    DNA     circular SYN 
FEATURES             Location/Qualifiers
     source          1..300
     promoter        10..60
                     /label="J23100_Promoter"
     RBS             61..80
                     /label="B0034_RBS"
     CDS             81..180
                     /label="mCherry"
     terminator      181..220
                     /label="T1_Terminator"
     rep_origin      230..260
                     /label="pUC_ori"
     CDS             261..295
                     /label="AmpR"
ORIGIN      
        1 atcgatcgat ttgacagcta gctcagtcct aggtataatg ctagctacta gagaaagagg
       61 agaaatacta gatggtgagc aagggcgagg aggataacat ggccatcatc aaggagttca
      121 tgcgcttcaa ggtgcacatg gagggctccg tgaacggcca cgagttcgag atcgagggcg
      181 actaggcctg gataaaaaat gagctgttga caattaatca tcggctcgta taatgtgtgg
      241 aattgtgagc ggataacaat tttcacacag gaaacagcta tgaccatgat tacgccaagc
//"""

def create_parts_library():
    print("--- Scraping Plasmids into Library ---")
    construct = parse_genbank(mock_content=mock_gb_file)
    library = [part for part, orientation in construct.parts_layout]    
    print(f"Successfully loaded {len(library)} parts into the library:")
    for part in library:
        role = part.roles[0]
        role_name = role.name if hasattr(role, 'name') else role
        print(f" - {part.name} ({role_name}) | {len(part.sequence)} bp")
    return library

library = create_parts_library()


def assemble_naive_plasmid(library, desired_roles):
    print("\n--- Assembling Naive Plasmid ---")
    assembled_sequence = ""
    plasmid_map = []
    current_bp = 1
    
    for desired_role in desired_roles:
        # ROBUST MATCHING: Compare the string names instead of the Enum objects
        # This completely bypasses the cross-file import mismatch!
        part = next((p for p in library if any(desired_role.name == (r.name if hasattr(r, 'name') else r) for r in p.roles)), None)
        
        if not part:
            print(f"Warning: No part found for role {desired_role.name}. Skipping.")
            continue
            
        part_length = len(part.sequence)
        end_bp = current_bp + part_length - 1
        
        plasmid_map.append({
            "name": part.name,
            "role": desired_role.name,
            "start": current_bp,
            "end": end_bp,
            "length": part_length
        })
        
        assembled_sequence += part.sequence
        current_bp = end_bp + 1
        
    return {
        "map": plasmid_map,
        "sequence": assembled_sequence,
        "total_length": len(assembled_sequence)
    }

# Define our ideal plasmid architecture
architecture = [
    Role.promoter, 
    Role.rbs, 
    Role.protein_tag,  # Changed from Role.cds so mCherry gets used!
    Role.terminator, 
    Role.origin, 
    Role.selection_marker
]

new_plasmid = assemble_naive_plasmid(library, architecture)

# Print the final report
print(f"\nFINAL PLASMID REPORT (Total Length: {new_plasmid['total_length']} bp)")
print(f"{'PART NAME':<20} | {'ROLE':<18} | {'LOCATION'}")
print("-" * 60)
for feature in new_plasmid['map']:
    loc_str = f"{feature['start']} - {feature['end']}"
    print(f"{feature['name']:<20} | {feature['role']:<18} | {loc_str}")

print("\nFULL SEQUENCE:")
# Wrap sequence to 80 characters for readability
for i in range(0, len(new_plasmid['sequence']), 80):
    print(new_plasmid['sequence'][i:i+80])