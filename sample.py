import os
from genbank_to_python import parse_genbank

def process_folder(folder_name):
    if not os.path.exists(folder_name):
        print(f"Error: Folder '{folder_name}' not found.")
        return

    # filter for GenBank files
    gbk_files = [f for f in os.listdir(folder_name) if f.endswith(('.gbk', '.gb'))]

    for file_name in gbk_files:
        file_path = os.path.join(folder_name, file_name)
        print(f"=== PROCESSING: {file_name} ===")
        
        try:
            plasmid = parse_genbank(file_path)
            
            # output the Lean4 code block
            print(plasmid.to_lean_definition())
            
            # count parts for summary
            print(f"Summary: Extracted {len(plasmid.parts_layout)} parts from {file_name}.")
            print("\n" + "="*60 + "\n")
            
        except Exception as e:
            print(f"Error processing {file_name}: {e}")

if __name__ == "__main__":
    process_folder("genbank files")