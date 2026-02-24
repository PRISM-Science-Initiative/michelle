import os
os.environ['KMP_DUPLICATE_LIB_OK'] = 'True'
from driver import InferenceEngine

def process_folder(folder_name):
    if not os.path.exists(folder_name):
        print(f"Error: Folder '{folder_name}' not found.")
        return
    # set your API key in your environment or here
    engine = InferenceEngine(api_key=None)

    # filter for GenBank files
    gbk_files = [f for f in os.listdir(folder_name) if f.endswith(('.gbk', '.gb'))]

    for file_name in gbk_files:
        file_path = os.path.join(folder_name, file_name)
        print(f"=== PROCESSING: {file_name} ===")
        
        try:
            construct = engine.process_construct(file_path)

            # --- VALIDATION AUDIT ---
            print("\nPART VALIDATION REPORT (WEIGHTED SPLIT):")
            header = f"{'Part Name':<18} | {'Seq (15bp)':<15} | {'Phys (50%)':<10} | {'Sem (50%)':<10} | {'Total'}"
            print(header)
            print("-" * len(header))

            for part, _ in construct.parts_layout:
                # Get a short snippet of the sequence for visual verification
                seq_snippet = f"{part.sequence[:12]}..." if len(part.sequence) > 12 else part.sequence
                
                # Access the weighted contributions stored by the engine
                print(f"{part.name[:18]:<18} | "
                      f"{seq_snippet:<15} | "
                      f"{part.phys_contrib:<10.3f} | "
                      f"{part.sem_contrib:<10.3f} | "
                      f"{part.confidence_score}")
            
            avg_score = sum(p.confidence_score for p, _ in construct.parts_layout) / len(construct.parts_layout)
            print(f"\nConstruct Effectiveness (Avg Confidence): {avg_score:.2f}")
            print("\n" + "="*80 + "\n")
            
        except Exception as e:
            print(f"Error processing {file_name}: {e}")

if __name__ == "__main__":
    process_folder("genbank files")