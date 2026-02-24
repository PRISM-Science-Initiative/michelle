from semantic_inferencer import SemanticInferencer
from llm_inferencer import LLMInferencer
from sequence_inferencer import SequenceInferencer
from genbank_to_python import parse_genbank
from felix_schema import Role, BioPart

class InferenceEngine:
    def __init__(self, api_key):
        self.semantic = SemanticInferencer()
        self.llm = LLMInferencer(api_key)
        self.physical = SequenceInferencer()
        # Define deterministic weights
        self.WEIGHTS = {"physical": 0.50, "semantic": 0.30, "llm": 0.20}

        SO_MAPPING = {
        Role.promoter: "SO:0000167",
        Role.cds: "SO:0000316",
        Role.rbs: "SO:0000139",
        Role.terminator: "SO:0000141",
        Role.recombinase_site: "SO:0000298"
    }

def process_construct(self, file_path):
        # 1. Extract non-dummy raw parts
        construct = parse_genbank(file_path)
        # 2. Iterate through real BioPart objects
        for part, orientation in construct.parts_layout:
            # Tier 1: Semantic Role Guess
            role_obj, sem_conf = self.semantic.get_score(part.metadata)
            
            # Tier 2: Physical Truth (Sequence Motifs)
            # Update part roles temporarily to validate the guess
            part.roles = [role_obj.name]
            phys_score = self.physical.validate_physical_truth(part)

            # Tier 3: LLM Refinement (Escalate if ambiguous)
            if sem_conf < 0.6:
                llm_roles, llm_conf = self.llm.get_role_list(part.metadata)
                final_roles = llm_roles
            else:
                final_roles = [role_obj.name]
                llm_conf = sem_conf

            # 3. Weighted Deterministic Scoring
            final_weight = (
                (phys_score * self.WEIGHTS["physical"]) +
                (sem_conf * self.WEIGHTS["semantic"]) +
                (llm_conf * self.WEIGHTS["llm"])
            )

            # Update the real part object with inferred data
            part.roles = final_roles
            part.confidence_score = round(final_weight, 2)
            
            # Audit Step: Map to SO Term for registry compliance
            so_term = self.SO_MAPPING.get(role_obj, "SO:0000110") # Default: misc_feature
            part.metadata += f" | SO_Term: {so_term}"

        return construct

# --- EXECUTION ---
if __name__ == "__main__":
    engine = InferenceEngine(api_key="sk-...")
    final_construct = engine.process_construct("addgene-plasmid-48138.gbk")
    
    # Export to Lean4 for formal verification
    print(final_construct.to_lean_definition())