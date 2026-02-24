from semantic_inferencer import SemanticInferencer
from llm_inferencer import LLMInferencer
from sequence_inferencer import SequenceInferencer
from genbank_to_python import parse_genbank
from felix_schema import Role, BioPart

class InferenceEngine:
    def __init__(self, api_key):
        self.semantic = SemanticInferencer()
        # self.llm = LLMInferencer(api_key)
        self.physical = SequenceInferencer()
        # Define deterministic weights
        self.WEIGHTS = {"physical": 0.50, "semantic": 0.50} # add the LLM part later

        self.SO_MAPPING = {
            Role.promoter: "SO:0000167",
            Role.cds: "SO:0000316",
            Role.rbs: "SO:0000139",
            Role.terminator: "SO:0000141",
            Role.recombinase_site: "SO:0000298"
    }

    def process_construct(self, file_path):
            construct = parse_genbank(file_path)
            for part, orientation in construct.parts_layout:
                # Tier 1: Semantic Role Guess
                clean_metadata = part.name # Just the label, e.g., "CMV promoter"
                role_obj, sem_conf = self.semantic.get_score(clean_metadata)
                # Tier 2: Physical Truth (Sequence Motifs)
                part.roles = [role_obj.name]
                phys_score = self.physical.validate_physical_truth(part)

                # Calculation of weighted contributions
                phys_contrib = phys_score * self.WEIGHTS["physical"]
                sem_contrib = sem_conf * self.WEIGHTS["semantic"]
                final_weight = phys_contrib + sem_contrib

                # Store the breakdown for transparency/debugging
                part.phys_score = phys_score
                part.sem_score = sem_conf
                part.phys_contrib = phys_contrib
                part.sem_contrib = sem_contrib
                
                # Finalize Part
                part.roles = [role_obj.name]
                part.confidence_score = round(final_weight, 3)
                # part.roles = [role_obj.name]
                # phys_score = self.physical.validate_physical_truth(part)
                # # Tier 3: LLM Refinement (Escalate if ambiguous)
                # # if sem_conf < 0.6:
                # #     llm_roles, llm_conf = self.llm.get_role_list(part.metadata + "\nSequence:" + part.sequence[:200])
                # #     # if this isn't enough of the sequence, we can try more
                # #     final_roles = llm_roles
                # # else:
                # final_roles = [role_obj.name]
                #     # llm_conf = sem_conf

                # # 3. Weighted Deterministic Scoring
                # final_weight = (
                #     (phys_score * self.WEIGHTS["physical"]) +
                #     (sem_conf * self.WEIGHTS["semantic"])
                # )

                # # Update the real part object with inferred data
                # part.roles = final_roles
                # part.confidence_score = round(final_weight, 3)
                
                # # map to SO Term for registry compliance
                # so_term = self.SO_MAPPING.get(role_obj, "SO:0000110") # Default: misc_feature
                # part.metadata += f" | SO_Term: {so_term}"

                if final_weight < 0.40:
                    # If we aren't sure, don't label it as a functional part
                    part.roles = ["unknown"]
                    part.metadata += " | INFERENCE_WARNING: Low Confidence"

            return construct