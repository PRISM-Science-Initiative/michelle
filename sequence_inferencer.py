import os
from Bio import motifs
from Bio.Seq import Seq
from felix_schema import Role

class SequenceInferencer:    
    def __init__(self, motif_dir="motifs/curated",
                 hash_file="data/hash_registry.tsv"):
        # Tier 1: Exact Hash Registry (Ground Truth Sequences)
        self.HASH_REGISTRY = self._load_hash_registry(hash_file) 
            # edit the hash file if you need to

        # Tier 2: Motif Library
        self.motif_library = {}
        if os.path.exists(motif_dir):
            self._load_motifs_recursive(motif_dir)
    
    # loads exact sequences from the tsv file
    def _load_hash_registry(self, path):
        registry = {}
        if not os.path.exists(path): return registry
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("sequence"): continue
                parts = line.split("\t")
                if len(parts) >= 2:
                    # STRIP EVERYTHING: spaces, newlines, tabs
                    clean_seq = parts[0].lower().replace(" ", "").strip()
                    registry[clean_seq] = Role[parts[1].strip()]
        return registry

    # load motifs and convert to PSSM (logs instead math tingz)
    def _load_motifs_recursive(self, directory):
        for root, dirs, files in os.walk(directory):
            for filename in files:
                if filename.endswith(".jaspar"):
                    path = os.path.join(root, filename)
                    with open(path) as handle:
                        for motif in motifs.parse(handle, "jaspar"):
                            # Normalize counts to a Log-Odds Position Specific Scoring Matrix
                            pssm = motif.counts.normalize().log_odds()
                            # Use both folder name and motif name for indexing
                            folder_name = os.path.basename(root).lower()
                            key = f"{folder_name}_{motif.name.lower()}"
                            self.motif_library[key] = pssm
    
    # search motif library for keyword matches & takes one with best matrix match
    def _best_motif_score(self, sequence, keyword):
        seq_obj = Seq(sequence)
        best_raw_score = -100.0 # Start extremely low
        target_pssm = None

        if not self.motif_library: return 0.5
        for name, pssm in self.motif_library.items():
            if keyword in name:
                scores = pssm.calculate(seq_obj)
                if len(scores) == 0: continue
                if max(scores) > best_raw_score:
                    best_raw_score = max(scores)
                    target_pssm = pssm

        if target_pssm is None:
            return 0.5
                    
        denom = target_pssm.max - target_pssm.min
        if denom == 0:
            # If the motif has no variance, any match is a 1.0
            return 1.0 if best_raw_score >= target_pssm.min else 0.0

        normalized = (best_raw_score - target_pssm.min) / denom
        return round(float(max(0.0, normalized)), 3)

    # actually does the validation check
    #    sequence of priority goes hash registry > CDS logic > promoter motif > hard coded terminator detection
    def validate_physical_truth(self, biopart):
        seq = biopart.sequence.lower().replace("\n", "").replace(" ", "")
        role_guess = biopart.primary_type 
        
        # --- TIER 1: HASH REGISTRY (EXACT MATCH) ---
        #  if the sequence exists in our known library
        if seq in self.HASH_REGISTRY:
            if self.HASH_REGISTRY[seq].name == role_guess:
                return 1.0 # 100% confidence on exact match
            
        # Logic for Origins and Markers (The "Meaningful" Rule)
        # These don't have motifs, so we verify them by length/composition
        if role_guess in [Role.origin.name, Role.selection_marker.name]:
            if len(seq) > 500: return 0.9  # High confidence for large structural genes
            if len(seq) > 100: return 0.7
            return 0.4

        # --- TIER 2: LOGICAL RULES (CDS) ---
        if role_guess == Role.cds.name:
            score = 0.0
            # 1. Length evidence
            if len(seq) < 50: return 0.01
            if len(seq) > 300: score += 0.4
            # 2. Start codons (not just ATG)
            if seq[:3] in ["atg", "gtg", "ttg"]: score += 0.15
            # 3. Stop codon near end (within last 6bp)
            if any(stop in seq[-6:] for stop in ["taa", "tag", "tga"]):
                score += 0.15
            # 4. Internal stop penalty
            first_80 = seq[:int(len(seq)*0.8)]
            if any(s in first_80 for s in ["taa", "tag", "tga"]):
                score -= 0.4 # Significant penalty
                
            return round(max(0.01, min(score, 1.0)), 3)

        # --- TIER 3: DYNAMIC MOTIF SCANNING (PROMOTERS/TERMINATORS) ---
        if role_guess == Role.promoter.name:
            score = 0.0
            tata = self._best_motif_score(seq, "tata")
            caat = self._best_motif_score(seq, "caat")
            gc = self._best_motif_score(seq, "gc")
            motif_score = max(tata, caat, gc)
            score += 0.5 * motif_score

            # AT richness (common in promoters)
            at_content = (seq.count("a") + seq.count("t")) / len(seq)
            if at_content > 0.6:
                score += 0.2
            if len(seq) < 500:
                score += 0.1  # promoters are usually shorter

            return round(min(score, 1.0), 3)

        if role_guess == Role.terminator.name:
            score = 0.0
            # PolyA signal
            if "aataaa" in seq or "attaaa" in seq:
                score += 0.5
            # Poly-T stretch (prokaryotic)
            if "tttt" in seq[-20:]:
                score += 0.3
            # GC hairpin proxy (simple heuristic)
            gc_content = (seq.count("g") + seq.count("c")) / len(seq)
            if gc_content > 0.5:
                score += 0.1
            poly = self._best_motif_score(seq, "poly")
            score += 0.3 * poly

            return round(min(score, 1.0), 3)
        motif_categories = {
            Role.rbs.name: ["rbs"],
            Role.insulator.name: ["insulator"],
            Role.non_coding_rna.name: ["non_coding_rna"],
            Role.structural.name: ["structural"],
            Role.enhancer.name: ["enhancer"]
        }

        if role_guess in motif_categories:
            best_cat_score = 0.0
            for keyword in motif_categories[role_guess]:
                best_cat_score = max(best_cat_score, self._best_motif_score(seq, keyword))
            return best_cat_score


        return 0.5 # Default fallback