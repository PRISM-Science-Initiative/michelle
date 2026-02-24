import os
import hashlib
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
            self._load_motifs(motif_dir)
    
    # loads exact sequences from the tsv file
    def _load_hash_registry(self, path):
        registry = {}
        if not os.path.exists(path):
            return registry
        with open(path) as f:
            for line in f:
                seq, role = line.strip().split("\t")
                registry[seq.lower()] = Role[role]
        return registry

    # load motifs and convert to PSSM (logs instead math tingz)
    def _load_motifs(self, directory):
        for filename in os.listdir(directory):
            if not filename.endswith(".jaspar"):
                continue
            path = os.path.join(directory, filename)
            with open(path) as handle:
                for motif in motifs.parse(handle, "jaspar"):
                    pwm = motif.counts.normalize()
                    pssm = pwm.log_odds()
                    key = motif.name.lower()
                    self.motif_library[key] = pssm
    
    # search motif library for keyword matches & takes one with best matrix match
    def _best_motif_score(self, sequence, keyword):
        seq_obj = Seq(sequence)
        best = 0.0
        for name, pssm in self.motif_library.items():
            if keyword not in name:
                continue
            scores = pssm.calculate(seq_obj)
            if len(scores) == 0:
                continue
            score = (max(scores) - pssm.min) / (pssm.max - pssm.min)
            best = max(best, score)
        return round(float(best), 3)

    # actually does the validation check
    #    sequence of priority goes hash registry > CDS logic > promoter motif > hard coded terminator detection
    def validate_physical_truth(self, biopart):
        seq = biopart.sequence.lower().replace("\n", "").replace(" ", "")
        role_guess = biopart.primary_type 
            # I will have to change this so that its inferencing the thing tiself instead of just using the baseline guess
        
        # --- TIER 1: HASH REGISTRY (EXACT MATCH) ---
        #  if the sequence exists in our known library
        if seq in self.HASH_REGISTRY:
            if self.HASH_REGISTRY[seq].name == role_guess:
                return 1.0 # 100% confidence on exact match

        # --- TIER 2: LOGICAL RULES (CDS) ---
            # dr naik didn't like these so... whoops
        if role_guess == Role.cds.name:
            has_start = seq.startswith("atg")
            has_stop = any(seq.endswith(s) for s in ["taa", "tag", "tga"])
            is_triplet = len(seq) % 3 == 0
            if has_start and has_stop and is_triplet: return 1.0
            return 0.3

        # --- TIER 3: DYNAMIC MOTIF SCANNING (PROMOTERS/TERMINATORS) ---
        if role_guess == Role.promoter.name:
            tata_score = self._best_motif_score(seq, "tata")
            caat_score = self._best_motif_score(seq, "caat")
            gc_score = self._best_motif_score(seq, "gc")
            return max(tata_score, caat_score, gc_score, 0.2)

        if role_guess == Role.terminator.name:
            # canonical polyA signals
            if "aataaa" in seq or "attaaa" in seq:
                return 0.95
            poly = self._best_motif_score(seq, "poly")
            return max(poly, 0.2)

        return 0.5 # Default fallback