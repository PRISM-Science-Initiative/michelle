import torch
import numpy as np
import faiss
from transformers import AutoTokenizer, AutoModel
from felix_schema import Role

# this is also chopped as hell and untested rn please don't use this,
#    it's just framework for what I will do later

class SemanticInferencer:
    def __init__(self):
        # Load BioBERT - specialized for biomedical text
        self.tokenizer = AutoTokenizer.from_pretrained("dmis-lab/biobert-v1.1-pubmed")
        self.model = AutoModel.from_pretrained("dmis-lab/biobert-v1.1-pubmed")
        
        # Ground Truth Reference Library
        self.references = [
            (Role.promoter, "DNA sequence that initiates transcription of a particular gene"),
            (Role.promoter, "CMV EF1a U6 CAG promoter high level expression"),
            (Role.terminator, "sequence that marks the end of a gene or operon in genomic DNA"),
            (Role.recombinase_site, "LoxP Lox2272 FRT site for site specific recombination"),
            (Role.structural, "Inverted Terminal Repeat ITR for viral packaging"),
            (Role.origin, "Origin of replication pBR322 ori f1 high copy number"),
            (Role.rbs, "Ribosome binding site Kozak sequence translation initiation"),
            (Role.verification, "PCR primer binding site sequencing verification probe"),
            (Role.cds, "Protein coding sequence gene open reading frame")
        ]
        self.roles, self.descriptions = zip(*self.references)
        self.index = self._build_index()

    def _get_embedding(self, text):
        inputs = self.tokenizer(text, return_tensors="pt", truncation=True, padding=True, max_length=128)
        with torch.no_grad():
            outputs = self.model(**inputs)
        return outputs.last_hidden_state.mean(dim=1).numpy()

    def _build_index(self):
        # vectorize all reference descriptions
        embeddings = np.vstack([self._get_embedding(d) for d in self.descriptions])
        # Initialize FAISS Index (L2 distance for similarity)
        index = faiss.IndexFlatL2(embeddings.shape[1])
        index.add(embeddings.astype('float32'))
        return index

    def get_score(self, text):
        query_vec = self._get_embedding(text).astype('float32')
        dist, indices = self.index.search(query_vec, k=1)
        confidence = 1.0 / (1.0 + dist[0][0])
        return self.roles[indices[0][0]], round(float(confidence), 2)