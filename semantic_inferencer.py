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
        self.tokenizer = AutoTokenizer.from_pretrained("monologg/biobert_v1.1_pubmed")
        self.model = AutoModel.from_pretrained("monologg/biobert_v1.1_pubmed")
        
        # Ground Truth Reference Library
        # semantic_inferencer.py
        # semantic_inferencer.py
        self.references = [
            (Role.promoter, "promoter CMV U6 EF1a CAG T7 SV40 transcription start"),
            (Role.terminator, "terminator PolyA signal BGH hGH transcription stop"),
            (Role.cds, "CDS coding sequence gene ORF EGFP mCherry Cre Cas9"),
            (Role.origin, "origin replication ori pBR322 f1 colE1 high copy"),
            (Role.verification, "primer sequencing forward reverse bind"),
            (Role.selection_marker, "selection marker resistance bla AmpR KanR NeoR"),
            (Role.recombinase_site, "LoxP FRT recombinase site Cre-mediated"),
            (Role.structural, "ITR LTR viral structural packaging signal")
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
        if dist[0][0] > 100.0: 
            return Role.selection_marker, 0.05
        confidence = 1.0 / (1.0 + dist[0][0])
        return self.roles[indices[0][0]], round(float(confidence), 2)