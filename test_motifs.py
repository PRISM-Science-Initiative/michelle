from Bio import motifs
import os

folder = "motifs/jaspar_raw"

motif_library = {}

for file in os.listdir(folder):
    if file.endswith(".jaspar"):
        path = os.path.join(folder, file)

        with open(path) as f:
            m = motifs.read(f, "jaspar")
            motif_library[m.name] = m

print("Loaded motifs:", len(motif_library))