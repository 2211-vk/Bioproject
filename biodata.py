from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import dendropy
import numpy as np
import sqlite3
import pandas as pd
import json

# Set email
Entrez.email = "akamipersona1111@gmail.com"

# Download antibiomatic genes
handle = Entrez.efetch(db="nucleotide", id="CP000944", rettype="fasta", retmode="text")  # Public ID
record = SeqIO.read(handle, "fasta")
seq = Seq(str(record.seq[:1000]))  # Cắt ngắn cho demo

# Features: One-hot encode
def one_hot_encode(seq):
    mapping = {'A': [1,0,0,0], 'C': [0,1,0,0], 'G': [0,0,1,0], 'T': [0,0,0,1]}
    return np.array([mapping.get(base, [0,0,0,0]) for base in seq])

features = one_hot_encode(seq)
print("Gen features shape:", features.shape)  # (1000, 4)

# Simulate cây tiến hóa đơn giản
tree = dendropy.Tree.get(path="MTBC_SNP_align_200.phy_phyml_tree.nexus", schema="nexus")  # Tải public tree hoặc generate
paths = []
for leaf in tree.leaf_node_iter():
    path = [node.taxon.label if node.taxon else "internal" for node in leaf.ancestor_iter(inclusive=True)]
    path.reverse()  # từ root → lá
    paths.append(path)
print("Evolutionary paths:", len(paths))
print(paths)

# Tính gc_content
gc_content = (seq.count('G') + seq.count('C')) / len(seq)

# Tạo DataFrame cho gen_metadata (thay vì SQLite)
df_gen = pd.DataFrame({
    'id': [1],  # Ví dụ ID
    'length': [len(seq)],
    'gc_content': [gc_content],
    'resistance_label': ["resistant"]
})

# Tạo DataFrame cho evol_paths (nodes lưu JSON để dễ import)
df_paths = pd.DataFrame({
    'path_id': range(1, len(paths) + 1),
    'nodes': [json.dumps(p) for p in paths],  # Lưu list thành JSON string
    'probability': np.random.rand(len(paths))  # Giả lập probs
})

# Export sang CSV cho Access import
df_gen.to_csv('gen_metadata.csv', index=False)
df_paths.to_csv('evol_paths.csv', index=False)
print("Exported to CSV for Access import.")

