from Bio import Entrez, SeqIO, Phylo
from Bio.Seq import Seq
import dendropy
import numpy as np
import pandas as pd
import json
from ete3 import Tree, TreeStyle

def setup_database(email: str = "akamipersona1111@gmail.com", id: str = "CP000944", phylo_tree: str = "MTBC_SNP_align_200.phy_phyml_tree.nexus"): # Data for M. tuberculosis H37Rv
    # Set email
    Entrez.email = email

    # Download antibiomatic genes
    handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    seq = Seq(str(record.seq[:1000]))

    # Features: One-hot encode
    def one_hot_encode(seq):
        mapping = {'A': [1,0,0,0], 'C': [0,1,0,0], 'G': [0,0,1,0], 'T': [0,0,0,1]}
        return np.array([mapping.get(base, [0,0,0,0]) for base in seq])

    features = one_hot_encode(seq)
    print("Gen features shape:", features.shape)  # (1000, 4)

    # Simulate evolutionary paths from a phylogenetic tree
    tree = dendropy.Tree.get(path=phylo_tree, schema=phylo_tree.split('.')[-1])
    paths = []
    for leaf in tree.leaf_node_iter():
        path = [node.taxon.label if node.taxon else "internal" for node in leaf.ancestor_iter(inclusive=True)]
        path.reverse()
        paths.append(path)
    print("Evolutionary paths:", len(paths))

    # Create DataFrame for gen_metadata
    gc_content = (seq.count('G') + seq.count('C')) / len(seq)

    df_gen = pd.DataFrame({
        'id': [1],
        'length': [len(seq)],
        'gc_content': [gc_content],
        'resistance_label': ["resistant"]
    })

    # Create DataFrame for evol_paths
    df_paths = pd.DataFrame({
        'path_id': range(1, len(paths) + 1),
        'nodes': [json.dumps(p) for p in paths],
        'probability': np.random.rand(len(paths))
    })

    df_gen.to_csv('gen_metadata.csv', index=False)
    df_paths.to_csv('evol_paths.csv', index=False)
    print("Exported to CSV for Access import.")

def visualize(file: str = "MTBC_SNP_align_200.phy_phyml_tree.nexus"):
    def build_tree(fname: str):
        """ build a tree from a file """
        tree = Phylo.read(fname, fname.split('.')[-1])
        return tree
    
    def build_tree_for_ete3(clade):
        """ Convert a biopython tree to an ete3 tree """
        builded_tree = Tree()
        builded_tree.name = clade.name if clade.name else ""
        
        for child in clade.clades:
            child_tree = build_tree_for_ete3(child)
            # Set branch length, default to 1.0 if None
            branch_len = child.branch_length if child.branch_length is not None else 3.0
            builded_tree.add_child(child=child_tree, name=child.name, dist=branch_len)
        
        return builded_tree
    
    if __name__ == '__main__':
        # Build tree from file
        tree = build_tree(file)
        
        # Optional: view with Phylo
        # Phylo.draw(tree)
        
        # Convert to ete3 format
        builded_tree = build_tree_for_ete3(tree.root)
        
        # Setup tree style
        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.mode = "c"  # circular mode
        ts.arc_start = -90  # 0 degrees = 3 o'clock
        ts.arc_span = 360
        ts.branch_vertical_margin = 50
        ts.scale = 800
        
   # Customize leaf name size and style
    from ete3 import TextFace, NodeStyle, AttrFace
    for node in builded_tree.traverse():
        # Tạo style cho node
        nstyle = NodeStyle()
        nstyle["size"] = 2  # Kích thước node (chấm tròn)
        nstyle["fgcolor"] = "#0066cc"  # Màu node (xanh dương)
        nstyle["hz_line_color"] = "#666666"  # Màu nhánh ngang
        nstyle["vt_line_color"] = "#666666"  # Màu nhánh dọc
        nstyle["hz_line_width"] = 2  # Độ dày nhánh ngang
        nstyle["vt_line_width"] = 2  # Độ dày nhánh dọc
        
        # Áp dụng style cho node
        node.set_style(nstyle)
    
    # Tùy chỉnh leaf names
    for leaf in builded_tree.iter_leaves():
        leaf.name = leaf.name if leaf.name else ""
        
        # Tạo text face với font nhỏ hơn
        name_face = TextFace(leaf.name, fsize=20, fgcolor="#ff6600")
        leaf.add_face(name_face, column=0, position="branch-right")
        
        # Style riêng cho leaf nodes
        leaf_style = NodeStyle()
        leaf_style["size"] = 6  # Leaf lớn hơn một chút
        leaf_style["fgcolor"] = "#ff6600"  # Màu cam cho leaf
        leaf_style["shape"] = "circle"  # Hình dạng: circle, square, sphere
        leaf_style["hz_line_color"] = "#ff6600"
        leaf_style["vt_line_color"] = "#37B2FF"
        leaf_style["hz_line_width"] = 2
        leaf_style["vt_line_width"] = 2
        leaf.set_style(leaf_style)
    
    # Display the tree
    builded_tree.show(tree_style=ts)
    output_file = "phylogenetic_tree.pdf"
    builded_tree.render(output_file, w=1200, h=1200, tree_style=ts, dpi=300)
    print(f"Tree saved to {output_file}")

setup_database()
visualize()
