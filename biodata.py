"""
BIODATA.PY - Complete biological data preparation pipeline

Functions:
1. Fetch genome from NCBI (with user input: email, id, phylo_tree)
2. Prepare gene mapping (H37Rv reference + TB-Profiler resistance mutations)
3. Visualize phylogenetic tree with resistance coloring and branch lengths
4. Generate output files for GA-based resistance prediction

Default values:
- email: akamipersona1111@gmail.com
- genome_id: NC_000962.3 (M. tuberculosis H37Rv)
- phylo_tree: MTBC_SNP_align_200.phy_phyml_tree.nexus
"""

from Bio import Entrez, SeqIO, Phylo
from Bio.Seq import Seq
import dendropy
import numpy as np
import pandas as pd
import json
import logging
import re
from pathlib import Path
from ete3 import Tree, TreeStyle, NodeStyle, TextFace
import sys

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

builded_tree = None
genomic_features = None
tree = None

# =================== TB RESISTANCE MUTATIONS DATABASE ===================
KNOWN_TB_RESISTANCE_MUTATIONS = {
    'katG': {
        'gene_name': 'catalase-peroxidase KatG',
        'drug': 'INH (Isoniazid)',
        'mutations': [
            {'mutation': 'S315T', 'effect': 'Loss of KatG activation', 'confidence': 0.95},
            {'mutation': 'S315N', 'effect': 'Reduced KatG activation', 'confidence': 0.85},
            {'mutation': 'S315A', 'effect': 'Reduced KatG activation', 'confidence': 0.80},
            {'mutation': 'R463L', 'effect': 'Loss of catalase activity', 'confidence': 0.75},
        ]
    },
    'inhA': {
        'gene_name': 'enoyl-CoA reductase InhA',
        'drug': 'INH, ETH (Ethionamide)',
        'mutations': [
            {'mutation': 'I21V', 'effect': 'INH resistance', 'confidence': 0.90},
            {'mutation': 'I21T', 'effect': 'INH resistance', 'confidence': 0.88},
            {'mutation': 'S94A', 'effect': 'INH/ETH resistance', 'confidence': 0.85},
            {'mutation': 'C61S', 'effect': 'INH resistance', 'confidence': 0.80},
        ]
    },
    'fabG1': {
        'gene_name': '3-ketoacyl-CoA thiolase',
        'drug': 'INH',
        'mutations': [
            {'mutation': 'C203S', 'effect': 'Reduced INH activation', 'confidence': 0.75},
        ]
    },
    'rpoB': {
        'gene_name': 'RNA polymerase beta subunit',
        'drug': 'RIF (Rifampicin), RFB',
        'mutations': [
            {'mutation': 'H445L', 'effect': 'High-level RIF resistance', 'confidence': 0.99},
            {'mutation': 'H445Y', 'effect': 'High-level RIF resistance', 'confidence': 0.99},
            {'mutation': 'H445D', 'effect': 'High-level RIF resistance', 'confidence': 0.99},
            {'mutation': 'S450F', 'effect': 'High-level RIF resistance', 'confidence': 0.98},
            {'mutation': 'S450L', 'effect': 'High-level RIF resistance', 'confidence': 0.98},
            {'mutation': 'D435V', 'effect': 'RIF resistance', 'confidence': 0.95},
            {'mutation': 'L430P', 'effect': 'RIF resistance', 'confidence': 0.92},
            {'mutation': 'Q432K', 'effect': 'RIF resistance', 'confidence': 0.85},
        ]
    },
    'gyrA': {
        'gene_name': 'DNA gyrase subunit A',
        'drug': 'FQ (Fluoroquinolone)',
        'mutations': [
            {'mutation': 'D94G', 'effect': 'High-level FQ resistance', 'confidence': 0.98},
            {'mutation': 'D94A', 'effect': 'High-level FQ resistance', 'confidence': 0.98},
            {'mutation': 'D94H', 'effect': 'High-level FQ resistance', 'confidence': 0.95},
            {'mutation': 'A94V', 'effect': 'FQ resistance', 'confidence': 0.90},
            {'mutation': 'S95T', 'effect': 'FQ resistance', 'confidence': 0.95},
            {'mutation': 'S95F', 'effect': 'FQ resistance', 'confidence': 0.92},
        ]
    },
    'gyrB': {
        'gene_name': 'DNA gyrase subunit B',
        'drug': 'FQ (Fluoroquinolone)',
        'mutations': [
            {'mutation': 'E21Q', 'effect': 'Low-level FQ resistance', 'confidence': 0.70},
            {'mutation': 'A74S', 'effect': 'Low-level FQ resistance', 'confidence': 0.70},
        ]
    },
    'etaC': {
        'gene_name': 'transcriptional regulator (ETH resistance)',
        'drug': 'ETH (Ethionamide)',
        'mutations': [
            {'mutation': 'G1A', 'effect': 'Loss of ETH activation', 'confidence': 0.85},
            {'mutation': '-10C->T', 'effect': 'Promoter mutation', 'confidence': 0.75},
        ]
    },
    'embA': {
        'gene_name': 'arabinosyl transferase EmbA',
        'drug': 'EMB (Ethambutol)',
        'mutations': [
            {'mutation': 'M1V', 'effect': 'EMB resistance', 'confidence': 0.85},
            {'mutation': 'G406D', 'effect': 'EMB resistance', 'confidence': 0.80},
        ]
    },
    'embB': {
        'gene_name': 'arabinosyl transferase EmbB',
        'drug': 'EMB (Ethambutol)',
        'mutations': [
            {'mutation': 'M306V', 'effect': 'High-level EMB resistance', 'confidence': 0.95},
            {'mutation': 'M306I', 'effect': 'High-level EMB resistance', 'confidence': 0.92},
            {'mutation': 'M306L', 'effect': 'High-level EMB resistance', 'confidence': 0.88},
            {'mutation': 'D1024A', 'effect': 'EMB resistance', 'confidence': 0.80},
        ]
    },
    'embC': {
        'gene_name': 'arabinosyl transferase EmbC',
        'drug': 'EMB (Ethambutol)',
        'mutations': [
            {'mutation': 'G1A', 'effect': 'EMB resistance', 'confidence': 0.78},
        ]
    },
    'pncA': {
        'gene_name': 'pyrazinamidase',
        'drug': 'PZA (Pyrazinamide)',
        'mutations': [
            {'mutation': 'H51R', 'effect': 'Loss of PZA activation', 'confidence': 0.92},
            {'mutation': 'T134P', 'effect': 'Loss of PZA activation', 'confidence': 0.88},
            {'mutation': 'H96R', 'effect': 'Loss of PZA activation', 'confidence': 0.85},
            {'mutation': 'K131E', 'effect': 'Reduced PZA activation', 'confidence': 0.80},
        ]
    },
    'rpsA': {
        'gene_name': '30S ribosomal protein S1',
        'drug': 'PZA (alternative target)',
        'mutations': [
            {'mutation': 'H51R', 'effect': 'PZA resistance', 'confidence': 0.70},
        ]
    },
    'rpsL': {
        'gene_name': '30S ribosomal protein S12',
        'drug': 'SM (Streptomycin)',
        'mutations': [
            {'mutation': 'K43R', 'effect': 'High-level SM resistance', 'confidence': 0.98},
            {'mutation': 'K43N', 'effect': 'High-level SM resistance', 'confidence': 0.95},
            {'mutation': 'K43E', 'effect': 'High-level SM resistance', 'confidence': 0.92},
        ]
    },
    'rrs': {
        'gene_name': '16S rRNA',
        'drug': 'SM, AMI (Amikacin), KAN (Kanamycin)',
        'mutations': [
            {'mutation': 'A1401G', 'effect': 'SM/AMI/KAN resistance', 'confidence': 0.98},
            {'mutation': 'G1484T', 'effect': 'SM resistance', 'confidence': 0.90},
            {'mutation': 'C1402T', 'effect': 'SM resistance', 'confidence': 0.85},
        ]
    },
    'rplC': {
        'gene_name': '50S ribosomal protein L3',
        'drug': 'LZD (Linezolid)',
        'mutations': [
            {'mutation': 'A154V', 'effect': 'Linezolid resistance', 'confidence': 0.85},
        ]
    },
    'rplD': {
        'gene_name': '50S ribosomal protein L4',
        'drug': 'LZD (Linezolid)',
        'mutations': [
            {'mutation': 'A194V', 'effect': 'Linezolid resistance', 'confidence': 0.80},
        ]
    },
}

# =================== GENOME FETCHING ===================
def fetch_genome_from_ncbi(email=None, genome_id=None):
    """
    Fetch genome sequence from NCBI with user input.
    
    Parameters:
    - email: NCBI Entrez email (default: prompt user)
    - genome_id: NCBI accession ID (default: prompt user or NC_000962.3)
    
    Returns: (record, full_seq, seq_1000bp)
    """
    
    if email is None:
        user_email = input("\nPlease provide your NCBI email (or Enter for default): ").strip()
        email = user_email if user_email else "akamipersona1111@gmail.com"
    
    if genome_id is None:
        user_id = input("Enter NCBI Genome ID (or Enter for NC_000962.3 - H37Rv): ").strip()
        genome_id = user_id if user_id else "NC_000962.3"
    
    logging.info(f"Loading genome {genome_id} from NCBI...")
    Entrez.email = email
    
    try:
        handle = Entrez.efetch(db="nuccore", id=genome_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        full_seq = str(record.seq)
        seq_1000bp = Seq(full_seq[:1000])
        
        logging.info(f"LOADING SUCCESS!")
        logging.info(f"  - ID: {record.id}")
        logging.info(f"  - Lenght: {len(full_seq)} bp")
        
        return record, full_seq, seq_1000bp
    except Exception as e:
        logging.error(f"ERROR: {e}")
        raise

# =================== GENE MAPPING PREPARATION ===================
def create_gene_mapping_from_manual_data():
    """Create gene mapping from H37Rv (NC_000962.3) annotation"""
    
    logging.info("create gene mapping from H37Rv annotation...")
    
    gene_coordinates = {
        'gyrA': {'start': 7380, 'end': 9818, 'strand': '+'},
        'gyrB': {'start': 9818, 'end': 11488, 'strand': '+'},
        'rpoB': {'start': 761155, 'end': 765422, 'strand': '+'},
        'katG': {'start': 1772200, 'end': 1775315, 'strand': '+'},
        'inhA': {'start': 1673425, 'end': 1674327, 'strand': '+'},
        'fabG1': {'start': 1674531, 'end': 1675343, 'strand': '+'},
        'embA': {'start': 3636303, 'end': 3637838, 'strand': '+'},
        'embB': {'start': 3637838, 'end': 3639378, 'strand': '+'},
        'embC': {'start': 3639378, 'end': 3640999, 'strand': '+'},
        'pncA': {'start': 2288681, 'end': 2289180, 'strand': '+'},
        'rpsA': {'start': 3871903, 'end': 3872796, 'strand': '+'},
        'rpsL': {'start': 1473920, 'end': 1474411, 'strand': '+'},
        'rrs': {'start': 1473197, 'end': 1474651, 'strand': '+'},
        'rplC': {'start': 1474644, 'end': 1476170, 'strand': '+'},
        'rplD': {'start': 1476206, 'end': 1477390, 'strand': '+'},
        'etaC': {'start': 920361, 'end': 920942, 'strand': '+'},
    }
    
    gene_mapping = []
    
    for gene_id, gene_info in KNOWN_TB_RESISTANCE_MUTATIONS.items():
        coord = gene_coordinates.get(gene_id, {'start': 0, 'end': 0, 'strand': '+'})
        
        for mutation in gene_info['mutations']:
            gene_mapping.append({
                'gene_id': gene_id,
                'gene_name': gene_info['gene_name'],
                'start': coord['start'],
                'end': coord['end'],
                'strand': coord['strand'],
                'mutation': mutation['mutation'],
                'effect': mutation['effect'],
                'drug': gene_info['drug'],
                'confidence': mutation['confidence']
            })
    
    df_mapping = pd.DataFrame(gene_mapping)
    logging.info(f"Create: {len(df_mapping)} mutations, {len(df_mapping['gene_id'].unique())} genes")
    
    return df_mapping

def save_gene_mapping_json(df_mapping, output_file='gene_mapping_tb.json'):
    """Save gene mapping as JSON"""
    
    genes_dict = {}
    
    for gene_id in df_mapping['gene_id'].unique():
        gene_df = df_mapping[df_mapping['gene_id'] == gene_id]
        first_row = gene_df.iloc[0]
        
        genes_dict[gene_id] = {
            'gene_name': first_row['gene_name'],
            'start': int(first_row['start']),
            'end': int(first_row['end']),
            'strand': first_row['strand'],
            'drug': first_row['drug'],
            'mutations': [
                {
                    'mutation': row['mutation'],
                    'effect': row['effect'],
                    'confidence': float(row['confidence'])
                }
                for _, row in gene_df.iterrows()
            ]
        }
    
    output = {
        'organism': 'Mycobacterium tuberculosis',
        'strain': 'H37Rv',
        'reference': 'NC_000962.3',
        'genes': genes_dict
    }
    
    with open(output_file, 'w') as f:
        json.dump(output, f, indent=2)
    
    logging.info(f"Save {output_file}")
    return output

def save_gene_mapping_csv(df_mapping, output_file='gene_mapping_tb.csv'):
    """Save gene mapping as CSV"""
    df_mapping.to_csv(output_file, index=False)
    logging.info(f"Save {output_file}")

# =================== GENOMIC DATA PREPARATION ===================
def prepare_genome_data(record, full_seq, seq_1000bp):
    """Prepare genomic features and save FASTA"""
    global genomic_features
    
    logging.info("Preparing genetic data...")
    
    # Save full FASTA
    output_fasta = 'h37rv_ref.fasta'
    with open(output_fasta, 'w') as f:
        f.write(f">{record.id}\n{full_seq}\n")
    logging.info(f"Save {output_fasta}")
    
    # One-hot encode first 1000bp
    def one_hot_encode(seq):
        mapping = {'A': [1,0,0,0], 'C': [0,1,0,0], 'G': [0,0,1,0], 'T': [0,0,0,1]}
        return np.array([mapping.get(base, [0,0,0,0]) for base in seq])
    
    features = one_hot_encode(seq_1000bp)
    genomic_features = features
    np.save('genomic_features.npy', features)
    logging.info(f"Save genomic_features.npy (shape: {features.shape})")
    
    # GC content
    gc_content = (seq_1000bp.count('G') + seq_1000bp.count('C')) / len(seq_1000bp)
    
    # Create metadata
    df_gen = pd.DataFrame({
        'genome_id': [record.id],
        'length': [len(full_seq)],
        'gc_content': [gc_content],
        'resistance_label': ['unknown']
    })
    
    df_gen.to_csv('gen_metadata.csv', index=False)
    logging.info(f"Save gen_metadata.csv")

def prepare_evolutionary_paths(phylo_tree):
    """Prepare evolutionary paths from phylogenetic tree"""
    logging.info("Preparing evolutionary paths...")
    global tree
    
    tree = dendropy.Tree.get(path=phylo_tree, schema=phylo_tree.split('.')[-1])
    paths = []
    for leaf in tree.leaf_node_iter():
        path = [node.taxon.label if node.taxon else "internal" for node in leaf.ancestor_iter(inclusive=True)]
        path.reverse()
        paths.append(path)
    
    df_paths = pd.DataFrame({
        'path_id': range(1, len(paths) + 1),
        'nodes': [json.dumps(p) for p in paths],
        'probability': np.random.rand(len(paths))
    })
    
    df_paths.to_csv('evol_paths.csv', index=False)
    logging.info(f"Save evol_paths.csv ({len(paths)} paths)")

# =================== TREE VISUALIZATION ===================
def visualize(file: str = "MTBC_SNP_align_200.phy_phyml_tree.nexus"):
    """Visualize phylogenetic tree with resistance coloring and branch lengths"""
    global builded_tree
    
    def build_tree(fname: str):
        """Build tree from file"""
        tree = Phylo.read(fname, fname.split('.')[-1])
        return tree
    
    def build_tree_for_ete3(clade):
        """Convert BioPython tree to ete3 tree"""
        builded_tree = Tree()
        builded_tree.name = clade.name if clade.name else ""
        
        for child in clade.clades:
            child_tree = build_tree_for_ete3(child)
            branch_len = child.branch_length if child.branch_length is not None else 0.01
            builded_tree.add_child(child=child_tree, name=child.name, dist=branch_len)
        
        return builded_tree
    
    if builded_tree is None:
        logging.info(f"Build tree from {file}...")
        tree_bio = build_tree(file)
        builded_tree = build_tree_for_ete3(tree_bio.root)
    
    # Setup tree style
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.mode = "c"
    ts.arc_start = -90
    ts.arc_span = 360
    ts.branch_vertical_margin = 50
    ts.scale = 800
    
    # Customize node styles
    leaf_color = "#ff6600"
    line_color = "#ff6600"
    for node in builded_tree.traverse():
        nstyle = NodeStyle()
        
        if node.is_leaf():
            # All leaf nodes use orange color
            nstyle["size"] = 8
            nstyle["fgcolor"] = leaf_color
            nstyle["shape"] = "circle"
            nstyle["hz_line_color"] = line_color
            nstyle["vt_line_color"] = line_color
            nstyle["hz_line_width"] = 2
            nstyle["vt_line_width"] = 2
        else:
            nstyle["size"] = 2
            nstyle["fgcolor"] = "#0066cc"
            nstyle["hz_line_color"] = "#888888"
            nstyle["vt_line_color"] = "#888888"
            nstyle["hz_line_width"] = 1
            nstyle["vt_line_width"] = 1

        node.set_style(nstyle)

        # Add leaf labels
        if node.is_leaf():
            leaf_name = node.name if node.name else ""
            lf = TextFace(leaf_name, fsize=15, fgcolor="#333333")
            node.add_face(lf, column=0, position="branch-right")
        
        # Add branch lengths
        if not node.is_leaf() and node.up:
            branch_len = node.dist if hasattr(node, 'dist') and node.dist else 0
            if branch_len > 0:
                dist_text = f"{branch_len:.3f}"
                df = TextFace(dist_text, fsize=8, fgcolor="#666666")
                node.add_face(df, column=1, position="branch-top")
        
        if node.is_leaf() and node.up:
            branch_len = node.dist if hasattr(node, 'dist') and node.dist else 0
            if branch_len > 0:
                dist_text = f"{branch_len:.4f}"
                df = TextFace(dist_text, fsize=12, fgcolor="#999999")
                node.add_face(df, column=1, position="branch-top")

    # Tweak TreeStyle
    ts.show_leaf_name = False
    ts.mode = "c"
    ts.arc_start = -90
    ts.arc_span = 360
    ts.branch_vertical_margin = 40
    ts.scale = 1200

    # Export PDF
    output_pdf = "phylogenetic_tree.pdf"
    builded_tree.render(output_pdf, w=2000, h=2000, tree_style=ts, dpi=300)
    logging.info(f"Save tree: {output_pdf}")

# =================== PER-SAMPLE SEQUENCE GENERATION & MUTATION DETECTION ===================
def generate_sample_sequences_with_mutations(full_seq: str, phylo_tree: str, df_mapping: pd.DataFrame, output_dir: str = "sample_sequences", mutation_data: str = 'mutations.csv'):
    """
    Generate per-sample genome sequences with simulated resistance mutations.
    
    Uses phylogenetic tree distances to:
    1. Create per-sample sequences by introducing mutations at resistance loci.
    2. Assign resistance labels based on mutation profile (resistant if >=1 mutation, else susceptible).
    3. Save per-sample FASTA files.
    4. Return updated metadata with true resistance labels.
    """
    
    logging.info("\nCreate sequences for each sample with mutations...")
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)

    mutation_df = pd.read_csv(mutation_data)
    R_mutations = mutation_df[mutation_df['confidence'].str.contains('Assoc w R', na=False)]['Mutation'].str.upper().tolist()
    KNOWN_RESISTANCE_GENES = [
        'KATG', 'INHA', 'RPOB', 'RPOC', 'GYRA', 'GYRB', 'PNCA', 'RPSL', 'RRS', 'EIS',
        'EMBB', 'ETHR', 'TLYA', 'RRB', 'FABG1'
    ]

    # Calculate leaf distances from root
    LEAF_DISTANCES = {}
    max_distance = 0.0
    
    for leaf in tree.leaf_node_iter():
        leaf_name = leaf.taxon.label.strip("'\"")
        distance = leaf.distance_from_root()
        LEAF_DISTANCES[leaf_name] = distance
    
        if distance > max_distance:
            max_distance = distance
            
    sample_names = list(LEAF_DISTANCES.keys())

    # Build mutation locus map: {gene_id: [(nuc_pos, mutation_str, confidence), ...]}
    mutation_loci = {}
    for gene_id in df_mapping['gene_id'].unique():
        gene_df = df_mapping[df_mapping['gene_id'] == gene_id]
        first_row = gene_df.iloc[0]
        start = first_row['start']
        strand = first_row['strand']
        
        mutation_loci[gene_id.upper()] = []
        
        for _, row in gene_df.iterrows():
            mut_str = row['mutation'].upper()
            confidence = row.get('confidence', 0.8)
            
            # Parse codon number from mutation string (e.g., S315T -> 315)
            match = re.search(r'(\d+)', mut_str)
            if match:
                codon_num = int(match.group(1))
                # Calculate genomic nucleotide position
                if strand == '+':
                    nuc_pos = start + (codon_num - 1) * 3
                else:
                    nuc_pos = start - (codon_num - 1) * 3
                
                # Validate position is within sequence
                if 0 <= nuc_pos < len(full_seq):
                    mutation_loci[gene_id.upper()].append((nuc_pos, mut_str, confidence))
    
    # Generate per-sample sequences
    sample_metadata = []
    np.random.seed(42)
    
    for sample_id in sample_names:
        sample_seq = list(full_seq)  # Copy reference sequence
        detected_mutations = []
        r_mutations = []
        
        # For each gene, probabilistically introduce mutations
        evo_dist = LEAF_DISTANCES.get(sample_id, 0.001)
        evo_factor = min(evo_dist / max_distance * 2.5, 1.0)
        for gene_id, mutations in mutation_loci.items():
            if gene_id not in KNOWN_RESISTANCE_GENES:
                continue
            for nuc_pos, mut_str, confidence in mutations:
                # Probability of mutation presence increases with confidence
                # and a small random component (simulating evolutionary history)
                prob_have_mutation = confidence * evo_factor * 0.8 + np.random.rand() * 0.05
                prob_have_mutation = min(prob_have_mutation, 0.90)
                
                if np.random.rand() < prob_have_mutation:
                    # Introduce mutation (simple SNP at position)
                    if nuc_pos < len(sample_seq):
                        bases = ['A', 'C', 'G', 'T']
                        sample_seq[nuc_pos] = np.random.choice([b for b in bases if b != sample_seq[nuc_pos]])
                        detected_mutations.append({'gene_id': gene_id, 'mutation': mut_str, 'position': nuc_pos})
                        if mut_str in R_mutations:
                            r_mutations.append(mut_str)
        # Determine resistance label: resistant if has any mutation
        if r_mutations:
            resistance_label = "resistant"
        else:
            resistance_prob = evo_factor * 0.85
            resistance_label = "resistant" if np.random.rand() < resistance_prob else "susceptible"
        # Save sample FASTA
        sample_fasta = Path(output_dir) / f"{sample_id}.fasta"
        with open(sample_fasta, 'w') as f:
            f.write(f">{sample_id}\n")
            f.write(''.join(sample_seq))
            f.write('\n')
        
        # Add to metadata
        sample_metadata.append({
            'genome_id': sample_id,
            'length': len(full_seq),
            'gc_content': (full_seq.count('G') + full_seq.count('C')) / len(full_seq),
            'resistance_label': resistance_label,
            'num_mutations': len(detected_mutations),
            'mutations': json.dumps([m['mutation'] for m in detected_mutations]),
            'evo_distance': round(evo_dist, 6)
        })
    
    # Save updated metadata
    df_meta = pd.DataFrame(sample_metadata)
    df_meta.to_csv('gen_metadata.csv', index=False)
    
    logging.info(f"Create {len(sample_names)} samples: {sum(df_meta['resistance_label']=='resistant')} resistant, {sum(df_meta['resistance_label']=='susceptible')} susceptible")
    
    return df_meta

def detect_mutations_in_sample(sample_fasta: str, ref_seq: str, gene_mapping: dict) -> dict:
    """
    Detect resistance mutations in a single sample by comparing to reference.
    
    Parameters:
    - sample_fasta: path to sample FASTA file
    - ref_seq: reference sequence (H37Rv)
    - gene_mapping: gene mapping dict with loci
    
    Returns:
    - detected_mutations: dict {gene_id: [mutation_str, ...]}
    """
    
    detected = {}
    
    try:
        # Read sample sequence
        record = SeqIO.read(sample_fasta, "fasta")
        sample_seq = str(record.seq)
    except:
        return detected
    
    # Compare at known resistance loci
    for gene_id, gene_info in gene_mapping.items():
        detected[gene_id] = []
        start = gene_info.get('start', 0)
        
        for mut_info in gene_info['mutations']:
            mut_str = mut_info['mutation']
            
            # Parse position from mutation string
            match = re.search(r'(\d+)', mut_str)
            if match:
                codon_num = int(match.group(1))
                nuc_pos = start + (codon_num - 1) * 3
                
                if 0 <= nuc_pos < len(sample_seq) and nuc_pos < len(ref_seq):
                    # Check if the position differs
                    if sample_seq[nuc_pos] != ref_seq[nuc_pos]:
                        detected[gene_id].append(mut_str)
    
    return detected

# =================== MAIN PIPELINE ===================
def setup_database(email: str = None, genome_id: str = None, phylo_tree: str = None):
    """
    Complete biological data preparation pipeline
    
    Parameters:
    - email: NCBI email (default: prompt user)
    - genome_id: NCBI genome ID (default: prompt user or NC_000962.3)
    - phylo_tree: Path to phylogenetic tree (default: MTBC_SNP_align_200.phy_phyml_tree.nexus)
    """
    
    print("\n" + "="*80)
    print("PREPARING BIOLOGICAL DATA FOR TB RESISTANCE PREDICTION")
    print("="*80)
    
    # If no phylo_tree provided, use default
    if phylo_tree is None:
        phylo_tree = "MTBC_SNP_align_200.phy_phyml_tree.nexus"
    
    # Fetch genome from NCBI
    record, full_seq, seq_1000bp = fetch_genome_from_ncbi(email, genome_id)
    
    # Prepare genomic data
    prepare_genome_data(record, full_seq, seq_1000bp)
    
    # Prepare evolutionary paths
    prepare_evolutionary_paths(phylo_tree)
    
    # Prepare gene mapping
    logging.info("\n--- Gene Mapping Preparation ---")
    df_mapping = create_gene_mapping_from_manual_data()
    save_gene_mapping_json(df_mapping)
    save_gene_mapping_csv(df_mapping)
    
    # Generate per-sample sequences with mutations
    logging.info("\n--- Per-Sample Sequence Generation ---")
    df_metadata = generate_sample_sequences_with_mutations(full_seq, phylo_tree, df_mapping)
    
    print("\n" + "="*80)
    print("BIOLOGICAL DATA PREPARATION COMPLETED!")
    print("="*80)
    print("\nFiles tao:")
    print("  + h37rv_ref.fasta")
    print("  + gen_metadata.csv")
    print("  + evol_paths.csv")
    print("  + genomic_features.npy")
    print("  + gene_mapping_tb.json")
    print("  + gene_mapping_tb.csv")
    print("="*80 + "\n")

if __name__ == "__main__":
    setup_database()
    visualize()
