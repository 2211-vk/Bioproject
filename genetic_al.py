"""
Workflow:
1. Load H37Rv reference genome + gene mapping
2. Align tree genomes (from phylogenetic tree) to H37Rv reference
3. Detect mutations at resistance-associated loci
4. Map mutations to genes → predict resistance phenotype
5. Output: genome_id | gene_name | mutation | drug | confidence
"""

import numpy as np
import pandas as pd
from deap import base, creator, tools, algorithms
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from ete3 import Tree
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_recall_fscore_support
import random
import logging
import json
import warnings
import re
from pathlib import Path
from typing import Dict, List
from collections import defaultdict

warnings.filterwarnings('ignore')
np.random.seed(42)
random.seed(42)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# =================== CONFIGURATION ===================
CONFIG = {
    'h37rv_fasta': 'h37rv_ref.fasta',
    'gene_mapping_file': 'gene_mapping_tb.json',
    'tree_file': 'MTBC_SNP_align_200.phy_phyml_tree.nexus',
    'metadata_file': 'gen_metadata.csv',
    'evol_paths_file': 'evol_paths.csv',
    'population_size': 100,
    'num_generations': 200,
    'crossover_prob': 0.7,
    'mutation_prob': 0.3,
    'tournament_size': 3,
    'test_size': 0.2
}
GENOME_CACHE = None
bounds = np.array([
    (0.001, 0.05),   # mutation_rate
    (0.0, 0.3),      # hgt_probability
    (0.5, 3.0),      # antibiotic_pressure
    (0.1, 2.0),      # efflux_pump_activity
    (0.0, 0.5),      # fitness_cost
    (0.0, 1.0),      # recombination_rate
    (0.0, 1.0),      # mobile_element_freq
    (0.0, 2.0),      # epistasis_coefficient
    (0.0, 1.0),      # population_bottleneck
    (0.0, 1.0),      # transmission_rate
    (0.0, 1.0),      # biofilm_formation
    (0.0, 2.0),      # stress_response
    (0.1, 2.0),       # threshold
    (0.0, 1.0),      # mutation_weight
    (0.0, 1.0)       # phylo_weight
])
dn_ds_results = None
# =================== DATA LOADER ===================
class DataLoader:
    """Load phylogenetic tree + metadata"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.tree = None
        self.tree_type = None
        self.metadata = None
        self.labels = None
        self.leaf_dict = {}
    
    def load_all(self):
        """Load all data"""
        logging.info("=== LOADING PHYLOGENETIC DATA ===")
        
        self._load_tree()
        self._load_metadata()
        
        return self
    
    def _load_tree(self):
        """Load phylogenetic tree"""
        try:
            try:
                self.tree = Tree(self.config['tree_file'], format=1)
                self.leaf_dict = {leaf.name: leaf for leaf in self.tree.get_leaves()}
                self.tree_type = 'ete3'
                logging.info(f"  ✓ Tree (ete3): {len(self.leaf_dict)} leaves")
            except:
                self.tree = Phylo.read(self.config['tree_file'], 'nexus')
                self.leaf_dict = {term.name: term for term in self.tree.get_terminals()}
                self.tree_type = 'phylo'
                logging.info(f"  ✓ Tree (Phylo): loaded")
        except Exception as e:
            logging.warning(f"  ⚠️ Load tree error: {e}")
            self._create_demo_tree()
            self.tree_type = 'ete3'
    
    def _create_demo_tree(self):
        """Create demo tree"""
        newick = "((A:0.1,B:0.2)AB:0.1,(C:0.15,D:0.25)CD:0.1)root;"
        self.tree = Tree(newick)
        self.leaf_dict = {leaf.name: leaf for leaf in self.tree.get_leaves()}
    
    def _load_metadata(self):
        """Load metadata"""
        try:
            self.metadata = pd.read_csv(self.config['metadata_file'])
            
            if 'resistance_label' in self.metadata.columns:
                self.labels = self.metadata['resistance_label'].map({
                    'resistant': 1, 'susceptible': 0
                }).values
            else:
                self.labels = np.random.randint(0, 2, len(self.metadata))
            
            if len(self.metadata) <= 1:
                self._expand_metadata_from_tree()
            
            logging.info(f"  ✓ Metadata: {len(self.metadata)} samples")
        except Exception as e:
            logging.warning(f"  ⚠️ Metadata error: {e}")
            self._create_metadata_from_tree()
    
    def _expand_metadata_from_tree(self):
        """Expand metadata from tree"""
        leaf_names = list(self.leaf_dict.keys())
        new_metadata = []
        for leaf_name in leaf_names:
            sample = {
                'genome_id': leaf_name,
                'resistance_label': np.random.choice(['resistant', 'susceptible'], p=[0.3, 0.7])
            }
            new_metadata.append(sample)
        
        self.metadata = pd.DataFrame(new_metadata)
        self.labels = self.metadata['resistance_label'].map({
            'resistant': 1, 'susceptible': 0
        }).values
        logging.info(f"  ✓ Expanded metadata: {len(self.metadata)} samples")
    
    def _create_metadata_from_tree(self):
        """Create metadata from tree"""
        leaf_names = list(self.leaf_dict.keys())
        metadata = []
        for leaf_name in leaf_names:
            sample = {
                'genome_id': leaf_name,
                'resistance_label': np.random.choice(['resistant', 'susceptible'], p=[0.3, 0.7])
            }
            metadata.append(sample)
        
        self.metadata = pd.DataFrame(metadata)
        self.labels = self.metadata['resistance_label'].map({
            'resistant': 1, 'susceptible': 0
        }).values
    
# =================== GENE REFERENCE LOADER ===================
class GeneReferenceLoader:
    """Load H37Rv reference genome + gene mapping"""
    
    def __init__(self, config: Dict, data_loader: DataLoader = None):
        self.config = config
        self.data_loader = data_loader
        self.h37rv_seq = None
        self.gene_mapping = None
        self.mutation_loci = {}
    
    def load_all(self):
        """Load all reference data"""
        logging.info("=== LOADING REFERENCE DATA ===")
        
        self._load_h37rv_fasta()
        self._load_gene_mapping()
        self._build_mutation_loci()
        
        logging.info(f"✓ Loaded H37Rv: {len(self.h37rv_seq)} bp")
        logging.info(f"✓ Loaded gene mapping: {len(self.mutation_loci)} resistance loci")
        
        return self
    
    def _load_h37rv_fasta(self):
        """Load H37Rv reference sequence"""
        try:
            for record in SeqIO.parse(self.config['h37rv_fasta'], 'fasta'):
                self.h37rv_seq = str(record.seq)
                logging.info(f"  ✓ H37Rv: {len(self.h37rv_seq)} bp")
                return
        except Exception as e:
            logging.error(f"  ✗ Error loading H37Rv: {e}")
    
    def _load_gene_mapping(self):
        """Load gene mapping from JSON"""
        try:
            with open(self.config['gene_mapping_file'], 'r') as f:
                data = json.load(f)
                self.gene_mapping = data['genes']
                logging.info(f"  ✓ Gene mapping: {len(self.gene_mapping)} genes")
        except Exception as e:
            logging.error(f"  ✗ Error loading gene mapping: {e}")
    
    def _build_mutation_loci(self):
        """Build index of mutation positions with proper genomic coordinate calculation"""
        
        for gene_id, gene_info in self.gene_mapping.items():
            start = gene_info['start']
            strand = gene_info.get('strand', '+')
            
            for mutation_info in gene_info['mutations']:
                mutation = mutation_info['mutation']
                confidence = mutation_info['confidence']
                
                try:
                    # Parse codon/position number from mutation string (e.g., S315T -> 315)
                    match = re.search(r'(\d+)', mutation)
                    if not match:
                        logging.debug(f"  Could not parse position from {mutation}")
                        continue
                    
                    pos = int(match.group(1))
                    
                    # Calculate genomic nucleotide position
                    if strand == '+':
                        codon_pos = start + (pos - 1) * 3
                    else:
                        # For reverse strand, adjust calculation
                        codon_pos = start - (pos - 1) * 3
                    
                    # Validate position
                    if codon_pos < 0 or codon_pos >= len(self.h37rv_seq):
                        logging.debug(f"  Position out of range for {mutation}: {codon_pos}")
                        continue
                    
                    if codon_pos not in self.mutation_loci:
                        self.mutation_loci[codon_pos] = []
                    
                    self.mutation_loci[codon_pos].append({
                        'gene_id': gene_id,
                        'gene_name': gene_info['gene_name'],
                        'mutation': mutation,
                        'drug': gene_info['drug'],
                        'confidence': confidence,
                        'position': codon_pos,
                        'strand': strand
                    })
                    
                except Exception as e:
                    logging.debug(f"  Error processing {mutation}: {e}")
                    continue
        
        logging.info(f"  ✓ Built mutation loci: {len(self.mutation_loci)} positions")
    
    def detect_mutations_in_sample(self, genome_id: str, sample_dir: str = 'sample_sequences') -> dict:
        """
        Detect resistance mutations in a sample by comparing to H37Rv reference.
        
        Parameters:
        - genome_id: sample identifier
        - sample_dir: directory containing per-sample FASTA files
        
        Returns:
        - detected_mutations: {gene_id: [mutation_str, ...]}
        """
        # detected = {}
        # sample_fasta = Path(sample_dir) / f"{genome_id}.fasta"
        
        # if not sample_fasta.exists():
        #     logging.debug(f"  Sample file not found: {sample_fasta}")
        #     return detected
        
        # try:
        #     # Read sample sequence
        #     for record in SeqIO.parse(str(sample_fasta), 'fasta'):
        #         sample_seq = str(record.seq)
        #         break
        # except Exception as e:
        #     logging.debug(f"  Error reading sample {genome_id}: {e}")
        #     return detected
        
        # # Initialize gene dict
        # for gene_id in self.gene_mapping.keys():
        #     detected[gene_id] = []
        GENOME_CACHE = self.load_sample_sequences()
        if genome_id in GENOME_CACHE:
            sample_seq = GENOME_CACHE[genome_id]
        else:
            return {} 

        detected = defaultdict(list)
        ref_seq = self.h37rv_seq
        
        # Check each known mutation locus
        for codon_pos, mutations in self.mutation_loci.items():
            if codon_pos >= len(sample_seq) or codon_pos >= len(ref_seq):
                continue
            
            # Check if position differs from reference
            if sample_seq[codon_pos] != ref_seq[codon_pos]:
                # Mutation is present
                for mut_info in mutations:
                    gene_id = mut_info['gene_id']
                    mutation_str = mut_info['mutation']
                    
                    if mutation_str not in detected[gene_id]:
                        detected[gene_id].append(mutation_str)
        
        return detected
    
    def load_sample_sequences(self):
        global GENOME_CACHE
        GENOME_CACHE = {}
        if not self.data_loader or self.data_loader.metadata is None:
            logging.warning("No metadata available for loading sequences")
            return GENOME_CACHE
        for idx, row in self.data_loader.metadata.iterrows():
            genome_id = row['genome_id']
            fasta_path = Path("sample_sequences") / f"{genome_id}.fasta"
            if fasta_path.exists():
                try:
                    record = next(SeqIO.parse(fasta_path, "fasta"))
                    GENOME_CACHE[genome_id] = str(record.seq)
                except:
                    GENOME_CACHE[genome_id] = self.h37rv_seq
            else:
                GENOME_CACHE[genome_id] = self.h37rv_seq
        return GENOME_CACHE

# =================== GA SETUP FOR GENE PREDICTION ===================
def setup_ga_toolbox():
    """Setup DEAP toolbox"""
    if hasattr(creator, "FitnessMax"):
        del creator.FitnessMax
    if hasattr(creator, "Individual"):
        del creator.Individual
    
    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)

def create_individual():
    """Create individual with 12 evolutionary parameters"""
    global bounds
    ind = [random.uniform(low, high) for low, high in bounds]
    return creator.Individual(ind)

# =================== GENE-LEVEL PREDICTOR ===================
class GeneLevelPredictor:
    """Predict resistance-causing mutations at gene level"""

    param_names = [
        "mutation_rate", "hgt_probability", "antibiotic_pressure", "efflux_pump",
        "fitness_cost", "recombination", "mobile_elements", "epistasis",
        "bottleneck", "transmission", "biofilm", "stress_response", 
        "threshold", "mutation_weight", "phylo_weight"
    ]
    
    def __init__(self, gene_loader: GeneReferenceLoader, data_loader: DataLoader):
        self.gene_loader = gene_loader
        self.data_loader = data_loader
        # Preload detected mutations for all samples
        self.detected_mutations_cache = {}
        self.convergent_mutations_cache = None
        self._precompute_convergent_mutations()

        if not hasattr(self, 'mutation_confidence_map'):
            self.mutation_confidence_map = {}
            for _, mut_infos in self.gene_loader.mutation_loci.items():
                for mut_info in mut_infos:
                    key = (mut_info['gene_id'], mut_info['mutation'])
                    self.mutation_confidence_map[key] = mut_info['confidence']

    
    def evaluate(self, individual):
        """
        Evaluate fitness by comparing observed mutation profiles to true labels.
        
        Fitness = accuracy of predicting resistance based on detected mutations
        vs. true resistance labels from metadata.
        """
        params = {name: val for name, val in zip(self.param_names, individual)}
        global dn_ds_results

        if not hasattr(self, 'dn_ds_cache'):
            self.dn_ds_cache = dn_ds_results or {}
        correct_predictions = 0
        total_predictions = 0
        
        for idx, row in self.data_loader.metadata.iterrows():
            if idx >= len(self.data_loader.labels):
                break
            
            genome_id = row.get('genome_id', f'sample_{idx}')
            true_label = self.data_loader.labels[idx]
            
            # Detect mutations in this sample (cache results)
            if genome_id not in self.detected_mutations_cache:
                self.detected_mutations_cache[genome_id] = self.gene_loader.detect_mutations_in_sample(genome_id)
            
            detected_muts = self.detected_mutations_cache[genome_id]
            
            # Calculate mutation score
            m_score = 0.0
            # convergent_muts = self.detect_convergent_mutations(self.data_loader.tree, self.detected_mutations_cache)
            convergent_muts = self.convergent_mutations_cache or {}
            for gene_id, mutations in detected_muts.items():
                # Get dN/dS ratio for this gene
                dn_ds_info = self.dn_ds_cache.get(gene_id, {'omega': 1.0})
                omega = dn_ds_info['omega']
                for mut_str in mutations:
                    key = (gene_id, mut_str)
                    base_confidence = self.mutation_confidence_map[(gene_id, mut_str)]
                    if key in convergent_muts:
                        base_confidence += 0.2*convergent_muts[key]['count']

                    # Apply GA parameters to modulate confidence
                    weight = (
                        1.0 +
                        params['antibiotic_pressure'] * 0.3 +
                        params['efflux_pump'] * 0.2 -
                        params['fitness_cost'] * 0.5
                    )
                    if omega > 1.0:
                        weight += min(0.3, (omega - 1.0) * 0.2)

                    m_score += base_confidence * weight
            
            # Calculate phylogenetic distance score
            leaf = self.data_loader.leaf_dict.get(genome_id)
            if leaf:
                p_score = self._calculate_evolutionary_distance(leaf, params)
            else:
                p_score = 0.0

            # Combine scores
            alpha = params.get('mutation_weight', 0.7)
            beta = params.get('phylo_weight', 0.3)
            total_confidence = m_score*alpha + p_score*beta

            # Predict resistant if sufficient evidence
            mutation_threshold = params.get('threshold', 0.5)  # At least one confident mutation
            predicted_resistance = 1 if total_confidence > mutation_threshold else 0
            
            if predicted_resistance == true_label:
                correct_predictions += 1
            
            total_predictions += 1
        
        accuracy = correct_predictions / total_predictions if total_predictions > 0 else 0.0
        return (accuracy,)
    
    def _calculate_evolutionary_distance(self, leaf, params: Dict) -> float:
        """Calculate evolutionary distance score"""
        score = 0.0
        node = leaf
        
        # Handle both ete3 Tree and Phylo Clade
        if hasattr(node, 'up'):  # ete3 Tree
            while node.up:
                dist = node.dist if hasattr(node, 'dist') and node.dist else 0.01
                mutations = params['mutation_rate'] * dist * 1000
                score += mutations
                
                if random.random() < params['hgt_probability']:
                    score += 5.0 * params['antibiotic_pressure']
                
                node = node.up
        else:  # Phylo Clade - use branch length
            branch_length = leaf.branch_length if hasattr(leaf, 'branch_length') else 0.01
            mutations = params['mutation_rate'] * branch_length * 1000
            score += mutations
            
            if random.random() < params['hgt_probability']:
                score += 5.0 * params['antibiotic_pressure']
        
        resistance_score = (
            score * params['antibiotic_pressure'] +
            params['efflux_pump'] * 2.0 +
            params['mobile_elements'] * 1.5
        ) / (1 + params['fitness_cost'])
        
        return resistance_score
    
    def count_synonymous_sites(self, codon: str):
        """
        Count number of synonymous sites in a codon.
        
        A synonymous site is a position where a substitution would not
        change the amino acid (silent mutation).
        
        Parameters:
        -----------
        codon : str
            3-nucleotide codon sequence (e.g., "ATG")
        
        Returns:
        --------
        float : number of synonymous sites (0-3)
        
        Example:
        --------
        Codon "ATG" (Methionine):
        - Position 1 (A): Change to C,G,T → all non-synonymous → 0 synonymous
        - Position 2 (T): Change to A,C,G → all non-synonymous → 0 synonymous  
        - Position 3 (G): Change to A → ATG→ATA (Met→Ile) non-syn
                          Change to C → ATG→ATC (Met→Ile) non-syn
                          Change to T → ATG→ATT (Met→Ile) non-syn
                          → 0 synonymous
        Total: 0 synonymous sites
        
        Codon "GCT" (Alanine):
        - Position 3 can change to A,C,G → GCA,GCC,GCG all code Alanine
        → ~1 synonymous site at position 3
        """
        if len(codon) != 3:
            return 0.0
        
        try:
            ref_aa = Seq(codon).translate(to_stop=False)
        except:
            return 0.0
        
        synonymous_count = 0.0
        bases = ['A', 'T', 'G', 'C']
        
        # Check each position in codon
        for pos in range(3):
            syn_changes = 0
            total_changes = 0
            
            original_base = codon[pos]
            
            # Try changing to each other base
            for base in bases:
                if base == original_base:
                    continue
                
                # Create mutated codon
                mutated_codon = codon[:pos] + base + codon[pos+1:]
                
                try:
                    mutated_aa = Seq(mutated_codon).translate(to_stop=False)
                    total_changes += 1
                    
                    # Check if amino acid unchanged (synonymous)
                    if mutated_aa == ref_aa:
                        syn_changes += 1
                except:
                    continue
            
            # Proportion of synonymous changes at this position
            if total_changes > 0:
                synonymous_count += syn_changes / total_changes
        
        return synonymous_count
    
    def calculate_dn_ds_for_gene(self, gene_id, sample_sequences, reference_seq):
        """
        Calculate dN/dS ratio for a gene
        dN/dS > 1 → positive selection
        """

        # Get gene region
        gene_info = self.gene_loader.gene_mapping[gene_id]
        start, end = gene_info['start'], gene_info['end']

        # Extract gene sequences
        ref_gene_seq = reference_seq[start:end]
        sample_gene_seqs = [seq[start:end] for seq in sample_sequences]

        # Count synonymous and non-synonymous substitutions
        synonymous_subs = 0
        non_synonymous_subs = 0
        synonymous_sites = 0
        non_synonymous_sites = 0

        for i in range(0, len(ref_gene_seq), 3):
            ref_codon = ref_gene_seq[i:i+3]
            if len(ref_codon) != 3:
                continue
            
            try:
                ref_aa = Seq(ref_codon).translate()
            except:
                continue
            
            # Count sites
            synonymous_sites += self.count_synonymous_sites(ref_codon)
            non_synonymous_sites += 3 - self.count_synonymous_sites(ref_codon)

            # Count substitutions across all samples
            for sample_seq in sample_gene_seqs:
                sample_codon = sample_seq[i:i+3]
                if len(sample_codon) != 3:
                    continue
                
                try:
                    sample_aa = Seq(sample_codon).translate()
                except:
                    continue
                
                if ref_codon != sample_codon:
                    if ref_aa == sample_aa:
                        synonymous_subs += 1
                    else:
                        non_synonymous_subs += 1

        # Calculate dN/dS
        dN = non_synonymous_subs / non_synonymous_sites if non_synonymous_sites > 0 else 0
        dS = synonymous_subs / synonymous_sites if synonymous_sites > 0 else 1e-10
        print(synonymous_subs, synonymous_sites)
        omega = np.divide(dN, dS, out=np.zeros_like(dN), where=dS!=0)
        return {
            'omega': omega,
            'dN': dN,
            'dS': dS,
            'non_syn_subs': non_synonymous_subs,
            'syn_subs': synonymous_subs,
            'non_syn_sites': non_synonymous_sites,
            'syn_sites': synonymous_sites,
            'selection_type': (
                'positive' if omega > 1.0 else
                'neutral' if omega == 1.0 else
                'purifying'
            )
        }
    
    def calculate_dn_ds_for_all_genes(self) -> Dict[str, Dict]:
        """
        Calculate dN/dS for all resistance genes across all samples.
        
        Returns:
        --------
        dict : {gene_id: {omega, dN, dS, ...}}
        """
        logging.info("\n=== CALCULATING dN/dS FOR ALL GENES ===")
        
        # Collect all sample sequences
        sample_sequences = []
        for idx, row in self.data_loader.metadata.iterrows():
            genome_id = row.get('genome_id', f'sample_{idx}')
            sample_fasta = Path('sample_sequences') / f"{genome_id}.fasta"
            if sample_fasta.exists():
                try:
                    for record in SeqIO.parse(str(sample_fasta), 'fasta'):
                        sample_sequences.append(str(record.seq))
                        break
                except Exception as e:
                    logging.debug(f"Error reading {genome_id}: {e}")
        
        if not sample_sequences:
            logging.warning("No sample sequences found for dN/dS calculation")
            return {}
        
        # Calculate dN/dS for each gene
        dn_ds_results = {}
        for gene_id in self.gene_loader.gene_mapping.keys():
            result = self.calculate_dn_ds_for_gene(
                gene_id,
                sample_sequences,
                self.gene_loader.h37rv_seq
            )
            dn_ds_results[gene_id] = result
            
            # Log interesting results
            if result['omega'] > 1.5:
                logging.info(f"  🔴 {gene_id}: omega={result['omega']:.2f} (STRONG positive selection)")
            elif result['omega'] > 1.0:
                logging.info(f"  🟡 {gene_id}: omega={result['omega']:.2f} (positive selection)")
        
        return dn_ds_results
    
    def count_synonymous_sites(self, codon: str) -> float:
        """
        Count number of synonymous sites in a codon.
        
        A synonymous site is a position where a substitution would not
        change the amino acid (silent mutation).
        
        Parameters:
        -----------
        codon : str
            3-nucleotide codon sequence (e.g., "ATG")
        
        Returns:
        --------
        float : number of synonymous sites (0-3)
        
        Example:
        --------
        Codon "ATG" (Methionine):
        - Position 1 (A): Change to C,G,T → all non-synonymous → 0 synonymous
        - Position 2 (T): Change to A,C,G → all non-synonymous → 0 synonymous  
        - Position 3 (G): Change to A → ATG→ATA (Met→Ile) non-syn
                          Change to C → ATG→ATC (Met→Ile) non-syn
                          Change to T → ATG→ATT (Met→Ile) non-syn
                          → 0 synonymous
        Total: 0 synonymous sites
        
        Codon "GCT" (Alanine):
        - Position 3 can change to A,C,G → GCA,GCC,GCG all code Alanine
        → ~1 synonymous site at position 3
        """
        if len(codon) != 3:
            return 0.0
        
        try:
            ref_aa = Seq(codon).translate(to_stop=False)
        except:
            return 0.0
        
        synonymous_count = 0.0
        bases = ['A', 'T', 'G', 'C']
        
        # Check each position in codon
        for pos in range(3):
            syn_changes = 0
            total_changes = 0
            
            original_base = codon[pos]
            
            # Try changing to each other base
            for base in bases:
                if base == original_base:
                    continue
                
                # Create mutated codon
                mutated_codon = codon[:pos] + base + codon[pos+1:]
                
                try:
                    mutated_aa = Seq(mutated_codon).translate(to_stop=False)
                    total_changes += 1
                    
                    # Check if amino acid unchanged (synonymous)
                    if mutated_aa == ref_aa:
                        syn_changes += 1
                except:
                    continue
            
            # Proportion of synonymous changes at this position
            if total_changes > 0:
                synonymous_count += syn_changes / total_changes
        
        return synonymous_count
    
    def export_dn_ds_results(self, dn_ds_results: Dict, output_file: str = "dn_ds_analysis.csv"):
        """
        Export dN/dS analysis results to CSV.

        Parameters:
        -----------
        dn_ds_results : dict
            Results from calculate_dn_ds_for_all_genes()
        output_file : str
            Output CSV filename
        """

        records = []
        for gene_id, result in dn_ds_results.items():
            records.append({
                'gene_id': gene_id,
                'omega': result['omega'],
                'dN': result['dN'],
                'dS': result['dS'],
                'non_syn_subs': result['non_syn_subs'],
                'syn_subs': result['syn_subs'],
                'non_syn_sites': result['non_syn_sites'],
                'syn_sites': result['syn_sites'],
                'selection_type': result['selection_type']
            })

        df = pd.DataFrame(records)
        df = df.sort_values('omega', ascending=False)
        df.to_csv(output_file, index=False)

        logging.info(f"✓ Saved dN/dS analysis: {output_file}")

        # Print summary
        print("\n" + "="*80)
        print("dN/dS ANALYSIS SUMMARY")
        print("="*80)

        print(f"\nTotal genes analyzed: {len(df)}")
        print(f"Genes under positive selection (omega > 1): {sum(df['omega'] > 1)}")
        print(f"Genes under purifying selection (omega < 1): {sum(df['omega'] < 1)}")
        print(f"Genes under neutral evolution (omega ≈ 1): {sum((df['omega'] >= 0.9) & (df['omega'] <= 1.1))}")

        print(f"\nTop 10 genes under POSITIVE selection:")
        print(df.head(10)[['gene_id', 'omega', 'selection_type']].to_string(index=False))

        print(f"\nTop 10 genes under PURIFYING selection:")
        print(df.tail(10)[['gene_id', 'omega', 'selection_type']].to_string(index=False))

        print("\n" + "="*80)

        return df
    
    def _precompute_convergent_mutations(self):
        """Pre-compute convergent mutations across all samples"""
        logging.info("Pre-computing convergent mutations...")

        if not self.detected_mutations_cache:
            for idx, row in self.data_loader.metadata.iterrows():
                genome_id = row.get('genome_id', f'sample_{idx}')
                if genome_id not in self.detected_mutations_cache:
                    self.detected_mutations_cache[genome_id] = \
                        self.gene_loader.detect_mutations_in_sample(genome_id)

        # Calculate convergent mutations
        self.convergent_mutations_cache = self.detect_convergent_mutations(
            self.data_loader.tree,
            self.detected_mutations_cache
        )

        count = len(self.convergent_mutations_cache)
        logging.info(f"Finished! Found {count} convergent mutations")

    def detect_convergent_mutations(self, tree, detected_mutations_all_samples):
        """
        Find mutations that evolved independently in multiple lineages
        → Signal of strong selection pressure (convergent evolution)
        """

        # Group mutations by type
        mutation_to_samples = defaultdict(list)
        for sample_id, muts in detected_mutations_all_samples.items():
            for gene_id, mutations in muts.items():
                for mut_str in mutations:
                    mutation_to_samples[(gene_id, mut_str)].append(sample_id)

        # Find mutations in multiple distant lineages
        convergent_mutations = {}
        for (gene_id, mut_str), sample_ids in mutation_to_samples.items():
            if len(sample_ids) < 2:
                continue
            
            # Check if samples are from different lineages (distant in tree)
            if self.data_loader.tree_type == 'ete3':
                leaves = [tree.search_nodes(name=sid)[0] for sid in sample_ids if tree.search_nodes(name=sid)]
            else:
                leaves = [list(tree.find_clades(name=sid))[0] for sid in sample_ids if list(tree.find_clades(name=sid))]
            if len(leaves) < 2:
                continue
            
            # Calculate pairwise distances
            min_distance = float('inf')
            for i in range(len(leaves)):
                for j in range(i+1, len(leaves)):
                    if self.data_loader.tree_type == 'ete3':
                        dist = leaves[i].get_distance(leaves[j])
                    else:
                        dist = tree.distance(leaves[i], leaves[j])
                    min_distance = min(min_distance, dist)

            # If samples are distant → convergent
            if min_distance > 0.1:  # Threshold for "distant"
                convergent_mutations[(gene_id, mut_str)] = {
                    'count': len(sample_ids),
                    'min_distance': min_distance,
                    'samples': sample_ids
                }

        return convergent_mutations

    def predict_genes_for_genome(self, genome_id: str, params: Dict) -> List[Dict]:
        """
        Predict resistance-causing genes for a genome.
        
        Only includes mutations that are actually detected in the sample.
        Uses mutation confidence and GA parameters to compute probability.
        """
        predictions = []
        
        # Get detected mutations for this sample
        if genome_id not in self.detected_mutations_cache:
            self.detected_mutations_cache[genome_id] = self.gene_loader.detect_mutations_in_sample(genome_id)
        
        detected_muts = self.detected_mutations_cache[genome_id]
        
        # Only predict for detected mutations
        for gene_id, mutations in detected_muts.items():
            if not mutations:  # No mutations detected in this gene
                continue
            
            gene_info = self.gene_loader.gene_mapping.get(gene_id)
            if not gene_info:
                continue
            
            for mutation_str in mutations:
                confidence = 0.0
                for mut_info in gene_info['mutations']:
                    if mut_info['mutation'] == mutation_str:
                        confidence = mut_info['confidence']
                        break
                
                prob = (
                    params['mutation_rate'] * 0.1 +
                    params['antibiotic_pressure'] * 0.2 +
                    (1 - params['fitness_cost']) * confidence * 0.7
                )
                
                prob = min(1.0, max(0.0, prob))
                
                predictions.append({
                    'genome_id': genome_id,
                    'gene_id': gene_id,
                    'gene_name': gene_info['gene_name'],
                    'mutation': mutation_str,
                    'drug': gene_info['drug'],
                    'confidence': confidence,
                    'probability': prob
                })
        
        return predictions

# =================== MAIN GA EXECUTION ===================
def run_genetic_al_gene_level():
    """Run GA for gene-level resistance prediction"""
    
    logging.info("\n" + "="*80)
    logging.info("GENETIC ALGORITHM v2 - GENE-LEVEL ANTIBIOTIC RESISTANCE PREDICTION")
    logging.info("="*80)
    
    data_loader = DataLoader(CONFIG).load_all()
    gene_loader = GeneReferenceLoader(CONFIG, data_loader).load_all()
    
    setup_ga_toolbox()
    predictor = GeneLevelPredictor(gene_loader, data_loader)
    predictor.detected_mutations_cache = {}
    for genome_id in data_loader.metadata['genome_id']:
        predictor.detected_mutations_cache[genome_id] = predictor.gene_loader.detect_mutations_in_sample(genome_id)
    
    # ===== CALCULATE dN/dS BEFORE GA =====
    logging.info("\n=== PRE-GA: dN/dS ANALYSIS ===")
    dn_ds_results = predictor.calculate_dn_ds_for_all_genes()
    df_dn_ds = predictor.export_dn_ds_results(dn_ds_results, "dn_ds_analysis.csv")
    
    n_samples = len(data_loader.labels)
    if n_samples > 10:
        try:
            train_idx, test_idx = train_test_split(
                range(n_samples),
                test_size=CONFIG['test_size'],
                random_state=42,
                stratify=data_loader.labels
            )
            logging.info(f"✓ Train/test split: {len(train_idx)} / {len(test_idx)}")
        except ValueError:
            train_idx = test_idx = list(range(n_samples))
            logging.warning("⚠️ Using full dataset for train/test")
    else:
        train_idx = test_idx = list(range(n_samples))
    
    toolbox = base.Toolbox()
    toolbox.register("individual", create_individual)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("mate", tools.cxBlend, alpha=0.5)
    toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=0.1, indpb=0.2)
    toolbox.register("select", tools.selTournament, tournsize=CONFIG['tournament_size'])
    toolbox.register("evaluate", predictor.evaluate)
    
    min_vals = bounds[:,0]
    max_vals = bounds[:,1]
    
    def check_bounds(min_v, max_v):
        def decorator(func):
            def wrapper(*args, **kwargs):
                offspring = func(*args, **kwargs)
                for child in offspring:
                    child[:] = np.clip(child, min_v, max_v)
                return offspring
            return wrapper
        return decorator
    
    toolbox.decorate("mate", check_bounds(min_vals, max_vals))
    toolbox.decorate("mutate", check_bounds(min_vals, max_vals))
    
    population = toolbox.population(n=CONFIG['population_size'])
    hof = tools.HallOfFame(10)
    
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("max", np.max)
    
    logging.info(f"\n=== RUNNING GA ({CONFIG['num_generations']} generations) ===")
    population, logbook = algorithms.eaSimple(
        population, toolbox,
        cxpb=CONFIG['crossover_prob'],
        mutpb=CONFIG['mutation_prob'],
        ngen=CONFIG['num_generations'],
        stats=stats,
        halloffame=hof,
        verbose=True
    )
    
    export_gene_predictions(hof, logbook, predictor, data_loader, gene_loader)
    
    return hof, logbook, df_dn_ds

# =================== EXPORT GENE-LEVEL PREDICTIONS ===================
def export_gene_predictions(hof, logbook, predictor, data_loader, gene_loader):
    """Export gene-level resistance predictions"""
    
    logging.info("\n=== EXPORTING GENE-LEVEL PREDICTIONS ===")
    
    best = hof[0]
    best_params = {name: float(best[i]) for i, name in enumerate(predictor.param_names)}
    
    logging.info("\n=== DETECTED MUTATIONS PER SAMPLE ===")
    for idx, row in data_loader.metadata.iterrows():
        genome_id = row.get('genome_id', f'sample_{idx}')
        true_label = 'resistant' if data_loader.labels[idx] == 1 else 'susceptible'
        
        detected = predictor.detected_mutations_cache.get(genome_id, {})
        num_mutations = sum(len(muts) for muts in detected.values())
        
        logging.info(f"  {genome_id}: {num_mutations} mutations detected (true: {true_label})")
        for gene_id, muts in detected.items():
            if muts:
                logging.info(f"    - {gene_id}: {', '.join(muts)}")
    
    all_predictions = []
    for idx, row in data_loader.metadata.iterrows():
        genome_id = row.get('genome_id', f'sample_{idx}')
        gene_preds = predictor.predict_genes_for_genome(genome_id, best_params)
        
        for pred in gene_preds:
            all_predictions.append({
                'genome_id': genome_id,
                'true_label': 'resistant' if data_loader.labels[idx] == 1 else 'susceptible',
                'gene_id': pred['gene_id'],
                'gene_name': pred['gene_name'],
                'mutation': pred['mutation'],
                'drug': pred['drug'],
                'mutation_confidence': f"{pred.get('confidence', 0.0):.4f}",
                'resistance_probability': f"{pred['probability']:.4f}"
            })
    
    df_predictions = pd.DataFrame(all_predictions)
    df_predictions.to_csv("genetic_al_gene_predictions.csv", index=False)
    logging.info(f"\n✓ Saved genetic_al_gene_predictions.csv ({len(df_predictions)} predictions)")
    
    with open("genetic_al_best_gene_params.json", "w") as f:
        json.dump(best_params, f, indent=2)
    logging.info("✓ Saved genetic_al_best_gene_params.json")
    
    df_stats = pd.DataFrame(logbook)
    df_stats.to_csv("genetic_al_gene_evolution_stats.csv", index=False)
    logging.info("✓ Saved genetic_al_gene_evolution_stats.csv")
    
    print_gene_prediction_summary(best_params, df_predictions, logbook, data_loader)

def print_gene_prediction_summary(params: Dict, predictions_df: pd.DataFrame, logbook, data_loader=None):
    """Print summary of gene-level predictions"""
    
    print("\n" + "="*80)
    print("GENE-LEVEL RESISTANCE PREDICTION RESULTS")
    print("="*80)
    
    print(f"\nPredictions Summary:")
    print(f"  Total predictions: {len(predictions_df)}")
    print(f"  Unique genomes: {predictions_df['genome_id'].nunique()}")
    print(f"  Resistance genes detected: {predictions_df['gene_id'].nunique()}")
    print(f"  Drug categories: {predictions_df['drug'].nunique()}")
    
    print(f"\nTop Resistance Genes:")
    if len(predictions_df) > 0:
        gene_counts = predictions_df['gene_name'].value_counts().head(5)
        for gene, count in gene_counts.items():
            print(f"  * {gene}: {count} times")
    
    print(f"\nDrugs Associated with Detected Mutations:")
    if len(predictions_df) > 0:
        drug_counts = predictions_df['drug'].value_counts()
        for drug, count in drug_counts.items():
            print(f"  * {drug}: {count} mutations")
    
    print(f"\nModel Performance:")
    print(f"  Initial avg fitness: {logbook[0]['avg']:.4f}")
    print(f"  Final avg fitness: {logbook[-1]['avg']:.4f}")
    print(f"  Best fitness achieved: {logbook[-1]['max']:.4f}")
    
    print(f"\nOutput Files:")
    print(f"  * genetic_al_gene_predictions.csv")
    print(f"  * genetic_al_best_gene_params.json")
    print(f"  * genetic_al_gene_evolution_stats.csv")
    
    print("\n" + "="*80)
