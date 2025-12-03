"""
GENETIC ALGORITHM v2 - Gene-Level Antibiotic Resistance Prediction
Dự đoán GEN ĐỘT BIẾN nào gây kháng thuốc ở từng chủng vi khuẩn
Sử dụng: Phylogenetic Tree + H37Rv reference genome + TB-Profiler resistance database

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
from Bio import Phylo, SeqIO, Align
from ete3 import Tree
from sklearn.metrics import accuracy_score, precision_recall_fscore_support
from sklearn.model_selection import train_test_split
import random
import logging
import json
import warnings
import re
from pathlib import Path
from typing import Dict, List, Tuple

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
    'num_generations': 10,
    'crossover_prob': 0.7,
    'mutation_prob': 0.2,
    'tournament_size': 3,
    'test_size': 0.2
}

# =================== GENE REFERENCE LOADER ===================
class GeneReferenceLoader:
    """Load H37Rv reference genome + gene mapping"""
    
    def __init__(self, config: Dict):
        self.config = config
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
        import re
        from pathlib import Path
        
        detected = {}
        sample_fasta = Path(sample_dir) / f"{genome_id}.fasta"
        
        if not sample_fasta.exists():
            logging.debug(f"  Sample file not found: {sample_fasta}")
            return detected
        
        try:
            # Read sample sequence
            for record in SeqIO.parse(str(sample_fasta), 'fasta'):
                sample_seq = str(record.seq)
                break
        except Exception as e:
            logging.debug(f"  Error reading sample {genome_id}: {e}")
            return detected
        
        # Initialize gene dict
        for gene_id in self.gene_mapping.keys():
            detected[gene_id] = []
        
        # Check each known mutation locus
        for codon_pos, mutations in self.mutation_loci.items():
            if codon_pos >= len(sample_seq):
                continue
            
            # Check if position differs from reference
            if sample_seq[codon_pos] != self.h37rv_seq[codon_pos]:
                # Mutation is present
                for mut_info in mutations:
                    gene_id = mut_info['gene_id']
                    mutation_str = mut_info['mutation']
                    
                    if mutation_str not in detected[gene_id]:
                        detected[gene_id].append(mutation_str)
        
        return detected

# =================== DATA LOADER ===================
class DataLoader:
    """Load phylogenetic tree + metadata"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.tree = None
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
                logging.info(f"  ✓ Tree (ete3): {len(self.leaf_dict)} leaves")
            except:
                self.tree = Phylo.read(self.config['tree_file'], 'nexus')
                self.leaf_dict = {term.name: term for term in self.tree.get_terminals()}
                logging.info(f"  ✓ Tree (Phylo): loaded")
        except Exception as e:
            logging.warning(f"  ⚠️ Load tree error: {e}")
            self._create_demo_tree()
    
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
    bounds = [
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
    ]
    ind = [random.uniform(low, high) for low, high in bounds]
    return creator.Individual(ind)

# =================== GENE-LEVEL PREDICTOR ===================
class GeneLevelPredictor:
    """Predict resistance-causing mutations at gene level"""
    
    param_names = [
        "mutation_rate", "hgt_probability", "antibiotic_pressure", "efflux_pump",
        "fitness_cost", "recombination", "mobile_elements", "epistasis",
        "bottleneck", "transmission", "biofilm", "stress_response"
    ]
    
    def __init__(self, gene_loader: GeneReferenceLoader, data_loader: DataLoader):
        self.gene_loader = gene_loader
        self.data_loader = data_loader
        # Preload detected mutations for all samples
        self.detected_mutations_cache = {}
    
    def evaluate(self, individual):
        """
        Evaluate fitness by comparing observed mutation profiles to true labels.
        
        Fitness = accuracy of predicting resistance based on detected mutations
        vs. true resistance labels from metadata.
        """
        params = {name: val for name, val in zip(self.param_names, individual)}
        
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
            
            # Predict resistance based on detected mutations
            num_mutations = sum(len(muts) for muts in detected_muts.values())
            
            # Confidence-weighted mutation count
            total_confidence = 0
            for gene_id, mutations in detected_muts.items():
                for mut_str in mutations:
                    # Find confidence score for this mutation
                    for codon_pos, mut_infos in self.gene_loader.mutation_loci.items():
                        for mut_info in mut_infos:
                            if mut_info['gene_id'] == gene_id and mut_info['mutation'] == mut_str:
                                total_confidence += mut_info['confidence']
            
            # Predict resistant if sufficient evidence
            mutation_threshold = 0.5  # At least one confident mutation
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
                # Find confidence score for this mutation
                confidence = 0.0
                for mut_info in gene_info['mutations']:
                    if mut_info['mutation'] == mutation_str:
                        confidence = mut_info['confidence']
                        break
                
                # Calculate probability using GA parameters
                prob = (
                    params['mutation_rate'] * 0.1 +
                    params['antibiotic_pressure'] * 0.2 +
                    (1 - params['fitness_cost']) * confidence * 0.7
                )
                
                # Clamp probability to [0, 1]
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
    
    gene_loader = GeneReferenceLoader(CONFIG).load_all()
    data_loader = DataLoader(CONFIG).load_all()
    
    setup_ga_toolbox()
    predictor = GeneLevelPredictor(gene_loader, data_loader)
    
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
    
    min_vals = [0.001, 0.0, 0.5, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    max_vals = [0.05, 0.3, 3.0, 2.0, 0.5, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 2.0]
    
    def check_bounds(min_v, max_v):
        def decorator(func):
            def wrapper(*args, **kwargs):
                offspring = func(*args, **kwargs)
                for child in offspring:
                    for i in range(len(child)):
                        child[i] = np.clip(child[i], min_v[i], max_v[i])
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
    
    return hof, logbook

# =================== EXPORT GENE-LEVEL PREDICTIONS ===================
def export_gene_predictions(hof, logbook, predictor, data_loader, gene_loader):
    """Export gene-level resistance predictions"""
    
    logging.info("\n=== EXPORTING GENE-LEVEL PREDICTIONS ===")
    
    best = hof[0]
    best_params = {name: float(best[i]) for i, name in enumerate(predictor.param_names)}
    
    # First, log detected mutations per sample for diagnostics
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

# Note: This module contains GA classes and run function only.
# The pipeline orchestration (fetch -> visualize -> quantum -> predict -> report)
# is handled by `main.py`. Import `run_genetic_al_gene_level` and call from there.
