import os
import logging
from pathlib import Path
from Bio import Entrez
import pandas as pd
import json
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PyPDF2 import PdfMerger

import biodata
import genetic_al

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

CONFIG = {
	'h37rv_id': 'NC_000962.3',
	'h37rv_fasta': 'h37rv_ref.fasta',
	'tree_file': 'MTBC_SNP_align_200.phy_phyml_tree.nexus',
	'gene_mapping_file': 'gene_mapping_tb.json'
}


def fetch_reference(email: str = 'akamipersona1111@gmail.com'):
	"""Ensure H37Rv reference FASTA exists locally; fetch from NCBI if missing."""
	fasta_path = Path(CONFIG['h37rv_fasta'])
	if fasta_path.exists():
		logging.info(f"Reference fasta already exists: {fasta_path}")
		return str(fasta_path)

	logging.info(f"Fetching H37Rv ({CONFIG['h37rv_id']}) from NCBI...")
	Entrez.email = email
	try:
		handle = Entrez.efetch(db='nuccore', id=CONFIG['h37rv_id'], rettype='fasta', retmode='text')
		seq = handle.read()
		with open(fasta_path, 'w', encoding='utf-8') as f:
			f.write(seq)
		logging.info(f"Saved reference to {fasta_path}")
		return str(fasta_path)
	except Exception as e:
		logging.error(f"Failed to fetch reference: {e}")
		raise


def prepare_gene_mapping():
	"""Run gene mapping preparation using consolidated biodata module."""
	logging.info("Preparing gene mapping (biodata.setup_database)...")
	# biodata.setup_database() handles all: fetch genome, create mapping, visualize tree
	biodata.setup_database()
	if not Path(CONFIG['gene_mapping_file']).exists():
		raise FileNotFoundError(f"Expected gene mapping not found: {CONFIG['gene_mapping_file']}")
	logging.info("Gene mapping ready.")


def visualize_tree():
	"""Call biodata.visualize to render tree PDF/PNG."""
	logging.info("Rendering phylogenetic tree PDF/PNG...")
	biodata.visualize(file=CONFIG['tree_file'])


def quantum_simulation_placeholder():
	"""Placeholder for future quantum simulation integration."""
	logging.info("Running quantum simulation placeholder (no-op)...")
	out = Path('quantum_placeholder.txt')
	out.write_text('Quantum simulation placeholder. Implement quantum pipeline here.')
	return str(out)


def run_prediction():
	"""Run GA-based prediction by invoking genetic_al.run_genetic_al_gene_level()."""
	logging.info("Starting gene-level GA prediction (genetic_al)...")
	hof, logbook = genetic_al.run_genetic_al_gene_level()
	logging.info("GA finished")
	return hof, logbook


def generate_report(pdf_path: str = 'genetic_al_report.pdf'):
	"""Create a multi-page PDF report with tree image, top genes bar chart, and GA fitness progression."""
	logging.info("Generating report PDF...")
	preds_file = Path('genetic_al_gene_predictions.csv')
	stats_file = Path('genetic_al_gene_evolution_stats.csv')
	params_file = Path('genetic_al_best_gene_params.json')
	# Create body pages into a temporary PDF, then merge tree PDF (if present) at front
	body_pdf = Path('report_body.pdf')
	with PdfPages(str(body_pdf)) as pdf:
		# Page 1 (body): Top genes bar chart
		if preds_file.exists():
			df = pd.read_csv(preds_file)
			if not df.empty:
				top_genes = df['gene_name'].value_counts().head(10)
				fig, ax = plt.subplots(figsize=(11, 6))
				top_genes.plot(kind='bar', ax=ax, color='tab:blue')
				ax.set_title('Top Detected Resistance Genes')
				ax.set_ylabel('Counts')
				ax.set_xlabel('Gene')
				plt.tight_layout()
				pdf.savefig(fig)
				plt.close(fig)

		# Page 2 (body): GA fitness progression
		if stats_file.exists():
			try:
				df_stats = pd.read_csv(stats_file)
				if 'avg' in df_stats.columns:
					fig, ax = plt.subplots(figsize=(11, 6))
					ax.plot(df_stats['avg'], label='avg')
					if 'max' in df_stats.columns:
						ax.plot(df_stats['max'], label='max')
					ax.set_title('GA Fitness Progression')
					ax.set_xlabel('Generation')
					ax.set_ylabel('Fitness')
					ax.legend()
					plt.tight_layout()
					pdf.savefig(fig)
					plt.close(fig)
			except Exception as e:
				logging.warning(f"Could not plot GA stats: {e}")

		# Page 3 (body): Best parameters summary
		if params_file.exists():
			with open(params_file, 'r') as f:
				params = json.load(f)
			fig, ax = plt.subplots(figsize=(11, 8))
			ax.axis('off')
			txt = '\n'.join([f"{k}: {v:.4f}" for k, v in params.items()])
			ax.text(0.01, 0.99, 'Best GA Parameters\n\n'+txt, va='top', ha='left', fontsize=10, family='monospace')
			pdf.savefig(fig)
			plt.close(fig)

	# Merge tree PDF (if exists) in front of the body PDF
	final_pdf = Path(pdf_path)
	merger = PdfMerger()
	try:
		merger.append(str(body_pdf))
		merger.write(str(final_pdf))
		merger.close()
		logging.info(f"Report saved to {final_pdf}")
	except Exception as e:
		logging.error(f"Failed to merge PDFs: {e}")
		# fallback: move body_pdf to final name
		body_pdf.replace(final_pdf)

	# cleanup
	try:
		if body_pdf.exists():
			body_pdf.unlink()
	except Exception:
		pass

	return str(final_pdf)


if __name__ == '__main__':
	# High-level orchestrator
	# try:
	fetch_reference()
	prepare_gene_mapping()
	visualize_tree()
	quantum_simulation_placeholder()
	hof, logbook = run_prediction()
	report = generate_report()
	logging.info('Pipeline complete. Report: %s', report)
	# except Exception as e:
		# logging.error('Pipeline failed: %s', e)

