#!/usr/bin/env python
"""
This script parses a single cell dataset in MTX format to dataframe (tab delimited), loompy and hd5 formats. 
The script will also place all the output files as well as the orignal MTX files in a output directory
given by the user. The script is useful to parse the output of tools such as STARSolo, Kallisto bustools
and Salmon Alevin. A set of basic QC plots will also be generated in the output folder. Users can use
an option to convert Ensembl ids to gene names if needed. 

Type:

postprocessing.py --help

To get extra information on the input parameters. 

How to use it:

postprocessing.py --input data.mtx.gz --features features.txt --barcodes barcodes.txt --output out 

After this you should find in the out directory the input matrix in the original format as well as the
other formats. The QC figures will be placed in a subfolder called figures.  

Author: Jose Fernandez Navarro (jfnavarro) <jc.fernandez.navarro@gmail.com>
"""
from __future__ import print_function
import os
import sys
import argparse
import scipy.io
import shutil
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def main(matrix_file, features_file, barcodes_file, transpose, gene_names, organism, output):

	if not os.path.isfile(matrix_file):
		sys.stderr.write("Error, input matrix {} is not a file\n".format(matrix_file))
		sys.exit(1)

	if not os.path.isfile(features_file):
		sys.stderr.write("Error, input features {} is not a file\n".format(features_file))
		sys.exit(1)

	if not os.path.isfile(barcodes_file):
		sys.stderr.write("Error, input barcodes {} is not a file\n".format(barcodes_file))
		sys.exit(1)
	
	if (organism is not None and not gene_names) or (gene_names and organism is None):
		sys.stderr.write("Error, --organism must be used together with --gene-names\n")
		sys.exit(1)

	# TODO check that organism is valid if gene_names is True

	if os.path.isdir(output):
		print("Ouput directory {} exists, files will be overwritten..".format(output))
	os.makedirs(os.path.abspath(output), exist_ok=True)

    # Parse the features
	features = pd.read_csv(features_file, sep='\t', header=None, index_col=None).iloc[:,0]

    # Parse the barcodes
	barcodes = pd.read_csv(barcodes_file, sep='\t', header=None, index_col=None).iloc[:,0]

    # Parse the matrix of counts
	matrix = scipy.io.mmread(matrix_file)

    # Create Pandas
	df = pd.DataFrame(matrix.toarray().transpose() if transpose else matrix.toarray(), 
	                  index=features, columns=barcodes)

    # Convert Ensembl ids to gene names
	if gene_names:
		annot = sc.queries.biomart_annotations(
        	organism, ["ensembl_gene_id", "external_gene_name", "chromosome_name"],
    	).set_index("ensembl_gene_id")
		df.index = [x.split('.')[0] for x in df.index]
		intersect = np.intersect1d(df.index, annot.index)
		df = df.loc[intersect, :]
		df.index = annot.loc[intersect, 'external_gene_name']	

	# Rename rows/columns
	df.index = df.index.rename("Genes")
	df.columns = df.columns.rename("Barcodes")

	# Create Anndata and rename duplicates (if any)
	adata = sc.AnnData(df.transpose())
	adata.var_names_make_unique()
	adata.obs_names_make_unique()
	df.index = adata.var.index
	df.columns = adata.obs.index

	# Write Pandas
	df.to_csv(os.path.join(output, 'matrix.tsv'), sep='\t', header=True, index=True)

	# Write hdf5
	df.to_hdf(os.path.join(output, 'matrix.h5'), key='df', mode='w')

	# Write mtx
	extension = os.path.splitext(matrix_file)[1] if not matrix_file.endswith('.mtx') else ''
	shutil.copy(matrix_file, os.path.join(output, 'matrix.mtx{}'.format(extension)))
	shutil.copy(features_file, os.path.join(output, 'features.mtx.txt'))
	shutil.copy(barcodes_file, os.path.join(output, 'barcodes.mtx.txt'))

	# Write h5ad
	adata.write(os.path.join(output, 'matrix.h5ad'))

	# Write loom
	adata.write_loom(os.path.join(output, 'matrix.loom'))

	# Some basic QC plots
	sc.settings.figdir = os.path.join(output, sc.settings.figdir)
	sc.set_figure_params(dpi=180, color_map='viridis_r')
	sc.settings.authoshow = False
	sc.settings.autosave = True
	sc.pl.highest_expr_genes(adata, n_top=20)
	adata.var['mt'] = adata.var_names.str.startswith('MT-')
	sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, multi_panel=True)
	sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='_mt.pdf')
	sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='_genes.pdf')
	fig, axs = plt.subplots(1, 2, figsize=(15, 4))
	sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
	sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[1])
	fig.savefig(os.path.join(sc.settings.figdir, 'distplot.pdf'))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True, 
    	description='This tool parses a matrix of counts in mtx format to other formats (hdf5, loompy and tsv)')
    parser.add_argument('--matrix', 
    	metavar="[FILE]", 
    	required=True,
    	help='Path to the matrix of counts in mtx format (gzipped or not) where genes are rows and barcodes are columns')
    parser.add_argument('--features',
    	metavar="[FILE]",
    	required=True,
    	help='Path to the column names (features) in tab delimited text format')
    parser.add_argument('--barcodes',
    	metavar="[FILE]",
    	required=True,
    	help='Path to the row names (barcodes) in tab delimited text format')
    parser.add_argument('--transpose',
	    action="store_true",
		default=False,
    	help='Use this if the input matrix has genes as columns and barcodes as rows')
    parser.add_argument('--gene-names',
	    action="store_true",
		default=False,
    	help='Use this if you want to convert the Ensembl ids to gene names')
    parser.add_argument('--organism',
    	metavar="[STR]",
    	required=False,
		default=None,
    	help='BioMart organism to use when using --gene-names. Options are: hsapiens, mmusculus, drerio, etc..')
    parser.add_argument('--output',
    	metavar="[STRING]",
    	required=True,
    	help='Path to the output folder')
    
    args = parser.parse_args()
    main(args.matrix, args.features, args.barcodes, args.transpose, args.gene_names, args.organism, args.output)
