#!/usr/bin/env python
"""
This script parses a single cell dataset in MTX format to dataframe (tab delimited) and hd5 formats. 
The script will also place all the output files as well as the orignal MTX files in a output directory
given by the user. The script is useful to parse the output of tools such as STARSolo, Kallisto bustools
and Salmon Alevin. 

How to use it:

python postprocessing.py --input data.mtx.gz --features features.txt --barcodes barcodes.txt --output out 

After this you should find in the out directory the input matrix in the original format as well as the
other formats. The input matrix in MTX format must be gzipped. 

Author: Jose Fernandez Navarro (jfnavarro) <jc.fernandez.navarro@gmail.com>
"""
from __future__ import print_function
import os
import sys
import argparse
import scipy.io
import shutil
import pandas as pd


def main(matrix_file, features_file, barcodes_file, output):

	if not os.path.isfile(matrix_file):
		sys.stderr.write("Error, input matrix {} is not a file\n".format(matrix_file))
		sys.exit(1)

	if not os.path.isfile(features_file):
		sys.stderr.write("Error, input features {} is not a file\n".format(features_file))
		sys.exit(1)

	if not os.path.isfile(barcodes_file):
		sys.stderr.write("Error, input barcodes {} is not a file\n".format(barcodes_file))
		sys.exit(1)

	if os.path.isdir(output):
		print("Ouput directory {} exists, files will be overwritten..".format(output))
	os.makedirs(os.path.abspath(output), exist_ok=True)

    # Parse the features
	features = pd.read_csv(features_file, sep='\t',
        header=None, index_col=None).values.flatten()

    # Parse the barcodes
	barcodes = pd.read_csv(barcodes_file, sep='\t',
    	header=None, index_col=None).values.flatten()

    # Parse the matrix of counts
	matrix = scipy.io.mmread(matrix_file)

    # Create Pandas
	df = pd.DataFrame(matrix.toarray().transpose(), index=features, columns=barcodes)

    # Write Pandas
	df.to_csv(os.path.join(output, 'matrix.tsv'), sep='\t', header=True, index=True)

    # Write hdf5
	df.to_hdf(os.path.join(output, 'matrix.h5'), key='df', mode='w')

    # Write mtx
	shutil.copy(matrix_file, os.path.join(output, 'matrix.mtx.gz' if matrix_file.endswith('gz') else 'matrix.mtx'))
	shutil.copy(features_file, os.path.join(output, 'features.mtx.txt'))
	shutil.copy(barcodes_file, os.path.join(output, 'barcodes.mtx.txt'))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=True, 
    	description='This tool parses a matrix of counts in mtx format to other formats (hdf5 and tsv)')
    parser.add_argument('--matrix', 
    	metavar="[FILE]", 
    	required=True,
    	help='Path to the matrix of counts in mtx format (gzipped or not)')
    parser.add_argument('--features',
    	metavar="[FILE]",
    	required=True,
    	help='Path to the column names (features) in text format')
    parser.add_argument('--barcodes',
    	metavar="[FILE]",
    	required=True,
    	help='Path to the row names (barcodes) in text format')
    parser.add_argument('--output',
    	metavar="[STRING]",
    	required=True,
    	help='Path to the output folder')
    
    args = parser.parse_args()
    main(args.matrix, args.features, args.barcodes, args.output)
