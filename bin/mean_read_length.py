#!/usr/bin/env python3
"""
This script returns the average read length over one or multiple fasta/fastq files using seqkit.

Author: Florian De Temmerman (fbdtemme) <floriandetemmerman@gmail.com>
"""

import argparse
import gzip
import multiprocessing
from concurrent.futures import ThreadPoolExecutor
from functools import partial
import mimetypes

import numpy as np
from Bio import SeqIO

N_CORES = multiprocessing.cpu_count()

def get_average_read_length_seqkit(fastq_file):
    encoding = mimetypes.guess_type(fastq_file)[1]  # uses file extension
    _open = partial(gzip.open, mode='rt') if encoding == 'gzip' else open

    with _open(fastq_file) as f:
        sizes = [len(rec) for rec in SeqIO.parse(f, "fastq")]

    result = np.mean(sizes)
    return result

def main(args, threads=N_CORES):
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_objs = []
        for read in args.reads:
            future_objs.append(executor.submit(get_average_read_length_seqkit, read))

    mean_length = np.mean([f.result() for f in future_objs])

    print(int(mean_length), end='')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=True, 
    	description='This tool returns the mean read length of a list of fastq/fasta files')
    parser.add_argument('reads', help='Path to the matrix of counts in mtx format (gzipped or not)', nargs='+')
    parser.add_argument('--threads', help="The number of threads to use.", type=int, default=4)

    args = parser.parse_args()
    main(args)


