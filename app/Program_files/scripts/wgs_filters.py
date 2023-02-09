#!/usr/bin/env python
# coding: utf-8


from glob import glob
import pandas as pd
import numpy as np
import argparse
import time
from io import StringIO

parser = argparse.ArgumentParser(description= 'FJD CNV database step1')

parser.add_argument('-f', '--file', help='wgs variant file', required=True)
parser.add_argument('-c', '--canonical', help='canonical transcripts', required=True, action='store_true')
parser.add_argument('-af', '--frequency', help='allele frequency', required=True)
parser.add_argument('-r', '--region', help='genomic region', required=True, action='store_true')
parser.add_argument('-p', '--panel', help='virtual panel of genes', required=False)

start_time = time.time()

args = parser.parse_args()

df = pd.read_csv(args.file, sep="\t")


if args.canonical == True:
	df = df[df['CANONICAL'] == 'YES']


if args.region == True:
	df = df[(df['Genomic_region'] == 'EXONIC' ) | (df['Genomic_region'] == 'SPLICING' )]


if args.frequency != None:
	df = df[(df['gnomADg_AF_popmax'] <= float(args.frequency) ) | (df['gnomADg_AF_popmax'].isna())]


if args.panel != None:

    with open(args.panel) as f:
        genes = f.read().splitlines()
        print(genes)
        df = df[df['SYMBOL'].isin(genes)]

#output = StringIO()
#df.to_csv(output)
#print(output.getvalue())
df.to_csv( "/home/raquel/PriorR/tmp/wgs.tsv", sep  ='\t')

print(time.time() - start_time)
