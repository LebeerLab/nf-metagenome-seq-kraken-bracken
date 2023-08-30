#!/usr/bin/env python3
# Extract sequences classified as specified taxa in Kraken
# Script by JK

# Modules
import argparse
from argparse import RawTextHelpFormatter
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import gzip

# Functions
def err(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# Usage
parser = argparse.ArgumentParser(
	formatter_class=RawTextHelpFormatter,
	description='extracts specified taxa sequences from Kraken classification - script by JK',
	usage='\n  %(prog)s --kraken [kraken.out] --taxids [file] READS > extracted.fastq'
	'\n  %(prog)s --kraken [kraken.out] --taxids [file] READS | pigz > extracted.fastq.gz')
parser.add_argument('reads', metavar='READS', help='single FASTQ file')
parser.add_argument('--kraken', metavar='FILE', help='output file from Kraken')
parser.add_argument('--taxids', metavar='FILE', help='file with list of taxa ID to extract')
parser.add_argument('--version', action='version', version='0.1')
args = parser.parse_args()

# Setup list of taxa IDs of relevance
taxids = [line.rstrip('\n') for line in open(args.taxids)]

# Obtain read IDs of reads classified as one of the taxa of relevance
count = 0
seqids = []
with open(args.kraken) as k:	
	for line in k:
		count += 1
		err('Checking reads ... {}'.format(count), end='\r')		
		row = line.split()
		if row[2] in taxids:
			seqids.append(row[1])

# Read counter
err('\nFound {} of {} reads.'.format(len(seqids), count))

# Print reads to stdout where read ID matches reads of relevance
count = 0

def extract_reads(handle, count):
    for title, seq, qual in FastqGeneralIterator(handle):
        id = title.split()[0]
        taxid = title.split("kraken:taxid|")[1]
        if id in seqids:
            count += 1
            err('Writing to file ... {}'.format(count), end='\r')
            print("%s\t%s" % (taxid, seq))
            #print("@%s\n%s\n+\n%s" % (title, seq, qual))

try: 		# Handle gzipped reads
	with gzip.open(args.reads, 'rt') as in_handle:
		extract_reads(in_handle, count)
		
except:		# Works for non-gzipped reads
	with open(args.reads, 'r') as in_handle:
		extract_reads(in_handle, count)
