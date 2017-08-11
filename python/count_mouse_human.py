#!/usr/bin/env python

from __future__ import division
import argparse
import math

def main():
    """Counts the number of spurious mouse-human contacts in a bead-mixing
    SPRITE experiment. Used to estimate noise.

    Output is six numbers printed to stdout:
      1. Total number of mouse-human contacts
      2. Total number of mouse reads
      3. Total number of human reads
      4. Total number of contacts
      5. Fraction of contacts that are mouse-human
      6. Expected fraction of contacts that are mouse-human
    """
    args = parse_arguments()

    total_contacts = 0
    mouse_human_contacts = 0
    mouse_reads = 0
    human_reads = 0

    with open(args.input, 'r') as f:
        for line in f:
            reads = get_reads_from_line(line)

            cluster_size = len(reads)
            if cluster_size > args.max_cluster_size:
                continue

            no_dupl_reads = set(reads)
            cluster_size = len(no_dupl_reads)
            if cluster_size > 1:
                mouse_count = 0
                human_count = 0
                for read in no_dupl_reads:
                    organism = get_organism_from_read(read)
                    if organism == "mouse":
                        mouse_count += 1
                    elif organism == "human":
                        human_count += 1
		total_contacts += nCr(cluster_size, 2)
		mouse_human_contacts += mouse_count * human_count
		mouse_reads += mouse_count
		human_reads += human_count

    total_reads = mouse_reads + human_reads

    print mouse_human_contacts,\
          mouse_reads,\
          human_reads,\
          total_contacts,\
          mouse_human_contacts / total_contacts,\
          mouse_reads / total_reads * human_reads / total_reads * 2


def nCr(n, r):
    f = math.factorial
    return f(n) // f(r) // f(n-r)

def get_reads_from_line(line):
    """Parses a cluster-file line into a list of reads. Ignores the 0th field,
    which contains the barcode string"""
    return line.rstrip().split()[1:]

def get_organism_from_read(read):
    """Assumes a read is named something like chr1_human:344102"""
    return read.split(":")[0].split("_")[-1]

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Calculates noise in a human/mouse SPRITE experiment')

    parser.add_argument('--input',
        metavar = "FILE",
        action = "store",
        help = "The input SPRITE clusters file.")

    parser.add_argument('--max_cluster_size',
        metavar = 'INT',
        type = int,
        action = 'store',
        default = 10000,
        help = "Maximum cluster size to consider (default 10000)")

    return parser.parse_args()

if __name__ == "__main__":
    main()
