#!/usr/bin/python
import argparse
import pysam

def main():
    args = parse_arguments()
    filter_reads(args)

def parse_arguments():
    parser = argparse.ArgumentParser(description =
            "This program removes reads from a BAM file if they do not " +
            "all tags. It does this by searching for the string " +
            "\"NOT_FOUND\" in the query name.")
    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'Input BAM file')
    parser.add_argument('-o', '--output', action = 'store', metavar = 'FILE',
                        help = 'Output BAM file containing reads that pass')
    return parser.parse_args()

def filter_reads(args):
    with pysam.AlignmentFile(args.input, "rb") as input_file,\
         pysam.AlignmentFile(args.output, "wb", template = input_file) as output_file:

        for read in input_file.fetch(until_eof = True):
            if not "NOT_FOUND" in read.query_name:
                output_file.write(read)

if __name__ == "__main__":
    main()
