#!/usr/bin/python
import argparse
import pysam

def main():
    args = parse_arguments()
    filter_reads(args)

def parse_arguments():

    parser = argparse.ArgumentParser(description =
            "This program removes reads from a BAM file according to the " +
            "number of mismatches and ambiguous bases they contain. This " +
            "program will not adjust the read flags for paired reads if one " +
            "read of the pair is removed and the other retained.")
    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'Input BAM file')
    parser.add_argument('-o', '--output', action = 'store', metavar = 'FILE',
                        help = 'Output BAM file containing reads that pass')
    parser.add_argument('--edit_max', action = 'store', metavar = 'Y',
                        type = int, default = 0,
                        help = ('If the edit distance between the read ' +
                                'sequence and the reference sequence falls ' +
                                'within [X, Y], keep the read. Otherwise ' +
                                'remove it. Default = 0'))
    parser.add_argument('--edit_min', action = 'store', metavar = 'X',
                        type = int, default = 0,
                        help = ('If the edit distance between the read ' +
                                'sequence and the reference sequence falls ' +
                                'within [X, Y], keep the read. Otherwise ' +
                                'remove it. Default = 0'))
    parser.add_argument('--mapq_min', action = 'store', metavar = 'M',
                        type = int, default = 0,
                        help = ('If the MAPQ score of this read falls within ' +
                                '[M, N], keep the read. Otherwise remove ' +
                                'it. Default = 0'))
    parser.add_argument('--mapq_max', action = 'store', metavar = 'N',
                        type = int, default = 255,
                        help = ('If the MAPQ score of this read falls within ' +
                                '[M, N], keep the read. Otherwise remove ' +
                                'it. Default = 255'))
    parser.add_argument('--paired', action = 'store_true',
                        help = ('Set if the BAM file contains a paired-end ' +
                                'alignment. Defaults to single-read alignment.'))
    args = parser.parse_args()

    if args.edit_max < args.edit_min:
        exit("Invalid edit-distance parameters. Max: " +
             str(args.edit_max) + ", min: " + str(args.edit_min))

    if args.mapq_max < args.mapq_min:
        exit("Invalid MAPQ parameters. Max: " +
             srt(args.mapq_max) + ", min: " + str(args.mapq_min))

    return args

def filter_reads(args):
    if args.paired:
        filter_paired_reads(args)
    else:
        filter_single_reads(args)

def filter_single_reads(args):
    with pysam.AlignmentFile(args.input, "rb") as input_file,\
         pysam.AlignmentFile(args.output, "wb", template = input_file) as output_file:

        for read in input_file.fetch(until_eof = True):
            is_valid_read = (not read.is_unmapped and
                             has_valid_edit_distance(read, args.edit_min, args.edit_max) and
                             has_valid_mapq_score(read, args.mapq_min, args.mapq_max))
            if is_valid_read:
                output_file.write(read)

def filter_paired_reads(args):
    with pysam.AlignmentFile(args.input, "rb") as input_file,\
         pysam.AlignmentFile(args.output, "wb", template = input_file) as output_file:

        reads = input_file.fetch(until_eof = True)
        for read in reads:
            mate = reads.next()
            if read.query_name != mate.query_name:
                exit("Successive read-names do not match. Is your BAM file " +
                     "sorted by read-name?")
            if all([not read.is_unmapped,
                    not mate.is_unmapped,
                    has_valid_edit_distance(read, args.edit_min, args.edit_max),
                    has_valid_edit_distance(mate, args.edit_min, args.edit_max),
                    has_valid_mapq_score(read, args.mapq_min, args.mapq_max),
                    has_valid_mapq_score(mate, args.mapq_min, args.mapq_max)]):
                
                output_file.write(read)
                output_file.write(mate)

def has_valid_edit_distance(read, lo, hi):
    score = int(read.get_tag("NM"))
    return lo <= score <= hi

def has_valid_mapq_score(read, lo, hi):
    return lo <= read.mapping_quality <= hi

if __name__ == "__main__":
    main()
