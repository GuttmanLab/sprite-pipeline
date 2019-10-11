#!/usr/bin/python
import argparse
import pysam
import assembly

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
    parser.add_argument('--assembly', metavar = "ASSEMBLY",
                        action = "store",
                        choices = ["mm9", "mm10", "hg19", "hg38", "none"],
                        default = "none",
                        help = "The genome assembly. (default mm9)")
    return parser.parse_args()


def add_chr_to_bam_header(bam_path, chroms):
    '''Add chr to alignments done with ensembl reference

        Args:
            bam_path (str): path to bam file
            chroms (list): list of core chromosomes to retain in the header

        Note:
            Assumes main chromosomes come first followed by alt, that will be removed 

        Return:
            header dict for pysam
    '''
    #get bam header
    with pysam.AlignmentFile(bam_path, 'rb') as input_file:
        bam_header = input_file.header.to_dict()
    #edit bam header
    to_delete = []
    for n, i in enumerate(bam_header['SQ']):
        chrom = i['SN']
        if 'chr' not in chrom:
            if 'chr' + chrom in chroms:
                bam_header['SQ'][n]['SN'] = 'chr' + chrom
            elif chrom == 'MT':
                bam_header['SQ'][n]['SN'] = 'chrM'
            else:
                to_delete.append(n)
    #need to delete largest index first
    for i in sorted(to_delete, reverse=True):
        del bam_header['SQ'][i]

    return(bam_header)



def filter_reads(args):
    '''Filter bam

    Params:
        args: argparse arguments

    Notes:
        Will filter reads without a full barcode, will remove alt chromosome and MT
        and if aligned with ensembl (1, 2, ..., X, Y) will change chromosomes to ucsc
        scheme (chr1, chr2, ..., chrX, chrY)

    '''

    if args.assembly != 'none':
        #filter out reads that do not map to main assembly
        chroms = list(assembly.build(args.assembly, 1)._chromsizes.keys())
        #add chr to bam header
        bam_header = add_chr_to_bam_header(args.input, chroms)

    # with pysam.AlignmentFile("/mnt/data/2p5-4.RNA.Aligned.sortedByCoord.out.bam.featureCounts.bam", "rb") as input_file,\
    #      pysam.AlignmentFile("/mnt/data/2p5-4.RNA.test.bam", "wb", header = bam_header) as output_file:
    filtered_count = 0
    out_count = 0
    with pysam.AlignmentFile(args.input, "rb") as input_file:
        if args.assembly == 'none':
            output_file = pysam.AlignmentFile(args.output, "wb", template = input_file)
            for read in input_file.fetch(until_eof = True):
                if not "NOT_FOUND" in read.query_name:
                    output_file.write(read)
                    out_count += 1
                else:
                    filtered_count += 1
        else:
            output_file = pysam.AlignmentFile(args.output, "wb", header = bam_header)
            for read in input_file.fetch(until_eof = True):
                ref_name = read.reference_name if 'chr' in read.reference_name else 'chr' + read.reference_name 
                if not "NOT_FOUND" in read.query_name and ref_name in chroms:
                    output_file.write(read)
                    out_count += 1
                else:
                    filtered_count += 1

        output_file.close()

    print('Filtered reads:', filtered_count)
    print('Written out reads:', out_count)



if __name__ == "__main__":
    main()
