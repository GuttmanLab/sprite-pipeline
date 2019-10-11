import argparse
import cluster as c

def main():
    args = parse_arguments()
    clusters = c.get_clusters_fastq(args.input, args.num_tags, args.orientation, args.umi_len)
    c.write_clusters_to_file(clusters, args.output, args.unique)
    print("done")

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generates a clusters file from a fastq file.')
    parser.add_argument('-i', '--input',
                        metavar = "FILE",
                        action = "store",
                        help = "The input fastq file.")
    parser.add_argument('-o', '--output',
                        metavar = "FILE",
                        action = "store",
                        help = "The output clusters file.")
    parser.add_argument('-n', '--num_tags',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        help = "The number of tags contained in the barcode " +
                               "of each fastq record.")
    parser.add_argument('-u', '--umi_len',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        help = "UMI length")
    parser.add_argument('-r', '--orientation',
                        type = str,
                        action = 'store',
                        help = "Orientation 'r1' or 'r2'")
    parser.add_argument('-q', '--unique',
                        action='store_true',
                        help = "Keep only unique UMIs in clusters")
    return parser.parse_args()

if __name__ == "__main__":
    main()
