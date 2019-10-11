import argparse
import cluster as c

def main():
    args = parse_arguments()
    clusters = c.get_clusters(args.input, args.num_tags, args.genome_1, 
                              args.genome_2)
    c.write_clusters_to_file(clusters, args.output)
    print("done")

def parse_arguments():
    parser = argparse.ArgumentParser(
        description = 'Generates a clusters file from a BAM file.')
    parser.add_argument('-i', '--input',
                        metavar = "FILE",
                        action = "store",
                        nargs='+',
                        required=True,
                        help = "The input BAM file(s).")
    parser.add_argument('-o', '--output',
                        metavar = "FILE",
                        action = "store",
                        required=True,
                        help = "The output clusters file.")
    parser.add_argument('-n', '--num_tags',
                        metavar = 'INT',
                        type = int,
                        action = 'store',
                        required=True,
                        help = "The number of tags contained in the barcode " +
                               "of each BAM record.")
    parser.add_argument('-g1', '--genome_1',
                        type = str,
                        action = 'store',
                        default = "None",
                        help = "Genome 1 strain used with SNPsplit")
    parser.add_argument('-g2', '--genome_2',
                        type = str,
                        action = 'store',
                        default = "None",
                        help = "Genome 2 strain used with SNPsplit")
    return parser.parse_args()

if __name__ == "__main__":
    main()
