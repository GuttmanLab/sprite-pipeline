import argparse
import assembly
import numpy as np
import numpy.ma as ma

def main():
    args = parse_arguments()

    contacts = np.loadtxt(args.heatmap)
    resolution = args.resolution
    assem = assembly.build(args.assembly, resolution)
    hubfile = args.hub

    compute_hub_avgs_and_print(contacts, hubfile, assem)


def compute_hub_avgs_and_print(contacts, hubfile, assem):

    with open(hubfile, 'r') as f:
        for line in f:
            fields = line.strip().split()
            chrom = fields[0]
            start_pos = fields[1]
            end_pos = fields[2]
            idx = assem.get_index(chrom, start_pos)

            chrom_start_idx = assem.get_offset(chrom)
            chrom_end_idx = chrom_start_idx + (assem.get_size(chrom) // assem._resolution)

            chrom_mask = ma.make_mask(np.zeros(contacts.shape[0], shrink = False)
            chrom_mask[chrom_start_idx, chrom_end_idx + 1] = True
            avg = np.mean(np.mean(ma.array(contacts[idx], mask = chrom_mask)))
            print "\t".join([chrom, start_pos, end_pos, avg])


def parse_arguments():
    parser = argparse.ArgumentParser(
        description = "Generates a bedgraph of hub-contact frequencies.")

    parser.add_argument('--heatmap',
                        metavar = "FILE",
                        action = "store",
                        help = "A heatmap file (raw contacts).")
    parser.add_argument('--resolution',
                        metavar = "INT",
                        type = int,
                        action = "store",
                        default = 1000000,
                        help = "Resolution in bp. Defaults to 1000000 (1 MB).")
    parser.add_argument('--hub',
                        metavar = "FILE",
                        action = "store",
                        help = "A BED file of hub regions.")
    parser.add_argument('--assembly',
                        metavar = "ASSEMBLY",
                        action = "store",
                        choices = ["mm9", "mm10", "hg19"],
                        default = "mm9",
                        help = "The genome assembly. (default mm9)")
    return parser.parse_args()


if __name__ == "__main__":
    main()
