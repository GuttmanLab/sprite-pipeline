from collections import Counter
import gzip
import operator
import pysam
import re
import sys

# Program for checking barcoding success rate

def main():
    lig = LigationEfficiency()
    lig.count_barcodes(sys.argv[1])
    lig.print_to_stdout()


class LigationEfficiency:

    def __init__(self):
        self._aggregate_count = Counter()
        self._position_count = Counter()
        self._pattern = re.compile('\[([a-zA-Z0-9_\-]+)\]')
        self._total = 0

    def count_barcodes(self, filename):
        if filename.lower().endswith(".bam"):
            self.count_barcodes_in_bam_file(filename)
        elif filename.lower().endswith((".fastq", ".fq")):
            self.count_barcodes_in_fastq_file(filename)
        elif filename.lower().endswith((".fastq.gz", ".fq.gz")):
            self.count_barcodes_in_fastqgz_file(filename)

    def count_barcodes_in_bam_file(self, bamfile):
        with pysam.AlignmentFile(bamfile, "rb") as f:
            for read in f.fetch(until_eof=True):
                self.count_barcodes_in_name(read.query_name)
                self._total += 1

    def count_barcodes_in_fastq_file(self, fastqfile):
        with open(fastqfile, "r") as f:
            for line in f:
                self.count_barcodes_in_name(line)
                next(f)
                next(f)
                next(f)
                self._total += 1

    def count_barcodes_in_fastqgz_file(self, fastqgzfile):
        with gzip.open(fastqgzfile, "r") as f:
            for line in f:
                self.count_barcodes_in_name(line)
                next(f)
                next(f)
                next(f)
                self._total += 1

    def count_barcodes_in_name(self, name):
        barcodes = self._pattern.findall(name)
        num_found = 0
        pos = 0
        for barcode in barcodes:
            if barcode != "NOT_FOUND":
                num_found += 1
                self._position_count[pos] += 1
            pos += 1
        self._aggregate_count[num_found] += 1

    def print_to_stdout(self):

        counts = sorted(self._aggregate_count.items(),
                        key=operator.itemgetter(0))

        for num_barcodes, count in counts:
            pct = "{0:.1f}%".format(100.0 * count / self._total)
            barcode = "barcode" if num_barcodes == 1 else "barcodes"
            print (str(count) + " (" + pct + ") reads found with " +
                str(num_barcodes) + " " + barcode + ".")

        print
        counts = sorted(self._position_count.items(),
                        key=operator.itemgetter(0))

        for position, count in counts:
            pct = "{0:.1f}%".format(100.0 * count / self._total)
            print (str(count) + " (" + pct + ") barcodes found in position " +
                str(position + 1) + ".")

if __name__ == "__main__":
    main()

