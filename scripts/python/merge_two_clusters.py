from collections import defaultdict
import sys

barcode_to_positions = defaultdict(list)

with open(sys.argv[1], 'r') as f:
    for line in f:
        fields = line.rstrip().split()
        barcode_to_positions[fields[0]] = fields[1:]

with open(sys.argv[2], 'r') as f:
    for line in f:
        fields = line.rstrip().split()
        print fields[0] + "\t" + "\t".join(set(barcode_to_positions[fields[0]] + fields[1:]))
        del barcode_to_positions[fields[0]]

for barcode, positions in barcode_to_positions.iteritems():
    print barcode + "\t" + "\t".join(positions)


