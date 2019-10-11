
from collections import defaultdict

cluster ='/mnt/data/combined.clusters'

count = 0
bedgraph = defaultdict()
with open(cluster, 'r') as in_clusters:
    for line in clusters:
        if count < 10:
            barcode, *cluster = line.rstrip().split('\t')
            count += 1
            for i in cluster:
                if 'RNA' in i:
                    #DPM(+)_chr1:3095318-3095358
                    #chr1:3086404-3086440
                    chr, start_end = i.split('_')[-1].split(':')
                    start, end = start_end
                    bedgraph[]
        else:
            break
