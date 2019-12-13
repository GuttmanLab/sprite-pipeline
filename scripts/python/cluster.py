import pysam
import re
import gzip
import os
import sys

__author__ = "Noah Ollikainen, Charlotte A Lai, Peter Chovanec"

class Position:
    """This class represents a genomic position, with type of nucleic acid (DNA)

    Methods:
    - to_string(): Returns a string representation of this position in the form
      "DNA(feature)_chrX:1000"
    """

    def __init__(self, type, feature, chromosome, start_coordinate, end_coordinate):
        self._type = type
        self._feature = feature
        self._chromosome = chromosome
        self._start_coordinate = start_coordinate
        self._end_coordinate = end_coordinate

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        return (self._type == other._type and
                self._feature == other._feature and
                self._chromosome == other._chromosome and
                self._start_coordinate == other._start_coordinate and
                self._end_coordinate == other._end_coordinate)

    def __hash__(self):
        return hash((self._type, self._feature, self._chromosome, 
                     self._start_coordinate, self._end_coordinate))

    def to_string(self):
        try:
            out = self._type + "[" + self._feature + "]" + "_" + \
                self._chromosome + ":" + \
                str(self._start_coordinate) + "-" + str(self._end_coordinate)
        except:
            print(self._type, self._feature, self._chromosome)
            print('Elements are not as expect!')
            sys.exit()
        return out

class Cluster:
    """This class represents a barcoding cluster as a collection of genomic
    positions.

    The underlying data structure is a set, so duplicate positions are
    discarded.

    Methods:
    - add_position(position): Adds a genomic position to this cluster

    - size(): Returns the number of reads or positions in this cluster

    - to_string(): Returns a string representation of this cluster as a
      tab-delimtited series of positions. See Position#to_string for how
      positions are represented as strings.
    """

    def __init__(self):
        self._positions = set()

    def add_position(self, position):
        self._positions.add(position)

    def size(self):
        return len(self._positions)

    def to_string(self):
        position_strings = [position.to_string() for position in self._positions]
        return "\t".join(position_strings)


class Clusters:
    """This class represents a collection of barcoding clusters.

    Methods:
    - get_cluster(barcode): Returns the cluster that corresponds to the given
      barcode. If the cluster does not exist, it is initialized (with zero
      positions), and this cluster.to_string()cluster.to_string()cluster.to_string()cluster.to_string()cluster.to_string()empty cluster is returned.

    - add_position(barcode, position): Adds the position to the cluster
      that corresponds with the given barcodes

    - to_strings(): Returns an iterator over the string representations of all
      of the contained clusters.
    """
    def __init__(self):
        self._clusters = {}

    def get_cluster(self, barcode):
        if barcode not in self._clusters:
            self._clusters[barcode] = Cluster()
        return self._clusters[barcode]

    def add_position(self, barcode, position):
        self.get_cluster(barcode).add_position(position)

    def to_strings(self):
        for barcode, cluster in self._clusters.items():
            yield barcode + '\t' + cluster.to_string()

    def remove_cluster(self, barcode):
        del self._clusters[barcode]

    def unique(self):
        for barcode, cluster in self._clusters.items():
            yield barcode + '\t' + cluster.unique()



def get_clusters(bamfile, num_tags):
    """Parses a BAM file, groups positions into clusters according to their
    barcodes, and returns the resulting structure.

    Each BAM record must have the barcode stored in the query name like so:

    ORIGINAL_READ_NAME::[Tag1][Tag2][Tag3]

    The tags should be enclosed in brackets and separated from
    the original read name with a double-colon.
    """
    #strip RPM DPM from barcode
    #TODO add file name as an additional barcode
    
    clusters = Clusters()
    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')

    for bam in bamfile:
        #get sample name from bamfile
        file_name = os.path.basename(bam)
        sample_name = file_name.split('.')[0]
        try:
            with pysam.AlignmentFile(bam, "rb") as f:
                for read in f.fetch(until_eof = True):
                    name = read.query_name
                    match = pattern.search(name)
                    barcode = list(match.groups())
                    strand = '+' if not read.is_reverse else '-'
                    position = Position('DNA', strand, read.reference_name,
                                        read.reference_start, read.reference_end)
                    barcode.append(sample_name)
                    barcode_str = ".".join(barcode)
                    clusters.add_position(barcode_str, position)
        except ValueError:
            print('BAM file provided is not a BAM or is empty!')

    return clusters



def write_clusters_to_file(clusters, outfile, unique=False):
    """Writes a Clusters object to a file"""

    count = 0
    with open(outfile, 'w') as f:

        if unique:
            for cluster_string in clusters.unique():
                f.write(cluster_string)
                f.write("\n")
                count += 1
        else:
            for cluster_string in clusters.to_strings():
                f.write(cluster_string)
                f.write("\n")
                count += 1
    print('Number of clusters written:',count)


def fastq_parse(fp):
    """
    Parse fastq file.
    """
    linecount = 0
    name, seq, thrd, qual = [None] * 4
    for line in fp:

        linecount += 1
        if linecount % 4 == 1:
            try:
                name = line.decode('UTF-8').rstrip()
            except AttributeError:
                name = line.rstrip()
            assert name.startswith('@'),\
                   "ERROR: The 1st line in fastq element does not start with '@'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 2:
            try:
                seq = line.decode('UTF-8').rstrip()
            except AttributeError:
                seq = line.rstrip()
        elif linecount % 4 == 3:
            try:
                thrd = line.decode('UTF-8').rstrip()
            except AttributeError:
                thrd = line.rstrip()
            assert thrd.startswith('+'),\
                   "ERROR: The 3st line in fastq element does not start with '+'.\n\
                   Please check FastQ file near line number %s" % (linecount)
        elif linecount % 4 == 0:
            try:
                qual = line.decode('UTF-8').rstrip()
            except AttributeError:
                qual = line.rstrip()
            assert len(seq) == len(qual),\
                    "ERROR: The length of Sequence and Quality aren't equal.\n\
                    Please check FastQ file near line number %s" % (linecount)

            yield name, seq, thrd, qual,
            name, seq, thrd, qual = [None] * 4



def file_open(filename):
    """
    Open as normal or as gzip
    Faster using zcat?
    """
    #does file exist?
    f = open(filename,'rb')
    if (f.read(2) == b'\x1f\x8b'): #compressed alsways start with these two bytes
        f.seek(0) #return to start of file
        return gzip.GzipFile(fileobj=f, mode='rb')
    else:
        f.seek(0)
        return f
