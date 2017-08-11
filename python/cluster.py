import pysam
import re

class Position:
    """This class represents a genomic position.

    Methods:
    - to_string(): Returns a string representation of this position in the form
      "chrX:1000"
    """

    def __init__(self, chromosome, coordinate):
        self._chromosome = chromosome
        self._coordinate = coordinate

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        return (self._chromosome == other._chromosome and
                self._coordinate == other._coordinate)

    def __hash__(self):
        return hash((self._chromosome, self._coordinate))

    def to_string(self):
        return self._chromosome + ":" + str(self._coordinate)


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
      positions), and this empty cluster is returned.

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
        for barcode, cluster in self._clusters.iteritems():
            yield barcode + "\t" + cluster.to_string()

    def remove_cluster(self, barcode):
        del self._clusters[barcode]


def get_clusters(bamfile, num_tags):
    """Parses a BAM file, groups positions into clusters according to their
    barcodes, and returns the resulting structure.

    Each BAM record must have the barcode stored in the query name like so:

    ORIGINAL_READ_NAME::[Tag1][Tag2][Tag3]

    The tags should be enclosed in brackets and separated from
    the original read name with a double-colon.
    """

    clusters = Clusters()
    pattern = re.compile('::' + num_tags * '\[([a-zA-Z0-9_\-]+)\]')

    with pysam.AlignmentFile(bamfile, "rb") as f:

        for read in f.fetch(until_eof = True):
            position = Position(read.reference_name, read.reference_start)
            name = read.query_name
            match = pattern.search(name)
            barcode = ".".join(match.groups())
            clusters.add_position(barcode, position)

    return clusters


def write_clusters_to_file(clusters, outfile):
    """Writes a Clusters object to a file"""

    with open(outfile, 'w') as f:
        for cluster_string in clusters.to_strings():
            f.write(cluster_string)
            f.write("\n")
