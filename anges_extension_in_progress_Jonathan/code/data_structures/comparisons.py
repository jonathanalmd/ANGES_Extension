import genomes
import markers
import intervals
import itertools
import collections

class Genome_pairs_index:

    # Constructor: minimalist, populate this datastructure with the
    # 'index_...' functions
    def __init__( self ):
        self.pairs = {}
        self.species = ''

    def __repr__(self):
        #Return the direct representation of Genome_pairs class
        return "Genome_pairs_index {}".format(self.species)


    # Function to process a genome and get a dictionary of adjacent pairs of markers.
    # Arguments:
    #    genome: a Genome object
    # Post-condition:
    #    'pairs' will have the following structure:
    #    dictionary, keyed with marker IDs, each element is a:
    #     -> dictionary, keyed with marker IDs, each element is a:
    #         -> list, where each element is a:
    #             -> tuple, containing the pair of markers corresponding to the
    #                dictionary keys for this list, followed by the locus of the
    #                pair.
    #
    #    The orientations of the loci indicate which marker comes first in the pair.
    #    Example: pair ('a', 'b') will appear twice: under key 'a' with positive
    #    orientation, and under key 'b' with negative orientation.
    def index_pairs( self, genome ):
        self.species = genome.species

        for key,chrom in genome.chromosomes.iteritems():
            # Is this chromosome long enough for an adjacency?
            if len( chrom ) < 2: continue

            # Loop over the chromosome
            for i in xrange( len(chrom) ):
                # Add new entry to pairs dictionary
                if not chrom[i].id in self.pairs:
                    self.pairs[ chrom[i].id ] = {}

                # Index the pair with marker to the right, if not head/tail pair.
                if i < len(chrom) - 1 and chrom[i].id != markers.sibling( chrom[i+1].id ):
                    if not chrom[i+1].id in self.pairs[ chrom[i].id ]:
                        self.pairs[ chrom[i].id ][ chrom[i+1].id ] = []
                    self.pairs[ chrom[i].id ][ chrom[i+1].id ].append(
                            ( markers.Locus( genome.species,
                                             key,
                                             chrom[i].locus.start,
                                             chrom[i+1].locus.end,
                                             1, '' ),
                              chrom[i],
                              chrom[i+1] ) )
                # Index the pair with marker to the left
                if i > 0 and chrom[i].id != markers.sibling( chrom[i-1].id ):
                    if not chrom[i-1].id in self.pairs[ chrom[i].id ]:
                        self.pairs[ chrom[i].id ][ chrom[i-1].id ] = []
                    self.pairs[ chrom[i].id ][ chrom[i-1].id ].append(
                            ( markers.Locus( genome.species,
                                             key,
                                             chrom[i-1].locus.start,
                                             chrom[i].locus.end,
                                             -1, '' ),
                              chrom[i-1],
                              chrom[i] ) )


    # Function to search the indexed genome for the give pair of markers.
    # Arguments:
    #    marker1, marker2: the ordered pair of markers to search for.
    # Output:
    #    A list of tuples, see format in comments for indexing function.
    def query( self, marker_id1, marker_id2 ):
        answer = []
        # Can we find the pair?
        if marker_id1 in self.pairs and marker_id2 in self.pairs[ marker_id1 ]:
            answer = self.pairs[ marker_id1 ][ marker_id2 ]

        return answer


# Function to find the adjacencies between two genomes.
# Arguments:
#    genome1, genome2: Genome objects
# Output:
#    An IntervalDict of adjacencies
def find_adjacencies( genome1, genome2, all_match ):
    adjs, _ = find_intervals( genome1, genome2, all_match, strip=False)
    return adjs

# Function to find the RSIs between two genomes.
# Arguments:
#    genome1, genome2: Genome objects
# Output:
#    An IntervalDict of RSIs
def find_RSIs( genome1, genome2, all_match):
    _, RSIs = find_intervals( genome1, genome2, all_match, strip=True )
    return RSIs

# Function to find the adjacencies and repeat spanning intervals between two
# genomes.
# Arguments:
#    genome1, genome2: Genome objects
#    strip: boolean - whether to strip given genomes or not.
# Output:
#    adjs, RSIs: IntervalDict objects
def find_intervals( genome1, genome2, all_match, strip):
    # NOTE: There might be trouble with RSIs that are the delimited at front
    #       and back by the same marker.

    # Only strip genomes (to find RSIs) if instructed in arguments
    stripped1 = genome1
    stripped2 = genome2
    if strip:
        stripped1 = strip_genome_unique( stripped1 )
        stripped2 = strip_genome_unique( stripped2 )

    # Index one genome, loop over the other.
    index = Genome_pairs_index()
    index.index_pairs( stripped1 )

    ints = intervals.IntervalDict()
    for key,chrom in stripped2.chromosomes.iteritems():
        for i in xrange( len(chrom) - 1 ):
            pair2 = ( chrom[i].id, chrom[i+1].id )
            if markers.are_siblings( *pair2 ):
                continue
            full_ids2 = [ marker.id for marker in
                         genome2.chromosomes[ key ] \
                             [ chrom[i].index : chrom[i+1].index+1 ] ]
            pairs = index.query( *pair2 )
            # Find all pairs in the index with a full match (or all match)
            matches = []
            match_ids = []
            for pair1 in pairs:
                # do we have a match between pair1 and pair2?
                
                full_ids1 = [ marker.id for marker in
                              genome1.chromosomes[ pair1[0].chromosome ] \
                                  [ pair1[1].index : pair1[2].index+1 ] ]

                cmp_result, _ = compare_marker_intervals(
                    full_ids1,
                    full_ids2, 
                    all_match
                    )
                if cmp_result:
                    matches.append( pair1 )
                    match_ids = full_ids1

            if matches:
                # Add a new interval to ints if not already there, and add the
                # loci from matches to it. If the interval is already there, it
                # must already have the loci from matches.
                if not full_ids2 in ints:
                    ints.add( intervals.Interval(
                        id=''.join(full_ids2),
                        marker_ids=full_ids2,
                        loci=[l for l, _, _ in matches],
                        order=intervals.Order( 1 ),
                        weight=1,
                        comment='',
                        ) )
                # Add the newly found locus in genome 2 to the interval.
                interval = ints[match_ids]
                locus2 = markers.Locus(
                    species=genome2.species,
                    chromosome=key,
                    start=chrom[i].locus.start,
                    end=chrom[i+1].locus.end,
                    orientation=1,
                    comment='',
                    )
                if interval.marker_ids == list( reversed( full_ids2 ) ):
                    locus2.orientation *= -1
                interval.loci.append( locus2 )

    # Filter the ints to distinguish repeat spanning intervals.
    adjs = intervals.IntervalDict()
    RSIs = intervals.IntervalDict()
    for interval in ints.itervalues():
        if len( interval.marker_ids ) == 2:
            adjs.add( interval )
        else:
            RSIs.add ( interval )

    return adjs, RSIs



# Function to add new adjacencies to an existing database of adjacencies,
# updating the weights of the adjacencies along the way.
# Arguments:
#    adjacencies: IntervalDict - database to update.
#    new_adjacencies: IntervalDict, intervals to add to existing database.
def add_intervals( ints, new_ints ):
    for interval in new_ints.itervalues():
        if interval.marker_ids not in ints:
            ints.add( interval )
        else:
            found = ints[interval.marker_ids]
            new_loci = interval.loci
            if found.marker_ids == list( reversed( interval.marker_ids ) ):
                for locus in new_loci:
                    locus.orientation *= -1
            found.loci += new_loci
            # Remove duplicate loci
            found.loci = sorted( found.loci )
            found.loci = [ locus for locus, _ in \
                           itertools.groupby( found.loci ) ]

# Function to set the weights of adjacencies.
# Currently just counts the number of occurences of the adjacencies.
# Arguments:
#   adjacencies: nested dict of Interval objects who's weights to set.
def set_interval_weights( ints ):
    for interval in ints.itervalues():
        interval.weight = len( interval.loci )


# Auxilary functions:

# Function to return a genome with only unique markers
# Arguments:
#     genome: Genome
# Output:
#     A new genome that is a copy of 'genome', but without markers with
#     a 'copy_number' greater than 1
def strip_genome_unique( genome ):
    stripped_genome = genomes.Genome( genome.species )
    for key,chrom in genome.chromosomes.iteritems():
        stripped_genome.chromosomes[ key ] = [ marker for marker in chrom if marker.copy_number == 1 ]

    return stripped_genome

# Function to check if two lists of doubled marker IDs are equivalent.
# This function checks for positive and negative orientation.
# Arguments:
#     markers1, markers2: the lists of marker IDs
# Output:
#     tuple of boolean (true of IDs are equivalent) and orientation, for example (True, -1)
def compare_marker_intervals( marker_ids1, marker_ids2, all_match):
    if all_match:
        # if set(marker_ids1) == set(marker_ids2):
        # collections.Counter(marker_ids1) == collections.Counter(marker_ids2)
        if sorted(marker_ids1) == sorted(marker_ids2):
            return (True, 1) # 1, -1 or 0?? 
        else:
            return (False, 0)
    else:
        if marker_ids1 == marker_ids2:
            return ( True, 1 )
        elif marker_ids1 == list( reversed( marker_ids2 ) ):
            return ( True, -1 )
        else:
            return ( False, 0 )

