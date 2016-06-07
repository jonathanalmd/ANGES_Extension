import markers

import sys
import bisect


# Genome data structure, kept ordered, allows in-order insertions.
# Represents the genome as the chromosomes dictionary of ordered
# lists of Marker objects.
class Genome:

    # Constructor
    # species: string, which species the genome is of.
    # chromosomes: a dictionary of list of Marker objects.
    def __init__( self, species ):
        self.species = species
        self.chromosomes = {}

    # Add locus in order to correct chromosome.
    # NOTE: this function does not set the 'index' field of inserted marker objects.
    def add_marker( self, id, locus, copy_number ):
        chrom = locus.chromosome
        if not chrom in self.chromosomes:
            self.chromosomes[chrom] = []

        # Insert the marker into correct chromosome in order.
        i = bisect.bisect_left( [x.locus.start for x in self.chromosomes[chrom]], locus.start )
        self.chromosomes[chrom].insert( i, Marker( id, locus, copy_number, 0 ) )

    # Find index of a marker in a chromosome, given its start location in the chromosome.
    # Returns negative index if not found.
    def find_marker( self, chrom, location ):
        current_chrom = self.chromosomes[chrom]
        i = bisect.bisect_left( [x.locus.start for x in current_chrom], location )
        if current_chrom[i].locus.start == location:
            return i
        else:
            return -sys.maxint-1


# Simplified data structure for a marker, with the following fields:
#    id: string - identity.
#    locus: Locus - position in genome.
#    copy_number: int - the number of markers with this ID in the ancestral genome.
#    index: int - index of this marker in chromosome in Genome object.
class Marker:
    def __init__( self, id, locus, copy_number, index ):
        self.id = id
        self.locus = locus
        self.copy_number = copy_number
        self.index = index

    def __repr__(self):
        #Return the direct representation of Genome_pairs_index class
        return "Marker id={} locus={} copy_number={} index={}".format(self.id, self.locus, self.copy_number, self.index)


# Get dictionary of sorted genomes
# Inpus:
# hom_fams: a list of HomFam objects
# species: a list of species whos genomes we want
# Output:
# genomes: a dictionary of genomes, with species names as keys
def get_genomes(hom_fams,species):
    # Dictionary of genomes
    genomes = {}
    # Create for each species in list
    for s in species:
        genomes[s] = Genome(s)
    # Add markers in place
    for hom_fam in hom_fams:
        for locus in hom_fam.loci:
            if locus.species in species:
                genomes[locus.species].add_marker( hom_fam.id, locus, hom_fam.copy_number )

    # Set marker indices
    for key1, genome in genomes.iteritems():
        for key2, chrom in genome.chromosomes.iteritems():
            for i in xrange( len( chrom ) ):
                chrom[i].index = i

    return genomes


# Function to return a list of the doubled counterparts of the given hom_fams.
# Arguments:
#   hom_fams: list of HomFam objects - the hom_fams to double.
# Output:
#   list of HomFam objects, doubled.
def double_oriented_markers( hom_fams ):
    doubled_hom_fams = []
    for hom_fam in hom_fams:
        # Get the new lists of loci:
        head_loci = []
        tail_loci = []
        for locus in hom_fam.loci:
            first =  markers.Locus( locus.species,
                                    locus.chromosome,
                                    locus.start,
                                    locus.start,
                                    0, "" )
            last = markers.Locus( locus.species,
                                  locus.chromosome,
                                  locus.end,
                                  locus.end,
                                  0, "" )
            if locus.orientation >= 0:
                head_loci.append( first )
                tail_loci.append( last )
            else:
                head_loci.append( last )
                tail_loci.append( first )

        doubled = markers.to_doubled( hom_fam.id )
        doubled_hom_fams.append( markers.HomFam( doubled[0],
                                                 head_loci,
                                                 hom_fam.copy_number,
                                                 hom_fam.comment ) )
        doubled_hom_fams.append( markers.HomFam( doubled[1],
                                                 tail_loci,
                                                 hom_fam.copy_number,
                                                 hom_fam.comment ) )
    return doubled_hom_fams
