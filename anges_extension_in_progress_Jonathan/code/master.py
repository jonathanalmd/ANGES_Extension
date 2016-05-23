#!/usr/bin/env python

# Master script - produces the ancestral genome of given species.
# Takes four file names as input. How to execute:
# python master.py <HomFam data> <species pairs> <output dir>


import sys
import time

from itertools import combinations
from collections import defaultdict

from data_structures import markers
from data_structures import intervals
from data_structures import genomes
from data_structures import comparisons

import optimization
import assembly
    
def getOverlappedPairs(hom_fams):
    """
    getOverlappedPairs: receives a list of hom_fams and returns a list of overlapping pairs
    hom_fams - HomFam: list of objects from HomFam class
    """
#    import pdb; pdb.set_trace()
    loci_dict = defaultdict(list)
    for family in hom_fams:
#        loci_list.append(sorted(family.loci,key = lambda x_loci: x_loci.species))
        for loci in family.loci:
            loci_dict[loci.species].append(loci)
    print loci_dict
    overlapping_pairs_list = []
    for species,value in loci_dict.items():
        print (species,value)
#    for family in combinations(hom_fams, 2):
#        overlapping_pairs_aux = family[0].overlappingPairs(family[1])
#        if overlapping_pairs_aux:
#            overlapping_pairs_list.append(overlapping_pairs_aux)
#    return overlapping_pairs_list

def filterByID(hom_fams, ids):
    """ 
    filterByID: Receives a list of IDs and remove them from the main markers list
    Returns a list 
    hom_fams - HomFam: list of HomFam
    ids - int: list of IDs to be removed from the hom_fams list
    """
    return filter(lambda fam: int(fam.id) not in ids, hom_fams)

def filterByCopyNumber(hom_fams, threshold):
    """
    filterByCopyNumber: receives a threshold and filter the main markers list by comparing the 
    copy_number with the threshold (if copy_number > threshold, remove from the list).
    Returns a list without the filtered markers.
    hom_fams - HomFam: List of markers
    theshold - int: filter using the threshold (filter if the copy_number is > threshold)
    """
    return filter(lambda fam: fam.copy_number <= threshold, hom_fams)

def main():
    # Parse arguments:
    if len( sys.argv ) != 4:
        print ( "%s  ERROR (master.py) - script called with incorrect number "
                "of arguments." %strtime() )
        sys.exit()
    hom_fams_file = sys.argv[1]
    pairs_file = sys.argv[2]
    output_dir = sys.argv[3]
    try:
        log = open( output_dir + "/log", 'w' )
    except IOError:
        print ( "%s  ERROR (master.py) - could not open log file: %s"
                %( strtime(), output_dir + "/log" ) )
        sys.exit()

    debug = None
    if __debug__:
        try:
            debug = open( output_dir + "/debug", 'w' )
        except IOError:
            log.write( "ERROR (master.py) - could not open debug file %s"
                       % (output_dir + "/debug") )

    # Parse the species pair file, put result in list.
    species_pairs = []
    try:
        pairs_stream = open( pairs_file, 'r' )
    except IOError:
        log.write( "%s  ERROR (master.py) - could not open pairs file: %s\n"
                   %( strtime(), pairs_file ) )
        sys.exit()
    for pair in pairs_stream:
        # Assume format of "<species1> <species2>", respect comments
        pair = pair.strip()
        pair = pair.split()
        if pair[0][0] != "#" and len( pair ) == 2:
            species_pairs.append( pair )
    pairs_stream.close()
    log.write( "%s  Read %s species pairs.\n"
               %( strtime(), len( species_pairs ) ) )
    log.flush()

    # Parse the hom fams file.
    hom_fams = markers.read_hom_families_file( hom_fams_file )
    log.write( "%s  Read homologous families from file.\n"
               %( strtime() ) )
    log.flush()
    
    #Get all overlapped pairs
    overlapped_pairs_list = getOverlappedPairs(hom_fams)
    print overlapped_pairs_list
        
    # Since the markers are all oriented, double them.
    hom_fams = genomes.double_oriented_markers( hom_fams )
    # Construct genome objects, based on hom_fams and a list of all species.
    gens = genomes.get_genomes(
        hom_fams,
        list( set( species for pair in species_pairs for species in pair ) ),
        )
    log.write( "%s  Constructed genomes of %s species.\n"
               %( strtime(), len( gens ) ) )
    log.flush()

    # For each pair of species, compare the species to find adjacencies.
    adjacencies = intervals.IntervalDict()
    for pair in species_pairs:
        new_adjacencies = comparisons.find_adjacencies( gens[ pair[0] ],
                                                        gens[ pair[1] ] )
        comparisons.add_intervals( adjacencies, new_adjacencies )
    comparisons.set_interval_weights( adjacencies )
    log.write( "%s  Found %s adjacencies with total weight of %s.\n"
               %( strtime(),
                  len( adjacencies ),
                  adjacencies.total_weight ) )
    log.flush()
    write_intervals( adjacencies, output_dir + "/adjacencies", log )

    # Do the same for repeat spanning intervals
    RSIs = intervals.IntervalDict()
    for pair in species_pairs:
        new_RSIs = comparisons.find_RSIs( gens[ pair[0] ], gens[ pair[1] ] )
        comparisons.add_intervals( RSIs, new_RSIs )
    comparisons.set_interval_weights( RSIs )
    log.write( "%s  Found %s repeat spanning intervals with total weight of"
    " %s.\n" %( strtime(), len( RSIs ), RSIs.total_weight ) )
    log.flush()
    write_intervals( RSIs, output_dir + "/RSIs", log )

    # Select maximal subsets of adjacencies that are realizable.
    realizable_adjacencies = optimization.opt_adjacencies(
        hom_fams, adjacencies
        )

    write_intervals( realizable_adjacencies,
                     output_dir + "/realizable_adjacencies",
                     log )
    log.write( "%s  Found %s realizable adjacencies with total weight of %s.\n"
               %( strtime(),
                  len( realizable_adjacencies ),
                  realizable_adjacencies.total_weight ) )
    log.flush()

    # Keep track of adjacencies that have been discarded.
    discarded_adjacencies = intervals.IntervalDict()
    for adj in adjacencies.itervalues():
        if not adj.marker_ids in realizable_adjacencies:
            discarded_adjacencies.add( adj )
    write_intervals( discarded_adjacencies,
                     output_dir + "/discarded_adjacencies",
                     log )

    # Select maximal subsets of RSIs that are realizable.
    realizable_RSIs = optimization.opt_RSIs_greedy(
        hom_fams,
        realizable_adjacencies,
        RSIs,
        "mixed",
        debug,
        )
    write_intervals( realizable_RSIs, output_dir + "/realizable_RSIs", log )
    log.write( "%s  Found %s realizable repeat spanning intervals with total "
               "weight of %s.\n"
               %( strtime(),
                  len( realizable_RSIs ),
                  realizable_RSIs.total_weight ) )
    log.flush()

    # Keep track of RSIs that have been discarded.
    discarded_RSIs = intervals.IntervalDict()
    for RSI in RSIs.itervalues():
        if not RSI.marker_ids in realizable_RSIs:
            discarded_RSIs.add( RSI )
    write_intervals( discarded_RSIs,
                     output_dir + "/discarded_RSIs",
                     log )

    # Construct ancestral genome based on realizable intervals.
    ancestor_name = "ANCESTOR"
    ancestor_hom_fams = assembly.assemble(
        hom_fams,
        realizable_adjacencies,
        realizable_RSIs,
        ancestor_name,
        )
    markers.write_hom_families_file(
        output_dir + "/ancestor_hom_fams",
        ancestor_hom_fams,
        )
    # To order the hom_fams in chromosomes, create a Genome object with
    # the new hom_fams.
    ancestor_genomes = genomes.get_genomes(
        ancestor_hom_fams,
        [ ancestor_name ]
        )
    ancestor_genome = next( ancestor_genomes.itervalues() )
    try:
        genome_output = open( output_dir + "/ancestor_genome", 'w' )
        genome_output.write( ">" + ancestor_name + "\n" )
        for chrom_id, chrom in ancestor_genome.chromosomes.iteritems():
            genome_output.write( "#" + chrom_id + "\n" )
            for marker in chrom:
                if marker.locus.orientation > 0:
                    orient = "+"
                elif marker.locus.orientation < 0:
                    orient = "-"
                else:
                    orient = "x"
                genome_output.write( marker.id + " " + orient + "\n" )
    except IOError:
        log.write( "%s  ERROR (master.py) - could not write ancestor genome to " "file: %s\n" %( strtime(), output_dir + "/ancestor_genome" ) )
        sys.exit()
    log.write( "%s  Assembled the ancestral genome, found a total of %s CARs.\n"
               %( strtime(), len( ancestor_genome.chromosomes ) ) )
    log.write( "%s  Done.\n"%strtime() )


# Function to format time to string
def strtime():
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


# Function to write intervals to file
# Arguments:
#   ints: IntervalDict
#   f: the file to write to
#   log: log stream to write errors to
def write_intervals( ints, f, log ):
    try:
        output = open( f, 'w' )
        for interval in ints.itervalues():
            output.write( str( interval ) + "\n" )
    except IOError:
        log.write( "%s  ERROR (master.py) - could not write intervals to"
                   "file: %s\n" %( strtime(), f ) )
        sys.exit()



main()
