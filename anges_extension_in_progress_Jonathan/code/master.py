#!/usr/bin/env python

# Master script - produces the ancestral genome of given species.
# Takes four file names as input. How to execute:
# python master.py <HomFam data> <species pairs> <output dir>
import sys

from data_structures import markers
from data_structures import intervals
from data_structures import genomes
from data_structures import comparisons
from data_structures import process

import optimization
import assembly

# -------------------------------------------- MAIN ---------------------------------------------------
def main():
    # Parse arguments:
    #hom_fams_file, pairs_file, output_dir
   
    #sys.argv[1] = ../data/configuration_file 

    master_script_obj = process.MasterScript()

    master_script_obj.readConfigFile(sys.argv[1], len(sys.argv))
    
    io_dict, markers_param_dict = master_script_obj.getIODictionary()

    hom_fams_file = io_dict["homologous_families"]
    pairs_file = io_dict["species_tree"]
    output_dir = io_dict["output_directory"]

    master_script_obj.setFileStreams()
    log, debug, hom_fams_stream, pairs_stream = master_script_obj.getFileStreams()

    # Parse the species pair file, put result in list.
    species_pairs = []
    for pair in pairs_stream:
        # Assume format of "<species1> <species2>", respect comments
        pair = pair.strip()
        pair = pair.split()
        if pair[0][0] != "#" and len( pair ) == 2:
            species_pairs.append( pair )
    pairs_stream.close()
    log.write( "%s  Read %s species pairs.\n"
               %( process.strtime(), len( species_pairs ) ) )
    log.flush()



    # Parse the hom fams file.
    hom_fams = markers.read_hom_families_file( hom_fams_file )
    log.write( "%s  Read homologous families from file.\n"
               %( process.strtime() ) )
    log.flush()



    
    #Get all overlapped pairs
    overlapped_pairs_list = process.getOverlappingPairs(hom_fams)
  
    # Since the markers are all oriented, double them.
    hom_fams = genomes.double_oriented_markers( hom_fams )
 
    # Construct genome objects, based on hom_fams and a list of all species.
    gens = genomes.get_genomes(
        hom_fams,
        list( set( species for pair in species_pairs for species in pair ) ),
        )
    log.write( "%s  Constructed genomes of %s species.\n"
               %( process.strtime(), len( gens ) ) )
    log.flush()

    # For each pair of species, compare the species to find adjacencies.
    adjacencies = intervals.IntervalDict()
    for pair in species_pairs:
        new_adjacencies = comparisons.find_adjacencies( gens[ pair[0] ],
                                                        gens[ pair[1] ] )
        comparisons.add_intervals( adjacencies, new_adjacencies )
    comparisons.set_interval_weights( adjacencies )
    log.write( "%s  Found %s adjacencies with total weight of %s.\n"
               %( process.strtime(),
                  len( adjacencies ),
                  adjacencies.total_weight ) )
    log.flush()
    process.write_intervals( adjacencies, output_dir + "/adjacencies", log )

    # Do the same for repeat spanning intervals
    RSIs = intervals.IntervalDict()
    for pair in species_pairs:
        new_RSIs = comparisons.find_RSIs( gens[ pair[0] ], gens[ pair[1] ] )
        comparisons.add_intervals( RSIs, new_RSIs )
    comparisons.set_interval_weights( RSIs )
    log.write( "%s  Found %s repeat spanning intervals with total weight of"
    " %s.\n" %( process.strtime(), len( RSIs ), RSIs.total_weight ) )
    log.flush()
    process.write_intervals( RSIs, output_dir + "/RSIs", log )

    # Select maximal subsets of adjacencies that are realizable.
    realizable_adjacencies = optimization.opt_adjacencies(
        hom_fams, adjacencies
        )

    process.write_intervals( realizable_adjacencies,
                     output_dir + "/realizable_adjacencies",
                     log )
    log.write( "%s  Found %s realizable adjacencies with total weight of %s.\n"
               %( process.strtime(),
                  len( realizable_adjacencies ),
                  realizable_adjacencies.total_weight ) )
    log.flush()

    # Keep track of adjacencies that have been discarded.
    discarded_adjacencies = intervals.IntervalDict()
    for adj in adjacencies.itervalues():
        if not adj.marker_ids in realizable_adjacencies:
            discarded_adjacencies.add( adj )
    process.write_intervals( discarded_adjacencies,
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
    process.write_intervals( realizable_RSIs, output_dir + "/realizable_RSIs", log )
    log.write( "%s  Found %s realizable repeat spanning intervals with total "
               "weight of %s.\n"
               %( process.strtime(),
                  len( realizable_RSIs ),
                  realizable_RSIs.total_weight ) )
    log.flush()

    # Keep track of RSIs that have been discarded.
    discarded_RSIs = intervals.IntervalDict()
    for RSI in RSIs.itervalues():
        if not RSI.marker_ids in realizable_RSIs:
            discarded_RSIs.add( RSI )
    process.write_intervals( discarded_RSIs,
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
        log.write( "%s  ERROR (master.py) - could not write ancestor genome to " "file: %s\n" %( process.strtime(), output_dir + "/ancestor_genome" ) )
        sys.exit()
    log.write( "%s  Assembled the ancestral genome, found a total of %s CARs.\n"
               %( process.strtime(), len( ancestor_genome.chromosomes ) ) )
    log.write( "%s  Done.\n"%process.strtime() )




main()
