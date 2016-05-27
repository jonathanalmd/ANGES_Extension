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
    master_script_obj.setConfigParams(sys.argv[1], len(sys.argv))
    
    io_dict, markers_param_dict = master_script_obj.getIODictionary() # REMOVE LATER

    hom_fams_file = io_dict["homologous_families"] # REMOVE LATER
    pairs_file = io_dict["species_tree"] # REMOVE LATER
    output_dir = io_dict["output_directory"] # REMOVE LATER

    master_script_obj.setFileStreams()
    log, debug, hom_fams_stream, pairs_stream = master_script_obj.getFileStreams() # REMOVE LATER


    # PARSE PHASE
    # Parse the species pair file, put result in list.
    master_script_obj.parseSpeciesPairs()
    species_pairs = master_script_obj.getSpeciesPairs() # REMOVE LATER

    # Parse the hom fams file.
    master_script_obj.read_hom_families_file()
    hom_fams = master_script_obj.getHomFamList() # REMOVE LATER
 


    # MARKERS PHASE
    #Get all overlapped pairs
    #overlapped_pairs_list = process.getOverlappingPairs(hom_fams) 
  
    # Since the markers are all oriented, double them.
    master_script_obj.doubleMarkers()
    hom_fams = master_script_obj.getHomFamList() # REMOVE LATTER
 

    # GENOME PHASE
    # Construct genome objects, based on hom_fams and a list of all species.
    master_script_obj.constructGenomes()
    gens = master_script_obj.getGenomes() # REMOVE LATER



    # FIND ADJACENCIES PHASE
    # For each pair of species, compare the species to find adjacencies.
    master_script_obj.solveAdjacencies()
    adjacencies = master_script_obj.getAdjacencies() # REMOVE LATTER

    # Do the same for repeat spanning intervals
    master_script_obj.solveRSIs()
    RSIs = master_script_obj.getRSIs() # REMOVE LATTER


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





    # ANCESTRAL GENOME CONSTRUCTION PHASE
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
