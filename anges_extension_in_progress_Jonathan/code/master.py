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
    """
    5 main phases:
        parse phase
        markers phase
        genome phase
        adjacencies phase
        genome construction phase
    """
    master_script_obj = process.MasterScript()

    # ------------------------ PARSE PHASE -------------------------------
    # Parse arguments: 
    #sys.argv[1] = ../data/configuration_file 
    master_script_obj.setConfigParams(sys.argv[1], len(sys.argv))
    master_script_obj.setFileStreams()

    # Parse the species pair file, put result in list.
    master_script_obj.parseSpeciesPairs()

    # Parse the hom fams file.
    master_script_obj.read_hom_families_file()
 

    # ------------------------- MARKERS PHASE --------------------------
    #Get all overlapped pairs
    #overlapped_pairs_list = process.getOverlappingPairs(hom_fams) 
  
    # Since the markers are all oriented, double them.
    master_script_obj.doubleMarkers()


    # -------------------------- GENOME PHASE ---------------------------
    # Construct genome objects, based on hom_fams and a list of all species.
    master_script_obj.constructGenomes()


    # -------------------- FIND ADJACENCIES PHASE -----------------------
    # For each pair of species, compare the species to find adjacencies.
    master_script_obj.solveAdjacencies()

    # Do the same for repeat spanning intervals
    master_script_obj.solveRSIs()

    # Select maximal subsets of adjacencies that are realizable.
    master_script_obj.selectMaxAdjacencies()

    # Keep track of adjacencies that have been discarded.
    master_script_obj.trackDiscardedAdjacencies()

    # Select maximal subsets of RSIs that are realizable.
    master_script_obj.selectMaxRSIs()

    # Keep track of RSIs that have been discarded.
    master_script_obj.trackDiscardedRSIs()


    # -------------- ANCESTRAL GENOME CONSTRUCTION PHASE -----------------
    # Construct ancestral genome based on realizable intervals.

    master_script_obj.genomeConstructionPhase() # TODO: SPLIT


main()
