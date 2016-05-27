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
    master_script_obj.parsePhase(sys.argv[1], len(sys.argv))

    # ------------------------- MARKERS PHASE --------------------------
    master_script_obj.markersPhase()

    # -------------------------- GENOME PHASE ---------------------------
    # Construct genome objects, based on hom_fams and a list of all species.
    master_script_obj.genomePhase()

    # -------------------- FIND ADJACENCIES PHASE -----------------------
    master_script_obj.adjacenciesPhase()

    # -------------- ANCESTRAL GENOME CONSTRUCTION PHASE -----------------
    # Construct ancestral genome based on realizable intervals.

    master_script_obj.genomeConstructionPhase() # TODO: SPLIT


main()
