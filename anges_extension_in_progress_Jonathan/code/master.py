#!/usr/bin/env python

# Master script - produces the ancestral genome of given species.
# Takes four file names as input. How to execute:
# python master.py <HomFam data> <species pairs> <output dir>
import sys

from data_structures import process

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

    # ------------------------ PARSE_MARKERS PHASE -------------------------------
    master_script_obj.parse_markersPhase(sys.argv[1], len(sys.argv))

    # -------------------- FIND ADJACENCIES PHASE -----------------------
    # Construct genome objects, based on hom_fams and a list of all species and deal with adjacencies
    master_script_obj.adjacenciesPhase()

    # -------------- ANCESTRAL GENOME CONSTRUCTION PHASE -----------------
    # Construct ancestral genome based on realizable intervals.
    master_script_obj.genomeConstructionPhase() # TODO: SPLIT

    master_script_obj.c1pPhase()


main()
