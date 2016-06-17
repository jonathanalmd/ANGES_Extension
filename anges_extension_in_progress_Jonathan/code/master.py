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
    master_script_obj.setConfigParams(sys.argv[1], len(sys.argv))


    # ------------------------ PARSE_MARKERS PHASE -------------------------------
    if not master_script_obj.receivedAcsFile():
        master_script_obj.parse_markersPhase()
    # -------------------- FIND ADJACENCIES PHASE -----------------------
        master_script_obj.adjacenciesPhase()

    # -------------- ANCESTRAL GENOME CONSTRUCTION PHASE -----------------
    # Construct ancestral genome based on realizable intervals.
    if not master_script_obj.doC1PorNot():
        master_script_obj.c1pPhase()
    else: # do MWM
        master_script_obj.genomeConstructionPhase() 

main()
