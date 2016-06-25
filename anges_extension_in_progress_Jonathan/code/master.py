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
    if not master_script_obj.receivedAcsFile(): # if received ACS, skip this part
        master_script_obj.parse_markersPhase()
    # -------------------- FIND ADJACENCIES PHASE -----------------------
        master_script_obj.adjacenciesPhase()

    # -------------- ANCESTRAL GENOME CONSTRUCTION PHASE -----------------
    if master_script_obj.doC1PorNot():
        # Construct ancestral genome based on realizable adjacencies.
        print("Running C1P")
        master_script_obj.c1pPhase()
    else: # do MWM
        print("Running MWM")
        # ----------------- INTERVALS PHASE ------------------
        master_script_obj.intervalsPhase()
        # Construct ancestral genome based on realiazable intervals
        master_script_obj.genomeConstructionPhase() 

main()
