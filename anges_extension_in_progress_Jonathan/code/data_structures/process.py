import sys
import time

from collections import defaultdict

from data_structures import markers
from data_structures import intervals
from data_structures import genomes
from data_structures import comparisons

import optimization
import assembly



import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

import parameters

import subprocess
import copy
import shutil
import threading
import time

global log_file

def print2(string):    

    try:
        print string
        log_file.write(string + '\n')
#        sys.__stdout__.write(string)
#        sys.__stdout__.flush()
    except:
        pass
    #endtry

    #endig
#enddef

def callprocess(params):
    process = subprocess.Popen(params, bufsize=0, executable=None, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=None, close_fds=False, shell=False, cwd=None, env=None, universal_newlines=False, startupinfo=None, creationflags=0)
    
    while process.poll() == None:
        s = process.stdout.read()
    
        if s and s != "":
            print2(s)
        #endif
        
        s = process.stderr.read()
    
        if s and s != "":
            print2(s)
        #endif
    
        time.sleep(.05)
    #endif
    
    s = process.stdout.read()

    if s and s != "":
        print2(s)
    #endif
    
    s = process.stderr.read()

    if s and s != "":
        print2(s)
    #endif
#enddef


def do_c1p(code_dir, params, c1p_dir, cars_dir, output_prefix, wacs, quiet, suffix, message, make_C1P, compute_PQRtree, m = False):
    acs_c1p          = c1p_dir + "/" + output_prefix + "ACS_C1P_" + suffix    # C1P matrix file
    acs_discarded = c1p_dir + "/" + output_prefix + "ACS_DISC_" + suffix    # matrix of removed rows file
    pq_tree = cars_dir + "/" + output_prefix + "PQTREE_" + suffix # PQ-tree file    
    if params.markers_doubled:
        pq_tree_doubled = cars_dir + "/" + output_prefix + "PQTREE_DOUBLED_" + suffix # PQ-tree file    
    #endif
    if not quiet:
        print2("----> Making (weighted)ACS matrix C1P (" + message + "): " + wacs + " " + acs_discarded)
    #endif

    if m:
        callprocess(["python", code_dir + make_C1P, wacs, 'max', acs_c1p, acs_discarded])
    else:
        callprocess(["python", code_dir + make_C1P, wacs, acs_c1p, acs_discarded])
    #endif

    if not quiet:
        print2("----> Creating a PQ-tree: " + pq_tree)
    #endif

    if params.markers_doubled:
        callprocess(["python", code_dir + compute_PQRtree, acs_c1p, pq_tree_doubled, params.output_ancestor])
        if not quiet:
            print2("----> Halving PQ-tree columns") 
        #endif
        callprocess(["python", code_dir + "/C1P/C1P_halve_PQRtree.py", pq_tree_doubled, pq_tree])
    else:
        callprocess(["python", code_dir + compute_PQRtree, acs_c1p, pq_tree, params.output_ancestor])
    #endif
#enddef



    
def strtime():
    """
    Function to format time to string
    """
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
#enddef

class MasterMarkers:
    """
    Parse and solve markers info
    """
    def __init___(self):
        self.hom_fams_file_stream = None
        self.pairs_file_stream = None
    #enddef

    def setInputStreams(self, species_tree_dir, hom_fam_dir):
        try:
            self.pairs_file_stream = open(species_tree_dir, 'r')
        except IOError:
            log.write( "{}  ERROR (master.py) - could not open pairs file: {}\n"
                       .format(strtime(), self.io_dict["species_tree"]))
            sys.exit()
        
        try:
            self.hom_fams_file_stream = open(hom_fam_dir, 'r')
        except IOError:
            log.write( "{}  ERROR (master.py) - could not open homologous families file: {}\n"
                       .format(strtime(), self.io_dict["homologous_families"]))
            sys.exit()
    #enddef

    def parseSpeciesPairs(self, species_set, species_tree_dir, log):
        """
        Populating the species_pairs list
        """
        line_number = 1
        species_pairs = []
        for pair in self.pairs_file_stream:
        # Assume format of "<species1> <species2>", respect comments
            pair = pair.strip()
            pair = pair.split()
            if pair[0][0] != "#":
                if len( pair ) == 2:
                    if pair[0] not in species_set:
                        print("Semantic Error at line {} (file: {}): '{}'\n\tInvalid species '{}'. Not listed in the Homologous Families file\n"
                                .format(line_number,species_tree_dir, pair[0]+' '+pair[1],pair[0]))
                    elif pair[1] not in species_set:
                        print("Semantic Error at line {} (file: {}): '{}'\n\tInvalid species '{}'. Not listed in the Homologous Families file\n"
                                .format(line_number,species_tree_dir, pair[0]+' '+pair[1],pair[1]))
                    else:
                        species_pairs.append( pair )
                else:
                    print("Syntatic Error at line {} (file: {})\n\tInvalid number of species. Each line must have only one pair of species\n"
                                .format(line_number,species_tree_dir))
            line_number = line_number + 1
        self.pairs_file_stream.close()
        log.write( "{}  Read {} species pairs.\n"
                   .format(strtime(), len(species_pairs)))
        log.flush()

        return species_pairs
    #enddef

    # reads hom. families from a file
    # file_name - str: the name of the file to read from
    # hom_fam_list - list of HomFam: the list to add to (Default = [])
    # Return - list of HomFam: the list of hom. familes read
    def parseHomFamilies(self, hom_fam_dir, log):
        """
        Populates hom_fam_list
        """
        hom_fam_list = []
        line_number = 1
        line = self.hom_fams_file_stream.readline()
        while len(line) > 0:
            trunc_line = line.strip()

            if len(trunc_line) > 0:
                if trunc_line[0] == '>':
                    # read first hom. family
                    hom_fam, line = markers.HomFam.from_file(self.hom_fams_file_stream, trunc_line, line_number, hom_fam_dir)
                    if hom_fam != None:
                        hom_fam_list.append(hom_fam)
                    #endif

                    # read the rest of the hom. families
                    while len(line) > 0:
                        hom_fam, line = markers.HomFam.from_file(self.hom_fams_file_stream, line, line_number, hom_fam_dir)
                        if hom_fam != None:
                            hom_fam_list.append(hom_fam)
                        #endif
                    #endwhile
                else:
                    line = self.hom_fams_file_stream.readline()
                    line_number = line_number + 1
                #endif
            else:
                line = self.hom_fams_file_stream.readline()
                line_number = line_number + 1
            #endif
        #endif

        self.hom_fams_file_stream.close()

        log.write( "{}  Read homologous families from file.\n"
               .format( strtime() ) )
        log.flush()

        return hom_fam_list
    #enddef

    def doubleMarkers(self, hom_fam_list):
        return genomes.double_oriented_markers(hom_fam_list)
    #enddef

    def getOverlappingPairs(self, hom_fam_list, log):
        """
        getOverlappingPairs: receives a list of hom_fams and returns a list of overlapping pairs (each pair is a tuple containing two Locus)
        hom_fams - HomFam: list of objects from HomFam class
        """
        #import pdb; pdb.set_trace()
        loci_dict = defaultdict(list)
        overlapping_pairs_list = []

        for hom_fam_index, marker_family in enumerate(hom_fam_list):
            for locus_index, locus in enumerate(marker_family.loci):
               loci_dict[locus.species].append((hom_fam_index, locus_index))
            #endfor
        #endfor
        for species, species_indexes in loci_dict.items():
            # species name and a list of tuples with (hom_fam_index, locus_index)
            # hom_fam_index => access to hom_fam ID and loci list
            #print (species, species_indexes)
            i = 1
            locus1 = hom_fam_list[species_indexes[0][0]].loci[species_indexes[0][1]]
            while i < len(species_indexes):
                locus2 = hom_fam_list[species_indexes[i][0]].loci[species_indexes[i][1]]
                #print (locus1, locus2)
                is_overlapping_pair = locus1.overlappingPairs(locus2)
                if is_overlapping_pair:
                    overlapping_pairs_list.append(is_overlapping_pair)
                #endif
                i = i + 1
            #endwhile
        #endfor
        log.write("{}  {} overlapping pairs have been found.\n"
                    .format(strtime(), len(overlapping_pairs_list)))
        return overlapping_pairs_list
    #enddef

    def filterByID(self, id_list, hom_fam_list, log):
        """ 
        filterByID: Receives a list of IDs and remove them from the main markers list
        Returns a list 
        hom_fams - HomFam: list of HomFam
        ids - int: list of IDs to be removed from the hom_fams list
        """
        id_list = []
        filtered_list = filter(lambda fam: int(fam.id) not in id_list, hom_fam_list)
        for hom_fam in filtered_list:
            id_list.append(int(hom_fam.id))
        log.write( "{}  IDs {} filtred from homologous families list.\n"
                .format(strtime(), id_list))
        return filtered_list
    #enddef

    def filterByCopyNumber(self, copy_number_threshold, hom_fam_list, log):
        """
        filterByCopyNumber: receives a threshold and filter the main markers list by comparing the 
        copy_number with the threshold (if copy_number > threshold, remove from the list).
        Returns a list without the filtered markers.
        hom_fams - HomFam: List of markers
        theshold - int: filter using the threshold (filter if the copy_number is > threshold)
        """
        filtered_list = filter(lambda fam: fam.copy_number <= copy_number_threshold, hom_fam_list)
        filtered_quant = len(hom_fam_list) - len(filtered_list)
        log.write("{}  Filtered {} homologous families with Copy Number greater than {}.\n"
                    .format(strtime(), filtered_quant, copy_number_threshold ))
        return filtered_list
     #enddef
 
    def closePairsFile(self):
        self.pairs_file_stream.close()
    #enddef

    def closeHomFamsFile(self):
        self.hom_fams_file_stream.close()
    #enddef

class MasterAdjacencies:
    def __init__(self):
        self.adjacencies = intervals.IntervalDict()
        self.realizable_adjacencies = intervals.IntervalDict()
        self.discarded_adjacencies = intervals.IntervalDict()
    #enddef

    def getAdjacencies(self):
        return self.adjacencies
    #enddef

    def getRealizableAdjacencies(self):
        return self.realizable_adjacencies
    #enddef

    def getDiscardedAdjacencies(self):
        return self.discarded_adjacencies
    #enddef

    def solveAdjacencies(self, species_pairs, gens, output_directory, log, all_match):
        for pair in species_pairs:
            new_adjacencies = comparisons.find_adjacencies( gens[ pair[0] ],
                                                            gens[ pair[1] ], all_match)
            comparisons.add_intervals( self.adjacencies, new_adjacencies )
        comparisons.set_interval_weights( self.adjacencies )
        log.write( "{}  Found {} adjacencies with total weight of {}.\n"
                   .format( strtime(),
                      len( self.adjacencies ),
                      self.adjacencies.total_weight ) )
        log.flush()
        intervals.write_intervals(log, self.adjacencies, output_directory + "/adjacencies")
    #enddef

    def selectMaxAdjacencies(self, hom_fam_list, output_directory, log):
        self.realizable_adjacencies = optimization.opt_adjacencies(hom_fam_list, self.adjacencies)
        intervals.write_intervals(log, self.realizable_adjacencies, 
                         output_directory + "/realizable_adjacencies",
                        )
        log.write( "{}  Found {} realizable adjacencies with total weight of {}.\n"
               .format(strtime(),len(self.realizable_adjacencies), self.realizable_adjacencies.total_weight ) )
        log.flush()
    #enddef

    def trackDiscardedAdjacencies(self, output_directory, log):
        for adj in self.adjacencies.itervalues():
            if not adj.marker_ids in self.realizable_adjacencies:
                self.discarded_adjacencies.add( adj )
        intervals.write_intervals(log, self.discarded_adjacencies, output_directory + "/discarded_adjacencies")
    #enddef
#endclass

class MasterRSI:
    def __init__(self):
        self.RSIs = intervals.IntervalDict()
        self.realizable_RSIs = intervals.IntervalDict()
        self.discarded_RSIs = intervals.IntervalDict()
    #enddef

    def getRSIs(self):
        return self.RSIs
    #enddef

    def getRealizableRSIs(self):
        return self.realizable_RSIs
    #enddef

    def getDiscardedRSIs(self):
        return self.discarded_RSIs
    #enddef

    def solveRSIs(self, species_pairs, gens, output_directory, log, all_match):
        for pair in species_pairs:
            new_RSIs = comparisons.find_RSIs(gens[ pair[0] ], gens[ pair[1] ], all_match)
            comparisons.add_intervals( self.RSIs, new_RSIs )
        comparisons.set_interval_weights( self.RSIs )
        log.write( "{}  Found {} repeat spanning intervals with total weight of"
        " {}.\n" .format( strtime(), len( self.RSIs ), self.RSIs.total_weight ) )
        log.flush()
        intervals.write_intervals(log, self.RSIs, output_directory + "/RSIs")
    #enddef
    
    def selectMaxRSIs(self, hom_fam_list, realizable_adjacencies, output_directory, log, debug):
        self.realizable_RSIs = optimization.opt_RSIs_greedy(
            hom_fam_list,
            realizable_adjacencies, 
            self.RSIs,
            "mixed",
            debug,
            )
        intervals.write_intervals(log, self.realizable_RSIs, output_directory + "/realizable_RSIs")
        log.write( "{}  Found {} realizable repeat spanning intervals with total weight of {}.\n"
                   .format(
                    strtime(),
                    len(self.realizable_RSIs ),
                    self.realizable_RSIs.total_weight 
                    )
                )
        log.flush()
    #enddef
    
    def trackDiscardedRSIs(self, output_directory, log):
        for RSI in self.RSIs.itervalues():
            if not RSI.marker_ids in self.realizable_RSIs:
                self.discarded_RSIs.add( RSI )
        intervals.write_intervals(log, self.discarded_RSIs, output_directory + "/discarded_RSIs")
    #enddef
#endclass

class MasterGenConstruction:
    def __init__(self): #deal with adjacencies and genome construction
        self.RSI = MasterRSI()
        self.adj = MasterAdjacencies()
        self.ancestor_name = "ANCESTOR"
        self.ancestor_hom_fams = []
        self.gens = {}
    #enddef

    def getGenomes(self):
        return self.gens
    #enddef

    # writes hom. families to a file
    # file_name - str: the name of the file to write to
    # hom_fam_list - list of HomFam: the list to write (Default = [])
    def writeAncestorHomFams(self,file_name):
        file_stream = open(file_name, 'w')

        for hom_fam in self.ancestor_hom_fams:
            hom_fam.to_file( file_stream )
            file_stream.write("\n")
            file_stream.flush()
        #endif

        file_stream.close()
    #enddef

    def constructGenomes(self, species_pairs, hom_fam_list, log):
        self.gens = genomes.get_genomes(
            hom_fam_list,
            list( set( species for pair in species_pairs for species in pair ) ),
            )
        log.write( "{}  Constructed genomes of {} species.\n" 
                        .format( strtime(), len( self.gens ) ) )
        log.flush()
    #enddef

    def dealWithAdjPhase(self, species_pairs, hom_fam_list, output_directory, log, debug, all_match):
        # For each pair of species, compare the species to find adjacencies.
        self.adj.solveAdjacencies(species_pairs, self.gens, output_directory, log, all_match)

        # Do the same for repeat spanning intervals
        self.RSI.solveRSIs(species_pairs, self.gens, output_directory, log, all_match)

        # Select maximal subsets of adjacencies that are realizable.
        self.adj.selectMaxAdjacencies(hom_fam_list, output_directory, log)

        # Keep track of adjacencies that have been discarded.
        self.adj.trackDiscardedAdjacencies(output_directory, log)

        # Select maximal subsets of RSIs that are realizable.
        self.RSI.selectMaxRSIs(hom_fam_list, self.adj.realizable_adjacencies, output_directory, log, debug)

        # Keep track of RSIs that have been discarded.
        self.RSI.trackDiscardedRSIs(output_directory, log)
    #enddef

    def dealWithConstructionPhase(self, hom_fam_list, output_directory, log):
        self.ancestor_hom_fams = assembly.assemble(
            hom_fam_list,
            self.adj.realizable_adjacencies,
            self.RSI.realizable_RSIs,
            self.ancestor_name,
            )
        self.writeAncestorHomFams(
            output_directory + "/ancestor_hom_fams",
            )
        # To order the hom_fams in chromosomes, create a Genome object with
        # the new hom_fams.
        self.ancestor_genomes = genomes.get_genomes( 
            self.ancestor_hom_fams,
            [ self.ancestor_name ]
            )
        ancestor_genome = next( self.ancestor_genomes.itervalues() )
        try:
            genome_output = open( output_directory + "/ancestor_genome", 'w' )
            genome_output.write( ">" + self.ancestor_name + "\n" )
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
            loglog.write( "{}  ERROR (master.py) - could not write ancestor genome to " "file: {}\n"
                            .format(strtime(), output_directory + "/ancestor_genome" ) )
            sys.exit()

        log.write( "{}  Assembled the ancestral genome, found a total of {} CARs.\n"
                   .format(strtime(), len(ancestor_genome.chromosomes) ) )
        log.write( "{}  Done.\n".format(strtime()) )
    #enddef


#endclass

class MasterScript:
    """
    MasterScript Class: coordinates the main process using some aux classes
    A MasterScript object must have all the most important information regarding homologous families and spcecies pairs and this object will coordinate all phases 
    io_dict - {string : string}: dictionary that will be populated based on the configuration.conf file 
                                this dictionary will have all the IO directories
    markers_param_dict - {string : int} or {string : [int]}: also will be populated based on the configuration.conf file (parse_markersPhase)
                                this dictionary will have all the markers running parameters (different execution modes)
    log - filestream: stream for the log file
    debug - filestream: stream for the debug file
    species_pairs: list of species_pairs -> populated based on the species_pairs file, during the parse_markersPhase
    hom_fam_list: list of homologous families -> populated based on the hom_fams_file, during the parse_markersPhase
    overlapped_pairs_list: list of overlapping pairs, populated based on the hom_fam_list, during the parse_markersPhase
    genome_construction_obj - MasterGenConstruction class object -> this object deal with adjacenciesPhase and genomeConstructionPhase and 
                                                                    has information regarding gens, adjacencies, ancestral genomes (including methods to manipulate these information)  
    """
    def __init__(self):
        self.io_dict = {}
        self.markers_param_dict = {}
        self.run_param_dict = {}

        self.log = None
        self.debug = None
    
        self.species_pairs = []
        self.hom_fam_list = []
        self.species_set = set()
        self.overlapping_pairs_list = []
        self.genome_construction_obj = MasterGenConstruction()
    #enddef

    def setOutputStreams(self):
        try:
            self.log = open(self.io_dict["output_directory"] + "/log", 'w' )
        except IOError:
            print ( "{}  ERROR (master.py) - could not open log file: {}"
                    .format(strtime(), self.io_dict["output_directory"] + "/log" ) )
            sys.exit()
        if __debug__:
            try:
                self.debug = open(self.io_dict["output_directory"] + "/debug", 'w' )
            except IOError:
                log.write( "ERROR (master.py) - could not open debug file {}"
                           .format(self.io_dict["output_directory"] + "/debug"))
    #enddef

    def setConfigParams(self, config_file, len_arguments):
        """
        Function to read the configuration file and return all the interpreted information
        Arguments:
        config_file - string: configuration file directory
        len_arguments - int: arguments length

        Returns: homologous families file directory
                 species tree file directory 
                 output directory
                 markers parameters (markers_doubled - int, markers_unique - int, markers_universal - int,
                                     markers_overlap - int, filter_copy_number - int, filter_by_id - list of int)
                where homologous families, species tree and output are keys for a dictionary and their respective 
                directories are the values for these keys;
                markers parameters is a dictionary using markers_doubled, ..., filter_by_id as keys and the values
                are the parameter values of each parameter (informed in the configuration file)
        """

        if len_arguments != 2:
            print ( "{}  ERROR (master.py -> process.py) - script called with incorrect number "
                     "of arguments.".format(strtime()))
            sys.exit()
        #endif
        
        try:
            config = {}
            execfile(config_file, config)
            #collect the information from config file
            self.io_dict["homologous_families"]           = config["homologous_families"]
            self.io_dict["species_tree"]                  = config["species_tree"]
            self.io_dict["output_directory"]              = config["output_directory"]
            self.io_dict["c1p_parameters"]                = config["c1p_parameters"]
            self.io_dict["wacs_file"]                     = config["wacs_directory"]

            self.markers_param_dict["markers_doubled"]    = config["markers_doubled"]
            self.markers_param_dict["markers_unique"]     = config["markers_unique"]
            self.markers_param_dict["markers_universal"]  = config["markers_universal"]
            self.markers_param_dict["markers_overlap"]    = config["markers_overlap"]
            self.markers_param_dict["filter_copy_number"] = config["filter_copy_number"]
            self.markers_param_dict["filter_by_id"]       = config["filter_by_id"]

            self.run_param_dict["all_match"]            = config["all_match"]

            self.debug = config["debug"]
            
        except IOError:
            print("{}  ERROR (master.py -> process.py) - could not open configuration file: {}\n"
                    .format(strtime(),config_file))
            sys.exit()
        config.clear()
    #enddef

    def parse_markersPhase(self, config_file_directory, len_input_arguments):
        """
        Deal with everything realated to input (configuration, markers and species pairs) in order to take information from these files

        """
        # Parse arguments: 
        #sys.argv[1] = ../data/configuration_file 
        self.setConfigParams(config_file_directory, len_input_arguments)
        self.setOutputStreams() #set Log and Debug
        # MasterMarkers class: methods used in order to deal with input files and take information from them (populate the species_pairs list and hom_fam_list)
        markers_phase_obj = MasterMarkers() 
        markers_phase_obj.setInputStreams(self.io_dict["species_tree"],self.io_dict["homologous_families"]) #set pairs file stream and hom fams file stream
        # Parse the hom fams file.
        self.hom_fam_list = markers_phase_obj.parseHomFamilies(self.io_dict["homologous_families"], self.log)
        self.getSpeciesList()
        # Parse the species pair file, put result in list.
        self.species_pairs = markers_phase_obj.parseSpeciesPairs(self.species_set, self.io_dict["species_tree"], self.log)
        # hom_fam_list and species_pairs list are now populated

        # Filter by ID
        if self.markers_param_dict["filter_by_id"]:
            self.hom_fam_list = markers_phase_obj.filterByID(self.markers_param_dict["filter_by_id"], self.hom_fam_list, self.log)

        # Filter by Copy Number
        if self.markers_param_dict["filter_copy_number"] != 0:
            self.hom_fam_list = markers_phase_obj.filterByCopyNumber(self.markers_param_dict["filter_copy_number"], self.hom_fam_list, self.log)

        #Get all overlapped pairs
        if self.markers_param_dict["markers_overlap"] == 1:
            self.overlapping_pairs_list = markers_phase_obj.getOverlappingPairs(self.hom_fam_list, self.log)
        
        # Since the markers are all oriented, double them.
        self.hom_fam_list = markers_phase_obj.doubleMarkers(self.hom_fam_list)


        del markers_phase_obj
    #enddef

    def adjacenciesPhase(self):
        # Genome construction
        self.genome_construction_obj.constructGenomes(self.species_pairs, self.hom_fam_list, self.log)
        self.genome_construction_obj.dealWithAdjPhase(self.species_pairs, self.hom_fam_list, self.io_dict["output_directory"], self.log, self.debug, self.run_param_dict["all_match"])
    #enddef










































    def c1pPhase(self):
        parameters_file = self.io_dict["c1p_parameters"]

        params = parameters.Parameters()
        params.from_file(parameters_file)
        
        code_dir = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + "/code/c1p_files" # directory where the code is stored
        
        event = threading.Event()
        
        working_dir = os.path.dirname(parameters_file)

        # set the working directory to the directory of the parameters file
        if working_dir and working_dir != '':
            os.chdir(working_dir)
        #endif
        #endif

        ## CREATING DIRECTORIES ------------------------------------------------------------------------------
        output_prefix = params.output_ancestor + "_"         # Prefix of all filenames
        output_dir = params.output_dir
        output_tmp = params.output_dir # directory for intermediate files

        # create output directory
        try:
            os.mkdir(output_dir)
        except:
            pass
        #endtry

        ## CREATING LOG FILE
        output_log = output_tmp + "/" + output_prefix + "LOG"    # captures code output

        log_file = open(output_log, 'w')
            
        print2("----> Started")

        # TRACKING
        quiet = False        # True if in quiet mode

        ## CHECKING INCOMPATIBILITIES ------------------------------------------------------------------------

        ## General errors

        # Missing input
        if not params.markers_provided:
            print2("ERROR: no markers file")
            sys.exit(-1)
        #endif
        if not params.species_tree_provided:
            print2("ERROR: no species tree file")
            sys.exit(-1)
        #endif

        # Conflicting models
        if params.c1p_circular and params.c1p_linear:
            print2("ERROR: chose either circular or linear C1P")
            sys.exit(-1)
        #endif

        if (not params.markers_unique) and (params.acs_ra):
            print2("ERROR: reliable adjacencies are not defined for repeated markers")
            sys.exit(-1)
        #endif

        if params.acs_correction >= 2 and (not params.c1p_spectral):
            print2("ERROR: corrected ACS with Xs require using the spectral seriation algorithm")
            sys.exit(-1)
        #endif    
        #if params.acs_correction < 2 and params.c1p_sandwich:
        #    print2("ERROR: the sandwich c1p requires using corrected ACS with X's ")
        #    sys.exit(-1)
        #endif    


        ## Current limitations

        if params.acs_ra and (not params.acs_sa):
            print2("ERROR: computing reliable adjacencies require computing supported adjacencies")
            # could be addressed by forcing params.acs_sa to be true if params.acs_ra is true
            sys.exit(-1)
        #endif

        if params.acs_weight!=1 and not params.acs_file_provided:
            print2("ERROR: weighting ACS by linear interpolation is mandatory currently ")
            sys.exit(-1)
        #endif    

        #if not params.markers_unique:
        #    print2("ERROR: data with non unique markers are currently not handled")
        #    sys.exit(-1)
        #endif

        #if params.c1p_sandwich:
        #    print2("ERROR: C1P sandwich not handled currently")
        #    sys.exit(-1)
        #endif
        #if params.c1p_spectral:
        #    print2("ERROR: spectral seriation not handled currently")
        #    sys.exit(-1)
        #endif

        if params.c1p_circular:
            if params.acs_aci or params.acs_sci or params.acs_mci:
                print2("ERROR: circular chromosomes can only be computed from adjacencies")
                sys.exit(-1)
                #endif
        #endif

        if params.c1p_spectral:
            if params.markers_doubled:
                print2("ERROR: spectral seriation can not be used with doubled markers")
                sys.exit(-1)
                #endif
            if params.c1p_circular:
                print2("ERROR: spectral seriation can not be used with circular chromosomes")
                sys.exit(-1)
                #endif
            if params.c1p_telomeres:
                print2("ERROR: spectral seriation can not be used with telomeric ACS")
                sys.exit(-1)
            #endif
        #endif


        ## TRACKING OPTIONS -----------------------------------------------------------------------

        #output = sys.stdout

        #if quiet:
        #    output = open(os.devnull, 'w')
        ##endif

        ## COPYING INPUT FILES -------------------------------------------------------------------------
        input_dir = output_tmp + "/INPUT"

        # create input directory
        try:
            os.mkdir(input_dir)
        except:
            pass
        #endtry

        if not quiet:
            print2("----> Copying input")
        #endif

        copied_params = copy.copy(params)

        copied_parameters = input_dir + "/" + output_prefix + "PARAMETERS"
        copied_markers = input_dir + "/" + output_prefix + "MARKERS"
        copied_tree = input_dir + "/" + output_prefix + "SPECIES_TREE"
        copied_acs = input_dir + "/" + output_prefix + "ACS"
        copied_pairs = input_dir + "/" + output_prefix + "SPECIES_PAIRS"

        copied_params.markers_input = output_prefix + "MARKERS"
        copied_params.species_tree = output_prefix + "SPECIES_TREE"
        copied_params.output_dir = ".."

        if os.path.abspath(params.markers_input) != os.path.abspath(copied_markers):
            shutil.copy(params.markers_input, copied_markers)
        #endif

        if os.path.abspath(params.species_tree) != os.path.abspath(copied_tree):
            shutil.copy(params.species_tree, copied_tree)
        #endif

        if params.acs_pairs_provided:
            copied_params.acs_pairs = output_prefix + "SPECIES_PAIRS"
            if os.path.abspath(params.acs_pairs) != os.path.abspath(copied_pairs):
                shutil.copy(params.acs_pairs, copied_pairs)
            #endif
        #endif
        if params.acs_file_provided:
            copied_params.acs_file = output_prefix + "ACS"
            if os.path.abspath(params.acs_file) != os.path.abspath(copied_acs):
                shutil.copy(params.acs_file, copied_acs)
            #endif
        #endif

        copied_params.to_file(copied_parameters)

        ## READING MARKERS ------------------------------------------------------------------------------

        # markers_dir = output_tmp + "/MARKERS"
        # #create markers directory
        # try:
        #     os.mkdir(markers_dir)
        # except:
        #     pass
        # #endtry

        # markers_file = markers_dir + "/" + output_prefix + "MARKERS"
        # if not quiet:
        #     print2("----> Filtering markers")
        # #endif
        # callprocess(["python", code_dir + "/MARKERS/MARKERS_filter_species.py", params.markers_input, params.species_tree, str(params.markers_unique), str(params.markers_universal), markers_file])
        # if params.markers_doubled:
        #     markers_file_doubled = markers_file + "_DOUBLED"
        #      if not quiet:
        #          print2("----> Doubling markers")
        #      #endif
        #      callprocess(["python", code_dir + "/MARKERS/MARKERS_double.py", markers_file, markers_file_doubled])
        #     markers_file=markers_file_doubled
        # #endif
        # if not quiet:
        #     print2("----> Markers: " + markers_file)
        # #endif

        markers_dir = output_tmp + "/MARKERS"
        markers_file = markers_dir + "/" + output_prefix + "MARKERS"


        ## COMPUTING ACS ------------------------------------------------------------------------------
        acs_dir = output_tmp + "/ACS"        # intermediate files from ACS code
        mci = acs_dir + "/" + output_prefix + "MCI"        # maximal common intervals file
        sci = acs_dir + "/" + output_prefix + "SCI"        # strong common intervals file
        aci = acs_dir + "/" + output_prefix + "ACI"        # strong common intervals file
        sa = acs_dir + "/" + output_prefix + "SA"        # supported adjacencies file
        ra = acs_dir + "/" + output_prefix + "RA"        # reliable adjacencies file
        acs = acs_dir + "/" + output_prefix + "ACS"        # ancestral contiguous sets file

        # # make ACS directory
        # try:
        #     os.mkdir(acs_dir)
        # except:
        #     pass
        # #endtry

        # # Clear ACS file
        # if not quiet:
        #     print2("----> Cleaning ACS file")
        # #endif
        # acsf = file(acs, 'w')
        # acsf.close()

        # # Computing species
        # if not quiet:
        #     print2("----> Computing species")
        # #endif

        # species=acs_dir + "/" + output_prefix + "SPECIES"
        # if not quiet:
        #     print2("------> All species: " + species)
        # #endif
        # callprocess(["python", code_dir + "/TREES/TREES_list_species.py", params.species_tree, 'ALL', species])

        # # collect ingroup and outgroup species
        # in_species = []
        # out_species = []
        # ingroup = False

        # for line in open(species):
        #     if line.rstrip() == '#ingroup':
        #         ingroup = True
                
        #         continue
        #     elif line.rstrip() == '#outgroup':
        #         ingroup = False
                
        #         continue
        #     #endif
            
        #     if ingroup:
        #         in_species.append(line.rstrip())
        #     else:
        #         out_species.append(line.rstrip())
        #     #endif
        # #endfor

        # # Computing species pairs
        # if not quiet:
        #     print2("----> Computing species pairs to compare")
        # #endif
        # if params.acs_pairs_provided:
        #     if not quiet:
        #         print2("------> File provided: " + params.acs_pairs)
        #         species_pairs=params.acs_pairs
        #     #endif
        # else:
        #     species_pairs=acs_dir + "/" + output_prefix + "PAIRS"
        #     if not quiet:
        #         print2("------> All informative pairs: " + species_pairs)
        #     #endif
        #     callprocess(["python", code_dir + "/TREES/TREES_list_species_pairs.py", params.species_tree, species_pairs])
        # #endif


        # ingroup = False

        # # Computing ACS
        # if not quiet:
        #     print2("----> Computing ancestral contiguous sets: " + acs)
        # #endif
        # for line in open(species_pairs):
        #     if line == None or line == '' or line == '\n':        continue
        #     #endif    
        #     sp = line.rstrip().split(' ')    
        #     if line.rstrip() == '#ingroup':
        #         ingroup = True
                
        #         continue
        #     elif line.rstrip() == '#outgroup':
        #         ingroup = False
                
        #         continue
        #     #endif
        #     if not quiet:
        #         print2("------> Species: " + sp[0] + " " + sp[1])
        #     #endif        
            
        #     acs_sp = acs + "_" + sp[0] + "_" + sp[1]
                
        #     file_list = []        # list of file to join to make the matrix
        #     if params.acs_mci:
        #         if not quiet:
        #             print2("--------> Computing maximum common intervals: " + mci + "_" + sp[0] + "_" + sp[1])
        #         #endif
        #         if params.markers_unique == 2:
        #             callprocess(["python", code_dir + "/ACS/ACS_compute_common_intervals_unique.py", markers_file, sp[0], sp[1], mci + "_" + sp[0] + "_" + sp[1], "MAX", str(int(params.c1p_circular)) ])
        #         else:
        #             callprocess(["python", code_dir + "/ACS/ACS_compute_common_intervals_nonunique.py", markers_file, sp[0], sp[1], mci + "_" + sp[0] + "_" + sp[1], "MAX", str(int(params.c1p_circular)) ])
        #         file_list.append(mci + "_" + sp[0] + "_" + sp[1]);
        #         #endif
        #     #endif
        #     if params.acs_sci:
        #         if not quiet:
        #             print2("--------> Computing strong common intervals: " + sci + "_" + sp[0] + "_" + sp[1])
        #         #endif    
        #         if params.markers_unique == 2:
        #             callprocess(["python", code_dir + "/ACS/ACS_compute_common_intervals_unique.py", markers_file, sp[0], sp[1], sci + "_" + sp[0] + "_" + sp[1], "STRONG", str(int(params.c1p_circular)) ])
        #         else:
        #             callprocess(["python", code_dir + "/ACS/ACS_compute_common_intervals_nonunique.py", markers_file, sp[0], sp[1], sci + "_" + sp[0] + "_" + sp[1], "STRONG", str(int(params.c1p_circular)) ])
        #         #endif
        #         file_list.append(sci + "_" + sp[0] + "_" + sp[1]);
        #     #endif
        #     if params.acs_aci:
        #         if not quiet:
        #             print2("--------> Computing all common intervals: " + aci + "_" + sp[0] + "_" + sp[1])
        #         #endif    
        #         if params.markers_unique == 2:
        #             callprocess(["python", code_dir + "/ACS/ACS_compute_common_intervals_unique.py", markers_file, sp[0], sp[1], aci + "_" + sp[0] + "_" + sp[1], "ALL", str(int(params.c1p_circular)) ])
        #         else:
        #             callprocess(["python", code_dir + "/ACS/ACS_compute_common_intervals_nonunique.py", markers_file, sp[0], sp[1], aci + "_" + sp[0] + "_" + sp[1], "ALL", str(int(params.c1p_circular)) ])
        #         #endif
        #          file_list.append(aci + "_" + sp[0] + "_" + sp[1]);
        #      #endif    
        #      if params.acs_sa:
        #          if not quiet:
        #              print2("--------> Computing supported adjacencies: " + sa + "_" + sp[0] + "_" + sp[1])
        #          #endif
        #         if params.markers_unique == 2:
        #             callprocess(["python", code_dir + "/ACS/ACS_compute_supported_adjacencies_unique.py", markers_file, sp[0], sp[1], sa + "_" + sp[0] + "_" + sp[1], str(int(params.c1p_circular))])        
        #         else:
        #             callprocess(["python", code_dir + "/ACS/ACS_compute_supported_adjacencies_nonunique.py", markers_file, sp[0], sp[1], sa + "_" + sp[0] + "_" + sp[1], str(int(params.c1p_circular))])        
        #         #endif
        #         file_list.append(sa + "_" + sp[0] + "_" + sp[1])
        #      #endif    
        #      if params.acs_ra:
        #          if ingroup:
        #              for sp3 in out_species:
        #                  if sp3 == sp[0] or sp3 == sp[1]:
        #                      continue
        #                  #endif
                     
        #                  if not quiet:
        #                      print2("--------> Computing reliable adjacencies: " + ra + "_" + sp[0] + "_" + sp[1] + "_" + sp3)
        #                  #endif
                     
        #                  callprocess(["python", code_dir + "/ACS/ACS_compute_reliable_adjacencies_unique.py", markers_file, sa + "_" + sp[0] + "_" + sp[1], sp[0], sp[1], sp3, ra + "_" + sp[0] + "_" + sp[1] + "_" + sp3, str(int(params.c1p_circular))])
        #                 file_list.append(ra + "_" + sp[0] + "_" + sp[1] + "_" + sp3)
                        
        #                 if not quiet:
        #                      print2("--------> Computing reliable adjacencies: " + ra + "_" + sp[1] + "_" + sp[0] + "_" + sp3)
        #                  #endif
                        
        #                 callprocess(["python", code_dir + "/ACS/ACS_compute_reliable_adjacencies_unique.py", markers_file, sa + "_" + sp[0] + "_" + sp[1], sp[1], sp[0], sp3, ra + "_" + sp[1] + "_" + sp[0] + "_" + sp3, str(int(params.c1p_circular))])
        #                 file_list.append(ra + "_" + sp[1] + "_" + sp[0] + "_" + sp3)
        #             #endfor
        #          else:
        #              if sp[0] in out_species:
        #                  for sp3 in in_species:
        #                      if sp3 == sp[0] or sp3 == sp[1]:
        #                          continue
        #                      #endif
                         
        #                      if not quiet:
        #                          print2("--------> Computing reliable adjacencies: " + ra + "_" + sp[0] + "_" + sp[1] + "_" + sp3)
        #                      #endif
                             
        #                      callprocess(["python", code_dir + "/ACS/ACS_compute_reliable_adjacencies_unique.py", markers_file, sa + "_" + sp[0] + "_" + sp[1], sp[0], sp[1], sp3, ra + "_" + sp[0] + "_" + sp[1] + "_" + sp3, str(int(params.c1p_circular))])
        #                     file_list.append(ra + "_" + sp[0] + "_" + sp[1] + "_" + sp3)
        #                 #endfor
        #             else:
        #                  for sp3 in in_species:
        #                      if sp3 == sp[0] or sp3 == sp[1]:
        #                          continue
        #                      #endif
                         
        #                      if not quiet:
        #                          print2("--------> Computing reliable adjacencies: " + ra + "_" + sp[1] + "_" + sp[0] + "_" + sp3)
        #                      #endif
                         
        #                     callprocess(["python", code_dir + "/ACS/ACS_compute_reliable_adjacencies_unique.py", markers_file, sa + "_" + sp[0] + "_" + sp[1], sp[1], sp[0], sp3, ra + "_" + sp[1] + "_" + sp[0] + "_" + sp3, str(int(params.c1p_circular))])
        #                     file_list.append(ra + "_" + sp[1] + "_" + sp[0] + "_" + sp3)
        #                 #endfor
        #             #endif
        #          #endif
        #      #endif
             
        #      if len(file_list) == 0:
        #          if not quiet:
        #              print2("------> Skipping ACS computation")
        #          #endif
                 
        #          break
        #      #endif
                
        #     prog = ["python", code_dir + "/ACS/ACS_join_files.py", "0"]        # join files program commandline
        #     prog.extend(file_list)
        #     prog.append(acs_sp)
        #     if not quiet:
        #         print2("------> Joining ACS files: " + acs_sp)
        #     #endif
             
        #      callprocess(prog)
                  
        #      # adding telomeres
        #     if params.c1p_telomeres > 0:
        #         acs_telomeres = acs_sp + "_TEL"
        #         if not quiet:
        #             print2("------> Adding telomeres")
        #         #endif
            
        #         callprocess(["python", code_dir + "/ACS/ACS_add_telomeres.py", markers_file, sp[0], sp[1], acs_sp, acs_telomeres])
                
        #         acs_sp = acs_telomeres
        #     #endif
            
        #     # Correcting computed ACS
        #     if (not ingroup and params.markers_universal == 1) or (params.markers_universal == 0) and params.acs_correction > 0:
        #         if params.acs_correction == 1:
        #             acs_corrected = acs_sp + "_MM1"
        #         elif params.acs_correction == 2:
        #             acs_corrected = acs_sp + "_MMX"
        #         elif params.acs_correction == 3:
        #             acs_corrected = acs_sp + "_MMX1"
        #         #endif
                
        #         if not quiet:
        #             print2("------> Correcting computed ancestral contiguous sets: " + acs_sp)
        #            #endif        
                    
        #         if params.acs_correction > 0:
        #             callprocess(["python", code_dir + "/ACS/ACS_add_missing_markers.py", markers_file, sp[0], sp[1], acs_sp, str(params.acs_correction), acs_corrected])
                
        #             acs_sp = acs_corrected
        #         #endif
        #     #endif
            
        #      if not quiet:
        #          print2("------> Joining ACS files: " + acs_sp + ", " + acs)
        #      #endif    
            
        #     callprocess(["python", code_dir + "/ACS/ACS_join_files.py", acs, acs_sp, acs])
        # #endfor

        # # Reading provided ACS (if any)
        # if params.acs_file_provided:
        #      if not quiet:
        #          print2("----> Adding provided ACS files: " + params.acs_file)
        #      #endif

        #     callprocess(["python", code_dir + "/ACS/ACS_join_files.py", acs, params.acs_file, acs])
        # #endif

        # ## WEIGHTING ACS ------------------------------------------------------------------------------
        # weight_dir = output_tmp + "/WEIGHT"
        # weight_input=acs
        # wacs = weight_dir + "/" + output_prefix + "WACS"    # weighted ancestral contiguous sets
        # # make WEIGHT directory
        # try:
        #     os.mkdir(weight_dir)
        # except:
        #     pass
        # #endtry
        # if params.acs_weight==1:
        #     if not quiet:
        #         print2("----> Weighting ACS: " + wacs)
        #         #endif
        #     if params.markers_doubled:
        #         callprocess(["python", code_dir +"/WEIGHT/WEIGHT_weight_linear_interpolation.py", weight_input, params.species_tree, wacs, "d"])
        #     else:
        #         callprocess(["python", code_dir +"/WEIGHT/WEIGHT_weight_linear_interpolation.py", weight_input, params.species_tree, wacs])
        #     #endif
        # elif params.acs_weight==0:        # weight already given
        #     wacs = acs
        # #endif

        weight_dir = output_tmp + "/WEIGHT"
        weight_input=acs
        # wacs = weight_dir + "/" + output_prefix + "WACS"    # weighted ancestral contiguous sets
        wacs = self.io_dict["wacs_file"]

        ## COMPUTING A C1P MATRIX ------------------------------------------------------------------------------
        c1p_dir = output_tmp + "/C1P"        # intermediate files from C1P code
        cars_dir = output_tmp + "/CARS"
        pqr_tree = cars_dir + "/" + output_prefix + "PQRTREE"        # PQR-tree file
        if params.markers_doubled:
            pqr_tree_doubled = cars_dir + "/" + output_prefix + "PQRTREE_DOUBLED"
        #endif
        # make C1P directory
        try:
             os.mkdir(c1p_dir)
        except:
             pass
        #endtry
        # make CARS directory
        try:
             os.mkdir(cars_dir)
        except:
             pass
        #endtry
        if params.c1p_linear or params.c1p_telomeres > 0:
            if not quiet:
                print2("----> Creating PQR-tree: " + pqr_tree)
            #endif

            if params.markers_doubled:
                callprocess(["python", code_dir +"/C1P/C1P_compute_PQRtree.py", wacs, pqr_tree_doubled, params.output_ancestor])
                if not quiet:
                    print2("----> Halving PQR-tree columns") 
                        #endif
                callprocess(["python", code_dir +"/C1P/C1P_halve_PQRtree.py", pqr_tree_doubled, pqr_tree])
            else:
                callprocess(["python", code_dir +"/C1P/C1P_compute_PQRtree.py", wacs, pqr_tree, params.output_ancestor])
            #endif
        elif params.c1p_circular:
            if not quiet:
                print2("----> Creating PQR-tree: " + pqr_tree)
            #endif

            if params.markers_doubled:
                callprocess(["python", code_dir +"/C1P/C1P_compute_PQCRtree.py", wacs, pqr_tree_doubled, params.output_ancestor])
                if not quiet:
                    print2("----> Halving PQR-tree columns") 
                        #endif
                callprocess(["python", code_dir +"/C1P/C1P_halve_PQRtree.py", pqr_tree_doubled, pqr_tree])
            else:
                callprocess(["python", code_dir +"/C1P/C1P_compute_PQCRtree.py", wacs, pqr_tree, params.output_ancestor])
            #endif
        #endif

        # heuristic
        if params.c1p_heuristic and params.c1p_telomeres == 0:
            if params.c1p_linear:
                do_c1p(code_dir, params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "HEUR", "heuristic", "/C1P/C1P_make_C1P_heuristic.py", "/C1P/C1P_compute_PQRtree.py")
            elif params.c1p_circular:
                do_c1p(code_dir, params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "HEUR", "heuristic", "/C1P/C1P_make_circC1P_heuristic.py", "/C1P/C1P_compute_PQCRtree.py", True)
            #endif
        #endif

        # branch and bound
        if params.c1p_bab and params.c1p_telomeres == 0:
            if params.c1p_linear:
                do_c1p(code_dir, params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "BAB", "branch and bound", "/C1P/C1P_make_C1P_branch_and_bound.py", "/C1P/C1P_compute_PQRtree.py")
            elif params.c1p_circular:
                do_c1p(code_dir, params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "BAB", "branch and bound", "/C1P/C1P_make_circC1P_branch_and_bound.py", "/C1P/C1P_compute_PQCRtree.py", True)
            #endif
        #endif

        # seriation
        if params.c1p_spectral:
            pq_tree = cars_dir + "/" + output_prefix + "PQTREE_SERIATION"    

            if params.acs_correction >= 2:
                if not quiet:
                    print2("----> Computing PQ-tree from correlation matrix of a ternary matrix using parameter alpha") 
                        #endif
                callprocess(["python", code_dir +"/SERIATION/SERIATION_compute_PQtree_dotproduct_correlation_matrix.py", wacs, str(params.c1p_spectral_alpha), pq_tree, params.output_ancestor])
            else:
                if not quiet:
                    print2("----> Computing PQ-tree using spectral seriation on correlation matrix") 
                        #endif
                callprocess(["python", code_dir +"/SERIATION/SERIATION_compute_PQtree_correlation_matrix.py", wacs, pq_tree, params.output_ancestor])

                #endif
        #endif

        # telomere heuristic
        if params.c1p_telomeres == 1:
            do_c1p(code_dir, params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "TEL_HEUR", "heuristic", "/C1P/C1P_make_mC1P_heuristic.py", "/C1P/C1P_compute_PQRtree.py")
        #endif

        # telomere branch and bound, add after
        if params.c1p_telomeres == 2:
            do_c1p(code_dir, params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "TEL_BAB1", "branch and bound, added after", "/C1P/C1P_make_mC1P_branch_and_bound_both.py", "/C1P/C1P_compute_PQRtree.py")
        #endif

        #telomere branch and bound, add during
        if params.c1p_telomeres == 3:
            do_c1p(code_dir, params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "TEL_BAB2", "branch and bound, added during", "/C1P/C1P_make_mC1P_branch_and_bound.py", "/C1P/C1P_compute_PQRtree.py")
        #endif

        if not quiet:
             print2("----> Done")
        #endif

        log_file.close()




























    def genomeConstructionPhase(self):
        self.genome_construction_obj.dealWithConstructionPhase(self.hom_fam_list, self.io_dict["output_directory"], self.log)
    #enddef

    # finds a HomFam in a list given a marker id
    # hom_fam_list - list of HomFam: the list to look in
    # marker_id - str: the marker_id to look for
    # Return - HomFam: the hom. family with marker id, marker_id
    def findHomFam(marker_id):
        # Filter the list
        matches = filter( lambda hom_fam: hom_fam.id == marker_id, self.hom_fam_list )
        # Return first element of the list
        if matches:
            return matches[0]
        else:
            return None
    #enddef

    def closeLogFile(self):
        self.log.close()
    #enddef

    def closeDebugFile(self):
        self.log.close()
    #enddef

    def closeAllFiles(self):
        self.closeLogFile()
        self.closeDebugFile()
    #enddef

    def getSpeciesList(self):
        for marker_family in self.hom_fam_list:
            for locus_index,locus in enumerate(marker_family.loci):
                self.species_set.add(marker_family.loci[locus_index].species)

    def getIODictionary(self):
        return self.io_dict, self.markers_param_dict
    #enddef

    def getFileStreams(self):
        return self.log, self.debug
    #enddef

    def getSpeciesPairs(self):
        return self.species_pairs
    #enddef

    def getHomFamList(self):
        return self.hom_fam_list
    #enddef
#endclass
