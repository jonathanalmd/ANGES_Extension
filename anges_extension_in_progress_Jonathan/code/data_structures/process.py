import sys
import time

from collections import defaultdict

from data_structures import markers
from data_structures import intervals
from data_structures import genomes
from data_structures import comparisons

import optimization
import assembly
    
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

    def parseSpeciesPairs(self, log):
        """
        Populating the species_pairs list
        """
        species_pairs = []
        for pair in self.pairs_file_stream:
        # Assume format of "<species1> <species2>", respect comments
            pair = pair.strip()
            pair = pair.split()
            if pair[0][0] != "#" and len( pair ) == 2:
                species_pairs.append( pair )
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

    def solveAdjacencies(self, species_pairs, gens, output_directory, log):
        for pair in species_pairs:
            new_adjacencies = comparisons.find_adjacencies( gens[ pair[0] ],
                                                            gens[ pair[1] ] )
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

    def solveRSIs(self, species_pairs, gens, output_directory, log):
        for pair in species_pairs:
            new_RSIs = comparisons.find_RSIs(gens[ pair[0] ], gens[ pair[1] ] )
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

    def dealWithAdjPhase(self, species_pairs, hom_fam_list, output_directory, log, debug):
        # For each pair of species, compare the species to find adjacencies.
        self.adj.solveAdjacencies(species_pairs, self.gens, output_directory, log)

        # Do the same for repeat spanning intervals
        self.RSI.solveRSIs(species_pairs, self.gens, output_directory, log)

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

        self.log = None
        self.debug = None
    
        self.species_pairs = []
        self.hom_fam_list = []
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

            self.markers_param_dict["markers_doubled"]    = config["markers_doubled"]
            self.markers_param_dict["markers_unique"]     = config["markers_unique"]
            self.markers_param_dict["markers_universal"]  = config["markers_universal"]
            self.markers_param_dict["markers_overlap"]    = config["markers_overlap"]
            self.markers_param_dict["filter_copy_number"] = config["filter_copy_number"]
            self.markers_param_dict["filter_by_id"]       = config["filter_by_id"]

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
        # Parse the species pair file, put result in list.
        self.species_pairs = markers_phase_obj.parseSpeciesPairs(self.log)
        # Parse the hom fams file.
        self.hom_fam_list = markers_phase_obj.parseHomFamilies(self.io_dict["homologous_families"], self.log)
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
        self.genome_construction_obj.dealWithAdjPhase(self.species_pairs, self.hom_fam_list, self.io_dict["output_directory"], self.log, self.debug)
    #enddef

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
