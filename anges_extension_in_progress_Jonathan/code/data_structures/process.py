
import sys
import time

from collections import defaultdict

from data_structures import markers
from data_structures import intervals
from data_structures import genomes
from data_structures import comparisons

import optimization
import assembly

    
def getOverlappingPairs(hom_fams):
    """
    getOverlappingPairs: receives a list of hom_fams and returns a list of overlapping pairs (each pair is a tuple containing two Locus)
    hom_fams - HomFam: list of objects from HomFam class
    """
        #import pdb; pdb.set_trace()
    loci_dict = defaultdict(list)
    overlapping_pairs_list = []
    for hom_fam_index, marker_family in enumerate(hom_fams):
        for locus_index, locus in enumerate(marker_family.loci):
           loci_dict[locus.species].append((hom_fam_index, locus_index))
        #endfor
    #endfor
    for species, species_indexes in loci_dict.items():
        # species name and a list of tuples with (hom_fam_index, locus_index)
        # hom_fam_index => access to hom_fam ID and loci list
        #print (species, species_indexes)
        i = 1
        locus1 = hom_fams[species_indexes[0][0]].loci[species_indexes[0][1]]
        while i < len(species_indexes):
            locus2 = hom_fams[species_indexes[i][0]].loci[species_indexes[i][1]]
            #print (locus1, locus2)
            is_overlapping_pair = locus1.overlappingPairs(locus2)
            if is_overlapping_pair:
                overlapping_pairs_list.append(is_overlapping_pair)
            #endif
            i = i + 1
        #endwhile
    #endfor
    return overlapping_pairs_list
#enddef

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

def strtime():
    """
    Function to format time to string
    """
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())


def write_intervals( ints, f, log ):
    """
    Function to write intervals to file
    Arguments:
    ints: IntervalDict
    f: the file to write to
    log: log stream to write errors to
    """
    try:
        output = open( f, 'w' )
        for interval in ints.itervalues():
            output.write( str( interval ) + "\n" )
    except IOError:
        log.write( "%s  ERROR (master.py) - could not write intervals to"
                   "file: %s\n" %( strtime(), f ) )
        sys.exit()


class MasterScript:
    def __init__(self):
        self.io_dict = {}
        self.markers_param_dict = {}
        self.log = None
        self.debug = None
        self.hom_fams_file_stream = None
        self.pairs_file_stream = None

        self.species_pairs = []
        self.hom_fam_list = []
        self.gens = {}
        self.adjacencies = intervals.IntervalDict()
        self.RSIs = intervals.IntervalDict()

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
        self.printDictInfo()
    #enddef

    def getIODictionary(self):
        return self.io_dict, self.markers_param_dict

    def printDictInfo(self):
        print("Configuration File info:\n")
        for key, value in self.io_dict.items():
            print(key, value)
        for key, value in self.markers_param_dict.items():
            print(key, value)
        print("\n")

    def setFileStreams(self):
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
        try:
            self.pairs_file_stream = open(self.io_dict["species_tree"], 'r')
        except IOError:
            log.write( "{}  ERROR (master.py) - could not open pairs file: {}\n"
                       .format(strtime(), self.io_dict["species_tree"]))
            sys.exit()
        
        try:
            self.hom_fams_file_stream = open(self.io_dict["homologous_families"], 'r')
        except IOError:
            log.write( "{}  ERROR (master.py) - could not open homologous families file: {}\n"
                       .format(strtime(), self.io_dict["homologous_families"]))
            sys.exit()

    def closeLogFile(self):
        self.log.close()

    def closeDebugFile(self):
        self.log.close()

    def closePairsFile(self):
        self.pairs_file_stream.close()

    def closeHomFamsFile(self):
        self.hom_fams_file_stream.close()

    def closeAllFiles(self):
        self.closeLogFile()
        self.closeDebugFile()
        self.closePairsFile()
        self.closeHomFamsFile()

    def getFileStreams(self):
        return self.log, self.debug, self.hom_fams_file_stream, self.pairs_file_stream

    def getSpeciesPairs(self):
        return self.species_pairs

    def getHomFamList(self):
        return self.hom_fam_list
    
    def getGenomes(self):
        return self.gens

    def getAdjacencies(self):
        return self.adjacencies

    def getRSIs(self):
        return self.RSIs

    def parseSpeciesPairs(self):
        for pair in self.pairs_file_stream:
        # Assume format of "<species1> <species2>", respect comments
            pair = pair.strip()
            pair = pair.split()
            if pair[0][0] != "#" and len( pair ) == 2:
                self.species_pairs.append( pair )
        self.pairs_file_stream.close()
        self.log.write( "{}  Read {} species pairs.\n"
                   .format(strtime(), len(self.species_pairs)))
        self.log.flush()

    # reads hom. families from a file
    # file_name - str: the name of the file to read from
    # hom_fam_list - list of HomFam: the list to add to (Default = [])
    # Return - list of HomFam: the list of hom. familes read
    def read_hom_families_file(self):
        line = self.hom_fams_file_stream.readline()

        while len(line) > 0:
            trunc_line = line.strip()

            if len(trunc_line) > 0:
                if trunc_line[0] == '>':
                    # read first hom. family
                    hom_fam, line = markers.HomFam.from_file(self.hom_fams_file_stream, trunc_line)
                    if hom_fam != None:
                        self.hom_fam_list.append(hom_fam)
                    #endif

                    # read the rest of the hom. families
                    while len(line) > 0:
                        hom_fam, line = markers.HomFam.from_file(self.hom_fams_file_stream, line)
                        if hom_fam != None:
                            self.hom_fam_list.append(hom_fam)
                        #endif
                    #endwhile
                else:
                    line = self.hom_fams_file_stream.readline()
                #endif
            else:
                line = self.hom_fams_file_stream.readline()
            #endif
        #endif

        self.hom_fams_file_stream.close()

        self.log.write( "{}  Read homologous families from file.\n"
               .format( strtime() ) )
        self.log.flush()
    #enddef

    def doubleMarkers(self):
        self.hom_fam_list = genomes.double_oriented_markers(self.hom_fam_list)

    def constructGenomes(self):
        self.gens = genomes.get_genomes(
            self.hom_fam_list,
            list( set( species for pair in self.species_pairs for species in pair ) ),
            )
        self.log.write( "{}  Constructed genomes of {} species.\n" 
                        .format( strtime(), len( self.gens ) ) )
        self.log.flush()

    def solveAdjacencies(self):
        for pair in self.species_pairs:
            new_adjacencies = comparisons.find_adjacencies( self.gens[ pair[0] ],
                                                            self.gens[ pair[1] ] )
            comparisons.add_intervals( self.adjacencies, new_adjacencies )
        comparisons.set_interval_weights( self.adjacencies )
        self.log.write( "{}  Found {} adjacencies with total weight of {}.\n"
                   .format( strtime(),
                      len( self.adjacencies ),
                      self.adjacencies.total_weight ) )
        self.log.flush()
        write_intervals( self.adjacencies, self.io_dict["output_directory"] + "/adjacencies", self.log )

    def solveRSIs(self):
        for pair in self.species_pairs:
            new_RSIs = comparisons.find_RSIs( self.gens[ pair[0] ], self.gens[ pair[1] ] )
            comparisons.add_intervals( self.RSIs, new_RSIs )
        comparisons.set_interval_weights( self.RSIs )
        self.log.write( "{}  Found {} repeat spanning intervals with total weight of"
        " {}.\n" .format( strtime(), len( self.RSIs ), self.RSIs.total_weight ) )
        self.log.flush()
        write_intervals( self.RSIs, self.io_dict["output_directory"] + "/RSIs", self.log )
#TODO:
# - finish the master.py integration with process.py (much work)
