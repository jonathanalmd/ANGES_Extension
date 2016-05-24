
# locus
class Locus:
    # Locus constructor.
    # Parameters
    # species: str - species name
    # chromosome: str - chromosome name
    # start: int - locus start
    # end: int - locus end
    # orientation: int - locus orientation (positive for forwards, negative for reversed and 0 for unoriented)
    # comment: str - comment
    def __init__(self, species, chromosome, start, end, orientation, comment):
        self.species = species
        self.chromosome = chromosome
        self.start = start
        self.end = end
        if self.start > self.end: raise ValueError("The Locus ends before it begins.")
        self.orientation = orientation
        self.comment = comment

    # Comparison functions so that loci can be checked for equality and ordered.
    def __eq__( self, other ):
        return ( self.species, self.chromosome, self.start, self.end, self.orientation ) == \
               ( other.species, other.chromosome, other.start, other.end, other.orientation )

    def __lt__( self, other ):
        return ( self.species, self.chromosome, self.start, self.end, self.orientation ) < \
               ( other.species, other.chromosome, other.start, other.end, other.orientation )

    # Returns a string representing the locus with given the separators. Will not print the comment if it is empty.
    # <species name>sp1<chromosome name>sp2<start>sp3<end>sp4<orientation>sp5<comment>
    # Parameters
    # sp1: str
    # sp2: str
    # sp3: str
    # sp4: str
    # sp5: str
    # Return: str - string representation of the locus
    def str_general(self,sp1,sp2,sp3,sp4,sp5):
        if self.orientation > 0:
            orient_str = "+"
        elif self.orientation < 0:
            orient_str = "-"
        else: # == 0
            orient_str = "X"
        #endif
        if len(self.comment) > 0:
            comment_str = sp5 +     self.comment
        else:
            comment_str = ""
        #endif
        return self.species + sp1 + self.chromosome + sp2 + str(self.start) + sp3 + str(self.end) + sp4 +orient_str + comment_str
    #enddef

    # Returns a string representing the locus. Will not print the comment is it is empty.
    # <species name>.<chromosome name>:<start>-<end> <orientation> #<comment>
    # or (if the species name contains '.')
    # <species name> <chromosome name> <start> <end> <orientation> #<comment>
    # Return: str - string representation of the locus
    def __str__(self):
        if self.species.find(".") < 0:
            return self.str_general(".",":","-"," "," #")
        else:
            return self.str_general(" "," "," "," "," #")
        #endif
    #enddef

    
    def __repr__(self):
        #Return the direct representation of Locus class
        return "Locus {}.{}:{}-{} {}".format(self.species, self.chromosome, self.start, self.end, self.orientation)

    # Static method which parses a locus from a string.
    # <species name>.<chromosome name>:<start>-<end> <orientation> #<comment>
    # or
    # <species name> <chromosome name> <start> <end> <orientation> #<comment>
    # Parameters
    # string - str: the string to parse
    # Return Locus - the locus in the file (None on failure)
    @staticmethod
    def from_string(string):
        #assert type(string) is str

        trunc_str = string.strip()
        full_str = trunc_str

        if len(trunc_str) == 0 or trunc_str[0] == '>' or trunc_str[0] == '#':
            return None
        #endif

        # Extracting the comments field
        comment_split = trunc_str.split('#', 1)
        if len(comment_split) > 1:
            comment = comment_split[1]
        else:
            comment = ""
        #endif
        trunc_str = comment_split[0].strip()

        # using space separated format
        split_str =     trunc_str.split()
        if len(split_str) >= 4:
            species = split_str[0]
            chromosome = split_str[1]

            try:
                start = int(split_str[2])
            except:
                print("Warning: start coordinate not an integer for locus. String : \'" + full_str + "\'")
                return None
            #endtry

            try:
                end = int(split_str[2])
            except:
                print("Warning: end coordinate not an integer for locus. String : \'" + full_str + "\'")
                return None
            #endtry

            if len(split_str) == 4:
                orientation = 0
            else:
                trunc_str = split_str[4].strip()

                if trunc_str[0] == '+':
                    orientation = 1
                elif trunc_str[0] == '-':
                    orientation = -1
                elif trunc_str[0].lower == 'x':
                    orientation = 0
                else:
                    print("Warning: unknown orientation for locus. String : \'" + full_str + "\'")
                    return None
                #endif
            #endif

            return Locus(species, chromosome, start, end, orientation, comment)
        #endif

        # using old format
        species, trunc_str = Locus.obj_split(trunc_str, ".", 2, full_str, "species")
        if trunc_str == None:
            return None
        #endif

        chromosome, trunc_str = Locus.obj_split(trunc_str, ":", 2, full_str, "chromosome")
        if trunc_str == None:
            return None
        #endif

        start, trunc_str = Locus.obj_split(trunc_str, "-", 2, full_str, "start")
        if trunc_str == None:
            return None
        #endif
        try:
            start = int(start)
        except:
            print("Warning: start coordinate not an integer for locus. String : \'" + full_str + "\'")
            return None
        #endtry

        end, trunc_str = Locus.obj_split(trunc_str, " ", 1, full_str, "end")
        if trunc_str == None:
            return None
        #endif
        try:
            end = int(end)
        except:
            print("Warning: end coordinate not an integer for locus. String : \'" + full_str + "\'")
            return None
        #endtry

        if len(trunc_str) == 0:
            orientation = 0
        else:
            trunc_str = trunc_str.strip()

            if trunc_str[0] == '+':
                orientation = 1
            elif trunc_str[0] == '-':
                orientation = -1
            elif trunc_str[0].lower == 'x':
                orientation = 0
            else:
                print("Warning: unknown orientation for locus. String : \'" + full_str + "\'")
                return None
            #endif
        #endif

        return Locus(species, chromosome, start, end, orientation, comment)
    #enddef

    # Static auxilary method for from_string.
    # Parameters:
    # string - str
    # delimeter - str
    # size - int
    # full_string - str
    # warning name - str
    # Return - str, str
    @staticmethod
    def obj_split(string, delimiter, size, full_string, warning_name):
        #assert type(string) is str
        #assert type(delimiter) is str
        #assert type(size) is int
        #assert type(full_string) is str
        #assert type(warning_name) is str

        split = string.split(delimiter, 1)

        if len(split) < size:
            print("Warning: could not process " + warning_name + " for locus. String: \'" + full_string + "\'")

            return None, None
        #endif

        if len(split) < 2:
            split.append(None)
        #endif

        return split[0], split[1]
    #enddef

    def isOverlap(self, other):
        """
        isOverlap: compares two locus and returns True if the two are overlapped
        self - Locus: the first locus to be compared
        other - Locus: the second locus to be compared
        """
        if not isinstance(other, Locus):
            raise TypeError("Trying to compare two instances with different types")
        #endif
        if self.species == other.species and self.chromosome == other.chromosome:
            if self.start <= other.start:
                return self.end >= other.start 
            else: # self.start >= other.start
                return other.end >= self.start
            #endif
        else:
            return False
        #endif
    #enddef

    # isInclusion compares two locus and returns True if one of the two markers is included
    # self - Locus
    # other - another Locus
    def isInclusion(self, other):
        """
        isInclusion: compares two locus and returns True if one of the two markers is included 
        to the other one.
        self - Locus: the first locus to be compared
        other - Locus: the secont locus to be compared with the first
        """
        if not isinstance(other, Locus):
            raise TypeError("Trying to compare a Locus vs non-Locus")
        if self.species == other.species and self.chromosome == other.chromosome:
            if self.start >= other.start:
                return self.end <= other.end
            else: # other.start >= self.start
                return other.end <= self.end
            #endif
        else:
            return False
        #endif
    #enddef
    
    def overlappingPairs(self, other):
        """
        overlappingPairs: receives two Locus and checks if there is a overlap/inclusion
        Returns a tuple with the two overlapped loci (if isOverlap) or a empty tuple
        self - Locus: the first Locus object to be compared
        other - Locus: the second Locus object to be compared with the first one
        """
        if not isinstance(other, Locus):
            raise TypeError("Trying to compare two instances with different types")
        if self.isOverlap(other):
            return (self,other)
        else:
            return ()

    #enddef

#endclass


# homologous family
class HomFam:
    # HomFam constructor
    # ident - str: marker identity
    # loci - list of Locus: hom. family loci
    # copy_number - int: copy number
    # comment - str: comment
    def __init__(self, ident, loci, copy_number, comment):
        #assert type(ident) is str
        #assert type(loci) is list
        #assert type(copy_number) is int
        #assert type(comment) is str

        self.id = ident
        self.loci = loci
        self.copy_number = copy_number
        self.comment = comment
    #enddef

    # tests equality of HomFams
    # other - HomFam
    # Return - boolean: compares ids
    def __eq__(self, other):
        if type(other) is HomFam:
            return self.id == other.id
        else:
            return 0
        #endif
    #enddef

    # returns the opening string of the HomFam in the form
    # \><id> <copy number> #<comment>
    # Return - str
    def opening_string(self):
        if len(self.comment) > 0:
            comment_str = " #" + self.comment
        else:
            comment_str = ""
        #endif

        if self.copy_number != 1:
            copy_str = " " +  str(self.copy_number)
        else:
            copy_str = ""
        #endif

        return ">" + self.id + copy_str + comment_str
    #enddef

    # returns a string representing itself in the format
    # <id> <copy number> #<comment>
    # locus1
    # locus2
    # ...
    # locusn
    # Return - str
    def __str__(self):
        string = self.opening_string() + "\n"

        for l in self.loci:
            string = string + str(l) + "\n"
        #endfor

        return string
    #enddef
    def __repr__(self):
        """
        __repr__: direct representation of HomFam class 
        """
        return  "{}\nTotal loci: {}\n".format(self.opening_string(), len(self.loci))
            

    # writes the string representation to file
    # see __str__ for format
    # file_stream - file: file stream to write to
    def to_file(self, file_stream):
        #assert type(file_stream) is file
        file_stream.write(str(self))
    #enddef

    # reads the HomFam from a file
    # file_stream - the file_stream to read from
    # first_line - the first line of the file_stream if already read from otherwise None
    # Return - HomFam
    @staticmethod
    def from_file(file_stream, first_line):
        #assert type(file_stream) is file
        #assert type(first_line) is str

        if first_line == None:
            line = file_stream.readline()
        else:
            line = first_line
        #endif

        if len(line) <= 0:
            return None, line
        #endif

        hom_fam = HomFam.opening_from_string(line)

        line = file_stream.readline()

        while len(line) > 0:
            line = line.strip()

            if len(line) > 0:
                if line[0]== '>':
                    break
                #endif

                if hom_fam != None:
                    locus = Locus.from_string(line)

                    if locus != None:
                        hom_fam.loci.append(locus)
                    #endif
                #endif
            #endif

            line = file_stream.readline()
        #endwhile

        return hom_fam, line
    #enddef

    # parses a HomFam from a string
    # see opening_str for format
    # string - str: the string to parse
    @staticmethod
    def opening_from_string(string):
        trunc_line = string.strip()

        if len(trunc_line) < 0 or trunc_line[0] != '>':
            print("Warning: not a hom. family. String: \'" + string + "\'")

            return None
        #endfor

        comment_split = trunc_line.split('#', 1) # 1 = number of lines to be made
        if len(comment_split) >= 2:
            comment = comment_split[1]
        else:
            comment = ""
        #endif
        trunc_line = comment_split[0].strip()

        split = trunc_line.split(' ', 1) # 1 = number of lines to be made

        if len(split[0]) > 0:
            ident = split[0][1:] 
        else:
            print("Warning: hom. family missing identity: \'" + string + "\'")

            return None
        #endif

        copy_number = 1
        if len(split) > 1:
            if len(split[1]) > 0:
                copy_number = split[1].strip()

                try:
                    copy_number = int(copy_number)
                except:
                    print("Warning: copy number not an integer for homology family. String : \'" + string + "\'")
                    return None
                #endtry
            else:
                print("Warning: unknown information in hom. family. String: \'" + string + "\'")
            #endif
        #endif

        return HomFam(ident, [], copy_number, comment)
    #enddef



#endclass


# reads hom. families from a file
# file_name - str: the name of the file to read from
# hom_fam_list - list of HomFam: the list to add to (Default = [])
# Return - list of HomFam: the list of hom. familes read
def read_hom_families_file(file_name, hom_fam_list = []):
    file_stream = open(file_name)

    line = file_stream.readline()

    while len(line) > 0:
        trunc_line = line.strip()

        if len(trunc_line) > 0:
            if trunc_line[0] == '>':
                # read first hom. family
                hom_fam, line = HomFam.from_file(file_stream, trunc_line)
                if hom_fam != None:
                    hom_fam_list.append(hom_fam)
                #endif

                # read the rest of the hom. families
                while len(line) > 0:
                    hom_fam, line = HomFam.from_file(file_stream, line)
                    if hom_fam != None:
                        hom_fam_list.append(hom_fam)
                    #endif
                #endwhile
            else:
                line = file_stream.readline()
            #endif
        else:
            line = file_stream.readline()
        #endif
    #endif

    file_stream.close()

    return hom_fam_list
#enddef


# writes hom. families to a file
# file_name - str: the name of the file to write to
# hom_fam_list - list of HomFam: the list to write (Default = [])
def write_hom_families_file(file_name, hom_fam_list):
    file_stream = open(file_name, 'w')

    for hom_fam in hom_fam_list:
        hom_fam.to_file( file_stream )
        file_stream.write("\n")
        file_stream.flush()
    #endif

    file_stream.close()
#enddef


# finds a HomFam in a list given a marker id
# hom_fam_list - list of HomFam: the list to look in
# marker_id - str: the marker_id to look for
# Return - HomFam: the hom. family with marker id, marker_id
def find_marker(hom_fam_list, marker_id):
    # Filter the list
    matches = filter( lambda hom_fam: hom_fam.id == marker_id, hom_fam_list )
    # Return first element of the list
    if matches:
        return matches[0]
    else:
        return None
#enddef


# Function that returns the sibling of the head/tail of an oriented, doubled marker.
# Arguments:
#    marker: string
# Output:
#    string - the 'sibling' of marker
def sibling( marker ):
    if marker[ -2 : ] == "_h":
        return marker[ : -2 ] + "_t"
    elif marker[ -2 : ] == "_t":
        return marker[ : -2 ] + "_h"
    return ""

# Function that tells if two markers are siblings.
# Arguments:
#   marker1, marker2: string
# Output:
#   boolean: if marker1 and marker2 are siblings.
def are_siblings( marker1, marker2 ):
    return marker1 == sibling( marker2 )

# Functions that translate between doubled / undoubled markers.
# Arguments:
#   marker: string
# Output:
#   (to) (string, string)
#   (from) string
def to_doubled( marker ):
    return ( marker + "_h", marker + "_t" )
def from_doubled( marker ):
    return marker[:-2]
