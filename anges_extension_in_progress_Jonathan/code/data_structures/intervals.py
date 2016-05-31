from decimal import *

import markers


# intervals of markers and adjacencies, syntenies
class Interval:
    # Interval constructor
    # ident - str: the interval identity
    # marker_ids - list of strings: the list of marker IDs in the interval
    # loci - list of Locus: the list of loci the interval is at
    # order - Order: the ordering of the intrerval (linear, unordered)
    # weight - Decimal: the weight of the interval (~probability of appearance in the ancestor)
    # comment - str: comment
    def __init__(self, id, marker_ids, loci, order, weight, comment):
        self.id = id
        self.marker_ids = marker_ids
        self.loci = loci
        self.order = order
        self.weight = weight
        self.comment = comment
    #enddef

    # returns the header of the interval in the form
    # \><id> <weight> <order>  <marker1> <marker2> ... <markern> # comment
    # Return - str
    def opening_string(self):
        if len(self.comment) > 0:
            comment_str = " #" + self.comment
        else:
            comment_str = ""
        #endif

        string = ">" + self.id +  " " + str(self.weight) + " " + str(self.order) + " "

        for m in self.marker_ids:
            string = string + " " + m
        #endfor

        return string + comment_str
    #enddef

    # returns the string of the interval in the form
    # \><id> <weight> <order>  <marker1> <marker2> ... <markern> # comment
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

    # writes the interval to file
    # see __str__ for format
    # file_stream - file: the file_stream to write to
    def to_file(self, file_stream):
        file_stream.write(str(self))
    #enddef

    # reads an interval from a file
    # file_stream - the file_stream to read from
    # first_line - the first line of the file_stream if already read from otherwise None
    # Return - Interval
    @staticmethod
    def from_file(file_stream, first_line):
        if first_line == None:
            line = file_stream.readline()
        else:
            line = first_line
        #endif

        if len(line) <= 0:
            return None, line
        #endif

        interval = Interval.opening_from_string(line)

        line = file_stream.readline()

        while len(line) > 0:
            line = line.strip()

            if len(line) > 0:
                if line[0]== '>':
                    break
                #endif

                if interval != None:
                    locus = markers.Locus.from_string(line)

                    if locus != None:
                        interval.loci.append(locus)
                    #endif
                #endif
            #endif

            line = file_stream.readline()
        #endwhile

        return interval, line
    #enddef

    # parses an Interval from a string
    # see opening_string for format
    # string - str: the string to parse
    @staticmethod
    def opening_from_string(string):
        #assert type(string) is str

        trunc_line = string.strip()

        if len(trunc_line) <= 0 or trunc_line[0] != '>':
            print("Warning: not an interval. String: \'" + string + "\'")

            return None
        #endfor

        trunc_line = trunc_line[1:]

        comment_split = trunc_line.split('#', 1)
        if len(comment_split) >= 2:
            comment = comment_split[1]
        else:
            comment = ""
        #endif
        trunc_line = comment_split[0].strip()

        split = trunc_line.split()

        if len(split) < 4:
            print("Warning: not enough entries for interval. String: \'" + string + "\'")

            return None
        #endfor

        ident = split[0].strip()

        try:
            weight = Decimal(split[1].strip())
        except:
            print("Warning: weight not a number for interval. String: \'" + string + "\'")

            return None
        #endfor

        try:
            order_type = int(split[2].strip())
        except:
            print("Warning: order type unknown for interval. String: \'" + string + "\'")

            return None
        #endfor

        marker_ids = []

        for i in xrange(3, len(split)):
            marker_ids.append( split[i].strip() )
        #endfor

        order = Order(order_type)

        return Interval(ident, marker_ids, [], order, weight, comment)
    #enddef

    # Function to check for equality between intervals. Based on the
    # marker_ids, accepts reversed order of the IDs.
    def __eq__( self, other ):
        return ( self.marker_ids == other.marker_ids
                 or reversed( self.marker_ids ) == other.marker_ids )

#endclass

# marker ordering type in an interval
class Order:
    # Order constructor
    # order_type - int: type of order, 0 if unordered and 1 if ordered
    #UNUSED: markers_list - list of IDs of HomFam: interval's list of markers
    def __init__(self, order_type):
        #UNUSED: assert type(markers_list) is list
        #UNUSED: self.markers = markers_list

        #assert type(order_type) is int
        self.order_type = order_type

        if order_type == 0:
            self.is_ordered = False
        else:
            self.is_ordered = True
        #endif
        #UNUSED: self.__iter__ = markers_list.__iter__
    #enddef

    # returns a string representation of the Order for file saving
    # return - str: '0' if unordered and '1' if ordered
    def __str__(self):
        return str(self.order_type)
    #enddef
#enddef

# Class to store intervals. Offers average constant time lookup, insertion and
# deletion based on marker ID list, and average constant time query of intervals
# that contain a specified ID as endpoint.
class IntervalDict:
    def __init__( self ):
        # Main interval dictionary.
        self.ints = {}
        # To look up intervals by endpoint, keep an additional dictionary that
        # maps each marker ID to a list keys for 'ints' (list of marker IDs).
        self.endpoints = {}

    # Custom hash function for list of marker IDs, hashes ids and reversed(ids)
    # to the same value. The results will be the keys in the internal dict.
    def hash_ids( self, ids ):
        ids1 = tuple( ids )
        ids2 = tuple( reversed( ids ) )
        if ''.join(ids1) < ''.join(ids2):
            return hash(ids1)
        else:
            return hash(ids2)

    def __getitem__( self, ids ):
        return self.ints[self.hash_ids(ids)]

    def __setitem__( self, ids, interval ):
        self.ints[self.hash_ids(ids)] = interval
        # Manage endpoint database:
        id1, id2 = ids[0], ids[-1]
        if id1 not in self.endpoints:
            self.endpoints[ id1 ] = []
        if id2 not in self.endpoints:
            self.endpoints[ id2 ] = []
        self.endpoints[ id1 ].append( ids )
        self.endpoints[ id2 ].append( ids )

    def __delitem__( self, ids ):
        del self.ints[self.hash_ids(ids)]
        # Manage endpoint database:
        id1, id2 = ids[0], ids[-1]
        for i, ids_other in enumerate( self.endpoints[ id1 ] ):
            if ids == ids_other:
                del self.endpoints[ id1 ][i]
        for i, ids_other in enumerate( self.endpoints[ id2 ] ):
            if ids == ids_other:
                del self.endpoints[ id2 ][i]
        if not self.endpoints[ id1 ]:
            del self.endpoints[ id1 ]
        if not self.endpoints[ id2 ]:
            del self.endpoints[ id2 ]


    def __iter__( self ):
        return iter(self.ints)

    def itervalues( self ):
        return self.ints.itervalues()

    def __contains__( self, ids ):
        return self.hash_ids(ids) in self.ints

    def __len__( self ):
        return len(self.ints)

    def __nonzero__( self ):
        return bool( len( self.ints ) )

    def add( self, interval ):
        self[ interval.marker_ids ] = interval

    def intervals_with( self, identity ):
        if identity in self.endpoints:
            return [ self[ ids ] for ids in self.endpoints[ identity ] ]
        else:
            return []

    @property
    def total_weight( self ):
        answer = 0
        for interval in self.itervalues():
            answer += interval.weight
        return answer

# reads intervals from a file
# file_name - str: the name of the file to read from
# interval_list - list of Interval: the list to add to (Default = [])
# Return - list of Interval: the list of intervals read
def read_intervals_file(file_name, interval_list = []):
    file_stream = open(file_name)

    line = file_stream.readline()

    while len(line) > 0:
        trunc_line = line.strip()

        if len(trunc_line) > 0:
            if trunc_line[0] == '>':
                # read first interval
                interval, line = Interval.from_file(file_stream, trunc_line)
                if interval != None:
                    interval_list.append(interval)
                #endif

                # read the rest of the intervals
                while len(line) > 0:
                    interval, line = Interval.from_file(file_stream, line)
                    if interval != None:
                        interval_list.append(interval)
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

    return interval_list
#enddef


# writes intervals to a file
# file_name - str: the name of the file to write to
# hom_fam_list - list of Interval: the list to write (Default = [])
def write_intervals_file(file_name, interval_list):
    #assert type(file_name) is str
    #assert type(interval_list) is list

    file_stream = open(file_name, 'w')

    for interval in interval_list:
        file_stream.write(str(interval) + '\n')
    #endif

    file_stream.close()
#enddef

def write_intervals(log, ints, f):
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
#enddef