from collections import deque

from data_structures import genomes
from data_structures import intervals
from data_structures import markers

DEBUG = False

# Function to decide if a genome instance is realizable, when the genome only
# has adjacencies and repeat spanning intervals, under the given genome model.
# Assumes that markers are doubled.
# Arguments:
#    hom_fams: list of HomFam objects - the markers to process
#    adjacencies: IntervalDict
#    RSIs: IntervalDict
#   genome_model: string - either "linear" or "mixed"
# Output:
#    Boolean - is genome realiable?
# NOTE: in the future there should be an argument for intervals containing only
#       unique markers.
def decision_RSIs( hom_fams, adjacencies, RSIs, genome_model, debug=None ):
    realizable = True
    # Create simplified copied of the arguments, for processing:
    # A set of all markers, a set of the unique markers, a set of the
    # unoriented markers and a dict for multiplicity. The 'unoriented' markers
    # is a subset of V_1. The markers in V_1 - unoriented are in fact oriented
    # marker extremities, so they technically don't have an orientation by
    # themselves.
    #NOTE: for backtracking, these datastructures should already be there
    #       (partially), only to be augmented with single additional RSI.
    V_1 = set([])
    multiplicity = {}
    unoriented = set([])
    for hom_fam in hom_fams:
        V_1.add( hom_fam.id )
        multiplicity[ hom_fam.id ] = hom_fam.copy_number

    # Create a simple, symmetric copy of adjacencies for processing, in the
    # form of a dictionary (keyed by IDs) of sets of IDs
    E_1 = {}
    for adj in adjacencies.itervalues():
        id1 = adj.marker_ids[0]
        id2 = adj.marker_ids[-1]
        add_edge( E_1, id1, id2 )

    # First stage:
    # Loop over RSIs, for each RSI create new unique unoriented markers and
    # edges.
    #NOTE: for backtracking, this should only be done for one RSI.
    for RSI in RSIs.itervalues():
        previous = RSI.marker_ids[0]
        # In order to give each marker in the repeat spanning interval a
        # unique ID, enumerate them and use index in ID. Need to make sure that
        # head/tail pair have the same ID except for the _h or _t suffix, this
        # is why we use 'i/2' as part of the new ID (which is rounded down to
        # and integer).
        for i, r_id in enumerate( RSI.marker_ids[1:-1] ):
            new_id = "RSI_internal_" + RSI.id + "_" + str(i/2) + "_" + r_id
            V_1.add( new_id )
            multiplicity[ new_id ] = 1
            # Add edge to previous marker in RSI if not siblings
            if not markers.are_siblings( new_id, previous ):
                add_edge( E_1, previous, new_id )
            # Decrease multiplicity, check if still non-negative:
            #NOTE: may turn repeat into unique marker.
            multiplicity[ r_id ] -= 1
            if multiplicity[ r_id ] < 0:
                realizable = False
            previous = new_id
        if not realizable: break
        # Handle endpoint adjacencies
        add_edge( E_1, RSI.marker_ids[-1], previous )
        remove_edge( E_1, RSI.marker_ids[0], RSI.marker_ids[1] )
        remove_edge( E_1, RSI.marker_ids[-1], RSI.marker_ids[-2] )

    if not realizable: return False

    # To finish off, loop over the repeats, and remove the ones who's
    # multiplicity has become zero.
    to_remove = set([])
    for marker in V_1:
        if multiplicity[ marker ] == 0:
            edges_to_remove = []
            if marker in E_1:
                for m in E_1[ marker ]:
                    edges_to_remove.append( m )
            for m in edges_to_remove:
                remove_edge( E_1, m, marker )
            to_remove.add( marker )
    V_1.difference_update( to_remove )

    # Second stage:
    # Create shallow copies of V_1, E_1 (for clarity)
    V_2 = V_1
    E_2 = E_1
    #NOTE: in the future, if we need to keep the data structures and pass them
    #      on, we would need to perform deep copies of V_1, E_1, multiplicity
    #      and unoriented. Then, only work with affected RSIs.

    # Do a DFS of the unique markers to find connected components, at the same
    # time check that the connected component will be realizable.
    # The connected component consists only of unique markers with adjacencies.
    to_explore = set( filter( lambda m: multiplicity[m] == 1, V_2 ) )
    search_queue = deque([])
    current_component = set([])
    current_frontier = []
    current_component_edges = {}
    while to_explore:

        # Add a unique marker to the queue
        search_queue.appendleft( next( iter( to_explore ) ) )

        while search_queue:
            current_vertex = search_queue.pop()
            to_explore.remove( current_vertex )
            # Get the neighbourhood of the current vertex, including "sibling"
            # of an oriented vertex.
            neighbours = set([])
            if current_vertex in E_2:
                neighbours = E_2[ current_vertex ]
            if not current_vertex in unoriented:
                other_half = markers.sibling( current_vertex )
                neighbours.add( other_half )
            unique_neighbours = filter( lambda m: multiplicity[m] == 1,
                                        neighbours )
            repeat_neighbours = filter( lambda m: multiplicity[m] > 1,
                                        neighbours )
            # Process neighbourhood.
            for neighbour in neighbours:
                if not neighbour == markers.sibling( current_vertex ):
                    add_edge( current_component_edges,
                              current_vertex,
                              neighbour )
            for neighbour in unique_neighbours:
                if neighbour in to_explore and not neighbour in search_queue:
                    search_queue.appendleft( neighbour )
            current_frontier += repeat_neighbours
            current_component.add( current_vertex )

        #NOTE: should use C1P in the future, when we have intervals as well.
        # Check if the current connected component is realizable.
        if not decision_adjacencies_internal(
                    current_component,
                    current_component_edges,
                    multiplicity,
                    unoriented,
                    genome_model ):
            realizable = False
            if debug:
                debug.write( "decision_RSIs: Unrealizable after connected "
                             "component adjacency check.\n" )
            break
        # Collapse the current connected component to a single vertex.
        V_2.difference_update( current_component )
        for key1,edges in current_component_edges.iteritems():
            for key2 in edges:
                remove_edge( E_2, key1, key2 )
        new_vertex = "connected_component_" + next( iter( current_component ) )
        V_2.add( new_vertex )
        multiplicity[ new_vertex ] = 1
        unoriented.add( new_vertex )
        for repeat in current_frontier:
            add_edge( E_2, new_vertex, repeat )

    if not realizable: return False

    # Check the realizability of the resulting instance: unique vertices should
    # now all have degree two at most, repeats no more than their multiplicity.
    realizable = decision_adjacencies_internal(
        V_2,
        adjacencies,
        multiplicity,
        unoriented,
        genome_model,
        debug
        )

    if not realizable:
        if debug:
            debug.write( "decision_RSIs: Unrealizable after collapsed "
                         "connected components adjacencency check.\n" )
        return False

    return True

# Function to decide if an instance of markers and only adjacencies is
# realizable under the given genome model.
# Intended 'internal' since it relies on data structures used only in this
# module.
# Arguments:
#   markers: set of strings - the markers to check.
#   adjacencies: dict of sets of strings - keyed by and contains marker IDs
#                 denoting adjacencies.
#   multiplicity: dict of ints - multiplicity of markers
#   unordered: set of strings - the markers that are unoriented (so not
#               doubled)
#   genome_model: string - either "linear" or "mixed"
# Output:
#   boolean, is the instance realizable?
def decision_adjacencies_internal( markers,
                                   adjacencies,
                                   multiplicity,
                                   unoriented,
                                   genome_model,
                                   debug=None ):
    realizable = True
    degree_sum = 0
    multiplicity_sum = 0
    for marker in markers:
        degree = 0
        if marker in adjacencies:
            degree = len( adjacencies[ marker ] )
        degree_sum += degree
        # For unoriented markers, degree can be at most 2 times multiplicity.
        # For oriented (doubled) markers, degree can be at most multiplicity.
        if marker in unoriented:
            if degree > 2 * multiplicity[ marker ]:
                realizable = False
                if debug:
                    debug.write( "decision_adjacencies_internal: Unrealizable,"
                    " unoriented marker \'%s\' with multiplicity %s has"
                    " degree %s, neighbours %s.\n"
                        % ( marker, multiplicity[marker],
                            degree, adjacencies[marker] ) )
                break
            multiplicity_sum += 2 * multiplicity[ marker ]
        else:
            if degree > multiplicity[ marker ]:
                realizable = False
                if debug:
                    debug.write( "decision_adjacencies_internal: Unrealizable,"
                    " oriented marker \'%s\' with multiplicity %s has degree"
                    " %s, neighbours %s.\n"
                        % ( marker, multiplicity[marker],
                            degree, adjacencies[marker] ) )
                break
            multiplicity_sum += multiplicity[ marker ]

    if genome_model == "linear":
        if degree_sum >= multiplicity_sum:
            realizable = False
    elif genome_model == "mixed":
        if degree_sum > multiplicity_sum:
            realizable = False
    else:
        print """ERROR (decision_adjacencies_internal):
        Unknown genome model!"""

    return realizable


# Auxiliary functions

# Simple function to add edge to symmetric nested edge dictionary:
# Arguments:
#    E: the nested edge dictionary
#    key1, key2: the 'vertices' to link with edges
def add_edge( E, key1, key2 ):
    if not key1 in E:
        E[ key1 ] = set([])
    if not key2 in E:
        E[ key2 ] = set([])
    E[ key1 ].add( key2 )
    E[ key2 ].add( key1 )


# Simple function to remove     edge to symmetric nested edge dictionary.
# Arguments:
#    E: the nested edge dictionary
#    key1, key2: the 'vertices' between which to remove an edge.
def remove_edge( E, key1, key2 ):
    if key1 in E:
        if key2 in E[ key1 ]:
            E[ key1 ].remove( key2 )
        if not    E[ key1 ]:
            del E[ key1 ]
    if key2 in E:
        if key1 in E[ key2 ]:
            E[ key2 ].remove( key1 )
        if not E[ key2 ]:
            del E[ key2 ]
