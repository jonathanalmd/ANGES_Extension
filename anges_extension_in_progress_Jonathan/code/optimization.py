from data_structures import intervals

import decisions
import networkx


# Function to find maximal subset (w.r.t weight) of adjacencies between given
# markers such that this subset is realizable under the mixed genome model.
# Arguments:
#   hom_fams: list of HomFam objects - the doubled markers to consider.
#   adjacencies: IntervalDict - the adjacencies to optimize
# Output:
#   max_adjacencies: IntervalDict - the maximal subset of adjacencies that are
#                                   realizable.
def opt_adjacencies( hom_fams, adjacencies ):


    #import pdb; pdb.set_trace()

    # We need to find a 2m-matching of the markers, but since we assume that we
    # are working with doubled markers, in this case we need an m-matching.
    # Create a dictionary that maps markers to their multiplicity:
    multiplicity = {}
    # multiplicity_high = []
    for hom_fam in hom_fams:
        multiplicity[ hom_fam.id ] = hom_fam.copy_number
        # if hom_fam.copy_number > 1:
            # multiplicity_high.append(hom_fam.id)
    # print multiplicity_high

    # First, create a networkx graph to encode the markers and adjacencies in
    # the arguments, to make further processing easier.
    G_0 = networkx.Graph()
    G_0.add_nodes_from( [ hom_fam.id for hom_fam in hom_fams ] )
    for adjacency in adjacencies.itervalues():
        id1, id2 = adjacency.marker_ids[0], adjacency.marker_ids[-1]
        G_0.add_edge( id1, id2, weight=adjacency.weight, adj=adjacency )
        # if id1 in multiplicity_high:
        #     if id2 in multiplicity_high:
        #         G_9.add_edge(id1,id2, weight=adjacency.weight, adj=adjacency)
        #     else:
        #         G_9.add_edge(id1,id1, weight=-1, adj=adjacency)
        # elif id2 in multiplicity_high:
        #     G_9.add_edge(id2,id2, weight=-1,adj=adjacency)
    # G_9.add_edge('179_t','16_t',weight=-1)
    # edges_list = G_9.edges()
    # print edges_list

    
    
    # Make a new networkx graph, and create structure required for the matching
    # algorithm to work.
    G = networkx.Graph()
    for m in G_0:
        # Create m 'marker_vertices' for a marker of multiplicity m
        marker_vertices = to_marker_vertices( m, multiplicity[ m ] )
        # Create d 'edge_vertices' for the same marker (of degree d)
        edge_vertices_weights = [
            ( to_edge_vertex( m, neighbour ),
              G_0[ m ][ neighbour ][ 'weight' ] ) \
            for neighbour in G_0[ m ]
            ]
        G.add_nodes_from( marker_vertices )
        G.add_nodes_from( [ v for (v,_) in edge_vertices_weights ] )
        # Create edges between the new vertice, construct complete bipartite
        # graph on marker_vertices, edge_vertices. Use the weights of the
        # edge_vertices.
        for edge_vertex, w in edge_vertices_weights:
            for marker_vertex in marker_vertices:
                G.add_edge( edge_vertex, marker_vertex, weight = w )
    # Add edges between the appropriate 'edge_vertices'
    for m1, m2, attr in G_0.edges( data = True ):
        G.add_edge( to_edge_vertex( m1, m2 ),
                    to_edge_vertex( m2, m1 ),
                    weight = attr[ 'weight' ] )

    # Find optimal set of adjacencies.
    # The matching is a dictionary such that dict[v]==u iff the nodes v,u in G
    # are connected in the matching.
    matching = networkx.max_weight_matching( G )

    # Translate the matching solution back to the m-matching we want on the
    # given markers.
    G_8 = networkx.Graph()
    max_adjacencies = intervals.IntervalDict()
    for u, v in matching.iteritems():
        # If an edge in G_0 is part of the m-matching, it will be 'caught' four
        # times in this loop, namely, twice (once for each direction) for each
        # of the two marker_vertex / edge_vertex pairs in G that represent one
        # edge in G_0. We only record the adjacency once.
        # Two cases: u,v form a vertex/edge pair or u,v are both edge_vertices.
        # Only the first case indicates an edge that is part of the m-matching
        # of G_0.
        if is_marker_vertex( u ) and not is_marker_vertex( v ):
            m1, m2 = from_edge_vertex( v )
            # Check that the adjacency between m1 and m2 is in the matching at
            # the other marker as well.
            marker_vertices_2 = to_marker_vertices( m2, multiplicity[ m2 ] )
            edge_vertex_2 = to_edge_vertex( m2, m1 )
            adj_at_other_marker = (
                edge_vertex_2 in matching and
                any( [ matching[ edge_vertex_2 ] == m
                     for m in marker_vertices_2 ] )
                )
            if edge_vertex_2 in matching and adj_at_other_marker:
                adjacency = G_0[ m1 ][ m2 ][ 'adj' ]
                max_adjacencies.add( adjacency )
                if multiplicity[m1] > 1:
                    if multiplicity[m2] > 1:
                        G_8.add_edge(m1,m2, weight=adjacency.weight, adj= adjacency)
                    else:
                        G_8.add_edge(m1,m1, weight=1, adj= adjacency)
                elif multiplicity[m2] > 1:
                    G_8.add_edge(m2,m2,weight=1, adj= adjacency)

    high_cp_graph = [c for c in sorted(networkx.connected_components(G_8), key=len, reverse=True)]
    # print high_cp_graph
    
    rc_total_list = []
    rc_total_list_int = []
    dont_print = []
    for rc in high_cp_graph:
        rc_elements = ""
        rc_elements_int = []
        for rc_elem in rc:
            if int(rc_elem[:-2]) not in dont_print:
                dont_print.append(int(rc_elem[:-2]))
                rc_elements = rc_elements + (rc_elem[:-2]) + " "
                rc_elements_int.append(int(rc_elem[:-2]))
        if rc_elements:
            rc_total_list.append(rc_elements[:-1])
            rc_total_list_int.append(rc_elements_int)

    return max_adjacencies, rc_total_list, rc_total_list_int


# Function to find maximal subset (w.r.t weight) of RSIs between given markers,
# and with a given realizable set of adjacencies, such that this subset is
# realizable under the given genome model.
# Arguments:
#   hom_fams: list of HomFam objects - the doubled markers to consider.
#   adjacencies: IntervalDict - the adjacencies between markers.
#   RSIs: IntervalDict - the RSIs to optimize over.
# Output:
#   realizable_RSIs: IntervalDict - the RSIs that are realizable.
def opt_RSIs_greedy( hom_fams, adjacencies, RSIs, genome_model, debug=None ):
    realizable_RSIs = intervals.IntervalDict()
    # Sort the given RSIs in reversed order by weight.
    sorted_RSIs = sorted( list( RSIs.itervalues() ),
                          key = lambda RSI: - RSI.weight )

    # Main greedy loop:
    for RSI in sorted_RSIs:
        # Add new RSI to the realizable RSIs.
        realizable_RSIs.add( RSI )
        # Check if the genome is still realizable.
        if debug:
            debug.write( "opt_RSIs_greedy: added RSI %s.\n" %RSI.marker_ids )
        if not decisions.decision_RSIs(
                hom_fams,
                adjacencies,
                realizable_RSIs,
                genome_model,
                debug ):
            del realizable_RSIs[ RSI.marker_ids ]
            if debug:
                debug.write( "opt_RSIs_greedy: removed last RSI.\n" )

    return realizable_RSIs





# Auxiliary functions for adjacency optimization

# Functions to map a pair of vertices to an edge_vertex and back.
# The edge_vertex will be at the first marker given/returned.
def to_edge_vertex( u, v ):
    return "edgevertex_" + u + "_to_" + v
def from_edge_vertex( e_v ):
    pair = e_v.split( "_to_" )
    u = pair[0][11:]
    v = pair[1]
    return u, v

# Functions to map a vertex to a marker_vertex and back.
def to_marker_vertices( v, n ):
    return [ "markervertex_" + v + "_" + str(i) for i in range( n ) ]
def from_marker_vertex( v ):
    return v[13:-2]

# Function to check if a vertex is a marker_vertex
def is_marker_vertex( v ):
    if v[:13] == "markervertex_":
        return True
    else:
        return False
