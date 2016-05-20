from data_structures import markers
from data_structures import intervals

from collections import deque


# Function to assembly a realizable genome, given a set of adjacencies and RSIs.
# Arguments:
#   hom_fams: list of HomFam objects.
#   adjacencies: IntervalDict
#   RSIs: IntervalDict
#   species: String - the species name of the assembled genome.
# Output:
#   hom_fams: a list of HomFam objects that define the assembled genome. The
#             markers are oriented and not doubled.
def assemble( hom_fams, adjacencies, RSIs, species ):
    # Preparation:
    # Create dict of new hom_fams (keyed by IDs) for easy access while adding
    # loci. To be converted to list at the end of this algorithm.
    hom_fam_dict = {}
    # Keep track of the multiplicity of the markers.
    multiplicity = {}
    for hom_fam in hom_fams:
        multiplicity[ hom_fam.id ] = hom_fam.copy_number
    # Create an IntervalDict for adjacencies to be explored, that only contain
    # unique markers.
    adjs_to_explore = intervals.IntervalDict()
    for adj in adjacencies.itervalues():
        if ( multiplicity[ adj.marker_ids[0] ] == 1 and
                multiplicity[ adj.marker_ids[-1] ] == 1 ):
            adjs_to_explore.add( adj )
    # Create a copy of the RSIs IntervalDict, so we can also keep track of
    # which RSIs have been used.
    RSIs_to_explore = intervals.IntervalDict()
    for RSI in RSIs.itervalues():
        RSIs_to_explore.add( RSI )

    # Main loop:
    current_chromosome = 1
    while adjs_to_explore or RSIs_to_explore:
        # Prepare to assemble a single CAR.
        Q = deque()
        if adjs_to_explore:
            adj = next( adjs_to_explore.itervalues() )
            Q.appendleft( adj )
            del adjs_to_explore[ adj.marker_ids ]
        else:
            RSI = next( RSIs_to_explore.itervalues() )
            Q.appendleft( RSI )
            del RSIs_to_explore[ RSI.marker_ids ]
        CAR = deque()
        while Q:
            I = Q.pop()
            # Take L to be the endpoint markers in I that are not already in
            # CAR. If an endpoint marker in I is already in CAR, it has to be
            # at one of the extremes.
            id1, id2 = I.marker_ids[0], I.marker_ids[-1]
            if CAR:
                L = [ m for m in [ id1, id2 ] if m != CAR[0] and m != CAR[-1] ]
            else:
                L = [ id1, id2 ]
            for endpoint in L:
                sibling = markers.sibling( endpoint )
                other_end = id2 if endpoint == id1 else id1
                # Step 1: extend the CAR.
                if CAR and other_end == CAR[0]:
                    if other_end == id1:
                        for m in I.marker_ids[1:]:
                            CAR.appendleft( m )
                    else:
                        for m in reversed( I.marker_ids[:-1] ):
                            CAR.appendleft( m )
                    CAR.appendleft( sibling )
                elif CAR and other_end == CAR[-1]:
                    if other_end == id1:
                        for m in I.marker_ids[1:]:
                            CAR.append( m )
                    else:
                        for m in reversed( I.marker_ids[:-1] ):
                            CAR.append( m )
                    CAR.append( sibling )
                else:
                    # Assert: CAR is empty
                    CAR.append( endpoint )
                    CAR.append( sibling )

                # Step 2: attempt to add another interval to Q.
                next_adjs = adjs_to_explore.intervals_with( sibling )
                if next_adjs:
                    # Assert: len( next_adjs ) == 1, or the genome wouldn't be
                    # realizable.
                    Q.appendleft( next_adjs[0] )
                    del adjs_to_explore[ next_adjs[0].marker_ids ]
                next_RSIs = RSIs_to_explore.intervals_with( sibling )
                if next_RSIs:
                    # Assert: len( next_RSIs ) == 1, or the genome wouldn't be
                    # realizable.
                    Q.appendleft( next_RSIs[0] )
                    del RSIs_to_explore[ next_RSIs[0].marker_ids ]
        # Finalize the construction of the current chromosome. Undouble the
        # markers on the spot.
        current_position = 0
        CAR_list = list( CAR )
        for i, marker in enumerate( CAR_list[:-1] ):
            # Catch all head / tail pairs exactly once.
            sibling = markers.sibling( marker )
            if sibling != CAR_list[ i + 1 ]:
                continue
            undoubled = markers.from_doubled( marker )
            if not undoubled in hom_fam_dict:
                hom_fam_dict[ undoubled ] = markers.HomFam(
                    undoubled,                     # ID
                    [],                         # Loci
                    multiplicity[ marker ],     # Copy number
                    ''                          # Comment
                    )
            if ( marker, sibling ) == markers.to_doubled( undoubled ):
                orientation = 1
            else:
                orientation = -1
            #NOTE: Currently, this doesn't not respect original marker length.
            locus = markers.Locus(
                species,                        # Species
                str( current_chromosome ),      # Chromosome
                current_position,               # Start
                current_position + 1,           # End
                orientation,                    # Orientation
                ''                              # Comment
                )
            hom_fam_dict[ undoubled ].loci.append( locus )
            current_position += 2
        current_chromosome += 1

    hom_fams = hom_fam_dict.values()
    return hom_fams
