



# class for sorting
class PartitionTree:
    # initializes a PartitionTree
    # p - part.Partition, l - list of set
    def __init__(self, p, l):
        self._part = p        # partition
        self._list = l        # list of sets that make up the partition via refinement
    #enddef
#enddef

class OverlapComponent:
    def __init__(self, param):
        if type(param) == list:
            self._head = param[0]
            self._rows = param[1:]
            self._all_rows = param
            self._support = set([i for r in param for i in r._set])
        else:
            self._head = param
            self._rows = []
            self._all_rows = [param]
            self._support = param._set
        #endif
        
        self._nsupport = len(self._support)
    #enddef
    
    def add_row(self, row):
        self._rows.append(row)
        self._all_rows.append(row)
        self._support = self._support | row._set
        
        self._nsupport = len(self._support)
    #enddef
    
    def join_components(self, comp, pos):
        self._rows.extend(comp._all_rows)
        self._all_rows.extend(comp._all_rows)
        self._support = self._support | comp._support
        
        self._nsupport = len(self._support)
    #enddef
    
    def __lt__(self, a):
        if self._nsupport == a._nsupport:
            return min(self._support) < min(a._support)
        #endif    
        
        return self._nsupport < a._nsupport
    #enddef
    
    def __le__(self, a):
        if self._nsupport == a._nsupport:            
            return min(self._support) <= min(a._support)
        #endif
    
        return self._nsupport < a._nsupport
    #enddef
#endclass

# represents a Telomere
class Telomere:
    def __init__(self):
        pass
    #enddif

    def __str__(self):
        return 'T'
    #enddef
#endclass

class C1P:
    def __init__(self):
        #PartitionTree partition_tree
        #OvarlapComponent
        #Telomere
        self.c1p_linear = 0   # 1 for working with linear chromosomes
		self.c1p_circular = 1 # 1 for working with a unique circular chromosomes
		self.c1p_telomeres = 0 # CURRENTLY NOT IMPLEMENTED
		self.c1p_heuristic = 1 # Using a greedy heuristic
		self.c1p_bab = 0 # Using a branch-and-bound
		self.c1p_spectral = 0 # Using spectral seriation CURRENTLY NOT IMPLEMENTED
		self.markers_doubled = 0

	def run():
		## COMPUTING A C1P MATRIX ------------------------------------------------------------------------------
		output_tmp = "../output" # directory for intermediate files
		c1p_dir = output_tmp + "/C1P"		# intermediate files from C1P code
		cars_dir = output_tmp + "/CARS"
		pqr_tree = cars_dir + "/" + output_prefix + "PQRTREE"		# PQR-tree file
		if self.markers_doubled:
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
		if self.c1p_linear or self.c1p_telomeres > 0:
			if not quiet:
				print2("----> Creating PQR-tree: " + pqr_tree)
			#endif

			if self.markers_doubled:
				callprocess(["python", code_dir +"/C1P/C1P_compute_PQRtree.py", wacs, pqr_tree_doubled, params.output_ancestor])
		 		if not quiet:
		 			print2("----> Halving PQR-tree columns") 
		                #endif
				callprocess(["python", code_dir +"/C1P/C1P_halve_PQRtree.py", pqr_tree_doubled, pqr_tree])
			else:
				callprocess(["python", code_dir +"/C1P/C1P_compute_PQRtree.py", wacs, pqr_tree, params.output_ancestor])
			#endif
		elif self.c1p_circular:
			if not quiet:
				print2("----> Creating PQR-tree: " + pqr_tree)
			#endif

			if self.markers_doubled:
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
		if self.c1p_heuristic and self.c1p_telomeres == 0:
			if self.c1p_linear:
				do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "HEUR", "heuristic", "/C1P/C1P_make_C1P_heuristic.py", "/C1P/C1P_compute_PQRtree.py")
			elif self.c1p_circular:
				do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "HEUR", "heuristic", "/C1P/C1P_make_circC1P_heuristic.py", "/C1P/C1P_compute_PQCRtree.py", True)
			#endif
		#endif

		# branch and bound
		if params.c1p_bab and self.c1p_telomeres == 0:
			if self.c1p_linear:
				do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "BAB", "branch and bound", "/C1P/C1P_make_C1P_branch_and_bound.py", "/C1P/C1P_compute_PQRtree.py")
			elif self.c1p_circular:
				do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "BAB", "branch and bound", "/C1P/C1P_make_circC1P_branch_and_bound.py", "/C1P/C1P_compute_PQCRtree.py", True)
			#endif
		#endif

		# seriation
		if self.c1p_spectral:
			pq_tree = cars_dir + "/" + output_prefix + "PQTREE_SERIATION"	

			if params.acs_correction >= 2:
		 		if not quiet:
		 			print2("----> Computing PQ-tree from correlation matrix of a ternary matrix using parameter alpha") 
		                #endif
				callprocess(["python", code_dir +"/SERIATION/SERIATION_compute_PQtree_dotproduct_correlation_matrix.py", wacs, str(self.c1p_spectral_alpha), pq_tree, params.output_ancestor])
			else:
		 		if not quiet:
		 			print2("----> Computing PQ-tree using spectral seriation on correlation matrix") 
		                #endif
				callprocess(["python", code_dir +"/SERIATION/SERIATION_compute_PQtree_correlation_matrix.py", wacs, pq_tree, params.output_ancestor])

		        #endif
		#endif

		# telomere heuristic
		if self.c1p_telomeres == 1:
			do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "TEL_HEUR", "heuristic", "/C1P/C1P_make_mC1P_heuristic.py", "/C1P/C1P_compute_PQRtree.py")
		#endif

		# telomere branch and bound, add after
		if self.c1p_telomeres == 2:
			do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "TEL_BAB1", "branch and bound, added after", "/C1P/C1P_make_mC1P_branch_and_bound_both.py", "/C1P/C1P_compute_PQRtree.py")
		#endif

		#telomere branch and bound, add during
		if self.c1p_telomeres == 3:
			do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, "TEL_BAB2", "branch and bound, added during", "/C1P/C1P_make_mC1P_branch_and_bound.py", "/C1P/C1P_compute_PQRtree.py")
		#endif


	def do_c1p(params, c1p_dir, cars_dir, output_prefix, wacs, quiet, suffix, message, make_C1P, compute_PQRtree, m = False):
		acs_c1p	      = c1p_dir + "/" + output_prefix + "ACS_C1P_" + suffix	# C1P matrix file
		acs_discarded = c1p_dir + "/" + output_prefix + "ACS_DISC_" + suffix	# matrix of removed rows file
		pq_tree = cars_dir + "/" + output_prefix + "PQTREE_" + suffix # PQ-tree file	
		if self.markers_doubled:
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

		if self.markers_doubled:
			callprocess(["python", code_dir + compute_PQRtree, acs_c1p, pq_tree_doubled, params.output_ancestor])
			if not quiet:
				print2("----> Halving PQ-tree columns") 
			#endif
			callprocess(["python", code_dir + "/C1P/C1P_halve_PQRtree.py", pq_tree_doubled, pq_tree])
		else:
			callprocess(["python", code_dir + compute_PQRtree, acs_c1p, pq_tree, params.output_ancestor])
		#endif
	#enddef

    # compares partition trees
    def comp_part(a, b):
        return len(a._part._support) > len(b._part._support)
    #enddef

    # tests whether two sets overlap
    # a - set, b - set
    # return - bool
    def overlap(a, b):
        if len(a & b) == 0:
            return False
        #endif
        
        if a < b or b < a:
            return False
        #endif
        
        return True
    #enddef

    # auxiliary 'f' method for traverse
    # returns true if v._value and params['v']._value do not overlap
    # v - vertex, params = {'og': list of tree.SpanningTree, 'j': int,
    # 'v': tree.Vertex, 'found': int}
    # return - bool
    def test_overlap(v, params):
        if (overlap(v._value._set, params['v']._value._set)):
            if params['found'] >= 0:        # join trees
                params['og'][params['found']].join(params['og'][params['j']], v, params['v'])
                
                del params['og'][params['j']]
            else:        # add to tree
                params['og'][params['j']].add_vertex(params['v'], v)
            #endif
            
            return False
        #endif
            
        return True
    #enddef

    # auxiliary 'g' method
    # flips a tree around so e is the head with child v
    # v - vertex, e - vertex, b - bool, params = {'found': int}
    def flip_tree_node(v, e, b, params):
        if not b and params['found'] >= 0:
            v.remove_edge(e)
            e.add_edge(v)
        #endif
        
        return b
    #enddef

    # auxilary 'f' method for traverse
    # refines part with v._value
    # v - vertex, part - part.Partition
    # return - bool
    def refine(v, part):
        return part.refine(v._value)
    #enddef

    # auxiliary 'g' method for traverse
    # does nothing but passthrough b
    # v - unused, e - unused, b - A, params - unused
    # return - A
    def do_nothing(v, e, b, params):
        return b
    #endef

    # auxiliary 'f' method for traverse
    # adds the row, v._value, to matrix
    # v - vertex, matrix - bm.BinaryMatrix
    # return - bool; True
    def make_matrix(v, matrix):
        matrix.add_row_info(v._value)

        return True
    #endef

    # auxiliary 'f' method for traverse
    # adds the row, v._value, to the list params
    # v - vertex
    # params - list 
    #
    # return - bool; True
    def add_to_list(v, params):
        params.append(v._value)
        
        return True
    #enddef

    # traverses a tree with head v
    # calling f(v, params) before traversing v's children
    # and g(v, child, child_return_code, params) after traversing each child
    # if f and g return False the traversal exits
    # params is used as parameters for both f and g
    # v - tree.Vertex
    # f - boolean(tree.Vertex, params)
    # g - boolean(tree.Vertex, tree.Vertex, boolean, params)
    # params - anything
    #
    # return - bool
    def traverse(v, f, g, params):
        stack = [[f, None, v]]
        go = True

        while len(stack) > 0:
            job = stack.pop()
        
            if job[0] == f:        # prefix
                if go:
                    go = f(job[2], params)
                            
                    if job[1] != None:
                        stack.append([g, job[1], job[2]])
                    #endif
                    
                    if go: 
                        for e in job[2]._edges:
                            stack.append([f, job[2], e])
                        #endfor
                    #endif
                #endif
            else:        # postfix
                g(job[1], job[2], go, params)
            #endif
        #endwhile
        
        return go
    #enddef

    # mem, sort
    # tests if the row, row, can be added to the sets of partitions p and
    # still be C1P
    # row - set, p - list of PartitionTree, mem - list of rows
    # return - bool
    def test_row(row, p):
        found = False        # True if first intersecting part is found
        failed = False        # True if not C1P
        j = 0        # current index
        mem = [[], []]        # altered rows
            
        # find overlaps
        while j < len(p):
            if len(p[j]._part._support & row._set) > 0:
                overlaps = False        # True if current component overlaps with row
            
                for s in p[j]._list:
                    if overlap(row._set, s._set):
                        overlaps = True
                        
                        break
                    #endif
                #endfor
                
                if overlaps:
                            
                    # refine the new found partition
                    if not found:
                        pNew = PartitionTree(p[j]._part.copy(), copy.copy(p[j]._list))
                        found = True
                        failed = not pNew._part.refine(row)
                            
                        if not failed:
                            pNew._list.append(row)
                            
                            mem[0].append(p[j])
                            mem[1].append(pNew)
                        #endif
                    else:        # if partition already found join the two partitions
                        # error (this should not happen)
                        if len(pNew._part._part._value) == 0:
                            print 'error1'
                            print row
                        
                            return
                        #endif
                        
                        pCopy = PartitionTree(p[j]._part.copy(), copy.copy(p[j]._list))
                        failed = not pNew._part.join(pCopy._part)
                            
                        if not failed:
                            pNew._list.extend(p[j]._list)
                            mem[0].append(p[j])
                        #endif
                    #endif
                    
                    # error (this should not happen)
                    if len(pNew._part._part._value) == 0:
                        print 'error2'
                        print row
                        
                        return
                    #endif
                    
                    if failed:
                        break
                    #endif
                #endif
            j = j + 1
        #endwhile
            
        # create new partition (new overlap component)
        if not found:
            par =  PartitionTree(part.Partition(row), [row])
        
            sort.insert(p, par, 0, len(p), comp_part)
            mem = [[], [par]]
        elif not failed:                    
            # remove refined rows
            for m in mem[0]:
                p.remove(m)
            #endfor

            # add refined row
            sort.insert(p, pNew, 0, found, comp_part)
        else:
            mem = [[], []]
        #endif
            
        return not failed, mem
    #endef

    # refines an overlapping matrix in the current order of rows
    # m - bm.BinaryMatrix
    # return - bool
    def refine_matrix(m):
        p = part.Partition(m.get_row_info(0))
        C1P = True
            
        for i in xrange(1, m._height):        
            if not p.refine(m.get_row_info(i)):
                return False                    
            #endif
        #endfor
        
        return True
    #enddef

    # creates the overlap graph (as spanning trees) for a given matrix, matrix
    #
    # ALGORITHM PREMISE:
    # OG = set of trees with vertex set a subset of the rows of matrix.
    # At any step in the alogithm each tree in OG will have the property that each vertex overlaps with each of its children.
    # At the end of the algorithm OG will be a set of spanning trees of the connect components of the overlap graph of matrix
    # and each tree's infix order will give a partitio refinement order.
    # 
    # matrix - bm.BinaryMatrix
    #
    # return - list of SpanningTree
    def create_overlap_graph(matrix):
        og = []        # overlap graph as a list of tree.SpanningTree

        for i in xrange(matrix._height):
            row = matrix.get_row_info(i)        # current row
            f = -1        # index of overlap graph found for row (-1 = undefined)
            v = tree.Vertex(row)        # tree.Vertex for row to
                                        # be put in spanning tree
            j = 0            # index of overlap graph
                
            # try to add to existing graph
            while j < len(og):        
                params = {'og': og, 'j': j, 'v': v, 'found': f}
                
                if len(og[j]._support & row._set) > 0 and not traverse(og[j]._head, test_overlap, flip_tree_node, params):
                    if f < 0:
                        f = j
                        j += 1        
                    #endif
                else:
                    j += 1
                #endif
            #endfor
            
            # make new spanning tree
            if f < 0:
                og.append(tree.SpanningTree(v))
            #endif
        #endfor
        
        rows = []
        
        for g in og:
            rows.append([])
        
            traverse(g._head, add_to_list, do_nothing, rows[-1])
        #endfor
        
        return [OverlapComponent(r) for r in rows]
    #endfor

    # checks whether the binary matrix, matrix, has the consectutive-ones property 
    # by using an algorithm akin to the algorithm of Section 5 of (McConell 2004)
    # matrix - bm.BinaryMatrix
    # return - bool
    def check_C1P(matrix):
        # create overlap graphs
        og = create_overlap_graph(matrix)    # overlap graph (as a spanning tree)
                    
        # traverse overlap graph and refine partitions
        for j in xrange(len(og)):
            p = part.Partition(og[j]._head)        # current partition to refine
                    
            for v in og[j]._rows:
                if not p.refine(v):
                    return False
                #endif
            #endfor
        #endfor
        
        return True
    #enddef

    def make_intersect_components(matrix):
        mats = []
        sups = []
        
        for r in matrix._rows:
            found = False
        
            for i in xrange(len(mats)):
                if len(sups[i] & r._set) > 0:
                    mats[i].add_row_info(r)
                    
                    sups[i] = sups[i] | r._set
                    
                    found = True
                    
                    for k in xrange(len(mats) - 1, i, -1):
                        if len(sups[k] & r._set) > 0:
                            for row in mats[k]._rows:
                                mats[i].add_row_info(row)
                            #endfor
                            
                            sups[i] = sups[i] | sups[k]
                            
                            del sups[k]
                            del mats[k]
                        #endif
                    #endfor
                    
                    break
                #endif
            #endfor
            
            if not found:
                mats.append(bm.BinaryMatrix())

                mats[-1].add_row_info(r)
                
                sups.append(r._set)
            #endif
        #endfor
        
        return mats
    #enddef

    # takes a matrix, m, and splits it into matrices that each only have one 
    # connected component
    # m - BinaryMatrix
    # return list of bm.BinaryMatrix
    def split_matrix(m):
        # create overlap graph
        g = create_overlap_graph(m)        # overlap graph
        ms = []        # split matrices
        
        for c in g:
            # make a new matrix for each overlap graph
            n = bm.BinaryMatrix()        # current overlap component
            
            for r in c._all_rows:
                n.add_row_info(r)
            #endfor
            
            ms.append(n)
        #endfor
        
        return ms
    #enddef

    # removes rows until the matrix is C1P
    # matrix - bm.BinaryMatrix
    # return - list of set
    def make_C1P(matrix):
        p = []        # list of partitions and spanning trees used
        i = 0        # iterator
        rows = []        # rows to remove to make C1P

        # sort the matrix by weight
        matrix.sort()
            
        # process each row
        while i < matrix._height:    
            row = matrix.get_row_info(i)        # current row
            yes, temp =  test_row(row, p)        # true if test_row succeeds, aux
            
            if not yes:
                # remove row from matrix
    #            print row._set
                            
                rows.append(i)
            #endif
            
            i += 1
        #endwhile
            
        return rows
    #enddef
  