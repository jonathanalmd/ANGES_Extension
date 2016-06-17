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

import copy
import shutil
import time

import bm
import c1p
import pqtree

import mc1p
import bab
import tree
import sort
from babtester import *
import math
import string

class MasterC1P:
    def __init__(self):
        self.markers_input = ""
        self.species_tree = ""
        self.acs_file = ""
        self.output_dir = ""
        self.output_ancestor = ""

        self.code_dir = ""
        self.pqr_tree_doubled = ""
        self.c1p_dir = ""        # intermediate files from C1P code
        self.cars_dir = ""
        self.pqr_tree = ""
        self.output_prefix = ""

        self.markers_doubled = 0
        self.markers_unique = 0

        self.c1p_circular = 0
        self.c1p_linear = 0
        self.acs_ra = 0
        self.acs_correction = 0
        self.c1p_spectral = 0
        self.acs_weight = 0
        self.acs_aci = 0
        self.acs_sci = 0
        self.acs_mci = 0
        self.c1p_telomeres = 0
        self.c1p_heuristic = 0
        self.c1p_bab = 0
        self.c1p_spectral_alpha = 0

        self.quiet = False

        self.markers_provided = False
        self.species_tree_provided = False
        self.acs_file_provided = False
        self.acs_pairs_provided = False

    def setConfigParams(self, config_file):
        try:
            config = {}
            execfile(config_file, config)
            #collect the information from config file
            self.markers_input      = config["homologous_families"]
            self.species_tree       = config["tree_file"]
            self.acs_file           = config["acs_file"]
            self.output_dir         = config["output_c1p"]
            self.output_ancestor    = config["output_ancestor"]

            self.markers_doubled    = config["markers_doubled"]
            self.markers_unique     = config["markers_unique"]

            self.c1p_circular       = config["c1p_circular"]
            self.c1p_linear         = config["c1p_linear"]
            self.acs_ra             = config["acs_ra"]
            self.acs_correction     = config["acs_correction"]
            self.c1p_spectral       = config["c1p_spectral"]
            self.acs_weight         = config["acs_weight"]
            self.acs_aci            = config["acs_aci"]
            self.acs_sci            = config["acs_sci"]
            self.acs_mci            = config["acs_mci"]
            self.c1p_telomeres      = config["c1p_telomeres"]
            self.c1p_heuristic      = config["c1p_heuristic"]
            self.c1p_bab            = config["c1p_bab"]
            self.c1p_spectral_alpha = config["c1p_spectral_alpha"]


            self.output_prefix = self.output_ancestor + "_"         # Prefix of all filenames
            self.markers_provided      = config["markers_provided"]
            self.species_tree_provided = config["species_tree_provided"]
            self.acs_file_provided     = config["acs_file_provided"]
            self.acs_pairs_provided    = config["acs_pairs_provided"]


            self.c1p_dir = self.output_dir + "/C1P"        # intermediate files from C1P code
            self.cars_dir = self.output_dir + "/CARS"
            self.pqr_tree = self.cars_dir + "/" + self.output_prefix + "PQRTREE"        # PQR-tree file

        except IOError:
            print("{}  ERROR (master.py -> process.py) - could not open configuration file: {}\n"
                    .format(strtime(),config_file))
            sys.exit()
        config.clear()
    #enddef

    def checkErrors(self):
          # Missing input
        if not self.markers_provided:
            print("ERROR: no markers file")
            sys.exit(-1)
        #endif
        if not self.species_tree_provided:
            print("ERROR: no species tree file")
            sys.exit(-1)
        #endif

        # Conflicting models
        if self.c1p_circular and self.c1p_linear:
            print("ERROR: chose either circular or linear C1P")
            sys.exit(-1)
        #endif

        if (not self.markers_unique) and (self.acs_ra):
            print("ERROR: reliable adjacencies are not defined for repeated markers")
            sys.exit(-1)
        #endif

        if self.acs_correction >= 2 and (not self.c1p_spectral):
            print("ERROR: corrected ACS with Xs require using the spectral seriation algorithm")
            sys.exit(-1)
        #endif    
        #if self.acs_correction < 2 and self.c1p_sandwich:
        #    print("ERROR: the sandwich c1p requires using corrected ACS with X's ")
        #    sys.exit(-1)
        #endif    


        ## Current limitations

        if self.acs_ra and (not self.acs_sa):
            print("ERROR: computing reliable adjacencies require computing supported adjacencies")
            # could be addressed by forcing self.acs_sa to be true if self.acs_ra is true
            sys.exit(-1)
        #endif

        if self.acs_weight!=1 and not self.acs_file_provided:
            print("ERROR: weighting ACS by linear interpolation is mandatory currently ")
            sys.exit(-1)
        #endif    


        if self.c1p_circular:
            if self.acs_aci or self.acs_sci or self.acs_mci:
                print("ERROR: circular chromosomes can only be computed from adjacencies")
                sys.exit(-1)
                #endif
        #endif

        if self.c1p_spectral:
            if self.markers_doubled:
                print("ERROR: spectral seriation can not be used with doubled markers")
                sys.exit(-1)
                #endif
            if self.c1p_circular:
                print("ERROR: spectral seriation can not be used with circular chromosomes")
                sys.exit(-1)
                #endif
            if self.c1p_telomeres:
                print("ERROR: spectral seriation can not be used with telomeric ACS")
                sys.exit(-1)
            #endif
        #endif

    # self.do_c1p("HEUR", "heuristic", "/C1P/C1P_make_C1P_heuristic.py", "/C1P/C1P_compute_PQRtree.py")

    def do_c1p(self, suffix, message, make_C1P, compute_PQRtree, m = False):
        acs_c1p          = self.c1p_dir + "/" + self.output_prefix + "ACS_C1P_" + suffix    # C1P matrix file
        acs_discarded = self.c1p_dir + "/" + self.output_prefix + "ACS_DISC_" + suffix    # matrix of removed rows file
        pq_tree = self.cars_dir + "/" + self.output_prefix + "PQTREE_" + suffix # PQ-tree file    
        if self.markers_doubled:
            pq_tree_doubled = self.cars_dir + "/" + self.output_prefix + "PQTREE_DOUBLED_" + suffix # PQ-tree file    
        #endif
        if not self.quiet:
            print("----> Making (weighted)ACS matrix C1P (" + message + "): " + self.acs_file + " " + acs_discarded)
        #endif

        if m:
            if make_C1P == "heuristic":
                m = bm.BinaryMatrix()       # matrix
                mat_rem = bm.BinaryMatrix() # rows discarded
                                        
                m.from_file(self.acs_file)

                rows = c1p.make_C1P(m)      # rows to remove to make C1P
                            
                for j in xrange(len(rows) - 1, -1, -1):
                    mat_rem.add_row_info(m.get_row_info(rows[j]))
                    
                    m.remove_row(rows[j])
                #endfor
                            
                f = file(acs_c1p, 'w')      # C1P matrix file
                            
                f.write(str(m))
                            
                f.close()

                f = open(acs_discarded, 'w')

                mat_rem.write(f.write)

                f.close()
  
            elif make_C1P == "circ_heuristic":
                m = bm.BinaryMatrix()

                m.from_file(self.acs_file)

                mat2, mat_rem = cc1p.make_circC1P(m, 'max')

                f = open(self.acs_c1p, 'w')

                mat2.write(f.write)

                f.close()

                f = open(self.acs_discarded, 'w')

                mat_rem.write(f.write)

                f.close()

            elif make_C1P == "branch_and_bound":
                prop = False     # True if using weights as probs
                mat = bm.BinaryMatrix()     # matrix
                            
                mat.from_file(self.acs_file)
                            
                ms = c1p.split_matrix(mat)      # split matrices
                matb = bm.BinaryMatrix()        # C1P matrix
                mat_rem = bm.BinaryMatrix()     # rows removed
                            
                j = 1       # iterator for tracing
                            
                del mat

                for m in ms:
                    print str(j) + '/' + str(len(ms))
                                
                    j += 1
                            
                    # sort matrix to get heuristic as first answer
                    m.sort()
                            
                    # branch and bound
                    # optimal rows to remove to make compnonent C1P
                    rows = bab.branch_and_bound(m, prop, BABTester(m._height))
                        
                    rows.sort()
                        
                    for i in xrange(len(rows) - 1, -1, -1):
                        row = m.get_row_info(rows[i])       # row to remove
                                    
                        mat_rem.add_row_info(row)
                            
                        m.remove_row(rows[i])
                    #endfor
                    
                    # collect usable rows into the C1P matrix
                    for r in m._rows:
                        matb.add_row_info(r)
                    #endfor
                    
                    print ''
                #endfor
                    
                f = file(acs_c1p, 'w')
                            
                f.write(str(matb))
                            
                f.close()

                f = file(acs_discarded, 'w')
                            
                f.write(str(mat_rem))
                            
                f.close()
            elif make_C1P == "circ_branch_and_bound":
                mat = bm.BinaryMatrix()     # matrix
            
                mat.from_file(self.acs_file)

                mode = 'max'

                if mode == 'whole':
                    matb, mat_rem = circC1P_bab(mat)
                elif mode == 'max':
                    maxs = c1p.make_intersect_components(mat)       # split matrices
                    matb = bm.BinaryMatrix()        # C1P matrix
                    mat_rem = bm.BinaryMatrix()     # rows removed
                    
                    j = 1       # iterator for tracing
                                
                    del mat
                    
                    for max_comp in maxs:
                        print 'Max:' + str(j) + '/' + str(len(maxs)) + ' '
                        
                        j += 1
                        
                        circC1P_bab(max_comp, matb, mat_rem)
                    #endfor
                #endif
                    
                f = file(self.acs_c1p, 'w')
                            
                f.write(str(matb))
                            
                f.close()

                f = file(self.acs_discarded, 'w')
                            
                f.write(str(mat_rem))
                            
                f.close()

            elif make_C1P == "m_heuristic":
                m = bm.BinaryMatrix()       # matrix
                mat_rem = bm.BinaryMatrix() # rows discarded
                                        
                m.from_file(self.acs_file)
                            
                rows = mc1p.make_mC1P(m)        # rows to remove to make C1P
                            
                for j in xrange(len(rows) - 1, -1, -1):
                    mat_rem.add_row_info(m.get_row_info(rows[j]))

                    m.remove_row(rows[j])
                #endfor
                            
                f = file(self.acs_c1p, 'w')      # C1P matrix file
                            
                f.write(str(m))
                            
                f.close()

                f = open(self.acs_discarded, 'w')

                mat_rem.write(f.write)

                f.close()

            elif make_C1P == "m_branch_and_bound":
                prop = False
                mat = bm.BinaryMatrix()     # matrix
                rows_rem = []       # rows removed

                mat.from_file(self.acs_file)

                ms = c1p.make_intersect_components(mat)     # intersect components
                matb = bm.BinaryMatrix()        # mC1P matrix
                mat_rem = bm.BinaryMatrix()     # rows removed
                            
                j = 1       # iterator for tracing

                for m in ms:
                    print str(j) + '/' + str(len(ms)) + ' ',
                    sys.stdout.flush()
                                
                    j += 1
                    
                    # sort matrix to get heuristic as first answer
                    m.sort()
                    
                    # branch and bound
                    # optimal rows to remove to make compnonent mC1P
                    rows = bab.branch_and_bound(m, prop, MC1PTester(m))
                        
                    rows.sort()

                    for i in xrange(len(rows) - 1, -1, -1):
                        row = m.get_row_info(rows[i])       # row to remove
                                    
                        mat_rem.add_row_info(row)
                            
                        m.remove_row(rows[i])
                    #endfor
                    
                    # collect usable rows into the C1P matrix
                    for r in m._rows:
                        matb.add_row_info(r)
                    #endfor
                    
                    print ''
                #endfor

                f = file(acs_c1p, 'w')

                f.write(str(matb))

                f.close()

                f = file(acs_discarded, 'w')
                            
                f.write(str(mat_rem))
                            
                f.close()

            elif make_C1P == "m_branch_and_bound_both":
                prop = False     # True if using weights as probs
                mat = bm.BinaryMatrix()     # matrix

                mat.from_file(self.acs_file)

                matb, mat_rem = C1P_and_mC1P_bab(mat)

                f = file(acs_c1p, 'w')

                f.write(str(matb))

                f.close()

                f = file(acs_discarded, 'w')
                            
                f.write(str(mat_rem))
                            
                f.close()

        else:
            if make_C1P == "heuristic":
                m = bm.BinaryMatrix()       # matrix
                mat_rem = bm.BinaryMatrix() # rows discarded
                                        
                m.from_file(self.acs_file)

                rows = c1p.make_C1P(m)      # rows to remove to make C1P
                            
                for j in xrange(len(rows) - 1, -1, -1):
                    mat_rem.add_row_info(m.get_row_info(rows[j]))
                    
                    m.remove_row(rows[j])
                #endfor
                            
                f = file(acs_c1p, 'w')      # C1P matrix file
                            
                f.write(str(m))
                            
                f.close()

                f = open(acs_discarded, 'w')

                mat_rem.write(f.write)

                f.close()
  
            elif make_C1P == "circ_heuristic":
                m = bm.BinaryMatrix()

                m.from_file(self.acs_file)

                mat2, mat_rem = cc1p.make_circC1P(m, 'max')

                f = open(self.acs_c1p, 'w')

                mat2.write(f.write)

                f.close()

                f = open(self.acs_discarded, 'w')

                mat_rem.write(f.write)

                f.close()

            elif make_C1P == "branch_and_bound":
                prop = False     # True if using weights as probs
                mat = bm.BinaryMatrix()     # matrix
                            
                mat.from_file(self.acs_file)
                            
                ms = c1p.split_matrix(mat)      # split matrices
                matb = bm.BinaryMatrix()        # C1P matrix
                mat_rem = bm.BinaryMatrix()     # rows removed
                            
                j = 1       # iterator for tracing
                            
                del mat

                for m in ms:
                    print str(j) + '/' + str(len(ms))
                                
                    j += 1
                            
                    # sort matrix to get heuristic as first answer
                    m.sort()
                            
                    # branch and bound
                    # optimal rows to remove to make compnonent C1P
                    rows = bab.branch_and_bound(m, prop, BABTester(m._height))
                        
                    rows.sort()
                        
                    for i in xrange(len(rows) - 1, -1, -1):
                        row = m.get_row_info(rows[i])       # row to remove
                                    
                        mat_rem.add_row_info(row)
                            
                        m.remove_row(rows[i])
                    #endfor
                    
                    # collect usable rows into the C1P matrix
                    for r in m._rows:
                        matb.add_row_info(r)
                    #endfor
                    
                    print ''
                #endfor
                    
                f = file(acs_c1p, 'w')
                            
                f.write(str(matb))
                            
                f.close()

                f = file(acs_discarded, 'w')
                            
                f.write(str(mat_rem))
                            
                f.close()
            elif make_C1P == "circ_branch_and_bound":
                mat = bm.BinaryMatrix()     # matrix
            
                mat.from_file(self.acs_file)

                mode = 'whole'

                if mode == 'whole':
                    matb, mat_rem = circC1P_bab(mat)
                elif mode == 'max':
                    maxs = c1p.make_intersect_components(mat)       # split matrices
                    matb = bm.BinaryMatrix()        # C1P matrix
                    mat_rem = bm.BinaryMatrix()     # rows removed
                    
                    j = 1       # iterator for tracing
                                
                    del mat
                    
                    for max_comp in maxs:
                        print 'Max:' + str(j) + '/' + str(len(maxs)) + ' '
                        
                        j += 1
                        
                        circC1P_bab(max_comp, matb, mat_rem)
                    #endfor
                #endif
                    
                f = file(self.acs_c1p, 'w')
                            
                f.write(str(matb))
                            
                f.close()

                f = file(self.acs_discarded, 'w')
                            
                f.write(str(mat_rem))
                            
                f.close()

            elif make_C1P == "m_heuristic":
                m = bm.BinaryMatrix()       # matrix
                mat_rem = bm.BinaryMatrix() # rows discarded
                                        
                m.from_file(self.acs_file)
                            
                rows = mc1p.make_mC1P(m)        # rows to remove to make C1P
                            
                for j in xrange(len(rows) - 1, -1, -1):
                    mat_rem.add_row_info(m.get_row_info(rows[j]))

                    m.remove_row(rows[j])
                #endfor
                            
                f = file(self.acs_c1p, 'w')      # C1P matrix file
                            
                f.write(str(m))
                            
                f.close()

                f = open(self.acs_discarded, 'w')

                mat_rem.write(f.write)

                f.close()

            elif make_C1P == "m_branch_and_bound":
                prop = False
                mat = bm.BinaryMatrix()     # matrix
                rows_rem = []       # rows removed

                mat.from_file(self.acs_file)

                ms = c1p.make_intersect_components(mat)     # intersect components
                matb = bm.BinaryMatrix()        # mC1P matrix
                mat_rem = bm.BinaryMatrix()     # rows removed
                            
                j = 1       # iterator for tracing

                for m in ms:
                    print str(j) + '/' + str(len(ms)) + ' ',
                    sys.stdout.flush()
                                
                    j += 1
                    
                    # sort matrix to get heuristic as first answer
                    m.sort()
                    
                    # branch and bound
                    # optimal rows to remove to make compnonent mC1P
                    rows = bab.branch_and_bound(m, prop, MC1PTester(m))
                        
                    rows.sort()

                    for i in xrange(len(rows) - 1, -1, -1):
                        row = m.get_row_info(rows[i])       # row to remove
                                    
                        mat_rem.add_row_info(row)
                            
                        m.remove_row(rows[i])
                    #endfor
                    
                    # collect usable rows into the C1P matrix
                    for r in m._rows:
                        matb.add_row_info(r)
                    #endfor
                    
                    print ''
                #endfor

                f = file(acs_c1p, 'w')

                f.write(str(matb))

                f.close()

                f = file(acs_discarded, 'w')
                            
                f.write(str(mat_rem))
                            
                f.close()

            elif make_C1P == "m_branch_and_bound_both":
                prop = False     # True if using weights as probs
                mat = bm.BinaryMatrix()     # matrix

                mat.from_file(self.acs_file)

                matb, mat_rem = C1P_and_mC1P_bab(mat)

                f = file(acs_c1p, 'w')

                f.write(str(matb))

                f.close()

                f = file(acs_discarded, 'w')
                            
                f.write(str(mat_rem))
                            
                f.close()

        #endif

        if not self.quiet:
            print("----> Creating a PQ-tree: " + pq_tree)
        #endif

        if self.markers_doubled:
            self.computePQRtree()
            if not self.quiet:
                print("----> Halving PQ-tree columns") 
            #endif
            self.halvePQRtree(pq_tree_doubled, pq_tree)
        else:
            self.computePQRtree()
        #endif
    #enddef

    def computePQRtree(self):
        if self.markers_doubled:
            m = bm.BinaryMatrix()       # matrix
            m.from_file(self.acs_file)
                        
            f = file(self.pqr_tree_doubled, 'w')      # PQR-tree file

            f.write(">" + self.output_ancestor + "\n")
            pqtree.make_PQR_tree(m).write(f.write)
                        
            f.close()
            if not self.quiet:
                print("----> Halving PQR-tree columns") 
                    #endif
            self.halvePQRtree(self.pqr_tree_doubled, self.pqr_tree)
        else:
            m = bm.BinaryMatrix()       # matrix
            m.from_file(self.acs_file)
                        
            f = file(self.pqr_tree, 'w')      # PQR-tree file

            f.write(">" + self.output_ancestor + "\n")
            pqtree.make_PQR_tree(m).write(f.write)
                        
            f.close()
        #endif

    def halvePQRtree(self, pqr_tree_doubled, pqr_tree):
        name_doubled=pqr_tree_doubled
        name_halved=pqr_tree

        input_file=open(name_doubled,"r").readlines()
        output_file=open(name_halved,"w")

        for i in range(len(input_file)):
            if input_file[i][0]==">" or input_file[i][0]=="#":
                output_file.write(input_file[i])
            else:
                mots=input_file[i].split()
                m=0
                while m<len(mots):
                    if mots[m].find("_")>=0 or mots[m] == "T":
                        output_file.write(mots[m]+" ")
                        m=m+1
                    else:
                        m1=int(mots[m])
                        try:
                            m2=int(mots[m+1])
                        except:
                            if mots[m-1].find("C")>=0:
                                m2=-100
                            else:
                                m=m+1
                                
                                continue
                            #endif
                        #endtry
                        
                        if abs(m1-m2)!=1:
                            if mots[m-1].find("C")>=0:
                                if m1%2==0:
                                    output_file.write("-"+str(m1/2)+" ")
                                else:
                                    output_file.write(str((m1+1)/2)+" ")
                                #endif
                                
                                m=m+1
                                
                                continue
                            else:
                                print("ERROR:  C1P_halve_PQRtree " + str(m1) + " " + str(m2))
                                sys.exit(0)
                            #endif
                        #endif
                        
                        if m1<m2 and m2%2==0:
                            output_file.write(str(m2/2)+" ")
                        elif m2<m1 and m1%2==0:
                            output_file.write("-"+str(m1/2)+" ")
                        else:
                            if mots[m-1].find("C")>=0:
                                if m1%2==0:
                                    output_file.write("-"+str(m1/2)+" ")
                                else:
                                    output_file.write(str((m1+1)/2)+" ")
                                #endif
                                
                                m=m+1
                                
                                continue
                            else:
                                print("ERROR:  C1P_halve_PQRtree " + str(m1) + " " + str(m2))
                                sys.exit(0)
                            #endif
                        #endif
                        
                        m=m+2
                    #endif
                #endwhile
                
                output_file.write("\n")
            #endif
        #endwhile

        output_file.close()



    def computePQCRtree(self):
        if self.markers_doubled:
            m = bm.BinaryMatrix()       # matrix
            m.from_file(self.acs_file)
                        
            f = file(self.pqr_tree_doubled, 'w')      # PQCR-tree file

            f.write(">" + self.output_ancestor + "\n")
            pqtree.make_PQCR_tree(m).write(f.write)
                        
            f.close()
            if not self.quiet:
                print("----> Halving PQR-tree columns") 
                    #endif
            self.halvePQRtree(self.pqr_tree_doubled, self.pq_tree)

        else:
            m = bm.BinaryMatrix()       # matrix
            m.from_file(self.acs_file)
                        
            f = file(self.pqr_tree, 'w')      # PQCR-tree file

            f.write(">" + self.output_ancestor + "\n")
            file_write = pqtree.make_PQCR_tree(m)
            file_write.write(f.write)
                        
            f.close()
        #endif

    def computePQtreeDotProductCorrelationMatrix(self, pq_tree):
        m = bm.BinaryMatrix()       # matrix
        m.from_file(self.acs_file)

        alpha=float(str(self.c1p_spectral_alpha))

        # Create numpy matrix

        indices={}
        i=0

        # indices is a dictionary using the column names as keys for the matrix indices.
        for col in m.get_support():
            indices[col]=i
            i+=1 

        # Check if there are columns that only appear as Xs.
        for row in m:
            try:
                row._Xs
            except:
                row._Xs=[]
            for col in row._Xs:
                if col not in indices.keys():
                    indices[col]=i
                    i+=1    

        A=numpy.zeros([len(indices),len(indices)])

        # Compute correlation matrix
        for col1 in m.get_support():
            for row in m:
                try:
                    row._Xs
                except:
                    row._Xs= []
                    if col1 in row._set:
                        for col2 in row._set:
                            A[indices[col1],indices[col2]]+=1
                        if row._Xs:
                            for col2 in row._Xs:
                                A[indices[col1],indices[col2]]+=alpha
                    elif row._Xs and col1 in row._Xs:
                        for col2 in row._set:
                            A[indices[col1],indices[col2]]+=alpha
                        for col2 in row._Xs:
                            A[indices[col1],indices[col2]]+=alpha**2
                    
        # Find PQ tree using seriation
        T=spectral.seriation(A,indices,0)

        o=open(pq_tree,'w')
        o.write('>'+self.output_ancestor.upper()+'\n')
        i=1

        for t in T:
            if '_P' in t.printTree().split(' ') or '_Q' in t.printTree().split(' '):
                o.write('#CAR'+str(i)+'\n')
            o.write(t.printTree()+'\n')
            i+=1

        o.close()

    def computePQtreeCorrelationMatrix(self, pq_tree):
        m = bm.BinaryMatrix()       # matrix
        m.from_file(self.acs_file)

        # Create numpy matrix

        indices={}
        i=0

        # indices is a dictionary using the column names as keys for the matrix indices.
        for col in m.get_support():
            indices[col]=i
            i+=1 

        M=numpy.zeros([m._height,len(indices)])

        # Populate binary matrix
        i=0
        for row in m:
            for col in row._set:
                M[i,indices[col]]=1
            i+=1

        del m

        # Compute correlation matrix as matrix product
        A=numpy.dot(M.T,M)
        # Find PQ tree using seriation

        T=spectral.seriation(A,indices,0)

        o=open(pq_tree,'w')
        o.write('>'+self.output_ancestor.upper()+'\n')
        i=1

        for t in T:
            o.write('#CAR'+str(i)+'\n')
            o.write(t.printTree()+'\n')
            i+=1

        o.close()

       

    def makeC1PSpectral(self):
        pq_tree = self.cars_dir + "/" + self.output_prefix + "PQTREE_SERIATION"    

        if self.acs_correction >= 2:
            if not self.quiet:
                print("----> Computing PQ-tree from correlation matrix of a ternary matrix using parameter alpha") 
                    #endif
            self.computePQtreeDotProductCorrelationMatrix(pq_tree)
        else:
            if not self.quiet:
                print("----> Computing PQ-tree using spectral seriation on correlation matrix") 
                    #endif
            self.computePQtreeCorrelationMatrix(pq_tree)

            #endif



    #weight
    def makeC1PHeuristic(self):
        self.do_c1p("HEUR", "heuristic", "/C1P/C1P_make_C1P_heuristic.py", "/C1P/C1P_compute_PQRtree.py")


    def makeCircC1PHeuristic(self):
        self.do_c1p("HEUR", "heuristic", "/C1P/C1P_make_circC1P_heuristic.py", "/C1P/C1P_compute_PQCRtree.py", True)

    def makeC1PBranchAndBound(self):
        self.do_c1p("BAB", "branch and bound", "/C1P/C1P_make_C1P_branch_and_bound.py", "/C1P/C1P_compute_PQRtree.py")

    def makeCircC1PBranchAndBound(self):
        self.do_c1p("BAB", "branch and bound", "/C1P/C1P_make_circC1P_branch_and_bound.py", "/C1P/C1P_compute_PQCRtree.py", True)



    def makeMC1PHeuristic(self):
        self.do_c1p("TEL_HEUR", "heuristic", "/C1P/C1P_make_mC1P_heuristic.py", "/C1P/C1P_compute_PQRtree.py")

    def makeMC1PBranchAndBoundBoth(self):
        self.do_c1p("TEL_BAB1", "branch and bound, added after", "/C1P/C1P_make_mC1P_branch_and_bound_both.py", "/C1P/C1P_compute_PQRtree.py")

    def makeMC1PBranchAndBound(self):
        self.do_c1p("TEL_BAB2", "branch and bound, added during", "/C1P/C1P_make_mC1P_branch_and_bound.py", "/C1P/C1P_compute_PQRtree.py")

    def run(self):
        self.code_dir = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + "/code/c1p_files" # directory where the code is stored
        working_dir = os.path.dirname(self.species_tree)

        # set the working directory to the directory of the parameters file
        if working_dir and working_dir != '':
            os.chdir(working_dir)

        ## CREATING DIRECTORIES ------------------------------------------------------------------------------
        # create output directory
        try:
            os.mkdir(output_dir)
        except:
            pass
        #endtry
    
        print("----> Started")

        ## CHECKING INCOMPATIBILITIES ------------------------------------------------------------------------
        ## General errors
        self.checkErrors()

        ## COMPUTING ACS ------------------------------------------------------------------------------

        ## COMPUTING A C1P MATRIX ------------------------------------------------------------------------------

        if self.markers_doubled:
            self.pqr_tree_doubled = self.cars_dir + "/" + self.output_prefix + "PQRTREE_DOUBLED"
        # make C1P directory
        try:
            os.mkdir(self.c1p_dir)
        except:
             pass
        # make CARS directory
        try:
            os.mkdir(self.cars_dir)
        except:
            pass

        if self.c1p_linear or self.c1p_telomeres > 0:
            if not self.quiet:
                print("----> Creating PQR-tree: " + self.pqr_tree)
            self.computePQRtree()

        elif self.c1p_circular:
            if not self.quiet:
                print("----> Creating PQR-tree: " + self.pqr_tree)
            self.computePQCRtree()

        # heuristic
        if self.c1p_heuristic and self.c1p_telomeres == 0:
            if self.c1p_linear:
                self.makeC1PHeuristic()
            elif self.c1p_circular:
                self.makeCircC1PHeuristic()
            
        # branch and bound
        if self.c1p_bab and self.c1p_telomeres == 0:
            if self.c1p_linear:
                self.makeC1PBranchAndBound()
            elif self.c1p_circular:
                self.makeCircC1PBranchAndBound()
            
        # seriation
        if self.c1p_spectral:
            self.makeC1PSpectral()
        
        # telomere heuristic
        if self.c1p_telomeres == 1:
            self.makeMC1PHeuristic()
        
        # telomere branch and bound, add after
        if self.c1p_telomeres == 2:
            self.makeMC1PBranchAndBound()
        
        #telomere branch and bound, add during
        if self.c1p_telomeres == 3:
            self.makeMC1PBranchAndBoundBoth()
        
        if not self.quiet:
            print("----> Done")
        




    
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

    def isCopyNumGreaterThanOne(self, hom_fam_list):
        filtered_list = filter(lambda fam: fam.copy_number <= 1, hom_fam_list)
        filtered_quant = len(hom_fam_list) - len(filtered_list)

        if filtered_quant != 0:
            return True
        else:
            return False
 
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
        self.config_file_directory = ""
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

        self.received_acs = False
        self.doC1P = False
        self.doMWM = False
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
        config = {}
        try:
            execfile(config_file, config)
        except IOError:
            print("{}  ERROR (master.py -> process.py) - could not open configuration file: {}\n"
                    .format(strtime(),config_file))
            sys.exit()
        
        #collect the information from config file
        self.io_dict["homologous_families"]           = config["homologous_families"]
        self.io_dict["species_tree"]                  = config["species_tree"]
        self.io_dict["output_directory"]              = config["output_directory"]
        self.io_dict["acs_file"]                     = config["acs_file"]

        self.markers_param_dict["markers_doubled"]    = config["markers_doubled"]
        self.markers_param_dict["markers_unique"]     = config["markers_unique"]
        self.markers_param_dict["markers_universal"]  = config["markers_universal"]
        self.markers_param_dict["markers_overlap"]    = config["markers_overlap"]
        self.markers_param_dict["filter_copy_number"] = config["filter_copy_number"]
        self.markers_param_dict["filter_by_id"]       = config["filter_by_id"]

        self.run_param_dict["all_match"]            = config["all_match"]

        self.debug = config["debug"]

        self.config_file_directory = config_file
        if config["acs_file"] != "":
            self.received_acs = True
        
        config.clear()
        self.setOutputStreams() #set Log and Debug

    #enddef

    def receivedAcsFile(self):
        return self.received_acs

    def parse_markersPhase(self):
        """
        Deal with everything realated to input (configuration, markers and species pairs) in order to take information from these files

        """
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

        if markers_phase_obj.isCopyNumGreaterThanOne(self.hom_fam_list):
            self.doC1P = False
            self.doMWM = True
        else:
            self.doC1P = True
            self.doMWM = False

        del markers_phase_obj
    #enddef

    def adjacenciesPhase(self):
        # Genome construction
        self.genome_construction_obj.constructGenomes(self.species_pairs, self.hom_fam_list, self.log)
        self.genome_construction_obj.dealWithAdjPhase(self.species_pairs, self.hom_fam_list, self.io_dict["output_directory"], self.log, self.debug, self.run_param_dict["all_match"])
    #enddef








    def c1pPhase(self):
        c1p_obj = MasterC1P()
        c1p_obj.setConfigParams(self.config_file_directory)
        c1p_obj.run()







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

    def doC1PorNot(self):
        return self.doC1P # if false, do MWM

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
