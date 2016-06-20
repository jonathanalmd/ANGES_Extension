
import os
import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))) + '/COMMON')

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

        self.m = bm.BinaryMatrix()

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
        acs_c1p = self.c1p_dir + "/" + self.output_prefix + "ACS_C1P_" + suffix    # C1P matrix file
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
                #m = bm.BinaryMatrix()       # matrix
                mat_rem = bm.BinaryMatrix() # rows discarded
                                        
                #m.from_file(self.acs_file)

                rows = c1p.make_C1P(self.m)      # rows to remove to make C1P
                            
                for j in xrange(len(rows) - 1, -1, -1):
                    mat_rem.add_row_info(m.get_row_info(rows[j]))
                    
                    self.m.remove_row(rows[j])
                #endfor
                            
                f = file(acs_c1p, 'w')      # C1P matrix file
                            
                f.write(str(m))
                            
                f.close()

                f = open(acs_discarded, 'w')

                mat_rem.write(f.write)

                f.close()
  
            elif make_C1P == "circ_heuristic":
                #m = bm.BinaryMatrix()

                #m.from_file(self.acs_file)

                mat2, mat_rem = cc1p.make_circC1P(self.m, 'max')

                f = open(self.acs_c1p, 'w')

                mat2.write(f.write)

                f.close()

                f = open(self.acs_discarded, 'w')

                mat_rem.write(f.write)

                f.close()

            elif make_C1P == "branch_and_bound":
                prop = False     # True if using weights as probs
                #mat = bm.BinaryMatrix()     # matrix
                            
                #m.from_file(self.acs_file)
                            
                ms = c1p.split_matrix(self.m)      # split matrices
                matb = bm.BinaryMatrix()        # C1P matrix
                mat_rem = bm.BinaryMatrix()     # rows removed
                            
                j = 1       # iterator for tracing
                            
                #del mat

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
                        row = self.m.get_row_info(rows[i])       # row to remove
                                    
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
                #mat = bm.BinaryMatrix()     # matrix
            
                #mat.from_file(self.acs_file)

                mode = 'max'

                if mode == 'whole':
                    matb, mat_rem = circC1P_bab(self.m)
                elif mode == 'max':
                    maxs = c1p.make_intersect_components(self.m)       # split matrices
                    matb = bm.BinaryMatrix()        # C1P matrix
                    mat_rem = bm.BinaryMatrix()     # rows removed
                    
                    j = 1       # iterator for tracing
                                
                    #del mat
                    
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
                #m = bm.BinaryMatrix()       # matrix
                mat_rem = bm.BinaryMatrix() # rows discarded
                                        
                #m.from_file(self.acs_file)
                            
                rows = mc1p.make_mC1P(m)        # rows to remove to make C1P
                            
                for j in xrange(len(rows) - 1, -1, -1):
                    mat_rem.add_row_info(self.m.get_row_info(rows[j]))

                    self.m.remove_row(rows[j])
                #endfor
                            
                f = file(self.acs_c1p, 'w')      # C1P matrix file
                            
                f.write(str(self.m))
                            
                f.close()

                f = open(self.acs_discarded, 'w')

                mat_rem.write(f.write)

                f.close()

            elif make_C1P == "m_branch_and_bound":
                prop = False
                #mat = bm.BinaryMatrix()     # matrix
                rows_rem = []       # rows removed

                #mat.from_file(self.acs_file)

                ms = c1p.make_intersect_components(self.m)     # intersect components
                matb = bm.BinaryMatrix()        # mC1P matrix
                mat_rem = bm.BinaryMatrix()     # rows removed
                            
                j = 1       # iterator for tracing

                for m in ms:
                    print str(j) + '/' + str(len(ms)) + ' ',
                    sys.stdout.flush()
                                
                    j += 1
                    
                    # sort matrix to get heuristic as first answer
                    self.m.sort()
                    
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
                    for r in self.m._rows:
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
                #mat = bm.BinaryMatrix()     # matrix

                #mat.from_file(self.acs_file)

                matb, mat_rem = C1P_and_mC1P_bab(self.m)

                f = file(acs_c1p, 'w')

                f.write(str(matb))

                f.close()

                f = file(acs_discarded, 'w')
                            
                f.write(str(mat_rem))
                            
                f.close()

        else:
            if make_C1P == "heuristic":
                #m = bm.BinaryMatrix()       # matrix
                mat_rem = bm.BinaryMatrix() # rows discarded
                                        
                #m.from_file(self.acs_file)

                rows = c1p.make_C1P(self.m)      # rows to remove to make C1P
                            
                for j in xrange(len(rows) - 1, -1, -1):
                    mat_rem.add_row_info(m.get_row_info(rows[j]))
                    
                    self.m.remove_row(rows[j])
                #endfor
                            
                f = file(acs_c1p, 'w')      # C1P matrix file
                            
                f.write(str(m))
                            
                f.close()

                f = open(acs_discarded, 'w')

                mat_rem.write(f.write)

                f.close()
  
            elif make_C1P == "circ_heuristic":
                #m = bm.BinaryMatrix()

                #m.from_file(self.acs_file)

                mat2, mat_rem = cc1p.make_circC1P(self.m, 'max')

                f = open(self.acs_c1p, 'w')

                mat2.write(f.write)

                f.close()

                f = open(self.acs_discarded, 'w')

                mat_rem.write(f.write)

                f.close()

            elif make_C1P == "branch_and_bound":
                prop = False     # True if using weights as probs
                # mat = bm.BinaryMatrix()     # matrix
                            
                # mat.from_file(self.acs_file)
                            
                ms = c1p.split_matrix(self.m)      # split matrices
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
                # mat = bm.BinaryMatrix()     # matrix
            
                # mat.from_file(self.acs_file)

                mode = 'whole'

                if mode == 'whole':
                    matb, mat_rem = circC1P_bab(self.m)
                elif mode == 'max':
                    maxs = c1p.make_intersect_components(self.m)       # split matrices
                    matb = bm.BinaryMatrix()        # C1P matrix
                    mat_rem = bm.BinaryMatrix()     # rows removed
                    
                    j = 1       # iterator for tracing
                                
                    # del mat
                    
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
                #m = bm.BinaryMatrix()       # matrix
                mat_rem = bm.BinaryMatrix() # rows discarded
                                        
                #m.from_file(self.acs_file)
                            
                rows = mc1p.make_mC1P(self.m)        # rows to remove to make C1P
                            
                for j in xrange(len(rows) - 1, -1, -1):
                    mat_rem.add_row_info(m.get_row_info(rows[j]))

                    self.m.remove_row(rows[j])
                #endfor
                            
                f = file(self.acs_c1p, 'w')      # C1P matrix file
                            
                f.write(str(self.m))
                            
                f.close()

                f = open(self.acs_discarded, 'w')

                mat_rem.write(f.write)

                f.close()

            elif make_C1P == "m_branch_and_bound":
                prop = False
                # mat = bm.BinaryMatrix()     # matrix
                rows_rem = []       # rows removed

                # mat.from_file(self.acs_file)

                ms = c1p.make_intersect_components(self.m)     # intersect components
                matb = bm.BinaryMatrix()        # mC1P matrix
                mat_rem = bm.BinaryMatrix()     # rows removed
                            
                j = 1       # iterator for tracing

                for m in ms:
                    print str(j) + '/' + str(len(ms)) + ' ',
                    sys.stdout.flush()
                                
                    j += 1
                    
                    # sort matrix to get heuristic as first answer
                    self.m.sort()
                    
                    # branch and bound
                    # optimal rows to remove to make compnonent mC1P
                    rows = bab.branch_and_bound(m, prop, MC1PTester(m))
                        
                    rows.sort()

                    for i in xrange(len(rows) - 1, -1, -1):
                        row = self.m.get_row_info(rows[i])       # row to remove
                                    
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


                matb, mat_rem = C1P_and_mC1P_bab(self.m)

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
            #m = bm.BinaryMatrix()       # matrix
            #m.from_file(self.acs_file)
                        
            f = file(self.pqr_tree_doubled, 'w')      # PQR-tree file

            f.write(">" + self.output_ancestor + "\n")
            pqtree.make_PQR_tree(self.m).write(f.write)
                        
            f.close()
            if not self.quiet:
                print("----> Halving PQR-tree columns") 
                    #endif
            self.halvePQRtree(self.pqr_tree_doubled, self.pqr_tree)
        else:
            #m = bm.BinaryMatrix()       # matrix
            #m.from_file(self.acs_file)
                        
            f = file(self.pqr_tree, 'w')      # PQR-tree file

            f.write(">" + self.output_ancestor + "\n")
            pqtree.make_PQR_tree(self.m).write(f.write)
                        
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
            #m = bm.BinaryMatrix()       # matrix
            #m.from_file(self.acs_file)
                        
            f = file(self.pqr_tree_doubled, 'w')      # PQCR-tree file

            f.write(">" + self.output_ancestor + "\n")
            pqtree.make_PQCR_tree(self.m).write(f.write)
                        
            f.close()
            if not self.quiet:
                print("----> Halving PQR-tree columns") 
                    #endif
            self.halvePQRtree(self.pqr_tree_doubled, self.pq_tree)

        else:
            #m = bm.BinaryMatrix()       # matrix
            #m.from_file(self.acs_file)
                        
            f = file(self.pqr_tree, 'w')      # PQCR-tree file

            f.write(">" + self.output_ancestor + "\n")
            file_write = pqtree.make_PQCR_tree(self.m)
            file_write.write(f.write)
                        
            f.close()
        #endif

    def computePQtreeDotProductCorrelationMatrix(self, pq_tree):
        #m = bm.BinaryMatrix()       # matrix
        #m.from_file(self.acs_file)

        alpha=float(str(self.c1p_spectral_alpha))

        # Create numpy matrix

        indices={}
        i=0

        # indices is a dictionary using the column names as keys for the matrix indices.
        for col in self.m.get_support():
            indices[col]=i
            i+=1 

        # Check if there are columns that only appear as Xs.
        for row in self.m:
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
        for col1 in self.m.get_support():
            for row in self.m:
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
        #m = bm.BinaryMatrix()       # matrix
        #m.from_file(self.acs_file)

        # Create numpy matrix

        indices={}
        i=0

        # indices is a dictionary using the column names as keys for the matrix indices.
        for col in self.m.get_support():
            indices[col]=i
            i+=1 

        M=numpy.zeros([self.m._height,len(indices)])

        # Populate binary matrix
        i=0
        for row in self.m:
            for col in row._set:
                M[i,indices[col]]=1
            i+=1

        # del m

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

    def bmFromACSFile(self):
        self.m.from_file(self.acs_file)

    def run(self):
        self.bmFromACSFile()
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