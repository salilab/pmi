#!/usr/bin/env python
import IMP
import IMP.algebra




############################
#####Analysis tools
############################



class Alignment():

    """
    This class performs alignment and RMSD calculation for two sets of coordinates
 
    Inputs:

      - query = {'p1':coords(L,3), 'p2':coords(L,3)}
      - template = {'p1':coords(L,3), 'p2':coords(L,3)}

    The class also takes into accout non-equal stoichiometry of the proteins. If this
    is the case, the protein names of protein in multiple copies should be delivered 
    in the following form: nameA..1, nameA..2 (note two dots).
    """

    def __init__(self, template, query):

        global array,argwhere,mgrid,shape,reshape,zeros,diagonal,argsort,deepcopy,cdist,sqrt
        global product, permutations
        from numpy import array,argwhere,mgrid,shape,reshape,zeros,diagonal,argsort
        from copy import deepcopy
        from scipy.spatial.distance import cdist
        from math import sqrt
        from itertools import permutations, product

        self.query = query
        self.template = template
        
        if len(self.query.keys()) != len(self.template.keys()): print '''ERROR: the number of proteins
                               in template and query does not match!''';exit()

    def permute(self):

        self.proteins = sorted(self.query.keys())
        prots_uniq = [i.split('..')[0] for i in self.proteins]
        P = {}
        for p in prots_uniq:
            np = prots_uniq.count(p)
            copies = [i for i in self.proteins if i.split('..')[0]==p]
            prmts = list(permutations( copies, len(copies) )) 
            P[p] = prmts
        self.P = P
        self.Product = list(product(*P.values()))

    def get_rmsd(self):

        self.permute()

        template_xyz = []
        torder = sum([list(i) for i in self.Product[0]],[])
        for t in torder:
           template_xyz += [i for i in self.template[t]]
        template_xyz = array(template_xyz)
        
        self.rmsd = 10000000000.
        for comb in self.Product:
            order = sum([list(i) for i in comb],[])
            query_xyz = []
            for p in order: query_xyz += [i for i in self.query[p]]
            query_xyz = array(query_xyz)
            if len(template_xyz) != len(query_xyz):
                print '''Alignment.get_rmsd: ERROR: the number of coordinates
                               in template and query does not match!''';exit()
            dist = sqrt(sum(diagonal(cdist(template_xyz,query_xyz)**2)) / len(template_xyz))
            if dist < self.rmsd: self.rmsd = dist
        return self.rmsd

    def align(self):

        self.permute()

        template_xyz = []
        torder = sum([list(i) for i in self.Product[0]],[])
        for t in torder:
           template_xyz += [IMP.algebra.Vector3D(i) for i in self.template[t]]
        #template_xyz = array(template_xyz)
        
        self.rmsd, Transformation = 10000000000.,''
        for comb in self.Product:
            order = sum([list(i) for i in comb],[])
            query_xyz = []
            for p in order: query_xyz += [IMP.algebra.Vector3D(i) for i in self.query[p]]
            #query_xyz = array(query_xyz)

            if len(template_xyz) != len(query_xyz):
                print '''ERROR: the number of coordinates
                               in template and query does not match!''';exit()

            transformation = IMP.algebra.get_transformation_aligning_first_to_second(query_xyz,template_xyz)
            query_xyz_tr = [transformation.get_transformed(n) for n in query_xyz]
            
            dist = sqrt(sum(diagonal(cdist(template_xyz,query_xyz_tr)**2)) / len(template_xyz))
            if dist < self.rmsd:
                self.rmsd = dist
                Transformation = transformation
            
        return (self.rmsd, Transformation)
            

### TEST for the alignment ###
"""
import numpy as np
Proteins = {'a..1':np.array([np.array([-1.,1.])]),
            'a..2':np.array([np.array([1.,1.,])]),
            'a..3':np.array([np.array([-2.,1.])]),
            'b':np.array([np.array([0.,-1.])]),
            'c..1':np.array([np.array([-1.,-1.])]),
            'c..2':np.array([np.array([1.,-1.])]),
            'd':np.array([np.array([0.,0.])]),
            'e':np.array([np.array([0.,1.])])}

Ali = Alignment(Proteins, Proteins)
Ali.permute()
if Ali.get_rmsd() == 0.0: print 'successful test!'
else: print 'ERROR!'; exit()
"""      
        
         

# ----------------------------------
class Violations():

    def __init__(self, filename):
        global impem,deepcopy,cdist,array,argwhere,mgrid,shape,reshape,zeros,sqrt,diagonal,argsort
        import IMP.em as impem
        from numpy import array,argwhere,mgrid,shape,reshape,zeros,diagonal,argsort
        from copy import deepcopy
        from scipy.spatial.distance import cdist
        from math import sqrt
        self.violation_thresholds = {}
        self.violation_counts = {}
   
        data = open(filename)
        D = data.readlines()
        data.close()

        for d in D:
            d = d.strip().split()
            self.violation_thresholds[d[0]] = float(d[1])

    def get_number_violated_restraints(self, rsts_dict):
        num_violated = 0
        for rst in self.violation_thresholds:
            if rst not in rsts_dict: continue #print rst; 
            if float(rsts_dict[rst]) > self.violation_thresholds[rst]:
                num_violated += 1
                if rst not in self.violation_counts: self.violation_counts[rst] = 1
                else: self.violation_counts[rst] += 1
        return num_violated



# ----------------------------------
class Clustering():

    def __init__(self):

        global impem,deepcopy,cdist,array,argwhere,mgrid,shape,reshape,zeros,sqrt,diagonal,argsort,npsum
        import IMP.em as impem
        from numpy import array,argwhere,mgrid,shape,reshape,zeros,diagonal,argsort
        from numpy import sum as npsum
        from copy import deepcopy
        from scipy.spatial.distance import cdist
        from math import sqrt
        self.all_coords = {}
        self.structure_cluster_ids=None

    def set_template(self, part_coords):

        self.tmpl_coords = part_coords

    def fill(self, frame, Coords):
        """
        fill stores coordinates of a model into a dictionary all_coords,
        containint coordinates for all models.
        """

        self.all_coords[frame]= Coords

       
    def dist_matrix(self,is_mpi=False):
        from itertools import combinations
        
        if is_mpi:
          from mpi4py import MPI
          comm = MPI.COMM_WORLD
          rank = comm.Get_rank()
          number_of_processes = comm.size
        else:
          number_of_processes = 1
          rank=0
        
        self.model_list_names= self.all_coords.keys()
        self.model_indexes=range(len(self.model_list_names))
        self.model_indexes_dict=dict(zip(self.model_list_names,self.model_indexes))
        model_indexes_unique_pairs=list(combinations(self.model_indexes,2))
        
        my_model_indexes_unique_pairs=IMP.pmi.tools.chunk_list_into_segments(model_indexes_unique_pairs,number_of_processes)[rank]

        print "process %s assigned with %s pairs" % (str(rank),str(len(my_model_indexes_unique_pairs)))
        
        (raw_distance_dict,self.transformation_distance_dict)=self.matrix_calculation(self.all_coords,
                                                           self.tmpl_coords,
                                                           my_model_indexes_unique_pairs,
                                                           do_alignment=True)
        
        if number_of_processes>1:
            raw_distance_dict=IMP.pmi.tools.scatter_and_gather(raw_distance_dict)     
            pickable_transformations=self.get_pickable_transformation_distance_dict()
            pickable_transformations=IMP.pmi.tools.scatter_and_gather(pickable_transformations) 
            self.set_transformation_distance_dict_from_pickable(pickable_transformations)   

        self.raw_distance_matrix = zeros((len(self.model_list_names), len(self.model_list_names)))
        for item in raw_distance_dict:
                (f1,f2)=item
                self.raw_distance_matrix[f1,f2]= raw_distance_dict[item]
                self.raw_distance_matrix[f2,f1]= raw_distance_dict[item]           
    
    def do_cluster(self,number_of_clusters):
        from sklearn.cluster import KMeans
        kmeans=KMeans(n_clusters=number_of_clusters)
        kmeans.fit_predict(self.raw_distance_matrix)
        
        self.structure_cluster_ids=kmeans.labels_

    def get_pickable_transformation_distance_dict(self):
        pickable_transformations={}
        for label in self.transformation_distance_dict:
            tr=self.transformation_distance_dict[label]
            trans=tuple(tr.get_translation())
            rot=tuple(tr.get_rotation().get_quaternion())
            pickable_transformations[label]=(rot,trans)
        return pickable_transformations        
    
    def set_transformation_distance_dict_from_pickable(self,pickable_transformations):
        self.transformation_distance_dict={}
        for label in pickable_transformations:
            tr=pickable_transformations[label]
            trans=IMP.algebra.Vector3D(tr[1])
            rot=IMP.algebra.Rotation3D(tr[0])
            self.transformation_distance_dict[label]=IMP.algebra.Transformation3D(rot,trans)


    def save_distance_matrix_file(self,file_name='cluster.rawmatrix.pkl'):
        import pickle    
        outf = open(file_name,'w')
        
        #to pickle the transformation dictionary 
        #you have to save the arrays correposnding to
        # the transformations
        
        
        pickable_transformations=self.get_pickable_transformation_distance_dict()
        pickle.dump((self.structure_cluster_ids,self.model_list_names,self.raw_distance_matrix,pickable_transformations),outf)            
        outf.close()        

    def load_distance_matrix_file(self,file_name='cluster.rawmatrix.pkl'):
        import pickle        
        inputf = open(file_name,'r')    
        (self.structure_cluster_ids,self.model_list_names,self.raw_distance_matrix,pickable_transformation)=pickle.load(inputf) 
        
        self.set_transformation_distance_dict_from_pickable(pickable_transformation)
        self.model_indexes=range(len(self.model_list_names))
        self.model_indexes_dict=dict(zip(self.model_list_names,self.model_indexes))
        
        inputf.close()      
   
    def plot_matrix(self):
        import pylab as pl    
        from scipy.cluster import hierarchy as hrc    
        
        fig = pl.figure()
        ax = fig.add_subplot(211)
        dendrogram = hrc.dendrogram(hrc.linkage(self.raw_distance_matrix),color_threshold=7)
        leaves_order = dendrogram['leaves']

        ax = fig.add_subplot(212)
        cax = ax.imshow(self.raw_distance_matrix[leaves_order,:][:,leaves_order], interpolation='nearest')
        ax.set_yticks(range(len(self.model_list_names)))
        ax.set_yticklabels( [self.model_list_names[i] for i in leaves_order] )
        fig.colorbar(cax)

        pl.show()
    
    def get_model_index_from_name(self,name):
        return self.model_indexes_dict[name]
    
    def get_cluster_labels(self):
        #this list 
        return list(set(self.structure_cluster_ids))
    
    def get_number_of_clusters(self):
        return len(self.get_cluster_labels())
    
    def get_cluster_label_indexes(self, label):
        return [i for i,l in enumerate(self.structure_cluster_ids) if l==label]
        
    def get_cluster_label_names(self, label):
        return [self.model_list_names[i] for i in self.get_cluster_label_indexes(label)]
    
    def get_cluster_label_average_rmsd(self,label):
        
        indexes=self.get_cluster_label_indexes(label)
        
        if len(indexes)>1:
           sub_distance_matrix= self.raw_distance_matrix[indexes,:][:,indexes]
           average_rmsd=npsum(sub_distance_matrix)/(len(sub_distance_matrix)**2 - len(sub_distance_matrix))
        else:
           average_rmsd=0.0
        return average_rmsd
    
    def get_cluster_label_size(self,label):
        return len(self.get_cluster_label_indexes(label))
        
    def get_transformation_to_first_member(self,cluster_label,structure_index):
        reference=self.get_cluster_label_indexes(cluster_label)[0]
        return self.transformation_distance_dict[(reference,structure_index)]
        
    
    def matrix_calculation(self,all_coords, template_coords, list_of_pairs, do_alignment=False):

            import IMP
            import IMP.pmi
            import IMP.pmi.analysis
            
            model_list_names = all_coords.keys()
            rmsd_protein_names = all_coords[model_list_names[0]].keys()
            alignment_template_protein_names=template_coords.keys()
            raw_distance_dict={}
            transformation_distance_dict={}
            for (f1,f2) in list_of_pairs:

                    if not do_alignment:
                        # here we only get the rmsd, 
                        # we need that for instance when you want to cluster conformations
                        # globally, eg the EM map is a reference
                        Ali = IMP.pmi.analysis.Alignment(all_coords[model_list_names[f1]], all_coords[model_list_names[f2]])
                        r= Ali.get_rmsd()
                        
                    elif do_alignment:
                        # here we actually align the conformations first
                        # and than calculate the rmsd. We need that when the 
                        # protein(s) is the reference
                        coords_f1 = dict([(pr,all_coords[model_list_names[f1]][pr]) for pr in alignment_template_protein_names])
                        coords_f2 = dict([(pr,all_coords[model_list_names[f2]][pr]) for pr in alignment_template_protein_names])
                        
                        Ali = IMP.pmi.analysis.Alignment(coords_f1, coords_f2)
                        template_rmsd, transformation = Ali.align()
                        
                        # here we calculate the rmsd
                        # we will align two models based n the nuber of subunits provided
                        # and transform coordinates of model 2 to model 1
                        coords_f1 = dict([(pr,all_coords[model_list_names[f1]][pr]) for pr in rmsd_protein_names])
                        coords_f2 = {}
                        for pr in rmsd_protein_names:
                            coords_f2[pr] = [transformation.get_transformed(i) for i in all_coords[model_list_names[f2]][pr]]
                        
                        Ali = IMP.pmi.analysis.Alignment(coords_f1, coords_f2)
                        rmsd= Ali.get_rmsd()
                       
                    raw_distance_dict[(f1,f2)]= rmsd
                    raw_distance_dict[(f2,f1)]= rmsd
                    transformation_distance_dict[(f1,f2)]= transformation
                    transformation_distance_dict[(f2,f1)]= transformation
            
            return raw_distance_dict,transformation_distance_dict         
      

# ----------------------------------

class GetModelDensity():

    def __init__(self, custom_ranges,representation=None, voxel=5.0):
        '''
        custom_ranges ={'kin28':[['kin28',1,-1]],
                'density_name_1' :[('ccl1')],
                'density_name_2' :[(1,142,'tfb3d1'),(143,700,'tfb3d2')],
                
       
        '''
        global impem
        import IMP.em as impem
        
        self.representation=representation           
        self.voxel=voxel
        self.densities= {}
        self.count_models=0.0
        self.custom_ranges=custom_ranges

    def add_subunits_density(self,hierarchy=None):
          # the hierarchy is optional, if passed
          self.count_models+=1.0
          for density_name in self.custom_ranges:
              parts=[]
              if hierarchy:
                 all_particles_by_segments=[]
                 
              for seg in self.custom_ranges[density_name]:
                  if not hierarchy:
                     parts+=IMP.tools.select_by_tuple(self.representation,seg,resolution=1,name_is_ambiguous=False)
                  else:
                     
                     if type(seg)==str:
                        children=[child for child in hierarchy.get_children() if child.get_name()==seg]
                        s=IMP.atom.Selection(children)
                     if type(seg)==tuple:   
                        children=[child for child in hierarchy.get_children() if child.get_name()==seg[2]]
                        s=IMP.atom.Selection(children,residue_indexes=range(seg[0],seg[1]+1))
                     all_particles_by_segments+=s.get_selected_particles()   
              
              if hierarchy:       
                 part_dict=get_particles_at_resolution_one(hierarchy)
                 all_particles_by_resolution=[]
                 for name in part_dict:
                     all_particles_by_resolution+=part_dict[name]
                 parts=list(set(all_particles_by_segments) & set(all_particles_by_resolution))
                 
                 
              self.create_density_from_particles(parts, density_name)
          
    def normalize_density(self):
        pass

    def create_density_from_particles(self,ps, name,
                                  resolution=1,
                                  kernel_type='GAUSSIAN'):
        '''pass XYZR particles with mass and create a density from them.
        kernel type options are GAUSSIAN, BINARIZED_SPHERE, and SPHERE.'''
      
        kd={'GAUSSIAN':IMP.em.GAUSSIAN, 'BINARIZED_SPHERE':IMP.em.BINARIZED_SPHERE,
        'SPHERE':IMP.em.SPHERE}
   
        dmap=impem.SampledDensityMap(ps,resolution,self.voxel)
        dmap.calcRMS()
        if name not in self.densities: self.densities[name] = dmap
        else: self.densities[name].add(dmap)

    def write_mrc(self, prefix):

        for density_name in self.densities:
            self.densities[density_name].multiply(1./self.count_models)
            impem.write_map(self.densities[density_name],prefix+"_"+density_name+".mrc",impem.MRCReaderWriter())

# ----------------------------------

class GetContactMap():
    def __init__(self, distance=15.):
        global impem,deepcopy,cdist,array,argwhere,mgrid,shape,reshape,zeros,sqrt,diagonal,argsort,log
        import IMP.em as impem
        from numpy import array,argwhere,mgrid,shape,reshape,zeros,diagonal,argsort,log
        from copy import deepcopy
        from scipy.spatial.distance import cdist
        global itemgetter
        from operator import itemgetter
        
        self.distance = distance
        self.contactmap = ''
        self.namelist = []
        self.xlinks = 0
        self.XL = {}
        self.expanded = {}
        self.resmap = {}

    def set_prot(self, prot):
        self.prot = prot

    def get_subunit_coords(self,frame, align=0):
        coords= []
        radii= []
        namelist = []
        test,testr = [],[]
        for part in self.prot.get_children():
            SortedSegments = []
            print part
            for chl in part.get_children():
                start = IMP.atom.get_leaves(chl)[0]
                end   = IMP.atom.get_leaves(chl)[-1]

                startres = IMP.atom.Fragment(start).get_residue_indexes()[0]
                endres   = IMP.atom.Fragment(end).get_residue_indexes()[-1]
                SortedSegments.append((chl,startres))
            SortedSegments = sorted(SortedSegments, key=itemgetter(1))

            for sgmnt in SortedSegments:
                for leaf in IMP.atom.get_leaves(sgmnt[0]):
                    p= IMP.core.XYZR(leaf)
                    crd = array([p.get_x(),p.get_y(),p.get_z()])

                    coords.append(crd)
                    radii.append(p.get_radius())
                   
                    new_name = part.get_name()+'_'+sgmnt[0].get_name()+\
                                    '_'+str(IMP.atom.Fragment(leaf).get_residue_indexes()[0])
                    namelist.append(new_name)
                    self.expanded[new_name] = len(IMP.atom.Fragment(leaf).get_residue_indexes())
                    if part.get_name() not in self.resmap: self.resmap[part.get_name()] = {}
                    for res in IMP.atom.Fragment(leaf).get_residue_indexes():
                        self.resmap[part.get_name()][res] = new_name

        coords = array(coords)
        radii = array(radii)
        if len(self.namelist)==0:
            self.namelist = namelist
            self.contactmap = zeros((len(coords), len(coords)))
        distances = cdist(coords, coords)
        distances = (distances-radii).T - radii
        distances = distances<=self.distance
        self.contactmap += distances


    def add_xlinks(self, filname):
        self.xlinks = 1
        data = open(filname)
        D = data.readlines()
        data.close()

        for d in D:
            d = d.strip().split()
            t1, t2 = (d[0],d[1]), (d[1],d[0])
            if t1 not in self.XL:
                self.XL[t1] = [(int(d[2])+1, int(d[3])+1)]
                self.XL[t2] = [(int(d[3])+1, int(d[2])+1)]
            else:
                self.XL[t1].append((int(d[2])+1, int(d[3])+1))
                self.XL[t2].append((int(d[3])+1, int(d[2])+1))

        


    def dist_matrix(self, skip_cmap=0, skip_xl=1):
        K= self.namelist
        M= self.contactmap
        C,R = [],[]
        L= sum(self.expanded.values())

        # exp new
        if skip_cmap==0:
            Matrices = {}
            proteins = [p.get_name() for p in self.prot.get_children()]
            missing = []
            for p1 in xrange(len(proteins)):
                for p2 in xrange(p1,len(proteins)):
                    pl1,pl2=max(self.resmap[proteins[p1]].keys()),max(self.resmap[proteins[p2]].keys())
                    pn1,pn2=proteins[p1],proteins[p2]
                    mtr=zeros((pl1+1,pl2+1))
                    print 'Creating matrix for: ',p1,p2,pn1,pn2,mtr.shape,pl1,pl2
                    for i1 in xrange(1,pl1+1):
                        for i2 in xrange(1,pl2+1):
                            try:
                                r1=K.index( self.resmap[pn1][i1] )
                                r2=K.index( self.resmap[pn2][i2] )
                                r=M[r1,r2]
                                mtr[i1-1,i2-1]=r
                            except KeyError: missing.append((pn1,pn2,i1,i2)); pass
                    Matrices[(pn1,pn2)]=mtr

        # add cross-links
        if skip_xl==0:
            if self.XL=={}: print "ERROR: cross-links were not provided, use add_xlinks function!"; exit()
            Matrices_xl = {}
            proteins = [p.get_name() for p in self.prot.get_children()]
            missing_xl = []
            for p1 in xrange(len(proteins)):
                for p2 in xrange(p1,len(proteins)):
                    pl1,pl2=max(self.resmap[proteins[p1]].keys()),max(self.resmap[proteins[p2]].keys())
                    pn1,pn2=proteins[p1],proteins[p2]
                    mtr=zeros((pl1+1,pl2+1))
                    flg=0
                    try: xls = self.XL[(pn1,pn2)]
                    except KeyError:
                        try: xls = self.XL[(pn2,pn1)]; flg=1
                        except KeyError: flg=2
                    if flg==0:
                        print 'Creating matrix for: ',p1,p2,pn1,pn2,mtr.shape,pl1,pl2
                        for xl1,xl2 in xls:
                            if xl1>pl1: print 'X'*10,xl1,xl2; xl1=pl1
                            if xl2>pl2: print 'X'*10,xl1,xl2; xl2=pl2
                            mtr[xl1-1,xl2-1]=100
                    elif flg==1:
                        print 'Creating matrix for: ',p1,p2,pn1,pn2,mtr.shape,pl1,pl2
                        for xl1,xl2 in xls:
                            if xl1>pl1: print 'X'*10,xl1,xl2; xl1=pl1
                            if xl2>pl2: print 'X'*10,xl1,xl2; xl2=pl2
                            mtr[xl2-1,xl1-1]=100
                    else: print 'WTF!'; exit()
                    Matrices_xl[(pn1,pn2)]=mtr                

        # expand the matrix to individual residues
        NewM = []
        for x1 in xrange(len(K)):
            lst = []
            for x2 in xrange(len(K)):
                lst += [M[x1,x2]]*self.expanded[K[x2]]
            for i in xrange(self.expanded[K[x1]]): NewM.append(array(lst))
        NewM = array(NewM)

        # make list of protein names and create coordinate lists  
        C = proteins
        W,R = [],[]
        for i,c in enumerate(C):
            cl = max(self.resmap[c].keys())
            W.append(cl)
            if i==0: R.append(cl)
            else: R.append(R[-1]+cl)
        
        # start plotting
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        import scipy.sparse as sparse

        f = plt.figure()
        gs = gridspec.GridSpec(len(W), len(W),
                       width_ratios=W,
                       height_ratios=W)

        cnt = 0
        for x1,r1 in enumerate(R):
            if x1==0: s1=0
            else: s1 = R[x1-1]
            for x2,r2 in enumerate(R):
                if x2==0: s2=0
                else: s2 = R[x2-1]

                ax = plt.subplot(gs[cnt])
                if skip_cmap==0:
                    try: mtr = Matrices[(C[x1],C[x2])]
                    except KeyError: mtr = Matrices[(C[x2],C[x1])].T
                    #cax = ax.imshow(log(NewM[s1:r1,s2:r2] / 1.), interpolation='nearest', vmin=0., vmax=log(NewM.max()))
                    cax = ax.imshow(log(mtr), interpolation='nearest', vmin=0., vmax=log(NewM.max()))
                    ax.set_xticks([])
                    ax.set_yticks([])
                if skip_xl==0:
                    try: mtr = Matrices_xl[(C[x1],C[x2])]
                    except KeyError: mtr = Matrices_xl[(C[x2],C[x1])].T
                    cax = ax.spy(sparse.csr_matrix(mtr), markersize=10, color='white', linewidth=100, alpha=0.5)
                    ax.set_xticks([])
                    ax.set_yticks([])
                                        
                cnt+=1
                if x2==0: ax.set_ylabel(C[x1], rotation=90)
        plt.show()
        
        
###############################################
# these are post production function analysis
###############################################



def get_hier_from_rmf(model,frame_number,rmf_file):
      import IMP.rmf
      import RMF    
      print "getting coordinates for frame %i rmf file %s" % (frame_number,rmf_file)
                  
      # load the frame
      rh= RMF.open_rmf_file_read_only(rmf_file)
      try:
        prots=IMP.rmf.create_hierarchies(rh, model)  
      except:
        print "Unable to open rmf file %s" % (rmf_file)
        prot=None
        return prot                   
      IMP.rmf.link_hierarchies(rh, prots)
      prot=prots[0]
      try: 
        IMP.rmf.load_frame(rh, frame_number)        
      except:
        print "Unable to open frame %i of file %s" % (frame_number,rmf_file)
        prot=None           
      model.update()
      return prot


def get_particles_at_resolution_one(prot):
    '''
    this fucntion get the particles by resolution, without a Representation class initialized
    it is mainly used when the hierarchy is read from an rmf file
    it returns a dictionary of component names and their particles
    '''
    particle_dict={}
    allparticles=[]
    for c in prot.get_children():
        name=c.get_name()
        particle_dict[name]=IMP.atom.get_leaves(c)
        for s in c.get_children():
            if "_Res:1" in s.get_name() and "_Res:10" not in s.get_name(): 
                allparticles+=IMP.atom.get_leaves(s)
            if "Beads" in s.get_name():
                allparticles+=IMP.atom.get_leaves(s)
    
    
    particle_align=[]
    for name in particle_dict:
        
        particle_dict[name]=IMP.pmi.tools.sort_by_residues(list(set(particle_dict[name]) & set(allparticles)))
    
    return particle_dict

def select_by_tuple(first_res_last_res_name_tuple):
    first_res=first_res_last_res_hier_tuple[0]
    last_res=first_res_last_res_hier_tuple[1]
    name=first_res_last_res_hier_tuple[2]

