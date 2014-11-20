#!/usr/bin/env python

"""@namespace IMP.pmi.analysis
   Post-clustering analysis
"""

import IMP
import IMP.algebra
import IMP.em
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.rmf
import RMF
from operator import itemgetter
from copy import deepcopy
from math import log,sqrt
import itertools
import numpy as np

class Precision(object):
    def __init__(self,model,resolution='one',selection_dictionary={}):

        ''' selection_dictionary is a dictionary where we store coordinates
            selection_dictionary = {"Selection_name_1":selection_tuple1,
                                    "Selection_name_2":selection_tuple2}'''

        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            self.rank = comm.Get_rank()
            self.number_of_processes = comm.size
        except ImportError:
            self.number_of_processes=1
            self.rank=0


        self.styles=['pairwise_rmsd','pairwise_drmsd_k','pairwise_drmsd_Q','pairwise_drms_k','pairwise_rmsd','drmsd_from_center']
        self.style='pairwise_drmsd_k'
        self.structures_dictionary={}
        self.reference_structures_dictionary={}
        self.prots=[]
        self.protein_names=None
        self.len_particles_resolution_one=None
        self.model=model
        self.rmf_names_frames={}
        self.reference_rmf_names_frames=None
        self.reference_structure=None
        self.reference_prot=None
        self.selection_dictionary=selection_dictionary
        self.threshold=40.0

        self.residue_particle_index_map=None

        if resolution in ["one", "ten"]:
            self.resolution=resolution
        else:
            raise KeyError("no such resolution")



    def get_structure(self,rmf_frame_index,rmf_name):

        rh= RMF.open_rmf_file_read_only(rmf_name)
        prots=IMP.rmf.create_hierarchies(rh, self.model)
        IMP.rmf.load_frame(rh, rmf_frame_index)
        print "getting coordinates for frame %i rmf file %s" % (rmf_frame_index, rmf_name)
        del rh

        if self.resolution=='one':
            particle_dict=get_particles_at_resolution_one(prots[0])
        elif self.resolution=='ten':
            particle_dict=get_particles_at_resolution_ten(prots[0])

        protein_names=particle_dict.keys()
        particles_resolution_one=[]
        for k in particle_dict:
            particles_resolution_one+=(particle_dict[k])

        if self.protein_names==None:
            self.protein_names=protein_names
        else:
            if self.protein_names!=protein_names:
                print "Error: the protein names of the new coordinate set is not compatible with the previous one"

        if self.len_particles_resolution_one==None:
            self.len_particles_resolution_one=len(particles_resolution_one)
        else:
            if self.len_particles_resolution_one!=len(particles_resolution_one):
                raise ValueError("the new coordinate set is not compatible with the previous one")

        return particles_resolution_one, prots


    def add_structures(self,rmf_name_frame_tuples,structure_set_name):
        '''this function helps reading a list of rmfs, the input is a list of
        tuples, containing the file name and the frame number'''

        my_rmf_name_frame_tuples=IMP.pmi.tools.chunk_list_into_segments(
            rmf_name_frame_tuples,self.number_of_processes)[self.rank]


        for nfr,tup in enumerate(my_rmf_name_frame_tuples):
            rmf_name=tup[0]
            rmf_frame_index=tup[1]
            # the first frame stores the map between residues and particles
            self.add_structure(rmf_name,rmf_frame_index,structure_set_name)
            try:
                (a,b)=self.get_structure(rmf_frame_index,rmf_name)
            except:
                print "something wrong with the rmf"
                continue
            if self.residue_particle_index_map is None:


                try:
                    (particles_resolution_one, prots)=self.get_structure(rmf_frame_index,rmf_name)
                except:
                    print "something wrong with the rmf"
                    continue

                self.residue_particle_index_map={}

                for prot_name in self.protein_names:

                    self.residue_particle_index_map[prot_name]=self.get_residue_particle_index_map(prot_name,particles_resolution_one,prots[0])


        if self.number_of_processes > 1:
            # synchronize the internal self.selection_dictionary data structure
            self.rmf_names_frames=IMP.pmi.tools.scatter_and_gather(self.rmf_names_frames)

            if self.rank != 0:
                comm.send(self.structures_dictionary, dest=0, tag=11)

            elif self.rank == 0:
                for i in range(1, self.number_of_processes):
                    data_tmp = comm.recv(source=i, tag=11)

                    for key in self.structures_dictionary:

                        self.structures_dictionary[key].update(data_tmp[key])

                for i in range(1, self.number_of_processes):
                    comm.send(self.structures_dictionary, dest=i, tag=11)

            if self.rank != 0:
                self.structures_dictionary = comm.recv(source=0, tag=11)




    def add_structure(self,rmf_name,rmf_frame_index,structure_set_name):
        '''
        @structure_set_name string specify a name for the structure set. This is used for the precision and
               cross-precision
        '''
        if structure_set_name in self.structures_dictionary:
            cdict=self.structures_dictionary[structure_set_name]
            rmflist=self.rmf_names_frames[structure_set_name]
        else:
            self.structures_dictionary[structure_set_name]={}
            self.rmf_names_frames[structure_set_name]=[]
            cdict=self.structures_dictionary[structure_set_name]
            rmflist=self.rmf_names_frames[structure_set_name]

        try:
            (particles_resolution_one, prots)=self.get_structure(rmf_frame_index,rmf_name)
        except:
            print "something wrong with the rmf"
            return 0

        self.selection_dictionary.update({"All":self.protein_names})

        for selection_name in self.selection_dictionary:
            selection_tuple=self.selection_dictionary[selection_name]
            coords=self.select_coordinates(selection_tuple,particles_resolution_one,prots[0])

            if selection_name not in cdict:
                cdict[selection_name]=[coords]
            else:
                cdict[selection_name].append(coords)

        rmflist.append((rmf_name,rmf_frame_index))

        for prot in prots: IMP.atom.destroy(prot)



    def get_residue_particle_index_map(self,prot_name,structure,hier):
        residue_particle_index_map=[]
        s=IMP.atom.Selection(hier,molecules=[prot_name])
        all_selected_particles=s.get_selected_particles()
        intersection=list(set(all_selected_particles) & set(structure))
        sorted_intersection=IMP.pmi.tools.sort_by_residues(intersection)
        for p in sorted_intersection:
            residue_particle_index_map.append(IMP.pmi.tools.get_residue_indexes(p))
        return residue_particle_index_map


    def select_coordinates(self,tuple_selections,structure,prot):
        selected_coordinates=[]
        for t in tuple_selections:
            if type(t)==tuple and len(t)==3:
                s=IMP.atom.Selection(prot,molecules=[t[2]],residue_indexes=range(t[0],t[1]+1))
                all_selected_particles=s.get_selected_particles()
                intersection=list(set(all_selected_particles) & set(structure))
                sorted_intersection=IMP.pmi.tools.sort_by_residues(intersection)
                cc=map(lambda p: tuple(IMP.core.XYZ(p).get_coordinates()), sorted_intersection)
                selected_coordinates+=cc

            elif type(t)==str:
                s=IMP.atom.Selection(prot,molecules=[t])
                all_selected_particles=s.get_selected_particles()
                intersection=list(set(all_selected_particles) & set(structure))
                sorted_intersection=IMP.pmi.tools.sort_by_residues(intersection)
                cc=map(lambda p: tuple(IMP.core.XYZ(p).get_coordinates()), sorted_intersection)
                selected_coordinates+=cc
            else:
                raise ValueError("Selection error")
        return selected_coordinates

    def set_threshold(self,threshold):
        self.threshold=threshold


    def get_distance(self,structure_set_name1,structure_set_name2,selection_name,index1,index2):
        c1=self.structures_dictionary[structure_set_name1][selection_name][index1]
        c2=self.structures_dictionary[structure_set_name2][selection_name][index2]

        coordinates1=map(lambda c: IMP.algebra.Vector3D(c), c1)
        coordinates2=map(lambda c: IMP.algebra.Vector3D(c), c2)

        if self.style=='pairwise_drmsd_k':
            distance=IMP.atom.get_drmsd(coordinates1,coordinates2)
        if self.style=='pairwise_drms_k':
            distance=IMP.atom.get_drms(coordinates1,coordinates2)
        if self.style=='pairwise_drmsd_Q':
            distance=IMP.atom.get_drmsd_Q(coordinates1,coordinates2,self.threshold)

        if self.style=='pairwise_rmsd':
            distance=IMP.atom.get_rmsd(coordinates1,coordinates2,IMP.algebra.get_identity_transformation_3d())

        return distance

    def get_particle_distances(self,structure_set_name1,structure_set_name2,
                               selection_name,index1,index2):
        import numpy as np
        c1=self.structures_dictionary[structure_set_name1][selection_name][index1]
        c2=self.structures_dictionary[structure_set_name2][selection_name][index2]

        coordinates1=map(lambda c: IMP.algebra.Vector3D(c), c1)
        coordinates2=map(lambda c: IMP.algebra.Vector3D(c), c2)

        distances=[np.linalg.norm(a-b) for (a,b) in zip(coordinates1,coordinates2)]

        return distances

    def get_precision(self,outfile,structure_set_name1,structure_set_name2,
                      skip=None,selection_keywords=None):

        '''
        when the structure_set_name1 is different from the structure_set_name2, you get the
        cross-precision
        structure_set_name1  string name of the first structure set
        structure_set_name2  string name of the second structure set
        skip Int analyse every skip structure for the distance matrix calculation
        outfile   str  the output file name
        skip      int  skip every n frames

        '''
        of=open(outfile,"w")

        if selection_keywords is None:
            sel_keys=self.selection_dictionary.keys()
        else:
            sel_keys=selection_keywords
        centroid_index=0
        for selection_name in sel_keys:

            number_of_structures_1=len(self.structures_dictionary[structure_set_name1][selection_name])
            number_of_structures_2=len(self.structures_dictionary[structure_set_name2][selection_name])

            distances={}

            structure_pointers_1=range(number_of_structures_1)
            structure_pointers_2=range(number_of_structures_2)

            if skip is not None:
                structure_pointers_1=structure_pointers_1[0::skip]
                structure_pointers_2=structure_pointers_2[0::skip]

            pair_combination_list=list(itertools.product(structure_pointers_1,structure_pointers_2))

            if len(pair_combination_list)==0:
                raise ValueError("no structure selected. Check the skip parameter.")

            my_pair_combination_list=IMP.pmi.tools.chunk_list_into_segments(
                pair_combination_list,self.number_of_processes)[self.rank]
            my_length=len(my_pair_combination_list)

            for n,pair in enumerate(my_pair_combination_list):
                progression=int(float(n)/my_length*100.0)
                #print self.rank,progression,len(pair_combination_list),my_length,n
                distances[pair]=self.get_distance(structure_set_name1,structure_set_name2,
                                                  selection_name,pair[0],pair[1])

            if self.number_of_processes > 1:
                distances = IMP.pmi.tools.scatter_and_gather(distances)

            if self.rank == 0:

                if structure_set_name1==structure_set_name2:
                    structure_pointers=structure_pointers_1
                    number_of_structures=number_of_structures_1

                    # calculate the distance from the first centroid
                    # and determine the centroid

                    distance=0.0
                    distances_to_structure={}
                    distances_to_structure_normalization={}

                    for n in structure_pointers:
                        distances_to_structure[n]=0.0
                        distances_to_structure_normalization[n]=0

                    for k in distances:
                        distance+=distances[k]
                        distances_to_structure[k[0]]+=distances[k]
                        distances_to_structure[k[1]]+=distances[k]
                        distances_to_structure_normalization[k[0]]+=1
                        distances_to_structure_normalization[k[1]]+=1

                    for n in structure_pointers:
                        distances_to_structure[n]=distances_to_structure[n]/distances_to_structure_normalization[n]

                    #for n,d in enumerate(distances_to_structure):
                    #    print self.rmf_names_frames[n],d
                    min_distance=min([distances_to_structure[n] for n in distances_to_structure])
                    centroid_index=[k for k, v in distances_to_structure.iteritems() if v == min_distance][0]
                    centroid_rmf_name=self.rmf_names_frames[structure_set_name1][centroid_index]

                    centroid_distance=0.0
                    for n in range(number_of_structures):
                        centroid_distance+=self.get_distance(structure_set_name1,structure_set_name1,
                                                             selection_name,centroid_index,n)


                    #pairwise_distance=distance/len(distances.keys())
                    centroid_distance/=number_of_structures
                    #average_centroid_distance=sum(distances_to_structure)/len(distances_to_structure)
                    of.write(str(selection_name)+" "+structure_set_name1+
                                    " average centroid distance "+str(centroid_distance)+"\n")
                    of.write(str(selection_name)+" "+structure_set_name1+
                                    " centroid index "+str(centroid_index)+"\n")
                    of.write(str(selection_name)+" "+structure_set_name1+
                                    " centroid rmf name "+str(centroid_rmf_name)+"\n")

                average_pairwise_distances=sum(distances.values())/len(distances.values())
                of.write(str(selection_name)+" "+structure_set_name1+" "+structure_set_name2+
                           " average pairwise distance "+str(average_pairwise_distances)+"\n")


        of.close()
        return centroid_index

    def get_rmsf(self,structure_set_name,outdir="./",
                 skip=None,make_plot=False,set_plot_yaxis_range=None):

        outfile=outdir+"/rmsf.dat"
        # get the centroid structure for the whole complex
        centroid_index=self.get_precision(outfile,structure_set_name,structure_set_name,
                                          skip=skip,selection_keywords=["All"])

        for p in self.protein_names:
            outfile=outdir+"/rmsf."+p+".dat"
            of=open(outfile,"w")
            selection_name=p
            self.selection_dictionary.update({selection_name:[p]})
            try:
                number_of_structures=len(self.structures_dictionary[structure_set_name][selection_name])
            except KeyError:
                # that protein was not included in the selection
                continue
            rpim=self.residue_particle_index_map[p]

            residue_distances={}
            residue_nblock={}
            for index in range(number_of_structures):


                distances=self.get_particle_distances(structure_set_name,
                                                      structure_set_name,
                                                      selection_name,
                                                      centroid_index,index)
                for nblock,block in enumerate(rpim):

                    for residue_number in block:
                        residue_nblock[residue_number]=nblock
                        if residue_number not in residue_distances:
                            residue_distances[residue_number]=[distances[nblock]]
                        else:
                            residue_distances[residue_number].append(distances[nblock])

            residues=[]
            rmsfs=[]
            for rn in residue_distances:
                residues.append(rn)
                rmsf=np.std(residue_distances[rn])
                rmsfs.append(rmsf)
                of.write(str(rn)+" "+str(residue_nblock[rn])+" "+str(rmsf)+"\n")

            IMP.pmi.output.plot_xy_data(residues,rmsfs,title=outdir+"/rmsf."+p,display=False,
                                       set_plot_yaxis_range=set_plot_yaxis_range)
            of.close()




    def set_reference_structure(self,rmf_name,rmf_frame_index):
        (particles_resolution_one, prot)=self.get_structure(rmf_frame_index,rmf_name)
        self.reference_rmf_names_frames=(rmf_name,rmf_frame_index)


        for selection_name in self.selection_dictionary:
            selection_tuple=self.selection_dictionary[selection_name]
            coords=self.select_coordinates(selection_tuple,
                                       particles_resolution_one,prot)
            self.reference_structures_dictionary[selection_name]=coords


    def get_average_distance_wrt_reference_structure(self,structure_set_name,tuple_selections=None):

        for selection_name in self.selection_dictionary:
            reference_coordinates=self.reference_structures_dictionary[selection_name]
            coordinates2=map(lambda c: IMP.algebra.Vector3D(c), reference_coordinates)
            distances=[]

            for sc in self.structures_dictionary[structure_set_name][selection_name]:

                coordinates1=map(lambda c: IMP.algebra.Vector3D(c), sc)


                if self.style=='pairwise_drmsd_k':
                    distance=IMP.atom.get_drmsd(coordinates1,coordinates2)
                if self.style=='pairwise_drms_k':
                    distance=IMP.atom.get_drms(coordinates1,coordinates2)
                if self.style=='pairwise_drmsd_Q':
                    distance=IMP.atom.get_drmsd_Q(coordinates1,coordinates2,self.threshold)
                if self.style=='pairwise_rmsd':
                    distance=IMP.atom.get_rmsd(coordinates1,coordinates2,IMP.algebra.get_identity_transformation_3d())

                distances.append(distance)


            print selection_name,"average distance",sum(distances)/len(distances),"minimum distance",min(distances)


    def get_coordinates(self):
        pass

    def set_precision_style(self, style):
        if style in self.styles:
            self.style=style
        else:
            raise ValueError("No such style")
