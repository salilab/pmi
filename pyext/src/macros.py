import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import os



class ReplicaExchange0():

    def __init__(self,model,
                      representation,
                      sample_objects,
                      output_objects,
                      crosslink_restraints=None,
                      monte_carlo_temperature=1.0,
                      replica_exchange_minimum_temperature=1.0,
                      replica_exchange_maximum_temperature=2.5,
                      number_of_best_scoring_models=500,
                      monte_carlo_steps=10,
                      number_of_frames=1000,
                      write_initial_rmf=True,
                      initial_rmf_name_suffix="initial",
                      stat_file_name_suffix="stat",
                      best_pdb_name_suffix="model",
                      do_clean_first=True,
                      do_create_directories=True,
                      global_output_directory="./",
                      rmf_dir="rmfs/",
                      best_pdb_dir="pdbs/",
                      replica_stat_file_suffix="stat_replica",
                      em_object_for_rmf=None,
                      replica_exchange_object=None):

       self.model=model
       self.representation=representation
       self.crosslink_restraints=crosslink_restraints
       self.em_object_for_rmf=em_object_for_rmf
       self.sample_objects=sample_objects
       self.output_objects=output_objects
       self.replica_exchange_object=replica_exchange_object
       self.vars={}
       self.vars["monte_carlo_temperature"]=monte_carlo_temperature
       self.vars["replica_exchange_minimum_temperature"]=replica_exchange_minimum_temperature
       self.vars["replica_exchange_maximum_temperature"]=replica_exchange_maximum_temperature
       self.vars["number_of_best_scoring_models"]=number_of_best_scoring_models
       self.vars["monte_carlo_steps"]=monte_carlo_steps
       self.vars["number_of_frames"]=number_of_frames
       self.vars["write_initial_rmf"]=write_initial_rmf
       self.vars["initial_rmf_name_suffix"]=initial_rmf_name_suffix
       self.vars["best_pdb_name_suffix"]=best_pdb_name_suffix
       self.vars["stat_file_name_suffix"]=stat_file_name_suffix
       self.vars["do_clean_first"]=do_clean_first
       self.vars["do_create_directories"]=do_create_directories
       self.vars["global_output_directory"]=global_output_directory
       self.vars["rmf_dir"]=rmf_dir
       self.vars["best_pdb_dir"]=best_pdb_dir
       self.vars["replica_stat_file_suffix"]=replica_stat_file_suffix

    def show_info(self):
      print "ReplicaExchange0: it generates initial.*.rmf3, stat.*.out, rmfs/*.rmf3 for each replica "
      print "--- it stores the best scoring pdb models in pdbs/"
      print "--- the stat.*.out and rmfs/*.rmf3 are saved only at the lowest temperature"
      print "--- variables:"
      keys=self.vars.keys()
      keys.sort()
      for v in keys:
          print "------",v.ljust(30), self.vars[v]

    def get_replica_exchange_object(self):
        return self.replica_exchange_object

    def execute_macro(self):

      temp_index_factor=100000.0




      print "Setting up MonteCarlo"
      mc = IMP.pmi.samplers.MonteCarlo(self.model,
                                       self.sample_objects,
                                       self.vars["monte_carlo_temperature"])
      self.output_objects.append(mc)

# -------------------------------------------------------------------------

      print "Setting up ReplicaExchange"
      rex= IMP.pmi.samplers.ReplicaExchange(self.model,
                                            self.vars["replica_exchange_minimum_temperature"],
                                            self.vars["replica_exchange_maximum_temperature"],mc,
                                            replica_exchange_object=self.replica_exchange_object)
      self.replica_exchange_object=rex.rem

      myindex=rex.get_my_index()
      self.output_objects.append(rex)

      # must reset the minimum temperature due to the
      # different binary length of rem.get_my_parameter double and python float
      min_temp_index=int(min(rex.get_temperatures())*temp_index_factor)

# -------------------------------------------------------------------------

      globaldir=self.vars["global_output_directory"]+"/"
      rmf_dir=globaldir+self.vars["rmf_dir"]
      pdb_dir=globaldir+self.vars["best_pdb_dir"]

      if self.vars["do_clean_first"]:
        #to write
        pass
      if self.vars["do_create_directories"]:
        try:
            os.makedirs(globaldir)
        except:
            pass

        try:
            os.makedirs(rmf_dir)
        except:
            pass
        
        try:
            os.makedirs(pdb_dir)
        except:
            pass


# -------------------------------------------------------------------------

      sw = IMP.pmi.tools.Stopwatch()
      self.output_objects.append(sw)

      print "Setting up stat file"
      output = IMP.pmi.output.Output()
      low_temp_stat_file=globaldir+self.vars["stat_file_name_suffix"]+"."+str(myindex)+".out"
      output.init_stat2(low_temp_stat_file,
                  self.output_objects,
                  extralabels=["rmf_file","rmf_frame_index"])

      print "Setting up replica stat file"
      replica_stat_file=globaldir+self.vars["replica_stat_file_suffix"]+"."+str(myindex)+".out"
      output.init_stat2(replica_stat_file,[rex],extralabels=["score"])

      print "Setting up best pdb files"
      output.init_pdb_best_scoring(pdb_dir+"/"+
                                   self.vars["best_pdb_name_suffix"],
                                   self.representation.prot,
                                   self.vars["number_of_best_scoring_models"],
                                   replica_exchange=True)

#----------------------------------------------

      print "Setting up and writing initial rmf coordinate file"
      init_suffix=globaldir+self.vars["initial_rmf_name_suffix"]
      output.init_rmf(init_suffix+"."+str(myindex)+".rmf3",
                      [self.representation.prot])
      if self.crosslink_restraints:
         output.add_restraints_to_rmf(init_suffix+"."+str(myindex)+".rmf3",
                                      self.crosslink_restraints)
      output.write_rmf(init_suffix+"."+str(myindex)+".rmf3")
      output.close_rmf(init_suffix+"."+str(myindex)+".rmf3")
      
#----------------------------------------------
      
      print "Setting up production rmf files"
      if self.em_object_for_rmf!=None:
         output_hierarchies=[self.representation.prot,self.em_object_for_rmf.get_density_as_hierarchy()]
      else:
         output_hierarchies=[self.representation.prot]
      
      rmfname=rmf_dir+"/"+str(myindex)+".rmf3"
      output.init_rmf(rmfname, output_hierarchies)

      if self.crosslink_restraints:
         output.add_restraints_to_rmf(rmfname,self.crosslink_restraints)

      ntimes_at_low_temp=0

      if myindex==0:
         self.show_info()

      for i in range(self.vars["number_of_frames"]):

        mc.optimize(self.vars["monte_carlo_steps"])
        score=self.model.evaluate(False)
        output.set_output_entry("score",score)

        my_temp_index=int(rex.get_my_temp()*temp_index_factor)

        if min_temp_index==my_temp_index:
           print "--- frame %s score %s " % (str(i),str(score))

           output.write_pdb_best_scoring(score)
           output.write_rmf(rmfname)

           output.set_output_entry("rmf_file",rmfname)
           output.set_output_entry("rmf_frame_index",ntimes_at_low_temp)
           output.write_stat2(low_temp_stat_file)
           ntimes_at_low_temp+=1

        output.write_stat2(replica_stat_file)
        rex.swap_temp(i,score)






class ReplicaExchange1():
    '''
    Multiple state variant of the above class
    '''

    def __init__(self,model,
                      representations,
                      sample_objects,
                      output_objects,
                      crosslink_restraints=None,
                      monte_carlo_temperature=1.0,
                      replica_exchange_minimum_temperature=1.0,
                      replica_exchange_maximum_temperature=2.5,
                      number_of_best_scoring_models=500,
                      monte_carlo_steps=10,
                      number_of_frames=1000,
                      write_initial_rmf=True,
                      initial_rmf_name_suffix="initial",
                      stat_file_name_suffix="stat",
                      best_pdb_name_suffix="model",
                      do_clean_first=True,
                      do_create_directories=True,
                      rmf_dir="rmfs/",
                      best_pdb_dir="pdbs/"):

       self.model=model
       self.representations=representations
       self.crosslink_restraints=crosslink_restraints
       self.sample_objects=sample_objects
       self.output_objects=output_objects
       self.vars={}
       self.vars["number_of_states"]=len(self.representations)
       self.vars["monte_carlo_temperature"]=monte_carlo_temperature
       self.vars["replica_exchange_minimum_temperature"]=replica_exchange_minimum_temperature
       self.vars["replica_exchange_maximum_temperature"]=replica_exchange_maximum_temperature
       self.vars["number_of_best_scoring_models"]=number_of_best_scoring_models
       self.vars["monte_carlo_steps"]=monte_carlo_steps
       self.vars["number_of_frames"]=number_of_frames
       self.vars["write_initial_rmf"]=write_initial_rmf
       self.vars["initial_rmf_name_suffix"]=initial_rmf_name_suffix
       self.vars["best_pdb_name_suffix"]=best_pdb_name_suffix
       self.vars["stat_file_name_suffix"]=stat_file_name_suffix
       self.vars["do_clean_first"]=do_clean_first
       self.vars["do_create_directories"]=do_create_directories
       self.vars["rmf_dir"]=rmf_dir
       self.vars["best_pdb_dir"]=best_pdb_dir


    def show_info(self):
      print "ReplicaExchange0: it generates initial.*.rmf3, stat.*.out, rmfs/*.rmf3 for each replica "
      print "--- it stores the best scoring pdb models in pdbs/"
      print "--- the stat.*.out and rmfs/*.rmf3 are saved only at the lowest temperature"
      print "--- variables:"
      keys=self.vars.keys()
      keys.sort()
      for v in keys:
          print "------",v.ljust(30), self.vars[v]

    def execute_macro(self):


      temp_index_factor=100000.0

      if self.vars["do_clean_first"]:
        #to write
        pass


      if self.vars["do_create_directories"]:
        for n in range(self.vars["number_of_states"]):
          if not os.path.exists(self.vars["rmf_dir"]+str(n)):
             os.makedirs(self.vars["rmf_dir"]+str(n))
          if not os.path.exists(self.vars["best_pdb_dir"]+str(n)):
             os.makedirs(self.vars["best_pdb_dir"]+str(n))


      print "Setting up MonteCarlo"
      mc = IMP.pmi.samplers.MonteCarlo(self.model,
                                       self.sample_objects,
                                       self.vars["monte_carlo_temperature"])
      self.output_objects.append(mc)


      print "Setting up ReplicaExchange"
      rex= IMP.pmi.samplers.ReplicaExchange(self.model,
         self.vars["replica_exchange_minimum_temperature"],
         self.vars["replica_exchange_maximum_temperature"],mc)

      myindex=rex.get_my_index()
      self.output_objects.append(rex)

      # must reset the minimum temperature due to the
      # different binary length of rem.get_my_parameter double and python float
      min_temp_index=int(min(rex.get_temperatures())*temp_index_factor)

      sw = IMP.pmi.tools.Stopwatch()
      self.output_objects.append(sw)

      print "Setting up stat file"
      output = IMP.pmi.output.Output()
      output.init_stat2(self.vars["stat_file_name_suffix"]+"."+str(myindex)+".out",
                  self.output_objects,
                  extralabels=["rmf_file","rmf_frame_index"])

      print "Setting up best pdb files"
      for n in range(self.vars["number_of_states"]):
         output.init_pdb_best_scoring(self.vars["best_pdb_dir"]+str(n)+"/"+
                                   self.vars["best_pdb_name_suffix"],
                                   self.representations[n].prot,
                                   self.vars["number_of_best_scoring_models"],
                                   replica_exchange=True)

      print "Setting up and writing initial rmf coordinate file"
      init_suffix=self.vars["initial_rmf_name_suffix"]
      output.init_rmf(init_suffix+"."+str(myindex)+".rmf3",
                      [r.prot for r in self.representations])
      if self.crosslink_restraints:
         output.add_restraints_to_rmf(init_suffix+"."+str(myindex)+".rmf3",
                                      self.crosslink_restraints)
      output.write_rmf(init_suffix+"."+str(myindex)+".rmf3")
      output.close_rmf(init_suffix+"."+str(myindex)+".rmf3")

      print "Setting up production rmf files"
      rmfdir=self.vars["rmf_dir"]
      rmfname=rmfdir+"/"+str(myindex)+".rmf3"
      output.init_rmf(rmfname, [r.prot for r in self.representations])

      if self.crosslink_restraints:
         output.add_restraints_to_rmf(rmfname,self.crosslink_restraints)

      ntimes_at_low_temp=0

      if myindex==0:
         self.show_info()

      for i in range(self.vars["number_of_frames"]):

        mc.optimize(self.vars["monte_carlo_steps"])
        score=self.model.evaluate(False)
        my_temp_index=int(rex.get_my_temp()*temp_index_factor)

        if min_temp_index==my_temp_index:
           print "--- frame %s score %s " % (str(i),str(score))

           output.write_pdb_best_scoring(score)
           output.write_rmf(rmfname)

           output.set_output_entry("rmf_file",rmfname)
           output.set_output_entry("rmf_frame_index",ntimes_at_low_temp)
           output.write_stats2()
           ntimes_at_low_temp+=1

        rex.swap_temp(i,score)


def BuildModel0(m,data,resolutions=[1,10],missing_bead_size=20):
      '''
      The macro construct a component for each subunit (no splitting, nothing fancy)
      You can pass the resolutions and the bead size for the missing residue regions.
      To use this macro, you must provide the following data structure:
      Component  pdbfile    chainid  rgb color     fastafile     sequence id
                                                                        in fastafile
data = [("Rpb1",     pdbfile,   "A",     0.00000000,  (fastafile,    0)),
        ("Rpb2",     pdbfile,   "B",     0.09090909,  (fastafile,    1)),
        ("Rpb3",     pdbfile,   "C",     0.18181818,  (fastafile,    2)),
        ("Rpb4",     pdbfile,   "D",     0.27272727,  (fastafile,    3)),
        ("Rpb5",     pdbfile,   "E",     0.36363636,  (fastafile,    4)),
        ("Rpb6",     pdbfile,   "F",     0.45454545,  (fastafile,    5)),
        ("Rpb7",     pdbfile,   "G",     0.54545455,  (fastafile,    6)),
        ("Rpb8",     pdbfile,   "H",     0.63636364,  (fastafile,    7)),
        ("Rpb9",     pdbfile,   "I",     0.72727273,  (fastafile,    8)),
        ("Rpb10",    pdbfile,   "L",     0.81818182,  (fastafile,    9)),
        ("Rpb11",    pdbfile,   "J",     0.90909091,  (fastafile,   10)),
        ("Rpb12",    pdbfile,   "K",     1.00000000,  (fastafile,   11))]
      '''

      r=IMP.pmi.representation.Representation(m)
      hierarchies={}

      for d in data:
            component_name=d[0]
            pdb_file=d[1]
            chain_id=d[2]
            color_id=d[3]
            fasta_file=d[4][0]
            fastids=IMP.pmi.tools.get_ids_from_fasta_file(fasta_file)
            fasta_file_id=d[4][1]
            #avoid to add a component with the same name
            r.create_component(component_name,
                               color=color_id)

            r.add_component_sequence(component_name,
                                     fasta_file,
                                     id=fastids[fasta_file_id])

            hierarchies=r.autobuild_model(component_name,
                                           pdb_file,
                                           chain_id,
                                           resolutions=resolutions,
                                           missingbeadsize=missing_bead_size)

            r.show_component_table(component_name)

            r.set_rigid_bodies([component_name])

            r.set_chain_of_super_rigid_bodies(hierarchies,min_length=2,max_length=2)

            r.setup_component_sequence_connectivity(component_name,resolution=1)
            r.setup_component_geometry(component_name)


      r.setup_bonds()
      #put it at the end of rigid bodies
      r.set_floppy_bodies()

      #set current coordinates as reference for RMSD calculation
      r.set_current_coordinates_as_reference_for_rmsd("Reference")

      return r





class AnalysisReplicaExchange0():


    
    def __init__(self,model,
                      stat_file_name_suffix="stat",
                      best_pdb_name_suffix="model",
                      do_clean_first=True,
                      do_create_directories=True,
                      global_output_directory="./",
                      rmf_dir="rmfs/",
                      best_pdb_dir="pdbs/",
                      replica_stat_file_suffix="stat_replica",
                      global_analysis_result_directory="./analysis/"):

                      import glob
                      self.model=model
                      rmf_dir=global_output_directory+rmf_dir
                      self.rmf_files=glob.glob(rmf_dir+"/*.rmf")
                      stat_dir=global_output_directory
                      self.stat_files=glob.glob(stat_dir+"/stat.*.out")
                      
                      
    def clustering(self,score_key,
                        rmf_file_key,
                        rmf_file_frame_key,
                        feature_keys=None,
                        alignment_components=None,
                        number_of_best_scoring_models=10,
                        rmsd_calculation_components=None,
                        distance_matrix_file=None,
                        load_distance_matrix_file=False,
                        is_mpi=False,
                        number_of_clusters=1,
                        display_plot=False,
                        get_every=1,
                        density_custom_ranges=None):     
        
        '''
        the features are keywords for which you want to calculate average, medians, etc,
        plots.
        
        density_calculation_components is a list of tuples or strings:
                                      ["Rpb1", (20,100,"Rpb2"), .....]
                                      if that is None, don't calculate the density.
        
        
        '''
        
        from operator import itemgetter
        import IMP.pmi.analysis
        import numpy as np
        
        if is_mpi:
           from mpi4py import MPI
           comm = MPI.COMM_WORLD
           rank = comm.Get_rank()
           number_of_processes=comm.size
        else:
           rank=0
           number_of_processes=1

    
        if alignment_components==None:
           alignment_flag=0
        else:
           alignment_flag=1 
        
        if feature_keys==None: feature_keys=[]          
        
        print "setup clustering class"
        Clusters = IMP.pmi.analysis.Clustering()
    
        if not load_distance_matrix_file:    
            my_stat_files=IMP.pmi.tools.chunk_list_into_segments(self.stat_files,number_of_processes)[rank]

            score_list=[]
            rmf_file_list=[]
            rmf_file_frame_list=[]
            if score_key!=None:
               feature_keyword_list_dict={}
               for sf in my_stat_files:
                  print "getting data from file %s" % sf
                  po=IMP.pmi.output.ProcessOutput(sf)
                  keywords=po.get_keys()
                  feature_keywords=[score_key,rmf_file_key,rmf_file_frame_key]
                  for k in keywords:
                      for fk in feature_keys:
                          if fk in k: feature_keywords.append(k)
                  fields=po.get_fields(feature_keywords,get_every=get_every)
                  #check that all lengths are all equal
                  length_set=set()
                  for f in fields: length_set.add(len(fields[f]))
                  
                  # if some of the fields are missing, truncate
                  # the feature files to the shortest one
                  
                  if len(length_set)>1: 
                     print "AnalysisReplicaExchange0.clustering: the statfile is not synchronous"
                     minlen=min(length_set)
                     for f in fields: 
                         fields[f]=fields[f][0:minlen]
                     
                  # append to the lists
                  score_list+=fields[score_key]
                  rmf_file_list+=fields[rmf_file_key]
                  rmf_file_frame_list+=fields[rmf_file_frame_key] 
                                    
                  for k in feature_keywords:
                      if k in feature_keyword_list_dict:
                         feature_keyword_list_dict[k]+=fields[k]
                      else:
                         feature_keyword_list_dict[k]=fields[k]
                    

# --------------------------------------------------------------------------------------------
# broadcast the conformation features
        
            if number_of_processes>1:
                    score_list=IMP.pmi.tools.scatter_and_gather(score_list)
                    rmf_file_list=IMP.pmi.tools.scatter_and_gather(rmf_file_list)
                    rmf_file_frame_list=IMP.pmi.tools.scatter_and_gather(rmf_file_frame_list)  
                    for k in feature_keyword_list_dict:
                        feature_keyword_list_dict[k]=IMP.pmi.tools.scatter_and_gather(feature_keyword_list_dict[k])
            
            # sort by score and ge the best scoring ones
            score_rmf_tuples=zip(score_list,rmf_file_list,rmf_file_frame_list,range(len(score_list)))
            #numerically sorting in ascending order
            sorted_score_rmf_tuples=sorted(score_rmf_tuples,key=lambda x: float(x[0]))
            best_score_rmf_tuples=sorted_score_rmf_tuples[0:number_of_best_scoring_models]
            # we have to define a best score rank when we process in parallel
            # to get the correct position in the best_score_rmf_tuples
            for n,tpl in enumarate(best_score_rmf_tuples):
                tpl_new=tuple(list(tpl).append(n))
                best_score_rmf_tuples[n]=tpl_new
            
# prune the feature_keyword_list_dict, keep only the best scoring models
# and sort it
            
            best_score_feature_keyword_list_dict={}
            for tpl in best_score_rmf_tuples:
                
                index=tpl[3]
                for f in feature_keyword_list_dict:
                    if f in best_score_feature_keyword_list_dict:
                       best_score_feature_keyword_list_dict[f].append(feature_keyword_list_dict[f][index])
                    else:
                       best_score_feature_keyword_list_dict[f]=[feature_keyword_list_dict[f][index]]
            
            del feature_keyword_list_dict
            
            
            # here I've tested that feature_keyword_list_dict is correct on 2 CPUs
            
# --------------------------------------------------------------------------------------------        
# read the coordinates
# --------------------------------------------------------------------------------------------


            my_best_score_rmf_tuples=IMP.pmi.tools.chunk_list_into_segments(best_score_rmf_tuples,
                                                                            number_of_processes)[rank]
            # reading the coordinates
            
            all_coordinates=[]
            all_rmf_file_names=[]
            # it will be used to extract the features
            rmf_file_name_index_dict={}
            for cnt,tpl in enumerate(my_best_score_rmf_tuples):
                rmf_file=tpl[1]
                frame_number=tpl[2]
                
                prot=IMP.pmi.analysis.get_hier_from_rmf(self.model,frame_number,rmf_file)
                if not prot: continue
                # getting the particles                
                part_dict=IMP.pmi.analysis.get_particles_at_resolution_one(prot)
                
                # getting the coordinates
                model_coordinate_dict={}
                for pr in part_dict:
                    coords = np.array([np.array(IMP.core.XYZ(i).get_coordinates()) for i in part_dict[pr]])
                    model_coordinate_dict[pr] = coords 

                all_coordinates.append(model_coordinate_dict)
                frame_name=rmf_file+'|'+str(frame_number)
                all_rmf_file_names.append(frame_name)
                rmf_file_name_index_dict[frame_name]=tpl[4]


# --------------------------------------------------------------------------------------------
# broadcast the coordinates

            if number_of_processes>1:
                all_coordinates=IMP.pmi.tools.scatter_and_gather(all_coordinates)
                all_rmf_file_names=IMP.pmi.tools.scatter_and_gather(all_rmf_file_names)
                rmf_file_name_index_dict=IMP.pmi.tools.scatter_and_gather(rmf_file_name_index_dict)

            if rank==0:
               # save needed informations in external files
               self.save_objects([best_score_feature_keyword_list_dict,rmf_file_name_index_dict],".macro.pkl")


# -----------------------------------------------------------------------------------------------              

            for n,model_coordinate_dict in enumerate(all_coordinates):
                
                template_coordinate_dict= {}
                # let's try to align
                if alignment_flag==1 and len(Clusters.all_coords)==0:
                   for pr in alignment_components: 
                       template_coordinate_dict[pr] = model_coordinate_dict[pr]
                # set the first model as template coordinates
                   Clusters.set_template(template_coordinate_dict)
      
                # set particles to calculate RMSDs on
                # (if global, list all proteins, or just a few for local) 
                # setting all the coordinates for rmsd calculation
                rmsd_coordinate_dict = {}
                if rmsd_calculation_components==None: rmsd_calculation_components=alignment_components
                
                for pr in rmsd_calculation_components:
                    rmsd_coordinate_dict[pr] = model_coordinate_dict[pr]

                Clusters.fill(all_rmf_file_names[n], rmsd_coordinate_dict)
    
            print "Global calculating the distance matrix"
            
            # calculate distance matrix, all against all
            
            Clusters.dist_matrix(is_mpi=is_mpi)

            Clusters.do_cluster(number_of_clusters)
            if display_plot:
               if rank==0:
                  Clusters.plot_matrix()
            Clusters.save_distance_matrix_file(file_name=distance_matrix_file)               
          
# -----------------------------------------------------------------------------------------------           
# load the distance matrix from file

     
        else:
            Clusters.load_distance_matrix_file(file_name=distance_matrix_file)
            Clusters.do_cluster(number_of_clusters)
            [best_score_feature_keyword_list_dict,rmf_file_name_index_dict]=self.load_objects(".macro.pkl")
            if display_plot:
               if rank==0:
                  Clusters.plot_matrix()           

# -----------------------------------------------------------------------------------------------  
# now save all informations about the clusters

# -----------------------------------------------------------------------------------------------
#
        if rank==0:

          o=IMP.pmi.output.Output()
        
          for n,cl in enumerate(Clusters.get_cluster_labels()):
              
              # first initialize the Density class if requested
              
              if density_custom_ranges:
                DensModule = IMP.pmi.analysis.GetModelDensity(density_custom_ranges)
              
              print Clusters.get_cluster_label_average_rmsd(cl)
              print Clusters.get_cluster_label_size(cl)
              print Clusters.get_cluster_label_names(cl)
              dircluster="cluster."+str(n)+"/"
              try:
                 os.mkdir(dircluster)
              except:
                 pass 
              
              clusstat=open(dircluster+"stat.out","w")   
                              
              for k,structure_name in enumerate(Clusters.get_cluster_label_names(cl)):
                  
                  
                  #extract the features
                  tmp_dict={}
                  index=rmf_file_name_index_dict[structure_name]
                  for key in best_score_feature_keyword_list_dict:
                      tmp_dict[key]=best_score_feature_keyword_list_dict[key][index]
                  
                  rmf_name=structure_name.split("|")[0]
                  rmf_frame_number=int(structure_name.split("|")[1])
                  
                  clusstat.write(str(tmp_dict)+"\n")                  
                  
                  prot=IMP.pmi.analysis.get_hier_from_rmf(self.model,rmf_frame_number,rmf_name)
                  if not prot: continue
                  
                  if k>0:
                      model_index=Clusters.get_model_index_from_name(structure_name)
                      transformation=Clusters.get_transformation_to_first_member(cl,model_index)
                      
                      rbs=[]           
                      for p in IMP.atom.get_leaves(prot):              
                          if IMP.core.RigidBody.particle_is_instance(p):                 
                             rb=IMP.core.RigidMember(p).get_rigid_body()                     
                             if rb not in rbs:                      
                                rbs.append(rb)              
                                IMP.core.transform(rb,transformation)
                                
                          else:                 
                             IMP.core.transform(IMP.core.XYZ(p),transformation)              
                        #IMP.em.add_to_map(dmap_dict[name],particle_dict[name])

                  # add the density                  
                  if density_custom_ranges:
                     DensModule.add_subunits_density(prot)

                  o.init_pdb(dircluster+str(k)+".pdb",prot)        
                  o.write_pdb(dircluster+str(k)+".pdb")
                  o.init_rmf(dircluster+str(k)+".rmf3",[prot])   
                  #IMP.rmf.add_restraints(o.dictionary_rmfs[dircluster+str(n)+".rmf3"],rs)    
                  o.write_rmf(dircluster+str(k)+".rmf3")
                  o.close_rmf(dircluster+str(k)+".rmf3")     
              
              if density_custom_ranges:
                 DensModule.write_mrc(path=dircluster)
                 del DensModule


    def save_objects(self,objects,file_name):
        import pickle    
        outf = open(file_name,'w')
        pickle.dump(objects,outf)            
        outf.close()        

    def load_objects(self,file_name):
        import pickle        
        inputf = open(file_name,'r')    
        objects=pickle.load(inputf)         
        inputf.close()
        return objects














