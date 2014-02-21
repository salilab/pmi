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

      if myindex==0:
        if self.vars["do_clean_first"]:
          #to write
          pass
        if self.vars["do_create_directories"]:
          if not os.path.exists(globaldir): os.makedirs(globaldir)
          if not os.path.exists(rmf_dir): os.makedirs(rmf_dir)
          if not os.path.exists(pdb_dir): os.makedirs(pdb_dir)

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
                        alignment_components=None,
                        number_of_best_scoring_models=10,
                        rmsd_calculation_components=None,
                        distance_matrix_file=None,
                        load_distance_matrix_file=False,
                        number_of_processes=1):
                        
        from operator import itemgetter
        import IMP.pmi.analysis
        import IMP.rmf
        import RMF
        import numpy as np


        if alignment_components==None:
           alignment_flag=0
        else:
           alignment_flag=1           
        
        print "setup clustering class"
        Clusters = IMP.pmi.analysis.Clustering()


      
        score_list=[]
        rmf_file_list=[]
        rmf_file_frame_list=[]
        if score_key!=None:
           for sf in self.stat_files[0:1]:
              print "getting data from file %s" % sf
              po=IMP.pmi.output.ProcessOutput(sf)
              fields=po.get_fields([score_key,rmf_file_key,rmf_file_frame_key])
              #check that all lengths are all equal
              length_set=set()
              for f in fields: length_set.add(len(fields[f]))
              if len(length_set)>1: print "AnalysisReplicaExchange0.clustering: the statfile is not synchronous"; exit()
              # append to the lists
              score_list+=fields[score_key]
              rmf_file_list+=fields[rmf_file_key]
              rmf_file_frame_list+=fields[rmf_file_frame_key]                            
        
        # sort by score and ge the best scoring ones
        score_rmf_tuples=zip(score_list,rmf_file_list,rmf_file_frame_list)
        sorted_score_rmf_tuples=sorted(score_rmf_tuples,key=itemgetter(0))
        best_score_rmf_tuples=sorted_score_rmf_tuples[0:number_of_best_scoring_models]
        
        
        
        if not load_distance_matrix_file:          
          for cnt,tpl in enumerate(best_score_rmf_tuples):

              rmf_file=tpl[1]
              frame_number=tpl[2]
              print "getting coordinates for frame %i rmf file %s" % (frame_number,rmf_file)
                          
              # load the frame
              rh= RMF.open_rmf_file_read_only(rmf_file)
              prot=IMP.rmf.create_hierarchies(rh, self.model)[0]     
              IMP.rmf.link_hierarchies(rh, [prot])
              try: 
                IMP.rmf.load_frame(rh, frame_number)        
              except:
                print "Unable to open frame %i of file %s" % (frame_number,rmf_file)
                continue           
              self.model.update()        
              
              Clusters.set_prot(prot)
              
              # let's get all particles at highest resolution
              part_dict=IMP.pmi.analysis.get_particles_at_resolution_one(prot)
              
              template_coordinate_dict= {}
              
              # let's try to align
              if alignment_flag==1:
                 for pr in alignment_components:
                     coords = np.array([np.array(IMP.core.XYZ(i).get_coordinates()) for i in part_dict[pr]]) 
                     template_coordinate_dict[pr] = coords

              # set the first model as template coordinates
              if len(Clusters.all_coords)==0:
                  Clusters.set_template(template_coordinate_dict)

              # set particles to calculate RMSDs on
              # (if global, list all proteins, or just a few for local) 
              # setting all the coordinates for rmsd calculation
                
              model_coordinate_dict = {}
              if rmsd_calculation_components==None: rmsd_calculation_components=alignment_components
              
              for pr in rmsd_calculation_components:
                  coords = np.array([np.array(IMP.core.XYZ(i).get_coordinates()) for i in part_dict[pr]])
                  model_coordinate_dict[pr] = coords 
                 
                  
              Clusters.fill(rmf_file+'|'+str(frame_number), model_coordinate_dict)

          print "Global calculating the distance matrix"
          
          # calculate distance matrix, all against all
          
          Clusters.dist_matrix(file_name=distance_matrix_file,number_of_processes=number_of_processes)
          Clusters.plot_matrix()
           
        else:
          Clusters.load_distance_matrix_file(file_name=distance_matrix_file)
          Clusters.plot_matrix()           
        
        #this list 
        for cl in Clusters.get_cluster_labels():
          print Clusters.get_cluster_label_average_rmsd(cl)
          print Clusters.get_cluster_label_size(cl)
          print Clusters.get_cluster_label_names(cl)
         
    
    
    
    
