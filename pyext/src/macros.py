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
                      replica_stat_file_suffix="stat_replica" ):

       self.model=model
       self.representation=representation
       self.crosslink_restraints=crosslink_restraints
       self.sample_objects=sample_objects
       self.output_objects=output_objects
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

    def execute_macro(self):




      if self.vars["do_clean_first"]:
        #to write
        pass

      globaldir=self.vars["global_output_directory"]+"/"
      rmf_dir=globaldir+self.vars["rmf_dir"]
      pdb_dir=globaldir+self.vars["best_pdb_dir"]


      if self.vars["do_create_directories"]:
        if not os.path.exists(globaldir): os.makedirs(globaldir)
        if not os.path.exists(rmf_dir): os.makedirs(rmf_dir)
        if not os.path.exists(pdb_dir): os.makedirs(pdb_dir)


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
      
      # must reset the minimum and maximum temperature due to the 
      # different binary length of rem.get_my_parameter and python float
      self.vars["replica_exchange_minimum_temperature"]=rex.rem.get_my_parameter("temp")[0]
      self.vars["replica_exchange_maximum_temperature"]=rex.rem.get_my_parameter("temp")[-1]

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

      print "Setting up and writing initial rmf coordinate file"
      init_suffix=globaldir+self.vars["initial_rmf_name_suffix"]
      output.init_rmf(init_suffix+"."+str(myindex)+".rmf3",
                      [self.representation.prot])
      if self.crosslink_restraints:
         output.add_restraints_to_rmf(init_suffix+"."+str(myindex)+".rmf3",
                                      self.crosslink_restraints)
      output.write_rmf(init_suffix+"."+str(myindex)+".rmf3")
      output.close_rmf(init_suffix+"."+str(myindex)+".rmf3")

      print "Setting up production rmf files"
      rmfname=rmf_dir+"/"+str(myindex)+".rmf3"
      output.init_rmf(rmfname, [self.representation.prot])

      if self.crosslink_restraints:
         output.add_restraints_to_rmf(rmfname,self.crosslink_restraints)

      ntimes_at_low_temp=0

      if myindex==0:
         self.show_info()

      for i in range(self.vars["number_of_frames"]):

        mc.optimize(self.vars["monte_carlo_steps"])
        score=self.model.evaluate(False)
        output.set_output_entry("score",score)
        if rex.get_my_temp()==self.vars["replica_exchange_minimum_temperature"]:
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
      
      # must reset the minimum and maximum temperature due to the 
      # different binary length of rem.get_my_parameter and python float
      self.vars["replica_exchange_minimum_temperature"]=rex.rem.get_my_parameter("temp")[0]
      self.vars["replica_exchange_maximum_temperature"]=rex.rem.get_my_parameter("temp")[-1]

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
        if rex.get_my_temp()==self.vars["replica_exchange_minimum_temperature"]:
           print "--- frame %s score %s " % (str(i),str(score))

           output.write_pdb_best_scoring(score)
           output.write_rmf(rmfname)

           output.set_output_entry("rmf_file",rmfname)
           output.set_output_entry("rmf_frame_index",ntimes_at_low_temp)
           output.write_stats2()
           ntimes_at_low_temp+=1

        rex.swap_temp(i,score)


class BuildModel0():
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


    def __init__(self,data,resolutions=[1,10],missing_bead_size=20):

      self.r=IMP.pmi.representation.Representation(m)
      self.data=data
      hierarchies={}

      for d in self.data:
            component_name=d[0]
            pdb_file=d[1]
            chain_id=d[2]
            color_id=d[3]
            fasta_file=d[4][0]
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
