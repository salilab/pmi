class ReplicaExchange_0():
    
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
                      rmf_dir="rmfs/",
                      best_pdb_dir="pdbs/"):

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
       self.vars["rmf_dir"]=rmf_dir
       self.vars["best_pdb_dir"]=best_pdb_dir           
       

    def show_info(self):
      print "ReplicaExchange_0: it generates initial.*.rmf3, stat.*.out, rmfs/*.rmf3 for each replica "
      print "--- it stores the best scoring pdb models in pdbs/"
      print "--- the stat.*.out and rmfs/*.rmf3 are saved only at the lowest temperature"      
      print "--- variables:"
      keys=self.vars.keys()
      keys.sort()
      for v in keys: 
          print "------",v, self.vars[v]
       
    def execute_macro(self):
      import os
      
      self.show_info()
      
      if self.vars["do_clean_first"]:
        #to write
        pass
      
      
      if self.vars["do_create_directories"]:
        if not os.path.exists(self.vars["rmf_dir"]):
           os.makedirs(self.vars["rmf_dir"])
        if not os.path.exists(self.vars["best_pdb_dir"]):
           os.makedirs(self.vars["best_pdb_dir"])           
      
          
      mc = IMP.pmi.samplers.MonteCarlo(self.model,
                                       self.sample_objects,
                                       self.vars["monte_carlo_temperature"])
      self.output_objects.append(mc)

      rex= IMP.pmi.samplers.ReplicaExchange(self.model,
         self.vars["replica_exchange_minimum_temperature"],
         self.vars["replica_exchange_maximum_temperature"],mc)
      
      myindex=rex.get_my_index()
      self.output_objects.append(rex)

      sw = IMP.pmi.tools.Stopwatch()
      self.output_objects.append(sw)

      output = IMP.pmi.output.Output()
      output.init_stat2(self.vars["stat_file_name_suffix"]+"."+str(myindex)+".out", 
                  self.output_objects, 
                  extralabels=["rmf_file","rmf_frame_index"])
      
      output.init_pdb_best_scoring(self.vars["best_pdb_dir"]+"/"+
                                   self.vars["best_pdb_name_suffix"],
                                   self.representation.prot,
                                   self.vars["number_of_best_scoring_models"],
                                   replica_exchange=True)
      
      init_suffix=self.vars["initial_rmf_name_suffix"]                         
      output.init_rmf(init_suffix+"."+str(myindex)+".rmf3", 
                      [self.representation.prot])
      output.write_rmf(init_suffix+"."+str(myindex)+".rmf3")
      if self.crosslink_restraints:
         output.add_restraints_to_rmf(init_suffix+"."+str(myindex)+".rmf3",
                                      self.crosslink_restraints)
      output.close_rmf(init_suffix+"."+str(myindex)+".rmf3")

      rmfdir=self.vars["rmf_dir"]
      rmfname=rmfdir+"/"+str(myindex)+".rmf3"
      output.init_rmf(rmfname, [self.representation.prot])
      
      if self.crosslink_restraints:
         output.add_restraints_to_rmf(rmfname,self.crosslink_restraints)

      ntimes_at_low_temp=0

      for i in range(self.vars["number_of_frames"]):

        mc.optimize(self.vars["monte_carlo_steps"])
        score=m.evaluate(False)
        if rex.get_my_temp()==1.0:
           output.write_pdb_best_scoring(score)
           output.write_rmf(rmfname)

           output.set_output_entry("rmf_file",rmfname)
           output.set_output_entry("rmf_frame_index",ntimes_at_low_temp)
           output.write_stats2()       
           ntimes_at_low_temp+=1    

        rex.swap_temp(i,score)
