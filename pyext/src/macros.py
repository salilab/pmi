"""@namespace IMP.pmi.macros
   Set up of protocols.
"""

import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import os


class ReplicaExchange0(object):

    def __init__(self, model,
                 representation=None,
                 root_hier=None,
                 sample_objects=None, # DEPRECATED
                 monte_carlo_sample_objects=None,
                 molecular_dynamics_sample_objects=None,
                 output_objects=None,
                 crosslink_restraints=None,
                 monte_carlo_temperature=1.0,
                 replica_exchange_minimum_temperature=1.0,
                 replica_exchange_maximum_temperature=2.5,
                 num_sample_rounds=1,
                 number_of_best_scoring_models=500,
                 monte_carlo_steps=10,
                 molecular_dynamics_steps=10,
                 number_of_frames=1000,
                 nframes_write_coordinates=1,
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
                 atomistic=False,
                 replica_exchange_object=None):
        '''
        representation    IMP.pmi.Representation()      can be a list of representations (enables the multi state modeling)
        '''

        self.model = model
        self.vars = {}

        ### add check hierarchy is multistate
        if representation:
            if type(representation) == list:
                self.is_multi_state = True
                self.root_hiers = [r.prot for r in representation]
                self.vars["number_of_states"] = len(representation)
            else:
                self.is_multi_state = False
                self.root_hier = representation.prot
                self.vars["number_of_states"] = 1
        elif root_hier:
            states = IMP.atom.get_by_type(root_hier,IMP.atom.STATE_TYPE)
            self.vars["number_of_states"] = len(states)
            if len(states)>1:
                self.root_hiers = states
                self.is_multi_state = True
            else:
                self.root_hier = root_hier
                self.is_multi_state = False
        else:
            print "ERROR: Must provide representation or root_hier"
            return

        self.crosslink_restraints = crosslink_restraints
        self.em_object_for_rmf = em_object_for_rmf
        self.monte_carlo_sample_objects = monte_carlo_sample_objects
        if sample_objects is not None:
            self.monte_carlo_sample_objects+=sample_objects
        self.molecular_dynamics_sample_objects=molecular_dynamics_sample_objects
        self.output_objects = output_objects
        self.replica_exchange_object = replica_exchange_object
        self.vars["monte_carlo_temperature"] = monte_carlo_temperature
        self.vars[
            "replica_exchange_minimum_temperature"] = replica_exchange_minimum_temperature
        self.vars[
            "replica_exchange_maximum_temperature"] = replica_exchange_maximum_temperature
        self.vars["num_sample_rounds"] = num_sample_rounds
        self.vars[
            "number_of_best_scoring_models"] = number_of_best_scoring_models
        self.vars["monte_carlo_steps"] = monte_carlo_steps
        self.vars["molecular_dynamics_steps"]=molecular_dynamics_steps
        self.vars["number_of_frames"] = number_of_frames
        self.vars["nframes_write_coordinates"] = nframes_write_coordinates
        self.vars["write_initial_rmf"] = write_initial_rmf
        self.vars["initial_rmf_name_suffix"] = initial_rmf_name_suffix
        self.vars["best_pdb_name_suffix"] = best_pdb_name_suffix
        self.vars["stat_file_name_suffix"] = stat_file_name_suffix
        self.vars["do_clean_first"] = do_clean_first
        self.vars["do_create_directories"] = do_create_directories
        self.vars["global_output_directory"] = global_output_directory
        self.vars["rmf_dir"] = rmf_dir
        self.vars["best_pdb_dir"] = best_pdb_dir
        self.vars["atomistic"] = atomistic
        self.vars["replica_stat_file_suffix"] = replica_stat_file_suffix

    def show_info(self):
        print "ReplicaExchange0: it generates initial.*.rmf3, stat.*.out, rmfs/*.rmf3 for each replica "
        print "--- it stores the best scoring pdb models in pdbs/"
        print "--- the stat.*.out and rmfs/*.rmf3 are saved only at the lowest temperature"
        print "--- variables:"
        keys = self.vars.keys()
        keys.sort()
        for v in keys:
            print "------", v.ljust(30), self.vars[v]

    def get_replica_exchange_object(self):
        return self.replica_exchange_object

    def execute_macro(self):

        temp_index_factor = 100000.0
        samplers=[]
        sampler_mc=None
        sampler_md=None
        if self.monte_carlo_sample_objects is not None:
            print "Setting up MonteCarlo"
            sampler_mc = IMP.pmi.samplers.MonteCarlo(self.model,
                                                 self.monte_carlo_sample_objects,
                                                 self.vars["monte_carlo_temperature"])
            self.output_objects.append(sampler_mc)
            samplers.append(sampler_mc)
        if self.molecular_dynamics_sample_objects is not None:
            print "Setting up MolecularDynamics"
            sampler_md = IMP.pmi.samplers.MolecularDynamics(self.model,
                                                       self.molecular_dynamics_sample_objects,
                                                       self.vars["monte_carlo_temperature"])
            self.output_objects.append(sampler_md)
            samplers.append(sampler_md)
# -------------------------------------------------------------------------

        print "Setting up ReplicaExchange"
        rex = IMP.pmi.samplers.ReplicaExchange(self.model,
                                               self.vars[
                                                   "replica_exchange_minimum_temperature"],
                                               self.vars[
                                                   "replica_exchange_maximum_temperature"],
                                               samplers,
                                               replica_exchange_object=self.replica_exchange_object)
        self.replica_exchange_object = rex.rem

        myindex = rex.get_my_index()
        self.output_objects.append(rex)

        # must reset the minimum temperature due to the
        # different binary length of rem.get_my_parameter double and python
        # float
        min_temp_index = int(min(rex.get_temperatures()) * temp_index_factor)

# -------------------------------------------------------------------------

        globaldir = self.vars["global_output_directory"] + "/"
        rmf_dir = globaldir + self.vars["rmf_dir"]
        pdb_dir = globaldir + self.vars["best_pdb_dir"]

        if self.vars["do_clean_first"]:
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

            if not self.is_multi_state:
                try:
                    os.makedirs(pdb_dir)
                except:
                    pass
            else:
                for n in range(self.vars["number_of_states"]):
                    try:
                        os.makedirs(pdb_dir + "/" + str(n))
                    except:
                        pass

# -------------------------------------------------------------------------

        sw = IMP.pmi.tools.Stopwatch()
        self.output_objects.append(sw)

        print "Setting up stat file"
        output = IMP.pmi.output.Output(atomistic=self.vars["atomistic"])
        low_temp_stat_file = globaldir + \
            self.vars["stat_file_name_suffix"] + "." + str(myindex) + ".out"
        output.init_stat2(low_temp_stat_file,
                          self.output_objects,
                          extralabels=["rmf_file", "rmf_frame_index"])

        print "Setting up replica stat file"
        replica_stat_file = globaldir + \
            self.vars["replica_stat_file_suffix"] + "." + str(myindex) + ".out"
        output.init_stat2(replica_stat_file, [rex], extralabels=["score"])

        print "Setting up best pdb files"
        if not self.is_multi_state:
            output.init_pdb_best_scoring(pdb_dir + "/" +
                                         self.vars["best_pdb_name_suffix"],
                                         self.root_hier,
                                         self.vars[
                                             "number_of_best_scoring_models"],
                                         replica_exchange=True)
        else:
            for n in range(self.vars["number_of_states"]):
                output.init_pdb_best_scoring(pdb_dir + "/" + str(n) + "/" +
                                             self.vars["best_pdb_name_suffix"],
                                             self.root_hiers[n],
                                             self.vars[
                                                 "number_of_best_scoring_models"],
                                             replica_exchange=True)

# ---------------------------------------------

        if not self.em_object_for_rmf is None:
            if not self.is_multi_state:
                output_hierarchies = [
                    self.root_hier,
                    self.em_object_for_rmf.get_density_as_hierarchy(
                    )]
            else:
                output_hierarchies = self.root_hiers
                output_hierarchies.append(
                    self.em_object_for_rmf.get_density_as_hierarchy())
        else:
            if not self.is_multi_state:
                output_hierarchies = [self.root_hier]
            else:
                output_hierarchies = self.root_hiers

#----------------------------------------------
        print "Setting up and writing initial rmf coordinate file"
        init_suffix = globaldir + self.vars["initial_rmf_name_suffix"]
        output.init_rmf(init_suffix + "." + str(myindex) + ".rmf3",
                        output_hierarchies)
        if self.crosslink_restraints:
            output.add_restraints_to_rmf(
                init_suffix + "." + str(myindex) + ".rmf3",
                self.crosslink_restraints)
        output.write_rmf(init_suffix + "." + str(myindex) + ".rmf3")
        output.close_rmf(init_suffix + "." + str(myindex) + ".rmf3")

#----------------------------------------------

        print "Setting up production rmf files"

        rmfname = rmf_dir + "/" + str(myindex) + ".rmf3"
        output.init_rmf(rmfname, output_hierarchies)

        if self.crosslink_restraints:
            output.add_restraints_to_rmf(rmfname, self.crosslink_restraints)

        ntimes_at_low_temp = 0

        if myindex == 0:
            self.show_info()

        for i in range(self.vars["number_of_frames"]):
            for nr in range(self.vars["num_sample_rounds"]):
                if sampler_mc is not None:
                    sampler_mc.optimize(self.vars["monte_carlo_steps"])
                if sampler_md is not None:
                    sampler_md.optimize(self.vars["molecular_dynamics_steps"])
            score = self.model.evaluate(False)
            output.set_output_entry("score", score)

            my_temp_index = int(rex.get_my_temp() * temp_index_factor)

            if min_temp_index == my_temp_index:
                print "--- frame %s score %s " % (str(i), str(score))

                if i % self.vars["nframes_write_coordinates"]==0:
                    print '--- writing coordinates'
                    output.write_pdb_best_scoring(score)
                    output.write_rmf(rmfname)
                    output.set_output_entry("rmf_file", rmfname)
                    output.set_output_entry("rmf_frame_index", ntimes_at_low_temp)
                else:
                    output.set_output_entry("rmf_file", rmfname)
                    output.set_output_entry("rmf_frame_index", '-1')
                output.write_stat2(low_temp_stat_file)
                ntimes_at_low_temp += 1

            output.write_stat2(replica_stat_file)
            rex.swap_temp(i, score)


# -----------------------------------------------------------------------


def BuildModel0(
    m,
    data,
    resolutions=[1,
                 10],
    missing_bead_size=20,
        residue_per_gaussian=None):
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

    r = IMP.pmi.representation.Representation(m)

    # the dictionary for the hierarchies,
    hierarchies = {}

    for d in data:
                # retrieve the information from the data structure
        component_name = d[0]
        pdb_file = d[1]
        chain_id = d[2]
        color_id = d[3]
        fasta_file = d[4][0]
        # this function
        fastids = IMP.pmi.tools.get_ids_from_fasta_file(fasta_file)
        fasta_file_id = d[4][1]
        # avoid to add a component with the same name
        r.create_component(component_name,
                           color=color_id)

        r.add_component_sequence(component_name,
                                 fasta_file,
                                 id=fastids[fasta_file_id])

        hierarchies = r.autobuild_model(component_name,
                                        pdb_file,
                                        chain_id,
                                        resolutions=resolutions,
                                        missingbeadsize=missing_bead_size)

        r.show_component_table(component_name)

        r.set_rigid_bodies([component_name])

        r.set_chain_of_super_rigid_bodies(
            hierarchies,
            min_length=2,
            max_length=2)

        r.setup_component_sequence_connectivity(component_name, resolution=1)
        r.setup_component_geometry(component_name)

    r.setup_bonds()
    # put it at the end of rigid bodies
    r.set_floppy_bodies()

    # set current coordinates as reference for RMSD calculation
    r.set_current_coordinates_as_reference_for_rmsd("Reference")

    return r


# ----------------------------------------------------------------------

class BuildModel1(object):

    ''' this building scheme needs a data structure with the following fields
          comp_name
          hier_name
          color
          fasta_file
          fasta_id
          pdb_name
          chain_id
          res_range
          read_em_files
          bead_size
          rb
          super_rb
          em_num_components
          em_txt_file_name
          em_mrc_file_name
    '''

    def __init__(self, representation):
      self.simo=representation
      self.gmm_models_directory="."

    def set_gmm_models_directory(self,directory_name):
      self.gmm_models_directory=directory_name

    def build_model(self,data_structure,sequence_connectivity_scale=4.0):

      self.domain_dict={}
      self.resdensities={}
      super_rigid_bodies={}
      chain_super_rigid_bodies={}
      rigid_bodies={}

      for d in data_structure:
          comp_name         = d[0]
          hier_name         = d[1]
          color             = d[2]
          fasta_file        = d[3]
          fasta_id          = d[4]
          pdb_name          = d[5]
          chain_id          = d[6]
          res_range         = d[7][0:2]
          try:
             offset         = d[7][2]
          except:
             offset         = 0
          read_em_files     = d[8]
          bead_size         = d[9]
          rb                = d[10]
          super_rb          = d[11]
          em_num_components = d[12]
          em_txt_file_name  = d[13]
          em_mrc_file_name  = d[14]
          chain_of_super_rb = d[15]

          if comp_name not in self.simo.get_component_names():
             self.simo.create_component(comp_name,color=0.0)
             self.simo.add_component_sequence(comp_name,fasta_file,fasta_id)
          outhier=self.autobuild(self.simo,comp_name,pdb_name,chain_id,res_range,read=read_em_files,beadsize=bead_size,color=color,offset=offset)

          if not read_em_files is None:
             if em_txt_file_name is " ": em_txt_file_name=self.gmm_models_directory+"/"+hier_name+".txt"
             if em_mrc_file_name is " ": em_mrc_file_name=self.gmm_models_directory+"/"+hier_name+".mrc"


             dens_hier,beads=self.create_density(self.simo,comp_name,outhier,em_txt_file_name,em_mrc_file_name,em_num_components,read_em_files)
             self.simo.add_all_atom_densities(comp_name, hierarchies=beads)
             dens_hier+=beads

          else:
             dens_hier=[]

          self.resdensities[hier_name]=dens_hier
          self.domain_dict[hier_name]=outhier+dens_hier

          if rb is not None:
              if rb not in rigid_bodies:
                  rigid_bodies[rb]=[h for h in self.domain_dict[hier_name]]
              else:
                  rigid_bodies[rb]+=[h for h in self.domain_dict[hier_name]]


          if super_rb is not None:
              for k in super_rb:
                  if k not in super_rigid_bodies:
                     super_rigid_bodies[k]=[h for h in self.domain_dict[hier_name]]
                  else:
                     super_rigid_bodies[k]+=[h for h in self.domain_dict[hier_name]]

          if  chain_of_super_rb is not None:
              for k in chain_of_super_rb:
                  if k not in chain_super_rigid_bodies:
                     chain_super_rigid_bodies[k]=[h for h in self.domain_dict[hier_name]]
                  else:
                     chain_super_rigid_bodies[k]+=[h for h in self.domain_dict[hier_name]]

      self.rigid_bodies=rigid_bodies


      for c in self.simo.get_component_names():
          self.simo.setup_component_sequence_connectivity(c,scale=sequence_connectivity_scale)
          self.simo.setup_component_geometry(c)


      for rb in rigid_bodies:
          self.simo.set_rigid_body_from_hierarchies(rigid_bodies[rb])

      for k in super_rigid_bodies:
          self.simo.set_super_rigid_body_from_hierarchies(super_rigid_bodies[k])

      for k in chain_super_rigid_bodies:
          self.simo.set_chain_of_super_rigid_bodies(chain_super_rigid_bodies[k],2,3)

      self.simo.set_floppy_bodies()
      self.simo.setup_bonds()

    def get_density_hierarchies(self,hier_name_list):
        # return a list of density hierarchies
        # specify the list of hierarchy names
        dens_hier_list=[]
        for hn in hier_name_list:
            print hn
            dens_hier_list+=self.resdensities[hn]
        return dens_hier_list

    def get_pdb_bead_bits(self,hierarchy):
        pdbbits=[]
        beadbits=[]
        helixbits=[]
        for h in hierarchy:
           if "_pdb" in h.get_name():pdbbits.append(h)
           if "_bead" in h.get_name():beadbits.append(h)
           if "_helix" in h.get_name():helixbits.append(h)
        return (pdbbits,beadbits,helixbits)

    def scale_bead_radii(self,nresidues,scale):
        scaled_beads=set()
        for h in self.domain_dict:
          (pdbbits,beadbits,helixbits)=self.get_pdb_bead_bits(self.domain_dict[h])
          slope=(1.0-scale)/(1.0-float(nresidues))

          for b in beadbits:
              # I have to do the following
              # because otherwise we'll scale more than once
              if b not in scaled_beads:
                 scaled_beads.add(b)
              else:
                 continue
              radius=IMP.core.XYZR(b).get_radius()
              num_residues=len(IMP.pmi.tools.get_residue_indexes(b))
              scale_factor=slope*float(num_residues)+1.0
              print scale_factor
              new_radius=scale_factor*radius
              IMP.core.XYZR(b).set_radius(new_radius)
              print b.get_name()
              print "particle with radius "+str(radius)+" and "+str(num_residues)+" residues scaled to a new radius "+str(new_radius)


    def create_density(self,simo,compname,comphier,txtfilename,mrcfilename,num_components,read=True):
       #density generation for the EM restraint
       (pdbbits,beadbits,helixbits)=self.get_pdb_bead_bits(comphier)


       outhier=[]
       if read:
          if len(pdbbits)!=0:
            outhier+=simo.add_component_density(compname,
                                     pdbbits,
                                     num_components=num_components, # number of gaussian into which the simulated density is approximated
                                     resolution=0,      # resolution that you want to calculate the simulated density
                                     inputfile=txtfilename) # read what it was calculated before
          if len(helixbits)!=0:
            outhier+=simo.add_component_density(compname,
                                     helixbits,
                                     num_components=num_components, # number of gaussian into which the simulated density is approximated
                                     resolution=1,      # resolution that you want to calculate the simulated density
                                     inputfile=txtfilename) # read what it was calculated before


       else:
          if len(pdbbits)!=0:
            outhier+=simo.add_component_density(compname,
                                     pdbbits,
                                     num_components=num_components, # number of gaussian into which the simulated density is approximated
                                     resolution=0,      # resolution that you want to calculate the simulated density
                                     outputfile=txtfilename, # do the calculation
                                     outputmap=mrcfilename,
                                     multiply_by_total_mass=True) # do the calculation and output the mrc

          if len(helixbits)!=0:
            outhier+=simo.add_component_density(compname,
                                     helixbits,
                                     num_components=num_components, # number of gaussian into which the simulated density is approximated
                                     resolution=1,      # resolution that you want to calculate the simulated density
                                     outputfile=txtfilename, # do the calculation
                                     outputmap=mrcfilename,
                                     multiply_by_total_mass=True) # do the calculation and output the mrc

       return outhier,beadbits

    def autobuild(self,simo,comname,pdbname,chain,resrange,read=True,beadsize=5,color=0.0,offset=0):

        if pdbname is not None and pdbname is not "IDEAL_HELIX" and pdbname is not "BEADS" :
          if resrange[-1]==-1: resrange=(resrange[0],len(simo.sequence_dict[comname]))
          if read==False:
             outhier=simo.autobuild_model(comname,
                              pdbname=pdbname,
                              chain=chain,
                              resrange=resrange,
                              resolutions=[0,1,10],
                              offset=offset,
                              color=color,
                              missingbeadsize=beadsize)
          else:
             outhier=simo.autobuild_model(comname,
                              pdbname=pdbname,
                              chain=chain,
                              resrange=resrange,
                              resolutions=[1,10],
                              offset=offset,
                              color=color,
                              missingbeadsize=beadsize)


        elif pdbname is not None and pdbname is "IDEAL_HELIX" and pdbname is not "BEADS" :

          outhier=simo.add_component_ideal_helix(comname,
                                              resolutions=[1,10],
                                              resrange=resrange,
                                              color=color,
                                              show=False)

        elif pdbname is not None and pdbname is not "IDEAL_HELIX" and pdbname is "BEADS" :
          outhier=simo.add_component_necklace(comname,resrange[0],resrange[1],beadsize,color=color)

        else:

          seq_len=len(simo.sequence_dict[comname])
          outhier=simo.add_component_necklace(comname,
                                begin=1,
                                end=seq_len,
                                length=beadsize)

        return outhier


# ----------------------------------------------------------------------

class AnalysisReplicaExchange0(object):

    def __init__(self, model,
                 stat_file_name_suffix="stat",
                 # if you want to merge two calculation directories
                 merge_directories=["./"],
                 best_pdb_name_suffix="model",
                 do_clean_first=True,
                 do_create_directories=True,
                 global_output_directory="./",
                 rmf_dir="rmfs/",
                 best_pdb_dir="pdbs/",
                 replica_stat_file_suffix="stat_replica",
                 global_analysis_result_directory="./analysis/"):

        import glob
        self.model = model
        rmf_dir = global_output_directory + rmf_dir
        self.rmf_files = glob.glob(rmf_dir + "/*.rmf")
        stat_dir = global_output_directory
        self.stat_files = []
        # it contains the position of the root directories
        self.root_directory_dict = {}
        for rd in merge_directories:
            stat_files = glob.glob(rd + "/" + stat_dir + "/stat.*.out")
            self.stat_files += stat_files
            for s in stat_files:
                self.root_directory_dict[s] = rd





    def clustering(self, score_key,
                   rmf_file_key,
                   rmf_file_frame_key,
                   prefiltervalue=None,
                   # give a float value, create filtered files for
                   # which the score key is below a certain
                   # threshold (faster and less
                   # memory greedy)
                   feature_keys=None,
                   outputdir="./",
                   alignment_components=None,
                   number_of_best_scoring_models=10,
                   rmsd_calculation_components=None,
                   distance_matrix_file=None,
                   load_distance_matrix_file=False,
                   is_mpi=False,
                   skip_clustering=False,
                   # if True, will just extract the best scoring models
                   # and save the pdbs
                   number_of_clusters=1,
                   display_plot=False,
                   exit_after_display=True,
                   get_every=1,
                   first_and_last_frames=None,
                   # pass here a tuple with the first and the last
                   # frames to be analysed. If None (default)
                   # get all frames.
                   density_custom_ranges=None,
                   write_pdb_with_centered_coordinates=False,
                   voxel_size=5.0):
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
        import RMF
        import IMP.rmf

        if is_mpi:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            number_of_processes = comm.size
        else:
            rank = 0
            number_of_processes = 1

        if alignment_components is None:
            alignment_flag = 0
        else:
            alignment_flag = 1

        if feature_keys is None:
            feature_keys = []

        print "setup clustering class"
        Clusters = IMP.pmi.analysis.Clustering()

        if not load_distance_matrix_file:

            if len(self.stat_files)==0: print "ERROR: no stat file found in the given path"; return
            my_stat_files=IMP.pmi.tools.chunk_list_into_segments(self.stat_files,number_of_processes)[rank]


            score_list=[]
            rmf_file_list=[]
            rmf_file_frame_list=[]

            if not score_key is None:
                feature_keyword_list_dict = {}
                for sf in my_stat_files:

                    root_directory_of_stat_file = self.root_directory_dict[sf]

                    print "getting data from file %s" % sf
                    po = IMP.pmi.output.ProcessOutput(sf)
                    keywords = po.get_keys()

                    feature_keywords = [
                        score_key,
                        rmf_file_key,
                        rmf_file_frame_key]
                    for k in keywords:
                        for fk in feature_keys:
                            if fk in k:
                                feature_keywords.append(k)

                    if prefiltervalue is None:
                        fields = po.get_fields(
                            feature_keywords,
                            get_every=get_every)
                    else:
                        fields = po.get_fields(
                            feature_keywords,
                            filtertuple=(
                                score_key,
                                "<",
                                prefiltervalue),
                            get_every=get_every)
                    # check that all lengths are all equal
                    length_set = set()
                    for f in fields:
                        length_set.add(len(fields[f]))

                    # if some of the fields are missing, truncate
                    # the feature files to the shortest one

                    if len(length_set) > 1:
                        print "AnalysisReplicaExchange0.clustering: the statfile is not synchronous"
                        minlen = min(length_set)
                        for f in fields:
                            fields[f] = fields[f][0:minlen]

                    # append to the lists
                    score_list += fields[score_key]
                    rmffilelist = []
                    for rmf in fields[rmf_file_key]:
                        rmffilelist.append(
                            root_directory_of_stat_file +
                            "/" +
                            rmf)

                    rmf_file_list += rmffilelist

                    rmf_file_frame_list += fields[rmf_file_frame_key]

                    for k in feature_keywords:
                        if k in feature_keyword_list_dict:
                            feature_keyword_list_dict[k] += fields[k]
                        else:
                            feature_keyword_list_dict[k] = fields[k]

# ------------------------------------------------------------------------
# broadcast the conformation features
            if number_of_processes > 1:
                score_list = IMP.pmi.tools.scatter_and_gather(score_list)
                rmf_file_list = IMP.pmi.tools.scatter_and_gather(rmf_file_list)
                rmf_file_frame_list = IMP.pmi.tools.scatter_and_gather(
                    rmf_file_frame_list)
                for k in feature_keyword_list_dict:
                    feature_keyword_list_dict[k] = IMP.pmi.tools.scatter_and_gather(
                        feature_keyword_list_dict[k])

            # sort by score and ge the best scoring ones
            score_rmf_tuples = zip(
                score_list,
                rmf_file_list,
                rmf_file_frame_list,
                range(len(score_list)))

            # first slice the ensemble accoridng to the user requirements
            if not first_and_last_frames is None:
                nframes = len(score_rmf_tuples)
                first_frame = int(first_and_last_frames[0] * nframes)
                last_frame = int(first_and_last_frames[1] * nframes)
                if last_frame > len(score_rmf_tuples):
                    last_frame = -1
                score_rmf_tuples = score_rmf_tuples[first_frame:last_frame]

            # numerically sorting in ascending order
            sorted_score_rmf_tuples = sorted(
                score_rmf_tuples,
                key=lambda x: float(x[0]))
            best_score_rmf_tuples = sorted_score_rmf_tuples[
                0:number_of_best_scoring_models]
            # we have to define a best score rank when we process in parallel
            # to get the correct position in the best_score_rmf_tuples
            for n, tpl in enumerate(best_score_rmf_tuples):
                tpl_new = tpl + tuple([n])
                best_score_rmf_tuples[n] = tpl_new

# prune the feature_keyword_list_dict, keep only the best scoring models
# and sort it

            best_score_feature_keyword_list_dict = {}
            for tpl in best_score_rmf_tuples:

                index = tpl[3]
                for f in feature_keyword_list_dict:
                    if f in best_score_feature_keyword_list_dict:
                        best_score_feature_keyword_list_dict[f].append(
                            feature_keyword_list_dict[f][index])
                    else:
                        best_score_feature_keyword_list_dict[f] = [
                            feature_keyword_list_dict[f][index]]

            del feature_keyword_list_dict

# ------------------------------------------------------------------------

            my_best_score_rmf_tuples = IMP.pmi.tools.chunk_list_into_segments(
                best_score_rmf_tuples,
                number_of_processes)[rank]

# ------------------------------------------------------------------------

            # extract all information without clustering

            if skip_clustering:

               dircluster=outputdir+"/all_models."+str(n)+"/"
               try:
                   os.mkdir(outputdir)
               except:
                   pass

               try:
                   os.mkdir(dircluster)
               except:
                   pass

               clusstat=open(dircluster+"stat."+str(rank)+".out","w")

               for cnt,tpl in enumerate(my_best_score_rmf_tuples):
                    rmf_name=tpl[1]
                    rmf_frame_number=tpl[2]

                    tmp_dict={}
                    index=tpl[4]

                    for key in best_score_feature_keyword_list_dict:
                        tmp_dict[key]=best_score_feature_keyword_list_dict[key][index]

                    prot=IMP.pmi.analysis.get_hier_from_rmf(self.model,rmf_frame_number,rmf_name)

                    if not prot: continue

                    o=IMP.pmi.output.Output()
                    o.init_pdb(dircluster+str(cnt)+"."+str(rank)+".pdb",prot)
                    o.write_pdb(dircluster+str(cnt)+"."+str(rank)+".pdb",
                    translate_to_geometric_center=write_pdb_with_centered_coordinates)


                    tmp_dict["local_pdb_file_name"]=str(cnt)+"."+str(rank)+".pdb"
                    tmp_dict["rmf_file_full_path"]=rmf_name
                    tmp_dict["local_rmf_file_name"]=str(cnt)+"."+str(rank)+".rmf3"
                    tmp_dict["local_rmf_frame_number"]=0

                    #IMP.atom.destroy(prot)

                    clusstat.write(str(tmp_dict)+"\n")

                    o.init_rmf(dircluster+str(cnt)+"."+str(rank)+".rmf3",[prot])
                    #IMP.rmf.add_restraints(o.dictionary_rmfs[dircluster+str(n)+".rmf3"],restraints)
                    o.write_rmf(dircluster+str(cnt)+"."+str(rank)+".rmf3")
                    o.close_rmf(dircluster+str(cnt)+"."+str(rank)+".rmf3")

               return


            # here I've tested that feature_keyword_list_dict is correct on 2 CPUs

# --------------------------------------------------------------------------------------------
# read the coordinates
# --------------------------------------------------------------------------------------------

            '''
            dircluster = outputdir + "/all_models." + str(n) + "/"
            try:
                os.mkdir(outputdir)
            except:
                pass

            try:
                os.mkdir(dircluster)
            except:
                pass

            clusstat = open(dircluster + "stat." + str(rank) + ".out", "w")

            for cnt, tpl in enumerate(my_best_score_rmf_tuples):
                rmf_name = tpl[1]
                rmf_frame_number = tpl[2]

                tmp_dict = {}
                index = tpl[4]

                for key in best_score_feature_keyword_list_dict:
                    tmp_dict[
                        key] = best_score_feature_keyword_list_dict[
                        key][
                        index]

                prot = IMP.pmi.analysis.get_hier_from_rmf(
                    self.model,
                    rmf_frame_number,
                    rmf_name)

                if not prot:
                    continue

                o = IMP.pmi.output.Output()
                o.init_pdb(
                    dircluster + str(cnt) + "." + str(rank) + ".pdb",
                    prot)
                o.write_pdb(
                    dircluster + str(cnt) + "." + str(rank) + ".pdb",
                    translate_to_geometric_center=write_pdb_with_centered_coordinates)

                tmp_dict["pdb_file_name"] = str(
                    cnt) + "." + str(rank) + ".pdb"

                # IMP.atom.destroy(prot)

                clusstat.write(str(tmp_dict) + "\n")

                o.init_rmf(
                    dircluster + str(cnt) + "." + str(rank) + ".rmf3",
                    [prot])
                # IMP.rmf.add_restraints(o.dictionary_rmfs[dircluster+str(n)+".rmf3"],restraints)
                o.write_rmf(
                    dircluster + str(cnt) + "." + str(rank) + ".rmf3")
                o.close_rmf(
                    dircluster + str(cnt) + "." + str(rank) + ".rmf3")
             '''


            # here I've tested that feature_keyword_list_dict is correct on 2
            # CPUs
# ------------------------------------------------------------------------
# read the coordinates
# ------------------------------------------------------------------------
            # reading the coordinates
            all_coordinates = []
            rmsd_coordinates = []
            alignment_coordinates = []
            all_rmf_file_names = []
            # it will be used to extract the features
            rmf_file_name_index_dict = {}
            for cnt, tpl in enumerate(my_best_score_rmf_tuples):
                rmf_file = tpl[1]
                frame_number = tpl[2]

                prot = IMP.pmi.analysis.get_hier_from_rmf(
                    self.model,
                    frame_number,
                    rmf_file)
                if not prot:
                    continue
                # getting the particles
                part_dict = IMP.pmi.analysis.get_particles_at_resolution_one(
                    prot)

                all_particles=[pp for key in part_dict for pp in part_dict[key]]


                # getting the coordinates
                model_coordinate_dict = {}
                template_coordinate_dict={}
                rmsd_coordinate_dict={}
                for pr in part_dict:
                    coords = np.array(
                        [np.array(IMP.core.XYZ(i).get_coordinates()) for i in part_dict[pr]])
                    model_coordinate_dict[pr] = coords

                if alignment_components is not None:
                  for pr in alignment_components:
                    if type(alignment_components[pr]) is str:
                       name=alignment_components[pr]
                       s=IMP.atom.Selection(prot,molecule=name)
                    elif type(alignment_components[pr]) is tuple:
                       name=alignment_components[pr][2]
                       rend=alignment_components[pr][1]
                       rbegin=alignment_components[pr][0]
                       s=IMP.atom.Selection(prot,molecule=name,residue_indexes=range(rbegin,rend+1))
                    ps=s.get_selected_particles()
                    filtered_particles=list(set(ps)&set(all_particles))
                    coords = np.array(
                        [np.array(IMP.core.XYZ(i).get_coordinates()) for i in filtered_particles])
                    template_coordinate_dict[pr] = coords

                if rmsd_calculation_components is not None:
                  for pr in rmsd_calculation_components:
                    if type(rmsd_calculation_components[pr]) is str:
                       name=rmsd_calculation_components[pr]
                       s=IMP.atom.Selection(prot,molecule=name)
                    elif type(rmsd_calculation_components[pr]) is tuple:
                       name=rmsd_calculation_components[pr][2]
                       rend=rmsd_calculation_components[pr][1]
                       rbegin=rmsd_calculation_components[pr][0]
                       s=IMP.atom.Selection(prot,molecule=name,residue_indexes=range(rbegin,rend+1))
                    ps=s.get_selected_particles()
                    filtered_particles=list(set(ps)&set(all_particles))
                    coords = np.array(
                        [np.array(IMP.core.XYZ(i).get_coordinates()) for i in filtered_particles])
                    rmsd_coordinate_dict[pr] = coords


                all_coordinates.append(model_coordinate_dict)
                alignment_coordinates.append(template_coordinate_dict)
                rmsd_coordinates.append(rmsd_coordinate_dict)
                frame_name = rmf_file + '|' + str(frame_number)
                all_rmf_file_names.append(frame_name)
                rmf_file_name_index_dict[frame_name] = tpl[4]

# ------------------------------------------------------------------------
# broadcast the coordinates
            if number_of_processes > 1:
                all_coordinates = IMP.pmi.tools.scatter_and_gather(
                    all_coordinates)
                all_rmf_file_names = IMP.pmi.tools.scatter_and_gather(
                    all_rmf_file_names)
                rmf_file_name_index_dict = IMP.pmi.tools.scatter_and_gather(
                    rmf_file_name_index_dict)
                alignment_coordinates=IMP.pmi.tools.scatter_and_gather(
                    alignment_coordinates)
                rmsd_coordinates=IMP.pmi.tools.scatter_and_gather(
                    rmsd_coordinates)

            if rank == 0:
                # save needed informations in external files
                self.save_objects(
                    [best_score_feature_keyword_list_dict,
                     rmf_file_name_index_dict],
                    ".macro.pkl")

# ------------------------------------------------------------------------
            for n, model_coordinate_dict in enumerate(all_coordinates):

                template_coordinate_dict = {}
                # let's try to align
                if alignment_flag == 1 and len(Clusters.all_coords) == 0:
                    #for pr in alignment_components:
                    #    template_coordinate_dict[pr] = model_coordinate_dict[pr]
                # set the first model as template coordinates
                    Clusters.set_template(alignment_coordinates[n])

                # set particles to calculate RMSDs on
                # (if global, list all proteins, or just a few for local)
                # setting all the coordinates for rmsd calculation
                #rmsd_coordinate_dict = {}
                #if rmsd_calculation_components is None:
                #    rmsd_calculation_components = alignment_components

                #for pr in rmsd_calculation_components:
                #    rmsd_coordinate_dict[pr] = model_coordinate_dict[pr]

                Clusters.fill(all_rmf_file_names[n], rmsd_coordinates[n])

            print "Global calculating the distance matrix"

            # calculate distance matrix, all against all

            Clusters.dist_matrix(is_mpi=is_mpi)

            Clusters.do_cluster(number_of_clusters)
            if display_plot:
                if rank == 0:
                    Clusters.plot_matrix()
                if number_of_processes > 1:
                    comm.Barrier()
                if exit_after_display:
                    exit()
            Clusters.save_distance_matrix_file(file_name=distance_matrix_file)

# ------------------------------------------------------------------------
# load the distance matrix from file

        else:
            Clusters.load_distance_matrix_file(file_name=distance_matrix_file)
            print "clustering with %s clusters" % str(number_of_clusters)
            Clusters.do_cluster(number_of_clusters)
            [best_score_feature_keyword_list_dict,
                rmf_file_name_index_dict] = self.load_objects(".macro.pkl")
            if display_plot:
                if rank == 0:
                    Clusters.plot_matrix()
                if number_of_processes > 1:
                    comm.Barrier()
                if exit_after_display:
                    exit()

# ------------------------------------------------------------------------
# now save all informations about the clusters

# ------------------------------------------------------------------------
#

        if rank == 0:

            print Clusters.get_cluster_labels()

            for n, cl in enumerate(Clusters.get_cluster_labels()):

                print "rank %s " % str(rank)
                print "cluster %s " % str(n)
                print "cluster label %s " % str(cl)
                print Clusters.get_cluster_label_names(cl)

                # first initialize the Density class if requested

                if density_custom_ranges:

                    DensModule = IMP.pmi.analysis.GetModelDensity(
                        density_custom_ranges,
                        voxel=voxel_size)

                dircluster = outputdir + "/cluster." + str(n) + "/"
                try:
                    os.mkdir(outputdir)
                except:
                    pass

                try:
                    os.mkdir(dircluster)
                except:
                    pass

                rmsd_dict = {"AVERAGE_RMSD":
                             str(Clusters.get_cluster_label_average_rmsd(cl))}

                clusstat = open(dircluster + "stat.out", "w")

                for k, structure_name in enumerate(Clusters.get_cluster_label_names(cl)):

                        # extract the features
                    tmp_dict = {}
                    tmp_dict.update(rmsd_dict)
                    index = rmf_file_name_index_dict[structure_name]
                    for key in best_score_feature_keyword_list_dict:
                        tmp_dict[
                            key] = best_score_feature_keyword_list_dict[
                            key][
                            index]

                    # get the rmf name and the frame number from the list of
                    # frame names
                    rmf_name = structure_name.split("|")[0]
                    rmf_frame_number = int(structure_name.split("|")[1])

                    clusstat.write(str(tmp_dict) + "\n")

                    prot = IMP.pmi.analysis.get_hier_from_rmf(
                        self.model,
                        rmf_frame_number,
                        rmf_name)

                    if not prot:
                        continue

                    #rh= RMF.open_rmf_file_read_only(rmf_name)
                    # restraints=IMP.rmf.create_restraints(rh,self.model)
                    #del rh


                    if k > 0:
                        model_index = Clusters.get_model_index_from_name(
                            structure_name)
                        transformation = Clusters.get_transformation_to_first_member(
                            cl,
                            model_index)

                        rbs = set()
                        for p in IMP.atom.get_leaves(prot):
                            if not IMP.core.XYZR.get_is_setup(p):
                                print 'does this ever happen?'
                                IMP.core.XYZR.setup_particle(p)
                                IMP.core.XYZR(p).set_radius(0.0001)
                                IMP.core.XYZR(p).set_coordinates((0, 0, 0))

                            if IMP.core.RigidBodyMember.get_is_setup(p):
                                rb = IMP.core.RigidBodyMember(p).get_rigid_body()
                                rbs.add(rb)
                            else:
                                IMP.core.transform(IMP.core.XYZ(p),
                                                   transformation)
                        for rb in rbs:
                            IMP.core.transform(rb,transformation)

                    # add the density
                    if density_custom_ranges:
                        DensModule.add_subunits_density(prot)

                    o = IMP.pmi.output.Output()
                    o.init_pdb(dircluster + str(k) + ".pdb", prot)
                    o.write_pdb(dircluster + str(k) + ".pdb")

                    o.init_rmf(dircluster + str(k) + ".rmf3", [prot])
                    # IMP.rmf.add_restraints(o.dictionary_rmfs[dircluster+str(n)+".rmf3"],restraints)
                    o.write_rmf(dircluster + str(k) + ".rmf3")
                    o.close_rmf(dircluster + str(k) + ".rmf3")

                    del o
                    # IMP.atom.destroy(prot)

                if density_custom_ranges:
                    DensModule.write_mrc(path=dircluster)
                    del DensModule

        if is_mpi:
            comm.Barrier()

    def save_objects(self, objects, file_name):
        import pickle
        outf = open(file_name, 'w')
        pickle.dump(objects, outf)
        outf.close()

    def load_objects(self, file_name):
        import pickle
        inputf = open(file_name, 'r')
        objects = pickle.load(inputf)
        inputf.close()
        return objects
