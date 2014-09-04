## \example pmi/crosslinking/atomic_xl.py

import IMP
import IMP.atom
import IMP.pmi
import IMP.pmi.representation_new
import IMP.pmi.sequence_tools
import IMP.pmi.restraints_new.stereochemistry
import IMP.pmi.restraints_new.atomic_xl
import IMP.pmi.hierarchy_tools as hierarchy_tools
import IMP.pmi.sampling_tools as sampling_tools
import IMP.pmi.data_tools as data_tools
import IMP.pmi.macros
import IMP.isd_emxl

### parameters ###

# scoring params
xl_length = 5.0
slope = 0.01
charmm_weight = 0.1
psi_init = 0.0
sigma_init = 2.0

# sampling params
num_md = 40
num_mc = 0
num_rounds = 1
out_dir = "."
nframes = 10000
nmodels = 10
min_temp = 1.0
max_temp = 1.0
mc_max_step = 0.1

### create system ###
mdl = IMP.Model()
system = IMP.pmi.representation_new.System(mdl)
seqs = IMP.pmi.representation_new.Sequences(fasta_fn="data/example.fasta",
                                            name_map={'example_protein':'prot'})
state = system.create_state()
prot = state.create_molecule("prot", sequence=seqs["prot"], chain_id='A')
atomic = prot.add_structure(pdb_fn="data/example_protein.pdb",chain_id='A',offset=1)
prot.add_representation(atomic,'balls',[0])

state2 = system.create_state()
prot2 = state2.create_molecule("prot2", sequence=seqs["prot"], chain_id='A')
atomic2 = prot2.add_structure(pdb_fn="data/example_protein.pdb",chain_id='A',offset=1)
prot2.add_representation(atomic2,'balls',[0])

hier = system.build()

### create restraints ###
output_objects=[]
# charmm 1
charmm1 = IMP.pmi.restraints_new.stereochemistry.CharmmForceFieldRestraint(state.get_hierarchy(),
                                                                           enable_nonbonded=True)
charmm1.set_weight(charmm_weight)
charmm1.add_to_model()
print 'EVAL - CHARMM 1',mdl.evaluate(False)
output_objects.append(charmm1)
charmm2 = IMP.pmi.restraints_new.stereochemistry.CharmmForceFieldRestraint(state2.get_hierarchy(),
                                                                           enable_nonbonded=True)
charmm2.add_to_model()
charmm2.set_weight(charmm_weight)
print 'EVAL - CHARMM 2',mdl.evaluate(False)
output_objects.append(charmm2)

# crosslinks
xls={}
xls[0]=[{'r1': {'residue_index':88},
         'r2': {'residue_index':146},
         'Score': 1.0}]
xls[1]=[{'r1': {'residue_index':141},
         'r2': {'residue_index':146},
         'Score': 1.0}]
xlrs = IMP.pmi.restraints_new.atomic_xl.AtomicCrossLinkMSRestraint(hier,
                                                                   xls,
                                                                   length=xl_length,
                                                                   slope=slope,
                                                                   nstates=2,
                                                                   psi_init=psi_init,
                                                                   sigma_init=sigma_init)
xlrs.add_to_model()
output_objects.append(xlrs)
print 'EVAL - XLINKS',mdl.evaluate(False)

### sample ###

### OUTPUT HACK ###
class mini_output:
    def __init__(self,mdl):
        self.mdl=mdl
    def get_output(self):
        output={}
        output["SimplifiedModel_Total_Score"] = str(self.mdl.evaluate(False))
        return output
output_objects.append(mini_output(mdl))
##################

#### SAMPLING HACK ####
for p in IMP.core.get_leaves(hier):
    IMP.core.XYZ(p).set_coordinates_are_optimized(False)

md_objs=sampling_tools.enable_md_sampling(mdl,
                                          particles=xlrs.get_particles(),
                                          include_siblings=True,
                                          exclude_backbone=True)

#################


rex = IMP.pmi.macros.ReplicaExchange0(mdl,
                            root_hier = hier,
                            molecular_dynamics_sample_objects=md_objs,
                            output_objects = output_objects,
                            replica_exchange_minimum_temperature=min_temp,
                            replica_exchange_maximum_temperature=max_temp,
                            num_sample_rounds=num_rounds,
                            molecular_dynamics_steps=num_md,
                            monte_carlo_steps=num_mc,
                            number_of_best_scoring_models=nmodels,
                            number_of_frames = nframes,
                            write_initial_rmf=True,
                            global_output_directory=out_dir,
                            atomistic=True,
                            crosslink_restraints=xlrs.create_restraints_for_rmf())
rex.execute_macro()
