## \example pmi/sandbox_example.py

import IMP
import IMP.rmf
import IMP.pmi
import IMP.pmi.representation_new
import IMP.pmi.restraints_new.stereochemistry
import IMP.pmi.restraints_new.em
import IMP.pmi.sequence_tools
import IMP.pmi.hierarchy_tools
from IMP.pmi.io.data_storage import SelectionDict,SubsequenceData,CrossLinkData
import IMP.pmi.io.data_parsers as data_parsers
import IMP.pmi.samplers
import IMP.pmi.sampling_tools as sampling_tools
import IMP.pmi.macros
import RMF


# files and settings

#native
#init_fn = 'data/1oel_A_fit4.pdb'
#seq_fn = 'data/1oel_pdb.fasta'
#out_dir='out2/'
#pdb_offset=-1

# model
init_fn = 'data/1oel_A_from_1we3.pdb'
seq_fn = 'data/1oel.fasta'
out_dir='out/'
pdb_offset=0

gmm_fn = 'gmms/1oel_A_200.txt'
sse_fn = 'data/1oel_A_from_1we3.dssp'
map_fn = 'data/1oel_A_4.mrc'


# restraint options
slope=0.0
elastic_strength=100.0
elastic_cutoff=7.0
em_weight=1.0
charmm_weight=1.0

#sampling options
num_md = 400
num_mc = 0
num_rounds = 1
min_temp=1.0
max_temp=1.2
nmodels=100
nframes=10000

### setup 1-state system
mdl = IMP.Model()
system = IMP.pmi.representation_new.System(mdl)

### read sequences into a little data structure
seqs = IMP.pmi.representation_new.Sequences(fasta_fn=seq_fn,
                                            name_map={'1oel_A':'1oel'})

# setup closed
print 'setting up state 1'
state = system.create_state()
groel = state.create_molecule("1oel", sequence=seqs["1oel"], chain_id='A')
atomic = groel.add_structure(pdb_fn=init_fn,chain_id='A',offset=pdb_offset)
groel.add_representation(atomic,'balls',[0])

### build system
print 'building system'
hier = system.build()
print 'system is built!'

### add restraints
output_objects=[]


### CHARMM
print 'adding charmm 1'
charmm = IMP.pmi.restraints_new.stereochemistry.CharmmForceFieldRestraint(hier)
charmm.set_weight(charmm_weight)
charmm.add_to_model()
print 'EVAL - CHARMM',mdl.evaluate(False)
output_objects.append(charmm)


### EM
#gaussian hack, should be done in representation...
mass=0.0
for p in IMP.atom.get_leaves(hier):
    center=IMP.core.XYZ(p).get_coordinates()
    rad=IMP.core.XYZR(p).get_radius()
    mass += IMP.atom.Mass(p).get_mass()
    trans=IMP.algebra.Transformation3D(IMP.algebra.get_identity_rotation_3d(),
                                       center)
    shape=IMP.algebra.Gaussian3D(IMP.algebra.ReferenceFrame3D(trans),[rad]*3)
    IMP.core.Gaussian.setup_particle(p,shape)
gem = IMP.pmi.restraints_new.em.GaussianEMRestraint(hier=hier,
                                        target_fn=gmm_fn,
                                        target_radii_scale=3.0,
                                        target_mass_scale=mass,
                                        spherical_gaussians=True,
                                        pointwise_restraint=True,
                                        slope=slope,
                                        mm_container=charmm.get_close_pair_container(),
                                        use_log_score=False,
                                        orig_map_fn=map_fn)
gem.set_weight(em_weight)
gem.add_to_model()
output_objects.append(gem)
print 'EVAL - EM',mdl.evaluate(False)

### elastic network for SSEs
sses = data_parsers.parse_dssp(mdl,sse_fn)
ers=[]
for ns,sse in enumerate(sses['helix']+sses['beta']):
    er=IMP.pmi.restraints_new.stereochemistry.ElasticNetworkRestraint(hier,
                        selection_dicts=sse,
                        label='sse',
                        add_info_to_label=True,
                        strength=elastic_strength,
                        dist_cutoff=elastic_cutoff,
                        atom_type=IMP.atom.AtomType("CA"))
    er.set_weight(1.0)
    er.add_to_model()
    output_objects.append(er)
    ers.append(er)


### OUTPUT HACK ###
md_objects = [sampling_tools.SampleObjects(
    'Floppy_Bodies_SimplifiedModel',[IMP.core.get_leaves(hier)])]
IMP.pmi.sampling_tools.enable_md_sampling(mdl,hier)
class mini_output:
    def __init__(self,mdl):
        self.mdl=mdl
    def get_output(self):
        output={}
        output["SimplifiedModel_Total_Score"] = str(self.mdl.evaluate(False))
        return output
output_objects.append(mini_output(mdl))
##################

rex = IMP.pmi.macros.ReplicaExchange0(mdl,
                            root_hier = hier,
                            molecular_dynamics_sample_objects=md_objects,
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
                            atomistic=True)
rex.execute_macro()
