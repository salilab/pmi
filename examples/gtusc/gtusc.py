## \example pmi/sandbox_example.py

import IMP
import IMP.rmf
import IMP.pmi
import IMP.pmi.representation_new
import IMP.pmi.restraints_new.stereochemistry
import IMP.pmi.restraints_new.atomic_xl
import IMP.pmi.sequence_tools
import IMP.pmi.hierarchy_tools
import IMP.pmi.data_tools
import IMP.pmi.samplers
import IMP.pmi.macros
import RMF


# files and settings
closed_fn = 'data/closed.pdb'
sse_fn='data/closed.dssp'
open_fn = 'data/open.pdb'
seq_fn = 'data/sequences.fasta'
xl_fn='data/xlinks.txt'
out_dir='out/'

# restraint options
slope=0.05
elastic_strength=100.0
elastic_cutoff=7.0
symmetry_strength=10.0

#sampling options
num_md = 20
num_mc = 1
num_rounds = 1
min_temp=1.0
max_temp=1.2
nmodels=100
nframes=10000
charmm_weight=1.0
mc_max_step = 0.1

### setup 2-state system
mdl = IMP.Model()
system = IMP.pmi.representation_new.System(mdl)

### read sequences into a little data structure
seqs = IMP.pmi.representation_new.Sequences(fasta_fn=seq_fn,
                                            name_map={'GCP2_YEAST':'Spc97',
                                                      'GCP3_YEAST':'Spc98',
                                                      'ytubulin':'ytubulin'})

# setup closed
print 'setting up state 1'
state1 = system.create_state()
spc97 = state1.create_molecule("Spc97", sequence=seqs["Spc97"], chain_id='A')
spc98 = state1.create_molecule("Spc98", sequence=seqs["Spc98"], chain_id='B')
ytub = state1.create_molecule("ytub", sequence=seqs["ytubulin"], chain_id='C')
s97_atomic = spc97.add_structure(pdb_fn=closed_fn,chain_id='A')
s98_atomic = spc98.add_structure(pdb_fn=closed_fn,chain_id='B')
ytub_atomic = ytub.add_structure(pdb_fn=closed_fn,chain_id='C')
ytub.add_copy(pdb_fn=closed_fn,chain_id='D')
spc97.add_copy(pdb_fn=closed_fn,chain_id='G')
spc98.add_copy(pdb_fn=closed_fn,chain_id='H')
ytub.add_copy(pdb_fn=closed_fn,chain_id='I')
ytub.add_copy(pdb_fn=closed_fn,chain_id='J')
spc97.add_representation(s97_atomic,'balls',[0])
spc98.add_representation(s98_atomic,'balls',[0])
ytub.add_representation(ytub_atomic,'balls',[0])


print 'setting up state 2'
state2 = system.create_state()
spc97_2 = state2.create_molecule("Spc97", sequence=seqs["Spc97"], chain_id='A')
spc98_2 = state2.create_molecule("Spc98", sequence=seqs["Spc98"], chain_id='B')
ytub_2 = state2.create_molecule("ytub", sequence=seqs["ytubulin"], chain_id='C')
s97_atomic_2 = spc97_2.add_structure(pdb_fn=open_fn,chain_id='A')
s98_atomic_2 = spc98_2.add_structure(pdb_fn=open_fn,chain_id='B')
ytub_atomic_2 = ytub_2.add_structure(pdb_fn=open_fn,chain_id='C')
ytub_2.add_copy(pdb_fn=open_fn,chain_id='D')
spc97_2.add_copy(pdb_fn=open_fn,chain_id='G')
spc98_2.add_copy(pdb_fn=open_fn,chain_id='H')
ytub_2.add_copy(pdb_fn=open_fn,chain_id='I')
ytub_2.add_copy(pdb_fn=open_fn,chain_id='J')
spc97_2.add_representation(s97_atomic_2,'balls',[0])
spc98_2.add_representation(s98_atomic_2,'balls',[0])
ytub_2.add_representation(ytub_atomic_2,'balls',[0])

### build system
print 'building system'
hier = system.build()

#### MD SAMPLING HACK ####
vxkey = IMP.FloatKey('vx')
vykey = IMP.FloatKey('vy')
vzkey = IMP.FloatKey('vz')
for pp in IMP.core.get_leaves(hier):
    IMP.core.XYZ(pp).set_coordinates_are_optimized(True)
    pp.add_attribute(vxkey, 0.0)
    pp.add_attribute(vykey, 0.0)
    pp.add_attribute(vzkey, 0.0)
md_objects = [IMP.pmi.restraints_new.atomic_xl.SampleObjects(
    'Floppy_Bodies_SimplifiedModel',[IMP.core.get_leaves(hier)])]
#############


### add restraints
output_objects=[]

### symmetry
ctrans = IMP.algebra.Transformation3D(IMP.algebra.get_rotation_from_matrix(0.568131,-0.822938,0.0,
                                                                           0.822938,0.56813,0.0,
                                                                           0.0,0.0,1.0),
                                      IMP.algebra.Vector3D(0,0,-18.8))
otrans = IMP.algebra.Transformation3D(IMP.algebra.get_rotation_from_matrix(0.568131,-0.822938,0.0,
                                                                           0.822938,0.56813,0.0,
                                                                           0.0,0.0,1.0),
                                      IMP.algebra.Vector3D(0,0,-22.2))

sym_dat=[['Spc97',0,'Spc97',1,'Spc97'],
         ['Spc98',0,'Spc98',1,'Spc98'],
         ['ytub',0,'ytub',2,'ytub97'],
         ['ytub',1,'ytub',3,'ytub98']]
for nstate in range(2):
    if nstate==0:
        trans = ctrans
    else:
        trans = otrans
    for mol1,copy1,mol2,copy2,name in sym_dat:
        name+='_state'+str(nstate)
        sel={"atom_type":IMP.atom.AtomType("CA"),"state_index":nstate}
        r=IMP.pmi.restraints_new.stereochemistry.SymmetryRestraint(hier,
                                                                   [trans],
                                                                   [{'molecule':mol1,'copy_index':copy1},
                                                                    {'molecule':mol2,'copy_index':copy2}],
                                                                   extra_sel=sel,
                                                                   label=name,
                                                                   strength=symmetry_strength)
        r.add_to_model()
        output_objects.append(r)
        print 'EVAL - symmetry',name,mdl.evaluate(False)


### elastic network on ALL SSEs together!
sse_selections=IMP.pmi.data_tools.parse_dssp(sse_fn,'ABCDGHIJ')
all_sses=sse_selections['helix']+sse_selections['beta']
ers=[]
for i in range(2):
    print 'creating SSEs',i
    er=IMP.pmi.restraints_new.stereochemistry.ElasticNetworkRestraint(hier,
                                selection_dicts=(sse for sselist in all_sses for sse in sselist),
                                extra_sel={'atom_type':IMP.atom.AtomType('CA'),'state_index':i},
                                label='sses_'+str(i),
                                add_info_to_label=False,
                                strength=elastic_strength,
                                dist_cutoff=elastic_cutoff)
    er.set_weight(1.0)
    er.add_to_model()
    output_objects.append(er)
    ers.append(er)

### xlinks
xls = IMP.pmi.data_tools.parse_xlinks_davis(xl_fn,
                                            max_num=-1,
                                            name_map={'His-TEV-Tub4':'ytub'},
                                            named_offsets={'ytub':-33})
print 'setting up xl'
xlrs = IMP.pmi.restraints_new.atomic_xl.AtomicCrossLinkMSRestraint(hier,
                                                                   xls,
                                                                   length=11.4,
                                                                   slope=slope,
                                                                   nstates=2)
xlrs.add_to_model()
output_objects.append(xlrs)
print 'EVAL - XLINKS',mdl.evaluate(False)

# CHARMM
print 'adding charmm 1'
charmm1 = IMP.pmi.restraints_new.stereochemistry.CharmmForceFieldRestraint(state1.get_hierarchy(),
                                                                           enable_nonbonded=True)
charmm1.set_label('_1')
charmm1.set_weight(charmm_weight)
#charmm1.add_to_model()
#print 'EVAL - CHARMM 1'
#mdl.evaluate(False)
print '\n'
output_objects.append(charmm1)

print 'adding charmm 2'
charmm2 = IMP.pmi.restraints_new.stereochemistry.CharmmForceFieldRestraint(state2.get_hierarchy(),
                                                                           enable_nonbonded=True)
charmm2.set_label('_2')
charmm2.set_weight(charmm_weight)
charmm2.add_to_model()
print 'EVAL - CHARMM 2'
mdl.evaluate(False)

output_objects.append(charmm2)

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

rex = IMP.pmi.macros.ReplicaExchange0(mdl,
                            root_hier = hier,
                            molecular_dynamics_sample_objects=md_objects,
                            monte_carlo_sample_objects=xlrs.get_mc_sample_objects(mc_max_step),
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
                            crosslink_restraints=ers,
                            atomistic=True)
rex.execute_macro()
