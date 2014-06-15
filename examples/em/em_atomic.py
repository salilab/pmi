## \example pmi/em/em_atomic.py

import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.em
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros
import sys
import os
import math
############ SETTINGS ###############
# Data options
m=IMP.Model()
#in_fn='data/helix3.pdb'
in_fn='data/1oel_A_fit4.pdb'
sse_fn='data/helix3.dssp'
model_gmm_dir='model_gmm'
target_gmm='target_gmm/helix3_c100.txt'
out_dir='out_atomic_md100/'
total_mass=5485.224

# Sampling options
lowertemp=1 #float(sys.argv[1])
highertemp=1.1 #float(sys.argv[2])
nframes=10000 #int(sys.argv[3])
mdsteps=100 #int(sys.argv[4])
mcsteps=10
rb_max_trans=0.1
rb_max_rot=0.01
nframes_write_coordinates = 1
num_sd=0
num_cg=0

# EM options
cutoff_dist_mm=0
cutoff_dist_md=0
overlap_threshold=1e-4
target_radii_scale=2.0
em_weight=1
slope=1e-6
pointwise=True
local_mm=True
shuffle=False

# other restraint options
elastic_strength=50.0
elastic_cutoff=8.0





try:
  os.stat(out_dir)
except:
  os.mkdir(out_dir)


def setup_particles(ps):
    vxkey = IMP.FloatKey('vx')
    vykey = IMP.FloatKey('vy')
    vzkey = IMP.FloatKey('vz')
    for p in ps:
        IMP.core.XYZ(p).set_coordinates_are_optimized(True)
        #IMP.core.RigidBody(p).set_coordinates_are_optimized(True)
        p.add_attribute(vxkey, 0.0)
        p.add_attribute(vykey, 0.0)
        p.add_attribute(vzkey, 0.0)

############ SETUP COMPONENTS  #############
outputobjects = []
mc_objects = []
md_objects = []
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)
simo.create_component("chainA",color=0.0)
chainA=simo.add_component_pdb("chainA",in_fn,'A',
                              resolutions=[0],
                              #resrange=[55,80],
                              color=None,
                              offset=0,
                              cacenters=False,
                              show=False)
for a in IMP.core.get_leaves(chainA[0]):
  simo.set_floppy_bodies_from_hierarchies([a])

simo.set_super_rigid_body_from_hierarchies(chainA)
md_objects.append(simo)
ptsf=IMP.pmi.tools.ParticleToSampleFilter([simo])
ptsf.add_filter("SR_Bodies")
mc_objects.append(ptsf)


p=chainA[0].get_particle()
IMP.atom.Chain.setup_particle(p,'A')


############ SETUP RIGID BODIES  ############
#,output_map='atomic.mrc')

simo.set_current_coordinates_as_reference_for_rmsd('structure_native')
outputobjects.append(simo)
sse_selections=IMP.pmi.tools.parse_dssp(sse_fn,'A')


######### SETUP RESTRAINTS  ##########
### CHARMM
charmm=IMP.pmi.restraints.stereochemistry.CharmmForceFieldRestraint(simo)
charmm.add_to_model()
outputobjects.append(charmm)
print 'EVAL - charmm'
print m.evaluate(False)

#total_mass=sum((IMP.atom.Mass(p).get_mass() for p in IMP.core.get_leaves(simo.prot)))
#print total_mass



print 'going to set up all atom gaussians'
simo.add_all_atom_densities('chainA',chainA,output_map='test_map_1oel.mrc',voxel_size=2.0)
exit()
simo.set_super_rigid_bodies_max_trans(rb_max_trans)
simo.set_super_rigid_bodies_max_rot(rb_max_rot)
### finish setting up particles for MD
setup_particles(IMP.core.get_leaves(simo.prot))



### SSE elastic network
excluded_pairs=[]
print 'setting up sse restraints'
for n,sse in enumerate(sse_selections['helix']+sse_selections['beta']):
  sse_tuple=tuple([tuple(i) for i in sse])
  print 'adding tuple',sse_tuple
  #simo.set_super_rigid_bodies(sse_tuple)
  er=IMP.pmi.restraints.stereochemistry.ElasticNetworkRestraint(simo,sse_tuple,
                                                                strength=elastic_strength,
                                                                dist_cutoff=elastic_cutoff,
                                                                resolution=0,
                                                                ca_only=True)
  print 'done'
  er.set_weight(1.0)
  er.set_label('elastic_helix_%i'%(n))
  er.add_to_model()
  excluded_pairs+=er.get_excluded_pairs()
  outputobjects.append(er)
print 'EVAL - elastic network'
print m.evaluate(False)

### shuffle things
atoms=IMP.core.get_leaves(simo.prot)
cent=IMP.algebra.get_centroid([IMP.core.XYZ(p).get_coordinates() for p in atoms])
trans=IMP.algebra.get_random_local_transformation(cent,7.0,math.pi)
if shuffle:
    for p in atoms:
        IMP.core.transform(IMP.core.RigidBody(p),trans)




'''
print 'shuffling'
rex_shuffle=macros.ReplicaExchange0(m,
                            simo,
                            sampleobjects,
                            outputobjects,
                            sampler_type="MD",
                            replica_exchange_minimum_temperature=10.0,
                            replica_exchange_maximum_temperature=15.0,
                            number_of_best_scoring_models=10,
                            monte_carlo_steps=5000,
                            number_of_frames=1,
                            write_initial_rmf=True,
                            initial_rmf_name_suffix="initial",
                            stat_file_name_suffix="stat",
                            best_pdb_name_suffix="model",
                            do_clean_first=False,
                            do_create_directories=True,
                            global_output_directory=shuffle_dir)
rex_shuffle.execute_macro()
rex_object=rex_shuffle.get_replica_exchange_object()
print 'EVAL - shuffle'
print m.evaluate(False)
'''




### Gaussian EM



gem = IMP.pmi.restraints.em.GaussianEMRestraint(atoms,
                                                target_gmm,
                                                target_mass_scale=total_mass,
                                                cutoff_dist_model_model=cutoff_dist_mm,
                                                cutoff_dist_model_data=cutoff_dist_md,
                                                overlap_threshold=overlap_threshold,
                                                target_radii_scale=target_radii_scale,
                                                model_radii_scale=1.0,
                                                slope=slope,
                                                spherical_gaussians=True,
                                                pointwise_restraint=pointwise,
                                                local_mm=local_mm,
                                                close_pair_container=charmm.get_close_pair_container())
gem.set_weight(em_weight)
gem.add_to_model()
outputobjects.append(gem)
print 'EVAL - EM'
print m.evaluate(False)


########### SAMPLE ###########
'''
trans=IMP.algebra.Transformation3D(IMP.algebra.Vector3D(1,0,0))
for i in range(20):
    print 'trans',i,
    for p in atoms:
        IMP.core.transform(IMP.core.RigidBody(p),trans)
    print gem.evaluate()
'''

rex=IMP.pmi.macros.ReplicaExchange0(m,
                            simo,
                            sample_objects=None,
                            monte_carlo_sample_objects=mc_objects,
                            molecular_dynamics_sample_objects=md_objects,
                            output_objects=outputobjects,
                            replica_exchange_minimum_temperature=lowertemp,
                            replica_exchange_maximum_temperature=highertemp,
                            number_of_best_scoring_models=300,
                            monte_carlo_steps=mcsteps,
                            molecular_dynamics_steps=mdsteps,
                            number_of_frames=nframes,
                            write_initial_rmf=True,
                            initial_rmf_name_suffix="initial",
                            stat_file_name_suffix="stat",
                            best_pdb_name_suffix="model",
                            do_clean_first=False,
                            do_create_directories=True,
                            global_output_directory=out_dir,
                            atomistic=True)
#replica_exchange_object=rex_object)
rex.execute_macro()
