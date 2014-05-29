import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.em
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.representation as representation
import IMP.pmi.tools as tools
import IMP.pmi.samplers as samplers
import IMP.pmi.output as output
import IMP.pmi.macros as macros
import sys
import os
from collections import deque
from math import pi
############ SETTINGS ###############
m = IMP.Model()
in_fn='data/helix3.pdb'
sse_fn='data/helix3.dssp'
target_gmm='target_gmm/helix3_c20.txt'
total_mass=5485.224

target_map='data/helix3_4.mrc'

#in_fn='data/1oel_A_fit4.pdb'
#sse_fn='data/1oel_A.dssp'
#target_gmm='target_gmm/1oel_A_c100.txt'
#total_mass=50990.8117

out_dir='out_atomic_md'
#shuffle_dir='out_atomic_md/shuffle/'
model_gmm_dir='model_gmm'
cutoff_dist_em=5.0
model_radii_scale=1.0
target_radii_scale=3.0
lowertemp=1.0
highertemp=1.3
nframes=1000
nsteps=200
pointwise=True

elastic_strength=50.0
elastic_cutoff=8.0
em_weight=10

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
sampleobjects = []
simo = representation.Representation(m,upperharmonic=True,disorderedlength=False)
simo.create_component("chainA",color=0.0)
chainA=simo.add_component_pdb("chainA",in_fn,'A',
                              resolutions=[0],
                              #resrange=[55,80],
                              color=None,
                              offset=0,
                              cacenters=False,
                              show=False)
# each atom is a rigid body
for a in IMP.core.get_leaves(chainA[0]):
  simo.set_floppy_bodies_from_hierarchies([a])

# each residue is a super rigid body...

p=chainA[0].get_particle()
IMP.atom.Chain.setup_particle(p,'A')


############ SETUP RIGID BODIES  ############
#,output_map='atomic.mrc')

#simo.set_super_rigid_body_from_hierarchies(chainA)
simo.set_current_coordinates_as_reference_for_rmsd('structure_native')
outputobjects.append(simo)
sampleobjects.append(simo)
sse_selections=tools.parse_dssp(sse_fn,'A')


######### SETUP RESTRAINTS  ##########
#all_densities=set(p for comp in res_densities for p in IMP.core.get_leaves(comp))
#total_mass=sum((IMP.atom.Mass(p).get_mass() for p in all_densities))

### CHARMM
charmm=IMP.pmi.restraints.stereochemistry.CharmmForceFieldRestraint(simo)
charmm.add_to_model()
outputobjects.append(charmm)
print 'EVAL - charmm'
print m.evaluate(False)


print 'going to set up all atom gaussians'
simo.add_all_atom_densities('chainA',chainA)
### finish setting up particles for MD
setup_particles(IMP.core.get_leaves(simo.prot))



### SSE elastic network
excluded_pairs=[]
print 'setting up sse restraints'
for n,sse in enumerate(sse_selections['helix']+sse_selections['beta']):
  sse_tuple=tuple([tuple(i) for i in sse])
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
trans=IMP.algebra.get_random_local_transformation(cent,5.0,pi/3)
for p in atoms:
    #print IMP.core.RigidBody.get_is_setup(p)
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
                                                cutoff_dist_for_container=cutoff_dist_em,
                                                target_radii_scale=target_radii_scale,
                                                model_radii_scale=model_radii_scale,
                                                spherical_gaussians=True,
                                                pointwise_restraint=pointwise)
gem.set_weight(em_weight)
gem.add_to_model()
outputobjects.append(gem)
print 'EVAL - EM'
print m.evaluate(False)


### Traditional EM
'''
em = IMP.pmi.restraints.em.EMRestraint(atoms,
                                       target_fn=target_map,
                                       map_resolution=4)
em.set_weight(em_weight)
em.add_to_model()
outputobjects.append(em)
print 'EVAL - EM'
print m.evaluate(False)
'''


########### SAMPLE ###########

rex=macros.ReplicaExchange0(m,
                            simo,
                            sampleobjects,
                            outputobjects,
                            sampler_type="MD",
                            replica_exchange_minimum_temperature=lowertemp,
                            replica_exchange_maximum_temperature=highertemp,
                            number_of_best_scoring_models=100,
                            monte_carlo_steps=nsteps,
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




########### SAMPLE ###########
'''
output = IMP.pmi.output.Output()
output.init_pdb_best_scoring(out_dir+'/pdbs/model',
                            simo.prot,
                            100,
                            replica_exchange=False)
output.init_rmf(out_dir+"/rmfs/0.rmf3",[simo.prot])

#bd = IMP.atom.BrownianDynamics(m)
#bd.set_time_step(10)
for i in range(0, 1000):
   md.optimize(200)
   #bd.optimize(100)

   score=m.evaluate(False)
   output.write_pdb_best_scoring(score)
   output.write_rmfs()
   print i,score
'''

'''
set_diffusion_factor(IMP.core.get_leaves(simo.prot),diffusion_factor)
output = IMP.pmi.output.Output()
cg = IMP.core.ConjugateGradients(m)
output.init_rmf(out_dir+"/trajectory.rmf3",[simo.prot])

output.init_pdb_best_scoring(out_dir+'/pdbs/model',
                             simo.prot,
                             100,
                             replica_exchange=False)

scores=deque(maxlen=20)
bd_nsteps=25
for i in range(0, 1000):
    bd.optimize(bd_nsteps)
    score=m.evaluate(False)
    print i
    print 'BD:',score
    cg.optimize(20)
    score=m.evaluate(False)
    print 'CG:',score
    scores.append(score)
    spread=abs(max(scores)-min(scores))/min(scores)
    print '  spread',spread
    if abs(spread)<spread_threshold:
      bd_nsteps=int(1.1*bd_nsteps)
      print 'increasing bd_nsteps to',bd_nsteps
    else:
      bd_nsteps=max(10,int(0.9*bd_nsteps))
      if bd_nsteps>10:
        print 'decreasing bd_nsteps to',bd_nsteps
    output.write_rmfs()
    output.write_pdb_best_scoring(score)
'''
