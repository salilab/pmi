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

############ SETTINGS ###############
in_fn='data/helix3.pdb'
sse_fn='data/helix3.dssp'
target_gmm='target_gmm/helix3_c20.txt'
out_dir='out_flex'
model_gmm_dir='model_gmm'
rbmaxtrans = 1.0
rbmaxrot=0.025
nrmffiles=1000
cutoff_dist=1000.0
radii_mult=3.0
lowertemp=1.0
highertemp=5.0
nframes=1000
nsteps=20

try:
  os.stat(out_dir)
except:
  os.mkdir(out_dir)

############ SETUP COMPONENTS  #############
m = IMP.Model()
outputobjects = []
sampleobjects = []
simo = representation.Representation(m,upperharmonic=True,disorderedlength=False)
simo.create_component("chainA",color=0.0)
chainA=simo.add_component_pdb("chainA",in_fn,'A',
                              resolutions=[1],
                              color=None,
                              offset=0,
                              cacenters=True,
                              show=False)
simo.setup_component_geometry("chainA")
simo.setup_bonds()


############ SETUP DENSITIES  ############
sse_selections=tools.parse_dssp(sse_fn,'A')
res_densities=[]

for key in ('helix','beta'):
  for sse_tuple in sse_selections[key]:
    ps=simo.get_particles_from_selection_tuples(sse_tuple)
    name=key
    for t in sse_tuple:
      name+='_%s-%i-%i'%(t[2].strip('chain'),t[0],t[1])
    res_density=simo.add_component_density(sse_tuple[0][2],particles=ps,
                                           num_components=4,resolution=0,
                                           out_hier_name=name,
                                           inputfile='%s/%s.txt'%(model_gmm_dir,name))
                                           #outputfile='%s/%s.txt'%(model_gmm_dir,name),
                                           #outputmap='%s/%s.mrc'%(model_gmm_dir,name),
                                           #multiply_by_total_mass=True)
    res_densities+=res_density
    simo.set_rigid_body_from_hierarchies(res_density,particles=ps)

# add density for each individual loop residue
for nl,l in enumerate(sse_selections['loop']):
  start,stop,chain=l[0]
  for r in range(start,stop+1):
    sse_tuple=[[r,r,chain]]
    ps=simo.get_particles_from_selection_tuples(sse_tuple)
    name='res'
    for t in sse_tuple:
      name+='_%s-%i'%(t[2].strip('chain'),t[0])
    res_density=simo.add_component_density(sse_tuple[0][2],particles=ps,
                                           num_components=2,resolution=0,
                                           out_hier_name=name,
                                           inputfile='%s/%s.txt'%(model_gmm_dir,name))
                                           #outputfile='%s/%s.txt'%(model_gmm_dir,name),
                                           #outputmap='%s/%s.mrc'%(model_gmm_dir,name),
                                           #multiply_by_total_mass=True)
    res_densities+=res_density
    simo.set_rigid_body_from_hierarchies(res_density,particles=ps)



# rigid body params
simo.set_super_rigid_body_from_hierarchies(chainA)
simo.set_current_coordinates_as_reference_for_rmsd('structure_native')
simo.set_rigid_bodies_max_rot(rbmaxrot)
simo.set_rigid_bodies_max_trans(rbmaxtrans)
simo.set_floppy_bodies()
outputobjects.append(simo)
sampleobjects.append(simo)
simo.shuffle_configuration(10)

######### SETUP RESTRAINTS  ##########
# connectivity
simo.setup_component_sequence_connectivity("chainA",resolution=1.0)
#all_densities=set(p for comp in res_densities for p in IMP.core.get_leaves(comp))
#total_mass=sum((IMP.atom.Mass(p).get_mass() for p in all_densities))
#print 'total mass is',total_mass
total_mass=11269.2246

# EM
print 'total mass',total_mass
gem = IMP.pmi.restraints.em.GaussianEMRestraint(res_densities,
                                                target_gmm,
                                                target_mass_scale=total_mass,
                                                cutoff_dist_for_container=cutoff_dist,
                                                target_radii_scale=radii_mult,
                                                model_radii_scale=radii_mult)
gem.set_weight(50)
gem.add_to_model()
outputobjects.append(gem)
print 'EVAL 1'
print m.evaluate(False)

# add bonds and angles
listofexcludedpairs = []
for l in sse_selections['loop']:
  this_loop=(l[0][0]-2,l[0][1]+2,l[0][2]) #extend loops
  rbr = IMP.pmi.restraints.stereochemistry.ResidueBondRestraint(simo, this_loop, strength=100.0,
                                                                jitter=0.7)
  rbr.set_label('BOND_'+str(this_loop))
  rbr.set_weight(1.0)
  rbr.add_to_model()
  listofexcludedpairs += rbr.get_excluded_pairs()
  outputobjects.append(rbr)

  rar = IMP.pmi.restraints.stereochemistry.ResidueAngleRestraint(simo, this_loop,strength=100.0,
                                                                 anglemin=100.0,anglemax=140.0)
  rar.set_label('ANGLE_'+str(this_loop))
  rar.set_weight(1.0)
  rar.add_to_model()
  listofexcludedpairs += rar.get_excluded_pairs()
  outputobjects.append(rar)
print 'EVAL 2'
print m.evaluate(False)

# excluded volume
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo, resolution=1.0)
ev.add_to_model()
ev.add_excluded_particle_pairs(listofexcludedpairs)
outputobjects.append(ev)
m.update()

########### SAMPLE ###########
rex=macros.ReplicaExchange0(m,
                            simo,
                            sampleobjects,
                            outputobjects,
                            replica_exchange_minimum_temperature=lowertemp,
                            replica_exchange_maximum_temperature=highertemp,
                            number_of_best_scoring_models=500,
                            monte_carlo_steps=nsteps,
                            number_of_frames=nframes,
                            write_initial_rmf=True,
                            initial_rmf_name_suffix="initial",
                            stat_file_name_suffix="stat",
                            best_pdb_name_suffix="model",
                            do_clean_first=False,
                            do_create_directories=True,
                            global_output_directory=out_dir)

rex.execute_macro()
