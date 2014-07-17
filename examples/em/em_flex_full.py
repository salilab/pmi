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
target_gmm10='target_gmm/helix3_c10.txt'
target_gmm5='target_gmm/helix3_c5.txt'
out_dir='out_flex'
model_gmm_dir='model_gmm'
rbmaxtrans = 1.0
rbmaxrot=0.025
nrmffiles=1000
cutoff_dist=0.0
radii_mult=2.0
lowertemp=1.0
highertemp=4.0
nframes=1000
nsteps=20

em_weight20=70
em_weight10=20
em_weight5=10
loop_bond_jitter=0.0
loop_angle_min=100.0 * 1.0
loop_angle_max=140 * 1.0
bond_strength=10.0
angle_strength=10.0
dihedral_strength=10.0
elastic_strength=100.0
elastic_cutoff=8.0

resolutions=[0,1]

############ SETUP COMPONENTS  #############
m = IMP.Model()
outputobjects = []
sampleobjects = []
simo = representation.Representation(m,upperharmonic=True,disorderedlength=False)
sse_selections=tools.parse_dssp(sse_fn,'A')
try:
  os.stat(out_dir)
except:
  os.mkdir(out_dir)

### read one residue at a time!
# read ahead to get residue numbering
mh=IMP.atom.read_pdb(in_fn,m,IMP.atom.CAlphaPDBSelector())
rnums=[IMP.atom.Residue(r).get_index() for r in IMP.atom.get_by_type(mh,IMP.atom.RESIDUE_TYPE)]
name='chainA'
simo.create_component(name,color=0.0)
all_hiers=[]
all_den=[]
each_res=[]
for idx,rnum in enumerate(rnums):
  rname='chainA_res%i'%rnum
  res=simo.add_component_pdb(name,in_fn,chain='A',resrange=(rnum,rnum),
                             resolutions=resolutions,
                             color=None,
                             offset=0,
                             cacenters=True,
                             show=False,
                             read_ca_cb_only=True)
  simo.setup_component_geometry(name)
  den=simo.add_component_density(name,hierarchies=res,
                                 simulation_res=0.5,
                                 num_components=2,
                                 inputfile=os.path.join(model_gmm_dir,rname+'.txt'),
                                 #outputfile=os.path.join(model_gmm_dir,name+'.txt'),
                                 #outputmap=os.path.join(model_gmm_dir,name+'.mrc'),
                                 sampled_points=100000)
  print [p.get_name() for p in IMP.core.get_leaves(res[0])]
  print [p.get_name() for p in IMP.core.get_leaves(den[0])]
  simo.set_rigid_body_from_hierarchies(res+den)
  simo.set_super_rigid_body_from_hierarchies(res+den)
  each_res.append(res+den)
  all_hiers+=res
  all_hiers+=den
  all_den+=den
simo.setup_bonds()

#for ssize in range(2,6):
#  for i in range(0,len(each_res)-ssize):
#    simo.set_super_rigid_body_from_hierarchies([j for r in each_res[i:i+ssize] for j in r])


### super rigid bodies
simo.set_super_rigid_body_from_hierarchies(all_hiers)
#simo.set_chain_of_super_rigid_bodies(all_hiers,6,6)
simo.set_current_coordinates_as_reference_for_rmsd('structure_native')
simo.set_rigid_bodies_max_rot(rbmaxrot)
simo.set_rigid_bodies_max_trans(rbmaxtrans)
simo.set_floppy_bodies()
outputobjects.append(simo)
sampleobjects.append(simo)

### gather secondary structure info
#ss=['-']*num_res
ss=['T']*len(rnums)
hnums=[]
for n,sse in enumerate(sse_selections['helix']):
  for ns in range(sse[0][0],sse[0][1]):
    #ss[ns]='H'
    ss[ns]='C'
ss=''.join(ss)
print ss

######### SETUP RESTRAINTS  ##########
excluded_pairs=[]

#all_densities=set(p for comp in res_densities for p in IMP.core.get_leaves(comp))
#total_mass=sum((IMP.atom.Mass(p).get_mass() for p in all_densities))
#print 'total mass is',total_mass
total_mass=11269.2246

### EM
print 'total mass',total_mass
gem = IMP.pmi.restraints.em.GaussianEMRestraint(all_den,
                                                target_gmm,
                                                target_mass_scale=total_mass,
                                                cutoff_dist_for_container=cutoff_dist,
                                                target_radii_scale=radii_mult,
                                                model_radii_scale=radii_mult)
gem.set_weight(em_weight20)
gem.add_to_model()
outputobjects.append(gem)
print 'EVAL - EM'
print m.evaluate(False)


### pseudoatomic angles, bonds
par=IMP.pmi.restraints.stereochemistry.PseudoAtomicRestraint(rnums,simo,'chainA',kappa=10.0,
                                                             jitter_angle=0.34,
                                                             jitter_improper=0.34)
par.set_label('PseudoAtomicRestraints')
par.set_weight(1.0)
par.add_to_model()
excluded_pairs += par.get_excluded_pairs()
outputobjects.append(par)
print 'EVAL STEREOCHEMISTRY'
print m.evaluate(False)

### resolution 1 restraints
# dihedrals
ldr = IMP.pmi.restraints.stereochemistry.ResidueDihedralRestraint(simo, 'chainA',ss,
                                                                    strength=100.0)
ldr.set_label('DIHEDRALS_')
ldr.set_weight(1.0)
ldr.add_to_model()
excluded_pairs+= ldr.get_excluded_pairs()
outputobjects.append(ldr)

# angles
rar = IMP.pmi.restraints.stereochemistry.ResidueAngleRestraint(simo, 'chainA',
                                                               strength=angle_strength,
                                                               anglemin=loop_angle_min,
                                                               anglemax=loop_angle_max)
rar.set_label('ANGLES')
rar.set_weight(1.0)
rar.add_to_model()
excluded_pairs += rar.get_excluded_pairs()
outputobjects.append(rar)

# bonds
rbr = IMP.pmi.restraints.stereochemistry.ResidueBondRestraint(simo, 'chainA',
                                                              strength=bond_strength)
rbr.set_label('BONDS')
rbr.set_weight(1.0)
rbr.add_to_model()
excluded_pairs += rbr.get_excluded_pairs()
outputobjects.append(rbr)



### SSE elastic network
print 'setting up sse restraints'
for n,sse in enumerate(sse_selections['helix']):
  sse_tuple=tuple([tuple(i) for i in sse])
  shiers=IMP.pmi.tools.select_by_tuple(simo,sse_tuple[0])
  simo.set_super_rigid_body_from_hierarchies(shiers)
  er=IMP.pmi.restraints.stereochemistry.ElasticNetworkRestraint(simo,sse_tuple,
                                                                strength=elastic_strength,
                                                                dist_cutoff=elastic_cutoff)
  print 'done'
  er.set_weight(1.0)
  er.set_label('helix_%i'%(n))
  er.add_to_model()
  excluded_pairs+=er.get_excluded_pairs()
  outputobjects.append(er)
print 'EVAL - elastic network'
print m.evaluate(False)


# excluded volume
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo, resolution=1.0)
ev.add_to_model()
ev.add_excluded_particle_pairs(excluded_pairs)
outputobjects.append(ev)
m.update()

simo.shuffle_configuration(5)
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
