## \example pmi/em/em_rigid.py

import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.em
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.representation
import sys
import os

# SETTINGS ###############
in_fn = 'data/helix27.pdb'
target_gmm = 'target_gmm/helix27_c6.txt'
out_dir = 'out'
model_gmm_dir = 'model_gmm'
rbmaxtrans = 1.0
rbmaxrot = 0.025
nrmffiles = 1000
cutoff_dist = 1000.0
radii_mult = 3.0
lowertemp = 1.0
highertemp = 5.0
nframes = 1000
nsteps = 20

try:
    os.stat(out_dir)
except:
    os.mkdir(out_dir)

# SETUP COMPONENTS  #############
m = IMP.Model()
outputobjects = []
sampleobjects = []
simo = IMP.pmi.representation.Representation(
    m,
    upperharmonic=True,
    disorderedlength=False)
residues = []
res_densities = []
simo.create_component("chainA", color=0.0)
chainA = simo.add_component_pdb("chainA", in_fn, 'A',
                                resolutions=[1],
                                color=None,
                                offset=0,
                                cacenters=True,
                                show=False)
simo.setup_component_geometry("chainA")
chainA_density = simo.add_component_density('chainA', chainA,
                                            num_components=2, resolution=0,
                                            out_hier_name='chainA',
                                            inputfile='%s/helix27_c2.txt' % (model_gmm_dir))
                                         # outputfile='%s/helix27_c2.txt'%(model_gmm_dir),
                                         # outputmap='%s/helix27_c2.mrc'%(model_gmm_dir),
                                         # multiply_by_total_mass=True)
res_densities += chainA_density
simo.set_rigid_body_from_hierarchies(chainA + chainA_density)
simo.set_current_coordinates_as_reference_for_rmsd('structure_native')
simo.set_rigid_bodies_max_rot(rbmaxrot)
simo.set_rigid_bodies_max_trans(rbmaxtrans)
# simo.set_floppy_bodies()
outputobjects.append(simo)
sampleobjects.append(simo)
simo.shuffle_configuration(10)
# SETUP RESTRAINTS  ##########
# connectivity
simo.setup_component_sequence_connectivity("chainA", resolution=1.0)
all_densities = set(
    p for comp in res_densities for p in IMP.core.get_leaves(comp))
total_mass = sum((IMP.atom.Mass(p).get_mass() for p in all_densities))

# EM
print 'total mass', total_mass
gem = IMP.pmi.restraints.em.GaussianEMRestraint(res_densities,
                                                target_gmm,
                                                target_mass_scale=total_mass,
                                                cutoff_dist_for_container=cutoff_dist,
                                                target_radii_scale=radii_mult,
                                                model_radii_scale=radii_mult)
gem.set_weight(100)
gem.add_to_model()
outputobjects.append(gem)
print 'eval 1'
print m.evaluate(False)

# excluded volume
ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(
    simo,
    resolution=1.0)
ev.add_to_model()
outputobjects.append(ev)
m.update()

# SAMPLE ###########
rex = IMP.pmi.macros.ReplicaExchange0(m,
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
