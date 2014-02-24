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
target_gmm='target_gmm/helix3_c20.txt'
out_dir='out_beads'
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

############ SETUP COMPONENTS  #############
m = IMP.Model()
outputobjects = []
sampleobjects = []
simo = representation.Representation(m,upperharmonic=True,disorderedlength=False)
simo.create_component("chainA",color=0.0)

############ SETUP DENSITIES  ############
ps=[]
xyzrs=IMP.core.create_xyzr_particles(m,10,5,10)
for xyzr in xyzrs:
    p=xyzr.get_particle()
    print 'creating',xyzr
    m=IMP.atom.Mass.setup_particle(p,10)
    ps.append(p)

bead_density=simo.add_component_density('chainA',particles=ps,
                                        num_components=40,resolution=0,
                                        out_hier_name='chainA',
                                        kernel_type=IMP.em.SPHERE,
                                        #inputfile='%s/%s.txt'%(model_gmm_dir,name))
                                        outputfile='%s/test_beads.txt'%(model_gmm_dir),
                                        outputmap='%s/test_beads.mrc'%(model_gmm_dir),
                                        multiply_by_total_mass=True,
                                        intermediate_map_fn='intermediate_beads.mrc')
