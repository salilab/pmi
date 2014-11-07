## \example running the sampling

import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.restraints.basic
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output
import IMP.pmi.macros

import os

# setting up parameters

rbmaxtrans = 2.00
fbmaxtrans = 3.00
rbmaxrot=0.03
outputobjects = []
sampleobjects = []
nbestscoringmodels=50
nframes=1000
nsteps=10

# setting up topology

m = IMP.Model()
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)
execfile("topology.py")
total_mass=sum((IMP.atom.Mass(p).get_mass() for h in resdensities for p in IMP.atom.get_leaves(h)))


# randomize the initial configuration

simo.shuffle_configuration(100)

# defines the movers

simo.set_rigid_bodies_max_rot(rbmaxrot)
simo.set_floppy_bodies_max_trans(fbmaxtrans)
simo.set_rigid_bodies_max_trans(rbmaxtrans)
simo.set_floppy_bodies()
simo.setup_bonds()

outputobjects.append(simo)
sampleobjects.append(simo)

# scoring function

ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo,resolution=10)
ev.add_to_model()
outputobjects.append(ev)


columnmap={}
columnmap["Protein1"]="pep1.accession"
columnmap["Protein2"]="pep2.accession"
columnmap["Residue1"]="pep1.xlinked_aa"
columnmap["Residue2"]="pep2.xlinked_aa"
columnmap["IDScore"]=None
columnmap["XLUniqueID"]=None

ids_map=IMP.pmi.tools.map()
ids_map.set_map_element(1.0,1.0)

xl1 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
                                   '../data/polii_xlinks.csv',
                                   length=21.0,
                                   slope=0.01,
                                   columnmapping=columnmap,
                                   ids_map=ids_map,
                                   resolution=1.0,
                                   label="Mike",
                                   csvfile=True)
xl1.add_to_model()
xl1.set_label("Mike")
sampleobjects.append(xl1)
outputobjects.append(xl1)
xl1.set_psi_is_sampled(False)
psi=xl1.get_psi(1.0)[0]
psi.set_scale(0.05)










columnmap={}
columnmap["Protein1"]="prot1"
columnmap["Protein2"]="prot2"
columnmap["Residue1"]="res1"
columnmap["Residue2"]="res2"
columnmap["IDScore"]=None
columnmap["XLUniqueID"]=None

ids_map=IMP.pmi.tools.map()
ids_map.set_map_element(1.0,1.0)

xl2 = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(simo,
                                   '../data/polii_juri.csv',
                                   length=21.0,
                                   slope=0.01,
                                   columnmapping=columnmap,
                                   ids_map=ids_map,
                                   resolution=1.0,
                                   label="Juri",
                                   csvfile=True)
xl2.add_to_model()
xl2.set_label("Juri")
sampleobjects.append(xl2)
outputobjects.append(xl2)
xl2.set_psi_is_sampled(False)
psi=xl2.get_psi(1.0)[0]
psi.set_scale(0.05)



####################################################
# Monte Carlo
####################################################

mc = IMP.pmi.samplers.MonteCarlo(m,sampleobjects, 0.5)
mc.set_simulated_annealing(1.0, 10.0, 100, 10)
outputobjects.append(mc)

####################################################
# Output
####################################################

sw = IMP.pmi.tools.Stopwatch()
outputobjects.append(sw)

output = IMP.pmi.output.Output()
output.init_stat2("stat.out", outputobjects,extralabels=["rmf_file","rmf_frame_index"])
output.init_rmf("initial.rmf", [simo.prot])
output.write_rmf("initial.rmf")

output.init_pdb_best_scoring("pdbs/model",simo.prot,nbestscoringmodels)

#####################################################
#running simulation
#####################################################

rmffile="models.rmf"
output.init_rmf(rmffile, [simo.prot])
output.add_restraints_to_rmf(rmffile,[xl1,xl2])

for i in range(nframes):
    mc.optimize(nsteps)
    score = m.evaluate(False)
    print mc.get_frame_number()
    output.set_output_entry("rmf_file",rmffile)
    output.set_output_entry("rmf_frame_index",i)
    output.write_stats2()
    output.write_pdb_best_scoring(score)
    output.write_rmf(rmffile)
