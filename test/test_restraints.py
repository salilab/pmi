import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints
import IMP.pmi.representation
import IMP.pmi.stereochemistry
import IMP.pmi.basic
import IMP.pmi.proteomics
import IMP.pmi.crosslinking

#input parameter

pdbfile=IMP.pmi.get_data_path("1WCM.pdb")
fastafile=IMP.pmi.get_data_path("1WCM.fasta.txt")

components=["Rpb3","Rpb3.copy","Rpb4"]
            
chains="CCD"

colors=[ 0.,0.5,1.0]

beadsize=20

fastids=IMP.pmi.tools.get_ids_from_fasta_file(fastafile)

m=IMP.Model()
simo = IMP.pmi.representation.SimplifiedModel(m) 

hierarchies={}

for n in range(len(components)):
    simo.add_component_name(components[n],color=colors[n])
    simo.add_component_sequence(components[n],fastafile,id=fastids[n+2])
    hierarchies[components[n]]=simo.autobuild_pdb_and_intervening_beads(components[n],pdbfile,chains[n],
                                             resolutions=[1,10,100],beadsize=beadsize)
    simo.setup_component_sequence_connectivity(components[n],1)

ev1 = IMP.pmi.stereochemistry.ExcludedVolumeSphere(simo)
ev1.add_to_model()
print ev1.get_output()
print m.evaluate(False)

ev2 = IMP.pmi.stereochemistry.ExcludedVolumeSphere(simo,[simo.hier_dict["Rpb3"]],[simo.hier_dict["Rpb4"]])
ev2.add_to_model()
print ev2.get_output()
print m.evaluate(False)

rb = IMP.pmi.stereochemistry.ResidueBondRestraint(simo,(30,40,"Rpb3"))
rb.add_to_model()
print rb.get_output()
print m.evaluate(False)

ra = IMP.pmi.stereochemistry.ResidueAngleRestraint(simo,(30,40,"Rpb3"))
ra.add_to_model()
print ra.get_output()
print m.evaluate(False)

rd = IMP.pmi.stereochemistry.ResidueDihedralRestraint(simo,(30,40,"Rpb3"))
rd.add_to_model()
print ra.get_output()
print m.evaluate(False)

ss = IMP.pmi.stereochemistry.SecondaryStructure(simo,(30,40,"Rpb3"),"HHHHHHHHHHH")
ss.add_to_model()
print ss.get_output()
print m.evaluate(False)

eb1 = IMP.pmi.basic.ExternalBarrier(simo,50)
eb1.add_to_model()
print eb1.get_output()
print m.evaluate(False)

cr1  = IMP.pmi.proteomics.ConnectivityRestraint(simo,[(1,100,"Rpb3"),(1,100,"Rpb4")],resolution=100)
cr1.add_to_model()
print cr1.get_output()
print m.evaluate(False)

cr2 = IMP.pmi.proteomics.CompositeRestraint(simo,[(1,100,"Rpb3"),(200,300,"Rpb3")],
                                            [[(1,100,"Rpb4")],[(200,500,"Rpb4")]],resolution=100)
cr2.add_to_model()
print cr2.get_output()
print m.evaluate(False)

restraints='''#
Rpb3 Rpb4 100 150
Rpb4 Rpb4 50 150'''

cr2 = IMP.pmi.proteomics.AmbiguousCompositeRestraint(simo,restraints,resolution=1)
cr2.add_to_model()
print cr2.get_output()
print m.evaluate(False)


restraints='''#
Rpb3 Rpb4 100 150 0.5
Rpb4 Rpb4 50 150 0.7'''

pm = IMP.pmi.proteomics.SimplifiedPEMAP(simo,restraints,20,1,resolution=1)
pm.add_to_model()
print pm.get_output()
print m.evaluate(False)


restraints='''#
Rpb3 Rpb4 100 150
Rpb4 Rpb4 50 150 '''

xl1 = IMP.pmi.crosslinking.ConnectivityCrossLinkMS(simo,restraints,20,1,resolution=1)
xl1.add_to_model()
print xl1.get_output()
print m.evaluate(False)

restraints='''#
Rpb3 Rpb4 100 150
Rpb4 Rpb4 50 150 '''

xl2 = IMP.pmi.crosslinking.SimplifiedCrossLinkMS(simo,restraints,25,1,resolution=1)
xl2.add_to_model()
print xl2.get_output()
print m.evaluate(False)

restraints='''#
Rpb3 Rpb4 100 150
Rpb4 Rpb4 50 150 '''

xl3 = IMP.pmi.crosslinking.SigmoidCrossLinkMS(simo,restraints,25,1,5,0.02,resolution=1)
xl3.add_to_model()
print xl3.get_output()
print m.evaluate(False)

restraints='''#
Rpb3 Rpb4 100 150 1
Rpb4 Rpb4 50 150 1'''

xl4 = IMP.pmi.crosslinking.ISDCrossLinkMS(simo,restraints,25,resolution=1)
xl4.add_to_model()
print xl4.get_output()
print m.evaluate(False)

restraints='''#
100 Rpb3 150 Rpb4 0.5 epsilon1
50 Rpb4 150 Rpb4  0.7 epsilon2'''

xl5 = IMP.pmi.crosslinking.CysteineCrossLinkRestraint([simo],restraints)
xl5.add_to_model()
xl5.set_output_level("high")
print xl5.get_output()
print m.evaluate(False)