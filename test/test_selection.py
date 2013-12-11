import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints as restraints
import IMP.pmi.representation as representation
import IMP.pmi.tools as tools
import IMP.pmi.samplers as samplers
import IMP.pmi.output as output

#input parameter

pdbfile=IMP.pmi.get_data_path("1WCM.pdb")
fastafile=IMP.pmi.get_data_path("1WCM.fasta.txt")

components=["Rpb3","Rpb3.copy","Rpb4"]
            
chains="CCD"

colors=[ 0.,0.5,1.0]

beadsize=20

fastids=tools.get_ids_from_fasta_file(fastafile)

m=IMP.Model()
simo = representation.SimplifiedModel(m) 

hierarchies={}

for n in range(len(components)):
    simo.add_component_name(components[n],color=colors[n])
    simo.add_component_sequence(components[n],fastafile,id=fastids[n+2])
    hierarchies[components[n]]=simo.autobuild_pdb_and_intervening_beads(components[n],pdbfile,chains[n],
                                             resolutions=[1,10,100],beadsize=beadsize)
    simo.setup_component_sequence_connectivity(components[n],1)
    
print "All",len(tools.select(simo))
print "resolution=1",len(tools.select(simo,resolution=1))
print "resolution=1,resid=10", len(tools.select(simo,resolution=1,residue=10))
print "resolution=1,resid=10,name=Rpb3",len(tools.select(simo,resolution=1,name="Rpb3",residue=10))
print "resolution=1,resid=10,name=Rpb3,ambiguous",len(tools.select(simo,resolution=1,name="Rpb3",name_is_ambiguous=True,residue=10))
print "resolution=1,resid=10,name=Rpb4,ambiguous",len(tools.select(simo,resolution=1,name="Rpb4",name_is_ambiguous=True,residue=10))
print "resolution=1,resrange=(10,20),name=Rpb3",len(tools.select(simo,resolution=1,name="Rpb3",first_residue=10,last_residue=20))
print "resolution=1,resrange=(10,20),name=Rpb3,ambiguous",len(tools.select(simo,resolution=1,name="Rpb3",name_is_ambiguous=True,first_residue=10,last_residue=20))
print "resolution=10,resrange=(10,20),name=Rpb3",len(tools.select(simo,resolution=10,name="Rpb3",first_residue=10,last_residue=20))
print "resolution=10,resrange=(10,20),name=Rpb3,ambiguous",len(tools.select(simo,resolution=10,name="Rpb3",name_is_ambiguous=True,first_residue=10,last_residue=20))
print "resolution=100,resrange=(10,20),name=Rpb3",len(tools.select(simo,resolution=100,name="Rpb3",first_residue=10,last_residue=20))
print "resolution=100,resrange=(10,20),name=Rpb3,ambiguous",len(tools.select(simo,resolution=100,name="Rpb3",name_is_ambiguous=True,first_residue=10,last_residue=20))

for key in hierarchies:
    for h in hierarchies[key]:
       print h.get_name(),"resolution=1", len(tools.select(simo,resolution=1,hierarchies=[h]))
       print h.get_name(),"resolution=10", len(tools.select(simo,resolution=10,hierarchies=[h]))
       print h.get_name(),"resolution=100", len(tools.select(simo,resolution=100,hierarchies=[h]))

print "Beads",len(tools.select(simo,representation_type="Beads"))
print "Molecule",len(tools.select(simo,representation_type="Molecule"))
print "resolution=1,Molecule",len(tools.select(simo,resolution=1,representation_type="Molecule"))
print "resolution=10,Molecule",len(tools.select(simo,resolution=10,representation_type="Molecule"))
print "resolution=100,Molecule",len(tools.select(simo,resolution=100,representation_type="Molecule"))
print "resolution=1,Beads",len(tools.select(simo,resolution=1,representation_type="Beads"))
print "resolution=10,Beads",len(tools.select(simo,resolution=10,representation_type="Beads"))
print "resolution=100,Beads",len(tools.select(simo,resolution=100,representation_type="Beads"))
print "resolution=2",len(tools.select(simo,resolution=2))

print "resolution=7",len(tools.select(simo,resolution=7))
print "resolution=10",len(tools.select(simo,resolution=10))
print "resolution=100",len(tools.select(simo,resolution=100))