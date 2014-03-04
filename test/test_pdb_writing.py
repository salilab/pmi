import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.stereochemistry
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.samplers
import IMP.pmi.output

#input parameter

pdbfile=IMP.pmi.get_data_path("1WCM.pdb")
fastafile=IMP.pmi.get_data_path("1WCM.fasta.txt")

components=["Rpb1","Rpb2","Rpb3", "Rpb4"]
            
chains="ABCDEFGHIJKL"

colors=[ 0. , 0.5,  0.75, 1.0]

beadsize=20

fastids=IMP.pmi.tools.get_ids_from_fasta_file(fastafile)



m=IMP.Model()
simo = IMP.pmi.representation.SimplifiedModel(m) 


simo.add_component_name("Rpb1",color=colors[0])
simo.add_component_sequence("Rpb1",fastafile,id=fastids[0])
simo.autobuild_model("Rpb1",pdbfile,"A",
                                         resolutions=[1],beadsize=beadsize)
simo.setup_component_sequence_connectivity("Rpb1",1)

simo.add_component_name("Rpb2",color=colors[1])
simo.add_component_sequence("Rpb2",fastafile,id=fastids[1])
simo.autobuild_model("Rpb2",pdbfile,"B",
                                         resolutions=[10],beadsize=beadsize)
simo.setup_component_sequence_connectivity("Rpb2",1)

simo.add_component_name("Rpb3",color=colors[2])
simo.add_component_sequence("Rpb3",fastafile,id=fastids[2])
simo.autobuild_model("Rpb3",pdbfile,"C",
                                         resolutions=[0],beadsize=beadsize)
simo.setup_component_sequence_connectivity("Rpb3",1)

simo.add_component_name("Rpb4",color=colors[3])
simo.add_component_sequence("Rpb4",fastafile,id=fastids[3])
simo.autobuild_model("Rpb4",pdbfile,"D",
                                         resolutions=[0,1],beadsize=beadsize)
simo.setup_component_sequence_connectivity("Rpb4",1)

output=IMP.pmi.output.Output()
output.init_pdb("test_pdb_writing.pdb",simo.prot)
output.write_pdbs()
output.init_pdb_best_scoring("test_pdb_writing",simo.prot,10)
for i in range(20):
    score=-float(i)
    output.write_pdb_best_scoring(score)   
    

