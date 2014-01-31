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

components=["Rpb1","Rpb2","Rpb3","Rpb4",
            "Rpb5","Rpb6","Rpb7","Rpb8",
            "Rpb9","Rpb10","Rpb11","Rpb12"]
            
chains="ABCDEFGHIJKL"

colors=[ 0.        ,  0.09090909,  0.18181818,  0.27272727,  0.36363636,
        0.45454545,  0.54545455,  0.63636364,  0.72727273,  0.81818182,
        0.90909091,  1.        ]

beadsize=20

fastids=IMP.pmi.tools.get_ids_from_fasta_file(fastafile)



m=IMP.Model()
simo = IMP.pmi.representation.SimplifiedModel(m) 

for n in range(len(components)):
    simo.add_component_name(components[n],color=colors[n])
    simo.add_component_sequence(components[n],fastafile,id=fastids[n])
    simo.autobuild_pdb_and_intervening_beads(components[n],pdbfile,chains[n],
                                             resolutions=[1,10],beadsize=beadsize)
    simo.setup_component_sequence_connectivity(components[n],1)
    
ev=IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(simo,resolution=10)
ev.add_to_model()

o = IMP.pmi.output.Output()
o.init_rmf("conformations.rmf",[simo.prot])
o.write_rmf("conformations.rmf")

simo.optimize_floppy_bodies(1000)

o.write_rmf("conformations.rmf")
o.close_rmf("conformations.rmf")
