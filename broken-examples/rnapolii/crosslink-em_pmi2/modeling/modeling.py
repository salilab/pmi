## \example running the sampling

import IMP
import IMP.rmf
import IMP.pmi
import IMP.pmi.representation_new
import IMP.pmi.restraints_new.stereochemistry
#import IMP.pmi.restraints_new.atomic_xl
import IMP.pmi.sequence_tools
import IMP.pmi.hierarchy_tools
import IMP.pmi.data_tools
import IMP.pmi.samplers
import IMP.pmi.macros
import RMF
import string

### setup 2-state system
model = IMP.Model()
system = IMP.pmi.representation_new.System(model)

pdb_file  ='../data/1WCM.pdb'
fasta_file="../data/1WCM.fasta.txt"
fasta_id_protein_name_map={'1WCM:A|PDBID|CHAIN|SEQUENCE':'Rpb1',
                           '1WCM:B|PDBID|CHAIN|SEQUENCE':'Rpb2',
                           '1WCM:C|PDBID|CHAIN|SEQUENCE':'Rpb3',
                           '1WCM:D|PDBID|CHAIN|SEQUENCE':'Rpb4',
                           '1WCM:E|PDBID|CHAIN|SEQUENCE':'Rpb5',
                           '1WCM:F|PDBID|CHAIN|SEQUENCE':'Rpb6',
                           '1WCM:G|PDBID|CHAIN|SEQUENCE':'Rpb7',
                           '1WCM:H|PDBID|CHAIN|SEQUENCE':'Rpb8',
                           '1WCM:I|PDBID|CHAIN|SEQUENCE':'Rpb9',
                           '1WCM:J|PDBID|CHAIN|SEQUENCE':'Rpb10',
                           '1WCM:K|PDBID|CHAIN|SEQUENCE':'Rpb11',
                           '1WCM:L|PDBID|CHAIN|SEQUENCE':'Rpb12'}


### read sequences into a little data structure
seqs = IMP.pmi.representation_new.Sequences(fasta_fn=fasta_file,
                                            name_map=fasta_id_protein_name_map)






# setup RNA Pol II complex

state1 = system.create_state()
molecules={}
offset_list=[0,0,0,-3,0,0,0,0,0,0,0,0]
for i in range(1,13):

    name='Rpb'+str(i)

    if name == "Rpb4":
        continue
    chain=string.uppercase[i-1]
    print name, chain
    molecules[name]   = state1.create_molecule(name,  sequence=seqs[name],  chain_id=chain)
    structured_part   = molecules[name].add_structure(pdb_fn=pdb_file,chain_id=chain,offset=offset_list[i-1])
    unstructured_part = molecules[name].get_residues()-structured_part
    molecules[name].add_representation(structured_part,  'balls',[1,10])
    molecules[name].add_representation(unstructured_part,'balls',[10])

hier=system.build()
IMP.atom.show_molecular_hierarchy(hier)


dof=IMP.pmi.DOF()
rb1=dof.create_rigid_body()
rb1.add_rigid_members({'molecule':'Rpb1','residue_indexes':range(1,1140+1)})
rb1.set_nonrigid_members({'molecule':'Rpb1','residue_indexes':range(1,1140+1)} & unstructured_part)
rb1.set_max_rotation(0.01)
rb1.set_max_translation(1.0)
rb1.set_max_non_rigid_member_translation(1.0)

sc1=dof.create_symmetry(symmetry_type)
sc1.add_target_members({'molecule':'Rpb1','residue_indexes':range(1,1140+1)})
sc1.add_copy_members({'molecule':'Rpb1','residue_indexes':range(1,1140+1)})

fl1=dof.create_flexible()
fl1.add_members({'molecule':'Rpb1','resolution':[1],'residue_indexes':range(1141,1200)})
fl1.set_translation(1.0)

srb1=dof.super_rigid_body()
srb1.add_rigid_body(rb1)
srb1.add_flexible(fl1)
srb1.set_max_rotation(0.01)
srb1.set_max_translation(1.0)

nuis1=dof.create_parameter()
nuis1.add_particles(restraint)
nuis1.set_max_translation(1.0)
