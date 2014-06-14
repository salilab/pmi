#A new representation module. It helps to construct the hierarchy
#and deal with multi-state, multi-scale, multi-copies

#Usage Example:

import IMP
import IMP.atom
import IMP.pmi.representationnew as r

# create Sequence handler objects for proteins
sequence_prot1=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot1")
sequence_prot2=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot2")


# initialize the system
s=r.System()

for i in range(10):
  # create a new state
  state=s.create_state()

  # create a new molecule
  p1=state.create_molecule("Prot1",sequence=sequence_prot1)
  
  # add two identical copies to the molecule
  p1.add_copy()
  p1.add_copy()  

  # create another molecule
  p2=state.create_molecule("Prot2",sequence=sequence_prot2)

  # add structures to prot1
  s_1_100=p1.add_structure(pdbfile="mypdb.pdb",    resrange=[1,100])
  s_101_200=p1.add_structure(abinitio_helix=True,resrange=[101,200])
  
  # set the representation
  p1.set_representation(resolutions=[1,10,100])

# build the system
s.build()

# show the hierarchy to debug
IMP.atom.show_molecular_hierarchy(s.system_hierarchy)

