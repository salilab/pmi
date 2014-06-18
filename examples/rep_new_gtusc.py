## \example pmi/sandbox_example.py

import IMP
import IMP.pmi
import IMP.pmi.representation_new

### setup 1-state system
mdl = IMP.Model()
system = IMP.pmi.representation_new.System(mdl)
state1 = system.create_state()

### read sequences into a little data structure
seqs = IMP.pmi.representation_new.Sequences(fasta_fn="data/tusc_sequences.fasta",
                                            name_map={'GCP2_YEAST':'Spc97',
                                                      'GCP3_YEAST':'Spc98'})


### add molecules to the state
spc97 = state1.create_molecule("Spc97", sequence=seqs["Spc97"])
spc98 = state1.create_molecule("Spc98", sequence=seqs["Spc98"])
spc97.add_copy()

### load structure data for each molecule. returns a list of contiguous fragments
s97_atomic = spc97.add_structure(pdb_fn='data/tusc.pdb',chain='A',res_range=None,offset=None)
s97_nonatomic = spc97[:] - spc97_atomic
s97_a = spc97[5:10]
s97_b = spc97[8:15]


spc97.set_representation(res=10,"Beads")
spc97.set_representation(s97_atomic,resolution,representation_type)
#spc97.set_representation(s97_nonatomic,resolution,resolution_type)



s98_fragments = spc98.add_structure(pdb_fn='data/tusc.pdb',chain='B')


### utility to create different representations for atomic/non-atomic parts
# note, can also do this with spc97.set_representations(fragment,spec)
#  where fragment is spc97.range(X,Y) or spc97.range('X','Y')
#spc97.add_representations(resolutions_for_atomic_regions=[1],
#                          resolutions_for_missing_regions=[10])
