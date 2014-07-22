## \example pmi/rep_new_gtusc.py

import IMP
import IMP.pmi
import IMP.pmi.representation_new
import IMP.pmi.restraints_new.stereochemistry
import IMP.pmi.sequence_tools
import IMP.pmi.hierarchy_tools

### setup 1-state system
mdl = IMP.Model()
system = IMP.pmi.representation_new.System(mdl)
stateA = system.create_state()
stateB = system.create_state()


### read sequences into a little data structure
seqs = IMP.pmi.representation_new.Sequences(fasta_fn="gtusc_data/tusc_sequences.fasta",
                                            name_map={'GCP2_YEAST':'Spc97',
                                                      'GCP3_YEAST':'Spc98'})

### add molecules to the state
spc97A = stateA.create_molecule("Spc97", sequence=seqs["Spc97"])
spc97B = stateB.create_molecule("Spc97", sequence=seqs["Spc97"])

### load structure data for each molecule. returns a list of contiguous fragments
s97A_atomic = spc97A.add_structure(pdb_fn='gtusc_data/tusc.pdb',chain='A',res_range=[55,95])
s97B_atomic = spc97B.add_structure(pdb_fn='gtusc_data/tusc.pdb',chain='A',res_range=[55,95])

s97A_nonatomic = spc97A[:] - s97A_atomic

### add representation
spc97A.add_representation(s97A_atomic,'balls',[0])
spc97B.add_representation(s97B_atomic,'balls',[0])

### build system
hier = system.build()
#IMP.atom.show_molecular_hierarchy(hier)
#IMP.pmi.structure_tools.show_representation(IMP.pmi.hierarchy_tools.select_node(hier,0).get_children()[0].get_children()[1])

### add some restraints
charmmA = IMP.pmi.restraints_new.stereochemistry.CharmmForceFieldRestraint(mdl,
                                        IMP.pmi.hierarchy_tools.select_node(hier,0))
charmmA.add_to_model()
print 'one state',mdl.evaluate(False)
charmmB = IMP.pmi.restraints_new.stereochemistry.CharmmForceFieldRestraint(mdl,
                                        IMP.pmi.hierarchy_tools.select_node(hier,1))
charmmB.add_to_model()
print 'two states',mdl.evaluate(False)
