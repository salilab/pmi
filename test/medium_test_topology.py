from __future__ import print_function, division
import IMP
import IMP.atom
import IMP.pmi1
import IMP.pmi1.topology
import IMP.pmi1.tools
import IMP.test
import RMF
import IMP.rmf
import IMP.pmi1.macros
import os


def get_atomic_residue_list(residues):
    r1 = []
    for r in residues:
        ps = r.get_hierarchy().get_children()
        if len(ps) == 0:
            r1.append('-')
        else:
            r1.append(IMP.atom.get_one_letter_code(r.get_residue_type()))
    return ''.join(r1)


class Tests(IMP.test.TestCase):

    def test_read_sequences(self):
        """Test if the sequence reader returns correct strings"""
        # test without name map
        seqs0 = IMP.pmi1.topology.Sequences(
            self.get_input_file_name('seqs.fasta'))
        self.assertEqual(len(seqs0), 3)
        self.assertEqual(seqs0['Protein_1'], 'QEALVVKDLL')
        self.assertEqual(seqs0['Protein_2'], 'PEEDILKYVSYTL')
        self.assertEqual(seqs0['Protein_3'], 'QEALVVKDLL')

        # test with name map
        seqs = IMP.pmi1.topology.Sequences(
            self.get_input_file_name('seqs.fasta'),
            name_map={'Protein_1': 'Prot1',
                      'Protein_2': 'Prot2',
                      'Protein_3': 'Prot3'})
        self.assertEqual(len(seqs), 3)
        self.assertEqual(seqs['Prot1'], 'QEALVVKDLL')
        self.assertEqual(seqs['Prot2'], 'PEEDILKYVSYTL')
        self.assertEqual(seqs['Prot3'], 'QEALVVKDLL')

    def test_pmi_molecule_hierarchy(self):
        model=IMP.Model()
        seqs = IMP.pmi1.topology.Sequences(self.get_input_file_name('seqs.fasta'))
        state=IMP.atom.State.setup_particle(IMP.Particle(model),0)
        for seq in seqs:
            mol=IMP.atom.Molecule.setup_particle(IMP.Particle(model))
            ch=IMP.atom.Chain.setup_particle(mol,"A")
            ch.set_sequence(seqs[seq])
            state.add_child(mol)
            mol.set_name(seq)
            IMP.atom.Copy.setup_particle(mol,0)
            pmimol=IMP.pmi1.topology.PMIMoleculeHierarchy(mol)
            self.assertEqual(pmimol.get_sequence(),seqs[seq])
            self.assertEqual(pmimol.get_residue_indexes(),[])

if __name__ == '__main__':
    IMP.test.main()
