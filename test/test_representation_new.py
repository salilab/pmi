import IMP
import IMP.atom
import IMP.pmi
import IMP.pmi.representation_new as r
import IMP.test

class RepresentationNewTest(IMP.test.TestCase):

    def test_read_sequences(self):
        # test without name map
        seqs0=r.Sequences(self.get_input_file_name('seqs.fasta'))
        self.assertEqual(len(seqs0),3)
        self.assertEqual(seqs0['Protein_1'],'MEIKEVDDRA')
        self.assertEqual(seqs0['Protein_2'],'MELEPTLFGIIEA')
        self.assertEqual(seqs0['Protein_3'],'LAPQLLSQSHLQTFVSDV')

        # test with name map
        seqs=r.Sequences(self.get_input_file_name('seqs.fasta'),
                         name_map={'Protein_1':'Prot1',
                                   'Protein_2':'Prot2',
                                   'Protein_3':'Prot3'})
        self.assertEqual(len(seqs),3)
        self.assertEqual(seqs['Prot1'],'MEIKEVDDRA')
        self.assertEqual(seqs['Prot2'],'MELEPTLFGIIEA')
        self.assertEqual(seqs['Prot3'],'LAPQLLSQSHLQTFVSDV')

    def test_system_base(self):
        '''Test systembase functions like create hierarchy and create child'''
        sb=r._SystemBase()
        root=sb._create_hierarchy()
        child=sb._create_child(root)
        self.assertEqual(child.get_parent(),root)
        sb.build()

    def test_create_states(self):
        '''Test State-creation from System'''
        s=r.System()
        for i in range(10):
            self.assertEqual(s.get_number_of_states(),i)
            st=s.create_state()
            self.assertEqual(st.state.get_parent(),s.system)
            self.assertEqual(st.mdl,s.mdl)
        self.assertEqual(s.get_number_of_states(),10)


    def test_create_molecules(self):
        '''Test Molecule creation from State'''
        s=r.System()
        seqs=r.Sequences(self.get_input_file_name('seqs.fasta'),
                         name_map={'Protein_1':'Prot1',
                                   'Protein_2':'Prot2',
                                   'Protein_3':'Prot3'})

        # create state 1 with 2 molecules
        st=s.create_state()
        p1=st.create_molecule("Prot1",sequence=seqs["Prot1"])
        p2=st.create_molecule("Prot2",sequence=seqs["Prot2"])
        self.assertEqual(p1.molecule.get_parent(),st.state)
        self.assertEqual(p2.molecule.get_parent(),st.state)
        self.assertEqual(p1.mdl,st.mdl)
        self.assertEqual(p2.mdl,st.mdl)
        self.assertEqual(p1.name,"Prot1")
        self.assertEqual(p2.name,"Prot2")
        self.assertEqual(p1.number_of_copies,1)
        self.assertEqual(len(st.state.get_children()),2)

        # create state 2 with one molecule
        st2=s.create_state()
        p3=st2.create_molecule("Prot3",sequence=seqs["Prot3"])
        self.assertEqual(p3.molecule.get_parent(),st2.state)
        self.assertEqual(p3.mdl,st2.mdl)
        self.assertEqual(p3.name,"Prot3")
        self.assertEqual(p3.number_of_copies,1)
        self.assertEqual(len(st2.state.get_children()),1)


    def test_create_copies(self):
        '''Test creation of Copies (does NOT check Copy content)'''
        s=r.System()
        seqs=r.Sequences(self.get_input_file_name('seqs.fasta'),
                         name_map={'Protein_1':'Prot1',
                                   'Protein_2':'Prot2',
                                   'Protein_3':'Prot3'})

        # create state 1 with 2 molecules and several copies
        st1=s.create_state()
        p1=st1.create_molecule("Prot1",sequence=seqs["Prot1"])
        p2=st1.create_molecule("Prot2",sequence=seqs["Prot2"])
        p1.add_copy()
        p1.add_copy()
        p1.add_copy()
        p2.add_copy()

        # create state 2 with one molecule
        st2=s.create_state()
        p3=st2.create_molecule("Prot3",sequence=seqs["Prot3"])

        # build system
        s.build()

        # check number of molecules per state
        self.assertEqual(p1.number_of_copies+p2.number_of_copies,len(st1.state.get_children()))
        self.assertEqual(p3.number_of_copies,len(st2.state.get_children()))

    '''
    def test_add_structure(self):
        # incomplete
        s=r.System()
        # create a new state
        st1=s.create_state()
        seqs=r.Sequences(fasta_fn="my_fasta_file.fasta")
        p1=st1.create_molecule("Prot1",sequence=seqs["Prot1"])
        p2=st1.create_molecule("Prot2",sequence=seqs["Prot2"])
        st2=s.create_state()
        p3=st2.create_molecule("Prot3",sequence=seqs["Prot3"])
        p1.add_copy()
        p1.add_structure()
        p2.add_structure()
        p3.add_structure()

    def test_set_representation(self):
        # incomplete
        s=r.System()
        seqs=r.Sequences(fasta_fn="my_fasta_file.fasta")
        # create a new state
        st1=s.create_state()
        p1=st1.create_molecule("Prot1",sequence=seqs["Prot1"])
        p2=st1.create_molecule("Prot2",sequence=seqs["Prot2"])
        st2=s.create_state()
        p3=st2.create_molecule("Prot3",sequence=seqs["Prot3"])
        p1.add_copy()
        p1.add_structure()
        p2.add_structure()
        p3.add_structure()
        p1.set_representation()
        p2.set_representation()
        p3.set_representation()

    def test_build_system(self):
        # incomplete
        s=r.System()
        seqs=r.Sequences(fasta_fn="my_fasta_file.fasta")

        # create a new state
        st1=s.create_state()
        p1=st1.create_molecule("Prot1",sequence=seqs["Prot1"])
        p2=st1.create_molecule("Prot2",sequence=seqs["Prot2"])
        st2=s.create_state()
        p3=st2.create_molecule("Prot3",sequence=seqs["Prot3"])
        p1.add_copy()
        p1.add_structure()
        p2.add_structure()
        p3.add_structure()
        p1.set_representation()
        p2.set_representation()
        p3.set_representation()
        s.build()

    def test_build_partial1(self):
        # incomplete
        s=r.System()
        sequence_prot1=r.Sequences(fastafile="my_fasta_file.fasta",id_fastafile="Prot1")
        sequence_prot2=r.Sequences(fastafile="my_fasta_file.fasta",id_fastafile="Prot2")
        sequence_prot3=r.Sequences(fastafile="my_fasta_file.fasta",id_fastafile="Prot3")
        # create a new state
        st1=s.create_state()
        p1=st1.create_molecule("Prot1",sequence=sequence_prot1)
        p2=st1.create_molecule("Prot2",sequence=sequence_prot2)
        st2=s.create_state()
        p3=st2.create_molecule("Prot3",sequence=sequence_prot3)
        p1.add_copy()
        p1.add_structure()
        p2.add_structure()
        p3.add_structure()
        p1.set_representation()
        p2.set_representation()
        p3.set_representation()
        st1.build()
        st2.build()

    def test_build_partial2(self):
        # incomplete
        s=r.System()
        sequence_prot1=r.Sequences(fastafile="my_fasta_file.fasta",id_fastafile="Prot1")
        sequence_prot2=r.Sequences(fastafile="my_fasta_file.fasta",id_fastafile="Prot2")
        sequence_prot3=r.Sequences(fastafile="my_fasta_file.fasta",id_fastafile="Prot3")
        # create a new state
        st1=s.create_state()
        p1=st1.create_molecule("Prot1",sequence=sequence_prot1)
        p2=st1.create_molecule("Prot2",sequence=sequence_prot2)
        st2=s.create_state()
        p3=st2.create_molecule("Prot3",sequence=sequence_prot3)
        p1.add_copy()
        p1.add_structure()
        p2.add_structure()
        p3.add_structure()
        p1.set_representation()
        p2.set_representation()
        p3.set_representation()
        p1.build()
        p2.build()
        p3.build()
    '''
if __name__ == '__main__':
    IMP.test.main()
