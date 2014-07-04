import IMP
import IMP.atom
import IMP.pmi
import IMP.pmi.representation_new as r
import IMP.test
import IMP.pmi.sequence_tools

def get_atomic_residue_list(residues):
    r1=[]
    for r in residues:
        ps=r.hier.get_children()
        if len(ps)==0:
            r1.append('-')
        else:
            r1.append(IMP.atom.get_one_letter_code(r.get_residue_type()))
    return ''.join(r1)



class RepresentationNewTest(IMP.test.TestCase):
    def test_read_sequences(self):
        '''Test if the sequence reader returns correct strings'''
        # test without name map
        seqs0=r.Sequences(self.get_input_file_name('seqs.fasta'))
        self.assertEqual(len(seqs0),3)
        self.assertEqual(seqs0['Protein_1'],'QEALVVKDLL')
        self.assertEqual(seqs0['Protein_2'],'PEEDILKYVSYTL')
        self.assertEqual(seqs0['Protein_3'],'NVLIGLEGTY')

        # test with name map
        seqs=r.Sequences(self.get_input_file_name('seqs.fasta'),
                         name_map={'Protein_1':'Prot1',
                                   'Protein_2':'Prot2',
                                   'Protein_3':'Prot3'})
        self.assertEqual(len(seqs),3)
        self.assertEqual(seqs['Prot1'],'QEALVVKDLL')
        self.assertEqual(seqs['Prot2'],'PEEDILKYVSYTL')
        self.assertEqual(seqs['Prot3'],'NVLIGLEGTY')

    def test_system_base(self):
        '''Test systembase functions like create hierarchy and create child'''
        sb=r._SystemBase()
        root=sb._create_hierarchy()
        child=sb._create_child(root)
        self.assertEqual(child.get_parent(),root)

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
        m1=st.create_molecule("Prot1",sequence=seqs["Prot1"])
        m2=st.create_molecule("Prot2",sequence=seqs["Prot2"])
        self.assertEqual(m1.molecule.get_parent(),st.state)
        self.assertEqual(m2.molecule.get_parent(),st.state)
        self.assertEqual(m1.mdl,st.mdl)
        self.assertEqual(m2.mdl,st.mdl)
        self.assertEqual(m1.name,"Prot1")
        self.assertEqual(m2.name,"Prot2")
        self.assertEqual(m1.number_of_copies,1)
        self.assertEqual(len(st.state.get_children()),2)

        # create state 2 with one molecule
        st2=s.create_state()
        m3=st2.create_molecule("Prot3",sequence=seqs["Prot3"])
        self.assertEqual(m3.molecule.get_parent(),st2.state)
        self.assertEqual(m3.mdl,st2.mdl)
        self.assertEqual(m3.name,"Prot3")
        self.assertEqual(m3.number_of_copies,1)
        self.assertEqual(len(st2.state.get_children()),1)

        # test if sequences are OK
        self.assertEqual(''.join(r.get_code() for r in m1.residues),seqs["Prot1"])
        self.assertEqual(''.join(r.get_code() for r in m2.residues),seqs["Prot2"])
        self.assertEqual(''.join(r.get_code() for r in m3.residues),seqs["Prot3"])

    """
    def test_create_copies(self):
        '''Test creation of Copies (does NOT check Copy content)'''
        s=r.System()
        seqs=r.Sequences(self.get_input_file_name('seqs.fasta'),
                         name_map={'Protein_1':'Prot1',
                                   'Protein_2':'Prot2',
                                   'Protein_3':'Prot3'})

        # create state 1 with 2 molecules and several copies
        st1=s.create_state()
        m1=st1.create_molecule("Prot1",sequence=seqs["Prot1"])
        m2=st1.create_molecule("Prot2",sequence=seqs["Prot2"])
        m1.add_copy()
        m1.add_copy()
        m1.add_copy()
        m2.add_copy()

        # create state 2 with one molecule
        st2=s.create_state()
        m3=st2.create_molecule("Prot3",sequence=seqs["Prot3"])

        # build system
        s.build()

        # check number of molecules per state
        self.assertEqual(m1.number_of_copies+m2.number_of_copies,len(st1.state.get_children()))
        self.assertEqual(m3.number_of_copies,len(st2.state.get_children()))
    """

    def test_add_structure(self):
        '''Test adding partial structure data to a molecule'''

        s=r.System()
        st1=s.create_state()
        seqs=r.Sequences(self.get_input_file_name('seqs.fasta'),
                         name_map={'Protein_1':'Prot1',
                                   'Protein_2':'Prot2',
                                   'Protein_3':'Prot3'})

        m1=st1.create_molecule("Prot1",sequence=seqs["Prot1"])
        m2=st1.create_molecule("Prot2",sequence=seqs["Prot2"])
        res1=m1.add_structure(self.get_input_file_name('prot.pdb'),chain='A',res_range=(1,10),offset=-54)
        res2=m2.add_structure(self.get_input_file_name('prot.pdb'),chain='B',res_range=(1,13),offset=-179)

        # check that the molecule residues have the right info
        rlist1=get_atomic_residue_list(m1.residues)
        rlist2=get_atomic_residue_list(m2.residues)
        self.assertEqual(rlist1,'QE--VVKDL-')
        self.assertEqual(rlist2,'PEEDILKYVSYTL')

        # check that the returned Residue index sets are correct
        self.assertEqual(res1,set([m1.residues[i] for i in (0,1,4,5,6,7,8)]))
        self.assertEqual(res2,set([m2.residues[i] for i in range(0,13)]))

    def test_residue_access(self):
        '''test functions to retrieve residues'''
        s=r.System()
        st1=s.create_state()
        seqs=r.Sequences(self.get_input_file_name('seqs.fasta'),
                         name_map={'Protein_1':'Prot1',
                                   'Protein_2':'Prot2',
                                   'Protein_3':'Prot3'})

        m1=st1.create_molecule("Prot1",sequence=seqs["Prot1"])
        self.assertEqual(m1[:],set(m1.residues[:]))
        self.assertEqual(m1[1],m1.residues[1])
        self.assertEqual(m1[1:5],set(m1.residues[1:5]))
        self.assertEqual(m1['1'],m1.residues[0])
        self.assertEqual(m1.residue_range(1,5),set(m1.residues[1:5]))
        self.assertEqual(m1.residue_range('2','6'),set(m1.residues[1:5]))
        inv=m1[:]-m1[1:5]
        self.assertEqual(inv,set([m1.residues[0]]+m1.residues[5:10]))

    def test_add_representation(self):
        '''test if add_representations propulates the correct Residues'''
        s=r.System()
        st1=s.create_state()
        seqs=r.Sequences(self.get_input_file_name('seqs.fasta'),
                         name_map={'Protein_1':'Prot1',
                                   'Protein_2':'Prot2',
                                   'Protein_3':'Prot3'})
        m1=st1.create_molecule("Prot1",sequence=seqs["Prot1"])
        atomic_res=m1.add_structure(self.get_input_file_name('prot.pdb'),chain='A',res_range=(1,10),offset=-54)
        m1.add_representation(resolutions=[1])
        m1.add_representation(atomic_res,resolutions=[0])
        for na in (0,1,4,5,6,7,8):
            self.assertEqual(m1[na].representations['balls'],set([0,1]))
        for nna in (2,3,9):
            self.assertEqual(m1[nna].representations['balls'],set([1]))


    def test_build_system(self):
        s=r.System()
        st1=s.create_state()
        seqs=r.Sequences(self.get_input_file_name('seqs.fasta'),
                         name_map={'Protein_1':'Prot1',
                                   'Protein_2':'Prot2',
                                   'Protein_3':'Prot3'})
        m1=st1.create_molecule("Prot1",sequence=seqs["Prot1"])
        atomic_res=m1.add_structure(self.get_input_file_name('prot.pdb'),chain='A',
                                    res_range=(1,10),offset=-54)
        #m1.add_representation(m1[:]-atomic_res,resolutions=[1])
        m1.add_representation(atomic_res,resolutions=[0,1])
        m1.build(merge_type="backbone")
        frags = m1.molecule.get_children()
        #IMP.atom.show_molecular_hierarchy(m1.molecule)
        for f in frags:
            self.assertTrue(IMP.atom.Fragment.get_is_setup(f))
            self.assertTrue(IMP.atom.Representation.get_is_setup(f))
            rep = IMP.atom.Representation(f)
            self.assertEquals(rep.get_resolutions(),[0,1])

            IMP.pmi.structure_tools.show_representation(rep)

        # check if atomic created correctly
        for rnum,rname,anums in zip((1,2,5,6,7,8,9),'QEVVKDL',(9,9,7,7,9,8,8)):
            res = IMP.atom.Selection(m1.molecule,residue_index=rnum,resolution=0).get_selected_particles()
            self.assertEquals(len(res),anums)
            self.assertEquals(IMP.atom.Residue(IMP.atom.Atom(res[0]).get_parent()).get_residue_type(),
                              IMP.pmi.sequence_tools.get_residue_type_from_one_letter_code(rname))
            res1 = IMP.atom.Selection(m1.molecule,residue_index=rnum,resolution=1).get_selected_particles()
            self.assertEquals(len(res1),1)
            self.assertEquals(IMP.atom.Residue(IMP.atom.Atom(res1[0])).get_residue_type(),
                              IMP.pmi.sequence_tools.get_residue_type_from_one_letter_code(rname))

        # check if non-atomic are correct

if __name__ == '__main__':
    IMP.test.main()
