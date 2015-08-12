import IMP
import IMP.pmi
import IMP.pmi.dof
import IMP.pmi.topology
import IMP.test


class TestDOF(IMP.test.TestCase):
    def init_topology1(self,mdl):
        s = IMP.pmi.topology.System(mdl)
        st1 = s.create_state()
        seqs = IMP.pmi.topology.Sequences(self.get_input_file_name('seqs.fasta'))

        m1 = st1.create_molecule("Prot1",sequence=seqs["Protein_1"])
        atomic_res = m1.add_structure(self.get_input_file_name('prot.pdb'),
                                      chain_id='A',res_range=(1,10),offset=-54)
        m1.add_representation(atomic_res,resolutions=[0])
        m1.add_representation(resolutions=[1])
        hier = m1.build()
        return m1
    def init_topology3(self,mdl):
        s = IMP.pmi.topology.System(mdl)
        st1 = s.create_state()
        seqs = IMP.pmi.topology.Sequences(self.get_input_file_name('seqs.fasta'))

        m1 = st1.create_molecule("Prot1",sequence=seqs["Protein_1"])
        m2 = st1.create_molecule("Prot2",sequence=seqs["Protein_2"])
        m3 = st1.create_molecule("Prot3",sequence=seqs["Protein_3"])
        a1 = m1.add_structure(self.get_input_file_name('prot.pdb'),
                              chain_id='A',res_range=(1,10),offset=-54)
        a2 = m2.add_structure(self.get_input_file_name('prot.pdb'),
                              chain_id='B',res_range=(1,13),offset=-179)
        a3 = m3.add_structure(self.get_input_file_name('prot.pdb'),
                              chain_id='G',res_range=(1,10),offset=-54)
        m1.add_representation(a1,resolutions=[0])
        m1.add_representation(resolutions=[1])
        m2.add_representation(a2,resolutions=[0])
        m2.add_representation(resolutions=[1])
        m3.add_representation(a3,resolutions=[0])
        m3.add_representation(resolutions=[1])
        hier = s.build()
        return m1,m2,m3

    def test_mc_rigid_body(self):
        """Test creation of rigid body and nonrigid members"""
        mdl = IMP.Model()
        molecule = self.init_topology1(mdl)
        dof = IMP.pmi.dof.xDegreesOfFreedom(mdl)
        rb_movers = dof.create_rigid_body(molecule,
                                          nonrigid_parts = molecule.get_non_atomic_residues())
        mvs = dof.get_movers()
        self.assertEqual(len(mvs),4)

    def test_mc_super_rigid_body(self):
        mdl = IMP.Model()
        mols = self.init_topology3(mdl)
        dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
        rb1_mov = dof.create_rigid_body(mols[0],
                                        nonrigid_parts = mols[0].get_non_atomic_residues())
        rb2_mov = dof.create_rigid_body(mols[1],
                                        nonrigid_parts = mols[1].get_non_atomic_residues())
        rb3_mov = dof.create_rigid_body(mols[2],
                                        nonrigid_parts = mols[2].get_non_atomic_residues())
        srb_mover = dof.create_super_rigid_body(mols,chain_length=1) #<--- pivot points
        ### rbX = dof.create_rigid_body([mols[0],mols[1]]) should fail
        '''
        >>>print dof
        - super-rigid "SRB1"
          - rigid "Mol1" (8 rigid, 3 nonrigid)
          - rigid "Mol2" (8 rigid, 3 nonrigid)
          - rigid "Mol3" (8 rigid, 3 nonrigid)
        '''
        self.assertEqual(len(dof.get_movers()),10)

    def test_constraint_symmetry(self):

        ### create representation
        mdl = IMP.Model()
        s = IMP.pmi.topology.System(mdl)
        st1 = s.create_state()
        seqs = IMP.pmi.topology.Sequences(self.get_input_file_name('seqs.fasta'))

        m1 = st1.create_molecule("Prot1",sequence=seqs["Prot1"])
        a1 = m1.add_structure(self.get_input_file_name('prot.pdb'),
                              chain_id='A',res_range=(1,10),offset=-54)
        m1.add_representation(a1,resolutions=[0])
        m1.add_representation(resolutions=[1])
        m3 = m1.create_clone()

        m2 = st1.create_molecule("Prot1",sequence=seqs["Prot2"])
        a2 = m2.add_structure(self.get_input_file_name('prot.pdb'),
                              chain_id='B',res_range=(1,13),offset=-179)
        m2.add_representation(atomic_res,resolutions=[0])
        m2.add_representation(resolutions=[1])
        m4 = m2.create_clone()
        root = s.build()

        ### create movers and constraints
        dof = IMP.pmi.dof.DegreesOfFreedom()
        rb1_movers = dof.create_rigid_body(m1,
                                           nonrigid_parts = m1.get_non_atomic_residues())
        rb2_movers = dof.create_rigid_body(m2,
                                           nonrigid_parts = m2.get_non_atomic_residues())
        dof.create_rigid_body(m3,
                              nonrigid_parts = m3.get_non_atomic_residues())
        dof.create_rigid_body(m4,
                              nonrigid_parts = m4.get_non_atomic_residues())


        trans = IMP.algebra.Transformation3D(IMP.algebra.Vector3D(0,0,-10))

        dof.create_symmetry([m1,m2],[m3,m4],trans)
        # checks:
        #   check same number of particles
        #   disable ANY movers involving symmetry copies
        #    (later may support moving them WITH references, but requires code to propagate constraint)

        srb = dof.create_super_rigid_body([m1,m2])   # OK
        srb = dof.create_super_rigid_body([m3,m4]) # should raise exception


    #def test_mc_with_densities(self):
    #    pass

    '''
    def test_mc_flexible_beads(self):
        mdl = IMP.Model()
        molecule = self.init_topology1(mdl)
        hier = molecule.get_hierarchy()
        dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
        sel_nonrigid = IMP.atom.Selection(hier,residue_indexes=[3,4,10])
        fbs = dof.create_flexible_beads(sel_nonrigid,max_trans=1.0)
        mvs = fbs.get_movers()
        self.assertEqual(len(mvs),3)

    def test_multi_select(self):
        mdl = IMP.Model()
        molecule = self.init_topology1(mdl)
        hier = molecule.get_hierarchy()
        res = IMP.pmi.tools.select_at_all_resolutions(hier,residue_index=1)
        self.assertEqual(len(res['BEADS'][0]),9)
        self.assertEqual(len(res['BEADS'][1]),1)

        ps = molecule.get_particles_at_all_resolutions()
        self.assertEqual(len(ps),70)
    '''
    '''
    '''
    '''
    def test_mc_compound_body(self):
        # compund body is a mix of rigid and flexible parts
        # Do we use System???? Do we need a new decorator that says
        # where the structure is coming from????
        # IMP.atom.Source(p)
        # IMP.atom.Source.get_pdb_id()
        # IMP.atom.Source.get_pdb_id()
        # source can be Modeller, PDB, emdb, Coarse-grained, No-source
        # Invent StringKeys
        hierarchy=self.init_topology()
        dof=IMP.pmi.dof.DegreesOfFreedom()
        s=IMP.atom.Selection(hierarchy,molecule=?,resid=range(1,10))
        structured_handle,unstructures_handle=dof.create_compound_body(s)
    '''

    def test_mc_kinematic(self):
        pass

    def test_md(self):
        pass

if __name__ == '__main__':
    IMP.test.main()
