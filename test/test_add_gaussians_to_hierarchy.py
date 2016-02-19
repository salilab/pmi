import os
import IMP
import IMP.test

import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.output

class Tests(IMP.test.TestCase):
    def test_add_em_gmms_to_hierarchy(self):
        """Test adding EM Restraint GMMs to PMI2 Hierarchies and RMF"""

        mdl = IMP.Model()
        s = IMP.pmi.topology.System(mdl)
        st1 = s.create_state()

        seqs = IMP.pmi.topology.Sequences(self.get_input_file_name('seqs.fasta'))

        m1 = st1.create_molecule("Prot1",sequence=seqs["Protein_1"])


        m1.add_representation(m1.get_residues(),resolutions=[1], setup_particles_as_densities=True)

        hier = m1.build()

        densities = IMP.atom.Selection(hier,representation_type=IMP.atom.DENSITIES).get_selected_particles() +\
                    [r.get_hierarchy() for r in m1.get_non_atomic_residues()]

        gem = IMP.pmi.restraints.em.GaussianEMRestraint(densities,
                                                        target_fn=self.get_input_file_name('prot_gmm.txt'),
                                                        target_is_rigid_body=True)

        gem.set_label("em_1")
        gem.add_to_model()

        gem.add_target_density_to_hierarchy(st1.get_hierarchy())

        # Add a second gmm, which should become a second chain, B

        gem2 = IMP.pmi.restraints.em.GaussianEMRestraint(densities,
                                                        target_fn=self.get_input_file_name('prot_gmm.txt'),
                                                        target_is_rigid_body=True)

        gem2.set_label("em_2")
        gem2.add_to_model()

        gem2.add_target_density_to_hierarchy(st1.get_hierarchy())

        # Test that a second child molecule was added to State
        self.assertEqual(len(st1.get_hierarchy().get_children()), 2)

        # Test that there is a molecule named EM_maps as a child of State
        self.assertTrue("EM_maps" in [m.get_name() for m in st1.get_hierarchy().get_children()])

        # Test that two child chains were added to the molecule EM_maps
        n_maps = 0
        for m in st1.get_hierarchy().get_children():
            if m.get_name() == 'EM_maps':
                n_maps = len(m.get_children())
        self.assertEqual(n_maps, 2)


        output = IMP.pmi.output.Output()

        output.init_pdb("test_psf_writing.pdb", s.get_hierarchy())
        #output.write_pdb("test_psf_writing.pdb")
        output.write_psf("test_psf_writing.psf","test_psf_writing.pdb")

        os.unlink('test_psf_writing.psf')
        os.unlink('test_psf_writing.pdb')

if __name__ == '__main__':
    IMP.test.main()
