import IMP
import IMP.test
import IMP.pmi1.representation
import IMP.pmi1.restraints.npc


class Tests(IMP.test.TestCase):

    def make_representation(self):
        pdbfile = self.get_input_file_name("nonbond.pdb")
        fastafile = self.get_input_file_name("nonbond.fasta")
        fastids = IMP.pmi1.tools.get_ids_from_fasta_file(fastafile)

        m = IMP.Model()
        r = IMP.pmi1.representation.Representation(m)

        r.create_component("Nup84", color=0.)
        r.add_component_sequence("Nup84", fastafile, id=fastids[0])
        r.autobuild_model("Nup84", pdbfile, "A", offset=1,
                          resolutions=[1, 10], missingbeadsize=1)
        return m, r

    def test_xy_radial_pos_restraint(self):
        """Test XYRadialPositionRestraint"""
        m, r = self.make_representation()
        r = IMP.pmi1.restraints.npc.XYRadialPositionRestraint(
            representation=r, protein="Nup84", lower_bound=5., upper_bound=15.)
        r.set_label('Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['XYRadialPositionRestraint_Test']),
            17.57, delta=1e-2)

    def test_xy_radial_pos_lower_restraint(self):
        """Test XYRadialPositionLowerRestraint"""
        m, r = self.make_representation()
        r = IMP.pmi1.restraints.npc.XYRadialPositionLowerRestraint(
            representation=r, protein="Nup84", lower_bound=25.)
        r.set_label('Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['XYRadialPositionLowerRestraint_Test']),
            33.74, delta=1e-2)

    def test_xy_radial_pos_upper_restraint(self):
        """Test XYRadialPositionUpperRestraint"""
        m, r = self.make_representation()
        r = IMP.pmi1.restraints.npc.XYRadialPositionUpperRestraint(
            representation=r, protein="Nup84", upper_bound=10.)
        r.set_label('Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['XYRadialPositionUpperRestraint_Test']),
            84.48, delta=1e-2)

    def test_z_axial_pos_restraint(self):
        """Test ZAxialPositionRestraint"""
        m, r = self.make_representation()
        r = IMP.pmi1.restraints.npc.ZAxialPositionRestraint(
            representation=r, protein="Nup84", lower_bound=-90.,
            upper_bound=-80.)
        r.set_label('Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['ZAxialPositionRestraint_Test']),
            145.88, delta=1e-2)

    def test_z_axial_pos_lower_restraint(self):
        """Test ZAxialPositionLowerRestraint"""
        m, r = self.make_representation()
        r = IMP.pmi1.restraints.npc.ZAxialPositionLowerRestraint(
            representation=r, protein="Nup84", lower_bound=-60.)
        r.set_label('Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['ZAxialPositionLowerRestraint_Test']),
            62.76, delta=1e-2)

    def test_z_axial_pos_upper_restraint(self):
        """Test ZAxialPositionUpperRestraint"""
        m, r = self.make_representation()
        r = IMP.pmi1.restraints.npc.ZAxialPositionUpperRestraint(
            representation=r, protein="Nup84", upper_bound=-70.)
        r.set_label('Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['ZAxialPositionUpperRestraint_Test']),
            4.32, delta=1e-2)

    def test_y_axial_pos_restraint(self):
        """Test YAxialPositionRestraint"""
        m, r = self.make_representation()
        r = IMP.pmi1.restraints.npc.YAxialPositionRestraint(
            representation=r, protein="Nup84", lower_bound=15.,
            upper_bound=20.)
        r.set_label('Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['YAxialPositionRestraint_Test']),
            7.03, delta=1e-2)

    def test_y_axial_pos_lower_restraint(self):
        """Test YAxialPositionLowerRestraint"""
        m, r = self.make_representation()
        r = IMP.pmi1.restraints.npc.YAxialPositionLowerRestraint(
            representation=r, protein="Nup84", lower_bound=17.)
        r.set_label('Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['YAxialPositionLowerRestraint_Test']),
            21.64, delta=1e-2)

    def test_y_axial_pos_upper_restraint(self):
        """Test YAxialPositionUpperRestraint"""
        m, r = self.make_representation()
        r = IMP.pmi1.restraints.npc.YAxialPositionUpperRestraint(
            representation=r, protein="Nup84", upper_bound=10.)
        r.set_label('Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['YAxialPositionUpperRestraint_Test']),
            5.51, delta=1e-2)

    def test_membrane_surface_restraint(self):
        """Test MembraneSurfaceLocationRestraint"""
        m, r = self.make_representation()
        r = IMP.pmi1.restraints.npc.MembraneSurfaceLocationRestraint(
            representation=r, protein=(11, 11, "Nup84"),
            tor_R=50., tor_r=30.)
        r.set_label('Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['MembraneSurfaceLocationRestraint_Test']),
            8646.77, delta=1e-2)

    def test_membrane_surface_cond_restraint(self):
        """Test MembraneSurfaceLocationConditionalRestraint"""
        m, r = self.make_representation()
        r = IMP.pmi1.restraints.npc.MembraneSurfaceLocationConditionalRestraint(
            representation=r, protein1=(11, 11, "Nup84"),
            protein2=(12, 12, "Nup84"), tor_R=55., tor_r=30.)
        r.set_label('Test')
        r.add_to_model()
        out = r.get_output()
        s = out['MembraneSurfaceLocationConditionalRestraint_Test']
        self.assertAlmostEqual(float(s), 9584.33, delta=1e-2)

    def test_membrane_exclusion_restraint(self):
        """Test MembraneExclusionRestraint"""
        m, r = self.make_representation()
        # Set z=0 for C terminus
        s = IMP.atom.Selection(
            r.prot, molecule='Nup84',
            residue_index=11).get_selected_particles()[0]
        self.assertTrue(IMP.core.XYZ.get_is_setup(s))
        xyz = IMP.core.XYZ(s)
        xyz.set_coordinate(2, 0.)

        r = IMP.pmi1.restraints.npc.MembraneExclusionRestraint(
            representation=r, protein=(11, 11, "Nup84"), tor_R=30., tor_r=10.)
        r.set_label('Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['MembraneExclusionRestraint_Test']),
            2661.07, delta=1e-2)


if __name__ == '__main__':
    IMP.test.main()
