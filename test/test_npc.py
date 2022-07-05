import IMP
import IMP.test
import IMP.pmi.topology
import IMP.pmi.restraints.npc

class Tests(IMP.test.TestCase):

    def make_representation(self):
        pdbfile = self.get_input_file_name("nonbond.pdb")

        m = IMP.Model()
        s = IMP.pmi.topology.System(m)
        state = s.create_state()

        c = state.create_molecule("Nup84", sequence='KF')
        struc = c.add_structure(pdbfile, chain_id="A", offset=-10)
        c.add_representation(struc, resolutions=[1, 10])
        root_hier = s.build()
        return m, root_hier

    def test_xy_radial_pos_restraint(self):
        """Test XYRadialPositionRestraint"""
        m, root_hier = self.make_representation()
        r = IMP.pmi.restraints.npc.XYRadialPositionRestraint(
            hier=root_hier, protein="Nup84", lower_bound=5., upper_bound=15.,
            label='Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['XYRadialPositionRestraint_Score_Test']),
            17.57, delta=1e-2)

    def test_xy_radial_pos_lower_restraint(self):
        """Test XYRadialPositionLowerRestraint"""
        m, root_hier = self.make_representation()
        r = IMP.pmi.restraints.npc.XYRadialPositionLowerRestraint(
            hier=root_hier, protein="Nup84", lower_bound=25.,
            label='Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['XYRadialPositionLowerRestraint_Score_Test']),
            33.74, delta=1e-2)

    def test_xy_radial_pos_upper_restraint(self):
        """Test XYRadialPositionUpperRestraint"""
        m, root_hier = self.make_representation()
        r = IMP.pmi.restraints.npc.XYRadialPositionUpperRestraint(
            hier=root_hier, protein="Nup84", upper_bound=10.,
            label='Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['XYRadialPositionUpperRestraint_Score_Test']),
            84.48, delta=1e-2)

    def test_z_axial_pos_restraint(self):
        """Test ZAxialPositionRestraint"""
        m, root_hier = self.make_representation()
        r = IMP.pmi.restraints.npc.ZAxialPositionRestraint(
            hier=root_hier, protein="Nup84", lower_bound=-90.,
            upper_bound=-80., label='Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['ZAxialPositionRestraint_Score_Test']),
            145.88, delta=1e-2)

    def test_z_axial_pos_lower_restraint(self):
        """Test ZAxialPositionLowerRestraint"""
        m, root_hier = self.make_representation()
        r = IMP.pmi.restraints.npc.ZAxialPositionLowerRestraint(
            hier=root_hier, protein="Nup84", lower_bound=-60., label='Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['ZAxialPositionLowerRestraint_Score_Test']),
            62.76, delta=1e-2)

    def test_z_axial_pos_upper_restraint(self):
        """Test ZAxialPositionUpperRestraint"""
        m, root_hier = self.make_representation()
        r = IMP.pmi.restraints.npc.ZAxialPositionUpperRestraint(
            hier=root_hier, protein="Nup84", upper_bound=-70., label='Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['ZAxialPositionUpperRestraint_Score_Test']),
            4.32, delta=1e-2)

    def test_y_axial_pos_restraint(self):
        """Test YAxialPositionRestraint"""
        m, root_hier = self.make_representation()
        r = IMP.pmi.restraints.npc.YAxialPositionRestraint(
            hier=root_hier, protein="Nup84", lower_bound=15.,
            upper_bound=20., label='Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['YAxialPositionRestraint_Score_Test']),
            7.03, delta=1e-2)

    def test_y_axial_pos_lower_restraint(self):
        """Test YAxialPositionLowerRestraint"""
        m, root_hier = self.make_representation()
        r = IMP.pmi.restraints.npc.YAxialPositionLowerRestraint(
            hier=root_hier, protein="Nup84", lower_bound=17., label='Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['YAxialPositionLowerRestraint_Score_Test']),
            21.64, delta=1e-2)

    def test_y_axial_pos_upper_restraint(self):
        """Test YAxialPositionUpperRestraint"""
        m, root_hier = self.make_representation()
        r = IMP.pmi.restraints.npc.YAxialPositionUpperRestraint(
            hier=root_hier, protein="Nup84", upper_bound=10., label='Test')
        r.add_to_model()
        out = r.get_output()
        self.assertAlmostEqual(
            float(out['YAxialPositionUpperRestraint_Score_Test']),
            5.51, delta=1e-2)


if __name__ == '__main__':
    IMP.test.main()
