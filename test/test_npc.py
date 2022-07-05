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


if __name__ == '__main__':
    IMP.test.main()
