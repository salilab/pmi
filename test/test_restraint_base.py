from __future__ import print_function, division
import IMP
import IMP.algebra
import IMP.core
import IMP.pmi
import IMP.pmi.restraints
import IMP.test


class DistanceRestraint(IMP.pmi.restraints.RestraintBase):

    def __init__(self, p1, p2, d, k, name=None, label=None, weight=1.):
        m = p1.get_model()
        super(DistanceRestraint, self).__init__(m, name=name, label=label,
                                                weight=weight)
        f = IMP.core.Harmonic(d, k)
        s = IMP.core.DistancePairScore(f)
        r = IMP.core.PairRestraint(self.m, s, (p1, p2))
        self.rs.add_restraint(r)


class Tests(IMP.test.TestCase):

    def test_setup(self):
        m = IMP.Model()
        p1 = IMP.Particle(m)
        IMP.core.XYZ.setup_particle(p1)
        p2 = IMP.Particle(m)
        IMP.core.XYZ.setup_particle(p2)

        r = DistanceRestraint(p1, p2, 0., 1.)
        self.assertAlmostEqual(r.evaluate(), 0.0, delta=1e-6)
        output = r.get_output()
        self.assertEqual(output["DistanceRestraint_Score"], str(0.0))
        self.assertEqual(output["_TotalScore"], str(0.0))
        self.assertIsInstance(r.get_restraint_set(), IMP.RestraintSet)

    def test_setup_with_label_weight(self):
        m = IMP.Model()
        p1 = IMP.Particle(m)
        IMP.core.XYZ.setup_particle(p1)
        p2 = IMP.Particle(m)
        IMP.core.XYZ.setup_particle(p2)
        IMP.core.XYZ(p2).set_coordinates(IMP.algebra.Vector3D(0, 0, 10))

        r = DistanceRestraint(p1, p2, 0., 1., name="DistanceRestraint2",
                              weight=10., label="Test")
        self.assertAlmostEqual(r.evaluate(), 500, delta=1e-6)
        output = r.get_output()
        self.assertEqual(output["DistanceRestraint2_Score_Test"], str(500.0))
        self.assertEqual(output["_TotalScore"], str(500.0))

        r.set_weight(1.)
        output = r.get_output()
        self.assertEqual(output["DistanceRestraint2_Score_Test"], str(50.0))
        self.assertEqual(output["_TotalScore"], str(50.0))

    def test_set_label_after_init_raises_error(self):
        m = IMP.Model()
        p1 = IMP.Particle(m)
        IMP.core.XYZ.setup_particle(p1)
        p2 = IMP.Particle(m)
        IMP.core.XYZ.setup_particle(p2)

        r = DistanceRestraint(p1, p2, 0., 1., label="Test")
        with self.assertRaises(ValueError):
            r.set_label("Test2")

        r = DistanceRestraint(p1, p2, 0., 1.)
        r.set_label("Test")
        with self.assertRaises(ValueError):
            r.set_label("Test2")

        r = DistanceRestraint(p1, p2, 0., 1.)
        r.get_output()
        with self.assertRaises(ValueError):
            r.set_label("Test")

        r = DistanceRestraint(p1, p2, 0., 1.)
        r.evaluate()
        with self.assertRaises(ValueError):
            r.set_label("Test")

        r = DistanceRestraint(p1, p2, 0., 1.)
        r.add_to_model()
        with self.assertRaises(ValueError):
            r.set_label("Test")

        r = DistanceRestraint(p1, p2, 0., 1.)
        r.get_restraint_set()
        with self.assertRaises(ValueError):
            r.set_label("Test")

        r = DistanceRestraint(p1, p2, 0., 1.)
        r.get_restraint_for_rmf()
        with self.assertRaises(ValueError):
            r.set_label("Test")

        r = DistanceRestraint(p1, p2, 0., 1.)
        r.get_particles_to_sample()
        with self.assertRaises(ValueError):
            r.set_label("Test")


if __name__ == '__main__':
    IMP.test.main()
