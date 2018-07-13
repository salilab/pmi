import math
import IMP
import IMP.core
import IMP.pmi1.topology
import IMP.pmi1.restraints.basic
import IMP.test


def _harmonic_prob(x, x0, sig):
    return math.exp(-(x-x0)**2 / 2. / sig**2)


class Tests(IMP.test.TestCase):
    def test_bistable_distance(self):
        m = IMP.Model()
        p1 = IMP.core.XYZ.setup_particle(IMP.Particle(m))
        p2 = IMP.core.XYZ.setup_particle(IMP.Particle(m))
        p1.set_coordinates([0., 0., 0.])
        p2.set_coordinates([0., 0., 2.])
        dists = [1., 4.]
        weights = [.5, .5]
        sigmas = [1., 1.]
        r = IMP.pmi1.restraints.basic.BiStableDistanceRestraint(
            m, p1, p2, dists[0], dists[1], sigmas[0], sigmas[1], weights[0],
            weights[1])
        r.set_was_used(True)
        rscore = r.unprotected_evaluate(None)

        tprob = (.5 * _harmonic_prob(2., dists[0], sigmas[0]) +
                 .5 * _harmonic_prob(2., dists[1], sigmas[1]))
        tscore = -math.log(tprob)
        self.assertAlmostEqual(rscore, tscore, delta=1e-6)
        self.assertEqual(len(r.do_get_inputs()), 2)
        self.assertListEqual(r.do_get_inputs(), [p1, p2])
        self.assertRaises(
            ValueError, IMP.pmi1.restraints.basic.BiStableDistanceRestraint, m,
            p1, p2, dists[0], dists[1], sigmas[0], sigmas[1], weights[0],
            weights[1] + 1e-5)


if __name__ == '__main__':
    IMP.test.main()
