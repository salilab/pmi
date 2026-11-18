import IMP
import IMP.test
import IMP.core
import IMP.container
import IMP.pmi
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.restraints
import IMP.pmi.restraints.crosslinking
from math import log, exp, pi
import numpy as np


def lens_volume(a, b, r):
    """Volume of the intersection of balls of radii a and b, r apart"""
    rmin = np.minimum(a, b)
    rs = np.where(r > 0., r, 1.)
    mid = (pi * (a + b - r) ** 2 * (r * r + 2. * r * (a + b)
                                     - 3. * (a - b) ** 2) / (12. * rs))
    return np.where(r >= a + b, 0.,
                    np.where(r <= abs(a - b), 4. / 3. * pi * rmin ** 3, mid))


def sphere_area_within(r, d, radius):
    """Area of the sphere of radius r, centered at the origin, that lies
       inside the ball of the given radius centered at distance d"""
    mid = pi * r * (radius * radius - (r - d) ** 2) / d
    return np.where(r + d <= radius, 4. * pi * r * r,
                    np.where((r <= d - radius) | (r >= d + radius),
                             0., mid))


def forward_model(dist, sigmai, sigmaj, length):
    """Probability that two points, uniform in spheres of radii sigmai and
       sigmaj whose centers are dist apart, are closer than length
       (integrated over the distance s of the first point from the center
       of the second sphere)"""
    smax = dist + sigmai
    points = sorted(set([0., smax] + [x for x in (abs(dist - sigmai),
                                                  abs(sigmaj - length),
                                                  sigmaj + length)
                                      if 0. < x < smax]))
    total = 0.
    for lo, hi in zip(points[:-1], points[1:]):
        # Composite Simpson's rule on each smooth piece
        s = np.linspace(lo, hi, 81)
        f = (sphere_area_within(s, dist, sigmai)
             * lens_volume(sigmaj, length, s))
        h = (hi - lo) / 80.
        total += h / 3. * (f[0] + f[-1] + 4. * f[1:-1:2].sum()
                           + 2. * f[2:-1:2].sum())
    voli = 4.0 / 3.0 * pi * sigmai ** 3
    volj = 4.0 / 3.0 * pi * sigmaj ** 3
    return total / voli / volj


def get_probability(xyz1s,xyz2s,sigma1s,sigma2s,psis,length,slope):
    onemprob = 1.0

    for n in range(len(xyz1s)):
        xyz1=xyz1s[n]
        xyz2=xyz2s[n]
        sigma1=sigma1s[n]
        sigma2=sigma2s[n]
        psi = psis[n]
        psi = psi.get_scale()
        dist=IMP.core.get_distance(xyz1, xyz2)

        sigmai = sigma1.get_scale()
        sigmaj = sigma2.get_scale()
        pofr = forward_model(max(dist, 0.0001), sigmai, sigmaj, length)

        factor = (1.0 - (psi * (1.0 - pofr) + pofr * (1 - psi))*exp(-slope*dist))
        onemprob = onemprob * factor
    prob = 1.0 - onemprob
    return prob


def log_evaluate(restraints):
    prob = 1.0
    score = 0.0

    for r in restraints:
        prob *= r.unprotected_evaluate(None)
        if prob<=0.0000000001:
            score=score-log(prob)
            prob=1.0

    score=score-log(prob)
    return score


class Tests(IMP.test.TestCase):

    def test_restraint_probability_beads(self):
        m = IMP.Model()
        simo = IMP.pmi.topology.System(m)
        st1 = simo.create_state()
        for molname in ("ProtA", "ProtB"):
            mol = st1.create_molecule(molname, "A" * 10)
            mol.add_representation(mol.get_non_atomic_residues(),
                                   resolutions=[10])
        root_hier = simo.build()

        xl = IMP.pmi.restraints.crosslinking.DisulfideCrossLinkRestraint(
            root_hier,
            (1,1,"ProtA"),
            (1,1,"ProtB"),
            label="DisulfideBond1",
            resolution=1,
            slope=0.01)
        xl.add_to_model()

        p1=xl.xl["Particle1"]
        p2=xl.xl["Particle2"]
        sigma=xl.xl["Sigma"]
        psi=xl.xl["Psi"]

        xyz1=IMP.core.XYZ(p1)
        xyz2=IMP.core.XYZ(p2)

        for d in range(100):
            xyz2.set_coordinates((0,d,0))
            test_prob=get_probability([xyz1],[xyz2],[sigma],[sigma],[psi],6.5,0.01)
            self.assertAlmostEqual(xl.rs.unprotected_evaluate(None),-log(test_prob),delta=0.00001)


if __name__ == '__main__':
    IMP.test.main()
