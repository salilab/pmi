import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.proteomics
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.em
import IMP.pmi.representation
import IMP.pmi.tools

def get_random_gaussian_3d(center):
    std=np.random.random_sample(3,) * 20
    var=[s**2 for s in std]
    rot=IMP.algebra.get_random_rotation_3d()
    trans=IMP.algebra.Transformation3D(center,rot)
    return IMP.algebra.Gaussian3D(IMP.algebra.ReferenceFrame3D(trans),var)

class TestEMRestraint(IMP.test.TestCase):
    def setUp(self):
         IMP.test.TestCase.setUp(self)

        self.m = IMP.Model()
        self.p0 = IMP.Particle(self.m)
        self.g0 = IMP.core.Gaussian.setup_particle(self.p0,get_random_gaussian_3d([0,0,0]))
        IMP.atom.Mass.setup_particle(self.p0,np.random.rand()*10)
        self.p1 = IMP.Particle(self.m)
        self.g1 = IMP.core.Gaussian.setup_particle(self.p1,get_random_gaussian_3d([0,0,0]))
        IMP.atom.Mass.setup_particle(self.p1,np.random.rand()*10)

    def test_move_to_center(self):
        gem=IMP.pmi.restraints.em(IMP.core.Hierarchy.setup_particle(self.p0),
                                  target_ps=[self.p1],
                                  cutoff_dist_for_container=0.0,
                                  target_radii_scale=10.0,
                                  model_radii_scale=10.0)
        gem.add_to_model()

        s0=self.m.evaluate(False)
        print 'init score',s0

        trans=IMP.algebra.Transformation3D(np.random.random_sample(3,)*5)
        IMP.atom.transform(IMP.atom.RigidBody(self.p0),trans)
        s1=self.m.evaluate(False)
        print 'transforming with',trans
        print 'score after translating',s1

        gem.center_model_on_target_density()
        s2=self.m.evaluate(False)
        print 'score after re-centering',s2

        self.assertAlmostEqual(s0,s2)

if __name__ == '__main__':
    IMP.test.main()
