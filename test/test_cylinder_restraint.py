from __future__ import print_function, division
import IMP
import IMP.atom
import IMP.pmi
import IMP.pmi.topology
import IMP.test
import IMP.pmi.restraints.basic
import math


class Tests(IMP.test.TestCase):


    def test_values(self):
        mdl = IMP.Model()
        p=IMP.Particle(mdl)
        d=IMP.core.XYZR.setup_particle(p)
        d.set_coordinates((0,0,0))
        d.set_radius(1.0)
        IMP.atom.Mass.setup_particle(p,1.0)
        h=IMP.atom.Hierarchy.setup_particle(p)
        cr=IMP.pmi.restraints.basic.CylinderRestraint(mdl,[h],10,20)
        cr.set_was_used(True)
        for r in range(100):
            d.set_coordinates((r,0,0))
            cr.unprotected_evaluate(None)


    def test_angles(self):
        mdl = IMP.Model()
        p=IMP.Particle(mdl)
        d=IMP.core.XYZR.setup_particle(p)
        d.set_coordinates((0,0,0))
        d.set_radius(1.0)
        IMP.atom.Mass.setup_particle(p,1.0)
        h=IMP.atom.Hierarchy.setup_particle(p)
        cr=IMP.pmi.restraints.basic.CylinderRestraint(mdl,[h],10,20,-72,72)
        cr.set_was_used(True)
        for angle in range(360):
            anglerad=float(angle)/180.0*math.pi
            d.set_coordinates((math.cos(anglerad),math.sin(anglerad),0))

if __name__ == '__main__':
    IMP.test.main()
