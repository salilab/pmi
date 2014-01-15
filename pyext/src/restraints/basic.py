#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

class ExternalBarrier():

    def __init__(self, representation, radius,hierarchies=None,resolution=None):
        self.m = representation.prot.get_model()
        self.rs = IMP.RestraintSet(self.m, 'barrier')

        self.radius = radius
        self.label = "None"

        c3 = IMP.algebra.Vector3D(0, 0, 0)
        ub3 = IMP.core.HarmonicUpperBound(radius, 10.0)
        ss3 = IMP.core.DistanceToSingletonScore(ub3, c3)
        lsc = IMP.container.ListSingletonContainer(self.m)
        # IMP.atom.get_by_type
        particles=IMP.pmi.tools.select(representation,resolution=resolution,hierarchies=hierarchies)
        lsc.add_particles(particles)
        r3 = IMP.container.SingletonsRestraint(ss3, lsc)
        self.rs.add_restraint(r3)

    def set_label(self, label):
        self.label = label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_restraint(self):
        return self.rs

    def get_output(self):
        self.m.update()
        output = {}
        score = self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["ExternalBarrier_" + self.label] = str(score)
        return output