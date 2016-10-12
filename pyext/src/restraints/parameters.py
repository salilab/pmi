#!/usr/bin/env python
import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.isd
import IMP.pmi.tools
import IMP.pmi.restraints


class WeightRestraint(IMP.pmi.restraints.RestraintBase):

    def __init__(self, w, lower, upper, kappa, label=None, weight=1.):
        self.w = w
        m = self.w.get_model()
        super(WeightRestraint, self).__init__(m, label=label, weight=weight)
        self.lower = lower
        self.upper = upper
        self.kappa = kappa
        self.rs.add_restraint(
            IMP.isd.WeightRestraint(
                self.w,
                self.lower,
                self.upper,
                self.kappa))


class JeffreysPrior(IMP.pmi.restraints.RestraintBase):

    def __init__(self, nuisance, label=None, weight=1.):
        m = nuisance.get_model()
        super(JeffreysPrior, self).__init__(m, label=label, weight=weight)
        jp = IMP.isd.JeffreysRestraint(self.m, nuisance)
        self.rs.add_restraint(jp)
