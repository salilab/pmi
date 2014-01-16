#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container


class GaussianEMRestraint():

    def __init__(self, densities, target_fn):
        global sys, impisd2, tools
        import sys
        import IMP.isd2
        import IMP.isd2.gmm_tools
        import IMP.pmi.tools as tools

        # some parameters
        self.label="None"
        self.sigmaissampled = True
        self.sigmamaxtrans = 0.3
        self.sigmamin = 1.0
        self.sigmamax = 100.0
        self.sigmainit = 2.0
        self.cutoff_dist_for_container = 10.0
        self.tabexp = True
        self.label="None"

        # setup target GMM
        self.m = densities[0].get_model()
        target_ps = []
        IMP.isd2.gmm_tools.read_gmm_txt(target_ps, target_fn, self.m)

        # model GMM
        model_ps = []
        for h in densities:
            model_ps += IMP.core.get_leaves(h)

        # sigma particle
        self.sigmaglobal = tools.SetupNuisance(self.m, self.sigmainit,
                                               self.sigmamin, self.sigmamax,
                                               self.sigmaissampled).get_particle()

        # create restraint
        print 'target', len(target_ps), 'model', len(model_ps)
        self.gaussianEM_restraint = IMP.isd2.GaussianEMRestraint(
            model_ps, target_ps, self.sigmaglobal,
            self.cutoff_dist_for_container, False, False)
        print 'done setup'
        self.rs = IMP.RestraintSet(self.m, 'GaussianEMRestraint')
        self.rs.add_restraint(self.gaussianEM_restraint)

    def set_label(self, label):
        self.label = label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_particles_to_sample(self):
        ps = {}
        if self.sigmaissampled:
            ps["Nuisances_GaussianEMRestraint_sigma_" +
                self.label] = ([self.sigmaglobal], self.sigmamaxtrans)
        return ps

    def get_hierarchy(self):
        return self.prot

    def get_restraint_set(self):
        return self.rs

    def get_output(self):
        self.m.update()
        output = {}
        score = self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["GaussianEMRestraint_" +
               self.label] = str(self.rs.unprotected_evaluate(None))
        output["GaussianEMRestraint_sigma_" +
               self.label] = str(self.sigmaglobal.get_scale())
        return output
