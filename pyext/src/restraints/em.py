#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container


class GaussianEMRestraint():

    def __init__(self, densities, target_fn,cutoff_dist_for_container=10.0,target_mass=1.0,radii_mult_factor=1.0):
        global sys, tools
        import sys
        import IMP.isd_emxl
        import IMP.isd_emxl.gmm_tools
        import IMP.pmi.tools as tools
        from math import sqrt

        # some parameters
        self.label="None"
        self.sigmaissampled = False
        self.sigmamaxtrans = 0.3
        self.sigmamin = 1.0
        self.sigmamax = 100.0
        self.sigmainit = 2.0
        self.tabexp = False
        self.label="None"

        # setup target GMM
        self.m = densities[0].get_model()
        target_ps = []
        IMP.isd_emxl.gmm_tools.decorate_gmm_from_text(target_fn, target_ps, self.m)
        for p in target_ps:
            rmax=sqrt(max(IMP.core.Gaussian(p).get_variances()))*radii_mult_factor
            IMP.core.XYZR.setup_particle(p,rmax)
            mp=IMP.atom.Mass(p)
            mp.set_mass(mp.get_mass()*target_mass)

        # model GMM
        model_ps = []
        for h in densities:
            model_ps += IMP.core.get_leaves(h)

        # sigma particle
        self.sigmaglobal = tools.SetupNuisance(self.m, self.sigmainit,
                                               self.sigmamin, self.sigmamax,
                                               self.sigmaissampled).get_particle()

        # create restraint
        print 'target num particles',len(target_ps),'total weight',sum([IMP.atom.Mass(p).get_mass() for p in target_ps])
        print 'model num particles',len(model_ps),'total weight',sum([IMP.atom.Mass(p).get_mass() for p in model_ps])
        self.gaussianEM_restraint = IMP.isd_emxl.GaussianEMRestraint(self.m,
                                                                     IMP.get_indexes(model_ps),
                                                                     IMP.get_indexes(target_ps),
                                                                     self.sigmaglobal.get_particle().get_index(),
                                                                     cutoff_dist_for_container,
                                                                     False, False)
        print 'done EM setup'
        self.rs = IMP.RestraintSet(self.m, 'GaussianEMRestraint')
        self.rs.add_restraint(self.gaussianEM_restraint)

    def set_weight(self,weight):
        self.rs.set_weight(weight)

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
