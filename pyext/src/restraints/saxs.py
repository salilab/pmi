#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container


class SAXSISDRestraint():

    import IMP.saxs 
    import IMP.isd 
    import IMP.isd2 
    import IMP.pmi.tools

    def __init__(self, representation,profile,resolution=0,weight=1,
            ff_type=IMP.saxs.HEAVY_ATOMS):

        from numpy import array,eye
        self.m = representation.prot.get_model()
        self.label = "None"
        self.rs = IMP.RestraintSet(self.m, 'saxs')

        self.sigmamaxtrans = 0.05
        self.gammamaxtrans = 0.05
        self.prof = IMP.saxs.Profile(profile)

        atoms=IMP.pmi.tools.select(representation,resolution=resolution)

        # sigma nuisance
        self.sigma = IMP.pmi.tools.SetupNuisance(
            self.m,
            10.0,
            0.00001,
            100,
            True).get_particle(
        )

        # gamma nuisance, initial value is ML estimate with diagonal covariance
        print "create profile"
        self.th = IMP.saxs.Profile(self.prof.get_min_q(),
                                            self.prof.get_max_q(), self.prof.get_delta_q())
        print "calculate profile"
        
        print len(atoms)
        
        print ff_type
        
        #for a in atoms: print IMP.atom.Residue(a).get_residue_type()
        
        self.th.calculate_profile(atoms, ff_type)
        print "setup gamma"
        gammahat = array([self.prof.get_intensity(i) / self.th.get_intensity(i)
                          for i in xrange(self.prof.size() - 1)]).mean()
        self.gamma = IMP.pmi.tools.SetupNuisance(
                self.m, gammahat, 1e-12, 1e8, True).get_particle()

        self.w = IMP.pmi.tools.SetupWeight(self.m).get_particle()

        # take identity covariance matrix for the start
        self.cov = eye(self.prof.size()).tolist()

        print "create restraint"
        self.saxs = IMP.isd2.SAXSRestraint(self.prof, self.sigma,
                                           self.gamma, self.w)
        self.saxs.add_scatterer(atoms, self.cov, ff_type)

        print "done"
        self.rs.add_restraint(self.saxs)
        self.rs.set_weight(weight)

        # self.saxs_stuff={'nuis':(sigma,gamma),'cov':cov,
        #        'exp':prof,'th':tmp}

        self.rs2 = IMP.RestraintSet(self.m, 'jeffreys')
        j2 = IMP.isd.JeffreysRestraint(self.m, self.gamma)
        self.rs2.add_restraint(j2)

    def update_covariance_matrices(self, tau):
        import numpy
        tau = 0.1
        self.th.calculate_varianced_profile(
            self.atoms,
            impsaxs.HEAVY_ATOMS,
            tau)
        prof = numpy.zeros(self.th.size())
        Sigma = numpy.zeros((self.th.size(), self.th.size()))
        # absolute variance matrix
        for i in xrange(self.th.size()):
            prof[i] = self.th.get_intensity(i)
            for j in xrange(i, self.th.size()):
                Sigma[i, j] = self.th.get_variance(i, j)
                if j != i:
                    Sigma[j, i] = Sigma[i, j]
        # relative matrix
        for i in xrange(prof.shape):
            for j in xrange(prof.shape):
                Sigma[i, j] = Sigma[i, j] / (prof[i] * prof[j])
        self.cov = Sigma.tolist()
        self.saxs.set_cov(0, self.cov)

    def get_gamma_value(self):
        return self.gamma.get_scale()

    def set_label(self, label):
        self.label = label

    def add_to_model(self):
        self.m.add_restraint(self.rs)
        self.m.add_restraint(self.rs2)

    def set_gammamaxtrans(self, gammamaxtrans):
        self.gammamaxtrans = gammamaxtrans

    def set_sigmamaxtrans(self, sigmamaxtrans):
        self.sigmamaxtrans = sigmamaxtrans

    def get_particles_to_sample(self):
        ps = {}
        ps["Nuisances_SAXSISDRestraint_Sigma_" +
            self.label] = ([self.sigma], self.sigmamaxtrans)
        ps["Nuisances_SAXSISDRestraint_Gamma_" +
            self.label] = ([self.gamma], self.gammamaxtrans)
        return ps

    def get_output(self):
        self.m.update()
        output = {}
        score = self.rs.unprotected_evaluate(None)
        score2 = self.rs2.unprotected_evaluate(None)
        output["_TotalScore"] = str(score + score2)

        output["SAXSISDRestraint_Likelihood_" + self.label] = str(score)
        output["SAXSISDRestraint_Prior_" + self.label] = str(score2)
        output["SAXSISDRestraint_Sigma_" +
               self.label] = str(self.sigma.get_scale())
        output["SAXSISDRestraint_Gamma_" +
               self.label] = str(self.gamma.get_scale())
        return output
