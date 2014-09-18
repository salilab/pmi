#!/usr/bin/env python
import IMP
import IMP.atom
import IMP.em
import IMP.isd_emxl
import IMP.isd_emxl.gmm_tools
import IMP.pmi.tools as tools
from math import sqrt

class EMRestraint(object):
    def __init__(self,
                 root,
                 map_fn,
                 resolution,
                 origin=None,
                 voxel_size=None,
                 weight=1.0,
                 label="",
                 selection_dict=None):
        """ create a FitRestraint. can provide rigid bodies instead of individual particles """

        print 'FitRestraint: setup'
        print '\tmap_fn',map_fn
        print '\tresolution',resolution
        print '\tvoxel_size',voxel_size
        print '\torigin',origin
        print '\tweight',weight

        # some parameters
        self.mdl = root.get_model()
        self.label = label
        self.dmap = IMP.em.read_map(map_fn,IMP.em.MRCReaderWriter())
        dh = self.dmap.get_header()
        dh.set_resolution(resolution)
        if voxel_size:
            self.dmap.update_voxel_size(voxel_size)
        if type(origin)==IMP.algebra.Vector3D:
            self.dmap.set_origin(origin)
        elif type(origin)==list:
            self.dmap.set_origin(*origin)
        else:
            print 'FitRestraint did not recognize format of origin'
            exit()
        if selection_dict:
            ps=IMP.atom.Selection(root,**selection_dict).get_selected_particles()
        else:
            ps=IMP.atom.get_leaves(root)
        fr = IMP.em.FitRestraint(ps,self.dmap)
        self.rs = IMP.RestraintSet(self.mdl,weight,"FitRestraint")
        self.rs.add_restraint(fr)
        self.set_weight(weight)

    def set_weight(self,weight):
        self.weight = weight
        self.rs.set_weight(weight)

    def set_label(self, label):
        self.label = label

    def add_to_model(self):
        self.mdl.add_restraint(self.rs)

    def get_restraint_set(self):
        return self.rs

    def get_output(self):
        self.mdl.update()
        output = {}
        score = self.weight * self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["EMRestraint_" + self.label] = str(score)
        return output

    def evaluate(self):
        return self.weight * self.rs.unprotected_evaluate(None)

class GaussianEMRestraint(object):

    def __init__(self,
                 model_ps,
                 target_fn='',
                 cutoff_dist_model_model=0.0,
                 cutoff_dist_model_data=0.0,
                 overlap_threshold=0.0,
                 target_mass_scale=1.0,
                 target_radii_scale=1.0,
                 model_radii_scale=1.0,
                 slope=0.0,
                 spherical_gaussians=False,
                 pointwise_restraint=False,
                 local_mm=False,
                 close_pair_container=None,
                 label=""):

        # some parameters
        self.label = label
        sigmaissampled = False
        sigmamaxtrans = 0.3
        sigmamin = 1.0
        sigmamax = 100.0
        sigmainit = 2.0
        self.weight=1

        # setup target GMM
        self.m = model_ps[0].get_model()
        print 'will scale target mass by', target_mass_scale
        target_ps=[]
        if target_fn != '':
            IMP.isd_emxl.gmm_tools.decorate_gmm_from_text(
                target_fn,
                target_ps,
                self.m,
                radius_scale=target_radii_scale,
                mass_scale=target_mass_scale)
        else:
            print 'Gaussian EM restraint: must provide target density file'
            return

        # setup model GMM
        if model_radii_scale != 1.0:
            for p in model_ps:
                rmax = sqrt(max(IMP.core.Gaussian(p).get_variances())) * \
                    model_radii_scale
                if not IMP.core.XYZR.get_is_setup(p):
                    IMP.core.XYZR.setup_particle(p, rmax)
                else:
                    IMP.core.XYZR(p).set_radius(rmax)

        # sigma particle
        sigmaglobal = tools.SetupNuisance(self.m, sigmainit,
                                               sigmamin, sigmamax,
                                               sigmaissampled).get_particle()

        # create restraint
        print 'target num particles', len(target_ps), \
            'total weight', sum([IMP.atom.Mass(p).get_mass()
                                for p in target_ps])
        print 'model num particles', len(model_ps), \
            'total weight', sum([IMP.atom.Mass(p).get_mass()
                                for p in model_ps])

        update_density=False
        update_model=not spherical_gaussians
        tabexp=False
        if not pointwise_restraint:
            self.gaussianEM_restraint = \
               IMP.isd_emxl.GaussianEMRestraint(self.m,
                                                IMP.get_indexes(model_ps),
                                                IMP.get_indexes(target_ps),
                                                sigmaglobal.get_particle().get_index(),
                                                cutoff_dist_model_model,
                                                cutoff_dist_model_data,
                                                slope,
                                                update_model, update_density, tabexp)
        else:
            print 'USING POINTWISE RESTRAINT'
            print 'update model?',update_model
            self.gaussianEM_restraint = \
               IMP.isd_emxl.PointwiseGaussianEMRestraint(self.m,
                                                IMP.get_indexes(model_ps),
                                                IMP.get_indexes(target_ps),
                                                sigmaglobal.get_particle().get_index(),
                                                cutoff_dist_model_model,
                                                cutoff_dist_model_data,
                                                overlap_threshold,
                                                slope,
                                                update_model, update_density, local_mm,
                                                close_pair_container)

        print 'done EM setup'
        self.rs = IMP.RestraintSet(self.m, 'GaussianEMRestraint')
        self.rs.add_restraint(self.gaussianEM_restraint)

    def set_weight(self,weight):
        self.weight = weight
        self.rs.set_weight(weight)

    def set_label(self, label):
        self.label = label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_restraint_set(self):
        return self.rs

    def get_output(self):
        self.m.update()
        output = {}
        score = self.weight * self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["GaussianEMRestraint_" +
               self.label] = str(score)
        #output["GaussianEMRestraint_sigma_" +
        #       self.label] = str(sigmaglobal.get_scale())
        return output

    def evaluate(self):
        return self.weight * self.rs.unprotected_evaluate(None)
