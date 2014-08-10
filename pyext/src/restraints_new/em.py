#!/usr/bin/env python
import IMP
import IMP.atom
import IMP.em

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

        # some parameters
        self.mdl = root.get_model()
        self.label = label
        self.weight=1
        self.dmap = IMP.em.read_map(map_fn,IMP.em.MRCReaderWriter())
        dh = self.dmap.get_header()
        dh.set_resolution(resolution)
        if voxel_size:
            self.dmap.update_voxel_size(voxel_size)
        if type(origin)==IMP.algebra.Vector3D:
            self.dmap.set_origin(origin)
        elif type(origin)==list:
            self.dmap.set_origin(*origin)

        if selection_dict:
            ps=IMP.atom.Selection(root,**selection_dict).get_selected_particles()
        else:
            ps=IMP.atom.get_leaves(root)
        fr = IMP.em.FitRestraint(ps,self.dmap)
        self.rs = IMP.RestraintSet(self.mdl,weight,"FitRestraint")
        self.rs.add_restraint(fr)

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
