#!/usr/bin/env python
"""useful tools for setting up sampling"""

import IMP
import IMP.atom

class SampleObjects(object):
    """ hack class to provide things to sample for PMI::samplers """
    def __init__(self,dict_name,pack_in_dict):
        self.d={dict_name:pack_in_dict}
    def get_particles_to_sample(self):
        return self.d

def enable_md_sampling(mdl,
                       hier=None,
                       particles=None,
                       include_siblings=False,
                       exclude_backbone=False):
    """ Adds necessary attributes to the selected residues for MD sampling
    @param model            The IMP model
    @param hier             Hierarchy to sample
    @param particles        Particles to sample
    @param selection        Single or multiple selections for enabling sampling
    @param include_siblings Get the siblings of the passed particles and sample them too
    @param exclude_backbone Don't sample backbone atoms
    """
    vxkey = IMP.FloatKey('vx')
    vykey = IMP.FloatKey('vy')
    vzkey = IMP.FloatKey('vz')
    backbone = [IMP.atom.AT_C,IMP.atom.AT_N, IMP.atom.AT_CA]
    if particles is None:
        particles=[]
    if hier is not None:
        particles+=IMP.core.get_leaves(hier)
    all_ps=[]
    for p in particles:
        if include_siblings:
            ps=[x.get_particle() for x in IMP.atom.Hierarchy(p).get_parent().get_children()]
        else:
            ps=[p]
        for pp in ps:
            if exclude_backbone and IMP.atom.Atom(mdl,pp.get_index()).get_atom_type() in backbone:
                continue
            IMP.core.XYZ(mdl,pp.get_index()).set_coordinates_are_optimized(True)
            mdl.add_attribute(vxkey,pp.get_index(),0.0)
            mdl.add_attribute(vykey,pp.get_index(),0.0)
            mdl.add_attribute(vzkey,pp.get_index(),0.0)
            all_ps.append(pp)
    return [SampleObjects('Floppy_Bodies_SimplifiedModel',[all_ps])]
