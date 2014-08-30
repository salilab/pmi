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

def enable_md_sampling(hier=None,
                       particles=None,
                       include_siblings=True):
    """ Adds necessary attributes to the selected residues for MD sampling
    @param root             The root node
    @param particles        Particles to sample
    @param selection        Single or multiple selections for enabling sampling
    @param include_siblings Get the siblings of the passed particles and sample them too """
    vxkey = IMP.FloatKey('vx')
    vykey = IMP.FloatKey('vy')
    vzkey = IMP.FloatKey('vz')
    if particles is None:
        particles=[]
    if hier is not None:
        particles+=IMP.core.get_leaves(hier)
    all_ps=[]
    for p in particles:
        if include_siblings:
            ps=IMP.atom.Hierarchy(p).get_parent().get_children()
        else:
            ps=[p]
        all_ps+=ps
        for pp in ps:
            IMP.core.XYZ(pp).set_coordinates_are_optimized(True)
            pp.add_attribute(vxkey, 0.0)
            pp.add_attribute(vykey, 0.0)
            pp.add_attribute(vzkey, 0.0)
    return [SampleObjects('Floppy_Bodies_SimplifiedModel',[all_ps])]
