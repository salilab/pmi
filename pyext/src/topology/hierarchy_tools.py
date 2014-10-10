#!/usr/bin/env python

"""@namespace IMP.pmi.hierarchy_tools
   Useful tools for finding things in an IMP hierarchy
"""

import IMP
import IMP.atom
from collections import defaultdict

def combine_dicts(a1,a2):
    """ Combine Selection dictionaries. The first one has higher priority.
    Keeps unique keys from each dictionary.
    If keys overlap, but both point to lists, the lists are combined.
    If keys overlap but they don't point to lists, only the first is kept.
    """
    final={}
    if a1 is None:
        return a2
    elif a2 is None:
        return a1
    elif a1 is None and a2 is None:
        return None
    for k in a1:
        final[k]=a1[k]
    for k in a2:
        if k in final:
            if type(final[k]) is not list:
                print "ERROR: you have overlapping keys that aren't lists"
            else:
                final[k]+=a2[k]
        else:
            final[k]=a2[k]
    return final

def get_particles_within_zone(root,ps,radius,
                              entire_residues=False,
                              exclude_backbone=False):
    """Get a list of selections of residues within a zone around some particles
    @param root             The root node to search in
    @param ps               The target particles
    @param radius           The search radius per particle
    @param entire_residues  Expand selection to include entire residues
    @param exclude_backbone Don't select backbone atoms
    \note Currently this only works for atomic problems, and grabs entire residues...
    """
    ret = []
    atoms = IMP.core.get_leaves(root)
    nn = IMP.algebra.NearestNeighbor3D([IMP.core.XYZ(p).get_coordinates() for p in atoms])
    backbone = [IMP.atom.AT_C,IMP.atom.AT_N,IMP.atom.AT_CA]
    for p in ps:
        ns = nn.get_in_ball(IMP.core.XYZ(p).get_coordinates(),radius)
        this_ps=[]
        for n in ns: #looping over each particle near p's position
            if entire_residues:
                r = IMP.atom.get_residue(IMP.atom.Atom(atoms[n]))
                this_ps += [l.get_particle() for l in IMP.core.get_leaves(r)]
            else:
                this_ps += atoms[n]
        if exclude_backbone:
            this_ps = [tp for tp in this_ps if IMP.atom.Atom(tp).get_atom_type() not in backbone]
        ret+=this_ps
    return list(set(ret))

def select_node(root,state_num,resolution=None,representation_type=None):
    """Eventually this should allow you to select any NODE in a hierarchy.
    E.g. I want state 5, molecule 3, copy 2.
    Selection() gets you all the leaves, but actually we want the highest? node matching.
    """
    for state in root.get_children():
        sn = IMP.atom.State(state).get_state_index()
        if sn == state_num:
            return state

def get_residue_gaps_in_hierarchy(residues,start=None,end=None,contiguous=False):
    '''
    returns the residue index gaps and contiguous segments as tuples given the IMP.pmi.residues, the first
    residue and the last residue indexes. The list is organized as
    [[1,100,"cont"],[101,120,"gap"],[121,200,"cont"]]
    /note Any residues assigned a resolution must have an IMP.atom.Residue hierarchy
              containing at least a CAlpha. For missing residues, these can be constructed
              from the PDB file

        @param residues   Input list of pmi residues
        @param start      optional. If given will set the first residue index
        @param end        optional. If given will set the last residue index
        @param contiguous If true, it will return the list of contiguous fragments

    '''

    residue_indexes=[r.get_index() for r in residues]
    if start==None:
       start=min(residue_indexes)
    if end==None:
       end=max(residue_indexes)

    stretches = []
    returned_stretches=[]
    for n, rindex in enumerate(range(start, end + 1)):
        rh=residues[n].get_hierarchy()
        rindex=residue_indexes[n]

        if len(rh.get_children()) == 0:
            if n == 0:
                # set the initial condition
                rindexgap = start
                rindexcont = start - 1
            if rindexgap == rindex - 1:
                # residue is contiguous with the previously discovered gap
                stretches[-1][1] += 1
            else:
                # residue is not contiguous with the previously discovered gap
                # hence create a new gap tuple
                stretches.append([rindex, rindex,"gap"])
            # update the index of the last residue gap
            rindexgap = rindex
        else:
            if n == 0:
                # set the initial condition
                rindexgap = start - 1
                rindexcont = start
            if rindexcont == rindex - 1:
                # residue is contiguous with the previously discovered
                # continuous part
                stretches[-1][1] += 1
            else:
                # residue is not contiguous with the previously discovered continuous part
                # hence create a new cont tuple
                stretches.append([rindex, rindex,"cont"])
            # update the index of the last residue gap
            rindexcont = rindex

    for s in stretches:
      if not contiguous and s[2]=="gap":
         returned_stretches.append((s[0],s[1]))
      if contiguous and s[2]=="cont":
         returned_stretches.append((s[0],s[1]))


    return returned_stretches


def list_chunks_iterator(input_list, length):
    """ Yield successive length-sized chunks from a list.
    """
    for i in xrange(0, len(input_list), length):
        yield input_list[i:i + length]