#!/usr/bin/env python
"""useful tools for finding things in an IMP hierarchy"""

import IMP
import IMP.atom

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


def select_node(root,state_num,resolution=None,representation_type=None):
    """Eventually this should allow you to select any NODE in a hierarchy.
    E.g. I want state 5, molecule 3, copy 2.
    Selection() gets you all the leaves, but actually we want the highest? node matching.
    """
    for state in root.get_children():
        sn = IMP.atom.State(state).get_state_index()
        if sn == state_num:
            return state
