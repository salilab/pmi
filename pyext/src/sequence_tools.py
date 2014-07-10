#!/usr/bin/env python
"""utilities for interpreting and analyzing sequence data"""

import IMP
import IMP.atom

def get_residue_type_from_one_letter_code(code):
    threetoone = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                  'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                  'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                  'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                  'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'UNK': 'X'}
    one_to_three={}
    for k in threetoone:
        one_to_three[threetoone[k]] = k
    return IMP.atom.ResidueType(one_to_three[code])
