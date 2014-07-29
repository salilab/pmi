#!/usr/bin/env python

"""Tools for reading data
"""

import IMP
import IMP.algebra
import IMP.atom
import re
from collections import defaultdict

def parse_dssp(dssp_fn, limit_to_chains=''):
    """read dssp file, get SSEs. values are all PDB residue numbering.
    Returns selection kwargs for each SSE type (stored in a dictionary).
    Example for a structure with helix A:5-7 and Beta strands A:1-3,A:9-11:
    helix : [ [ {'chain':'A','residue_indexes': [5,6,7]} ] ]
    beta  : [ [ {'chain':'A','residue_indexes': [1,2,3]},
                {'chain':'A','residue_indexes': [9,10,11]} ] ]
    loop  : same format as helix
    """
    # setup
    sses = {'helix': [],
            'beta': [],
            'loop': []}
    helix_classes = 'GHI'
    strand_classes = 'EB'
    loop_classes = [' ', '', 'T', 'S']
    sse_dict = {}
    for h in helix_classes:
        sse_dict[h] = 'helix'
    for s in strand_classes:
        sse_dict[s] = 'beta'
    for l in loop_classes:
        sse_dict[l] = 'loop'

    # read file and parse
    start = False
    # temporary beta dictionary indexed by DSSP's ID
    beta_dict = defaultdict(list)
    prev_sstype = None
    cur_sse = {'chain':'','residue_indexes':[]}
    prev_beta_id = None
    for line in open(dssp_fn, 'r'):
        fields = line.split()
        chain_break = False
        if len(fields) < 2:
            continue
        if fields[1] == "RESIDUE":
            # Start parsing from here
            start = True
            continue
        if not start:
            continue
        if line[9] == " ":
            chain_break = True
        elif limit_to_chains != '' and line[11] not in limit_to_chains:
            continue

        # gather line info
        if not chain_break:
            pdb_res_num = int(line[5:10])
            chain = line[11]
            sstype = sse_dict[line[16]]
            beta_id = line[33]

        # decide whether to extend or store the SSE
        if prev_sstype is None:
            cur_sse = {'chain':chain,'residue_indexes':[pdb_res_num]}
        elif sstype != prev_sstype or chain_break:
            # add cur_sse to the right place
            if prev_sstype in ['helix', 'loop']:
                sses[prev_sstype].append([cur_sse])
            elif prev_sstype == 'beta':
                beta_dict[prev_beta_id].append(cur_sse)
            cur_sse = {'chain':chain,'residue_indexes':[pdb_res_num]}
        else:
            cur_sse['residue_indexes'].append(pdb_res_num)
        if chain_break:
            prev_sstype = None
            prev_beta_id = None
        else:
            prev_sstype = sstype
            prev_beta_id = beta_id
    # final SSE processing
    if not prev_sstype is None:
        if prev_sstype in ['helix', 'loop']:
            sses[prev_sstype].append([cur_sse])
        elif prev_sstype == 'beta':
            beta_dict[prev_beta_id].append(cur_sse)
    # gather betas
    for beta_sheet in beta_dict:
        sses['beta'].append(beta_dict[beta_sheet])
    return sses

def parse_xlinks_davis(data_fn,
                       max_num=-1,
                       name_map={},
                       named_offsets={}):
    """ Format from Trisha Davis. Lines are:
    ignore ignore seq1 seq2 >Name(res) >Name(res) score
    @param data_fn       The data file name
    @param max_num       Maximum number of XL to read (-1 is all)
    @param name_map      Dictionary mapping text file names to the molecule name
    @param named_offsets Integer offsets to apply to the indexing in the file
    Output is dictionary of selection commands:
    data[unique_id] =
              [ { 'r1': {'molecule':'A','residue_index':5},
                  'r2': {'molecule':'B','residue_index':100},
                  'score': 123 },
                { 'r1': {'molecule':'C','residue_index':63},
                  'r2': {'molecule':'D','residue_index':94},
                  'score': 600 }
              ]
    """

    inf=open(data_fn,'r')
    data=defaultdict(list)
    for nl,l in enumerate(inf):
        if max_num==-1 or nl<max_num:
            ig1,ig2,seq1,seq2,s1,s2,score=l.split()
            score=float(score)
            m1=re.search(r'>([\w-]+)\((\w+)\)',s1)
            m2=re.search(r'>([\w-]+)\((\w+)\)',s2)

            n1=m1.group(1)
            r1=int(m1.group(2))
            if n1 in name_map:
                n1=name_map[n1]
            if n1 in named_offsets:
                r1+=named_offsets[n1]

            n2=m2.group(1)
            r2=int(m2.group(2))
            if n2 in name_map:
                n2=name_map[n2]
            if n2 in named_offsets:
                r2+=named_offsets[n2]
            xl={}
            xl['r1']={'molecule':n1,'residue_index':r1}
            xl['r2']={'molecule':n2,'residue_index':r2}
            xl['score']=score
            data[nl].append(xl)

    inf.close()
    return data
