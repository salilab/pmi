# tools for reading data

import IMP
import IMP.algebra
import IMP.atom

def parse_dssp(dssp_fn, limit_to_chains=''):
    ''' read dssp file, get SSEs. values are all PDB residue numbering. returns dict of sel tuples
helix : [ [ ['A',5,7] ] , [['B',15,17]] , ...] two helices A:5-7,B:15-17
beta  : [ [ ['A',1,3] , ['A',100,102] ] , ...] one sheet: A:1-3 & A:100-102
loop  : same format as helix, it's the contiguous loops
'''

    from collections import defaultdict

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
            break

        # gather line info
        if not chain_break:
            pdb_res_num = int(line[5:10])
            chain = line[11]
            sstype = line[16]
            beta_id = line[33]

        # decide whether to extend or store the SSE
        if prev_sstype is None:
            cur_sse = {'chain':chain,'residue_indexes':[pdb_res_num]}
        elif sstype != prev_sstype or chain_break:
            # add cur_sse to the right place
            if sse_dict[prev_sstype] in ['helix', 'loop']:
                sses[sse_dict[prev_sstype]].append([cur_sse])
            if sse_dict[prev_sstype] == 'beta':
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
        if sse_dict[prev_sstype] in ['helix', 'loop']:
            sses[sse_dict[prev_sstype]].append([cur_sse])
        if sse_dict[prev_sstype] == 'beta':
            beta_dict[prev_beta_id].append(cur_sse)

    # gather betas
    for beta_sheet in beta_dict:
        sses['beta'].append(beta_dict[beta_sheet])

    return sses
