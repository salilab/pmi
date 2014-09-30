#!/usr/bin/env python
import IMP
import IMP.atom
from collections import defaultdict

class Subsequence(object):
    """ A light class to store multiple not-necessarily-contiguous residue ranges."""
    def __init__(self,chain=None,molecule=None,residue_tuple=None,subsequences=None):
        """Create subsequence and optionally pass the first contiguous range.
        @param chain The chain ID
        @param molecule The molecule name
        @param residue_tuple PDB-style inclusive residue range
        @param subsequences A list of other subsequences to combine (not implemented)
        """
        self.seqs=[]
        if chain or molecule or residue_tuple:
            self.add_range(chain,molecule,residue_tuple)
        if subsequences:
            pass
    def add_range(self,chain=None,molecule=None,residue_tuple=None):
        """Add some stuff to this subsequence
        @param chain The chain ID
        @param molecule The molecule name
        @param residue_tuple PDB-style inclusive residue range
        """
        self.seqs.append({'chain':chain,'molecule':molecule,'residue_tuple':residue_tuple})
    def get_selection(self,hier,**kwargs):
        """Create an IMP Selection from this subsequence
        @param hier An IMP hierarchy or list of them
        \note any additional keyword arguments will be appended to the selection
        """
        for nseq,seq in enumerate(self.seqs):
            args=kwargs
            if seq['chain']:
                args['chain']=seq['chain']
            if seq['molecule']:
                args['molecule']=seq['molecule']
            if seq['residue_tuple']:
                args['residue_indexes']=range(seq['residue_tuple'][0],
                                              seq['residue_tuple'][1]+1)
            sel = IMP.atom.Selection(hier,**args)
            if nseq==0:
                ret=sel
            else:
                ret|=sel
        return ret
    def __repr__(self):
        rep=''
        for nseq,seq in enumerate(self.seqs):
            this_str=[]
            if seq['residue_tuple'] is not None:
                this_str.append('%i-%i'%(seq['residue_tuple'][0],seq['residue_tuple'][1]))
            if seq['molecule'] is not None:
                this_str.append('%s'%seq['molecule'])
            if seq['chain'] is not None:
                this_str.append('%s'%seq['chain'])
            rep+='.'.join(this_str)
            if nseq < len(self.seqs)-1:
                rep+='_'
        return rep

class SubsequenceData(object):
    """ Group a bunch of subsequences with certain labels
    Use cases: storing lists of secondary structures from DSSP or PSIPRED
               storing lists of molecules that should be made symmetric
    """
    def __init__(self):
        """Setup groups of subsequences
        """
        self.data=defaultdict(list)

    def add_subsequence(self,label,subsequence):
        """ Append a new cross-link to a given unique ID
        @param label         label for this subsequence (e.g., 'helix')
        @param subsequence   a Subsequence object to store under that label
                             e.g. Subsequence(chain='A',residue_tuple=(10,15))
        """
        if type(subsequence) is not Subsequence:
            raise InputError('must provide a subsequence object')
        self.data[label].append(subsequence)

    def __getitem__(self,key):
        return self.data[key]

    def __repr__(self):
        return self.data.__repr__()

    def keys(self):
        return self.data.keys()
