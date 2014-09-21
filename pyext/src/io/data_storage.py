#!/usr/bin/env python
import IMP
import IMP.atom
from collections import defaultdict

class SelectionDict(object):
    """A class for storing kwargs to atom::Selection.
    This is useful when you want to read some experimental data but
    don't yet want to make a selection
    \note Has a check to ensures it only stores acceptable Selection keywords
    """
    def __init__(self,model,data={},**kwargs):
        """Setup a SelectionDict with an initial dictionary or just use kwargs
        @param model The IMP model you are using (needed for testing purposes)
        @param data  A dictionary containing all the selection keywords you want
        """
        self.model = model
        self.mh = IMP.atom.Hierarchy(IMP.Particle(self.model))
        self.data={}
        self.add_fields(**data)
        self.add_fields(**kwargs)

    def get_data(self):
        """Get the Selection arguments"""
        return self.data

    def add_fields(self,**kwargs):
        """Add some new keys and values to the dictionary.
        Checks to make sure they are valid IMP selection criteria.
        If those keys already exist, the new values will overrwrite old values!
        """
        for key,value in kwargs.iteritems():
            IMP.atom.Selection(self.mh,**{key:value})
            try:
                IMP.atom.Selection(self.mh,**{key:value})
            except:
                raise Exception("ERROR: you provided invalid key "+key)
        self.data.update(kwargs)

    def update(self,sdict):
        """Append another SelectionDict to this one."""
        self.add_fields(**sdict.get_data())

    def __getitem__(self,key):
        return self.data[key]

    def __setitem__(self,key,value):
        self.add_fields(key=value)

    def __iter__(self):
        return iter(self.data)

    def __repr__(self):
        return self.data.__repr__()

    def __cmp__(self,other):
        return self.data.__cmp__(other.get_data())

    def keys(self):
        return self.data.keys()

    def select(self,mh):
        return IMP.atom.Selection(mh,**self.data)

class CrossLinkData(object):
    """ This class handles and stores cross-link data sets.
    To be used with cross-link restraints. """

    def __init__(self,mdl,global_sdict=None):
        """Class to store info about XL's, conveniently as kwargs
        \param global_sdict Additional kwargs to add to all selections
                            (e.g., SelectionDict(atom_type=IMP.atom.AtomType('CA')))
        """
        self.mdl = mdl
        self.data=defaultdict(list)
        self.global_sdict=global_sdict

    def add_cross_link(self,unique_id,sdict1,sdict2,**kwargs):
        """ Append a new cross-link to a given unique ID
        @param unique_id the unique ID of the identified precursor ion
        @param sdict1    SelectionDict for the first atom
                         e.g. SelectionDict(molecule='A',residue_index=5)
        @param sdict2    SelectionDict for the second atom
        \note can provide additional kwargs if you have more descriptors for this XL (e.g., score)
        """
        if type(sdict1) is dict:
            sdict1 = SelectionDict(self.mdl,sdict1)
        if type(sdict2) is dict:
            sdict2 = SelectionDict(self.mdl,sdict2)
        if self.global_sdict is not None:
            sdict1.update(self.global_sdict)
            sdict2.update(self.global_sdict)
        this_xl = {'r1':sdict1,'r2':sdict2}
        this_xl.update(kwargs)
        self.data[unique_id].append(this_xl)

    def __getitem__(self,key):
        return self.data[key]

    def __len__(self):
        return len(self.data)

class SubsequenceData(object):
    """ A convenient way to store labeled lists of subsets of your molecule.
    Each label points to a list, in which each item is a list of SelectionDictionaries.
    Use cases: storing lists of secondary structures from DSSP or PSIPRED
               storing lists of molecules that should be made symmetric
    """
    def __init__(self,mdl,global_sdict=None):
        """Class to store info about XL's, conveniently as kwargs
        @param global_sdict Additional SelectionDict to add to all selections
                            (e.g., SelectionDict(atom_type=IMP.atom.AtomType('CA')) )
        """
        self.mdl = mdl
        self.data=defaultdict(list)
        self.global_sdict = global_sdict

    def add_subsequence(self,label,sdicts):
        """ Append a new cross-link to a given unique ID
        @param label   label for this subsequence (e.g., 'helix')
        @param sdicts  list of SelectionDicts (or just dicts) for this subset
                       e.g. you could add a beta sheet with:
                       [SelectionDict(molecule='A',residue_indexes=range(5,10))
                        SelectionDict(molecule='B',residue_indexes=range(90,95))]
        """
        if self.global_sdict is not None:
            for sd in sdicts:
                sd.update(self.global_sdict)
        tmp=[]
        for sdict in sdicts:
            if type(sdict) is SelectionDict:
                tmp.append(sdict)
            elif type(sdict) is dict:
                tmp.append(SelectionDict(self.mdl,sdict))
            else:
                print 'subsequence must consist of list of dict or SelectionDict'
        self.data[label].append(tmp)

    def get_data(self):
        return self.data

    def __getitem__(self,key):
        return self.data[key]

    def __repr__(self):
        return self.data.__repr__()

    def keys(self):
        return self.data.keys()
