#!/usr/bin/env python
from collections import defaultdict

class SelectionDictionary(object):
    """A class for storing kwargs to atom::Selection.
    Ensures it only stores accetable Selection keywords
    """
    def __init__(self,model,data={},**kwargs):
        """Setup a SelectionDictionary with an initial dictionary or just use kwargs
        @param model The IMP model you are using (needed for testing purposes)
        @param data  A dictionary containing all the selection keywords you want
        """
        self.add_fields(**data)
        self.add_fields(kwargs)
        self.model = model
        self.mh = IMP.atom.Hierarchy(IMP.Particle(self.model))
        self.data={}

    def get_data(self):
        """Get the Selection arguments"""
        return data

    def add_fields(self,**kwargs):
        """Add some new keys and values to the dictionary.
        Checks to make sure they are valid IMP selection criteria.
        If those keys already exist, the new values will overrwrite old values!
        """
        for key,value in kwargs:
            try:
                IMP.atom.Selection(self.mh,key=value)
            except:
                raise Exception("ERROR: you provided invalid key "+key)
        self.data.update(kwargs)

    def update(self,sdict):
        """Append another SelectionDictionary to this one."""
        self.add_fields(**sdict.get_data())

    def __getitem__(self,key):
        return self.data[key]

    def __setitem__(self,key,value):
        self.add_fields(key=value)

    def __iter__(self):
        return iter(self.data)

class CrossLinkData(object):
    """ This class handles and stores cross-link data sets.
    To be used with cross-link restraints. """

    def __init__(self,global_sdict=None):
        """Class to store info about XL's, conveniently as kwargs
        \param global_sdict Additional kwargs to add to all selections
                            (e.g., SelectionDictionary(atom_type=IMP.atom.AtomType('CA')))
        """
        self.data=defaultdict(list)
        self.global_sdict=global_sdict

    def add_cross_link(self,unique_id,sdict1,sdict2,**kwargs):
        """ Append a new cross-link to a given unique ID
        @param unique_id the unique ID of the identified precursor ion
        @param sdict1    SelectionDictionary for the first atom
                         e.g. SelectionDictionary(molecule='A',residue_index=5)
        @param sdict2    SelectionDictionary for the second atom
        \note can provide additional kwargs if you have more descriptors for this XL (e.g., score)
        """
        if self.global_sdict is not None:
            sdict1.update(self.global_sdict)
            sdict2.update(self.global_sdict)
        this_xl = {'r1':sdict1,'r2':sdict2}
        this_xl.update(kwargs)
        self.data[unique_id].append(this_xl)

    def get_data(self):
        return self.data


class SubsequenceData(object):
    """ A convenient way to store lists of subsets of your molecule.
    These can be labeled for easier retrieval
    Use cases: storing lists of secondary structures from DSSP or PSIPRED
               storing lists of molecules that should be made symmetric
    """
    def __init__(self,global_sdict=None):
        """Class to store info about XL's, conveniently as kwargs
        @param global_sdict Additional SelectionDictionary to add to all selections
                            (e.g., SelectionDictionary(atom_type=IMP.atom.AtomType('CA')) )
        """
        self.data=defaultdict(list)
        self.global_sdict = global_sdict


    def add_subsequence(self,label,sdicts):
        """ Append a new cross-link to a given unique ID
        @param label   label for this subsequence (e.g., 'helix')
        @param sdicts  list of Selection Dictionaries for this subset
                       e.g. [SelectionDictionary(molecule='A',residue_indexes=range(5,10))
                             SelectionDictionary(molecule='B',residue_indexes=range(90,100))]
        """
        tmp=[]
        if self.global_sdict is not None:
            for sd in sdicts:
                sd.update(self.global_sdict)
        self.data[label].append(tmp)

    def get_data(self):
        return self.data
