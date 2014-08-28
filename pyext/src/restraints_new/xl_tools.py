#!/usr/bin/env python


class CrossLinkDataSet(object):
    """ This class handles and stores cross-link data sets. 
    To be used with cross-link restraints. """
    def __init__(self):
        self.data={}
    
    def read_data_set(self,filename,format='csv'):
       """ Read a data set from a file. Not yet implemented. """
       pass
       
    def add_unique_id(self,unique_id):
       """ Create a new unique ID 
        @param unique_id the unique_id of the identified precursor ion       
       """    
       self.data[unique_id]=[]
    
    def add_cross_link(self,unique_id,sel1,sel2,**kwargs):
       """ Append a new cross-link to a given unique ID    
        @param unique_id the unique ID of the identified precursor ion
        @param sel1      IMP.Selection of the first atom
        @param sel2      IMP.Selection of the second atom
        @param **kwargs  additional descriptors
       """ 
       self.data[unique_id].append[{"Selection1":sel1,"Selection2":sel2}]
       self.data[unique_id][-1].update(kwargs)
    
    def get_data(self):
        return self.data