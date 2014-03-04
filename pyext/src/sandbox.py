#!/usr/bin/env python
import IMP
import IMP.pmi
import IMP.pmi.representation


class Representation():
    def __init__(self,m):
       self.repr=IMP.pmi.Representation(m)
    
    def create_component(self,"name"):
       return Component(self,"name")

class Component():
    def __init__(self,representation,name):
       representation.repr.create_component(name)
    
    def add_sequence(self,fastafile,filename,id=None):
       representation.repr.add_component_sequence(name,fastafile,id=id,format="FASTA")
    
    def create_structure(self,residuerange=None):
       pass
    
    def add_ideal_helix(1,100)
       pass