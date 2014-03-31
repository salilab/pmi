#!/usr/bin/env python

''' PMI2

USER ACTIONS:
add sequences (can have multiple copies and states)
add structure to those sequences
specify rigid bodies
specify what to model and what to ignore
specify resolutions (default, per-component, and per-body)
add restraints (including what resolution they act upon)

INTERNAL MODEL CONSTRUCTION (hopefully in C++):
reads coordinates, figures out gaps, creates "necklace" for missing chunks
creates all the different resolutions
stores it in a root hierarchy containing all the different resolutions in different branches of the tree.
create actual restraints

key: must store a session file
no redundant numbers
'''


import IMP
import IMP.pmi
import IMP.pmi.representation

class System():
    def __init__(self):
        self.mdl = IMP.Model()
        self.states=[]

    def create_state(self,state=None,state_num=None):
        if state_num is None:
            state_num=len(self.states)
        if state is None:
            new_state=State(self,state_num)
        else:
            new_state=state.clone()
            new_state.state_num=state_num

        self.h=IMP.atom.Hierarchy()
        IMP.atom.State.setup_particle(self.h,state_num)

        self.states.append(new_state)
        return new_state


class State():

    ''' A class to store and search various components'''

    def __init__(self,system,state_num):
        self.rep=IMP.pmi.Representation(system.mdl)
        self.components={}
        self.bodies=[]
        self.state_num=state_num

    def create_component(self,component_name,fasta_fn,copy_num=0):
        c=Component(self,component_name,fasta_fn,copy_num,state_num)
        self.components[(component_name,copy_num)]=c
        return c

    def create_body(self,rigid_segments=[],nonrigid_segments=[]):
        ''' makes this piece move together, though flexible components are allowed '''

    def clone(self):
        pass


class Component():

    ''' A part of a System '''

    def __init__(self,representation,name,fasta_fn,copy_num,state_num,id=None):
        representation.rep.create_component(name)
        self.name=name
        self.copy_num=copy_num
        self.state_num=state_num
        self.add_sequence(fasta_fn)

    def add_sequence(self,fasta_fn,id=None):
        representation.rep.add_component_sequence(self.name,fastafile,
                                                  id=id,format="FASTA")

    def add_coordinates(self,chain,res_range,offset=0):
        ''' read coordinates for this component'''

    def get_segment(self,**kwargs):
        ''' return a dictionary for creating a Selection '''

    def add_ideal_helix(self,residue_range=None):
        pass


''' notes from representation.py

add_component_pdb:  reads the pdb file then passes the particles to coarse_hierarchy()
coarse_hierarchy:   runs create simplified along backbone, which averages particle positions
                     (presumably handles gaps well)


'''
