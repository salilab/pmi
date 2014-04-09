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
import re

# some utility functions (to be placed in a tools file)
def read_sequence(fasta_fn,fasta_id):
    pass

def get_number_of_copies_of_a_component(system,state_num,component_name):
    ct=0
    pattern=r'%i:%s:(\d+)'%(state_num,component_name)
    for c in system.states[state_num]:
        if re.search(pattern,c.get_name()):
            ct+=1
    return ct

def clone_component(component,new_copy_num):
    '''probably implement this in C++. copies all fragments, and all levels of resolution'''

class System():
    '''
    a class to maintain an IMP hierarchy of a system.
    contains dictionaries to easily access states and components.
    all the functions concern adding or editing components of the root hierachy.

    data structure:
    root ->
      atom::State ->
        component copy 0 (==atom::Copy) ->
          Fragment (==atom::Resolution) ->
            res0,
            res1,
            res10,
            ...
        component copy 1
           (identical to copy 0)

    '''

    def __init__(self):
        self.mdl=IMP.Model()
        self.root=IMP.atom.Hierarchy.setup_particle(self.mdl,
                                                    self.mdl.add_particle("root"))
        # internal data
        self.sequences={}    # key: (state_num, component_name)

        # the values of these dictionaries are members of the hierarchy
        self.states={}       # key: state_num
        self.components={}   # key: (state_num,component_name,copy_num)

    def add_state(self,state_num=None,state_to_copy=None):
        '''insert a new State in the root hierarchy'''
        if state_num is None:
            state_num=len(self.states)
        if state_to_copy is None:
            new_state=IMP.atom.State.setup_particle(self.mdl,
                      self.mdl.add_particle("state%i"%state_num),state_num)
        else:
            pass
        self.root.add_child(new_state)
        self.states[state_num]=new_state
        return new_state

    def add_component(self,state_num,component_name,fasta_fn=None,fasta_id=None):
        '''insert a new component within a particular State.
        If you mean to add a copy of an existing component, use add_component_copies()'''

        if get_number_of_copies_of_a_component(self,state_num,component_name)!=0:
            print "error: you are trying to create a component that already exists"
            return
        copy_num=0
        name='%i:%s:%i'%(state_num,component_name,copy_num)
        new_component=IMP.atom.Copy.setup_particle(self.mdl,self.mdl.add_particle(name),
                                                   copy_num)

        self.states[state_num].add_child(new_component)
        self.components[(state_num,component_name,copy_num)]=new_component
        self.sequences[(state_num,component_name)]=read_sequence(fasta_fn,fasta_id)

    def add_coordinates(self,state_num,component_name,pdb_fn,residue_range=None,offset=None):
        '''read a pdb file and assign its coordinates to a particular component.
        NOTE: only edits copy 0. run this before you create additional copies'''

    def add_resolution_to_component(self,state_num,component_name,resolution,
                                    residue_range=None):
        '''simplify coordinates along backbone and add to the hierarchy.
        NOTE: only edits copy 0. run this before you create additional copies'''

    def add_beads_for_missing_residues(self,state_num,component_name,resolution,
                                       residue_range=None):
        '''search component for missing residues (those with no coordinates) and
        add a string of beads.
        NOTE: only edits copy 0. run this before you create additional copies'''

    def add_component_copies(self,component_name,state_num,num_copies,transform=None):
        '''make copies if a component and add them to the hierarchy.
        NOTE: do this AFTER you complete adding coordinates, resolutions, beads'''
        n=get_number_of_copies_of_a_component(self,state_num,component_name)
        for i in range(num_copies):
            copy_num=i+n
            new_component=clone_component(self.components[(state_num,component_name,0)],
                                          copy_num)
            self.states[state_num].add_child(new_component)
            self.components[(state_num,component_name,copy_num)]=new_component
            self.sequences[(state_num,component_name)]=read_sequence(fasta_fn,fasta_id)
