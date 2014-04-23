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


options:
python "Component" class stores coordinates, reference to the node in the hierarchy,
reference to a sequence, maybe info about breaks (otherwise all you have is a hierarchy...)


'''


import IMP
import IMP.pmi
import IMP.pmi.representation
import re

# some utility functions (to be placed in a tools file)


def read_sequence(fasta_fn, fasta_id):
    pass


def get_number_of_copies_of_a_component(system, state_num, component_name):
    ct = 0
    pattern = r'%i:%s:(\d+)' % (state_num, component_name)
    for c in system.states[state_num]:
        if re.search(pattern, c.get_name()):
            ct += 1
    return ct


def clone_component(component, new_copy_num):
    '''probably implement this in C++. copies all fragments, and all levels of resolution'''


class Handle():

    '''a class for accessing part of a Sequence.
    define: and, or, -, [] '''


class Sequence():

    ''' might simply be a string'''


class Component():

    ''' a class to store the sequence and coordinates of a component.
    it should act as a handle to itself, or at least return one easily.'''

    def __init__(self, state_num, name, copy_num,
                 sequence, hierarchy):
        self.state_num = state_num
        self.name = name
        self.copy_num = copy_num
        self.sequence = sequence
        self.hierarchy = hierarchy

    def set_coordinates(self, pdb_fn, chain, residue_range=None,
                        offset=None, read_non_water_atoms=False):
        '''read a pdb file and store the coordinates.
        does NOT add them to the hierarchy'''
        # read pdb file into an internal hierarchy

    def residue_range(self, start, stop):
        return Handle(self, start, stop)


class Body():

    ''' a class to store mixtures of subsets of components, with optional settings'''


class System():

    '''
    a class to make it easy to setup an IMP hierarchy of a system.
    The user can perform the following operations:

    1. read in sequences
    2. create Components using those sequences
    3. add known coordinates to Components
    4. create Bodies from arbitrary slices of the Components
    5. decide how to model the Bodies:
         - defining rigid and flexible parts
         - setting resolutions
         - adding beads for missing coordinates
         - add multiple copies
    6. call create(), which will set up all resolutions and bodies
       (also sets up some dictionaries so you can access the hierarchy using those classes)


    the final hierarchy always looks like:
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

    def __init__(self, model=None):
        if model is None:
            self.mdl = IMP.Model()
        else:
            self.mdl = model

        # data storage
        self.states = {}       # value: State, key: state_num
        # value: sequence key: (state_num,component_name)
        self.sequences = {}
        # value PMI.Component, key: (state_num,component_name,copy_num)
        self.components = {}
        self.bodies = []

    def add_states(self, num_states, state_to_copy=None):
        '''insert a new State in the root hierarchy'''
        n_prev = len(self.states)
        for ns in range(num_states):
            state_num = n_prev + ns
            if state_to_copy is None:
                new_state = IMP.atom.State.setup_particle(self.mdl,
                                                          self.mdl.add_particle("state%i" % state_num), state_num)
            else:
                # insert code for cloning a state
                pass
            self.root.add_child(new_state)
            self.states[state_num] = new_state

    def add_component(
        self,
        state_num,
        component_name,
        fasta_fn,
            fasta_id=None):
        '''insert a new component within a particular State.
        If you mean to add a copy of an existing component, use add_component_copies()'''

        if get_number_of_copies_of_a_component(self, state_num, component_name) != 0:
            print "error: you are trying to create a component that already exists"
            return

        # store the sequence
        seq = read_sequence(fasta_fn, fasta_id)
        self.sequences[(state_num, component_name)] = seq

        # create new node in the hierarchy and decorate it
        copy_num = 0
        name = '%i:%s:%i' % (state_num, component_name, copy_num)
        c = IMP.atom.Copy.setup_particle(
            self.mdl,
            self.mdl.add_particle(name),
            copy_num)
        self.states[state_num].add_child(c)

        # create a PMI.Component and store it
        comp = Component(self, state_num, component_name, copy_num,
                         sequence=seq,
                         hierarchy=c)
        self.components[(state_num, component_name, copy_num)] = comp
        return comp

    def add_component_copies(self, components, num_copies, transform=None):
        '''make copies if a component and add them to the hierarchy.
        NOTE: do this AFTER you complete adding coordinates, resolutions, beads'''
        n = get_number_of_copies_of_a_component(
            self,
            state_num,
            component_name)
        for i in range(num_copies):
            copy_num = i + n
            new_component = clone_component(
                self.components[(state_num, component_name, 0)],
                copy_num)
            self.states[state_num].add_child(new_component)
            self.components[(
                state_num,
                component_name,
                copy_num)] = new_component
            self.sequences[(
                state_num,
                component_name)] = read_sequence(fasta_fn,
                                                 fasta_id)

    def create_body(
        self, rigid_stuff=None, flexible_stuff=None, resolutions=[1],
            unstructured_resolution=None, representation_type='Beads'):
        '''gather arbitrary handles from various components. basically an intermediary to create()
'''

    def create(self):
        ''' loop through the states, components, and bodies and
        create the IMP hierarchy, resolutions, etc
        by default, you shouldn't have much to do -
        just add components and there should be simple defaults
        (including if you don't make rigid bodies or split into fragments)
'''

        self.root = IMP.atom.Hierarchy.setup_particle(self.mdl,
                                                      self.mdl.add_particle("root"))
