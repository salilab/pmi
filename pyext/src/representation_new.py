Simport IMP
import IMP.atom
import IMP.pmi


'''
A new representation module. It helps to construct the hierarchy
and deal with multi-state, multi-scale, multi-copies

Usage Example:

see representation_new_test.py
'''

#------------------------




class _SystemBase(object):

  '''This is the base class for System, _State and _Molecule
  classes. It contains shared functions in common to these classes'''
  
  def _create_hierarchy(self):
    '''create a new hierarchy'''
    tmp_part=IMP.kernel.Particle(self.model)    
    return IMP.atom.Hierarchy.setup_particle(tmp_part)  
  
  def _create_child(self,parent_hierarchy):
    '''create a new hierarchy, set it as child of the input
    one, and return it'''  
    child_hierarchy=self._create_hierarchy()
    parent_hierarchy.add_child(child_hierarchy)    
    return child_hierarchy
  
  def build(self):
    '''Build the coordinates of the system'''
    pass


class System(_SystemBase):
  '''This class initializes the root node of the global IMP.atom.Hierarchy.'''
  
  def __init__(self):
    self.model=IMP.Model()
    self.number_of_states=1
    # the root hierarchy node  
    self.system_hierarchy=self._create_hierarchy()
    self.system_hierarchy.set_name("System")
  
  def create_state(self): 
    '''returns a new IMP.pmi.representation_new._State(), increment the state index'''    
    self.number_of_states+=1
    return _State(self.system_hierarchy,self.number_of_states-1)
  
#------------------------

class _State(_SystemBase):
  '''This private class is called from the System class. 
  It is decorated with IMP.atom.State'''
  def __init__(self,system_hierarchy,state_index):
    '''
    Arguments:
    IMP.atom.Hierarchy                   system_hierarchy     the parent root hierarchy
    int                                  state_index          the index of the new state
    '''  
    self.model=system_hierarchy.get_model()
    self.state_hierarchy=self._create_child(system_hierarchy)
    self.state_hierarchy.set_name("State_"+str(state_index))    
    IMP.atom.State.setup_particle(self.state_hierarchy,state_index)
    
  
  def create_molecule(self,molecule_name,sequence=None):
    '''returns a new IMP.pmi.representation_new._Molecule(),
    Arguments:
    str                                  molecule_name       the name of the molecule
    IMP.pmi.representation_new.Sequence  sequence            sequence handler
    '''    
    return _Molecule(self.state_hierarchy,molecule_name,sequence)


#------------------------

class _Molecule(_SystemBase):
  '''This private class is called from the State class. 
  It is decorated with IMP.atom.Molecule and IMP.atom.Copy'''
  
  def __init__(self,state_hierarchy,molecule_name,sequence):
    '''
    Arguments:
    IMP.atom.Hierarchy                   state_hierarchy     the parent State-decorated hierarchy
    str                                  molecule_name       the name of the molecule
    IMP.pmi.representation_new.Sequence  sequence            sequence handler
    '''  
    self.model=state_hierarchy.get_model()  
    self.molecule_name=molecule_name
    self.sequence=sequence
    self.number_of_copies=1
    self.molecule_hierarchies=[]
    self.add_copy(parent_hierarchy=state_hierarchy)

  def add_copy(self,parent_hierarchy=None):
    '''Create a new copy of the Molecule. All copies are identical
    IMP.atom.Hierarchy  parent_hierarchy   if None, just use the previous copy in the list
                                           to get the parent, otherwise get give the parent hierarchy 
                                           directly.'''
    
    copy_index=self.number_of_copies-1
    if parent_hierarchy is None:
      parent_hierarchy=self.molecule_hierarchies[copy_index-1].get_parent()
    self.molecule_hierarchies.append(self._create_child(parent_hierarchy))
    self.molecule_hierarchies[copy_index].set_name(self.molecule_name+"_"+str(copy_index))    
    molecule=IMP.atom.Molecule.setup_particle(self.molecule_hierarchies[copy_index])  
    IMP.atom.Copy.setup_particle(self.molecule_hierarchies[copy_index],copy_index)    
    self.number_of_copies+=1

  def add_structure(self,**kwargs):
    '''handles the structure, such as a pdb, or
    ab-initio structure generators'''  
    pass
    
  def set_representation(self,**kwargs):  
    '''handles the IMP.atom.Representation decorators, such as multi-scale,
    density, etc.'''  
    pass


#------------------------

class Sequence(object):
  '''this class handles a sequence as a string or as a fasta file'''
  def __init__(self,**kwargs):
    pass




