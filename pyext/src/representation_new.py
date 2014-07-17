import IMP
import IMP.atom
import IMP.pmi
from collections import defaultdict
import IMP.pmi.structure_tools as structure_tools
from Bio import SeqIO

"""
A new representation module. It helps to construct the hierarchy
and deal with multi-state, multi-scale, multi-copies

Usage Example:

see representation_new_test.py

For each of the classes System, State, and Molecule, you store the root node and
  references to child classes (System>State>Molecule).
When you call build() on any of these classes, build() is also called for each of the child classes,
and the root IMP hierarchy is returned.
"""

#------------------------


class _SystemBase(object):
    """This is the base class for System, _State and _Molecule
    classes. It contains shared functions in common to these classes
    """

    def __init__(self,mdl=None):
        if mdl is None:
            self.mdl=IMP.Model()
        else:
            self.mdl=mdl

    def _create_hierarchy(self):
        """create a new hierarchy"""
        tmp_part=IMP.kernel.Particle(self.mdl)
        return IMP.atom.Hierarchy.setup_particle(tmp_part)

    def _create_child(self,parent_hierarchy):
        """create a new hierarchy, set it as child of the input
        one, and return it"""
        child_hierarchy=self._create_hierarchy()
        parent_hierarchy.add_child(child_hierarchy)
        return child_hierarchy

    def build(self):
        """Build the coordinates of the system.
        Loop through stored(?) hierarchies and set up coordinates!"""
        pass

#------------------------

class System(_SystemBase):
    """This class initializes the root node of the global IMP.atom.Hierarchy."""
    def __init__(self,mdl=None):
        _SystemBase.__init__(self,mdl)
        self._number_of_states = 0
        self.states = []
        self.built=False

        # the root hierarchy node
        self.hier=self._create_hierarchy()
        self.hier.set_name("System")

    def create_state(self):
        """returns a new IMP.pmi.representation_new._State(), increment the state index"""
        self._number_of_states+=1
        state = _State(self,self._number_of_states-1)
        self.states.append(state)
        return state

    def get_number_of_states(self):
        """returns the total number of states generated"""
        return self._number_of_states

    def get_hierarchy(self):
        return self.hier

    def build(self,**kwargs):
        """call build on all states"""
        if not self.built:
            for state in self.states:
                state.build(**kwargs)
            self.built=True
        return self.hier

#------------------------

class _State(_SystemBase):
    """This private class is constructed from within the System class.
    It wraps an IMP.atom.State
    """
    def __init__(self,system,state_index):
        """Define a new state
        @param system        the PMI System
        @param state_index   the index of the new state
        """
        self.mdl = system.get_hierarchy().get_model()
        self.system = system
        self.hier = self._create_child(system.get_hierarchy())
        self.hier.set_name("State_"+str(state_index))
        self.molecules = []
        IMP.atom.State.setup_particle(self.hier,state_index)
        self.built=False

    def create_molecule(self,name,sequence=None,chain_id='',molecule_to_copy=None):
        """Create a new Molecule within this State
        @param name                the name of the molecule (string) it must not
                                   contain underscores characters "_" and must not
                                   be already used
        @param sequence            sequence (string)
        @param chain_id            Chain id to assign to this molecule
        @param molecule_to_copy    Copy everything from an existing molecule. NOT IMPLEMENTED

        """
        # check the presence of underscores
        if "_" in name:
           raise WrongMoleculeName('A molecule name should not contain underscores characters')

        # check whether the molecule name is already assigned
        if name in [mol.name for mol in self.molecules]:
           raise WrongMoleculeName('Cannot use a molecule name already used')

        mol = _Molecule(self,name,sequence,chain_id,copy_num=0)
        self.molecules.append(mol)
        return mol

    def get_hierarchy(self):
        return self.hier

    def build(self,**kwargs):
        """call build on all molecules (automatically makes copies)"""
        if not self.built:
            for mol in self.molecules:
                mol.build(**kwargs)
            self.built=True
        return self.hier

#------------------------

class _Molecule(_SystemBase):
    """This private class is constructed from within the State class.
    It wraps an IMP.atom.Molecule and IMP.atom.Copy
    Structure is read using this class
    Resolutions and copies can be registered, but are only created when build() is called
    """

    def __init__(self,state,name,sequence,chain_id,copy_num):
        """Create copy 0 of this molecule.
        Arguments:
        @param state      the parent PMI State
        @param name       the name of the molecule (string)
        @param sequence   sequence (string)
        """
        # internal data storage
        self.mdl=state.get_hierarchy().get_model()
        self.state=state
        self.name=name
        self.sequence=sequence
        self.copies=[]
        self.built=False

        # create root node and set it as child to passed parent hierarchy
        self.hier = self._create_child(self.state.get_hierarchy())
        self.hier.set_name(self.name)
        IMP.atom.Copy.setup_particle(self.hier,copy_num)
        IMP.atom.Chain.setup_particle(self.hier,chain_id)

        # create Residues from the sequence
        self.residues=[]
        for ns,s in enumerate(sequence):
            r=_Residue(self,s,ns+1)
            self.residues.append(r)

    def __repr__(self):
        return self.state.get_hierarchy().get_name()+'_'+self.name+ \
            IMP.atom.Copy(self.hier).get_copy_number()

    def __getitem__(self,val):
        if isinstance(val,int):
            return self.residues[val]
        elif isinstance(val,str):
            return self.residues[int(val)-1]
        elif isinstance(val,slice):
            return set(self.residues[val])
        else:
            print "ERROR: range ends must be int or str. Stride must be int."

    def get_hierarchy(self):
        return self.hier

    def residue_range(self,a,b,stride=1):
        """get residue range. Use integers to get 0-indexing, or strings to get PDB-indexing"""
        if isinstance(a,int) and isinstance(b,int) and isinstance(stride,int):
            return set(self.residues[a:b:stride])
        elif isinstance(a,str) and isinstance(b,str) and isinstance(stride,int):
            return set(self.residues[int(a)-1:int(b)-1:stride])
        else:
            print "ERROR: range ends must be int or str. Stride must be int."

    def get_atomic_residues(self):
        """ Return a set of Residues that have associated structure coordinates """
        atomic_res=set()
        for res in self.residues:
            if len(res.hier.get_children())>0:
                atomic_res.add(res)
        return atomic_res

    def add_copy(self,pdb_fn,chain_id,res_range=[],offset=0,new_chain_id=None):
        """Create a new Molecule storing the new coordinates.
        Ensures that representations are identical to original molecule
        Will verify that the sequence is the same as that of the first structure.
        @param pdb_fn       The file to read
        @param chain_id     Chain ID to read
        @param res_range    Add only a specific set of residues
        @param offset       Apply an offset to the residue indexes of the PDB file
        @param new_chain_id If you want to set the chain ID of the copy to something
                            (defaults to what you extract from the PDB file)
        """
        if new_chain_id is None:
            new_chain_id=chain_id
        mol=_Molecule(self.state,self.name,self.sequence,new_chain_id,copy_num=len(self.copies)+1)
        self.copies.append(mol)
        new_atomic_res = mol.add_structure(pdb_fn,chain_id,res_range,offset)
        new_idxs = set([r.get_index() for r in new_atomic_res])
        orig_idxs = set([r.get_index() for r in self.get_atomic_residues()])
        if new_idxs!=orig_idxs:
            raise MoleculeCopyError("You added a structure for the copy with different residues")
        for orig,new in zip(self.residues,mol.residues):
            new.representations=orig.representations

    def add_structure(self,pdb_fn,chain_id,res_range=[],offset=0,model_num=None):
        """Read a structure and store the coordinates.
        Returns the atomic residues (as a set)
        @param pdb_fn    The file to read
        @param chain_id  Chain ID to read
        @param res_range Add only a specific set of residues
        @param offset    Apply an offset to the residue indexes of the PDB file
        @param model_num Read multi-model PDB and return that model
        \note After offset, we expect the PDB residue numbering to match the FASTA file
        """
        # get IMP.atom.Residues from the pdb file
        rhs=structure_tools.get_structure(self.mdl,pdb_fn,chain_id,res_range,offset)
        if len(rhs)>len(self.residues):
            print 'ERROR: You are loading',len(rhs), \
                'pdb residues for a sequence of length',len(self.residues),'(too many)'

        # load those into the existing pmi Residue objects, and return contiguous regions
        atomic_res=set() # collect integer indexes of atomic residues!
        for nrh,rh in enumerate(rhs):
            idx=rh.get_index()
            internal_res=self.residues[idx-1]
            if internal_res.get_code()!=IMP.atom.get_one_letter_code(rh.get_residue_type()):
                raise StructureError('ERROR: PDB residue is',
                                     IMP.atom.get_one_letter_code(rh.get_residue_type()),
                                     'and sequence residue is',internal_res.get_code())
            internal_res.set_structure(rh)
            atomic_res.add(internal_res)
        return atomic_res

    def add_representation(self,res_set=None,representation_type="balls",resolutions=[]):
        """handles the IMP.atom.Representation decorators, such as multi-scale,
        density, etc.
        @param res_set             set of PMI residues for adding the representation
        @param representation_type currently supports only balls
        @param resolutions         what resolutions to add to the residues
        """
        allowed_types=["balls"]
        if representation_type not in allowed_types:
            print "ERROR: Allowed representation types:",allowed_types
            return
        if res_set is None:
            res_set=set(self.residues)
        for res in res_set:
            res.add_representation(representation_type,resolutions)

    def build(self,merge_type="backbone",ca_centers=True,fill_in_missing_residues=True):
        """Create all parts of the IMP hierarchy
        including Atoms, Residues, and Fragments/Representations and, finally, Copies
        /note Any residues assigned a resolution must have an IMP.atom.Residue hierarchy
              containing at least a CAlpha. For missing residues, these can be constructed
              from the PDB file

        @param merge_type Principle for grouping into fragments.
                          "backbone": linear sequences along backbone are grouped
                          into fragments if they have identical sets of representations.
                          "volume": at each resolution, groups are made based on
                          spatial distance (not currently implemented)
        @param ca_centers For single-bead-per-residue only. Set the center over the CA position.
        """
        allowed_types=("backbone")
        if merge_type not in allowed_types:
            print "ERROR: Allowed merge types:",allowed_types
            return
        if not self.built:
            # for every Residue with tagged representation, build a CAlpha
            # NOT IMPLEMENTED
            structure_tools.fill_in_missing_backbone(self.residues)

            # group into Fragments along backbone
            if merge_type=="backbone":
                structure_tools.build_along_backbone(self.mdl,self.hier,self.residues,
                                                     IMP.atom.BALLS,ca_centers)


            # group into Fragments by volume
            elif merge_type=="volume":
                pass

            # build copies
            for copy in self.copies:
                copy.build()

            self.built=True

        return self.hier


#------------------------

class Sequences(object):
    """A dictionary-like wrapper for reading and storing sequence data"""
    def __init__(self,fasta_fn,name_map=None):
        """read a fasta file and extract all the requested sequences
        @param fasta_fn sequence file
        @param name_map dictionary mapping the fasta name to the stored name
        """
        self.sequences={}
        self.read_sequences(fasta_fn,name_map)
    def __len__(self):
        return len(self.sequences)
    def __contains__(self,x):
        return x in self.sequences
    def __getitem__(self,key):
        return self.sequences[key]
    def __repr__(self):
        ret=''
        for s in self.sequences:
            ret+='%s\t%s\n'%(s,self.sequences[s])
        return ret
    def read_sequences(self,fasta_fn,name_map=None):
        # read all sequences
        handle = open(fasta_fn, "rU")
        record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        handle.close()
        if name_map is None:
            for pn in record_dict:
                self.sequences[pn]=str(record_dict[pn].seq).replace("*", "")
        else:
            for pn in name_map:
                try:
                    self.sequences[name_map[pn]]=str(record_dict[pn].seq).replace("*", "")
                except:
                    print "tried to add sequence but: id %s not found in fasta file" % pn
                    exit()

#------------------------


class _Residue(object):
    """Stores basic residue information, even without structure available."""
    # Consider implementing __hash__ so you can select.
    def __init__(self,molecule,code,index):
        """setup a Residue
        @param molecule PMI Molecule to which this residue belongs
        @param code     one-letter residue type code
        @param index    PDB index
        """
        self.molecule = molecule
        self.hier = IMP.atom.Residue.setup_particle(IMP.Particle(molecule.mdl),
                                IMP.pmi.sequence_tools.get_residue_type_from_one_letter_code(code),
                                index)
        self.representations = defaultdict(set)
    def __str__(self):
        return self.get_code()
    def __repr__(self):
        return self.__str__()
    def __key(self):
        return (self.molecule,self.hier,
                frozenset((k,tuple(self.representations[k])) for k in self.representations))
    def __eq__(self,other):
        return type(other)==type(self) and self.__key() == other.__key()
    def __hash__(self):
        return hash(self.__key())
    def get_index(self):
        return self.hier.get_index()
    def get_code(self):
        return IMP.atom.get_one_letter_code(self.hier.get_residue_type())
    def get_residue_type(self):
        return self.hier.get_residue_type()
    def get_hierarchy(self):
        return self.hier
    def set_structure(self,res):
        if res.get_residue_type()!=self.hier.get_residue_type():
            raise StructureError("Adding structure to this residue, but it's the wrong type!")
        for a in res.get_children():
            self.hier.add_child(a)
            atype=IMP.atom.Atom(a).get_atom_type()
            a.get_particle().set_name('Atom %s of residue %i'%(atype.__str__().strip('"'),
                                                               self.hier.get_index()))
    def add_representation(self,rep_type,resolutions):
        self.representations[rep_type] |= set(resolutions)
