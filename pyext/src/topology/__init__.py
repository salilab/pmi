"""@namespace IMP.pmi.topology
   Set up of system representation from topology files.

   * Class for storing topology elements of PMI components
   * Functions for reading these elementsfrom a formatted PMI topology file
   * Functions for converting an existing IMP hierarchy into PMI topology
   * TopologyWriter for writing PMI topology files
"""

import IMP
import IMP.atom
import IMP.pmi
import IMP.base
import csv
import os
from collections import defaultdict
import system_tools
from Bio import SeqIO

def get_residue_type_from_one_letter_code(code):
    threetoone = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                  'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                  'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                  'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                  'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'UNK': 'X'}
    one_to_three={}
    for k in threetoone:
        one_to_three[threetoone[k]] = k
    return IMP.atom.ResidueType(one_to_three[code])

class StructureError(Exception):
    pass

#------------------------

class SystemBase(object):
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

class System(SystemBase):
    """This class initializes the root node of the global IMP.atom.Hierarchy."""
    def __init__(self,mdl=None,name="System"):
        SystemBase.__init__(self,mdl)
        self._number_of_states = 0
        self.states = []
        self.built=False

        # the root hierarchy node
        self.hier=self._create_hierarchy()
        self.hier.set_name(name)

    def create_state(self):
        """returns a new IMP.pmi.representation_new._State(), increment the state index"""
        self._number_of_states+=1
        state = _State(self,self._number_of_states-1)
        self.states.append(state)
        return state

    def __repr__(self):
        return self.hier.get_name()

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

class _State(SystemBase):
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

    def __repr__(self):
        return self.system.__repr__()+'.'+self.hier.get_name()

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
        if name in [mol.get_name() for mol in self.molecules]:
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

class _Molecule(SystemBase):
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
        self.sequence=sequence
        self.copies=[]
        self.built=False

        # create root node and set it as child to passed parent hierarchy
        self.hier = self._create_child(self.state.get_hierarchy())
        self.hier.set_name(name)
        IMP.atom.Copy.setup_particle(self.hier,copy_num)
        IMP.atom.Chain.setup_particle(self.hier,chain_id)

        # create Residues from the sequence
        self.residues=[]
        for ns,s in enumerate(sequence):
            r=_Residue(self,s,ns+1)
            self.residues.append(r)

    def __repr__(self):
        return self.state.__repr__()+'.'+self.hier.get_name()+'.'+ \
            str(IMP.atom.Copy(self.hier).get_copy_index())


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

    def get_name(self):
        return self.hier.get_name()

    def residue_range(self,a,b,stride=1):
        """get residue range. Use integers to get 0-indexing, or strings to get PDB-indexing"""
        if isinstance(a,int) and isinstance(b,int) and isinstance(stride,int):
            return set(self.residues[a:b:stride])
        elif isinstance(a,str) and isinstance(b,str) and isinstance(stride,int):
            return set(self.residues[int(a)-1:int(b)-1:stride])
        else:
            print "ERROR: range ends must be int or str. Stride must be int."

    def get_residues(self):
        """ Return all Residues as a set"""
        all_res=set()
        for res in self.residues:
            all_res.add(res)
        return all_res

    def get_atomic_residues(self):
        """ Return a set of Residues that have associated structure coordinates """
        atomic_res=set()
        for res in self.residues:
            if len(res.hier.get_children())>0:
                atomic_res.add(res)
        return atomic_res

    def get_non_atomic_residues(self):
        """ Return a set of Residues that don't have associated structure coordinates """
        non_atomic_res=set()
        for res in self.residues:
            if len(res.hier.get_children())==0:
                non_atomic_res.add(res)
        return non_atomic_res

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
        mol=_Molecule(self.state,self.get_hierarchy().get_name(),
                      self.sequence,new_chain_id,copy_num=len(self.copies)+1)
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
        rhs=system_tools.get_structure(self.mdl,pdb_fn,chain_id,res_range,offset)
        if len(rhs)>len(self.residues):
            print 'ERROR: You are loading',len(rhs), \
                'pdb residues for a sequence of length',len(self.residues),'(too many)'

        # load those into the existing pmi Residue objects, and return contiguous regions
        atomic_res=set() # collect integer indexes of atomic residues!
        for nrh,rh in enumerate(rhs):
            idx=rh.get_index()
            internal_res=self.residues[idx-1]
            if internal_res.get_code()!=IMP.atom.get_one_letter_code(rh.get_residue_type()):
                raise StructureError('ERROR: PDB residue index',idx,'is',
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

            # group into Fragments along backbone
            if merge_type=="backbone":
                system_tools.build_along_backbone(self.mdl,self.hier,self.residues,
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
                                get_residue_type_from_one_letter_code(code),
                                index)
        self.representations = defaultdict(set)
    def __str__(self):
        return self.get_code()+str(self.get_index())
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
    def get_has_coordinates(self):
        if len(self.hier.get_children())>0:
            return True
        else:
            return False
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


class TopologyReader(object):
    '''
    This class reads in a standard pipe-delimited PMI topology file
    and stores the items as a ComponentTopology class for input into IMP.pmi.autobuild_model()
    '''
    def __init__(self, topology_file):
        self.topology_file=topology_file
        self.component_list=[]
        self.defaults={'bead_size'                : 10,
                       'residue_range'            : 'all',
                       'pdb_offset'               : 0,
                       'em_residues_per_gaussian' : 0};
        self.component_list=self.import_topology_file(topology_file)


    def write_topology_file(self,outfile):
        f=open(outfile, "w")
        f.write("|directories|\n")
        #print self.defaults
        for key, value in self.defaults.iteritems():
            output="|"+str(key)+"|"+str(value)+"|\n"
            f.write(output)
        f.write("\n\n")
        f.write("|topology_dictionary|\n")
        f.write("|component_name|domain_name|fasta_fn|fasta_id|pdb_fn|chain|residue_range|pdb_offset|bead_size|em_residues_per_gaussian|\n")
        for c in self.component_list:
            output="|"+str(c.name)+"|"+str(c.domain_name)+"|"+str(c.fasta_file)+"|"+str(c.fasta_id)+"|"+str(c.pdb_file)+"|"+str(c.chain)+"|"+str(c.residue_range).strip("(").strip(")")+"|"+str(c.pdb_offset)+"|"+str(c.bead_size)+"|"+str(c.em_residues_per_gaussian)+"|\n"
            f.write(output)
        return outfile

    def get_component_topologies(self, topology_list = "all"):
        """ Return list of ComponentTopologies for selected components given a list of indices"""
        if topology_list == "all":
            topologies = self.component_list
        else:
            topologies=[]
            for i in topology_list:
                topologies.append(self.component_list[i])
        return topologies

    def set_dir(self, default_dir, new_dir):
        """ Changes the default directories and renames the files for each ComponentTopology object """
        if default_dir in self.defaults.keys():
            self.defaults[default_dir]=new_dir
        else:
            print default_dir, "is not a correct directory key"
            exit()
        for c in self.component_list:
            pdb_file=c.pdb_file.split("/")[-1]
            c.pdb_file=self._make_path(self.defaults['pdb_dir'],
                                        pdb_file)
            fasta_file=c.fasta_file.split("/")[-1]
            c.fasta_file=self._make_path(self.defaults['fasta_dir'],
                                        fasta_file)
            if c.gmm_file is not None:
                gmm_file=c.gmm_file.split("/")[-1]
                c.gmm_file=self._make_path(self.defaults['gmm_dir'],
                                        gmm_file)
                mrc_file=c.mrc_file.split("/")[-1]
                c.mrc_file=self._make_path(self.defaults['gmm_dir'],
                                        mrc_file)


    def import_topology_file(self, topology_file, append=False):
        """ Import system components from topology file. append=False will erase current topology and overwrite with new """
        is_defaults=False
        is_topology=False
        defaults_dict={}
        linenum=1

        if append==False:
            self.component_list=[]

        with open(topology_file) as infile:
            for line in infile:

                if line.lstrip()=="" or line[0]=="#":
                    continue

                elif line.split('|')[1]=="topology_dictionary":
                    is_topology=True

                elif is_topology==True and is_defaults==True:
                # Store the field names for this topology grid
                    topology_fields=line
                    is_defaults=False

                elif is_topology==True:
                # create a component_topology from this line
                    new_component=self.create_component_topology(line, topology_fields, self.defaults, linenum)
                    self.component_list.append(new_component)

                elif is_defaults==True:
                # grab value for default and put into class attribute
                    self.add_default_parameter(line, linenum)

                elif line.split('|')[1]=="directories":
                    is_defaults=True

                #print line, is_defaults, is_topology
                linenum=linenum+1
                #print self.defaults
        return self.component_list

    def _make_path(self, dirname, fname):
        "Get the full path to a file, possibly relative to the topology file"
        dirname = IMP.base.get_relative_path(self.topology_file, dirname)
        return os.path.join(dirname, fname)

    def create_component_topology(self, component_line, topology_fields, defaults, linenum, color="0.1"):

    #Reads a grid of topology values and matches them to their key.
    #Checks each value for correct syntax
    #Returns a list of ComponentTopology objects

        fields=topology_fields.split('|')
        values=component_line.split('|')
        c=ComponentTopology()
        no_error=True
    ##### Required fields
        c.name          = values[fields.index("component_name")].strip()
        c.domain_name   = values[fields.index("domain_name")].strip()
        c.fasta_file    = self._make_path(defaults['fasta_dir'],
                                          values[fields.index("fasta_fn")])
        c.fasta_id      = values[fields.index("fasta_id")].strip()
        c.pdb_file      = self._make_path(defaults['pdb_dir'],
                                          values[fields.index("pdb_fn")])
        # Need to find a way to define color
        c.color         = 0.1

        # PDB Chain
        if len(values[fields.index("chain")])==1 and values[fields.index("chain")].isupper()==True:
            c.chain = values[fields.index("chain")]
        else:
            print "PDB Chain format for component ", c.name, ", line ", linenum, " is not correct"
            print "Correct syntax is a single uppercase letter. |", values[fields.index("chain")], "| was given."
            no_error=False

    ##### Optional fields
        # Residue Range
        if "residue_range" in fields:
            f=values[fields.index("residue_range")]
            if f.strip()=='all' or str(f)=="":
                c.residue_range=(1,-1)
            # Make sure that is residue range is given, there are only two values and they are integers
            elif len(f.split(','))==2 and self.is_int(f.split(',')[0]) and self.is_int(f.split(',')[1]):
                c.residue_range=(int(f.split(',')[0]), int(f.split(',')[1]))
            else:
                print "Residue Range format for component ", c.name, ", line ", linenum, " is not correct"
                print "Correct syntax is two comma separated integers:  |start_res, end_res|. |", f, "| was given."
                print "To select all residues, indicate |\"all\"|"
                no_error=False
        else:
            c.residue_range=defaults["residue_range"]


        # PDB Offset
        if "pdb_offset" in fields:
            f=values[fields.index("pdb_offset")]
            if self.is_int(f):
                c.pdb_offset=int(f)
            else:
                print "PDB Offset format for component ", c.name, ", line ", linenum, " is not correct"
                print "The value must be a single integer. |", f, "| was given."
                no_error=False
        else:
            c.pdb_offset=defaults["pdb_offset"]

        # Bead Size
        if "bead_size" in fields:
            f=values[fields.index("bead_size")]
            if self.is_int(f):
                c.bead_size=int(f)
            else:
                print "Bead Size format for component ", c.name, ", line ", linenum, " is not correct"
                print "The value must be a single integer. |", f, "| was given."
                no_error=False
        else:
            c.bead_size=defaults["bead_size"]

        # EM Residues Per Gaussian
        if "em_residues_per_gaussian" in fields:
            f=values[fields.index("em_residues_per_gaussian")]
            if self.is_int(f):
                if int(f) > 0:
                    c.gmm_file=self._make_path(defaults['gmm_dir'],
                                           c.domain_name.strip() + ".txt")
                    c.mrc_file=self._make_path(defaults['gmm_dir'],
                                           c.domain_name.strip() + ".mrc")
                c.em_residues_per_gaussian=int(f)
            else:
                print "em_residues_per_gaussian format for component ", c.name, ", line ", linenum, " is not correct"
                print "The value must be a single integer. |", f, "| was given."
                no_error=False
        else:
            c.em_residues_per_gaussian=defaults["em_residues_per_gaussian"]

        if no_error==True:
            return c
        else:
            print "Fix Topology File syntax errors and rerun.  Exiting..."
            exit()


    def is_int(self, s):
       # is this string an integer?
        try:
            float(s)
            return float(s).is_integer()
        except ValueError:
            return False


    def add_default_parameter(self,line, linenum):
    #Separates a line into a key:value pair.

        f=line.split('|')
        if len(f) != 4:
            print "Default value syntax not correct for ", line
            print "Line number", linenum," contains ", len(f)-2, " fields."
            print "Please reformat to |KEY|VALUE|"
        self.defaults[f[1]]=f[2]



class ComponentTopology(object):
    '''
    Topology class stores the components required to build a standard IMP hierarchy
    using IMP.pmi.autobuild_model()
    '''
    def __init__(self):
        self.name=None
        self.domain_name=None
        self.fasta_file=None
        self.fasta_id=None
        self.pdb_file=None
        self.chain=None
        self.residue_range=None
        self.pdb_offset=None
        self.bead_size=None
        self.em_residues_per_gaussian=None
        self.gmm_file=None
        self.mrc_file=None
        self.color=None

    def recompute_default_dirs(self, topology):
        pdb_filename=self.pdb_file.split("/")[-1]
        self.pdb_filename=IMP.base.get_relative_path(topology.topology_file, topology.defaults)
