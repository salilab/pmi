"""@namespace IMP.pmi.topology
Set of python classes to create a multi-state, multi-resolution IMP hierarchy.
* Start by creating a System with `model = IMP.Model(); s = IMP.pmi.topology.System(model)`. The System will store all the states.
* Then call System.create_state(). You can easily create a multistate system by calling this function multiples times.
* For each State, call State.create_molecule() to add a Molecule (a uniquely named polymer). This function returns the Molecule object which can be passed to various PMI functions.
* Some useful functions to help you set up your Molecules:
 * Access the sequence residues with slicing (Molecule[a:b]) or functions like Molecule.get_atomic_residues() and Molecule.get_non_atomic_residues(). These functions all return python sets for easy set arithmetic using & (and), | (or), - (difference)
 * Molecule.add_structure() to add structural information from a PDB file.
 * Molecule.add_representation() to create a representation unit - here you can choose bead resolutions as well as alternate representations like densities or ideal helices.
 * Molecule.create_clone() lets you set up a molecule with identical representations, just a different chain ID. Use Molecule.create_copy() if you want a molecule with the same sequence but that allows custom representations.
* Once data has been added and representations chosen, call System.build() to create a canonical IMP hierarchy.
* Following hierarchy construction, setup rigid bodies, flexible beads, etc in IMP::pmi::dof.
* Check your representation with a nice printout: IMP::atom::show_with_representation()
See a [comprehensive example](https://integrativemodeling.org/nightly/doc/ref/pmi_2multiscale_8py-example.html) for using these classes.

Alternatively one can construct the entire topology and degrees of freedom via formatted text file with TopologyReader and IMP::pmi::macros::BuildSystem(). This is used in the [PMI tutorial](@ref rnapolii_stalk).
Note that this only allows a limited set of the full options available to PMI users (rigid bodies only, fixed resolutions).
"""

from __future__ import print_function
import IMP
import IMP.atom
import IMP.algebra
import IMP.pmi
import IMP.pmi.tools
import csv
import os
from collections import defaultdict
from bisect import bisect_left
from math import pi,cos,sin
from operator import itemgetter

def _build_ideal_helix(model, residues, coord_finder):
    """Creates an ideal helix from the specified residue range
    Residues MUST be contiguous.
    This function actually adds them to the TempResidue hierarchy
    """
    created_hiers = []

    # this function creates a CAlpha helix structure (which can be used for coarsening)
    for n, tempres in enumerate(residues):
        if tempres.get_has_structure():
            raise Exception("You tried to build ideal_helix for a residue "
                            "that already has structure:",tempres)
        if n>0 and (not tempres.get_index()==prev_idx+1):
            raise Exception("Passed non-contiguous segment to build_ideal_helix for",tempres.get_molecule())

        # New residue particle will replace the TempResidue's existing (empty) hierarchy
        rp = IMP.Particle(model)
        rp.set_name("Residue_%i" % tempres.get_index())

        # Copy the original residue type and index
        this_res = IMP.atom.Residue.setup_particle(rp,tempres.get_hierarchy())

        # Create the CAlpha
        ap = IMP.Particle(model)
        d = IMP.core.XYZR.setup_particle(ap)
        x = 2.3 * cos(n * 2 * pi / 3.6)
        y = 2.3 * sin(n * 2 * pi / 3.6)
        z = 6.2 / 3.6 / 2 * n * 2 * pi / 3.6
        d.set_coordinates(IMP.algebra.Vector3D(x, y, z))
        d.set_radius(1.0)
        a = IMP.atom.Atom.setup_particle(ap, IMP.atom.AT_CA)  # Decorating as Atom also decorates as Mass
        IMP.atom.Mass(ap).set_mass(110.0)
        this_res.add_child(a)

        # Add this structure to the TempResidue
        tempres.set_structure(this_res)
        created_hiers.append(this_res)
        prev_idx = tempres.get_index()
    coord_finder.add_residues(created_hiers) #the coord finder is for placing beads (later)


class Sequences(object):
    """A dictionary-like wrapper for reading and storing sequence data"""
    def __init__(self,fasta_fn,name_map=None):
        """read a fasta file and extract all the requested sequences
        @param fasta_fn sequence file
        @param name_map dictionary mapping the fasta name to final stored name
        """
        self.sequences = IMP.pmi.tools.OrderedDict()
        self.read_sequences(fasta_fn,name_map)
    def __len__(self):
        return len(self.sequences)
    def __contains__(self,x):
        return x in self.sequences
    def __getitem__(self,key):
        if type(key) is int:
            try:
                allseqs = list(self.sequences.keys())
                return self.sequences[allseqs[key]]
            except:
                raise Exception("You tried to access sequence num",key,"but there's only",len(self.sequences.keys()))
        else:
            return self.sequences[key]
    def __iter__(self):
        return self.sequences.__iter__()
    def __repr__(self):
        ret=''
        for s in self.sequences:
            ret += '%s\t%s\n'%(s,self.sequences[s])
        return ret
    def read_sequences(self,fasta_fn,name_map=None):
        code = None
        seq = None
        with open(fasta_fn,'r') as fh:
            for (num, line) in enumerate(fh):
                if line.startswith('>'):
                    if seq is not None:
                        self.sequences[code] = seq.strip('*')
                    code = line.rstrip()[1:]
                    if name_map is not None:
                        try:
                            code = name_map[code]
                        except:
                            pass
                    seq = ''
                else:
                    line = line.rstrip()
                    if line: # Skip blank lines
                        if seq is None:
                            raise Exception( \
    "Found FASTA sequence before first header at line %d: %s" % (num + 1, line))
                        seq += line
        if seq is not None:
            self.sequences[code] = seq.strip('*')


class TopologyReader(object):
    """Automatically setup Sytem and Degrees of Freedom with a formatted text file.
    The file is read in and each part of the topology is stored as a
    ComponentTopology object for input into IMP::pmi::macros::BuildSystem.
    The topology file should be in a simple pipe-delimited format:
    @code{.txt}
|molecule_name|color|fasta_fn|fasta_id|pdb_fn|chain|residue_range|pdb_offset|bead_size|em_residues_per_gaussian|rigid_body|super_rigid_body|chain_of_super_rigid_bodies|flags|
|Rpb1   |blue   |1WCM.fasta|1WCM:A|1WCM.pdb|A|1,1140   |0|10|0|1|1,3|1||
|Rpb1   |blue   |1WCM.fasta|1WCM:A|1WCM.pdb|A|1141,1274|0|10|0|2|1,3|1||
|Rpb1   |blue   |1WCM.fasta|1WCM:A|1WCM.pdb|A|1275,END |0|10|0|3|1,3|1||
|Rpb2   |red    |1WCM.fasta|1WCM:B|1WCM.pdb|B|all      |0|10|0|4|2,3|2||
|Rpb2.1 |green  |1WCM.fasta|1WCM:B|1WCM.pdb|B|all      |0|10|0|4|2,3|2||

    @endcode

    These are the fields you can enter:
    - `component_name`: Name of the component (chain). Serves as the parent
      hierarchy for this structure.
    - `color`: The color used in the output RMF file. Uses chimera names or R,G,B values
    - `fasta_fn`: Name of FASTA file containing this component.
    - `fasta_id`: String found in FASTA sequence header line. The sequence read
      from the file is assumed to be a protein sequence. If it should instead
      be treated as RNA or DNA, add an ',RNA' or ',DNA' suffix. For example,
      a `fasta_id` of 'myseq,RNA' will read the sequence 'myseq' from the
      FASTA file and treat it as RNA.
    - `pdb_fn`: Name of PDB file with coordinates (if available).
       If left empty, will set up as BEADS (you can also specify "BEADS")
       Can also write "IDEAL_HELIX".
    - `chain`: Chain ID of this domain in the PDB file.
    - `residue_range`: Comma delimited pair defining range.
       Can leave empty or use 'all' for entire sequence from PDB file.
       The second item in the pair can be END to select the last residue in the
       PDB chain.
    - `pdb_offset`: Offset to sync PDB residue numbering with FASTA numbering.
    - `bead_size`: The size (in residues) of beads used to model areas not
      covered by PDB coordinates. These will be automatically built.
    - `em_residues`: The number of Gaussians used to model the electron
      density of this domain. Set to zero if no EM fitting will be done.
      The GMM files will be written to <gmm_dir>/<component_name>_<em_res>.txt (and .mrc)
    - `rigid_body`: Leave empty if this object is not in a rigid body.
       Otherwise, this is a number corresponding to the rigid body containing
       this object. The number itself is just used for grouping things.
    - `super_rigid_body`: Like a rigid_body, except things are only occasionally rigid
    - `chain_of_super_rigid_bodies` For a polymer, create SRBs from groups.
    - `flags` additional flags for advanced options
    \note All filenames are relative to the paths specified in the constructor.

    """
    def __init__(self,
                 topology_file,
                 pdb_dir='./',
                 fasta_dir='./',
                 gmm_dir='./'):
        """Constructor.
        @param topology_file Pipe-delimited file specifying the topology
        @param pdb_dir Relative path to the pdb directory
        @param fasta_dir Relative path to the fasta directory
        @param gmm_dir Relative path to the GMM directory
        """
        self.topology_file = topology_file
        self.molecules = IMP.pmi.tools.OrderedDict() # key=molname, value=TempMolecule
        self.pdb_dir = pdb_dir
        self.fasta_dir = fasta_dir
        self.gmm_dir = gmm_dir
        self._components = self.read(topology_file)

    # Preserve old self.component_list for backwards compatibility
    @IMP.deprecated_method("2.7",
                       "Use 'get_components()' instead of 'component_list'.")
    def __get_component_list(self): return self._components
    component_list = property(__get_component_list)

    def write_topology_file(self,outfile):
        with open(outfile, "w") as f:
            f.write("|molecule_name|color|fasta_fn|fasta_id|pdb_fn|chain|"
                    "residue_range|pdb_offset|bead_size|em_residues_per_gaussian|"
                    "rigid_body|super_rigid_body|chain_of_super_rigid_bodies|\n")
            for c in self._components:
                output = c.get_str()+'\n'
                f.write(output)
        return outfile

    def get_components(self, topology_list = "all"):
        """ Return list of ComponentTopologies for selected components
        @param topology_list List of indices to return"""
        if topology_list == "all":
            topologies = self._components
        else:
            topologies=[]
            for i in topology_list:
                topologies.append(self._components[i])
        return topologies

    def get_molecules(self):
        return self.molecules

    def read(self, topology_file, append=False):
        """Read system components from topology file. append=False will erase
        current topology and overwrite with new
        """
        is_topology = False
        is_directories = False
        linenum = 1
        if append==False:
            self._components=[]

        with open(topology_file) as infile:
            for line in infile:
                if line.lstrip()=="" or line[0]=="#":
                    continue
                elif line.split('|')[1].strip() in ("molecule_name"):
                    is_topology=True
                    is_directories = False
                    old_format = False
                    continue
                elif line.split('|')[1] == "component_name":
                    is_topology = True
                    IMP.handle_use_deprecated(
                          "Old-style topology format (using "
                          "|component_name|) is deprecated. Please switch to "
                          "the new-style format (using |molecule_name|)\n")
                    old_format = True
                    is_directories = False
                    continue
                elif line.split('|')[1] == "directories":
                    IMP.handle_use_deprecated(
                          "Setting directories in the topology file "
                          "is deprecated. Please do so through the "
                          "TopologyReader constructor. Note that new-style "
                          "paths are relative to the current working "
                          "directory, not the topology file.\n")
                    is_directories = True
                elif is_directories:
                    fields = line.split('|')
                    setattr(self, fields[1],
                            IMP.get_relative_path(topology_file, fields[2]))
                if is_topology:
                    new_component = self._parse_line(line, linenum, old_format)
                    self._components.append(new_component)
                    linenum += 1
        return self._components

    def _parse_line(self, component_line, linenum, old_format):
        """Parse a line of topology values and matches them to their key.
        Checks each value for correct syntax
        Returns a list of Component objects
        fields:
        """
        c = _Component()
        values = [s.strip() for s in component_line.split('|')]
        errors = []

        ### Required fields
        if old_format:
            c.molname = values[1]
            c.copyname = ''
            c._domain_name = values[2]
            c.color = 'blue'
        else:
            names = values[1].split('.')
            if len(names)==1:
                c.molname = names[0]
                c.copyname = ''
            elif len(names)==2:
                c.molname = names[0]
                c.copyname = names[1]
            else:
                c.molname = names[0]
                c.copyname = names[1]
                errors.append("Molecule name should be <molecule.copyID>")
                errors.append("For component %s line %d " % (c.molname,linenum))
            c._domain_name = c.molname + '.' + c.copyname
            colorfields = values[2].split(',')
            if len(colorfields)==3:
                c.color = [float(x) for x in colorfields]
                if any([x>1 for x in c.color]):
                    c.color=[x/255 for x in c.color]
            else:
                c.color = values[2]
        c._orig_fasta_file = values[3]
        c.fasta_file = values[3]
        fasta_field = values[4].split(",")
        c.fasta_id  = fasta_field[0]
        c.fasta_flag = None
        if len(fasta_field) > 1:
            c.fasta_flag = fasta_field[1]
        c._orig_pdb_input = values[5]
        pdb_input = values[5]
        tmp_chain = values[6]
        rr = values[7]
        offset = values[8]
        bead_size = values[9]
        emg = values[10]
        if old_format:
            rbs = srbs = csrbs = ''
        else:
            rbs = values[11]
            srbs = values[12]
            csrbs = values[13]

        if c.molname not in self.molecules:
            self.molecules[c.molname] = _TempMolecule(c)
        else:
            # COPY OR DOMAIN
            c._orig_fasta_file = self.molecules[c.molname].orig_component.fasta_file
            c.fasta_id = self.molecules[c.molname].orig_component.fasta_id
            self.molecules[c.molname].add_component(c,c.copyname)

        # now cleanup input
        c.fasta_file = os.path.join(self.fasta_dir,c._orig_fasta_file)
        if pdb_input=="":
            errors.append("PDB must have BEADS, IDEAL_HELIX, or filename")
            errors.append("For component %s line %d is not correct"
                          "|%s| was given." % (c.molname,linenum,pdb_input))
        elif pdb_input in ("IDEAL_HELIX","BEADS"):
            c.pdb_file = pdb_input
        else:
            c.pdb_file = os.path.join(self.pdb_dir,pdb_input)

            # PDB chain must be one or two characters
            if len(tmp_chain)==1 or len(tmp_chain)==2:
                c.chain = tmp_chain
            else:
                errors.append("PDB Chain identifier must be one or two characters.")
                errors.append("For component %s line %d is not correct"
                              "|%s| was given." % (c.molname,linenum,tmp_chain))

        ### Optional fields
        # Residue Range
        if rr.strip()=='all' or str(rr)=="":
            c.residue_range = None
        elif len(rr.split(','))==2 and self._is_int(rr.split(',')[0]) and (self._is_int(rr.split(',')[1]) or rr.split(',')[1] == 'END'):
            # Make sure that is residue range is given, there are only two values and they are integers
            c.residue_range = (int(rr.split(',')[0]), rr.split(',')[1])
            if c.residue_range[1] != 'END':
                c.residue_range = (c.residue_range[0], int(c.residue_range[1]))
            # Old format used -1 for the last residue
            if old_format and c.residue_range[1] == -1:
                c.residue_range = (c.residue_range[0], 'END')
        else:
            errors.append("Residue Range format for component %s line %d is not correct" % (c.molname, linenum))
            errors.append("Correct syntax is two comma separated integers:  |start_res, end_res|. end_res can also be END to select the last residue in the chain. |%s| was given." % rr)
            errors.append("To select all residues, indicate |\"all\"|")

        # PDB Offset
        if self._is_int(offset):
            c.pdb_offset=int(offset)
        elif len(offset)==0:
            c.pdb_offset = 0
        else:
            errors.append("PDB Offset format for component %s line %d is not correct" % (c.molname, linenum))
            errors.append("The value must be a single integer. |%s| was given." % offset)

        # Bead Size
        if self._is_int(bead_size):
            c.bead_size=int(bead_size)
        elif len(bead_size)==0:
            c.bead_size = 0
        else:
            errors.append("Bead Size format for component %s line %d is not correct" % (c.molname, linenum))
            errors.append("The value must be a single integer. |%s| was given." % bead_size)

        # EM Residues Per Gaussian
        if self._is_int(emg):
            if int(emg) > 0:
                c.density_prefix = os.path.join(self.gmm_dir,c.get_unique_name())
                c.gmm_file = c.density_prefix+'.txt'
                c.mrc_file = c.density_prefix+'.gmm'

                c.em_residues_per_gaussian=int(emg)
            else:
                c.em_residues_per_gaussian = 0
        elif len(emg)==0:
            c.em_residues_per_gaussian = 0
        else:
            errors.append("em_residues_per_gaussian format for component "
                          "%s line %d is not correct" % (c.molname, linenum))
            errors.append("The value must be a single integer. |%s| was given." % emg)

        # rigid bodies
        if len(rbs)>0:
            if not self._is_int(rbs):
                errors.append("rigid bodies format for component "
                              "%s line %d is not correct" % (c.molname, linenum))
                errors.append("Each RB must be a single integer, or empty. "
                              "|%s| was given." % rbs)
            c.rigid_body = int(rbs)

        # super rigid bodies
        if len(srbs)>0:
            srbs = srbs.split(',')
            for i in srbs:
                if not self._is_int(i):
                    errors.append("super rigid bodies format for component "
                                  "%s line %d is not correct" % (c.molname, linenum))
                    errors.append("Each SRB must be a single integer. |%s| was given." % srbs)
            c.super_rigid_bodies = srbs

        # chain of super rigid bodies
        if len(csrbs)>0:
            if not self._is_int(csrbs):
                errors.append("em_residues_per_gaussian format for component "
                              "%s line %d is not correct" % (c.molname, linenum))
                errors.append("Each CSRB must be a single integer. |%s| was given." % csrbs)
            c.chain_of_super_rigid_bodies = csrbs

        # done
        if errors:
            raise ValueError("Fix Topology File syntax errors and rerun: " \
                             + "\n".join(errors))
        else:
            return c


    def set_gmm_dir(self,gmm_dir):
        """Change the GMM dir"""
        self.gmm_dir = gmm_dir
        for c in self._components:
            c.gmm_file = os.path.join(self.gmm_dir,c.get_unique_name()+".txt")
            c.mrc_file = os.path.join(self.gmm_dir,c.get_unique_name()+".mrc")
            print('new gmm',c.gmm_file)

    def set_pdb_dir(self,pdb_dir):
        """Change the PDB dir"""
        self.pdb_dir = pdb_dir
        for c in self._components:
            if not c._orig_pdb_input in ("","None","IDEAL_HELIX","BEADS"):
                c.pdb_file = os.path.join(self.pdb_dir,c._orig_pdb_input)

    def set_fasta_dir(self,fasta_dir):
        """Change the FASTA dir"""
        self.fasta_dir = fasta_dir
        for c in self._components:
            c.fasta_file = os.path.join(self.fasta_dir,c._orig_fasta_file)

    def _is_int(self, s):
       # is this string an integer?
        try:
            float(s)
            return float(s).is_integer()
        except ValueError:
            return False

    def get_rigid_bodies(self):
        """Return list of lists of rigid bodies (as domain name)"""
        rbl = defaultdict(list)
        for c in self._components:
            if c.rigid_body:
                rbl[c.rigid_body].append(c.get_unique_name())
        return rbl.values()

    def get_super_rigid_bodies(self):
        """Return list of lists of super rigid bodies (as domain name)"""
        rbl = defaultdict(list)
        for c in self._components:
            for rbnum in c.super_rigid_bodies:
                rbl[rbnum].append(c.get_unique_name())
        return rbl.values()

    def get_chains_of_super_rigid_bodies(self):
        """Return list of lists of chains of super rigid bodies (as domain name)"""
        rbl = defaultdict(list)
        for c in self._components:
            for rbnum in c.chain_of_super_rigid_bodies:
                rbl[rbnum].append(c.get_unique_name())
        return rbl.values()

class _TempMolecule(object):
    """Store the Components and any requests for copies"""
    def __init__(self,init_c):
        self.molname = init_c.molname
         # key=copy ID, value = list of domains
        self.domains = IMP.pmi.tools.OrderedDefaultDict(list)
        self.add_component(init_c,init_c.copyname)
        self.orig_copyname = init_c.copyname
        self.orig_component = self.domains[init_c.copyname][0]
    def add_component(self,component,copy_id):
        self.domains[copy_id].append(component)
        component.domainnum = len(self.domains[copy_id])-1
    def __repr__(self):
        return ','.join('%s:%i'%(k,len(self.domains[k])) for k in self.domains)

class _Component(object):
    """Stores the components required to build a standard IMP hierarchy
    using IMP.pmi.BuildModel()
    """
    def __init__(self):
        self.molname = None
        self.copyname = None
        self.domainnum = 0
        self.fasta_file = None
        self._orig_fasta_file = None
        self.fasta_id = None
        self.fasta_flag = None
        self.pdb_file = None
        self._orig_pdb_input = None
        self.chain = None
        self.residue_range = None
        self.pdb_offset = 0
        self.bead_size = 10
        self.em_residues_per_gaussian = 0
        self.gmm_file = ''
        self.mrc_file = ''
        self.density_prefix = ''
        self.color = 0.1
        self.rigid_body = None
        self.super_rigid_bodies = []
        self.chain_of_super_rigid_bodies = []

    def _l2s(self,l):
        return ",".join("%s" % x for x in l)

    def __repr__(self):
        return self.get_str()

    def get_unique_name(self):
        return "%s.%s.%i"%(self.molname,self.copyname,self.domainnum)

    def get_str(self):
        res_range = self.residue_range
        if self.residue_range is None:
            res_range = []
        name = self.molname
        if self.copyname!='':
            name += '.'+self.copyname
        if self.chain is None:
            chain = ' '
        else:
            chain = self.chain
        color=self.color
        if isinstance(color, list):
            color=','.join([str(x) for x in color])
        a= '|'+'|'.join([name,color,self._orig_fasta_file,self.fasta_id,
                         self._orig_pdb_input,chain,self._l2s(list(res_range)),
                             str(self.pdb_offset),str(self.bead_size),
                             str(self.em_residues_per_gaussian),
                             str(self.rigid_body) if self.rigid_body else '',
                             self._l2s(self.super_rigid_bodies),
                             self._l2s(self.chain_of_super_rigid_bodies)])+'|'
        return a

    # Preserve old self.name for backwards compatibility
    @IMP.deprecated_method("2.7", "Use 'molname' instead of 'name'.")
    def __get_name(self): return self.molname
    name = property(__get_name)

    # Preserve old self.domain_name for backwards compatibility
    @IMP.deprecated_method("2.7",
                           "Use 'get_unique_name()' instead of 'domain_name'.")
    def __get_domain_name(self): return self._domain_name
    domain_name = property(__get_domain_name)


class PMIMoleculeHierarchy(IMP.atom.Molecule):
    '''Extends the functionality of IMP.atom.Molecule'''

    def __init__(self,hierarchy):
        IMP.atom.Molecule.__init__(self,hierarchy)

    def get_state_index(self):
        state = self.get_parent()
        return IMP.atom.State(state).get_state_index()

    def get_copy_index(self):
        return IMP.atom.Copy(self).get_copy_index()

    def get_extended_name(self):
        return self.get_name()+"."+\
               str(self.get_copy_index())+\
               "."+str(self.get_state_index())

    def get_sequence(self):
        return IMP.atom.Chain(self).get_sequence()

    def get_residue_indexes(self):
        return IMP.pmi.tools.get_residue_indexes(self)

    def get_residue_segments(self):
        return IMP.pmi.tools.Segments(self.get_residue_indexes())

    def get_chain_id(self):
        return IMP.atom.Chain(self).get_id()

    def __repr__(self):
        s='PMIMoleculeHierarchy '
        s+=self.get_name()
        s+=" "+"Copy  "+str(IMP.atom.Copy(self).get_copy_index())
        s+=" "+"State "+str(self.get_state_index())
        s+=" "+"N residues "+str(len(self.get_sequence()))
        return s
