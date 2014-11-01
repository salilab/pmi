#!/usr/bin/env python

"""@namespace IMP.pmi.topology_io
   Utilities for interpreting and analyzing sequence data
"""

import IMP
import IMP.pmi
import csv

class TopologyReader(object):
    '''
    This class reads in a standard pipe-delimited PMI topology file
    and stores the items as a topology class for input into IMP.pmi.autobuild_model()
    '''
    def __init__(self, topology_file):
        is_defaults=False
        is_topology=False
        defaults_dict={}
        linenum=1
        self.component_list=[]
        self.defaults={'bead_size'                : 10,
                       'residue_range'            : 'all',
                       'pdb_offset'               : 0,
                       'em_residues_per_gaussian' : 0};


        with open(topology_file) as infile:
            for line in infile:

                if line.lstrip()=="":
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

    def create_component_topology(self, component_line, topology_fields, defaults, linenum, color="0.1"):

    #Reads a grid of topology values and matches them to their key.
    #Checks each value for correct syntax
    #Returns a list of ComponentTopology objects

        fields=topology_fields.split('|')
        values=component_line.split('|')
        c=ComponentTopology()   
        no_error=True
    ##### Required fields
        c.name          = values[fields.index("component_name")]
        c.domain_name   = values[fields.index("domain_name")]
        c.fasta_file    = defaults['fasta_dir'] + values[fields.index("fasta_fn")]
        c.fasta_id      = values[fields.index("fasta_id")]
        c.pdb_file      = defaults['pdb_dir'] + values[fields.index("pdb_fn")]
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
                    c.read_gmm_files=True
                    c.gmm_file=defaults['gmm_dir'] + c.domain_name.strip() + ".txt"
                    c.mrc_file=defaults['gmm_dir'] + c.domain_name.strip() + ".mrc"
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
        self.read_gmm_files=None
        self.gmm_file=None
        self.mrc_file=None
        self.read_em_files=None
        self.color=None


