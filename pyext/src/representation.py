#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.display
import IMP.pmi

class SimplifiedModel():
#Peter Cimermancic and Riccardo Pellarin
    '''
    This class creates the molecular hierarchies, representation, 
    sequence connectivity for the various involved proteins and 
    nucleic acid macromolecules:
    
    Create a protein, DNA or RNA, represent it as a set of connected balls of appropriate
    radii and number of residues, pdb at given resolution(s), or ideal helices.
    
    To initialize the class:
    
    m                the model
    upperharmonic    Bool This flag uses either harmonic (False) 
                     or upperharmonic (True) in the intra-pair 
                     connectivity restraint. Default is True.
    disorderedlength Bool This flag uses either disordered length 
                     calculated for random coil peptides (True) or zero 
                     surface-to-surface distance between beads (False)
                     as optimal distance for the sequence connectivity 
                     restraint. Default is False.
    
    How to use the SimplifiedModel class (typical use):
    

    
    m = IMP.Model()
    simo = representation.SimplifiedModel(m,upperharmonic=
                               True,disorderedlength=False)
    
    simo.add_component_name("prot1")
    simo.add_component_beads("prot1",[(1,31)],colors=[0.0])
    simo.add_component_pdb("prot1",'prot1.1.pdb', "A",
                 resolutions=[1,10], color=0.0, offset=-39)
    simo.add_component_beads("prot1",[(131,145)],colors=[0.0])
    simo.add_component_pdb("prot1",'prot1.2.pdb', "A",
                             resolutions=[1,10], color=0.0)
    simo.add_component_beads("prot1",[(626,635)],colors=[0.0])
 
    simo.setup_component_sequence_connectivity("prot1")    
    simo.set_rigid_bodies([("prot1",(30,130))])
    simo.set_rigid_bodies([("prot1",(146,625))]) 
    simo.set_floppy_bodies()
    simo.shuffle_configuration()
    
    '''

    def __init__(self,m,upperharmonic=True,disorderedlength=False):
        global random, itemgetter,tools,nrrand,array,imprmf,RMF,sqrt,output
        import random
        from math import sqrt as sqrt
        from operator import itemgetter
        import IMP.pmi.tools as tools
        import IMP.pmi.output as output
        from numpy.random import rand as nrrand
        from numpy import array
        import IMP.rmf as imprmf
        import RMF

        
        # this flag uses either harmonic (False) or upperharmonic (True)
        # in the intra-pair connectivity restraint. Harmonic is used whe you want to
        # remove the intra-ev term from energy calculations, e.g.:
        # upperharmonic=False
        # ip=simo.get_connected_intra_pairs()
        # ev.add_excluded_particle_pairs(ip)
        
        self.upperharmonic=upperharmonic
        self.disorderedlength=disorderedlength
        self.rigid_bodies=[]
        self.floppy_bodies=[]
        self.super_rigid_bodies=[]
        self.output_level="low"
        self.label="None"
  
        self.maxtrans_rb=0.15
        self.maxrot_rb=0.03
        self.maxtrans_srb=1.0
        self.maxrot_srb=0.025
        self.maxtrans_fb=0.15
        self.resolution=10.0
        self.bblenght=100.0
        self.kappa=100.0
        self.m = m
        
        self.unmodeledregions_cr_dict={}
        self.sortedsegments_cr_dict={}
        self.prot=IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.m))
        self.connected_intra_pairs=[]
        self.hier_dict={}
        self.color_dict={}
        self.sequence_dict={}
        self.hier_geometry_pairs={}
        self.elements={}
        self.linker_restraints=IMP.RestraintSet("linker_restraints")
        self.linker_restraints_dict={}
        self.threetoone={'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
                         'CYS':'C','GLU':'E','GLN':'Q','GLY':'G',
                         'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
                         'MET':'M','PHE':'F','PRO':'P','SER':'S',
                         'THR':'T','TRP':'W','TYR':'Y','VAL':'V','UNK':'X'}
        self.onetothree = {v:k for k, v in self.threetoone.items()}
        self.residuenamekey = IMP.kernel.StringKey("ResidueName")

    def shuffle_configuration(self,bounding_box_length=300.,translate=True):
        "shuffle configuration, used to restart the optimization"
        "it only works if rigid bodies were initialized"
        if len(self.rigid_bodies)==0:
            print "MultipleStates: rigid bodies were not intialized"
        hbbl=bounding_box_length/2
        if 1:
            ub = IMP.algebra.Vector3D(-hbbl,-hbbl,-hbbl)
            lb = IMP.algebra.Vector3D( hbbl, hbbl, hbbl)
            bb = IMP.algebra.BoundingBox3D(ub, lb)
            for rb in self.rigid_bodies:

                if translate==True: translation = IMP.algebra.get_random_vector_in(bb)
                else: translation = (rb.get_x(), rb.get_y(), rb.get_z())
                rotation = IMP.algebra.get_random_rotation_3d()
                transformation = IMP.algebra.Transformation3D(rotation, translation)
                rb.set_reference_frame(IMP.algebra.ReferenceFrame3D(transformation))
            for fb in self.floppy_bodies:
                translation = IMP.algebra.get_random_vector_in(bb)
                IMP.core.XYZ(fb).set_coordinates(translation)


    def add_component_name(self,name,color=None):
        protein_h = IMP.atom.Molecule.setup_particle(IMP.Particle(self.m))
        protein_h.set_name(name)
        self.hier_dict[name]=protein_h
        self.prot.add_child(protein_h)
        self.color_dict[name]=color
        self.elements[name]=[]
    
    def get_component_names(self):
        return self.hier_dict.keys()
    
    def add_component_sequence(self,name,filename,format="FASTA"):
        from Bio import SeqIO
        handle = open(filename, "rU")
        record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        handle.close()
        length=len(record_dict[name].seq)
        self.sequence_dict[name]=str(record_dict[name].seq)
        self.elements[name].append((length,length," ","end"))

    def add_pdb_and_intervening_beads(self,name,pdbname,chain,resolutions,resrange,beadsize,
                                      color=None,pdbresrange=None,offset=0,show=False,isnucleicacid=False,
                                      attachbeads=False):
        outhiers=[]

        if color==None:
           color=self.color_dict[name]
        
        #get the initial and end residues of the pdb
        t=IMP.atom.read_pdb( pdbname, self.m, 
        IMP.atom.AndPDBSelector(IMP.atom.ChainPDBSelector(chain),IMP.atom.ATOMPDBSelector()))

        #find start and end indexes
        
        start = IMP.atom.Residue(t.get_children()[0].get_children()[0]).get_index()
        end   = IMP.atom.Residue(t.get_children()[0].get_children()[-1]).get_index()

        if pdbresrange!=None:
           if pdbresrange[0]>start: start=pdbresrange[0]
           if pdbresrange[1]<end:   end=pdbresrange[1]  
       
        start=start+offset
        end  =end+offset    

        #attach floppy beads to the start or end of PDB structure
        initcoords=[None,None]
        if attachbeads==True:
            firstca= IMP.core.XYZ(t.get_children()[0].get_children()[0].get_children()[0]).get_coordinates()
            lastca=  IMP.core.XYZ(t.get_children()[0].get_children()[-1].get_children()[0]).get_coordinates()     
            initcoords=[firstca,lastca]

        #add pre-beads        
        for i in range(resrange[0],start-1,beadsize)[0:-1]:
            outhiers+=self.add_component_beads(name,[(i,i+beadsize-1)],colors=[color],incoord=initcoords[0])
        
        if resrange[0]<start-1:
           j=range(resrange[0],start-1,beadsize)[-1]
           outhiers+=self.add_component_beads(name,[(j,start-1)],colors=[color],incoord=initcoords[0])
        
        outhiers+=self.add_component_pdb(name,pdbname,chain,resolutions=resolutions, color=color,
                               resrange=resrange,offset=offset,isnucleicacid=isnucleicacid)

        #add after-beads
        for i in range(end+1,resrange[1],beadsize)[0:-1]:
            print "adding sphere from %d to %d , beadsize %d" % (i,i+beadsize-1,beadsize)
            outhiers+=self.add_component_beads(name,[(i,i+beadsize-1)],colors=[color],incoord=initcoords[1])
        
        if end+1<resrange[1]:
           j=range(end+1,resrange[1],beadsize)[-1]
           self.add_component_beads(name,[(j,resrange[1])],colors=[color],incoord=initcoords[1])
        
        #IMP.atom.show_molecular_hierarchy(self.hier_dict[name])
        return outhiers



    def add_component_pdb(self,name,pdbname,chain,resolutions,color=None,resrange=None,offset=0,
                                   cacenters=False,show=False,isnucleicacid=False,readnonwateratoms=False):
                                   
        '''
        resrange specify the residue range to extract from the pdb
        it is a tuple (beg,end). If not specified, it takes all residues belonging to
         the specified chain.
         
        chain can be either a string (eg, A, B, C) or an integer (0,1,2) in case you want
        to get the corresponding chain number in the pdb.
        '''
        
        if color==None:
           color=self.color_dict[name]         
        protein_h=self.hier_dict[name]
        outhiers=[]

        if type(chain)==str:
           if not readnonwateratoms:
              t=IMP.atom.read_pdb( pdbname, self.m, 
              IMP.atom.AndPDBSelector(IMP.atom.ChainPDBSelector(chain),IMP.atom.ATOMPDBSelector()))
           else:
              t=IMP.atom.read_pdb( pdbname, self.m, 
              IMP.atom.AndPDBSelector(IMP.atom.ChainPDBSelector(chain),IMP.atom.NonWaterPDBSelector()))  
           #get the first and last residue
           start = IMP.atom.Residue(t.get_children()[0].get_children()[0]).get_index()
           end   = IMP.atom.Residue(t.get_children()[0].get_children()[-1]).get_index()
           c=IMP.atom.Chain(IMP.atom.get_by_type(t, IMP.atom.CHAIN_TYPE)[0])
               
        elif type(chain)==int:
           if not readnonwateratoms:        
              s=IMP.atom.read_pdb( pdbname, self.m, IMP.atom.ATOMPDBSelector())
           else:
              s=IMP.atom.read_pdb( pdbname, self.m, IMP.atom.NonWaterPDBSelector())        
           t=IMP.atom.Chain(IMP.atom.get_by_type(s, IMP.atom.CHAIN_TYPE)[chain])
           #get the first and last residue
           start = IMP.atom.Residue(t.get_children()[0]).get_index()
           end   = IMP.atom.Residue(t.get_children()[-1]).get_index()           
           c=t           
           chain=t.get_id()
           del s,t
           
        
        if resrange!=None:
           if resrange[0]>start: start=resrange[0]
           if resrange[1]<end:   end=resrange[1] 
        
        
        if not isnucleicacid:
           #do what you have to do for proteins
           sel=IMP.atom.Selection(c,residue_indexes=range(start,end+1),atom_type=IMP.atom.AT_CA)

        else:
           #do what you have to do for nucleic-acids
           sel=IMP.atom.Selection(c,residue_indexes=range(start,end+1),atom_type=IMP.atom.AT_P)
           
         
        ps=sel.get_selected_particles()
        c0=IMP.atom.Chain.setup_particle(IMP.Particle(self.m),"X")

        for p in ps:
            par=IMP.atom.Atom(p).get_parent()
            ri=IMP.atom.Residue(par).get_index()
            IMP.atom.Residue(par).set_index(ri+offset)
            c0.add_child(par)
        start=start+offset
        end=end+offset
        
        self.elements[name].append((start,end,pdbname.split("/")[-1]+":"+chain,"pdb"))
        
        for r in resolutions:
            s=IMP.atom.create_simplified_along_backbone(c0, r)
                           
            chil=s.get_children()
            s0=IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.m))
            
            s0.set_name(name+'_%i-%i_pdb' % (start,end)+"_Res:"+str(r))
            for ch in chil: s0.add_child(ch)            
            protein_h.add_child(s0)
            outhiers+=[s0]
            del s
            for prt in IMP.atom.get_leaves(s0):
                ri=IMP.atom.Fragment(prt).get_residue_indexes()
                first=ri[0]
                last=ri[-1]
                if first==last:
                   prt.set_name(name+'_%i_pdb' % (first))
                else:                   
                   prt.set_name(name+'_%i-%i_pdb' % (first,last))
                
                if r==1:
                   #all the code below is to set the appropriate residue name

                   if not isnucleicacid:
                      sel=IMP.atom.Selection(c0,residue_index=first, atom_type=IMP.atom.AT_CA)
                   else:
                      sel=IMP.atom.Selection(c0,residue_index=first, atom_type=IMP.atom.AT_P)
                   
                   p=sel.get_selected_particles()[0]
                   rtobject=IMP.atom.Residue(IMP.atom.Atom(p).get_parent()).get_residue_type()
                   #the problem of the following function is that it does not 
                   #cover nucleic acid codes, therefore you'll get UNK for RNA and DNA pdbs.
                   rt=self.onetothree[IMP.atom.get_one_letter_code(rtobject)]
                   if cacenters:
                       #this specially sets the radii and centers on the 
                       #cartesian coordinates of the ca atoms when the resolution is 1
                       #beacuse simplified along backbone gets the COM positions
                       
                       xyz=IMP.core.XYZ(p).get_coordinates()
                       IMP.core.XYZ(prt).set_coordinates(xyz)
                       
                       try:
                          vol=IMP.atom.get_volume_from_residue_type(rtobject)
                          radius=IMP.algebra.get_ball_radius_from_volume_3d(vol)
                       except IMP.base.ValueException:
                          radius=IMP.core.XYZR(prt).get_radius()
                       IMP.core.XYZR(prt).set_radius(radius)
                    
                   #setting the residue type to the particle
                   
                   if name in self.sequence_dict:
                      try:
                         rt_final=self.onetothree[self.sequence_dict[name][first-1]]
                      except IndexError:
                         print "add_component_pdb> WARNING: Name: %s residue number: %d the input fasta sequence is shorter \
than the one read from the pdb, setting residue type from pdb. Check the fasta file or the pdb offset" % (name, first)
                         rt_final=rt
                   else:
                      rt_final=rt
                      
                   if rt_final!=rt:
                      print "add_component_pdb> WARNING: Name: %s residue number: %d Expected \
residue: %s residue found in pdb: %s residue types don't match between \
sequence and pdb, setting it to the one in the sequence, check \
the pdb offset" % (name,first,rt_final,rt)
                      
                   #I'm using particle attributes rather that Residue decorators
                   #because I don't want to mess up the hierarchy levels
                   prt.add_attribute(self.residuenamekey, rt_final)
                                      
                radius=IMP.core.XYZR(prt).get_radius()
                IMP.pmi.Uncertainty.setup_particle(prt,radius)
                IMP.pmi.Resolution.setup_particle(prt,r)
                
                #setting up color for each particle in the hierarchy, if colors missing in the colors list set it to red
                try:
                    clr=IMP.display.get_rgb_color(color)
                except:
                    colors.append(1.0)
                    clr=IMP.display.get_rgb_color(colors[pdb_part_count])
                IMP.display.Colored.setup_particle(prt,clr)

        if show:
           IMP.atom.show_molecular_hierarchy(protein_h)
        
        return outhiers


    def add_component_ideal_helix(self,name,resolutions,resrange,color=None,show=False):
    
        from math import pi,cos,sin
        
        protein_h=self.hier_dict[name]   
        outhiers=[]
        if color==None:
           color=self.color_dict[name]
        
        start=resrange[0]
        end=resrange[1]
        self.elements[name].append((start,end," ","helix"))
        c0=IMP.atom.Chain.setup_particle(IMP.Particle(self.m),"X")         
        for n,res in enumerate(range(start,end+1)):
            r=IMP.atom.Residue.setup_particle(IMP.Particle(self.m),IMP.atom.ALA,res)
            p=IMP.Particle(self.m)
            d=IMP.core.XYZR.setup_particle(p)
            x=2.3*cos(n*2*pi/3.6)
            y=2.3*sin(n*2*pi/3.6) 
            z=5.4/3.6/2*n*2*pi/3.6   
            d.set_coordinates(IMP.algebra.Vector3D(x,y,z))
            d.set_radius(2.9)
            #print d
            a=IMP.atom.Atom.setup_particle(p,IMP.atom.AT_CA)
            r.add_child(a)
            c0.add_child(r)
            
        for r in resolutions:

            s=IMP.atom.create_simplified_along_backbone(c0, r)
            chil=s.get_children()
            s0=IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.m))
            
            s0.set_name(name+'_%i-%i_helix' % (start,end)+"_Res:"+str(r))
            for ch in chil: s0.add_child(ch)            
            protein_h.add_child(s0)
            outhiers+=[s0]
            del s
            for prt in IMP.atom.get_leaves(s0):
                ri=IMP.atom.Fragment(prt).get_residue_indexes()
                first=ri[0]
                last=ri[-1]
                if first==last:
                   prt.set_name(name+'_%i_helix' % (first))
                else:                   
                   prt.set_name(name+'_%i-%i_helix' % (first,last))
                radius=IMP.core.XYZR(prt).get_radius()
                if r==1: 
                   if name in self.sequence_dict:
                      rt_final=self.onetothree[self.sequence_dict[name][first-1]]
                      rtobject=IMP.atom.ResidueType(rt_final)
                      vol=IMP.atom.get_volume_from_residue_type(rtobject)
                      radius=IMP.algebra.get_ball_radius_from_volume_3d(vol)
                   else:
                      rt_final="X"
                      rtobject=IMP.atom.ResidueType("ALA")
                      vol=IMP.atom.get_volume_from_residue_type(rtobject)
                      radius=IMP.algebra.get_ball_radius_from_volume_3d(vol)
                   prt.add_attribute(self.residuenamekey, rt_final)
                   IMP.core.XYZR(prt).set_radius(radius)

                IMP.pmi.Uncertainty.setup_particle(prt,radius)
                IMP.pmi.Resolution.setup_particle(prt,r)
                #setting up color for each particle in the hierarchy, if colors missing in the colors list set it to red
                try:
                    clr=IMP.display.get_rgb_color(color)
                except:
                    colors.append(1.0)
                    clr=IMP.display.get_rgb_color(colors[pdb_part_count])
                IMP.display.Colored.setup_particle(prt,clr)

        if show:
           IMP.atom.show_molecular_hierarchy(protein_h)
        
        return outhiers 
            

    def add_component_beads(self,name,ds,colors=None,incoord=None):
        from math import pi
        protein_h=self.hier_dict[name]
        outhiers=[]
        if colors==None:
           colors=[self.color_dict[name]]
            
        for n,dss in enumerate(ds):
            ds_frag=(dss[0],dss[1])
            self.elements[name].append((dss[0],dss[1]," ","bead"))
            prt=IMP.Particle(self.m)
            h=IMP.atom.Fragment.setup_particle(prt)
            h.set_residue_indexes(range(ds_frag[0],ds_frag[1]+1))
            
            if ds_frag[0]==ds_frag[1]:
               h.set_name(name+'_%i_bead' % (ds_frag[0]))
            else:
               h.set_name(name+'_%i-%i_bead' % (ds_frag[0],ds_frag[1]))               
            resolution=len(h.get_residue_indexes())
            try:
                clr=IMP.display.get_rgb_color(colors[n])
            except:
                clr=IMP.display.get_rgb_color(colors[0])
            
            IMP.display.Colored.setup_particle(prt,clr)
            
            #decorate particles according to their resolution
            IMP.pmi.Resolution.setup_particle(prt,1000000)
            p=IMP.atom.get_leaves(h)[0]
            IMP.core.XYZR.setup_particle(p)
            ptem=IMP.core.XYZR(p)
            mass =IMP.atom.get_mass_from_number_of_residues(resolution)
            volume=IMP.atom.get_volume_from_mass(mass)
            radius=0.8*(3.0/4.0/pi*volume)**(1.0/3.0)
            ptem.set_radius(radius)
            try:
                if tuple(incoord)!=None: ptem.set_coordinates(incoord)
            except TypeError: pass 
            IMP.pmi.Uncertainty.setup_particle(ptem,radius)
            self.floppy_bodies.append(p)
            outhiers+=[h]
            protein_h.add_child(h)
        
        return outhiers

    def add_component_necklace(self,name,begin,end,length):
        
        outhiers=[]
        #nbeads=len(range(begin,end,length))
        #lastend=range(begin,end,length)[-2]
        #if float(end-lastend+length)<length/2:
        #   length=length+int(float(end-i+length)/(nbeads-1))

        for i in range(begin,end,length)[0:-1]:
           outhiers+=self.add_component_beads(name,[(i,i+length-1)])
        outhiers+=self.add_component_beads(name,[(i+length,end)])
        
        return outhiers
        
 
    def setup_component_geometry(self,name,color=None):
        if color==None:
           color=self.color_dict[name]
        #this function stores all particle pairs
        #ordered by residue number, to be used 
        #to construct backbone traces
        self.hier_geometry_pairs[name]=[]
        protein_h=protein_h=self.hier_dict[name]
        pbr=tools.get_particles_by_resolution(protein_h,1.0)
        
        sortedparticles=[]
        
        for p in pbr:
            startres = IMP.atom.Fragment(p).get_residue_indexes()[0]
            sortedparticles.append((p,startres))
            sortedparticles = sorted(sortedparticles, key=itemgetter(1))
            
        for n in range(len(sortedparticles)-1):
            self.hier_geometry_pairs[name].append((sortedparticles[n][0],sortedparticles[n+1][0],color))
      
    def setup_component_sequence_connectivity(self,name,resolution=10):
        unmodeledregions_cr=IMP.RestraintSet("unmodeledregions")
        sortedsegments_cr=IMP.RestraintSet("sortedsegments")   
         
        protein_h=protein_h=self.hier_dict[name]    
        SortedSegments = []
        pbr=tools.get_particles_by_resolution(protein_h,resolution)
        
        for chl in protein_h.get_children():
            start = IMP.atom.get_leaves(chl)[0]
            end   = IMP.atom.get_leaves(chl)[-1]
            if not start in pbr: continue
            if not end   in pbr: continue
            startres = IMP.atom.Fragment(start).get_residue_indexes()[0]
            endres   = IMP.atom.Fragment(end).get_residue_indexes()[-1]
            SortedSegments.append((chl,startres))
        SortedSegments = sorted(SortedSegments, key=itemgetter(1))

        #connect the particles
        for x in xrange(len(SortedSegments)-1):
            last = IMP.atom.get_leaves(SortedSegments[x][0])[-1]
            first= IMP.atom.get_leaves(SortedSegments[x+1][0])[0]
            nreslast=len(IMP.atom.Fragment(last).get_residue_indexes())
            lastresn=IMP.atom.Fragment(last).get_residue_indexes()[-1]
            nresfirst=len(IMP.atom.Fragment(first).get_residue_indexes())
            firstresn=IMP.atom.Fragment(first).get_residue_indexes()[0]
            
            residuegap=firstresn-lastresn-1
            
            if self.disorderedlength and (nreslast/2+nresfirst/2+residuegap)>20.0 :
               #calculate the distance between the sphere centers using Kohn PNAS 2004               
               optdist=sqrt(5/3)*1.93*(nreslast/2+nresfirst/2+residuegap)**0.6
               #optdist2=sqrt(5/3)*1.93*((nreslast)**0.6+(nresfirst)**0.6)/2
               if self.upperharmonic:
                  hu=IMP.core.HarmonicUpperBound(optdist, self.kappa)
               else:
                  hu=IMP.core.Harmonic(optdist, self.kappa) 
               dps=IMP.core.DistancePairScore(hu)            
            else: #default
               optdist=0.0+residuegap*3.6
               if self.upperharmonic: #default
                  hu=IMP.core.HarmonicUpperBound(optdist, self.kappa)
               else:
                  hu=IMP.core.Harmonic(optdist, self.kappa)               
               dps=IMP.core.SphereDistancePairScore(hu)
            
            
            
            pt0=IMP.atom.Selection(last).get_selected_particles()[0]          
            pt1=IMP.atom.Selection(first).get_selected_particles()[0]           
            r=IMP.core.PairRestraint(dps,IMP.ParticlePair(pt0,pt1))
            
            print "Adding sequence connectivity restraint between", pt0.get_name(), " and ", pt1.get_name()
            sortedsegments_cr.add_restraint(r)
            self.linker_restraints_dict["LinkerRestraint-"+pt0.get_name()+"-"+pt1.get_name()]=r
            self.connected_intra_pairs.append((pt0,pt1))
            self.connected_intra_pairs.append((pt1,pt0))

        self.m.add_restraint(sortedsegments_cr)
        self.m.add_restraint(unmodeledregions_cr)
        self.linker_restraints.add_restraint(sortedsegments_cr)
        self.linker_restraints.add_restraint(unmodeledregions_cr)
        self.sortedsegments_cr_dict[name]=sortedsegments_cr
        self.unmodeledregions_cr_dict[name]=unmodeledregions_cr


    def add_component(self,name, chainnames, length, pdbs, 
                      init_coords=None, simplepdb=1, ds=None, colors=None, 
                      resolutions=None):
        '''
        Using multiple resolution (resolutions!=None)
        Using this representation, a crystallographic structure or a homology model
        can be represented by a set of models, coarse-grained at different levels (giving as input a list of sizes 
        for the beads, in residue number). All particles belonging to the models are constrained in the same rigid body. 
        Restraints constructors have an optional argument that corresponds to the resolution you want to apply it.
        If the resolution of the model is hybrid (i.e., some parts are multi resolution, and some 
        other are single resolution because they are flexible) the restraint picks up the most convenient particle.

        In this way cross-links will be applied directly to individual residues, excluded 
        volume is applied to intermediate resolutions (eg. 10 residues per bead) and 
        connectivity restraint from domain mapping is applied to larger beads (eg. 100 residues per bead), 
        saving a lot of computational time. I've benchmarked it and it is generally faster 
        (from twice, up to 2 orders of magnitude faster), and of course, we are not limited 
        by the usage of a single compromise resolution (eg, 5 residues per beads).
        '''
    
        if ds==None: ds=[]
        if colors==None: colors=[]
        if init_coords==None: init_coords=()
        if resolutions==None: resolutions=[]
        
        
        protein_h = IMP.atom.Molecule.setup_particle(IMP.Particle(self.m))
        protein_h.set_name(name)
        bb=IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-self.bblenght,-self.bblenght,-self.bblenght),
                                      IMP.algebra.Vector3D(self.bblenght, self.bblenght, self.bblenght))

        unmodeledregions_cr=IMP.RestraintSet("unmodeledregions")
        sortedsegments_cr=IMP.RestraintSet("sortedsegments")

        RigiParticles=[]
        FlexParticles=[]

        # work on PDB structures first
        #ds = [] #store regions without structures
        if len(pdbs)>0:
            bounds = []
            for pdb_part_count,pdb in enumerate(pdbs):
                sls=IMP.base.SetLogState(IMP.NONE)
                t=IMP.atom.read_pdb( pdb, self.m, 
                  IMP.atom.AndPDBSelector(IMP.atom.ChainPDBSelector(chainnames[pdb_part_count]), 
                                                                      IMP.atom.ATOMPDBSelector()))
                            
                del sls

                #find start and end indexes
                start = IMP.atom.Residue(t.get_children()[0].get_children()[0]).get_index()
                end   = IMP.atom.Residue(t.get_children()[0].get_children()[-1]).get_index()
                bounds.append(( start,end ))

                #IMP.atom.show_molecular_hierarchy(t)
                c=IMP.atom.Chain(IMP.atom.get_by_type(t, IMP.atom.CHAIN_TYPE)[0])
                if c.get_number_of_children()==0:
                    IMP.atom.show_molecular_hierarchy(t)
                # there is no reason to use all atoms, just approximate the pdb shape instead
                # add coarse level hierarchy
                
                if len(resolutions)==0:
                   #default
                   s=IMP.atom.create_simplified_along_backbone(c, self.resolution/2.0)
                   if simplepdb==1: s=IMP.atom.create_simplified_along_backbone(c, self.resolution/2.0)
                   else: s=IMP.atom.create_simplified_along_backbone(c, 1.)

                   s.set_name(chainnames[pdb_part_count]+str(pdb_part_count))
                   IMP.atom.destroy(t)
                   '''
                   # make the simplified structure rigid
                   if init_coords!=() and rbo==0:
                      print name,pdb
                      rb=IMP.atom.create_rigid_body(s)
                      rb.set_coordinates(init_coords)
                      rb.set_coordinates_are_optimized(True)
                      RigiParticles.append(rb)
                   #if rbo==1: rb.set_coordinates_are_optimized(False)
                   '''
                   #s0.add_child(s)                          
                   #protein_h.add_child(s0)
                
                   protein_h.add_child(s) 
               
                   for prt in IMP.atom.get_leaves(s):
                     IMP.pmi.Resolution.setup_particle(prt,self.resolution)               
                     #setting up color for each particle in the hierarchy, if colors missing in the colors list set it to red
                     try:
                       clr=IMP.display.get_rgb_color(colors[pdb_part_count])
                     except:
                       colors.append(1.0)
                       clr=IMP.display.get_rgb_color(colors[pdb_part_count])
                     IMP.display.Colored.setup_particle(prt,clr)
                   s.set_name(pdb)

                else:
                   #multiple resolutions
                   for r in resolutions:
                     s0=IMP.atom.create_simplified_along_backbone(c, r)
                     s0.set_name(chainnames[pdb_part_count]+str(pdb_part_count)+str(r))            
                     protein_h.add_child(s0)
                     for prt in IMP.atom.get_leaves(s0):
                       IMP.pmi.Resolution.setup_particle(prt,r)
                       #setting up color for each particle in the hierarchy, if colors missing in the colors list set it to red
                       try:
                         clr=IMP.display.get_rgb_color(colors[pdb_part_count])
                       except:
                         colors.append(1.0)
                         clr=IMP.display.get_rgb_color(colors[pdb_part_count])
                       IMP.display.Colored.setup_particle(prt,clr)
                   

            #calculate regions without structrue (un-modelled regions)
            for i,bnd in enumerate(bounds):
                if i==0:
                    if bnd[0]>1:
                       ds.append((1,bnd[0]-1))
                       if len(ds)>len(colors): colors.append(colors[pdb_part_count])
                    if len(bounds)==1:
                        if bnd[1]<length:
                           ds.append((bnd[1]+1,length))                          
                           if len(ds)>len(colors): colors.append(colors[pdb_part_count])
                else:
                    ds.append((bounds[i-1][1]+1, bnd[0]-1))                     
                    if len(ds)>len(colors): colors.append(colors[pdb_part_count])
                    if i==len(bounds)-1:
                        if bnd[1]<length:
                           ds.append((bnd[1]+1,length))                         
                           if len(ds)>len(colors): colors.append(colors[pdb_part_count])

        else:
            ds.append((1,length))



        #work on un-modelled regions
        randomize_coords = lambda c: tuple(10.*(nrrand(3)-0.5)+array(c))

        for n,ds_frag in enumerate(ds):
            if ds_frag[1]-ds_frag[0]==0: ds_frag=(ds_frag[0],ds_frag[1]+1)
            
            '''
            #alternatively: define fragments
            h=IMP.atom.Fragment.setup_particle(IMP.Particle(self.m))
            h.set_residue_indexes(range(ds_frag[0],ds_frag[1]))
            h.set_name(name+'_%i-%i' % (ds_frag[0],ds_frag[1]))
            '''
            
            h=IMP.atom.create_protein(self.m, name+'_%i-%i' % (ds_frag[0],ds_frag[1]), self.resolution, \
                          [ds_frag[0],ds_frag[1]])
                          
            for prt in IMP.atom.get_leaves(h):
                #setting up color for the particle, if colors missing in the colors list set it to red
                try:
                    clr=IMP.display.get_rgb_color(colors[n])
                except:
                    clr=IMP.display.get_rgb_color(1.0)
                IMP.display.Colored.setup_particle(prt,clr)
                
                #decorate particles according to their resolution
                IMP.pmi.Resolution.setup_particle(prt,self.resolution)

                ptem= prt.get_as_xyzr()
                ptem.set_radius(ptem.get_radius()*0.8)

                #if init_coords==(): ptem.set_coordinates((random.uniform(bb[0][0],bb[1][0]),\
                #      random.uniform(bb[0][1],bb[1][1]),random.uniform(bb[0][2],bb[1][2])))
                #else: ptem.set_coordinates(randomize_coords(init_coords))
                if len(pdbs)==0: ptem.set_coordinates((random.uniform(bb[0][0],bb[1][0]),\
                       random.uniform(bb[0][1],bb[1][1]),random.uniform(bb[0][2],bb[1][2])))
                else:
                  if len(resolutions)==0:
                    #default:
                    crd=IMP.atom.get_leaves(s)[0].get_as_xyz().get_coordinates()
                    #comment: why is that done in this way?
                    ptem.set_coordinates(randomize_coords(crd))

                FlexParticles.append(ptem.get_particle())

            Spheres = IMP.atom.get_leaves(h)
            
            for x,prt in enumerate(Spheres):
                if x==len(Spheres)-1: break
                #r=IMP.atom.create_connectivity_restraint([IMP.atom.Selection(prt),\
                #                          IMP.atom.Selection(Spheres[x+1])],-3.0,self.kappa)
                
                if self.disorderedlength:
                   nreslast=len(IMP.atom.Fragment(prt).get_residue_indexes())
                   nresfirst=len(IMP.atom.Fragment(Spheres[x+1]).get_residue_indexes())
                   #calculate the distance between the sphere centers using Kohn PNAS 2004
                   optdist=sqrt(5/3)*1.93*(nreslast/2+nresfirst/2)**0.6
                   #optdist2=sqrt(5/3)*1.93*((nreslast)**0.6+(nresfirst)**0.6)/2               
                   if self.upperharmonic:
                      hu=IMP.core.HarmonicUpperBound(optdist, self.kappa)
                   else:
                      hu=IMP.core.Harmonic(optdist, self.kappa) 
                   dps=IMP.core.DistancePairScore(hu)            
                else: #default
                   optdist=0.0
                   if self.upperharmonic: #default
                      hu=IMP.core.HarmonicUpperBound(optdist, self.kappa)
                   else:
                      hu=IMP.core.Harmonic(optdist, self.kappa)               
                   dps=IMP.core.SphereDistancePairScore(hu)

                pt0=IMP.atom.Selection(prt).get_selected_particles()[0]
                pt1=IMP.atom.Selection(Spheres[x+1]).get_selected_particles()[0]
                #print IMP.core.XYZR(pt0).get_radius()+IMP.core.XYZR(pt1).get_radius()
                r=IMP.core.PairRestraint(dps,IMP.ParticlePair(pt0,pt1))
                unmodeledregions_cr.add_restraint(r)
                self.connected_intra_pairs.append((pt0,pt1))
                self.connected_intra_pairs.append((pt1,pt0))                
                        # only allow the particles to separate by 1 angstrom
                        # self.m.set_maximum_score(r, self.kappa)
            protein_h.add_child(h)


        SortedSegments = []
        
        if len(resolutions)!=0: 
           pbr=tools.get_particles_by_resolution(protein_h,1.0)
        
        for chl in protein_h.get_children():
            
            start = IMP.atom.get_leaves(chl)[0]
            end   = IMP.atom.get_leaves(chl)[-1]
            
            if len(resolutions)!=0: 
               #skip particles            
               if not start in pbr: continue
               if not end   in pbr: continue
            
            startres = IMP.atom.Fragment(start).get_residue_indexes()[0]
            endres   = IMP.atom.Fragment(end).get_residue_indexes()[-1]
            
            
            SortedSegments.append((chl,startres))
        SortedSegments = sorted(SortedSegments, key=itemgetter(1))


        #connect the particles
        for x in xrange(len(SortedSegments)-1):
            last = IMP.atom.get_leaves(SortedSegments[x][0])[-1]
            first= IMP.atom.get_leaves(SortedSegments[x+1][0])[0]
            #r=IMP.atom.create_connectivity_restraint([IMP.atom.Selection(last),\
            #                                          IMP.atom.Selection(first)],-3.0,self.kappa)

            if self.disorderedlength:
               nreslast=len(IMP.atom.Fragment(last).get_residue_indexes())
               nresfirst=len(IMP.atom.Fragment(first).get_residue_indexes())
               #calculate the distance between the sphere centers using Kohn PNAS 2004               
               optdist=sqrt(5/3)*1.93*(nreslast/2+nresfirst/2)**0.6
               #optdist2=sqrt(5/3)*1.93*((nreslast)**0.6+(nresfirst)**0.6)/2
               if self.upperharmonic:
                  hu=IMP.core.HarmonicUpperBound(optdist, self.kappa)
               else:
                  hu=IMP.core.Harmonic(optdist, self.kappa) 
               dps=IMP.core.DistancePairScore(hu)            
            else: #default
               optdist=0.0
               if self.upperharmonic: #default
                  hu=IMP.core.HarmonicUpperBound(optdist, self.kappa)
               else:
                  hu=IMP.core.Harmonic(optdist, self.kappa)               
               dps=IMP.core.SphereDistancePairScore(hu)
            
            pt0=IMP.atom.Selection(last).get_selected_particles()[0]          
            pt1=IMP.atom.Selection(first).get_selected_particles()[0]
            #print IMP.core.XYZR(pt0).get_radius()+IMP.core.XYZR(pt1).get_radius()            
            r=IMP.core.PairRestraint(dps,IMP.ParticlePair(pt0,pt1))
            print pt0.get_name(),pt1.get_name()
            sortedsegments_cr.add_restraint(r)
            self.connected_intra_pairs.append((pt0,pt1))
            self.connected_intra_pairs.append((pt1,pt0))
            # only allow the particles to separate by 1 angstrom
            # self.m.set_maximum_score(r, self.kappa)

        self.m.add_restraint(sortedsegments_cr)
        self.m.add_restraint(unmodeledregions_cr)
        self.sortedsegments_cr_dict[name]=sortedsegments_cr
        self.unmodeledregions_cr_dict[name]=unmodeledregions_cr
        self.prot.add_child(protein_h)
        #self.rigid_bodies+=RigiParticles
        self.floppy_bodies+=FlexParticles

    def create_rotational_symmetry(self,maincopy,copies):
        #still working on it!
        from math import pi
        ncopies=len(copies)+1

        sel=IMP.atom.Selection(self.prot,molecule=maincopy)              
        mainparticles=sel.get_selected_particles()
        
        for k in range(len(copies)):
          rotation3D=IMP.algebra.get_rotation_about_axis(IMP.algebra.Vector3D(0,0,1), 2*pi/ncopies*(k+1))
          sm=IMP.core.TransformationSymmetry(rotation3D)

        
          sel=IMP.atom.Selection(self.prot,molecule=copies[k])              
          copyparticles=sel.get_selected_particles()

           
          mainpurged=[]
          copypurged=[]
          for n,p in enumerate(mainparticles):
           
            pc=copyparticles[n]
            
            #print p.get_name(),pc.get_name()


            #if (p.get_name()!=pc.get_name()): print p.get_name()+" "+pc.get_name()+" not the same particle"

            '''
            if IMP.core.RigidMember.particle_is_instance(p):
               rb=IMP.core.RigidMember(p).get_rigid_body()
               if not rb in mainpurged: 
                  mainpurged.append(rb)
                  IMP.pmi.Symmetric.setup_particle(rb,0)
                  print "added Rigid Body", type(rb),rb.get_name()
            elif IMP.core.NonRigidMember.particle_is_instance(p):   
                  mainpurged.append(p)
                  IMP.pmi.Symmetric.setup_particle(p,0)
                  print "added NonRigidMember", type(p),p.get_name()    
                  print "warning: NonRigidMember particles might not work with symmetries"
                  print "create rigid body with nonrigidmembers=False"                                        
            else:
               print "added Particle", type(p),p.get_name()
               mainpurged.append(p)
               IMP.pmi.Symmetric.setup_particle(p,0)
            '''
            #print "added Particle", type(p),p.get_name()
            mainpurged.append(p)
            IMP.pmi.Symmetric.setup_particle(p,0)           
            
            '''
            if IMP.core.RigidMember.particle_is_instance(pc):
               rbc=IMP.core.RigidMember(pc).get_rigid_body()
               if not rbc in copypurged: 
                  copypurged.append(rbc)
                  IMP.pmi.Symmetric.setup_particle(rbc,1)
                  print "added Rigid Body", type(rbc),rbc.get_name()
            elif IMP.core.NonRigidMember.particle_is_instance(pc):   
                  copypurged.append(pc)
                  print "added NonRigidMember", type(pc),pc.get_name()  
                  print "warning: NonRigidMember particles might not work with symmetries"
                  print "create rigid body with nonrigidmembers=False" 
                  IMP.pmi.Symmetric.setup_particle(pc,1)   
            else:
               print "added Particle", type(pc),pc.get_name()
               copypurged.append(pc)
               IMP.pmi.Symmetric.setup_particle(pc,1)
            '''

            #print "added Particle", type(pc),pc.get_name()
            copypurged.append(pc)
            IMP.pmi.Symmetric.setup_particle(pc,1)
           
                       
          lc=IMP.container.ListSingletonContainer(self.m)        
          for n,p in enumerate(mainpurged):
            
            pc=copypurged[n]
            print "setting "+p.get_name()+" as reference for "+pc.get_name()    
            
            
            IMP.core.Reference.setup_particle(pc,p)
            lc.add_particle(pc)
          
          c=IMP.container.SingletonsConstraint(sm,None,lc)
          self.m.add_score_state(c)        
        
        self.m.update()



    def link_components_to_rmf(self,rmfname,frameindex):
        '''
        load coordinates in the current representation
        this should be done only if the hierarchy self.prot is identical to the one
        i.e. all components were added
        as stored in the rmf
        '''
        rh= RMF.open_rmf_file(rmfname)
        imprmf.link_hierarchies(rh, [self.prot])
        imprmf.load_frame(rh, frameindex)
        del rh

    def create_components_from_rmf(self,rmfname,frameindex):
        '''
        still not working.
        create the representation (i.e. hierarchies) from the rmf file.
        it will be stored in self.prot, which will be overwritten.
        load the coordinates from the rmf file at frameindex.
        '''
        rh= RMF.open_rmf_file(rmfname)
        self.prot=imprmf.create_hierarchies(rh, self.m)[0]
        IMP.atom.show_molecular_hierarchy(self.prot)
        imprmf.link_hierarchies(rh, [self.prot])
        imprmf.load_frame(rh, frameindex)
        del rh
        for p in self.prot.get_children():
            self.add_component_name(p.get_name())
            self.hier_dict[p.get_name()]=p
        '''
        still missing: save rigid bodies contained in the rmf in self.rigid_bodies
        save floppy bodies in self.floppy_bodies
        get the connectivity restraints
        '''

    def set_rigid_body_from_hierarchies(self,hiers):
        rigid_parts=set()
        name=""
        print "set_rigid_body_from_hierarchies> setting up a new rigid body"
        for hier in hiers:
            ps=IMP.atom.get_leaves(hier)
            for p in ps:
              if IMP.core.RigidMember.particle_is_instance(p):
                 rb=IMP.core.RigidMember(p).get_rigid_body()
                 print "set_rigid_body_from_hierarchies> WARNING particle %s already belongs to rigid body %s" % (p.get_name(),rb.get_name())
              else:
                 rigid_parts.add(p)
            name+=hier.get_name()+"-"
            print "set_rigid_body_from_hierarchies> adding %s to the rigid body" % hier.get_name()
        rb=IMP.atom.create_rigid_body(list(rigid_parts))
        rb.set_coordinates_are_optimized(True)
        rb.set_name(name+"rigid_body")
        self.rigid_bodies.append(rb)

    def set_super_rigid_body_from_hierarchies(self,hiers):
        super_rigid_xyzs=set()
        super_rigid_rbs=set()
        name=""
        print "set_super_rigid_body_from_hierarchies> setting up a new SUPER rigid body"
        for hier in hiers:
            ps=IMP.atom.get_leaves(hier)
            for p in ps:
              if IMP.core.RigidMember.particle_is_instance(p):
                 rb=IMP.core.RigidMember(p).get_rigid_body()
                 super_rigid_rbs.add(rb)
              else:
                 super_rigid_xyzs.add(p)
            print "set_rigid_body_from_hierarchies> adding %s to the rigid body" % hier.get_name()
        self.super_rigid_bodies.append((super_rigid_xyzs,super_rigid_rbs))


    def set_rigid_bodies(self,subunits,coords=None,nonrigidmembers=True):
        if coords==None: coords=()
        #sometimes, we know about structure of an interaction
        #and here we make such PPIs rigid
        randomize_coords = lambda c: tuple(1.*(nrrand(3)-0.5)+array(c))
        
        rigid_parts=set()
        for s in subunits:
            if type(s)==type(tuple()) and len(s)==2:
               sel=IMP.atom.Selection(self.prot,molecule=s[0],residue_indexes=range(s[1][0],s[1][1]+1))
               if len(sel.get_selected_particles())==0: 
                  print "set_rigid_bodies: selected particle does not exists"
               for p in sel.get_selected_particles():
                  #if not p in self.floppy_bodies:
                      if IMP.core.RigidMember.particle_is_instance(p):
                         rb=IMP.core.RigidMember(p).get_rigid_body()
                         print "set_rigid_body_from_hierarchies> WARNING particle %s already belongs to rigid body %s" % (p.get_name(),rb.get_name())
                      else:
                         rigid_parts.add(p)
                     
               
            elif type(s)==type(str()):
               sel=IMP.atom.Selection(self.prot,molecule=s)
               if len(sel.get_selected_particles())==0: 
                  print "set_rigid_bodies: selected particle does not exists"
               for p in sel.get_selected_particles():
                  #if not p in self.floppy_bodies:
                      if IMP.core.RigidMember.particle_is_instance(p):
                         rb=IMP.core.RigidMember(p).get_rigid_body()
                         print "set_rigid_body_from_hierarchies> WARNING particle %s already belongs to rigid body %s" % (p.get_name(),rb.get_name())
                      else:
                         rigid_parts.add(p)
        
        rb=IMP.atom.create_rigid_body(list(rigid_parts))
        rb.set_coordinates_are_optimized(True)
        rb.set_name(''.join(str(subunits))+"_rigid_body")        
        if type(coords)==tuple and len(coords)==3: rb.set_coordinates(randomize_coords(coords))
        self.rigid_bodies.append(rb)
        
    def set_super_rigid_bodies(self,subunits,coords=None):
        super_rigid_xyzs=set()
        super_rigid_rbs=set()
        
        for s in subunits:
            if type(s)==type(tuple()) and len(s)==2:
               sel=IMP.atom.Selection(self.prot,molecule=s[0],residue_indexes=range(s[1][0],s[1][1]+1))
               if len(sel.get_selected_particles())==0: 
                  print "set_rigid_bodies: selected particle does not exists"
               for p in sel.get_selected_particles():
                      if IMP.core.RigidMember.particle_is_instance(p):
                         rb=IMP.core.RigidMember(p).get_rigid_body()
                         super_rigid_rbs.add(rb)
                      else:
                         super_rigid_xyzs.add(p)
            elif type(s)==type(str()):
               sel=IMP.atom.Selection(self.prot,molecule=s)
               if len(sel.get_selected_particles())==0: 
                  print "set_rigid_bodies: selected particle does not exists"
               for p in sel.get_selected_particles():
                  #if not p in self.floppy_bodies:
                      if IMP.core.RigidMember.particle_is_instance(p):
                         rb=IMP.core.RigidMember(p).get_rigid_body()
                         super_rigid_rbs.add(rb)
                      else:
                         super_rigid_xyzs.add(p)
        self.super_rigid_bodies.append((super_rigid_xyzs,super_rigid_rbs))       

    def set_floppy_bodies(self):
        for p in self.floppy_bodies:
            name=p.get_name()
            p.set_name(name+"_floppy_body")
            if IMP.core.RigidMember.particle_is_instance(p):
                print "I'm trying to make this particle flexible although it was assigned to a rigid body", p.get_name()
                rb=IMP.core.RigidMember(p).get_rigid_body()
                rb.set_is_rigid_member(p.get_particle_index(),False)
                p.set_name(p.get_name()+"_rigid_body_member")
    
    def get_particles_from_selection(self,selection_tuples):
        #to be used for instance by CompositeRestraint
        #selection tuples must be [(r1,r2,"name1"),(r1,r2,"name2"),....]
        particles=[]
        
        for s in selection_tuples:
            if type(s)==tuple and len(s)==3:
              sel=IMP.atom.Selection(self.prot,molecule=s[2],residue_indexes=range(s[0],s[1]+1))
            elif type(s)==str:
              sel=IMP.atom.Selection(self.prot,molecule=s)              
            ps=sel.get_selected_particles()
            print "get_particles_from_selection: "+str(s)+" selected "+str(len(ps))+" particles"
            particles+=ps
            
        return particles

    def get_connected_intra_pairs(self):
        return self.connected_intra_pairs

    def set_rigid_bodies_max_trans(self,maxtrans):
        self.maxtrans_rb=maxtrans

    def set_rigid_bodies_max_rot(self,maxrot):
        self.maxrot_rb=maxrot

    def set_floppy_bodies_max_trans(self,maxtrans):
        self.maxtrans_fb=maxtrans
        
    def set_rigid_bodies_as_fixed(self,rigidbodiesarefixed=True):
        '''
        this function will fix rigid bodies in their actual
        position. the get_particles_to_sample function will return
        just the floppy bodies.
        '''
        self.rigidbodiesarefixed=rigidbodiesarefixed
        

    def get_particles_to_sample(self):
        #get the list of samplable particles with their type
        #and the mover displacement. Everything wrapped in a dictionary,
        #to be used by samplers modules
        ps={}
        
        #remove symmetric particles: they are not sampled
        rbtmp=[]
        fbtmp=[]
        if not self.rigidbodiesarefixed:
            for rb in self.rigid_bodies:
               if IMP.pmi.Symmetric.particle_is_instance(rb):
                  if IMP.pmi.Symmetric(rb).get_symmetric()!=1:
                     rbtmp.append(rb)
               else: 
                  rbtmp.append(rb)
    
        for fb in self.floppy_bodies:
           if IMP.pmi.Symmetric.particle_is_instance(fb):
              if IMP.pmi.Symmetric(fb).get_symmetric()!=1:
                 fbtmp.append(fb)
           else: 
              fbtmp.append(fb)   
              
        self.rigid_bodies=rbtmp
        self.floppy_bodies=fbtmp

        ps["Rigid_Bodies_SimplifiedModel"]=(self.rigid_bodies,self.maxtrans_rb,self.maxrot_rb)
        ps["Floppy_Bodies_SimplifiedModel"]=(self.floppy_bodies,self.maxtrans_fb)
        ps["SR_Bodies_SimplifiedModel"]=(self.super_rigid_bodies,self.maxtrans_srb,self.maxrot_srb)
        print ps
        return ps
    
    def set_output_level(self,level):
        self.output_level=level
    
    def get_output(self):
        output={}
        score=0.0
        
        output["SimplifiedModel_Total_Score_"+self.label]=str(self.m.evaluate(False))  
        output["SimplifiedModel_Linker_Score_"+self.label]=str(self.linker_restraints.unprotected_evaluate(None))
        for name in self.sortedsegments_cr_dict:
            partialscore=self.sortedsegments_cr_dict[name].evaluate(False)
            score+=partialscore
            output["SimplifiedModel_Link_SortedSegments_"+name+"_"+self.label]=str(partialscore)
            partialscore=self.unmodeledregions_cr_dict[name].evaluate(False)
            score+=partialscore            
            output["SimplifiedModel_Link_UnmodeledRegions_"+name+"_"+self.label]=str(partialscore)
        for name in self.linker_restraints_dict:
            output[name+"_"+self.label]=str(self.linker_restraints_dict[name].unprotected_evaluate(None))
        if self.output_level=="high":
            #print coordinates
            for p in IMP.atom.get_leaves(self.prot):
                d=IMP.core.XYZR(p)
                output["Coordinates_"+p.get_name()+"_"+self.label]=str(d)
                
        output["_TotalScore"]=str(score)
        return output

    def get_hierarchy(self):
        return  self.prot

    def draw_hierarchy_graph(self):
        for c in IMP.atom.Hierarchy(self.prot).get_children():
            print "Drawing hierarchy graph for "+c.get_name()
            output.get_graph_from_hierarchy(c)
            

    def get_geometries(self):
        #create segments at the lowest resolution
        seggeos=[]
        for name in self.hier_geometry_pairs:
            for pt in self.hier_geometry_pairs[name]:
                p1=pt[0]
                p2=pt[1]
                color=pt[2]
                try:
                    clr=IMP.display.get_rgb_color(color)
                except:
                    clr=IMP.display.get_rgb_color(1.0)
                coor1=IMP.core.XYZ(p1).get_coordinates()
                coor2=IMP.core.XYZ(p2).get_coordinates()
                seg=IMP.algebra.Segment3D(coor1,coor2)
                seggeos.append(IMP.display.SegmentGeometry(seg,clr))
        return seggeos

    def setup_bonds(self):
        #create segments at the lowest resolution
        seggeos=[]
        for name in self.hier_geometry_pairs:
            for pt in self.hier_geometry_pairs[name]:
                p1=pt[0]
                p2=pt[1]
                IMP.atom.create_bond(IMP.atom.Bonded.setup_particle(p1),IMP.atom.Bonded.setup_particle(p2),1)

    def show_component_table(self,name):
        residues=set()
        prot=self.hier_dict[name]

        for p in IMP.atom.get_leaves(prot):
            residues.update(IMP.atom.Fragment(p).get_residue_indexes())
        
        if name in self.sequence_dict:
           lastresn_fromsequence=len(self.sequence_dict[name])
           residues.update([1,lastresn_fromsequence])
        
        firstresn=min(residues)
        lastresn=max(residues)
        
        for nres in range(firstresn,lastresn+1):
            s=IMP.atom.Selection(prot,residue_index=nres)
            ps=s.get_selected_particles()
            if len(ps)==0:  
               print "%20s %20s" % (name,nres), "**** not represented ****"
            else:
               print "%20s %7s" % (name,nres), " ".join(["%20s %7s" % (str(p.get_name()),
                     str(IMP.pmi.Resolution(p).get_resolution())) for p in ps])

                
    def draw_hierarchy_composition(self):

        ks=self.elements.keys()
        ks.sort()
 
        max=0
        for k in self.elements:
            for l in self.elements[k]:
                if l[1]>max: max=l[1]
        
        for k in ks:
            self.draw_component_composition(k,max)
        
    def draw_component_composition(self,name,max=1000):
            from matplotlib import pyplot
            import matplotlib as mpl    
            k=name
            list=sorted(self.elements[k], key=itemgetter(0))
            endres=list[-1][1]
            fig = pyplot.figure(figsize=(26.0*float(endres)/max+2,2))
            ax = fig.add_axes([0.05, 0.475, 0.9, 0.15])
            
            # Set the colormap and norm to correspond to the data for which
            # the colorbar will be used.
            cmap = mpl.cm.cool
            norm = mpl.colors.Normalize(vmin=5, vmax=10)
            bounds=[1]
            colors=[]
            
            print k
            
            for n,l in enumerate(list):
                firstres=l[0]
                lastres=l[1]
                if l[3]!="end":
                    if bounds[-1]!=l[0]:
                       colors.append("white")
                       bounds.append(l[0])
                       if l[3]=="pdb": colors.append("#99CCFF")
                       if l[3]=="bead": colors.append("#FFFF99")
                       if l[3]=="helix": colors.append("#33CCCC")
                       if l[3]!="end":
                          bounds.append(l[1]+1)
                    else:
                       if l[3]=="pdb": colors.append("#99CCFF")
                       if l[3]=="bead": colors.append("#FFFF99")
                       if l[3]=="helix": colors.append("#33CCCC")
                       if l[3]!="end":                   
                          bounds.append(l[1]+1)
                else:
                   if bounds[-1]-1==l[0]:
                      bounds.pop()
                      bounds.append(l[0])
                   else:
                      colors.append("white")
                      bounds.append(l[0])                     
                            
            bounds.append(bounds[-1])
            colors.append("white")
            
            '''
            for n,l in enumerate(list):
                firstres=l[0]
                lastres=l[1]
  
                if l[3]=="end" and bounds[-1]!=list[n-1][1]: 
                   colors.append("white")
                   bounds.append(list[n-1][1])
                elif bounds[-1]!=l[0]:
                   print bounds[-1],l[0]          
                   bounds.append(l[0])
                   if l[3]=="pdb": colors.append("#99CCFF")
                   if l[3]=="bead": colors.append("#FFFF99")
            
            if bounds[-1]!=endres:
               bounds.append(endres)
            '''
            cmap = mpl.colors.ListedColormap(colors)
            cmap.set_over('0.25')
            cmap.set_under('0.75')
            
            norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
            cb2 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                     norm=norm,
                                     # to use 'extend', you must
                                     # specify two extra boundaries:
                                     boundaries=bounds,
                                     ticks=bounds, # optional
                                     spacing='proportional',
                                     orientation='horizontal')
            
            extra_artists=[]
            npdb=0
            for l in list:  
                if l[3]=="pdb": 
                   npdb+=1                   
                   mid=1.0/endres*float(l[0])
                   #t =ax.text(mid, float(npdb-1)/2.0+1.5, l[2], ha="left", va="center", rotation=0,
                   #size=10)
                   #t=ax.annotate(l[0],2)
                   t=ax.annotate(l[2], xy=(mid, 1),  xycoords='axes fraction',
                   xytext=(mid+0.025, float(npdb-1)/2.0+1.5), textcoords='axes fraction',
                   arrowprops=dict(arrowstyle="->",
                                connectionstyle="angle,angleA=0,angleB=90,rad=10"),
                   )
                   extra_artists.append(t)
            
            #set the title of the bar
            title=ax.text(-0.005, 0.5, k, ha="right", va="center", rotation=90,
                size=15)
            
            extra_artists.append(title)
            #changing the xticks labels
            labels=len(bounds)*[" "]
            ax.set_xticklabels(labels)
            mid=1.0/endres*float(bounds[0])
            t=ax.annotate(bounds[0], xy=(mid, 0),  xycoords='axes fraction',
                   xytext=(mid-0.01, -0.5), textcoords='axes fraction',)            
            extra_artists.append(t)
            offsets=[0,-0.5,-1.0]
            nclashes=0
            for n in range(1,len(bounds)):
                if bounds[n]==bounds[n-1]: continue
                mid=1.0/endres*float(bounds[n])
                if (float(bounds[n])-float(bounds[n-1]))/max<=0.01:
                   nclashes+=1
                   offset=offsets[nclashes%3]
                else:
                   nclashes=0
                   offset=offsets[0]
                if offset>-0.75: 
                  t=ax.annotate(bounds[n], xy=(mid, 0),  xycoords='axes fraction',
                     xytext=(mid, -0.5+offset), textcoords='axes fraction')
                else:
                  t=ax.annotate(bounds[n], xy=(mid, 0),  xycoords='axes fraction',
                     xytext=(mid, -0.5+offset), textcoords='axes fraction',arrowprops=dict(arrowstyle="-"))
                extra_artists.append(t)
                            
            cb2.add_lines(bounds,["black"]*len(bounds),[1]*len(bounds))
            #cb2.set_label(k)   
            
            pyplot.savefig(k+"structure.png",dpi=150,transparent="True",bbox_extra_artists=(extra_artists), bbox_inches='tight')
            pyplot.show()

    def get_prot_name_from_particle(self,particle):
            names=self.get_component_names()
            particle0=particle
            name=None
            while not name in names:
                 h=IMP.atom.Hierarchy(particle0).get_parent()
                 name=h.get_name()
                 particle0=h.get_particle()
            return name

    def get_random_residue_pairs(self,names,resolution,number):
            from random import choice
            particles=[]
            for name in names:
                prot=self.hier_dict[name]
                particles+=tools.get_particles_by_resolution(prot,resolution)
            
            random_residue_pairs=[]
            for i in range(number):
                p1=choice(particles)
                p2=choice(particles)
                r1=choice(IMP.atom.Fragment(p1).get_residue_indexes())
                r2=choice(IMP.atom.Fragment(p2).get_residue_indexes())            
                name1=self.get_prot_name_from_particle(p1)
                name2=self.get_prot_name_from_particle(p2)
                random_residue_pairs.append((name1,r1,name2,r2))
            return random_residue_pairs           