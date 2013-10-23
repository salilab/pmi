#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.display
import IMP.pmi

class Rods():
    def __init__(self,m):
        self.m=m
        self.hier=IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.m))
        self.rigid_bodies=[]
        self.floppy_bodies=[]        
        self.maxtrans_rb=2.0
        self.maxtrans_fb=2.0  
        self.maxrot_rb=0.15
    
    def add_protein(self,name,(firstres,lastres)):
        from math import pi,cos,sin
        h = IMP.atom.Molecule.setup_particle(IMP.Particle(self.m))
        h.set_name(name)
        nres=lastres-firstres
        radius=(nres)*5/2/pi
        
        for res in range(firstres,lastres):
            alpha=2*pi/nres*(res-firstres)
            x=radius*cos(alpha)
            y=radius*sin(alpha)
            p=IMP.Particle(self.m)
            r=IMP.atom.Residue.setup_particle(p,IMP.atom.ALA,res)
            d=IMP.core.XYZR.setup_particle(p,5.0)
            d.set_coordinates(IMP.algebra.Vector3D((x,y,0)))   
            d.set_coordinates_are_optimized(True)                     
            h.add_child(r)
        self.hier.add_child(h)
    
    def get_hierarchy(self):
        return self.hier
    
    def set_rod(self,chainname,(firstres,lastres)):
        prb=IMP.Particle(self.m)
        sel=IMP.atom.Selection(self.hier,molecule=chainname,residue_indexes=range(firstres,lastres+1))
        ps=sel.get_selected_particles()
        rb=IMP.core.RigidBody.setup_particle(prb,ps)
        self.rigid_bodies.append(rb)
    
    def get_particles_to_sample(self):
        #get the list of samplable particles with their type
        #and the mover displacement. Everything wrapped in a dictionary,
        #to be used by samplers modules
        ps={}
        ps["Floppy_Bodies_Rods"]=(self.floppy_bodies,self.maxtrans_fb)
        ps["Rigid_Bodies_Rods"]=(self.rigid_bodies,self.maxtrans_rb,self.maxrot_rb)
        return ps                


class Beads():
    def __init__(self,m):
        self.m=m
        self.beads=[]
        self.nresidues=0
        self.hier=IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.m))
        self.floppy_bodies=[]
        self.maxtrans_fb=0.2
        self.particle_database={}

    def add_bead(self,radius,label="None",color=None):

        p=IMP.Particle(self.m)
        p.set_name(label)
        self.particle_database[label]=p
        self.floppy_bodies.append(p)
        #set default coordinates 0,0,0
        d=IMP.core.XYZ.setup_particle(p)
        d=IMP.core.XYZR.setup_particle(p,radius)
        d.set_coordinates_are_optimized(True)
        #a=IMP.atom.Atom.setup_particle(p,IMP.atom.AT_CA)
        #p=IMP.Particle(self.m)
        self.nresidues+=1
        #r=IMP.atom.Residue.setup_particle(p,IMP.atom.ALA,self.nresidues)
        #r.add_child(a)
        #self.hier.add_child(r)
        self.hier.add_child(p)
        if color!=None: self.set_color(label,color)
        return self.particle_database[label]
    
    def set_color(self,label,value):
        p=self.particle_database[label]
        clr=IMP.display.get_rgb_color(value)
        IMP.display.Colored.setup_particle(p,clr) 
    
    def set_floppy_bodies_max_trans(self,maxtrans):
        self.maxtrans_fb=maxtrans

    def get_hierarchy(self):
        return self.hier

    def get_bead(self,label):
        return self.particle_database[label]

    def set_maxtrans_fb(self,maxtrans_fb):
        self.maxtrans_fb=maxtrans_fb

    def get_particles_to_sample(self):
        #get the list of samplable particles with their type
        #and the mover displacement. Everything wrapped in a dictionary,
        #to be used by samplers modules
        ps={}
        ps["Floppy_Bodies_Beads"]=(self.floppy_bodies,self.maxtrans_fb)
        return ps

class MultipleStates():
    def __init__(self,nstates,m):
        global itertools, tools, restraints
        
        import itertools
        import IMP.pmi.tools as tools
        import IMP.pmi.restraints as restraints

        self.floppy_bodies=[]
        self.rigid_bodies=[]
        for ncopy in range(nstates):
            self.floppy_bodies.append([])
            self.rigid_bodies.append([])
        
        self.rigid_bodies_are_sampled=True
        self.floppy_bodies_are_sampled=True
        self.prot=[]
        self.refprot=[]
        self.prot_lowres={}
        self.nstates=nstates
        self.label="None"

        #model decorator list
        self.xyzmodellist=[]
        self.xyzreflist=[]
        self.maxtrans_rb=0.15
        self.maxrot_rb  =0.03
        self.maxtrans_fb=0.15
        self.m = m

    def get_model(self):
        return self.m

    def set_label(self,label):
        self.label=label
    
    def set_rigid_bodies_are_sampled(self,input=True):
        self.rigid_bodies_are_sampled=input

    def set_floppy_bodies_are_sampled(self,input=True):
        self.floppy_bodies_are_sampled=input
            
    def get_rigid_bodies(self):
        return  self.rigid_bodies

    def set_rigid_bodies_max_trans(self,maxtrans):
        self.maxtrans_rb=maxtrans

    def set_rigid_bodies_max_rot(self,maxrot):
        self.maxrot_rb=maxrot

    def set_floppy_bodies_max_trans(self,maxtrans):
        self.maxtrans_fb=maxtrans

    def get_hierarchies(self):
        return  self.prot

    def destroy_residues(self,segments):
        #segments are defined as a list of tuples ex [(res1,res2,chain),....]
        #this function must be called before the rigid body definition!
        for prot in self.prot:
            for segment in segments:
                #rinterval=[(segment[0],segment[1]+1)]
                if (segment[0]==-1 or segment[1]==-1):
                    s=IMP.atom.Selection(prot,chains=segment[2])
                else:
                    s=IMP.atom.Selection(prot,chains=segment[2],residue_indexes=range(segment[0],segment[1]+1))
                for p in s.get_selected_particles():
                    if IMP.core.RigidMember.particle_is_instance(p):
                        print "MultipleStates: one particle was not destroied because it was a RigidMember."
                    else:
                        #destroy the residue and the associated atom
                        a=IMP.atom.Atom(p)
                        r=IMP.atom.Residue(IMP.atom.Atom(p).get_parent())
                        IMP.atom.destroy(r)
                        #IMP.atom.destroy(a)
                        #IMP.atom.destroy(p)
            IMP.atom.show_molecular_hierarchy(prot)

    def add_residues_to_chains(self,residuechainlist,residue_type=IMP.atom.LYS):
        #add a list of residues to the corresponding list
        #for instance residuechainlist=[(35,"A"),(100,"B")] will add
        #residue 35 to chain A and residue 100 to chain B
        for rc in residuechainlist:

            s=IMP.atom.Selection(self.prot[0],chains=rc[1],residue_index=rc[0],atom_type=IMP.atom.AT_CA)

            print s.get_selected_particles()
            if len(s.get_selected_particles())==0:
                for prot in self.prot:
                    print "adding " + str(rc)
                    p=IMP.Particle(self.m)
                    #set default coordinates 0,0,0
                    d=IMP.core.XYZ.setup_particle(p)
                    IMP.core.XYZR.setup_particle(p,0.0)
                    d.set_coordinates_are_optimized(True)
                    a=IMP.atom.Atom.setup_particle(p,IMP.atom.AT_CA)
                    p=IMP.Particle(self.m)
                    r=IMP.atom.Residue.setup_particle(p,residue_type,rc[0])
                    r.add_child(a)
                    p=IMP.Particle(self.m)
                    c=IMP.atom.Chain.setup_particle(p,rc[1])
                    c.add_child(r)
                    prot.add_child(c)
                    print tools.get_residue_index_and_chain_from_particle(a)


            else:
                p=s.get_selected_particles()[0]
                print rc, s.get_selected_particles()[0] #, tools.get_residue_index_and_chain_from_particle(s.get_selected_particles()[0])



            #test that that was indeed added:

            s=IMP.atom.Selection(self.prot[0],chains=rc[1],residue_index=rc[0],atom_type=IMP.atom.AT_CA)

            print s.get_selected_particles()

    def add_beads(self,segments,xyzs=None,radii=None,colors=None):
        '''
        this method generate beads in missing portions.
        The segments argument must be a list of selections 
        in the form [(firstres,lastres,chain)]
        each selection will generate a bead
        '''
        if xyzs==None: xyzs=[]
        if radii==None: radii=[]
        if colors==None: colors=[]
        
        from math import pi
        
        for n,s in enumerate(segments):
          firstres=s[0]
          lastres=s[1]
          chainid=s[2]
          nres=s[1]-s[0]
          for prot in self.prot:         
            for prot in self.prot: 
                cps=IMP.atom.get_by_type(prot, IMP.atom.CHAIN_TYPE)
                for c in cps:
                   chain=IMP.atom.Chain(c)
                   if chain.get_id()==chainid:
                       p=IMP.Particle(self.m)
                       f=IMP.atom.Fragment.setup_particle(p)
                       rindexes=range(firstres,lasteres+1)
                       f.set_residue_indexes(rindexes)
                       f.set_name("Fragment_"+'%i-%i' % (firstres,lastres))
                       chain.add_child(f)
                       mass=len(rindexes)*110.0
                       vol=IMP.atom.get_volume_from_mass(mass)
                       if n+1>len(radii):
                          mass=len(rindexes)*110.0
                          vol=IMP.atom.get_volume_from_mass(mass)
                          radius=(3*vol/math.pi)**(1/3)
                       else:
                          radius=radii[n]

                       if n+1>len(xyzs):
                          x=0
                          y=0
                          z=0
                       else:
                          x=xyzs[n][0]
                          y=xyzs[n][1]         
                          z=xyzs[n][2]
                       
                       if n+1<=len(colors):
                          clr=IMP.display.get_rgb_color(colors[n])
                          IMP.display.Colored.setup_particle(prt,clr)
                                                                  
                       d=IMP.atom.XYZR.setup_particle(p,IMP.algebra.Sphere3D(x,y,z,radius))
                       
    
    def renumber_residues(self,chainid,newfirstresiduenumber):
            for prot in self.prot: 
                cps=IMP.atom.get_by_type(prot, IMP.atom.CHAIN_TYPE)
                for c in cps:
                    if IMP.atom.Chain(c).get_id()==chainid:
                       ps=c.get_children()
                       r=IMP.atom.Residue(ps[0])
                       ri=r.get_index()
                       offs=newfirstresiduenumber-ri
                       for p in ps:
                           r=IMP.atom.Residue(p)
                           ri=r.get_index()                            
                           r.set_index(ri+offs)
                
    def destroy_everything_but_the_residues(self,segments):
        #segments are defined as a list of tuples ex [(res1,res2,chain),....]
        for prot in self.prot:
            pstokeep=[]
            for segment in segments:

                #rinterval=[(segment[0],segment[1]+1)]
                if (segment[0]==-1 or segment[1]==-1):
                    s=IMP.atom.Selection(prot,chains=segment[2])
                else:
                    s=IMP.atom.Selection(prot,chains=segment[2],residue_indexes=range(segment[0],segment[1]+1))
                pstokeep+=s.get_selected_particles()

            for p in IMP.atom.get_leaves(prot):
                if p not in pstokeep:
                    if IMP.core.RigidMember.particle_is_instance(p):
                        print "MultipleStates: one particle was not destroied because it was a RigidMember."
                    else:
                        #destroy the residue and the associated atom
                        a=IMP.atom.Atom(p).get_parent()
                        r=IMP.atom.Residue(IMP.atom.Atom(p).get_parent())
                        #IMP.atom.destroy(a)
                        IMP.atom.destroy(r)
                        #self.m.remove_particle(p)

    def generate_linkers_restraint_and_floppy_bodies(self,segment):
        '''
        this methods automatically links the particles consecutively
        according to the sequence. The restraint applied is a harmonic upper bound,
        with a distance that is proportional to the number of residues
        in the gap.
        '''
        #this function will create floppy bodies where there are not
        #rigid bodies and moreover create a linker restraint between them
        linker_restraint_objects=[]
        for ncopy,prot in enumerate(self.prot):
            if (segment[0]==-1 or segment[1]==-1):
                s=IMP.atom.Selection(prot,chains=segment[2])
            else:
                s=IMP.atom.Selection(prot,chains=segment[2],residue_indexes=range(segment[0],segment[1]+1))
            residue_indexes=[]
            for p in s.get_selected_particles():

                (r,c)=tools.get_residue_index_and_chain_from_particle(p)

                if IMP.core.RigidMember.particle_is_instance(p):
                    Floppy=False
                else:
                    (r,c)=tools.get_residue_index_and_chain_from_particle(p)
                    p.set_name(str(r)+":"+str(c))
                    tools.set_floppy_body(p)
                    self.floppy_bodies[ncopy].append(p)
                    Floppy=True
                residue_indexes.append((r,Floppy,c,p))

            residue_indexes.sort()


            pruned_residue_list=[]
            r0=residue_indexes[0]
            pruned_residue_list.append(r0)

            #generate the list of residues that define the intervals
            #between rigid bodies and floppy bodies
            for i in range(1,len(residue_indexes)):
                r=residue_indexes[i]
                if r[1]==r0[1] and r[1]==False and IMP.core.RigidMember(r[3]).get_rigid_body() == IMP.core.RigidMember(r0[3]).get_rigid_body():
                    r0=r
                elif r[1]==r0[1] and r[1]==False and IMP.core.RigidMember(r[3]).get_rigid_body() != IMP.core.RigidMember(r0[3]).get_rigid_body():
                    pruned_residue_list.append(r0)
                    pruned_residue_list.append(r)
                    r0=r
                elif r[1]!=r0[1] and r0[1]==False:
                    pruned_residue_list.append(r0)
                    pruned_residue_list.append(r)
                    r0=r
                elif r[1]==r0[1] and r0[1]==True:
                    pruned_residue_list.append(r)
                    r0=r
                elif r[1]!=r0[1] and r[1]==False:
                    pruned_residue_list.append(r)
                    r0=r

            
            r0=pruned_residue_list[0]
            linkdomaindef=[]
            
            for i in range(1,len(pruned_residue_list)):
                r=pruned_residue_list[i]
                if r[1]==r0[1] and r[1]==False and IMP.core.RigidMember(r[3]).get_rigid_body() == IMP.core.RigidMember(r0[3]).get_rigid_body():
                    r0=r
                else:
                    linkdomaindef.append((r0[0],r[0],r[2]))
                    r0=r
            
            print " creating linker between atoms defined by: "+str(linkdomaindef)
            
            ld=restraints.LinkDomains(prot,linkdomaindef,1.0,3.0)
            ld.set_label(str(ncopy))
            ld.add_to_model()
            linker_restraint_objects.append(ld)
            prs=ld.get_pairs()
        
        return linker_restraint_objects





    def get_ref_hierarchies(self):
        return  self.refprot

    def get_number_of_states(self):
        return  self.nstates

    def get_rigid_bodies(self):
        rblist=[]
        for rbl in self.rigid_bodies:
            for rb in rbl:
                rblist.append(rb)
        return rblist

    def get_floppy_bodies(self):
        fblist=[]
        for fbl in self.floppy_bodies:
            for fb in fbl:
                fblist.append(fb)
        return fblist

    def set_rigid_bodies(self,rigid_body_list):
        if len(self.prot)==0:
            print "MultipleStates.set_rigid_bodies: hierarchy was not initialized"
            exit()
        for ncopy,prot in enumerate(self.prot):
            rbl=[]
            for element in rigid_body_list:
                atoms=[]
                for interval in element:
                #rinterval upper bound is incremented by one because the
                #residue_indexes attribute cuts the upper edge
                    #rinterval=[(interval[0],interval[1]+1)]
                    if (interval[0]==-1 or interval[1]==-1):
                        s=IMP.atom.Selection(prot,chains=interval[2])
                    else:
                        s=IMP.atom.Selection(prot,chains=interval[2],residue_indexes=range(interval[0],interval[1]+1))
                    for p in s.get_selected_particles():
                        atoms.append(IMP.core.XYZR(p))

                    #add low resolution representation to the rigid bodies
                    for key in self.prot_lowres:
                        if (interval[0]==-1 or interval[1]==-1):
                            s=IMP.atom.Selection(self.prot_lowres[key][ncopy],chains=interval[2])
                        else:
                            s=IMP.atom.Selection(self.prot_lowres[key][ncopy],chains=interval[2],
                                                residue_indexes=range(interval[0],interval[1]+1))
                        for p in s.get_selected_particles():
                            atoms.append(IMP.core.XYZR(p))


                if len(atoms)>0:
                    prb=IMP.Particle(self.m)
                    rb=IMP.core.RigidBody.setup_particle(prb,atoms)
                    rb.set_name(str(element))
                    rbl.append(rb)
                else:
                    print "MultipleStates.set_rigid_bodies: selection " + str(interval) + "  has zero elements"
            self.rigid_bodies[ncopy]+=rbl

    def set_floppy_bodies(self,floppy_body_list):
        #define flexible regions within rigid bodies

        if len(self.prot)==0:
            print "MultipleStates: hierarchy was not initialized"
            exit()

        for ncopy,prot in enumerate(self.prot):
            atoms=[]
            for element in floppy_body_list:

                for interval in element:
                #rinterval upper bound is incremented by one because the
                #residue_indexes attribute cuts the upper edge
                    #rinterval=[(interval[0],interval[1]+1)]
                    if (interval[0]==-1 or interval[1]==-1):
                        s=IMP.atom.Selection(prot,chains=interval[2])
                    else:
                        s=IMP.atom.Selection(prot,chains=interval[2],residue_indexes=range(interval[0],interval[1]+1))
                    for p in s.get_selected_particles():
                        (r,c)=tools.get_residue_index_and_chain_from_particle(p)
                        tools.set_floppy_body(p)
                        p.set_name(str(r)+":"+str(c))
                        atoms.append(IMP.core.XYZR(p))
            self.floppy_bodies[ncopy]+=atoms

    def get_particles_to_sample(self):
        #get the list of samplable particles with their type
        #and the mover displacement. Everything wrapped in a dictionary,
        #to be used by samplers modules
        ps={}
        rblist=self.get_rigid_bodies()
        fblist=self.get_floppy_bodies()
        if self.rigid_bodies_are_sampled:
           ps["Rigid_Bodies_MultipleStates"]=(rblist,self.maxtrans_rb,self.maxrot_rb)
        if self.floppy_bodies_are_sampled:
           ps["Floppy_Bodies_MultipleStates"]=(fblist,self.maxtrans_fb)
        return ps

    def set_hierarchy_from_pdb(self,pdblistoflist):
        "the input is a list of list of pdbs"
        "one list for each copy"
        #eg [["pdb1_copy0","pdb2_copy0"],["pdb1_copy1","pdb2_copy1"]]"
        for copy in range(0,self.nstates):
            prot=self.read_pdbs(pdblistoflist[copy])
            self.prot.append(prot)
            xyz=IMP.core.XYZs(IMP.atom.get_leaves(prot))
            self.xyzmodellist.append(xyz)

    def set_ref_hierarchy_from_pdb(self,pdblistoflist):
        "the input is a list of list of pdbs"
        "one list for each copy"
        #eg [["pdb1_copy0","pdb2_copy0"],["pdb1_copy1","pdb2_copy1"]]"
        for copy in range(0,self.nstates):
            prot=self.read_pdbs(pdblistoflist[copy])
            self.refprot.append(prot)
            xyz=IMP.core.XYZs(IMP.atom.get_leaves(prot))
            self.xyzreflist.append(xyz)


    def read_pdbs(self,list_pdb_file):
        """read pdbs from an external list file
        create a simplified representation
        if the pdbs are given a individual strings, it will read the
        pdbs and give the chain name as specified in the pdb
        If it s a tuple like (filename,chainname) it will read
        the pdb and assing a name chainname
        to the chain"""

        hier=IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.m)) #create an empty hierarchy
        for pdb in list_pdb_file:
            if type(pdb)==str:
                h=IMP.atom.read_pdb(pdb, self.m,IMP.atom.AndPDBSelector(IMP.atom.CAlphaPDBSelector(), 
                                                                      IMP.atom.ATOMPDBSelector()))

                '''
                #destroy CA atoms, for the future
                for p in IMP.atom.get_leaves(h):
                    coor=IMP.core.XYZ(p).get_coordinates()
                    r=IMP.atom.Hierarchy(p).get_parent()
                    IMP.core.XYZ.setup_particle(r,coor)
                    IMP.atom.destroy(p)
                '''

                cps=IMP.atom.get_by_type(h, IMP.atom.CHAIN_TYPE)

                '''
                #consolidate the chains
                for c in cps:
                    cid=c.get_id()
                    s0=IMP.atom.Selection(hier, chains=cid)
                    try:
                      p=s0.get_selected_particles()[0]
                      re=IMP.atom.Residue(IMP.atom.Atom(p).get_parent()
                      ch=IMP.atom.Chain(re).get_parent())

                    except:
                      continue
                '''

                hier.add_child(h) #add read chains into hierarchy
            if type(pdb)==tuple:
                h=IMP.atom.read_pdb(pdb[0], self.m,IMP.atom.AndPDBSelector(IMP.atom.CAlphaPDBSelector(), 
                                                                      IMP.atom.ATOMPDBSelector()))

                '''
                #destroy CA atoms, for the future
                for p in IMP.atom.get_leaves(h):
                    coor=IMP.core.XYZ(p).get_coordinates()
                    r=IMP.atom.Hierarchy(p).get_parent()
                    IMP.core.XYZ.setup_particle(r,coor)
                    IMP.atom.destroy(p)
                '''

                cps=IMP.atom.get_by_type(h, IMP.atom.CHAIN_TYPE)
                for cp in cps:
                    IMP.atom.Chain(cp).set_id(pdb[1])
                hier.add_child(h) #add read chains into hierarchy

        return hier

    def recenter(self,prot):
        "recenter the hierarchy"
        ps=IMP.atom.get_leaves(prot)
        center = IMP.algebra.get_zero_vector_3d()
        for l in ps:
            center += IMP.core.XYZ(l).get_coordinates()
        center /= len(ps)
        for l in ps:
            d = IMP.core.XYZ(l)
            d.set_coordinates(d.get_coordinates() - center)
            d.set_coordinates_are_optimized(True)

        '''
        # bug generating code: keeping it for history

        rb=IMP.atom.create_rigid_body(prot)
        rbcoord=rb.get_coordinates()
        rot=IMP.algebra.get_identity_rotation_3d()
        tmptrans=IMP.algebra.Transformation3D(rot,rbcoord)
        trans=tmptrans.get_inverse()
        IMP.core.transform(rb,trans)
        IMP.core.RigidBody.teardown_particle(rb)
        self.m.remove_particle(rb)
        '''


    def shuffle_configuration(self,bounding_box_length):
        "shuffle configuration, used to restart the optimization"
        "it only works if rigid bodies were initialized"
        if len(self.rigid_bodies)==0:
            print "MultipleStates: rigid bodies were not intialized"
        hbbl=bounding_box_length/2
        for rbl in self.rigid_bodies:
            for rb in rbl:
                ub = IMP.algebra.Vector3D(-hbbl,-hbbl,-hbbl)
                lb = IMP.algebra.Vector3D( hbbl, hbbl, hbbl)
                bb = IMP.algebra.BoundingBox3D(ub, lb)
                translation = IMP.algebra.get_random_vector_in(bb)
                rotation = IMP.algebra.get_random_rotation_3d()
                transformation = IMP.algebra.Transformation3D(rotation, translation)
                rb.set_reference_frame(IMP.algebra.ReferenceFrame3D(transformation))

    def generate_simplified_hierarchy(self,nres):
        #generate a new multistate hierarchy
        self.prot_lowres[nres]=[]
        for prot in self.prot:
            sh=IMP.atom.create_simplified_along_backbone(prot, nres, False)
            print IMP.atom.get_leaves(sh)
            #for p in IMP.atom.get_leaves(sh):
            #    IMP.atom.Atom.setup_particle(p,IMP.atom.AT_CA)
            #s=IMP.atom.Selection(sh, chains="A",
            #              residue_index=958)
            #print s.get_selected_particles()[0]
            self.prot_lowres[nres].append(sh)


    def get_simplified_hierarchy(self,nres):
        return self.prot_lowres[nres]

    def calculate_drms(self):
        # calculate DRMSD matrix

        if len(self.xyzmodellist)==0:
            print "MultipleStates: hierarchies were not intialized"
            
        if len(self.xyzreflist)==0:
            print "MultipleStates: reference hierarchies were not intialized"

        drmsd={}
        for i in range(len(self.xyzreflist)):
            for j in range(len(self.xyzmodellist)):
                try:
                    drmsdval=IMP.atom.get_drmsd(self.xyzmodellist[j],self.xyzreflist[i])
                except:
                    drmsdval=tools.get_drmsd(self.xyzmodellist[j],self.xyzreflist[i])
                drmsd["MultipleStates_DRMSD_"+str(i)+"-Model_"+str(j)]=drmsdval

        # calculate model-template assignment that gives minimum total drmsd
        min_drmsd=[]
        for assign in itertools.permutations(range(len(self.xyzreflist))):
            s=0.
            for i,j in enumerate(assign):
                s+=drmsd["MultipleStates_DRMSD_"+str(j)+"-Model_"+str(i)]
            min_drmsd.append(s)

        drmsd["MultipleStates_Total_DRMSD"]=min(min_drmsd)
        return drmsd

    def get_output(self):
        output={}
        if len(self.refprot)!=0:
            drms=self.calculate_drms()
            output.update(drms)          
        output["MultipleStates_Total_Score_"+self.label]=str(self.m.evaluate(False))
        return output





class SimplifiedModel():
#Peter Cimermancic and Riccardo Pellarin
    '''
    This class creates the molecular hierarchies for the various involved proteins.    
    '''

    def __init__(self,m,upperharmonic=True,disorderedlength=False):
        global random, itemgetter,tools,nrrand,array,imprmf,RMF,sqrt
        import random
        from math import sqrt as sqrt
        from operator import itemgetter
        import IMP.pmi.tools as tools
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
        self.output_level="low"
        self.label="None"
  
        self.maxtrans_rb=0.15
        self.maxrot_rb=0.03
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
        self.hier_geometry_pairs={}
        self.elements={}

    # create a protein, represented as a set of connected balls of appropriate
    # radii and number, chose by the resolution parameter and the number of
    # amino acids.


    def shuffle_configuration(self,bounding_box_length=300.,translate=True):
        "shuffle configuration, used to restart the optimization"
        "it only works if rigid bodies were initialized"
        if len(self.rigid_bodies)==0:
            print "MultipleStates: rigid bodies were not intialized"
        hbbl=bounding_box_length/2
        if 1:
            for rb in self.rigid_bodies:
                ub = IMP.algebra.Vector3D(-hbbl,-hbbl,-hbbl)
                lb = IMP.algebra.Vector3D( hbbl, hbbl, hbbl)
                bb = IMP.algebra.BoundingBox3D(ub, lb)
                if translate==True: translation = IMP.algebra.get_random_vector_in(bb)
                else: translation = (rb.get_x(), rb.get_y(), rb.get_z())
                rotation = IMP.algebra.get_random_rotation_3d()
                transformation = IMP.algebra.Transformation3D(rotation, translation)
                rb.set_reference_frame(IMP.algebra.ReferenceFrame3D(transformation))


    def add_component_name(self,name):
        protein_h = IMP.atom.Molecule.setup_particle(IMP.Particle(self.m))
        protein_h.set_name(name)
        self.hier_dict[name]=protein_h
        self.prot.add_child(protein_h)
        self.elements[name]=[]

    
    def add_component_pdb(self,name,pdbname,chain,resolutions,color,resrange=None,offset=0,show=False,isnucleicacid=False):
        #resrange specify the residue range to extract from the pdb
        #it is a tuple (beg,end). If not specified, it takes all residues belonging to
        # the specified chain.
         
        protein_h=self.hier_dict[name]
        
        t=IMP.atom.read_pdb( pdbname, self.m, 
        IMP.atom.AndPDBSelector(IMP.atom.ChainPDBSelector(chain),IMP.atom.ATOMPDBSelector()))

        #find start and end indexes
        
        start = IMP.atom.Residue(t.get_children()[0].get_children()[0]).get_index()
        end   = IMP.atom.Residue(t.get_children()[0].get_children()[-1]).get_index()

        c=IMP.atom.Chain(IMP.atom.get_by_type(t, IMP.atom.CHAIN_TYPE)[0])

        
        if resrange!=None:
           if resrange[0]>start: start=resrange[0]
           if resrange[1]<end:   end=resrange[1] 
        
        sel=IMP.atom.Selection(c,residue_indexes=range(start,end+1))
        
        '''
        if not isnucleicacid:
           #do what you have to do for proteins
           ,atom_type=IMP.atom.AT_CA)

        else:
           #do what you have to do for nucleic-acids
           sel=IMP.atom.Selection(c,residue_indexes=range(start,end+1),atom_type=IMP.atom.AT_P)
        '''
        
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
        
        
        if show:
           IMP.atom.show_molecular_hierarchy(c0)
        
        
        for r in resolutions:
            if r==0:
               #use atomistic representation
               s=c0
            else:
               s=IMP.atom.create_simplified_along_backbone(c0, r)
            chil=s.get_children()
            s0=IMP.atom.Hierarchy.setup_particle(IMP.Particle(self.m))
            s0.set_name(name+'_%i-%i' % (start,end)+"_Res:"+str(r))
            for ch in chil: s0.add_child(ch)            
            protein_h.add_child(s0)
            del s
            for prt in IMP.atom.get_leaves(s0):
                IMP.pmi.Resolution.setup_particle(prt,r)
                #setting up color for each particle in the hierarchy, if colors missing in the colors list set it to red
                try:
                    clr=IMP.display.get_rgb_color(color)
                except:
                    colors.append(1.0)
                    clr=IMP.display.get_rgb_color(colors[pdb_part_count])
                IMP.display.Colored.setup_particle(prt,clr)


    def add_component_beads(self,name,ds,colors):
        from math import pi
        protein_h=self.hier_dict[name]    
        for n,dss in enumerate(ds):
            ds_frag=(dss[0],dss[1])
            self.elements[name].append((dss[0],dss[1]," ","bead"))
            prt=IMP.Particle(self.m)
            h=IMP.atom.Fragment.setup_particle(prt)
            h.set_residue_indexes(range(ds_frag[0],ds_frag[1]+1))
            h.set_name(name+'_%i-%i' % (ds_frag[0],ds_frag[1]))
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
            self.floppy_bodies.append(p)
            protein_h.add_child(h)

    def setup_component_geometry(self,name,color=0):
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
            self.hier_geometry_pairs[name].append((sortedparticles[n][0],sortedparticles[n+1][0]))
        
    def setup_component_sequence_connectivity(self,name):
        unmodeledregions_cr=IMP.RestraintSet("unmodeledregions")
        sortedsegments_cr=IMP.RestraintSet("sortedsegments")    
        protein_h=protein_h=self.hier_dict[name]    
        SortedSegments = []
        pbr=tools.get_particles_by_resolution(protein_h,10.0)
        
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
            r=IMP.core.PairRestraint(dps,IMP.ParticlePair(pt0,pt1))
            print "Adding sequence connectivity restraint between", pt0.get_name(), " and ", pt1.get_name()
            sortedsegments_cr.add_restraint(r)
            self.connected_intra_pairs.append((pt0,pt1))
            self.connected_intra_pairs.append((pt1,pt0))

        self.m.add_restraint(sortedsegments_cr)
        self.m.add_restraint(unmodeledregions_cr)
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

    def create_components_from_rmf(self,rmfname,frameindex):
        '''
        still not working.
        create the representation (i.e. hierarchies) from the rmf file.
        it will be stored in self.prot, which will be overwritten.
        load the coordinates from the rmf file at frameindex.
        '''
        rh= RMF.open_rmf_file(rmfname)
        self.prot=imprmf.create_hierarchies(rh, self.m)[0]
        imprmf.link_hierarchies(rh, [self.prot])
        imprmf.load_frame(rh, frameindex)
        '''
        still missing: save rigid bodies contained in the rmf in self.rigid_bodies
        save floppy bodies in self.floppy_bodies
        get the connectivity restraints
        '''


    def set_rigid_bodies(self,subunits,coords=None,nonrigidmembers=True):
        if coords==None: coords=()
        #sometimes, we know about structure of an interaction
        #and here we make such PPIs rigid
        randomize_coords = lambda c: tuple(1.*(nrrand(3)-0.5)+array(c))
        
        rigid_parts=[]
        for s in subunits:
            if type(s)==type(tuple()) and len(s)==2:
               sel=IMP.atom.Selection(self.prot,molecule=s[0],residue_indexes=range(s[1][0],s[1][1]+1))
               
               ps=sel.get_selected_particles()
               for p in ps: 
                   print s,p
               rigid_parts+=sel.get_selected_particles()
               
               
            elif type(s)==type(str()):
               sel=IMP.atom.Selection(self.prot,molecule=s)
               rigid_parts+=sel.get_selected_particles()
               ps=sel.get_selected_particles()
               for p in ps: 
                   print s,p        
        
        
        rb=IMP.atom.create_rigid_body(rigid_parts)
        rb.set_coordinates_are_optimized(True)
        rb.set_name(''.join(str(subunits))+"_rigid_body")        
        if type(coords)==tuple and len(coords)==3: rb.set_coordinates(randomize_coords(coords))
        self.rigid_bodies.append(rb)
            
        '''
        if type(subunits[0])==str:
            rigid_parts = []
            for prt in self.prot.get_children():
                if prt.get_name() in subunits:
                    for frag in prt.get_children(): rigid_parts += IMP.atom.get_leaves(frag)
            
            if not nonrigidmembers:
               for p in rigid_parts:
                   if p in self.floppy_bodies:
                      rigid_parts.remove(p)
                    
            rb=IMP.atom.create_rigid_body(rigid_parts)
            rb.set_coordinates_are_optimized(True)
            rb.set_name(''.join(subunits)+"_rigid_body")
            if type(coords)==tuple and len(coords)==3: rb.set_coordinates(randomize_coords(coords))
            self.rigid_bodies.append(rb)

        elif type(subunits[0])==tuple or type(subunits[0])==list():
            rigid_parts = []
            #print '#####',subunits, [name[0] for name in subunits]
            for subunit in subunits:
                    print "WWWWW",subunit
                    name=subunit[0]
                    bounds=subunit[1]
                    s= IMP.atom.Selection(self.prot,molecule=name, residue_indexes=range(bounds[0],bounds[1]+1))
                    rigid_parts += s.get_selected_particles()
                    
                    for p in s.get_selected_particles():
                        print p, IMP.atom.Fragment(p).get_parent()
                    
                    #print prt,s,'\n\t',s.get_selected_particles()
                    #for f in prt.get_children(): print '\t\t',f
                    #print
            
            if not nonrigidmembers:
               for p in rigid_parts:
                   if p in self.floppy_bodies:
                      rigid_parts.remove(p)        
            
            rb=IMP.atom.create_rigid_body(rigid_parts)
            rb.set_coordinates_are_optimized(True)
            rb.set_name(''.join(str(subunits))+"_rigid_body")
            if type(coords)==tuple and len(coords)==3: rb.set_coordinates(randomize_coords(coords))
            self.rigid_bodies.append(rb)
         '''

    def set_floppy_bodies(self):
        for p in self.floppy_bodies:
            name=p.get_name()
            p.set_name(name+"_floppy_body")
            tools.set_floppy_body(p)
    
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

    def get_particles_to_sample(self):
        #get the list of samplable particles with their type
        #and the mover displacement. Everything wrapped in a dictionary,
        #to be used by samplers modules
        ps={}
        
        #remove symmetric particles: they are not sampled
        rbtmp=[]
        fbtmp=[]
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
        return ps
    
    def set_output_level(self,level):
        self.output_level=level
    
    def get_output(self):
        output={}
        score=0.0
        
        output["SimplifiedModel_Total_Score_"+self.label]=str(self.m.evaluate(False))        
        for name in self.sortedsegments_cr_dict:
            partialscore=self.sortedsegments_cr_dict[name].evaluate(False)
            score+=partialscore
            output["SimplifiedModel_Link_SortedSegments_"+name+"_"+self.label]=str(partialscore)
            partialscore=self.unmodeledregions_cr_dict[name].evaluate(False)
            score+=partialscore            
            output["SimplifiedModel_Link_UnmodeledRegions_"+name+"_"+self.label]=str(partialscore)
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
            tools.get_graph_from_hierarchy(c)
            

    def get_geometries(self):
        #create segments at the lowest resolution
        seggeos=[]
        for name in self.hier_geometry_pairs:
            for pt in self.hier_geometry_pairs[name]:
                p1=pt[0]
                p2=pt[1]
                coor1=IMP.core.XYZ(p1).get_coordinates()
                coor2=IMP.core.XYZ(p2).get_coordinates()
                seg=IMP.algebra.Segment3D(coor1,coor2)
                seggeos.append(IMP.display.SegmentGeometry(seg))
        return seggeos
        
    def draw_hierarchy_composition(self):
        from matplotlib import pyplot
        import matplotlib as mpl
        


        
        ks=self.elements.keys()
        ks.sort()
 
        max=0
        for k in self.elements:
            for l in self.elements[k]:
                if l[1]>max: max=l[1]
        
        
        for k in ks:
            list=sorted(self.elements[k], key=itemgetter(0))
            endres=list[-1][1]
            fig = pyplot.figure(figsize=(26.0*float(endres)/max+2,2))
            ax = fig.add_axes([0.05, 0.475, 0.9, 0.15])
            #ax = fig.add_axes([0.05, 0.475, 0.9, 0.15])
            
            # Set the colormap and norm to correspond to the data for which
            # the colorbar will be used.
            cmap = mpl.cm.cool
            norm = mpl.colors.Normalize(vmin=5, vmax=10)
    

    
            
            
            bounds=[1]
            colors=['white']
            
            
            for l in list:
                firstres=l[0]
                lastres=l[1]
                if l[3]=="pdb": colors.append("#99CCFF")
                if l[3]=="bead": colors.append("#FFFF99")                
                bounds.append(l[0])

            bounds.append(endres)
            
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