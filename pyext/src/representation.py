#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.display

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
        ps["Rigid_Bodies_MultipleStates"]=(rblist,self.maxtrans_rb,self.maxrot_rb)
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
#Peter Cimermancic
    '''
    This class creates the molecular hierarchies for the various involved proteins.
    '''

    def __init__(self,m,upperharmonic=True):
        global random, itemgetter,tools,nrrand,array
        import random
        from operator import itemgetter
        import IMP.pmi.tools as tools
        from numpy.random import rand as nrrand
        from numpy import array
        
        # this flag uses either harmonic (False) or upperharmonic (True)
        # in the intra-pair connectivity restraint. Harmonic is used whe you want to
        # remove the intra-ev term from energy calculations, e.g.:
        # upperharmonic=False
        # ip=simo.get_connected_intra_pairs()
        # ev.add_excluded_particle_pairs(ip)
        
        self.upperharmonic=upperharmonic
        self.rigid_bodies=[]
        self.floppy_bodies=[]
                 
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

    # create a protein, represented as a set of connected balls of appropriate
    # radii and number, chose by the resolution parameter and the number of
    # amino acids.


    def shuffle_configuration(self,bounding_box_length=100.,translate=True):
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



    def add_component(self,name,chainnames, length, pdbs, init_coords=None,simplepdb=1,ds=None,colors=None):
        
        if ds==None: ds=[]
        if colors==None: colors=[]
        if init_coords==None: init_coords=()
        
        
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

                #find start and end indeces
                start = IMP.atom.Residue(t.get_children()[0].get_children()[0]).get_index()
                end   = IMP.atom.Residue(t.get_children()[0].get_children()[-1]).get_index()
                bounds.append(( start,end ))

                #IMP.atom.show_molecular_hierarchy(t)
                c=IMP.atom.Chain(IMP.atom.get_by_type(t, IMP.atom.CHAIN_TYPE)[0])
                if c.get_number_of_children()==0:
                    IMP.atom.show_molecular_hierarchy(t)
                # there is no reason to use all atoms, just approximate the pdb shape instead
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
                protein_h.add_child(s)
                
                for prt in IMP.atom.get_leaves(s):
                    #setting up color for each particle in the hierarchy, if colors missing in the colors list set it to red
                    try:
                       clr=IMP.display.get_rgb_color(colors[pdb_part_count])
                    except:
                       clr=IMP.display.get_rgb_color(1.0)
                    IMP.display.Colored.setup_particle(prt,clr)  
                
                s.set_name(pdb)



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

        print ds
        #work on un-modelled regions
        randomize_coords = lambda c: tuple(10.*(nrrand(3)-0.5)+array(c))

        for n,ds_frag in enumerate(ds):
            if ds_frag[1]-ds_frag[0]==0: ds_frag=(ds_frag[0],ds_frag[1]+1)
            h=IMP.atom.create_protein(self.m, name+'_%i-%i' % (ds_frag[0],ds_frag[1]), self.resolution, \
                          [ds_frag[0],ds_frag[1]])
            for prt in IMP.atom.get_leaves(h):
                #setting up color for the particle, if colors missing in the colors list set it to red
                try:
                    clr=IMP.display.get_rgb_color(colors[n])
                except:
                    clr=IMP.display.get_rgb_color(1.0)
                IMP.display.Colored.setup_particle(prt,clr)


                ptem= prt.get_as_xyzr()
                ptem.set_radius(ptem.get_radius()*0.8)

                #if init_coords==(): ptem.set_coordinates((random.uniform(bb[0][0],bb[1][0]),\
                #      random.uniform(bb[0][1],bb[1][1]),random.uniform(bb[0][2],bb[1][2])))
                #else: ptem.set_coordinates(randomize_coords(init_coords))
                if len(pdbs)==0: ptem.set_coordinates((random.uniform(bb[0][0],bb[1][0]),\
                       random.uniform(bb[0][1],bb[1][1]),random.uniform(bb[0][2],bb[1][2])))
                else:
                    crd=IMP.atom.get_leaves(s)[0].get_as_xyz().get_coordinates()
                    ptem.set_coordinates(randomize_coords(crd))

                FlexParticles.append(ptem.get_particle())

            Spheres = IMP.atom.get_leaves(h)
            for x,prt in enumerate(Spheres):
                if x==len(Spheres)-1: break
                #r=IMP.atom.create_connectivity_restraint([IMP.atom.Selection(prt),\
                #                          IMP.atom.Selection(Spheres[x+1])],-3.0,self.kappa)
                if self.upperharmonic:
                   hu=IMP.core.HarmonicUpperBound(0., self.kappa)
                else:
                   hu=IMP.core.Harmonic(0., self.kappa)                   
                dps=IMP.core.SphereDistancePairScore(hu)
                pt0=IMP.atom.Selection(prt).get_selected_particles()[0]
                pt1=IMP.atom.Selection(Spheres[x+1]).get_selected_particles()[0]
                r=IMP.core.PairRestraint(dps,IMP.ParticlePair(pt0,pt1))
                unmodeledregions_cr.add_restraint(r)
                self.connected_intra_pairs.append((pt0,pt1))
                self.connected_intra_pairs.append((pt1,pt0))                
                        # only allow the particles to separate by 1 angstrom
                        # self.m.set_maximum_score(r, self.kappa)
            protein_h.add_child(h)


        SortedSegments = []
        for chl in protein_h.get_children():
            start = IMP.atom.get_leaves(chl)[0]
            end   = IMP.atom.get_leaves(chl)[-1]

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
            if self.upperharmonic:
               hu=IMP.core.HarmonicUpperBound(0., self.kappa)
            else:
               hu=IMP.core.Harmonic(0., self.kappa)    
            dps=IMP.core.SphereDistancePairScore(hu)
            pt0=IMP.atom.Selection(last).get_selected_particles()[0]
            pt1=IMP.atom.Selection(first).get_selected_particles()[0]
            r=IMP.core.PairRestraint(dps,IMP.ParticlePair(pt0,pt1))
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

    def set_rigid_bodies(self,subunits,coords=()):
        #sometimes, we know about structure of an interaction
        #and here we make such PPIs rigid
        randomize_coords = lambda c: tuple(1.*(nrrand(3)-0.5)+array(c))

        if type(subunits[0])==str:
            rigid_parts = []
            for prt in self.prot.get_children():
                if prt.get_name() in subunits:
                    for frag in prt.get_children(): rigid_parts += IMP.atom.get_leaves(frag)
            rb=IMP.atom.create_rigid_body(rigid_parts)
            rb.set_coordinates_are_optimized(True)
            if coords!=(): rb.set_coordinates(randomize_coords(coords))
            else: rb.set_coordinates(randomize_coords((0.,0.,0.)))
            self.rigid_bodies.append(rb)

        elif type(subunits[0])==tuple or type(subunits[0])==list():
            rigid_parts = []
            print '#####',subunits, [name[0] for name in subunits]
            for prt in self.prot.get_children():
                if prt.get_name() in [name[0] for name in subunits]:
                    prt_index=[name[0] for name in subunits].index(prt.get_name())
                    bounds=subunits[prt_index][1]
                    s= IMP.atom.Selection(prt, residue_indexes=range(bounds[0],bounds[1]+1))
                    rigid_parts += s.get_selected_particles()
                    print prt,s,'\n\t',s.get_selected_particles()
                    for f in prt.get_children(): print '\t\t',f
                    print
            rb=IMP.atom.create_rigid_body(rigid_parts)
            rb.set_coordinates_are_optimized(True)
            if coords!=(): rb.set_coordinates(randomize_coords(coords))
            else: rb.set_coordinates(randomize_coords((0.,0.,0.)))
            self.rigid_bodies.append(rb)


    def set_floppy_bodies(self):
        for p in self.floppy_bodies:
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
        ps["Rigid_Bodies_SimplifiedModel"]=(self.rigid_bodies,self.maxtrans_rb,self.maxrot_rb)
        ps["Floppy_Bodies_SimplifiedModel"]=(self.floppy_bodies,self.maxtrans_fb)
        return ps

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
        output["_TotalScore"]=str(score)
        return output

    def get_hierarchy(self):
        return  self.prot
