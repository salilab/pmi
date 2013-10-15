#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container


class LinkDomains():

    def __init__(self,prot,resrangelist,kappa,length=5.0):
        #generate a linker between residues using HarmonicUpperBound
        #restraints. Define a list of linked residues,
        #e.g.   [(3,5,"A"),(9,10,"B")]
        # will link residues 3 and 5 belonging to chain A and
        #residues 9 and 10 belonging to chain B
        self.rs = IMP.RestraintSet('linker')
        self.prot=prot
        self.kappa=kappa
        self.resrangelist=resrangelist
        self.length=length
        self.label="None"

        self.m=self.prot.get_model()
        self.pairs=[]

        for pair in self.resrangelist:
            c0=pair[2]
            r0=pair[0]
            c1=pair[2]
            r1=pair[1]
            try:
                s0=IMP.atom.Selection(self.prot, chains=c0, residue_index=r0,atom_type=IMP.atom.AT_CA)
                p0=s0.get_selected_particles()[0]
            except:
                "LinkDomains: error"
                continue
            try:
                s1=IMP.atom.Selection(self.prot, chains=c1, residue_index=r1,atom_type=IMP.atom.AT_CA)
                p1=s1.get_selected_particles()[0]
            except:
                "LinkDomains: error"
                continue
            #check this is the residue length (should be 4, but I use a larger length)
            dist0=float(pair[1]-pair[0])*self.length
            h=IMP.core.HarmonicUpperBound(dist0, self.kappa)
            dps=IMP.core.DistancePairScore(h)
            pr=IMP.core.PairRestraint(dps,IMP.ParticlePair(p0,p1))
            pr.set_name("LinkDomains_"+str(pair[0])+"-"+str(pair[1])+"_"+str(pair[2]))
            self.rs.add_restraint(pr)
            self.pairs.append((p0,p1,r0,c0,r1,c1,pr))

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_pairs(self):
        return self.pairs

    def get_hierarchy(self):
        return self.prot

    def get_kappa(self):
        return self.kappa

    def get_residue_ranges(self):
        return self.resrangelist

    def get_length(self):
        return self.length

    def get_restraint(self):
        return self.rs

    def get_restraints(self):
        rlist=[]
        for r in self.rs.get_restraints():
            rlist.append(IMP.core.PairRestraint.get_from(r))
        return rlist

    def print_chimera_pseudobonds(self,filesuffix,model=0):
        f=open(filesuffix+".chimera","w")
        atype="ca"
        for p in self.pairs:
            s="#"+str(model)+":"+str(p[2])+"."+p[3]+"@"+atype+" #"+str(model)+":"+str(p[4])+"."+p[5]+"@"+atype
            f.write(s+"\n")
        f.close()

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)
        output["LinkDomains_"+self.label]=str(score)
        for rst in self.rs.get_restraints():
            #output["LinkDomains_"+
            #        IMP.core.PairRestraint.get_from(rst).get_name()+
            #           "_"+self.label]=IMP.core.PairRestraint.get_from(rst).evaluate(False)
            output["LinkDomains_"+rst.get_name()+
                       "_"+self.label]=rst.unprotected_evaluate(None)

        for i in range(len(self.pairs)):

            p0=self.pairs[i][0]
            p1=self.pairs[i][1]
            r0=self.pairs[i][2]
            c0=self.pairs[i][3]
            r1=self.pairs[i][4]
            c1=self.pairs[i][5]

            label=str(r0)+":"+c0+"_"+str(r1)+":"+c1

            d0=IMP.core.XYZ(p0)
            d1=IMP.core.XYZ(p1)
            output["LinkDomains_Distance_"+label+"_"+self.label]=str(IMP.core.get_distance(d0,d1))

        return output


###########################################################################

class ConnectivityRestraint():

    def __init__(self,hier,selection_tuples,kappa=10.0,resolution=None,label="None"):
        '''
        generate a connectivity restraint between domains
        setting up the composite restraint
        example:
        cr=restraints.ConnectivityRestraint(prot,["CCC",(1,100,"TTT"),(100,150,"AAA")])
        cr.add_to_model()
        cr.set_label("CR1")
        '''
        self.weight=1.0
        self.hier=hier
        self.kappa=kappa
        self.label=label
        if self.label=="None": self.label=str(selection_tuples)
        self.rs = IMP.RestraintSet(label)
        
        self.m=self.hier.get_model()
        
        sels=[]
        

        for s in selection_tuples:

                   
           if type(s)==tuple and len(s)==3:
              hiers=[]
              for h in self.hier.get_children():
                if s[2] in h.get_name():
                   hiers.append(h)
              sel=IMP.atom.Selection(hierarchies=hiers,residue_indexes=range(s[0],s[1]+1))
           elif type(s)==str:
              hiers=[]
              for h in self.hier.get_children():
                if s in h.get_name():
                   hiers.append(h)
              sel=IMP.atom.Selection(hierarchies=hiers)  

           if resolution!=None:
              particles=[]
              for prot in self.hier.get_children():
                  particles+=IMP.pmi.tools.get_particles_by_resolution(prot,resolution)              
              selectedp=sel.get_selected_particles()
              #get the intersection to remove redundant particles
              sel=IMP.atom.Selection(list(set(selectedp) & set(particles)))
           sels.append(sel)
                        
        cr = IMP.atom.create_connectivity_restraint(sels, self.kappa, self.label)
        self.rs.add_restraint(cr)

    def set_label(self,label):
        self.label=label
        self.rs.set_name(label)
        for r in self.rs.get_restraints():
            r.set_name(label)

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_restraint(self):
        return self.rs

    def get_restraints(self):
        rlist=[]
        for r in self.rs.get_restraints():
            rlist.append(IMP.core.PairRestraint.get_from(r))
        return rlist

    def set_weight(self,weight):
        self.weight=weight
        self.rs.set_weight(weight)
 
    def get_output(self):
        self.m.update()
        output={}
        score=self.weight*self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)        
        output["ConnectivityRestraint_"+self.label]=str(score)
        return output


###########################################################################


class UpperBound():

    def __init__(self,prot,respairs,kappa,length=5.0):
        #generate a restraint between residues using HarmonicUpperBound
        #restraints. Define a list of linked residues,
        #e.g.   [(3,"A",5,"B"),(9,"B",10,"C")]
        # will link residues 3 and 5 belonging to chain A and B and
        #residues 9 and 10 belonging to chain B and C
        self.rs = IMP.RestraintSet('upperbound')

        self.prot=prot
        self.kappa=kappa
        self.respairs=respairs
        self.length=length
        self.label="None"

        self.m=self.prot.get_model()

        for pair in self.respairs:
            try:
                s0=IMP.atom.Selection(self.prot, chains=pair[1], residue_index=pair[0],atom_type=IMP.atom.AT_CA)
                p0=s0.get_selected_particles()[0]
            except:
                "UpperBound: error"
                continue
            try:
                s1=IMP.atom.Selection(self.prot, chains=pair[3], residue_index=pair[2],atom_type=IMP.atom.AT_CA)
                p1=s1.get_selected_particles()[0]
            except:
                "UpperBound: error"
                continue
            #check this is the residue length (should be 4, but I use a larger length)

            h=IMP.core.HarmonicUpperBound(self.length, self.kappa)
            dps=IMP.core.DistancePairScore(h)
            pr=IMP.core.PairRestraint(dps,IMP.ParticlePair(p0,p1))
            pr.set_name("UpperBound_"+str(pair[0])+"-"+str(pair[1])+"_"+str(pair[2]))
            self.rs.add_restraint(pr)

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_hierarchy(self):
        return self.prot

    def get_kappa(self):
        return self.kappa

    def get_residue_ranges(self):
        return self.respairs

    def get_length(self):
        return self.length

    def get_restraint(self):
        return self.rs

    def get_output(self):
        output={}
        self.m.update()
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)        
        output["UpperBound_"+self.label]=str(score)
        return output

###########################################################################


class ExcludedVolumeResidue():

    def __init__(self,prot,kappa):
        self.rs = IMP.RestraintSet('excluded_volume')
        self.prot=prot
        self.kappa=kappa
        self.label="None"
        self.m=self.prot.get_model()

        atoms=IMP.atom.get_by_type(self.prot, IMP.atom.ATOM_TYPE)
        for atom in atoms:
            restype=IMP.atom.Residue(IMP.atom.Atom(atom).get_parent()).get_residue_type()
            vol=IMP.atom.get_volume_from_residue_type(restype)
            radius=IMP.algebra.get_ball_radius_from_volume_3d(vol)
            IMP.core.XYZR(atom).set_radius(radius)
        lsa=IMP.container.ListSingletonContainer(self.m)
        lsa.add_particles(atoms)

        evr=IMP.core.ExcludedVolumeRestraint(lsa,self.kappa)
        self.rs.add_restraint(evr)

    def set_label(self,label):
        self.label=label

    def add_excluded_particle_pairs(self,excluded_particle_pairs):
        # add pairs to be filtered when calculating  the score
        lpc=IMP.container.ListPairContainer(self.m)
        lpc.add_particle_pairs(excluded_particle_pairs)
        icpf=IMP.container.InContainerPairFilter(lpc)
        IMP.core.ExcludedVolumeRestraint.get_from(self.rs.get_restraints()[0]).add_pair_filter(icpf)

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_hierarchy(self):
        return self.prot

    def get_kappa(self):
        return self.kappa

    def get_restraint(self):
        return self.rs

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)        
        output["ExcludedVolumeResidue_"+self.label]=str(score)
        return output

###########################################################################

class ExcludedVolumeSphere():

    def __init__(self,prot,resolution=None,kappa=1.0):
        self.rs = IMP.RestraintSet('excluded_volume')
        self.weight=1.0
        self.kappa=kappa
        self.prot=prot
        self.label="None"
        self.m=self.prot.get_model()
        
        if resolution==None:
          #default
          evr=IMP.atom.create_excluded_volume_restraint([prot])
        else:
           particles=[]
           lsa=IMP.container.ListSingletonContainer(self.m)
           for hier in prot.get_children():
              particles=IMP.pmi.tools.get_particles_by_resolution(hier,resolution)
              lsa.add_particles(particles)
           evr=IMP.core.ExcludedVolumeRestraint(lsa,self.kappa)   
           
                  
        self.rs.add_restraint(evr)

    def add_excluded_particle_pairs(self,excluded_particle_pairs):
        # add pairs to be filtered when calculating  the score
        lpc=IMP.container.ListPairContainer(self.m)
        lpc.add_particle_pairs(excluded_particle_pairs)
        icpf=IMP.container.InContainerPairFilter(lpc)
        IMP.core.ExcludedVolumeRestraint.get_from(self.rs.get_restraints()[0]).add_pair_filter(icpf)

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_hierarchy(self):
        return self.prot

    def get_kappa(self):
        return self.kappa

    def get_restraint(self):
        return self.rs

    def set_weight(self,weight):
        self.weight=weight
        self.rs.set_weight(weight)

    def get_output(self):
        self.m.update()
        output={}
        score=self.weight*self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)        
        output["ExcludedVolumeSphere_"+self.label]=str(score)
        return output

###########################################################################


class BipartiteExcludedVolumeResidue():

    def __init__(self,prot1,prot2,kappa):
        self.rs = IMP.RestraintSet('bipartite_excluded_volume')
        self.prot1=prot1
        self.prot2=prot2
        self.kappa=kappa
        self.label="None"
        self.m=self.prot.get_model()

        atoms1=IMP.atom.get_by_type(prot1, IMP.atom.ATOM_TYPE)
        ls1=IMP.container.ListSingletonContainer(self.m)
        ls1.add_particles(atoms1)
        for atom in atoms1:
            restype=IMP.atom.Residue(IMP.atom.Atom(atom).get_parent()).get_residue_type()
            vol=IMP.atom.get_volume_from_residue_type(restype)
            radius=IMP.algebra.get_ball_radius_from_volume_3d(vol)
            IMP.core.XYZR(atom).set_radius(radius)

        atoms2=IMP.atom.get_by_type(prot2, IMP.atom.ATOM_TYPE)
        ls2=IMP.container.ListSingletonContainer(self.m)
        ls2.add_particles(atoms2)
        for atom in atoms2:
            restype=IMP.atom.Residue(IMP.atom.Atom(atom).get_parent()).get_residue_type()
            vol=IMP.atom.get_volume_from_residue_type(restype)
            radius=IMP.algebra.get_ball_radius_from_volume_3d(vol)
            IMP.core.XYZR(atom).set_radius(radius)

        cbpc=IMP.container.CloseBipartitePairContainer(ls_ref,ls_symm,kappa,10.0)
        ssps=IMP.core.SoftSpherePairScore(kappa)
        evr3=IMP.container.PairsRestraint(ssps,cbpc)
        self.rs.add_restraint(evr3)
        self.m.add_restraint(rs)

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_hierarchy(self):
        return self.prot1, self.prot2

    def get_kappa(self):
        return self.kappa

    def get_restraint(self):
        return self.rs

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)        
        output["BipartiteExcludedVolumeResidue_"+self.label]=str(score)
        return output

###########################################################################

class TemplateRestraint():
    def __init__(self,ps1,ps2,cutoff=6.5,kappa=1.0,forcerb=False):
        self.m=ps1[0].get_model()
        self.label="None"
        self.cutoff=cutoff
        self.kappa=kappa
        #this parameter ovverides the rigid body filter below
        self.forcerb=forcerb
        self.rs=IMP.RestraintSet('template_restraint')
        for p1 in  ps1:
            for p2 in ps2:
                #check that the two particles are not in the same rigid body
                if(IMP.core.RigidMember.particle_is_instance(p1) and IMP.core.RigidMember.particle_is_instance(p2) and
                IMP.core.RigidMember(p1).get_rigid_body() == IMP.core.RigidMember(p2).get_rigid_body()) and not self.forcerb: continue
                d0=IMP.core.XYZ(p1)
                d1=IMP.core.XYZ(p2)
                dist=IMP.core.get_distance(d0,d1)
                if dist <= self.cutoff:
                    hf=IMP.core.Harmonic(dist,self.kappa)
                    dps=IMP.core.DistancePairScore(hf)
                    pr=IMP.core.PairRestraint(dps,IMP.ParticlePair(p1,p2))
                    self.rs.add_restraint(pr)

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_cutoff(self):
        return self.cutoff

    def get_kappa(self):
        return self.kappa

    def get_restraint(self):
        return self.rs

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)        
        output["TemplateRestraint_"+self.label]=str(score)
        return output


###########################################################################

class CompositeRestraint():

    def __init__(self,handleparticles,compositeparticles,cut_off=5.0,lam=1.0,label="None"):
        
        global imppmi
        import IMP.pmi as imppmi 
        #composite particles: all particles beside the handle
        self.label=label
        self.rs = IMP.RestraintSet('cr')
        self.m=handleparticles[0].get_model()       
        
        print handleparticles
        
        ln=imppmi.CompositeRestraint(self.m,handleparticles,cut_off,lam,True)
        for ps in compositeparticles:
            #composite particles is a list of list of particles
            ln.add_composite_particle(ps)
        
        self.rs.add_restraint(ln)
      
    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)        
        output["CompositeRestraint_"+self.label]=str(score)
        return output

###########################################################################

class MarginalChi3Restraint():

    def __init__(self,part1,part2):
        global impisd2, tools
        import IMP.isd2 as impisd2
        import IMP.pmi.tools as tools

        self.m=part1.get_model()
        self.label="None"
        self.rs=IMP.RestraintSet('chi3_restraint')
        self.sigmamaxtrans=0.1

        self.ps1=IMP.atom.get_leaves(part1)
        self.ps2=IMP.atom.get_leaves(part2)
        self.sigma=tools.SetupNuisance(self.m,1.0,0.1,100.0,True).get_particle()

        for i in range(len(self.ps1)):
            mc=impisd2.MarginalChi3Restraint(self.ps1[i],self.ps2[i],self.sigma)
            self.rs.add_restraint(mc)

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_restraint(self):
        return self.rs

    def get_particles_to_sample(self):
        ps={}
        ps["Nuisances_MarginalChi3Restraint_Sigma_"+self.label]=([self.sigma],self.sigmamaxtrans)
        return ps

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)        
        output["MarginalChi3Restraint_"+self.label]=str(score)
        output["MarginalChi3Restraint_Sigma_"+self.label]=str(self.sigma.get_scale())
        return output

###########################################################################



class ExternalBarrier():

    def __init__(self,prot,radius):
        self.rs = IMP.RestraintSet('barrier')
        self.prot=prot
        self.radius=radius
        self.label="None"

        self.m=self.prot.get_model()
        c3= IMP.algebra.Vector3D(0,0,0)
        ub3= IMP.core.HarmonicUpperBound(radius, 10.0)
        ss3= IMP.core.DistanceToSingletonScore(ub3, c3)
        lsc= IMP.container.ListSingletonContainer(self.m)
        #IMP.atom.get_by_type

        lsc.add_particles(IMP.atom.get_leaves(self.prot))
        r3= IMP.container.SingletonsRestraint(ss3, lsc)
        self.rs.add_restraint(r3)

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_hierarchy(self):
        return self.prot

    def get_radius(self):
        return self.radius

    def get_restraint(self):
        return self.rs

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)        
        output["ExternalBarrier_"+self.label]=str(score)
        return output

###########################################################################

class CrossLinkMS():
    '''
    this class initialize a CrossLinkMS restraint and contains
    all useful informations, such as the cross-link database, contained in self.pairs
    If restraint_file=None, it will proceed creating simulated data
    '''
    def __init__(self,prots,
                 listofxlinkertypes=["BS3","BS2G","EGS"],map_between_protein_names_and_chains=None,
                 sigmamin=1.0,sigmamax=1.0,sigmagrid=1,sigmaissampled=False,typeofprofile="gofr"):
        global impisd2,tools
        import IMP.isd2 as impisd2
        import IMP.pmi.tools as tools
        
        if map_between_protein_names_and_chains==None:
           map_between_protein_names_and_chains={}

        self.rs=IMP.RestraintSet('data')
        self.rs2=IMP.RestraintSet('prior')
        self.prots=prots
        self.label="None"
        self.pairs=[]
        self.m=self.prots[0].get_model()
        self.sigmamin=sigmamin
        self.sigmamax=sigmamax
        self.sigmagrid=sigmagrid
        self.sigmaissampled=sigmaissampled

        self.sigmatrans=0.1
        self.sigmaglobal=tools.SetupNuisance(self.m,self.sigmamin,
                 self.sigmamin,self.sigmamax,self.sigmaissampled).get_particle()
        self.outputlevel="low"
        self.listofxlinkertypes=listofxlinkertypes
        #simulated reaction rates
        self.reaction_rates=None
        self.allpairs_database=None
        self.residue_list=None

        self.crosslinker_dict=self.get_crosslinker_dict(typeofprofile)
        #this map is used in the new cross-link file reader
        #{"Nsp1":"A","Nup82":"B"} etc.
        self.mbpnc=map_between_protein_names_and_chains
        #check whther the file was initialized



    #-------------------------------

    def get_crosslinker_dict(self,typeofprofile="gofr"):
        #fill the cross-linker pmfs
        #to accelerate the init the list listofxlinkertypes might contain only yht needed crosslinks
        #type of profile can be gofr or pfes

        disttuple=(0.0, 200.0, 1000)
        omegatuple=(1.0, 1000.0, 30)
        sigmatuple=(self.sigmamin, self.sigmamax, self.sigmagrid)

        crosslinker_dict={}
        if "BS3" in self.listofxlinkertypes:
            crosslinker_dict["BS3"]=tools.get_cross_link_data("bs3l",
               "pmf_bs3l_tip3p.txt.standard",disttuple,omegatuple,sigmatuple,
                                            don=None,doff=None,prior=0,type_of_profile=typeofprofile)
        if "BS2G" in self.listofxlinkertypes:
            crosslinker_dict["BS2G"]=tools.get_cross_link_data("bs2gl",
              "pmf_bs2gl_tip3p.txt.standard",disttuple,omegatuple,sigmatuple,
                                            don=None,doff=None,prior=0,type_of_profile=typeofprofile)
        if "EGS" in self.listofxlinkertypes:
            crosslinker_dict["EGS"]=tools.get_cross_link_data("egl",
              "pmf_egl_tip3p.txt.standard",disttuple,omegatuple,sigmatuple,
                                            don=None,doff=None,prior=0,type_of_profile=typeofprofile)
        if "Short" in self.listofxlinkertypes:
            #setup a "short" xl with an half length of 10 Ang
            crosslinker_dict["Short"]=tools.get_cross_link_data_from_length(10.0,disttuple,omegatuple,sigmatuple)
        return crosslinker_dict

    #-------------------------------

    def add_restraints(self,restraint_file=None,oldfile=True):

        if restraint_file==None:
            #get the restraints from simulated data
            restraint_list=self.allpairs_database

        else:
            #get the restraints from external file
            f = open(restraint_file)
            restraint_list = f.readlines()

        self.index=0


        self.added_pairs_list=[]
        self.missing_residues=[]
        for line in restraint_list:
            #force_restraint=True makes all intra rigid body restraint to be accepted
            force_restraint=False

            if restraint_file==None:
                if line["Is_Detected"]==True:
                    crosslinker=line["Crosslinker"]
                    (r1,c1)=line["Identified_Pair1"]
                    (r2,c2)=line["Identified_Pair2"]
                    index+=1
                else:
                    continue

            elif oldfile:
                tokens=line.split()
                #skip character
                if (tokens[0]=="#"): continue
                r1=int(tokens[0])
                c1=tokens[1]
                r2=int(tokens[2])
                c2=tokens[3]
                crosslinker=tokens[4]

                #two restraints with the same index will be ambiguous
                self.index=int(tokens[5])

                #force restraint even if it belong to the same rigid body, use it for ambiguous restraints
                if (tokens[len(tokens)-1]=="F"): force_restraint=True

            else:
                #read with the new file parser
                totallist=eval(line)
                self.add_crosslink_according_to_new_file(totallist)
                #skip the rest
                continue


            print '''CrossLinkMS: attempting to add restraint between
                     residue %d of chain %s and residue %d of chain %s''' % (r1,c1,r2,c2)

            p1s=[]
            p2s=[]

            try:
                s1=IMP.atom.Selection(self.prots[0], residue_index=r1, chains=c1, atom_type=IMP.atom.AT_CA)
                p1=(s1.get_selected_particles()[0])
            except:
                print "CrossLinkMS: WARNING> residue %d of chain %s is not there" % (r1,c1)
                continue
            try:
                s2=IMP.atom.Selection(self.prots[0], residue_index=r2, chains=c2, atom_type=IMP.atom.AT_CA)
                p2=(s2.get_selected_particles()[0])
            except:
                print "CrossLinkMS: WARNING> residue %d of chain %s is not there" % (r2,c2)
                continue


            for copy in self.prots:
                s1=IMP.atom.Selection(copy, residue_index=r1, chains=c1, atom_type=IMP.atom.AT_CA)
                p1s.append(s1.get_selected_particles()[0])
                s2=IMP.atom.Selection(copy, residue_index=r2, chains=c2, atom_type=IMP.atom.AT_CA)
                p2s.append(s2.get_selected_particles()[0])

            #check whether the atom pair belongs to the same rigid body
            if(IMP.core.RigidMember.particle_is_instance(p1s[0]) and
               IMP.core.RigidMember.particle_is_instance(p2s[0]) and
               IMP.core.RigidMember(p1s[0]).get_rigid_body() ==
               IMP.core.RigidMember(p2s[0]).get_rigid_body() and not force_restraint):
                print '''CrossLinkMS: WARNING> residue %d of chain %s and
                       residue %d of chain %s belong to the same rigid body''' % (r1,c1,r2,c2)
                continue

            #this list contains the list of symmetric pairs to avoid duplicates
            if (p1s[0],p2s[0],crosslinker) in self.added_pairs_list:
                print "CrossLinkMS: WARNING> pair %d %s %d %s already there" % (r1,c1,r2,c2)
                continue
            if (p2s[0],p1s[0],crosslinker) in self.added_pairs_list:
                print "CrossLinkMS: WARNING> pair %d %s %d %s already there" % (r1,c1,r2,c2)
                continue

            print "CrossLinkMS: added pair %d %s %d %s" % (r1,c1,r2,c2)
            self.added_pairs_list.append((p1s[0],p2s[0],crosslinker))

            rs_name='restraint_'+str(index)

            ln=impisd2.CrossLinkMSRestraint(self.sigmaglobal,self.crosslinker_dict[crosslinker])
            for i in range(len(p1s)):
                ln.add_contribution(p1s[i],p2s[i])

            for i in range(len(self.prots)):
                self.pairs.append((p1s[i], p2s[i],
                                   crosslinker, rs_name,
                                   100, 100, (r1,c1,i), (r2,c2,i),
                                   crosslinker, i, ln))

            self.rs.add_restraint(ln)


        #self.rs2.add_restraint(impisd2.JeffreysRestraint(self.sigmaglobal))
        self.rs2.add_restraint(impisd2.UniformPrior(self.sigmaglobal,1000.0,self.sigmaglobal.get_upper()-1.0,self.sigmaglobal.get_lower()+0.1))
        print "CrossLinkMS: missing residues"
        for ms in self.missing_residues:
            print "CrossLinkMS:missing "+str(ms)








#---------------------------------
    def add_crosslink_according_to_new_file(self,totallist):
        force_restraint=False
        ambiguous_list=totallist[0]
        crosslinker=totallist[1]
        if (totallist[2]=="F"): force_restraint=True

        p1s=[]
        p2s=[]
        r1s=[]
        r2s=[]
        c1s=[]
        c2s=[]
        self.index+=1
        for pair in ambiguous_list:
            error=False

            try:
                c1=self.mbpnc[pair[0][0]]
            except:
                "CrossLinkMS: WARNING> protein name "+pair[0][0]+" was not defined"
                continue
            try:
                c2=self.mbpnc[pair[1][0]]
            except:
                "CrossLinkMS: WARNING> protein name "+pair[1][0]+" was not defined"
                continue
            r1=int(pair[0][1])
            r2=int(pair[1][1])

            print '''CrossLinkMS: attempting to add restraint between
                     residue %d of chain %s and residue %d of chain %s''' % (r1,c1,r2,c2)


            try:
                s1=IMP.atom.Selection(self.prots[0], residue_index=r1, chains=c1, atom_type=IMP.atom.AT_CA)
                p1=(s1.get_selected_particles()[0])
            except:
                print "CrossLinkMS: WARNING> residue %d of chain %s is not there" % (r1,c1)
                error=True
                self.missing_residues.append((r1,c1))
            try:
                s2=IMP.atom.Selection(self.prots[0], residue_index=r2, chains=c2, atom_type=IMP.atom.AT_CA)
                p2=(s2.get_selected_particles()[0])
            except:
                print "CrossLinkMS: WARNING> residue %d of chain %s is not there" % (r2,c2)
                error=True
                self.missing_residues.append((r2,c2))
            if error: continue


            s1=IMP.atom.Selection(self.prots[0], residue_index=r1, chains=c1, atom_type=IMP.atom.AT_CA)
            p1=s1.get_selected_particles()[0]
            s2=IMP.atom.Selection(self.prots[0], residue_index=r2, chains=c2, atom_type=IMP.atom.AT_CA)
            p2=s2.get_selected_particles()[0]
            #this list contains the list of symmetric pairs to avoid duplicates
            if (p1,p2,crosslinker) in self.added_pairs_list:
                print "CrossLinkMS: WARNING> pair %d %s %d %s already there" % (r1,c1,r2,c2)
                continue
            if (p2,p1,crosslinker) in self.added_pairs_list:
                print "CrossLinkMS: WARNING> pair %d %s %d %s already there" % (r1,c1,r2,c2)
                continue

            #check whether the atom pair belongs to the same rigid body
            if(IMP.core.RigidMember.particle_is_instance(p1) and
               IMP.core.RigidMember.particle_is_instance(p2) and
               IMP.core.RigidMember(p1).get_rigid_body() ==
               IMP.core.RigidMember(p2).get_rigid_body() and not force_restraint):
                print '''CrossLinkMS: WARNING> residue %d of chain %s and
                       residue %d of chain %s belong to the same rigid body''' % (r1,c1,r2,c2)
                continue

            p1s.append(p1)
            p2s.append(p2)
            r1s.append(r1)
            r2s.append(r2)
            c1s.append(c1)
            c2s.append(c2)

            print "CrossLinkMS: added pair %d %s %d %s" % (r1,c1,r2,c2)

            self.added_pairs_list.append((p1,p2,crosslinker))


        if len(p1s)>0:
            rs_name= '{:05}'.format(self.index % 100000)

            ln=impisd2.CrossLinkMSRestraint(self.sigmaglobal,self.crosslinker_dict[crosslinker])
            for i in range(len(p1s)):
                print rs_name,i
                ln.add_contribution(p1s[i],p2s[i])
                self.pairs.append((p1s[i], p2s[i], crosslinker, rs_name,
                     100, 100, (r1s[i],c1s[i],i), (r2s[i],c2s[i],i), crosslinker, i, ln))
            self.rs.add_restraint(ln)

#---------------------------------

    def simulate_data(self,crosslinker,weights,sensitivity_threshold=0.1,
                         false_positive_half=0.02,elapsed_time=0.01,
                         ratemin=100,ratemax=100):

        from random import choice
        from random import randrange
        from itertools import combinations
        from numpy.random import binomial

        self.weights=weights
        self.sensitivity_threshold=sensitivity_threshold
        self.false_positive_half=false_positive_half
        self.elapsed_time=elapsed_time
        self.ratemin=ratemin
        self.ratemax=ratemax

        #dictionary of random reaction rates
        #check if they were already initialized
        if self.reaction_rates==None:
            self.reaction_rates={}

            s0=IMP.atom.Selection(self.prots[0],
                                  residue_type=IMP.atom.get_residue_type("K"),
                                  atom_type=IMP.atom.AT_CA)
            self.residue_list=[]

            for p1 in s0.get_selected_particles():
                (r1,c1)=tools.get_residue_index_and_chain_from_particle(p1)
                self.residue_list.append((r1,c1))
                if self.ratemin!=self.ratemax:
                    self.reaction_rates[(r1,c1)]=randrange(self.ratemin,self.ratemax,1)
                else:
                    self.reaction_rates[(r1,c1)]=self.ratemax

        if self.allpairs_database==None:
            self.allpairs_database=[]

        #generate the restraints
        allcomb=list(combinations(self.residue_list,2))
        for ((r1,c1),(r2,c2)) in allcomb:

            p1s=[]
            p2s=[]

            for copy in self.prots:
                s1=IMP.atom.Selection(copy, residue_index=r1,
                                      chains=c1, atom_type=IMP.atom.AT_CA)
                p1s.append(s1.get_selected_particles()[0])
                s2=IMP.atom.Selection(copy, residue_index=r2,
                                      chains=c2, atom_type=IMP.atom.AT_CA)
                p2s.append(s2.get_selected_particles()[0])

            ln=impisd2.CrossLinkMSRestraint(self.sigmaglobal,self.crosslinker_dict[crosslinker])
            for i in range(len(p1s)):
                ln.add_contribution(p1s[i],p2s[i])
                d1=IMP.core.XYZ(p1s[i])
                d2=IMP.core.XYZ(p2s[i])
                dist=IMP.core.get_distance(d1,d2)
                reactionrate1=self.reaction_rates[(r1,c1)]
                reactionrate2=self.reaction_rates[(r2,c2)]
                prob=ln.get_marginal_probabilities()[i]
                effrate=float(reactionrate1*reactionrate2)/(reactionrate1+reactionrate2)
                probt=self.weights[i]*(1-exp(-effrate*prob*elapsed_time))
                falsepositiveprob=exp(-probt/false_positive_half)
                falsepositivebool=False
                falsepositive=binomial(n=1, p=falsepositiveprob)
                if (falsepositive==1):
                    falsepositivebool=True
                    randompair=choice(allcomb)
                    randpair1=randompair[0]
                    randpair2=randompair[1]
                else:
                    randpair1=(r1,c1)
                    randpair2=(r2,c2)
                if (probt>sensitivity_threshold):
                    detectedbool=True
                else:
                    detectedbool=False

                self.allpairs_database.append({})
                self.allpairs_database[-1]["Particle1"]=p1s[i]
                self.allpairs_database[-1]["Particle2"]=p2s[i]
                self.allpairs_database[-1]["Distance"]=dist
                self.allpairs_database[-1]["Crosslinker"]=crosslinker
                self.allpairs_database[-1]["IMPRestraint"]=ln
                self.allpairs_database[-1]["IMPRestraint_Probability"]=prob
                self.allpairs_database[-1]["Reaction_Rate1"]=reactionrate1
                self.allpairs_database[-1]["Reaction_Rate2"]=reactionrate2
                self.allpairs_database[-1]["Effective_Rate"]=effrate
                self.allpairs_database[-1]["CrossLink_Fraction"]=probt
                self.allpairs_database[-1]["Resid1_Chainid1_Copy1"]=(r1,c1,i)
                self.allpairs_database[-1]["Resid2_Chainid2_Copy2"]=(r2,c2,i)
                self.allpairs_database[-1]["Is_False_Positive"]=falsepositivebool
                self.allpairs_database[-1]["Identified_Pair1"]=randpair1
                self.allpairs_database[-1]["Identified_Pair2"]=randpair2
                self.allpairs_database[-1]["Is_Detected"]=detectedbool

    def set_hierarchy(self,prots):
        #we use it to change the hierarchy
        self.prots=prots

    def initialize_simulated_database(self):
        #we use it to restart the simulation
        self.allpairs_database=None

    def get_number_detected_inter(self,xl_type):
        #we use it to see ho many xls of a give type (eg. BS3) were detected as inter-chain
        ndetected=0
        for el in self.allpairs_database:
            if el["Is_Detected"]==True and \
                 ( el["Identified_Pair1"][1]!=el["Identified_Pair2"][1] ) and \
                 el["Crosslinker"]==xl_type:
                ndetected+=1
        return ndetected

    def get_number_detected_inter_false_positive(self,xl_type):
        #we use it to see ho many xls of a give type (eg. BS3) were detected as inter-chain
        ndetectedfp=0
        for el in self.allpairs_database:
            if el["Is_Detected"]==True and \
                 ( el["Identified_Pair1"][1]!=el["Identified_Pair2"][1] ) and \
                 el["Crosslinker"]==xl_type and el["Is_False_Positive"]==True:
                ndetectedfp+=1
        return ndetectedfp


    def show_simulated_data(self,what="Inter"):
        #"what" can be "All", "Detected", "FalsePositive", "TruePositive", "Intra", "Inter"
        if self.allpairs_database!=None:
            detectedlist=[]
            for el in self.allpairs_database:
                printbool=False
                if el["Is_Detected"]==True:
                    p1=el["Identified_Pair1"]
                    p2=el["Identified_Pair2"]
                    isf=el["Is_False_Positive"]
                    isinter=(el["Identified_Pair1"][1]!=el["Identified_Pair2"][1])
                    cl=el["Crosslinker"]
                    detectedlist.append((p1,p2,isf,cl,isinter))

                if el["Is_Detected"]==True and what=="Detected": printbool=True
                if el["Is_Detected"]==True and el["Is_False_Positive"]==True and what=="FalsePositive": printbool=True
                if el["Is_Detected"]==True and el["Is_False_Positive"]==False and what=="TruePositive": printbool=True
                if el["Is_Detected"]==True and what=="Intra" and \
                    ( el["Identified_Pair1"][1]==el["Identified_Pair2"][1] ) : printbool=True
                if el["Is_Detected"]==True and what=="Inter" and \
                    ( el["Identified_Pair1"][1]!=el["Identified_Pair2"][1] ) : printbool=True
                if what=="All": printbool=True

                if printbool:
                    print "Residue1: %6s, chainid1: %6s, copy1: %6d" % el["Resid1_Chainid1_Copy1"]
                    print "Residue2: %6s, chainid2: %6s, copy2: %6d" % el["Resid2_Chainid2_Copy2"]
                    keylist=el.keys()
                    keylist.sort()
                    for k in keylist:
                        print "----", k, el[k]
            print "r1 c1 r2 c2 FP XL Inter"
            for d in detectedlist:
                print d[0][0],d[0][1],d[1][0],d[1][1],d[2],d[3],d[4]
        else:
            print "CrossLinkMS: Simulated data not initialized"
            exit()

    def dump_simulated_data(self,filename="simulated_cross_link.dat"):
            #dump the whole simulated xl database on a file
        sclf=open(filename,"w")
        for el in self.allpairs_database:
            sclf.write(str(el))
            sclf.write("\n")
        sclf.close()

    def write_simulated_data(self,filename="detected_cross_link.dat"):
        #dump the whole simulated xl database on a file
        sclf=open(filename,"w")
        index=0
        for el in self.allpairs_database:
            if el["Is_Detected"]==True:
                index+=1
                p1=el["Identified_Pair1"]
                p2=el["Identified_Pair2"]
                isf=el["Is_False_Positive"]
                cl=el["Crosslinker"]
                s=" ".join([str(p1[0]), p1[1], str(p2[0]), p2[1], cl, str(index),str(isf), "T"])
                sclf.write(s)
                sclf.write("\n")
        sclf.close()

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)
        self.m.add_restraint(self.rs2)

    def get_hierarchies(self):
        return self.prots

    def get_particles(self):
        return self.sigmaglobal

    def get_restraint_sets(self):
        return self.rs,self.rs2

    def get_restraint(self):
        tmprs=IMP.RestraintSet('xlms')
        tmprs.add_restraint(self.rs)
        tmprs.add_restraint(self.rs2)
        return tmprs

    def set_output_level(self,level="low"):
        #this might be "low" or "high"
        self.outputlevel=level

    def print_chimera_pseudobonds(self,filesuffix,model=0):
        f=open(filesuffix+".chimera","w")
        atype="ca"
        for p in self.pairs:
            s="#"+str(model)+":"+str(p[6][0])+"."+p[6][1]+"@"+atype+" #"+str(model)+":"+str(p[7][0])+"."+p[7][1]+"@"+atype
            f.write(s+"\n")
        f.close()

    def get_particles_to_sample(self):
        ps={}
        if self.sigmaissampled:
            ps["Nuisances_CrossLinkMS_Sigma_"+self.label]=([self.sigmaglobal],self.sigmatrans)
        return ps

    def get_output(self):
        #content of the crosslink database pairs
        #self.pairs.append((p1s[i], p2s[i], crosslinker, rs_name, 100, 100, (r1,c1,i),  (r2,c2,i), crosslinker, i, ln))
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        score2=self.rs2.unprotected_evaluate(None)
        output["_TotalScore"]=str(score+score2)        
        
        output["CrossLinkMS_Likelihood_"+self.label]=str(score)
        output["CrossLinkMS_Prior_"+self.label]=str(score2)
        output["CrossLinkMS_Sigma"]=str(self.sigmaglobal.get_scale())


        if self.outputlevel=="high":
            for i in range(len(self.pairs)):

                p0=self.pairs[i][0]
                p1=self.pairs[i][1]
                crosslinker=self.pairs[i][2]
                ln=self.pairs[i][10]
                index=self.pairs[i][9]
                rsname=self.pairs[i][3]
                resid1=self.pairs[i][6][0]
                chain1=self.pairs[i][6][1]
                copy1=self.pairs[i][6][2]
                resid2=self.pairs[i][7][0]
                chain2=self.pairs[i][7][1]
                copy2=self.pairs[i][7][2]
                label_copy=str(rsname)+":"+str(index)+"-"+str(resid1)+":"+chain1+"_"+"-"+str(resid2)+":"+chain2+"_"+crosslinker
                output["CrossLinkMS_Partial_Probability_"+label_copy]=str(ln.get_marginal_probabilities()[index])

                if copy1==0:
                    label=str(resid1)+":"+chain1+"_"+str(resid2)+":"+chain2
                    #output["CrossLinkMS_Combined_Probability_"+str(rsname)+"_"+crosslinker+"_"+label]=str(ln.get_probability())
                    output["CrossLinkMS_Score_"+str(rsname)+"_"+crosslinker+"_"+label]=str(ln.unprotected_evaluate(None))


                d0=IMP.core.XYZ(p0)
                d1=IMP.core.XYZ(p1)
                output["CrossLinkMS_Distance_"+label_copy]=str(IMP.core.get_distance(d0,d1))

        return output

###############################################################

class BinomialXLMSRestraint():



    def __init__(self,m,prots,
                 listofxlinkertypes=["BS3","BS2G","EGS"],map_between_protein_names_and_chains=None,typeofprofile='pfes'):
        
        if map_between_protein_names_and_chains==None:
           map_between_protein_names_and_chains={}
        
        global impisd2, tools, exp
        import IMP.isd2 as impisd2
        import IMP.pmi.tools as tools
        
        self.setup=0
           
        self.label="None"
        self.rs = IMP.RestraintSet('xlms')
        self.rs2 = IMP.RestraintSet('jeffreys') 
        self.m=m       
        self.prots=prots
        self.pairs=[]
        
        self.weightmaxtrans=0.05
        self.weightissampled=False
       
        self.sigmainit=5.0
        self.sigmamin=1.0
        self.sigmaminnuis=0.0
        self.sigmamax=10.0
        self.sigmamaxnuis=11.0 
        self.nsigma=100
        self.sigmaissampled=True
        self.sigmamaxtrans=0.1
        
        self.betainit=1.0
        self.betamin=1.0
        self.betamax=4.0
        if self.setup==1:
           self.betaissampled=True
           print "BinomialXLMSRestraint: beta is sampled"
        if self.setup==0:
           self.betaissampled=False      
           print "BinomialXLMSRestraint: beta is NOT sampled"     
        self.betamaxtrans=0.01
        
        '''
        self.deltainit=0.001
        self.deltamin=0.001
        self.deltamax=0.1
        self.deltaissampled=False
        self.deltamaxtrans=0.001
        
        self.laminit=5.0
        self.lammin=0.01
        self.lamminnuis=0.00001
        self.lammax=10.0
        self.lammaxnuis=100.0
        self.lamissampled=False        
        self.lammaxtrans=0.1        
        '''
        
        self.epsilon=0.01
        self.psi_dictionary={}

        self.sigma=tools.SetupNuisance(self.m,self.sigmainit,
             self.sigmaminnuis,self.sigmamaxnuis,self.sigmaissampled).get_particle()

        self.beta=tools.SetupNuisance(self.m,self.betainit,
             self.betamin,self.betamax,self.betaissampled).get_particle()    
        
        '''
        self.delta=tools.SetupNuisance(self.m,self.deltainit,
             self.deltamin,self.deltamax,self.deltaissampled).get_particle()

        self.lam=tools.SetupNuisance(self.m,self.laminit,
             self.lamminnuis,self.lammaxnuis,self.lamissampled).get_particle()

        self.weight=tools.SetupWeight(m,False).get_particle()
        
        for n in range(len(self.prots)):
            self.weight.add_weight()
        '''
        self.outputlevel="low"
        self.listofxlinkertypes=listofxlinkertypes
        #simulated reaction rates
        self.reaction_rates=None
        self.allpairs_database=None
        self.residue_list=None

        self.crosslinker_dict=self.get_crosslinker_dict(typeofprofile)
        #this map is used in the new cross-link file reader
        #{"Nsp1":"A","Nup82":"B"} etc.
        self.mbpnc=map_between_protein_names_and_chains
        #check whther the file was initialized
    
    def create_psi(self,index,value):
        if value==None:
           self.psiinit=0.01
           self.psiissampled=True
           print "BinomialXLMSRestraint: psi "+str(index)+" is sampled"           
        else:
           self.psiinit=value 
           self.psiissampled=False
           print "BinomialXLMSRestraint: psi "+str(index)+" is NOT sampled"               
        self.psiminnuis=0.0000001
        self.psimaxnuis=0.4999999
        self.psimin=    0.01
        self.psimax=    0.49
        self.psitrans=  0.01 
        self.psi=tools.SetupNuisance(self.m,self.psiinit,
             self.psiminnuis,self.psimaxnuis,self.psiissampled).get_particle()
        self.psi_dictionary[index]=(self.psi,self.psitrans,self.psiissampled)    
        
    def get_psi(self,index,value):
        if not index in self.psi_dictionary:
           self.create_psi(index,value)
        return self.psi_dictionary[index]

    def get_crosslinker_dict(self,typeofprofile):
        #fill the cross-linker pmfs
        #to accelerate the init the list listofxlinkertypes might contain only yht needed crosslinks

        disttuple=(0.0, 200.0, 500)
        omegatuple=(0.01, 1000.0, 30)
        sigmatuple=(self.sigmamin, self.sigmamax, self.nsigma)

        crosslinker_dict={}
        if "BS3" in self.listofxlinkertypes:
            crosslinker_dict["BS3"]=tools.get_cross_link_data("bs3l",
               "pmf_bs3l_tip3p.txt.standard",disttuple,omegatuple,sigmatuple,don=None,doff=None,prior=1,type_of_profile=typeofprofile)
        if "BS2G" in self.listofxlinkertypes:
            crosslinker_dict["BS2G"]=tools.get_cross_link_data("bs2gl",
              "pmf_bs2gl_tip3p.txt.standard",disttuple,omegatuple,sigmatuple,don=None,doff=None,prior=1,type_of_profile=typeofprofile)
        if "EGS" in self.listofxlinkertypes:
            crosslinker_dict["EGS"]=tools.get_cross_link_data("egl",
              "pmf_egl_tip3p.txt.standard",disttuple,omegatuple,sigmatuple,don=None,doff=None,prior=1,type_of_profile=typeofprofile)
        if "Short" in self.listofxlinkertypes:
            #setup a "short" xl with an half length of 10 Ang
            crosslinker_dict["Short"]=tools.get_cross_link_data_from_length(10.0,disttuple,omegatuple,sigmatuple)
        return crosslinker_dict        

    def add_restraints(self,restraint_file=None,oldfile=False):

        #get the restraints from external file
        f = open(restraint_file)
        restraint_list = f.readlines()

        self.index=0
            
        self.added_pairs_list=[]
        self.missing_residues=[]
        for line in restraint_list:
            #force_restraint=True makes all intra rigid body restraint to be accepted
            force_restraint=False

            #read with the new file parser
            totallist=eval(line)
            self.add_crosslink_according_to_new_file(totallist)

        self.rs2.add_restraint(impisd2.UniformPrior(self.sigma,1000000000.0,self.sigmamax,self.sigmamin))
        #self.rs2.add_restraint(impisd2.JeffreysRestraint(self.sigma))
        
        for psiindex in self.psi_dictionary:
          if self.psi_dictionary[psiindex][2]:
            psip=self.psi_dictionary[psiindex][0]
            
            if self.setup==0:
               print "BinomialXLMSRestraint: setup 0, adding BinomialJeffreysPrior to psi particle "+str(psiindex) 
               self.rs2.add_restraint(impisd2.BinomialJeffreysPrior(psip))
               #self.rs2.add_restraint(impisd2.JeffreysRestraint(psip))
            
            self.rs2.add_restraint(impisd2.UniformPrior(psip,1000000000.0,self.psimax,self.psimin))
        
    def add_crosslink_according_to_new_file(self,totallist,constructor=0):
        #the constructor variable specify what constroctor to use in 
        #the restraint   0: use predetermined f.p.r. (psi) 
        #1: infer the f.p.r. defining a psi nuisance defined by a class
        force_restraint=False
        ambiguous_list=totallist[0]
        crosslinker=totallist[1]
        if (totallist[2]=="F"): force_restraint=True

        p1s=[]
        p2s=[]
        r1s=[]
        r2s=[]
        c1s=[]
        c2s=[]
        psis=[]
        psivalues=[]
        self.index+=1 
        for pair in ambiguous_list:
            error=False
           
            
            
            try:
               c1=self.mbpnc[pair[0][0]]
            except:
               print "BinomialXLMSRestraint: WARNING> protein name "+pair[0][0]+" was not defined"
               continue 
            try:                         
               c2=self.mbpnc[pair[1][0]]   
            except:
               print "BinomialXLMSRestraint: WARNING> protein name "+pair[1][0]+" was not defined"
               continue                         
            
            r1=int(pair[0][1])
            r2=int(pair[1][1])
            psi=float(pair[2])
            try:
               psivalue=float(pair[3])
            except:
               psivalue=None
            
            print '''CrossLinkMS: attempting to add restraint between
                     residue %d of chain %s and residue %d of chain %s''' % (r1,c1,r2,c2)


            try:
                s1=IMP.atom.Selection(self.prots[0], residue_index=r1, chains=c1, atom_type=IMP.atom.AT_CA)
                p1=(s1.get_selected_particles()[0])
            except:
                print "BinomialXLMSRestraint: WARNING> residue %d of chain %s is not there" % (r1,c1)
                error=True
                self.missing_residues.append((r1,c1))
            try:
                s2=IMP.atom.Selection(self.prots[0], residue_index=r2, chains=c2, atom_type=IMP.atom.AT_CA)
                p2=(s2.get_selected_particles()[0])
            except:
                print "BinomialXLMSRestraint: WARNING> residue %d of chain %s is not there" % (r2,c2)
                error=True
                self.missing_residues.append((r2,c2))            
            if error: continue
            

            s1=IMP.atom.Selection(self.prots[0], residue_index=r1, chains=c1, atom_type=IMP.atom.AT_CA)
            p1=s1.get_selected_particles()[0]
            s2=IMP.atom.Selection(self.prots[0], residue_index=r2, chains=c2, atom_type=IMP.atom.AT_CA)
            p2=s2.get_selected_particles()[0]
            #this list contains the list of symmetric pairs to avoid duplicates
            if (p1,p2,crosslinker) in self.added_pairs_list:
                print "BinomialXLMSRestraint: WARNING> pair %d %s %d %s already there" % (r1,c1,r2,c2)
                continue
            if (p2,p1,crosslinker) in self.added_pairs_list:
                print "BinomialXLMSRestraint: WARNING> pair %d %s %d %s already there" % (r1,c1,r2,c2)
                continue

            #check whether the atom pair belongs to the same rigid body
            if(IMP.core.RigidMember.particle_is_instance(p1) and
               IMP.core.RigidMember.particle_is_instance(p2) and
               IMP.core.RigidMember(p1).get_rigid_body() ==
               IMP.core.RigidMember(p2).get_rigid_body() and not force_restraint):
                print '''BinomialXLMSRestraint: WARNING> residue %d of chain %s and
                       residue %d of chain %s belong to the same rigid body''' % (r1,c1,r2,c2)
                continue

            p1s.append(p1)
            p2s.append(p2)
            r1s.append(r1)
            r2s.append(r2)
            c1s.append(c1)
            c2s.append(c2)
            psis.append(psi)
            psivalues.append(psivalue)
            
            print "BinomialXLMSRestraint: added pair %d %s %d %s" % (r1,c1,r2,c2)
                       
            self.added_pairs_list.append((p1,p2,crosslinker))


        if len(p1s)>0:
            rs_name= '{:05}'.format(self.index % 100000)
            
            if self.setup==0:  
               print "BinomialXLMSRestraint: constructor 0" 
               ln=impisd2.BinomialCrossLinkMSRestraint(self.m,self.sigma,self.epsilon,self.crosslinker_dict[crosslinker])

            if self.setup==1:  
               print "BinomialXLMSRestraint: constructor 1" 
               ln=impisd2.BinomialCrossLinkMSRestraint(self.m,self.sigma,self.beta,self.epsilon,self.crosslinker_dict[crosslinker])

            for i in range(len(p1s)):
                ln.add_contribution()

                psi=self.get_psi(psis[i],psivalues[i]) 
                                
                ln.add_particle_pair(i,(p1s[i].get_index(),p2s[i].get_index()),psi[0].get_particle().get_index())
                self.pairs.append((p1s[i], p2s[i], crosslinker, rs_name,
                     self.index, 100, (r1s[i],c1s[i],i), (r2s[i],c2s[i],i), crosslinker, i, ln))
                     
                h=IMP.core.Linear(0,0.03)
                dps=IMP.core.DistancePairScore(h)
                pr=IMP.core.PairRestraint(dps,IMP.ParticlePair(p1s[i],p2s[i]))  
                self.rs2.add_restraint(pr)
                
            self.rs.add_restraint(ln)

        print "BinomialXLMSRestraint: missing residues"
        for ms in self.missing_residues:
            print "BinomialXLMSRestraint:missing "+str(ms)

            
        #self.rs2.add_restraint(impisd2.IntensityThresholdRestraint(self.delta))
        #self.rs2.add_restraint(impisd2.UniformPrior(self.delta,1000000000.0,self.delta.get_upper(),self.delta.get_lower()))
        

        
        #exit()
            #self.rs2.add_restraint(impisd2.UniformPrior(psip,1000000000.0,0.5,0.01))
            
            
    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)
        self.m.add_restraint(self.rs2)

    def get_hierarchies(self):
        return self.prots

    def get_particles(self):
        return self.sigmaglobal

    def get_restraint_sets(self):
        return self.rs,self.rs2

    def get_restraint(self):
        tmprs=IMP.RestraintSet('xlms')
        tmprs.add_restraint(self.rs)
        tmprs.add_restraint(self.rs2)
        return tmprs

    def set_output_level(self,level="low"):
        #this might be "low" or "high"
        self.outputlevel=level

    def print_chimera_pseudobonds(self,filesuffix,model=0):
        f=open(filesuffix+".chimera","w")
        atype="ca"
        for p in self.pairs:
            s="#"+str(model)+":"+str(p[6][0])+"."+p[6][1]+"@"+atype+" #"+str(model)+":"+str(p[7][0])+"."+p[7][1]+"@"+atype
            f.write(s+"\n")
        f.close()

    def print_chimera_pseudobonds_with_psiindexes(self,filesuffix,model=0):
        
        f=open(filesuffix+".chimera","w")
        atype="ca"
        for p in self.pairs:
            s="#"+str(model)+":"+str(p[6][0])+"."+p[6][1]+"@"+atype+" #"+str(model)+":"+str(p[7][0])+"."+p[7][1]+"@"+atype
            f.write(s+"\n")
        f.close()

    def get_output(self):
        #content of the crosslink database pairs
        #self.pairs.append((p1s[i], p2s[i], crosslinker, rs_name, psi, 100, (r1,c1,i),  (r2,c2,i), crosslinker, i, ln))
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        score2=self.rs2.unprotected_evaluate(None)
        output["_TotalScore"]=str(score+score2)        
        output["CrossLinkMS_Likelihood_"+self.label]=str(score)
        output["CrossLinkMS_Prior_"+self.label]=str(score2)


        if self.outputlevel=="high":
            for i in range(len(self.pairs)):

                p0=self.pairs[i][0]
                p1=self.pairs[i][1]
                crosslinker=self.pairs[i][2]
                psiindex=self.pairs[i][4]
                ln=self.pairs[i][10]
                index=self.pairs[i][9]
                rsname=self.pairs[i][3]
                resid1=self.pairs[i][6][0]
                chain1=self.pairs[i][6][1]
                copy1=self.pairs[i][6][2]
                resid2=self.pairs[i][7][0]
                chain2=self.pairs[i][7][1]
                copy2=self.pairs[i][7][2]
                label_copy=str(rsname)+":"+str(index)+"-"+str(resid1)+":"+chain1+"_"+"-"+str(resid2)+":"+chain2+"_"+crosslinker

                if copy1==0:
                    label=str(resid1)+":"+chain1+"_"+str(resid2)+":"+chain2
                    #output["CrossLinkMS_Combined_Probability_"+str(rsname)+"_"+crosslinker+"_"+label]=str(ln.get_probability())
                    output["CrossLinkMS_Score_"+str(rsname)+"_"+str(index)+"_"+crosslinker+"_"+label]=str(ln.unprotected_evaluate(None))


                d0=IMP.core.XYZ(p0)
                d1=IMP.core.XYZ(p1)
                output["CrossLinkMS_Distance_"+str(index)+"_"+label_copy]=str(IMP.core.get_distance(d0,d1))
        
        output["CrossLinkMS_Sigma_"+self.label]=str(self.sigma.get_scale())
        '''
        output["CrossLinkMS_Delta_"+self.label]=str(self.delta.get_scale())
        output["CrossLinkMS_Lambda_"+self.label]=str(self.lam.get_scale())
        '''
        output["CrossLinkMS_Beta_"+self.label]=str(self.beta.get_scale())
        for psiindex in self.psi_dictionary:
            output["CrossLinkMS_Psi_"+str(psiindex)+"_"+self.label]=str(self.psi_dictionary[psiindex][0].get_scale())
        '''
        for n in range(self.weight.get_number_of_states()):
           output["CrossLinkMS_weights_"+str(n)+"_"+self.label]=str(self.weight.get_weight(n))
        '''
        return output

    def get_particles_to_sample(self):
        ps={}
        if self.sigmaissampled:
           ps["Nuisances_CrossLinkMS_Sigma_"+self.label]=([self.sigma],self.sigmamaxtrans)
        '''
        if self.deltaissampled:
           ps["Nuisances_CrossLinkMS_Delta_"+self.label]=([self.delta],self.deltamaxtrans)           
        if self.lamissampled:
           ps["Nuisances_CrossLinkMS_Lambda_"+self.label]=([self.lam],self.lammaxtrans)  
        '''
        if self.betaissampled:
           ps["Nuisances_CrossLinkMS_Beta_"+self.label]=([self.beta],self.betamaxtrans)
        
        for psiindex in self.psi_dictionary:
          if self.psi_dictionary[psiindex][2]:
            ps["Nuisances_CrossLinkMS_Psi_"+str(psiindex)+"_"+self.label]=([self.psi_dictionary[psiindex][0]],self.psi_dictionary[psiindex][1])
        
        '''
        if self.weightissampled:
           ps["Weights_CrossLinkMS_"+self.label]=([self.weight],self.weightmaxtrans)
        '''
        return ps    

##############################################################


class ConnectivityCrossLinkMS():
    '''
    this restraint allows ambiguous crosslinking between multiple copies
    it is a variant of the SimplifiedCrossLinkMS
    '''

    def __init__(self,prot,restraints_file,expdistance,strength,resolution=None):

        self.rs=IMP.RestraintSet('data')
        self.weight=1.0
        self.prot=prot
        self.label="None"
        self.pairs=[]
        self.m=self.prot.get_model()

        self.outputlevel="low"
        self.expdistance=expdistance
        self.strength=strength
        
        if resolution!=None:
          particles=[]
          for prot in self.prot.get_children():
             particles+=IMP.pmi.tools.get_particles_by_resolution(prot,resolution) 
        
        #fill the cross-linker pmfs
        #to accelerate the init the list listofxlinkertypes might contain only yht needed crosslinks


        for line in open(restraints_file):

            tokens=line.split()
            #skip character
            if (tokens[0]=="#"): continue
            r1=int(tokens[2])
            c1=tokens[0]
            r2=int(tokens[3])
            c2=tokens[1]

            
            hrc1=[]
            for h in self.prot.get_children():
                if c1 in h.get_name():
                   hrc1.append(h)
            
            #hrc1 = [h for h in self.prot.get_children() if c1 in h.get_name()][0]
            #print line
            
            s1=IMP.atom.Selection(hierarchies=hrc1, residue_index=r1)
            ps1=s1.get_selected_particles()
            

            
            if len(ps1)==0:
                print "ConnectivityCrossLinkMS: WARNING> residue %d of chain %s is not there" % (r1,c1)
                continue
            #hrc2 = [h for h in self.prot.get_children() if c2 in h.get_name()][0]
            
            hrc2=[]   
            for h in self.prot.get_children():
                if c2 in h.get_name():
                   hrc2.append(h)
            
            
            s2=IMP.atom.Selection(hierarchies=hrc2, residue_index=r2)
            ps2=s2.get_selected_particles()
            if len(ps2)==0:
                print "ConnectivityCrossLinkMS: WARNING> residue %d of chain %s is not there" % (r2,c2)
                continue

            if resolution!=None:
             
              #get the intersection to remove redundant particles
              ps1=(list(set(ps1) & set(particles)))
              ps2=(list(set(ps2) & set(particles)))
              s1=IMP.atom.Selection(ps1)
              s2=IMP.atom.Selection(ps2)
            
            
            sels=[s1,s2]
            cr = IMP.atom.create_connectivity_restraint(sels, self.expdistance,self.strength)
    
            self.rs.add_restraint(cr)
            self.pairs.append((ps1,hrc1,c1,r1,ps2,hrc2,c2,r2,cr))
            

    def set_label(self,label):
        self.label=label
        self.rs.set_name(label)
        for r in self.rs.get_restraints():
            r.set_name(label)

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_hierarchies(self):
        return self.prot

    def get_restraint_sets(self):
        return self.rs
        
    def get_restraint(self):
        return self.rs        

    def set_output_level(self,level="low"):
            #this might be "low" or "high"
        self.outputlevel=level

    def set_weight(self,weight):
        self.weight=weight
        self.rs.set_weight(weight)

    def get_output(self):
        #content of the crosslink database pairs
        #self.pairs.append((p1,p2,dr,r1,c1,r2,c2))
        self.m.update()

        output={}
        score=self.weight*self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)            
        output["SimplifiedCrossLinkMS_Score_"+self.label]=str(score)
        for n,p in enumerate(self.pairs):

            ps1=p[0]
            hrc1=p[1]
            c1=p[2]
            r1=p[3]
            ps2=p[4]
            hrc2=p[5]
            c2=p[6]
            r2=p[7]
            cr=p[8]
            for n1,p1 in enumerate(ps1):        
                name1=hrc1[n1].get_name()

                for n2,p2 in enumerate(ps2):
                  name2=hrc2[n2].get_name()                    
                  d1=IMP.core.XYZR(p1) 
                  d2=IMP.core.XYZR(p2)      
                  label=str(r1)+":"+name1+"_"+str(r2)+":"+name2    
                  output["ConnectivityCrossLinkMS_Distance_"+label]=str(IMP.core.get_distance(d1,d2))

            label=str(r1)+":"+c1+"_"+str(r2)+":"+c2
            output["ConnectivityCrossLinkMS_Score_"+label]=str(self.weight*cr.unprotected_evaluate(None))

        return output


###############################################################


class SimplifiedCrossLinkMS():

    def __init__(self,prot,restraints_file,expdistance,strength,resolution=None, columnmapping=None, truncatedharmonic=False):
        #columnindexes is a list of column indexes for protein1, protein2, residue1, residue2
        #by default column 0 = protein1; column 1 = protein2; column 2 = residue1; column 3 = residue2
        
        if columnmapping==None:
           columnmapping={}
           columnmapping["Protein1"]=0
           columnmapping["Protein2"]=1
           columnmapping["Residue1"]=2
           columnmapping["Residue2"]=3
        
        self.rs=IMP.RestraintSet('data')
        self.weight=1.0
        self.prot=prot
        self.label="None"
        self.pairs=[]
        self.m=self.prot.get_model()

        self.outputlevel="low"
        self.expdistance=expdistance
        self.strength=strength


        #fill the cross-linker pmfs
        #to accelerate the init the list listofxlinkertypes might contain only yht needed crosslinks
        protein1=columnmapping["Protein1"]
        protein2=columnmapping["Protein2"]
        residue1=columnmapping["Residue1"]
        residue2=columnmapping["Residue2"]        
        
        if resolution!=None:
          particles=[]
          for prot in self.prot.get_children():
             particles+=IMP.pmi.tools.get_particles_by_resolution(prot,resolution) 
        
        
        
        for line in open(restraints_file):

            tokens=line.split()
            #skip character
            if (tokens[0]=="#"): continue
            r1=int(tokens[residue1])
            c1=tokens[protein1]
            r2=int(tokens[residue2])
            c2=tokens[protein2]


            hrc1 = [h for h in self.prot.get_children() if h.get_name()==c1][0]
            s1=IMP.atom.Selection(hrc1, residue_index=r1)
            ps1=s1.get_selected_particles()

            hrc2 = [h for h in self.prot.get_children() if h.get_name()==c2][0]
            s2=IMP.atom.Selection(hrc2, residue_index=r2)
            ps2=s2.get_selected_particles()
       
            if resolution!=None: 
              #get the intersection to remove redundant particles
              ps1=(list(set(ps1) & set(particles)))
              ps2=(list(set(ps2) & set(particles)))
              
            if len(ps1)>1:
               print "SimplifiedCrossLinkMS: ERROR> residue %d of chain %s selects multiple particles"  % (r1,c1)
               exit()
            elif len(ps1)==0:
               print "SimplifiedCrossLinkMS: WARNING> residue %d of chain %s is not there" % (r1,c1)
               continue               

            if len(ps2)>1:
               print "SimplifiedCrossLinkMS: ERROR> residue %d of chain %s selects multiple particles"  % (r2,c2)
               exit()
            elif len(ps2)==0:
               print "SimplifiedCrossLinkMS: WARNING> residue %d of chain %s is not there" % (r2,c2)
               continue 

            
            p1=ps1[0]
            p2=ps2[0]

            if truncatedharmonic:
                #use truncated harmonic to account for outliers
                hub=IMP.core.TruncatedHarmonicBound(12.0,1.0/25.0,15.0,5)
            else: 
                hub= IMP.core.HarmonicUpperBound(self.expdistance,self.strength)
            df= IMP.core.SphereDistancePairScore(hub)
            dr= IMP.core.PairRestraint(df, (p1, p2))
            self.rs.add_restraint(dr)
            self.pairs.append((p1,p2,dr,r1,c1,r2,c2))


    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_hierarchies(self):
        return self.prot

    def get_restraint_sets(self):
        return self.rs
        
    def get_restraint(self):
        return self.rs        

    def get_restraints(self):
        rlist=[]
        for r in self.rs.get_restraints():
            rlist.append(IMP.core.PairRestraint.get_from(r))
        return rlist
    
    def get_particle_pairs(self):
        ppairs=[]
        for i in range(len(self.pairs)):        
            p0=self.pairs[i][0]
            p1=self.pairs[i][1]
            ppairs.append((p0,p1))
        return ppairs            
            
    def set_output_level(self,level="low"):
            #this might be "low" or "high"
        self.outputlevel=level

    def set_weight(self,weight):
        self.weight=weight
        self.rs.set_weight(weight)

    def get_output(self):
        #content of the crosslink database pairs
        #self.pairs.append((p1,p2,dr,r1,c1,r2,c2))
        self.m.update()

        output={}
        score=self.weight*self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)            
        output["SimplifiedCrossLinkMS_Score_"+self.label]=str(score)
        for i in range(len(self.pairs)):

            p0=self.pairs[i][0]
            p1=self.pairs[i][1]
            crosslinker='standard'
            ln=self.pairs[i][2]
            resid1=self.pairs[i][3]
            chain1=self.pairs[i][4]
            resid2=self.pairs[i][5]
            chain2=self.pairs[i][6]

            label=str(resid1)+":"+chain1+"_"+str(resid2)+":"+chain2
            output["SimplifiedCrossLinkMS_Score_"+crosslinker+"_"+label]=str(self.weight*ln.unprotected_evaluate(None))

            d0=IMP.core.XYZ(p0)
            d1=IMP.core.XYZ(p1)
            output["SimplifiedCrossLinkMS_Distance_"+label]=str(IMP.core.get_distance(d0,d1))

        return output

##############################################################

class SimplifiedPEMAP():

    def __init__(self,prot,restraints_file,expdistance,strength):

        self.rs=IMP.RestraintSet('data')
        self.prot=prot
        self.label="None"
        self.pairs=[]
        self.m=self.prot.get_model()

        self.outputlevel="low"
        self.expdistance=expdistance
        self.strength=strength

        #fill the cross-linker pmfs
        #to accelerate the init the list listofxlinkertypes might contain only yht needed crosslinks


        for line in open(restraints_file):

            tokens=line.split()
            #skip character
            if (tokens[0]=="#"): continue
            r1=int(tokens[2])
            c1=tokens[0]
            r2=int(tokens[3])
            c2=tokens[1]
            pcc=float(tokens[4])

            try:
                #hrc1 = [h for h in self.prot.get_children() if h.get_name()==c1][0]
                #s1=IMP.atom.Selection(hrc1, residue_index=r1)
                s1=IMP.atom.Selection(self.prot, residue_index=r1,chains=c1)
                p1=s1.get_selected_particles()[0]
                #print "SimplifiedPEMAP: HEREHERE> residue %d of chain %s is not there" % (r1,c1)
            except:
                print "SimplifiedPEMAP: WARNING> residue %d of chain %s is not there (w/ %d %s)" % (r1,c1,r2,c2)
                continue
            try:
                #hrc2 = [h for h in self.prot.get_children() if h.get_name()==c2][0]
                #s2=IMP.atom.Selection(hrc2, residue_index=r2)
                s2=IMP.atom.Selection(self.prot, residue_index=r2,chains=c2)
                p2=s2.get_selected_particles()[0]
                #print "SimplifiedPEMAP: HEREHERE> residue %d of chain %s is not there" % (r2,c2)
            except:
                print "SimplifiedPEMAP: WARNING> residue %d of chain %s is not there (w/ %d %s)" % (r2,c2,r1,c1)
                continue

            #This is harmonic potential for the pE-MAP data
            upperdist = self.get_upper_bond(pcc)
            limit=self.strength*(upperdist+15)**2+10.0
            hub= IMP.core.TruncatedHarmonicUpperBound(upperdist,self.strength,upperdist+15,limit)

            #This is harmonic for the X-link
            #hub= IMP.core.TruncatedHarmonicBound(17.0,self.strength,upperdist+15,limit)

            df= IMP.core.SphereDistancePairScore(hub)
            dr= IMP.core.PairRestraint(df, (p1, p2))
            self.rs.add_restraint(dr)
            self.pairs.append((p1,p2,dr,r1,c1,r2,c2))

            #Lower-bound restraint
            lowerdist = self.get_lower_bond(pcc)
            limit=self.strength*(lowerdist-15)**2+10.0
            hub2= IMP.core.TruncatedHarmonicLowerBound(lowerdist,self.strength,lowerdist+15,limit)

            #This is harmonic for the X-link
            #hub2= IMP.core.TruncatedHarmonicBound(17.0,self.strength,upperdist+15,limit)

            df2= IMP.core.SphereDistancePairScore(hub2)
            dr2= IMP.core.PairRestraint(df2, (p1, p2))
            self.rs.add_restraint(dr2)
            self.pairs.append((p1,p2,dr2,r1,c1,r2,c2))


    def get_upper_bond(self,pearsoncc):
        #return (pearsoncc-1.)/-0.0075
        return (pearsoncc-.5)/(-0.005415)

    def get_lower_bond(self,pearsoncc):
        return (pearsoncc-1.)/-0.0551

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_hierarchies(self):
        return self.prot

    def get_restraint_sets(self):
        return self.rs

    def set_output_level(self,level="low"):
            #this might be "low" or "high"
        self.outputlevel=level

    def get_output(self):
        #content of the crosslink database pairs
        #self.pairs.append((p1,p2,dr,r1,c1,r2,c2))
        self.m.update()

        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)            
        output["SimplifiedPEMAP_Score_"+self.label]=str(score)
        for i in range(len(self.pairs)):

            p0=self.pairs[i][0]
            p1=self.pairs[i][1]
            crosslinker='standard'
            ln=self.pairs[i][2]
            resid1=self.pairs[i][3]
            chain1=self.pairs[i][4]
            resid2=self.pairs[i][5]
            chain2=self.pairs[i][6]

            label=str(resid1)+":"+chain1+"_"+str(resid2)+":"+chain2
            output["SimplifiedPEMAP_Score_"+crosslinker+"_"+label]=str(ln.unprotected_evaluate(None))

            d0=IMP.core.XYZ(p0)
            d1=IMP.core.XYZ(p1)
            output["SimplifiedPEMAP_Distance_"+label]=str(IMP.core.get_distance(d0,d1))

        return output



##############################################################

class CrossLinkMSSimple():

    def __init__(self,prot,restraints_file,TruncatedHarmonic=True):
        """read crosslink restraints between two residue
        of different chains from an external text file
        syntax: part_name_1 part_name_2 distance error
        example:     0 1 1.0 0.1"""
        self.prot=prot
        self.rs=IMP.RestraintSet('xlms')
        self.pairs=[]
        self.m=self.prot.get_model()

        #crosslinker="BS3"
        #this is an harmonic potential with mean 12 and sigma=25, therefore
        #k=1/sigma^2

        if TruncatedHarmonic:
            #use truncated harmonic to account for outliers
            hf=IMP.core.TruncatedHarmonicBound(12.0,1.0/25.0,15.0,5)
        else:
            hf=IMP.core.Harmonic(12.0,1.0/25.0)
        dps=IMP.core.DistancePairScore(hf)
        
        #small linear contribution for long range
        h=IMP.core.Linear(0,0.03)
        dps2=IMP.core.DistancePairScore(h)

        index=0

        addedd_pairs_list=[]
        for line in open(restraints_file):

        #force_restraint=True makes all intra rigid body restraint to be accepted
            force_restraint=True
            tokens=line.split()

            #skip character
            if (tokens[0]=="#"): continue
            r1=int(tokens[0])
            c1=tokens[1]
            r2=int(tokens[2])
            c2=tokens[3]
            crosslinker=tokens[4]

            #two restraints with the same index will be ambiguous
            index=int(tokens[5])

            #force restraint even if it belong to the same rigid body, use it for ambiguous restraints
            if (tokens[len(tokens)-1]=="F"): force_restraint=True

            print "attempting to add restraint between residue %d of chain %s and residue %d of chain %s" % (r1,c1,r2,c2)

            p1s=[]
            p2s=[]

            #apply the cross-link to the main copy
            try:
                s1=IMP.atom.Selection(self.prot, residue_index=r1, chains=c1, atom_type=IMP.atom.AT_CA)
                p1=(s1.get_selected_particles()[0])
            except:
                print "WARNING> residue %d of chain %s is not there" % (r1,c1)
                continue
            try:
                s2=IMP.atom.Selection(self.prot, residue_index=r2, chains=c2, atom_type=IMP.atom.AT_CA)
                p2=(s2.get_selected_particles()[0])
            except:
                print "WARNING> residue %d of chain %s is not there" % (r2,c2)
                continue

            #part1=[]
            #part2=[]
            #this list contain the list of simmetric pairs to avoid duplications

            #this is the map between particle pairs and the restraints (there
            #might be more than one restraint per particle pair if you have
            #more than one cross-link type

            print "attempting to add restraint between residue %d of chain %s and residue %d of chain %s" % (r1,c1,r2,c2)

            #check whether the atom pair belongs to the same rigid body
            if(IMP.core.RigidMember.particle_is_instance(p1) and IMP.core.RigidMember.particle_is_instance(p2) and
                   IMP.core.RigidMember(p1).get_rigid_body() == IMP.core.RigidMember(p2).get_rigid_body() and not force_restraint):
                print "WARNING> residue %d of chain %s and residue %d of chain %s belong to the same rigid body" % (r1,c1,r2,c2)
                continue

                #this list contain the list of simmetric pairs to avoid duplications
            if (p1,p2,crosslinker) in addedd_pairs_list:
                print "WARNING> pair %d %s %d %s already there" % (r1,c1,r2,c2)
                continue
            if (p2,p1,crosslinker) in addedd_pairs_list:
                print "WARNING> pair %d %s %d %s already there" % (r1,c1,r2,c2)
                continue

            print "added pair %d %s %d %s" % (r1,c1,r2,c2)
            index+=1
            addedd_pairs_list.append((p1,p2,crosslinker))

            rs_name='restraint_'+str(index)

            ln=IMP.core.PairRestraint(dps,IMP.ParticlePair(p1,p2))
            ln.set_name("CrossLinkMSSimple_"+str(r1)+":"+str(c1)+"-"+str(r2)+":"+str(c2))
            ln.set_weight(1.0)

            self.rs.add_restraint(ln)
            

            pr=IMP.core.PairRestraint(dps2,IMP.ParticlePair(p1,p2))  

            self.rs.add_restraint(pr)


            self.pairs.append((p1,  p2,  crosslinker,  rs_name,  100,  100,  (r1,c1),  (r2,c2), crosslinker, ln,pr))


    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_restraint(self):
        return self.rs

    def print_chimera_pseudobonds(self,filesuffix,model=0):
        f=open(filesuffix+".chimera","w")
        atype="ca"
        for p in self.pairs:
            s="#"+str(model)+":"+str(p[6][0])+"."+p[6][1]+"@"+atype+" #"+str(model)+":"+str(p[7][0])+"."+p[7][1]+"@"+atype
            f.write(s+"\n")
        f.close()


    def get_output(self):
        #content of the crosslink database pairs
        #self.pairs.append((p1s[i], p2s[i], crosslinker, rs_name, 100, 100, (r1,c1),  (r2,c2), crosslinker, ln,pr))
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)            
        output["CrossLinkMSSimple_Score_"+self.label]=str(score)
        for i in range(len(self.pairs)):

            p0=self.pairs[i][0]
            p1=self.pairs[i][1]
            crosslinker=self.pairs[i][2]
            ln=self.pairs[i][9]
            pr=self.pairs[i][10]
            resid1=self.pairs[i][6][0]
            chain1=self.pairs[i][6][1]
            resid2=self.pairs[i][7][0]
            chain2=self.pairs[i][7][1]

            label=str(resid1)+":"+chain1+"_"+str(resid2)+":"+chain2
            output["CrossLinkMSSimple_Score_"+crosslinker+"_"+label]=str(ln.unprotected_evaluate(None))
            output["CrossLinkMSSimple_Score_Linear_"+crosslinker+"_"+label]=str(pr.unprotected_evaluate(None))
            d0=IMP.core.XYZ(p0)
            d1=IMP.core.XYZ(p1)
            output["CrossLinkMSSimple_Distance_"+label]=str(IMP.core.get_distance(d0,d1))

        return output


class SecondaryStructure():
    

    def __init__(self,prot,resrangetuple,ssstring,mixture=False,nativeness=1.0,kt_caff=0.1):
        #check that the secondary structure string
        #is compatible with the ssstring
        global impisd2
        import IMP.isd2 as impisd2
        
        self.prot=prot
        self.m=self.prot.get_model()
        self.dihe_dict={}
        self.ang_dict={}
        self.do_mix={}
        self.anglfilename=impisd2.get_data_path("CAAngleRestraint.dat")
        self.dihefilename=impisd2.get_data_path("CADihedralRestraint.dat")
        self.nativeness=nativeness
        self.kt_caff=kt_caff
        self.anglrs=IMP.RestraintSet("Angles")
        self.dihers=IMP.RestraintSet("Dihedrals")
        self.bondrs=IMP.RestraintSet("Bonds")
        self.label="None"
        if resrangetuple[1]-resrangetuple[0]+1+4!=len(ssstring):
            print "SecondaryStructure: residue range and SS string incompatible"
            exit()
        for i in range(resrangetuple[0]-2,resrangetuple[1]+3):
            self.dihe_dict[(i,resrangetuple[2])]=ssstring[i-resrangetuple[0]]
            self.ang_dict[(i,resrangetuple[2])]=ssstring[i-resrangetuple[0]]
            self.do_mix[(i,resrangetuple[2])]=mixture

        (bondrslist,anglrslist,diherslist,pairslist)=self.get_CA_force_field(resrangetuple)
        self.pairslist=pairslist
        self.anglrs.add_restraints(anglrslist)
        self.dihers.add_restraints(diherslist)
        self.bondrs.add_restraints(bondrslist)

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.anglrs)
        self.m.add_restraint(self.dihers)
        self.m.add_restraint(self.bondrs)

    def get_CA_force_field(self,resrange):
        bondrslist=[]
        anglrslist=[]
        diherslist=[]
        pairslist=[]
        # add bonds
        for res in range(resrange[0]-1,resrange[1]+1):
            ps=[]
            for delta in range(0,2):
                try:
                    s=IMP.atom.Selection(self.prot, chains=resrange[2], residue_index=res+delta, atom_type=IMP.atom.AT_CA)
                    ps.append(s.get_selected_particles()[0])
                except:
                    print "SecondaryStructure: Bond: didn't find "+resrange[2]+str(res+delta)
                    continue

            if (len(ps)!=2): continue
            pairslist.append(IMP.ParticlePair(ps[0],ps[1]))
            pairslist.append(IMP.ParticlePair(ps[1],ps[0]))
            br=self.get_distance_restraint(ps[0],ps[1],3.78,416.0)
            br.set_name('Bond_restraint')
            bondrslist.append(br)
        # add dihedrals
        for res in range(resrange[0]-2,resrange[1]+3):

            #if res not in dihe_dict: continue
            # get the appropriate parameters
            # get the particles
            ps=[]
            ssstring=""
            for delta in range(-2,+3):
                try:
                    s=IMP.atom.Selection(self.prot, chains=resrange[2],residue_index=res+delta, atom_type=IMP.atom.AT_CA)
                    ps.append(s.get_selected_particles()[0])
                    ssstring+=self.dihe_dict[(res+delta,resrange[2])]
                except:
                    print "SecondaryStructure: Dihedral: didn't find "+resrange[2]+str(res+delta)
                    continue

            if (len(ps)!=5): continue
            [phi0,phi1,score_dih]=self.read_potential_dihedral(ssstring,self.do_mix[(res,resrange[2])])
            pairslist.append(IMP.ParticlePair(ps[0],ps[3]))
            pairslist.append(IMP.ParticlePair(ps[3],ps[0]))
            pairslist.append(IMP.ParticlePair(ps[1],ps[4]))
            pairslist.append(IMP.ParticlePair(ps[4],ps[1]))
            dr=impisd2.CADihedralRestraint(ps[0],ps[1],ps[2],ps[3],ps[4],phi0,phi1,score_dih)
            dr.set_name('Dihedral restraint')
            diherslist.append(dr)
        # add angles
        for res in range(resrange[0]-1,resrange[1]+2):
            ps=[]
            ssstring=""
            for delta in range(-1,+2):
                try:
                    s=IMP.atom.Selection(self.prot, chains=resrange[2],residue_index=res+delta, atom_type=IMP.atom.AT_CA)
                    ps.append(s.get_selected_particles()[0])
                    ssstring+=self.ang_dict[(res+delta,resrange[2])]
                except:
                    print "SecondaryStructure: Angle: didn't find "+resrange[2]+str(res+delta)
                    continue

            if (len(ps)!=3): continue
            [psi,score_ang]=self.read_potential_angle(ssstring,self.do_mix[(res,resrange[2])])
            pairslist.append(IMP.ParticlePair(ps[0],ps[2]))
            pairslist.append(IMP.ParticlePair(ps[2],ps[0]))
            dr=impisd2.CAAngleRestraint(ps[0],ps[1],ps[2],psi,score_ang)
            dr.set_name('Angle restraint')
            anglrslist.append(dr)
        return (bondrslist,anglrslist,diherslist,pairslist)

    def read_potential_dihedral(self,string,mix=False):
    # read potentials for dihedral
        score_dih=[]
        phi0=[]; phi1=[]
        for i in range(0,36):
            phi0.append(i*10.0/180.0*math.pi)
            phi1.append(i*10.0/180.0*math.pi)
            for j in range(0,36):
                score_dih.append(0.0)
        # open file
        if not mix:
            f = open(self.dihefilename, 'r')
            for line in f.readlines():
                riga=(line.strip()).split()
                if (len(riga)==4 and riga[0]==string):
                    ii=int(float(riga[1])/10.0)
                    jj=int(float(riga[2])/10.0)
                    score_dih[ii*len(phi0)+jj]=-self.kt_caff*math.log(float(riga[3]))
            f.close()
        if mix:
            #mix random coil and native secondary structure
            counts=[]
            for i in range(0,36):
                for j in range(0,36):
                    counts.append(1.0)
            f = open(self.dihefilename, 'r')
            for line in f.readlines():
                riga=(line.strip()).split()
                if (len(riga)==4 and riga[0]==string):
                    ii=int(float(riga[1])/10.0)
                    jj=int(float(riga[2])/10.0)
                    counts[ii*len(phi0)+jj]+=self.nativeness*float(riga[3])
                if (len(riga)==4 and riga[0]=="-----"):
                    ii=int(float(riga[1])/10.0)
                    jj=int(float(riga[2])/10.0)
                    counts[ii*len(phi0)+jj]+=(1.0-self.nativeness)*float(riga[3])
            f.close()
            for i in range(len(counts)):
                score_dih[i]=-self.kt_caff*math.log(counts[i])
        return [phi0,phi1,score_dih]

    def read_potential_angle(self,string,mix=False):
    # read potentials for angles
        score_ang=[]
        psi=[]
        for i in range(0,180):
            psi.append(i/180.0*math.pi)
            score_ang.append(0.0)
        # read file
        if not mix:
            f = open(self.anglfilename, 'r')
            for line in f.readlines():
                riga=(line.strip()).split()
                if (len(riga)==3 and riga[0]==string):
                    ii=int(riga[1])
                    score_ang[ii]=-self.kt_caff*math.log(float(riga[2]))
            f.close()
        if mix:
            #mix random coil and native secondary structure
            counts=[]
            for i in range(0,180):
                counts.append(1.0)

            f = open(self.anglfilename, 'r')
            for line in f.readlines():
                riga=(line.strip()).split()
                if (len(riga)==3 and riga[0]==string):
                    ii=int(riga[1])
                    counts[ii]+=self.nativeness*float(riga[2])
                if (len(riga)==3 and riga[0]=="---"):
                    ii=int(riga[1])
                    counts[ii]+=(1.0-self.nativeness)*float(riga[2])
            f.close()
            for i in range(0,180):
                score_ang[i]=-self.kt_caff*math.log(counts[i])
        return [psi,score_ang]

    def get_excluded_pairs(self):
        return self.pairslist

    def get_restraint(self):
        tmprs=IMP.RestraintSet('tmp')
        tmprs.add_restraint(self.anglrs)
        tmprs.add_restraint(self.dihers)
        tmprs.add_restraint(self.bondrs)
        return tmprs

    def get_distance_restraint(self,p0, p1, d0, kappa):
        h=IMP.core.Harmonic(d0,kappa)
        dps=IMP.core.DistancePairScore(h)
        pr=IMP.core.PairRestraint(dps,IMP.ParticlePair(p0,p1))
        return pr

    def get_output(self):
        output={}
        self.m.update()
        score_angle=self.anglrs.unprotected_evaluate(None)
        score_dihers=self.dihers.unprotected_evaluate(None)
        score_bondrs=self.bondrs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score_angle+score_dihers+score_bondrs)    
        
        output["SecondaryStructure_Angles_"+self.label]=str(score_angle)
        output["SecondaryStructure_Dihedrals_"+self.label]=str(score_dihers)
        output["SecondaryStructure_Bonds_"+self.label]=str(score_bondrs)
        return output

###########################################################################

class WeightRestraint():
    def __init__(self,weight,lower,upper,kappa):
        global impisd2
        import IMP.isd2 as impisd2

        self.weight=weight
        self.m=self.weight.get_model()
        self.label="None"
        self.rs = IMP.RestraintSet('weight_restraint')
        self.lower =lower
        self.upper=upper
        self.kappa=kappa
        self.rs.add_restraint(impisd2.WeightRestraint(self.weight,self.lower,self.upper,self.kappa))

    def get_restraint(self,label):
        return self.rs

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def set_label(self,label):
        self.label=label

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)         
        output["WeightRestraint_"+self.label]=str(score)
        return output


class JeffreysPrior():
    def __init__(self,nuisance):
        global impisd2
        import IMP.isd2 as impisd2

        self.m=nuisance.get_model()
        self.label="None"
        self.rs = IMP.RestraintSet('jeffrey_prior')
        jp=impisd2.JeffreysRestraint(nuisance)
        self.rs.add_restraint(jp)

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def set_label(self,label):
        self.label=label

    def get_output(self):
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)         
        output["JeffreyPrior_"+self.label]=str(score)
        return output

###########################################################################

class SAXSISDRestraint():



    def __init__(self,prot,profile,weight=1):
        global impsaxs, impisd2, tools
        import IMP.saxs as impsaxs
        import IMP.isd2 as impisd2
        import IMP.pmi.tools as tools

        self.prot=prot
        self.m=self.prot.get_model()
        self.label="None"
        self.rs = IMP.RestraintSet('saxs')

        self.sigmamaxtrans=0.05
        self.gammamaxtrans=0.05
        self.prof = impsaxs.Profile(profile)

        atoms = IMP.atom.get_by_type(self.prot, IMP.atom.ATOM_TYPE)

        self.sigma = tools.SetupNuisance(self.m,10.0,0.00001,100,True).get_particle()

        self.th = impsaxs.Profile(self.prof.get_min_q(),
                self.prof.get_max_q(), self.prof.get_delta_q())

        self.th.calculate_profile(atoms)
        gammahat = array([self.prof.get_intensity(i)/self.th.get_intensity(i)
                            for i in xrange(self.prof.size()-1) ]).mean()

        self.gamma = tools.SetupNuisance(self.m,gammahat,0.01,20,True).get_particle()

        self.cov = eye(self.prof.size()).tolist()

        self.saxs = impisd2.SAXSRestraint(atoms, self.prof, self.sigma,
                                        self.gamma, self.cov, impsaxs.CA_ATOMS)

        self.rs.add_restraint(self.saxs)
        self.rs.set_weight(weight)

        #self.saxs_stuff={'nuis':(sigma,gamma),'cov':cov,
        #        'exp':prof,'th':tmp}

        self.rs2 = IMP.RestraintSet('jeffreys')
        j1 = impisd2.JeffreysRestraint(self.sigma)
        self.rs2.add_restraint(j1)
        j2 = impisd2.JeffreysRestraint(self.gamma)
        self.rs2.add_restraint(j2)

    def gen_covariance_matrices(self, tau, nsigma):
        v = tools.Variance(self.m, tau, nsigma, self.prot, self.th)
        v.run()
        self.cov = v.get_cov(relative=True)

    def get_gamma_value(self):
        return self.gamma.get_scale()

    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)
        self.m.add_restraint(self.rs2)

    def get_restraint(self):
        tmprs=IMP.RestraintSet('tmp')
        tmprs.add_restraint(self.rs)
        tmprs.add_restraint(self.rs2)
        return tmprs

    def set_gammamaxtrans(self,gammamaxtrans):
        self.gammamaxtrans=gammamaxtrans

    def set_sigmamaxtrans(self,sigmamaxtrans):
        self.sigmamaxtrans=sigmamaxtrans

    def get_particles_to_sample(self):
        ps={}
        ps["Nuisances_SAXSISDRestraint_Sigma_"+self.label]=([self.sigma],self.sigmamaxtrans)
        ps["Nuisances_SAXSISDRestraint_Gamma_"+self.label]=([self.gamma],self.gammamaxtrans)
        return ps

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        score2=self.rs2.unprotected_evaluate(None)
        output["_TotalScore"]=str(score+score2) 
        
        output["SAXSISDRestraint_Likelihood_"+self.label]=str(score)
        output["SAXSISDRestraint_Prior_"+self.label]=str(score2)
        output["SAXSISDRestraint_Sigma_"+self.label]=str(self.sigma.get_scale())
        output["SAXSISDRestraint_Gamma_"+self.label]=str(self.gamma.get_scale())
        return output


class CysteineCrossLinkRestraint():
    def __init__(self,prots,filename,cbeta=False,
                    betatuple=(0.03, 0.1),
                    disttuple=(0.0,25.0, 1000),
                    omegatuple=(1.0, 1000.0, 50),
                    sigmatuple=(0.3,0.3,1),
                    betaissampled=True,
                    sigmaissampled=False,
                    weightissampled=True,
                    epsilonissampled=True
                    ):
    #the file must have residue1 chain1 residue2 chain2 fractionvalue epsilonname
    #epsilonname is a name for the epsilon particle that must be used for that particular
    #residue pair, eg, "Epsilon-Intra-Solvent", or "Epsilon-Solvent-Membrane", etc.
        global impisd2, tools
        import IMP.isd2 as impisd2
        import IMP.pmi.tools as tools

        self.prots=prots
        self.rs=IMP.RestraintSet('Cysteine_Crosslink')
        self.m=self.prots[0].get_model()
        self.cbeta=cbeta
        self.epsilonmaxtrans=0.01
        self.sigmamaxtrans=0.1
        self.betamaxtrans=0.01
        self.weightmaxtrans=0.1
        self.label="None"
        self.outputlevel="low"
        self.betaissampled=betaissampled
        self.sigmaissampled=sigmaissampled
        self.weightissampled=weightissampled
        self.epsilonissampled=epsilonissampled

        betalower=betatuple[0]
        betaupper=betatuple[1]
        beta=0.04
        sigma=10.0
        betangrid=30
        crossdataprior=1

        # beta
        self.beta=tools.SetupNuisance(self.m,beta,betalower,betaupper,betaissampled).get_particle()
        # sigma
        self.sigma=tools.SetupNuisance(self.m,sigma,sigmatuple[0],sigmatuple[1],sigmaissampled).get_particle()
        #population particle
        self.weight=tools.SetupWeight(self.m,weightissampled).get_particle()

        #read the file
        fl=open(filename,"r")
        #determine the upperlimit for epsilon
        #and also how many epsilon are needed
        self.epsilons={}
        data=[]
        for line in fl:

            t=line.split()


            if t[5] in self.epsilons:
                if 1.0-float(t[4])<=self.epsilons[t[5]].get_upper():
                    self.epsilons[t[5]].set_upper(1.0-float(t[4]))
            else:
                self.epsilons[t[5]]=tools.SetupNuisance(self.m,
                               0.01,0.01,1.0-float(t[4]),epsilonissampled).get_particle()
            up=self.epsilons[t[5]].get_upper()
            low=self.epsilons[t[5]].get_lower()
            if up<low: self.epsilons[t[5]].set_upper(low)

            data.append((int(t[0]),t[1],int(t[2]),t[3],float(t[4]),t[5]))
        fl.close()


        # create CrossLinkData



        if not self.cbeta:
            crossdata=tools.get_cross_link_data("cysteine","cysteine_CA_FES.txt.standard",
                                                   disttuple,omegatuple,sigmatuple,disttuple[1],disttuple[1],1)
        else:
            crossdata=tools.get_cross_link_data("cysteine","cysteine_CB_FES.txt.standard",
                                                   disttuple,omegatuple,sigmatuple,disttuple[1],disttuple[1],1)


        # create grids needed by CysteineCrossLinkData
        fmod_grid=tools.get_grid(0.0, 1.0, 300, True)
        omega2_grid=tools.get_log_grid(0.001, 10000.0, 100)
        beta_grid=tools.get_log_grid(betalower,betaupper,betangrid)

        for d in data:
            resid1=d[0]
            chain1=d[1]
            resid2=d[2]
            chain2=d[3]
            fexp=d[4]
            epslabel=d[5]

            # CysteineCrossLinkData

            ccldata=impisd2.CysteineCrossLinkData(fexp,fmod_grid,omega2_grid,beta_grid)

            ccl=impisd2.CysteineCrossLinkRestraint(self.beta,self.sigma,self.epsilons[epslabel],self.weight,crossdata,ccldata)

            failed=False
            for i,prot in enumerate(self.prots):

                if not self.cbeta:
                    p1=None
                    p2=None


                    p1=tools.select_calpha_or_residue(prot=prot,chain=chain1,
                                    resid=resid1,ObjectName="CysteineCrossLink:",SelectResidue=True)
                    if p1==None: failed=True

                    p2=tools.select_calpha_or_residue(prot=prot,chain=chain2,
                                    resid=resid2,ObjectName="CysteineCrossLink:",SelectResidue=True)
                    if p2==None: failed=True

                else:
                    #use cbetas
                    p1=[]
                    p2=[]
                    for t in range(-1,2):
                        p=tools.select_calpha_or_residue(prot=prot,chain=chain1,
                                    resid=resid1+t,ObjectName="CysteineCrossLink:",SelectResidue=False)
                        if p!=None:
                            p1.append(p)
                        else:
                            failed=True

                        p=tools.select_calpha_or_residue(prot=prot,chain=chain2,
                                    resid=resid2+t,ObjectName="CysteineCrossLink:",SelectResidue=False)
                        if p!=None:
                            p2.append(p)
                        else:
                            failed=True

                if not self.cbeta:
                    if (p1!=None and p2!=None):
                        ccl.add_contribution(p1,p2)
                        d1=IMP.core.XYZ(p1)
                        d2=IMP.core.XYZ(p2)

                        print "Distance_"+str(resid1)+"_"+chain1+":"+str(resid2)+"_"+chain2, IMP.core.get_distance(d1,d2)

                else:
                    if (len(p1)==3 and len(p2)==3):
                        ccl.add_contribution(p1,p2)

            if not failed:
                self.rs.add_restraint(ccl)
                ccl.set_name("CysteineCrossLink_"+str(resid1)+"_"+chain1+":"+str(resid2)+"_"+chain2)



    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_particles_to_sample(self):
        ps={}
        if self.epsilonissampled:
            for eps in self.epsilons.keys():
                ps["Nuisances_CysteineCrossLinkRestraint_epsilon_"+eps+"_"+self.label]=([self.epsilons[eps]],self.epsilonmaxtrans)
        if self.betaissampled:
            ps["Nuisances_CysteineCrossLinkRestraint_beta_"+self.label]=([self.beta],self.betamaxtrans)
        if self.weightissampled:
            ps["Weights_CysteineCrossLinkRestraint_"+self.label]=([self.weight],self.weightmaxtrans)
        if self.sigmaissampled:
            ps["Nuisances_CysteineCrossLinkRestraint_"+self.label]=([self.sigma],self.sigmamaxtrans)
        return ps

    def set_output_level(self,level="low"):
                #this might be "low" or "high"
        self.outputlevel=level

    def get_restraint(self):
        return self.rs

    def get_sigma(self):
        return self.sigma

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)         
        output["CysteineCrossLinkRestraint_Score_"+self.label]=str(score)
        output["CysteineCrossLinkRestraint_sigma_"+self.label]=str(self.sigma.get_scale())
        for eps in self.epsilons.keys():
            output["CysteineCrossLinkRestraint_epsilon_"+eps+"_"+self.label]=str(self.epsilons[eps].get_scale())
        output["CysteineCrossLinkRestraint_beta_"+self.label]=str(self.beta.get_scale())
        for n in range(self.weight.get_number_of_states()):
            output["CysteineCrossLinkRestraint_weights_"+str(n)+"_"+self.label]=str(self.weight.get_weight(n))

        if self.outputlevel=="high":
            for rst in self.rs.get_restraints():
                output["CysteineCrossLinkRestraint_Total_Frequency_"+
                       impisd2.CysteineCrossLinkRestraint.get_from(rst).get_name()+
                       "_"+self.label]=impisd2.CysteineCrossLinkRestraint.get_from(rst).get_model_frequency()
                output["CysteineCrossLinkRestraint_Standard_Error_"+
                       impisd2.CysteineCrossLinkRestraint.get_from(rst).get_name()+"_"
                       +self.label]=impisd2.CysteineCrossLinkRestraint.get_from(rst).get_standard_error()
                if len(self.prots)>1:
                    for i in range(len(self.prots)):
                        output["CysteineCrossLinkRestraint_Frequency_Contribution_"+
                        impisd2.CysteineCrossLinkRestraint.get_from(rst).get_name()+
                        "_State_"+str(i)+"_"+self.label]=impisd2.CysteineCrossLinkRestraint.get_from(rst).get_frequencies()[i]

        return output

class GaussianEMRestraint():

    def __init__(self,prot,map_anchors_fn,segment_anchors=None,segment_parts=None,rigid=True):

        global sys, impisd2, tools
        import sys
        import IMP.isd2 as impisd2
        import IMP.pmi.tools as tools
        #import IMP.multifit
        
        if segment_anchors==None: segment_anchors=[]
        if segment_parts==None: segment_parts=[]        
        
        #dcoords=IMP.multifit.read_anchors_data(map_anchors_fn).points_
        self.prot=prot
        sel=IMP.atom.Selection(self.prot)
        #sel.set_atom_type(IMP.atom.AT_CA)
        self.model_anchors=sel.get_selected_particles()


        self.m=self.prot.get_model()

        data = open(map_anchors_fn)
        D = data.readlines()
        data.close()
        dcoords={}
        for d in D:
            d=d.strip().split('|')
            if len(d)==6: dcoords[int(d[1])] = IMP.algebra.Vector3D(float(d[2]),float(d[3]),float(d[4]))

        # parameters
        self.model_sigmas=[15.0]*len(self.model_anchors)

        #self.model_sigmas=[float(anch.get_as_xyzr().get_radius()) for anch in self.segment_parts]
        self.model_weights=[1.0]*len(self.model_anchors)
        self.density_sigmas=[15.0]*len(segment_anchors)
        self.density_weights=[5.0]*len(segment_anchors)
        self.sigmamaxtrans=0.1
        self.sigmamin=1.
        self.sigmamax=100.0
        self.sigmainit=10.0
        self.cutoff_dist_for_container=10.0
        self.rigid=rigid
        self.segment_anchors=segment_anchors
        self.segment_parts=segment_parts
        self.tabexp=True


        self.density_anchors=[]
        for d in dcoords:
            if self.segment_anchors==[]:
                p=IMP.Particle(self.m)
                self.density_anchors.append(p)
                IMP.core.XYZR.setup_particle(p,\
                                         IMP.algebra.Sphere3D(d,\
                                         self.density_sigmas[nd]*1.5))
            else:
                if d in self.segment_anchors:
                    p=IMP.Particle(self.m)
                    self.density_anchors.append(p)
                    IMP.core.XYZR.setup_particle(p,\
                                         IMP.algebra.Sphere3D(dcoords[d],\
                                         self.density_sigmas[0]*1.5))

        for np,p in enumerate(self.model_anchors):
            self.model_sigmas[np]=IMP.core.XYZR(p).get_radius()/1.5

        self.sigmaglobal=tools.SetupNuisance(self.m,self.sigmainit,
                 self.sigmamin,self.sigmamax,True).get_particle()
        print 'setting up restraint'

        self.gaussianEM_restraint=impisd2.GaussianEMRestraint(
            self.model_anchors,self.model_sigmas,self.model_weights,
            self.density_anchors,self.density_sigmas,self.density_weights,
            self.sigmaglobal.get_particle(),self.cutoff_dist_for_container,
            self.rigid,self.tabexp)
        print 'done setup'
        self.rs = IMP.RestraintSet('GaussianEMRestraint')
        self.rs.add_restraint(self.gaussianEM_restraint)





    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_particles_to_sample(self):
        ps={}
        ps["Nuisances_GaussianEMRestraint_sigma_"+self.label]=([self.sigmaglobal],self.sigmamaxtrans)
        return ps

    def get_hierarchy(self):
        return self.prot

    def get_restraint_set(self):
        return self.rs

    def get_output(self):
        self.m.update()
        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)         
        output["GaussianEMRestraint_"+self.label]=str(self.rs.unprotected_evaluate(None))
        output["GaussianEMRestraint_sigma_"+self.label]=str(self.sigmaglobal.get_scale())
        return output
