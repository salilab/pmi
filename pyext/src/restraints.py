#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container



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


##############################################################


class ConnectivityCrossLinkMS():
    '''
    this restraint allows ambiguous crosslinking between multiple copies
    it is a variant of the SimplifiedCrossLinkMS
    '''

    def __init__(self,prot,restraints_file,expdistance,strength=None,resolution=None):

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

            #calculate the radii to estimate the slope of the restraint
            if self.strength==None:
              rad1=0
              rad2=0
              for p in ps1:
                rad1+=IMP.pmi.Uncertainty(p).get_uncertainty()

              for p in ps2:
                rad2+=IMP.pmi.Uncertainty(p).get_uncertainty()

              rad1=rad1/len(ps1)
              rad2=rad2/len(ps2)

              self.strength=1/(rad1**2+rad2**2)

            
            sels=[s1,s2]
            cr = IMP.atom.create_connectivity_restraint(sels, self.expdistance,self.strength)
    
            self.rs.add_restraint(cr)
            self.pairs.append((ps1,hrc1,c1,r1,ps2,hrc2,c2,r2,cr))

    def plot_restraint(self,uncertainty1,uncertainty2,maxdist=50,npoints=10):
        import IMP.pmi.output as output
        
        p1=IMP.Particle(self.m)
        p2=IMP.Particle(self.m)
        d1=IMP.core.XYZR.setup_particle(p1)
        d2=IMP.core.XYZR.setup_particle(p2)
        d1.set_radius(uncertainty1)
        d2.set_radius(uncertainty2)                
        s1=IMP.atom.Selection(p1)
        s2=IMP.atom.Selection(p2)
        sels=[s1,s2]
        strength=1/(uncertainty1**2+uncertainty2**2)
        cr = IMP.atom.create_connectivity_restraint(sels, self.expdistance,strength)
        dists=[]
        scores=[]
        for i in range(npoints):
            d2.set_coordinates(IMP.algebra.Vector3D(maxdist/npoints*float(i),0,0))
            dists.append(IMP.core.get_distance(d1,d2))
            scores.append(cr.unprotected_evaluate(None))
        output.plot_xy_data(dists,scores)
            
            

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
        output["ConnectivityCrossLinkMS_Score_"+self.label]=str(score)
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

    def __init__(self,prot,restraints_file,expdistance,strength,resolution=None, columnmapping=None):
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
        self.already_added_pairs={}
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

            names = [h.get_name() for h in self.prot.get_children()]
            if c1 not in names or c2 not in names:
                print c1, ' OR ', c2, ' not in the model!'; continue
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
               print "particles are: ", [p.get_name() for p in ps1]
               exit()
            elif len(ps1)==0:
               print "SimplifiedCrossLinkMS: WARNING> residue %d of chain %s is not there" % (r1,c1)
               continue               

            if len(ps2)>1:
               print "SimplifiedCrossLinkMS: ERROR> residue %d of chain %s selects multiple particles"  % (r2,c2)
               print "particles are: ", [p.get_name() for p in ps2]
               exit()
            elif len(ps2)==0:
               print "SimplifiedCrossLinkMS: WARNING> residue %d of chain %s is not there" % (r2,c2)
               continue 

            
            p1=ps1[0]
            p2=ps2[0]
            
            if (p1,p2) in self.already_added_pairs:
               dr=self.already_added_pairs[(p1,p2)]
               weight=dr.get_weight()
               dr.set_weight(weight+1.0)
               print "SimplifiedCrossLinkMS> crosslink %d %s %d %s was already found, adding 1.0 to the weight, weight is now %d"  % (r1,c1,r2,c2,weight+1.0)
               continue
            
            else:
                        
              limit=self.strength*(self.expdistance+15)**2+10.0
              hub= IMP.core.TruncatedHarmonicUpperBound(self.expdistance,self.strength,self.expdistance+15.,limit)
              df= IMP.core.SphereDistancePairScore(hub)
              dr= IMP.core.PairRestraint(df, (p1, p2))
              dr.set_name(c1+":"+str(r1)+"-"+c2+":"+str(r2))
            
              self.rs.add_restraint(dr)
              self.pairs.append((p1,p2,dr,r1,c1,r2,c2))
              self.already_added_pairs[(p1,p2)]=dr
              self.already_added_pairs[(p2,p1)]=dr

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

    def plot_restraint(self,radius1,radius2,maxdist=50,npoints=10):
        import IMP.pmi.output as output
        
        p1=IMP.Particle(self.m)
        p2=IMP.Particle(self.m)
        d1=IMP.core.XYZR.setup_particle(p1)
        d2=IMP.core.XYZR.setup_particle(p2)
        d1.set_radius(radius1)
        d2.set_radius(radius2)                
        s1=IMP.atom.Selection(p1)
        s2=IMP.atom.Selection(p2)
        #limit=self.strength*(self.expdistance+1)**2+10.0
        #hub= IMP.core.TruncatedHarmonicUpperBound(self.expdistance,self.strength,self.expdistance+5.,limit)
        #df= IMP.core.SphereDistancePairScore(hub)
        #dr= IMP.core.PairRestraint(df, (p1, p2))
        dr=IMP.pmi.SigmoidRestraintSphere(self.m, p1, p2, self.expdistance, 1.0, 5.0) 
        dists=[]
        scores=[]
        for i in range(npoints):
            d2.set_coordinates(IMP.algebra.Vector3D(maxdist/npoints*float(i),0,0))
            dists.append(IMP.core.get_distance(d1,d2))
            scores.append(dr.unprotected_evaluate(None))
        output.plot_xy_data(dists,scores)

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


###############################################################


class SigmoidCrossLinkMS():

    def __init__(self,prot,restraints_file,inflection,slope,amplitude,resolution=None, columnmapping=None):
        #columnindexes is a list of column indexes for protein1, protein2, residue1, residue2
        #by default column 0 = protein1; column 1 = protein2; column 2 = residue1; column 3 = residue2


        
        if columnmapping==None:
           columnmapping={}
           columnmapping["Protein1"]=0
           columnmapping["Protein2"]=1
           columnmapping["Residue1"]=2
           columnmapping["Residue2"]=3
        
        self.rs=IMP.RestraintSet('data')
        self.rslin=IMP.RestraintSet('data')
        self.rssig=IMP.RestraintSet('data')
        self.weight=1.0
        self.prot=prot
        self.label="None"
        self.pairs=[]
        self.already_added_pairs={}
        self.m=self.prot.get_model()

        self.outputlevel="low"

        #small linear contribution for long range
        h=IMP.core.Linear(0,0.07)
        dps2=IMP.core.DistancePairScore(h)


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


            #hrc1 = [h for h in self.prot.get_children() if h.get_name()==c1][0]
            s1=IMP.atom.Selection(self.prot,molecule=c1, residue_index=r1)
            ps1=s1.get_selected_particles()

            #hrc2 = [h for h in self.prot.get_children() if h.get_name()==c2][0]
            s2=IMP.atom.Selection(self.prot,molecule=c2, residue_index=r2)
            ps2=s2.get_selected_particles()
       
            if resolution!=None: 
              #get the intersection to remove redundant particles
              ps1=(list(set(ps1) & set(particles)))
              ps2=(list(set(ps2) & set(particles)))
              
            if len(ps1)>1:
               print "SigmoidCrossLinkMS: ERROR> residue %d of chain %s selects multiple particles %s"  % (r1,c1,str(ps1))
               exit()
            elif len(ps1)==0:
               print "SigmoidCrossLinkMS: WARNING> residue %d of chain %s is not there" % (r1,c1)
               continue               

            if len(ps2)>1:
               print "SigmoidCrossLinkMS: ERROR> residue %d of chain %s selects multiple particles %s"  % (r2,c2,str(ps2))
               exit()
            elif len(ps2)==0:
               print "SigmoidCrossLinkMS: WARNING> residue %d of chain %s is not there" % (r2,c2)
               continue 

            
            p1=ps1[0]
            p2=ps2[0]
            
            if (p1,p2) in self.already_added_pairs:
               dr=self.already_added_pairs[(p1,p2)]
               weight=dr.get_weight()
               dr.increment_amplitude(amplitude)
               print "SigmoidCrossLinkMS> crosslink %d %s %d %s was already found, adding %d to the amplitude, amplitude is now %d"  % (r1,c1,r2,c2,amplitude,dr.get_amplitude())
               dr.set_name(c1+":"+str(r1)+"-"+c2+":"+str(r2)+"-ampl:"+str(dr.get_amplitude()))
               continue
            
            else:
                        
 
              dr= IMP.pmi.SigmoidRestraintSphere(self.m, p1, p2, inflection, slope, amplitude)
              dr.set_name(c1+":"+str(r1)+"-"+c2+":"+str(r2)+"-ampl:"+str(dr.get_amplitude()))
            
              self.rs.add_restraint(dr)
              self.rssig.add_restraint(dr)
              
              pr=IMP.core.PairRestraint(dps2,IMP.ParticlePair(p1,p2)) 

              self.rs.add_restraint(pr)
              self.rslin.add_restraint(pr)
              
              self.pairs.append((p1,p2,dr,pr,r1,c1,r2,c2))
              self.already_added_pairs[(p1,p2)]=dr
              self.already_added_pairs[(p2,p1)]=dr

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
        output["SigmoidCrossLinkMS_Score_"+self.label]=str(score)
        output["SigmoidCrossLinkMS_Sigmoidal_Score_"+self.label]=self.rssig.unprotected_evaluate(None)   
        output["SigmoidCrossLinkMS_Linear_Score_"+self.label]=self.rslin.unprotected_evaluate(None)
        for i in range(len(self.pairs)):

            p0=self.pairs[i][0]
            p1=self.pairs[i][1]
            crosslinker='standard'
            ln=self.pairs[i][2]
            pr=self.pairs[i][3]
            resid1=self.pairs[i][4]
            chain1=self.pairs[i][5]
            resid2=self.pairs[i][6]
            chain2=self.pairs[i][7]

            label=str(resid1)+":"+chain1+"_"+str(resid2)+":"+chain2
            output["SigmoidCrossLinkMS_Score_"+crosslinker+"_"+label]=str(self.weight*ln.unprotected_evaluate(None))
            output["SigmoidCrossLinkMS_Linear_Score_"+crosslinker+"_"+label]=str(self.weight*pr.unprotected_evaluate(None))
            d0=IMP.core.XYZR(p0)
            d1=IMP.core.XYZR(p1)
            output["SigmoidCrossLinkMS_Distance_"+label]=str(IMP.core.get_distance(d0,d1))

        return output

###############################################################

class ISDCrossLinkMS():

    def __init__(self,prot,restraints_file,length,resolution=None, columnmapping=None):
        #columnindexes is a list of column indexes for protein1, protein2, residue1, residue2
        #by default column 0 = protein1; column 1 = protein2; column 2 = residue1; column 3 = residue2;
        # column 4 = idscores

        global impisd2, tools
        import IMP.isd2 as impisd2
        import IMP.pmi.tools as tools
        
        if columnmapping==None:
           columnmapping={}
           columnmapping["Protein1"]=0
           columnmapping["Protein2"]=1
           columnmapping["Residue1"]=2
           columnmapping["Residue2"]=3
           columnmapping["IDScore"]=4
        
        self.rs=IMP.RestraintSet('data')
        self.rspsi=IMP.RestraintSet('prior_psi')
        self.rssig=IMP.RestraintSet('prior_sigmas')        
        self.rslin=IMP.RestraintSet('prior_linear')
        self.prot=prot
        self.label="None"
        self.pairs=[]
        self.m=self.prot.get_model()
        self.sigma_dictionary={}
        self.psi_dictionary={}
        
        self.ids_map=tools.map()
        #self.ids_map.set_map_element(20.0,0.05)      
        self.ids_map.set_map_element(35.0,0.01)   

        self.radius_map=tools.map()
        self.radius_map.set_map_element(2.5,2.5)
        self.radius_map.set_map_element(5.0,5.0)      
        self.radius_map.set_map_element(7.5,7.5)  
        self.radius_map.set_map_element(10,10)
        self.radius_map.set_map_element(15,15)      
        self.radius_map.set_map_element(20,20) 
        
        self.outputlevel="low"

        #small linear contribution for long range
        h=IMP.core.Linear(0,0.07)
        dps2=IMP.core.DistancePairScore(h)

        #fill the cross-linker pmfs
        #to accelerate the init the list listofxlinkertypes might contain only yht needed crosslinks
        protein1=columnmapping["Protein1"]
        protein2=columnmapping["Protein2"]
        residue1=columnmapping["Residue1"]
        residue2=columnmapping["Residue2"]        
        idscore=columnmapping["IDScore"]
        
        print idscore
        
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
            
            if tokens[idscore]=="High" : ids=1
            elif tokens[idscore]=="Low" : ids=0
            else: ids=float(tokens[idscore])


            s1=IMP.atom.Selection(self.prot,molecule=c1, residue_index=r1)
            ps1=s1.get_selected_particles()
 
            s2=IMP.atom.Selection(self.prot,molecule=c2, residue_index=r2)
            ps2=s2.get_selected_particles()
       
            if resolution!=None: 
              #get the intersection to remove redundant particles
              ps1=(list(set(ps1) & set(particles)))
              ps2=(list(set(ps2) & set(particles)))
              
            if len(ps1)>1:
               print "ISDCrossLinkMS: ERROR> residue %d of chain %s selects multiple particles %s"  % (r1,c1,str(ps1))
               exit()
            elif len(ps1)==0:
               print "ISDCrossLinkMS: WARNING> residue %d of chain %s is not there" % (r1,c1)
               continue               

            if len(ps2)>1:
               print "ISDCrossLinkMS: ERROR> residue %d of chain %s selects multiple particles %s"  % (r2,c2,str(ps2))
               exit()
            elif len(ps2)==0:
               print "ISDCrossLinkMS: WARNING> residue %d of chain %s is not there" % (r2,c2)
               continue 

            p1=ps1[0]
            p2=ps2[0]
            
            #remove in the future!!!
            if p1==p2: continue
            
            dr= impisd2.CrossLinkMSRestraint(self.m, length)

            
            #sigma1=self.get_sigma(IMP.pmi.Resolution(p1))[0]
            #sigma2=self.get_sigma(IMP.pmi.Resolution(p2))[0]    
            mappedr1=self.radius_map.get_map_element(IMP.pmi.Uncertainty(p1).get_uncertainty())
            sigma1=self.get_sigma(mappedr1)[0]
            mappedr2=self.radius_map.get_map_element(IMP.pmi.Uncertainty(p2).get_uncertainty())
            sigma2=self.get_sigma(mappedr2)[0]
            psival=self.ids_map.get_map_element(ids)
            psi=self.get_psi(psival)[0]
            p1i=p1.get_index()
            p2i=p2.get_index()
            s1i=sigma1.get_particle().get_index()
            s2i=sigma2.get_particle().get_index()
            psii=psi.get_particle().get_index()
            dr.add_contribution((p1i,p2i),(s1i,s2i),psii)
            
            print "--------------"
            print "ISDCrossLinkMS: generating cross-link restraint between"
            print "ISDCrossLinkMS: residue %d of chain %s and residue %d of chain %s" % (r1,c1,r2,c2)
            print "ISDCrossLinkMS: with sigma1 %f  sigma2 %f psi %s" % (mappedr1,mappedr2,psival) 
            print "ISDCrossLinkMS: between particles %s and %s" % (p1.get_name(),p2.get_name())
            
            self.rs.add_restraint(dr)
            self.rssig.add_restraint(dr)

            pr=IMP.core.PairRestraint(dps2,IMP.ParticlePair(p1,p2)) 

            self.rs.add_restraint(pr)
            self.rslin.add_restraint(pr)

            #check if the two residues belong to the same rigid body
            

            if(IMP.core.RigidMember.particle_is_instance(p1) and
               IMP.core.RigidMember.particle_is_instance(p2) and
               IMP.core.RigidMember(p1).get_rigid_body() ==
               IMP.core.RigidMember(p2).get_rigid_body()):
               xlattribute="intrarb"
            else:
               xlattribute="interrb"


            dr.set_name(xlattribute+"-"+c1+":"+str(r1)+"-"+c2+":"+str(r2))
          
            self.pairs.append((p1,p2,dr,r1,c1,r2,c2,xlattribute,mappedr1,mappedr2,psival))
        
    def create_sigma(self,resolution):
        self.sigmainit=resolution+2.0
        self.sigmaissampled=True             
        self.sigmaminnuis=0.0000001
        self.sigmamaxnuis=1000.0
        self.sigmamin=    resolution-resolution/2
        self.sigmamax=    500.0
        self.sigmatrans=  0.2
        self.sigma=tools.SetupNuisance(self.m,self.sigmainit,
             self.sigmaminnuis,self.sigmamaxnuis,self.sigmaissampled).get_particle()
        self.sigma_dictionary[resolution]=(self.sigma,self.sigmatrans,self.sigmaissampled)    
        self.rssig.add_restraint(impisd2.UniformPrior(self.sigma,1000000000.0,self.sigmamax,self.sigmamin))
        #self.rssig.add_restraint(impisd2.JeffreysRestraint(self.sigma))        
        
    def get_sigma(self,resolution):
        if not resolution in self.sigma_dictionary:
           self.create_sigma(resolution)
        return self.sigma_dictionary[resolution]


    def create_psi(self,value):
        self.psiinit=value
        self.psiissampled=True                     
        self.psiminnuis=0.0000001
        self.psimaxnuis=0.4999999
        self.psimin=    0.01
        self.psimax=    1.0
        self.psitrans=  0.01 
        self.psi=tools.SetupNuisance(self.m,self.psiinit,
             self.psiminnuis,self.psimaxnuis,self.psiissampled).get_particle()
        self.psi_dictionary[value]=(self.psi,self.psitrans,self.psiissampled)    
        self.rspsi.add_restraint(impisd2.UniformPrior(self.psi,1000000000.0,self.psimax,self.psimin))
        self.rspsi.add_restraint(impisd2.JeffreysRestraint(self.psi))
        
    def get_psi(self,value):
        if not value in self.psi_dictionary:
           self.create_psi(value)
        return self.psi_dictionary[value]


    def set_label(self,label):
        self.label=label

    def add_to_model(self):
        self.m.add_restraint(self.rs)
        self.m.add_restraint(self.rspsi)
        self.m.add_restraint(self.rssig)

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


    def get_output(self):
        #content of the crosslink database pairs
        #self.pairs.append((p1,p2,dr,r1,c1,r2,c2))
        self.m.update()

        output={}
        score=self.rs.unprotected_evaluate(None)
        output["_TotalScore"]=str(score)            
        output["ISDCrossLinkMS_Data_Score_"+self.label]=str(score)
        output["ISDCrossLinkMS_PriorSig_Score_"+self.label]=self.rssig.unprotected_evaluate(None)   
        output["ISDCrossLinkMS_PriorPsi_Score_"+self.label]=self.rspsi.unprotected_evaluate(None)  
        for i in range(len(self.pairs)):

            p0=self.pairs[i][0]
            p1=self.pairs[i][1]
            crosslinker='standard'
            ln=self.pairs[i][2]
            resid1=self.pairs[i][3]
            chain1=self.pairs[i][4]
            resid2=self.pairs[i][5]
            chain2=self.pairs[i][6]
            attribute=self.pairs[i][7]
            rad1=self.pairs[i][8]
            rad2=self.pairs[i][9]
            psi=self.pairs[i][10]
            
            label=attribute+"-"+str(resid1)+":"+chain1+"_"+str(resid2)+":"+chain2+"-"+str(rad1)+"-"+str(rad2)+"-"+str(psi)
            output["ISDCrossLinkMS_Score_"+crosslinker+"_"+label]=str(ln.unprotected_evaluate(None))
            d0=IMP.core.XYZ(p0)
            d1=IMP.core.XYZ(p1)
            output["ISDCrossLinkMS_Distance_"+label]=str(IMP.core.get_distance(d0,d1))

        for psiindex in self.psi_dictionary:
            output["ISDCrossLinkMS_Psi_"+str(psiindex)+"_"+self.label]=str(self.psi_dictionary[psiindex][0].get_scale())

        for resolution in self.sigma_dictionary:
            output["ISDCrossLinkMS_Sigma_"+str(resolution)+"_"+self.label]=str(self.sigma_dictionary[resolution][0].get_scale())

        return output

    def get_particles_to_sample(self):
        ps={}

        for resolution in self.sigma_dictionary:
          if self.sigma_dictionary[resolution][2]:
            ps["Nuisances_ISDCrossLinkMS_Sigma_"+str(resolution)+"_"+self.label]=\
                      ([self.sigma_dictionary[resolution][0]],self.sigma_dictionary[resolution][1])
        
        for psiindex in self.psi_dictionary:
          if self.psi_dictionary[psiindex][2]:
            ps["Nuisances_ISDCrossLinkMS_Psi_"+str(psiindex)+"_"+self.label]=([self.psi_dictionary[psiindex][0]],self.psi_dictionary[psiindex][1])
        
        return ps  


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



############################################################

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

    def __init__(self,prot,map_anchors_fn,segment_anchors=None,segment_parts=None,rigid=True,resolution=None):
        #segment parts should be a list containing protein names or a tuple with protein names and residue ranges:
        # [("ABC",1,100),"CYT","CDR"]

        global sys, impisd2, tools
        import sys
        import IMP.isd2 as impisd2
        import IMP.pmi.tools as tools
        #import IMP.multifit
        
        if segment_anchors==None: segment_anchors=[]
        if segment_parts==None: segment_parts=[]        
        
        #dcoords=IMP.multifit.read_anchors_data(map_anchors_fn).points_
        self.prot=prot
        self.m=self.prot.get_model()

        

        if resolution!=None:
           hierarchy_anchors=[]
           for prot in self.prot.get_children():
             hierarchy_anchors+=IMP.pmi.tools.get_particles_by_resolution(prot,resolution)         

        
        if segment_parts!=None:
            model_anchors=[]
            for seg in segment_parts:
                if type(seg)==str:
                   s=IMP.atom.Selection(self.prot,molecule=seg)
                   ps=s.get_selected_particles()
                   model_anchors+=ps
                elif type(seg)==tuple:
                   s=IMP.atom.Selection(self.prot,molecule=seg[0],residue_indexes=range(seg[1],seg[2]+1))
                   ps=s.get_selected_particles()
                   model_anchors+=ps                
            if resolution!=None:
              #get the intersection to remove redundant particles
              self.model_anchors=(list(set(model_anchors) & set(hierarchy_anchors)))

        for p in self.model_anchors:
            print p.get_name()  



        data = open(map_anchors_fn)
        D = data.readlines()
        data.close()
        dcoords={}
        for d in D:
            d=d.strip().split('|')
            if len(d)==6: dcoords[int(d[1])] = IMP.algebra.Vector3D(float(d[2]),float(d[3]),float(d[4]))

        # parameters
        self.model_sigmas=[15.0]*len(self.model_anchors)
        
        self.model_sigmas=[]
        self.model_weights=[]
        for p in self.model_anchors: self.model_sigmas.append(IMP.core.XYZR(p).get_radius())
        for p in self.model_weights: self.model_weights.append(len(IMP.atom.Fragment(p).get_residue_indexes()))



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
