#!/usr/bin/env python
import IMP
import IMP.algebra


class Stopwatch():

    def __init__(self,isdelta=True):
        global time
        import time

        self.starttime=time.clock()
        self.label="None"
        self.isdelta=isdelta

    def set_label(self,labelstr):
        self.label=labelstr

    def get_output(self):
        output={}
        if self.isdelta:
            newtime=time.clock()
            output["Stopwatch_"+self.label+"_delta_seconds"]=str(newtime-self.starttime)
            self.starttime=newtime
        else:
            output["Stopwatch_"+self.label+"_elapsed_seconds"]=str(time.clock()-self.starttime)
        return output


class SetupNuisance():
    def __init__(self,m,initialvalue,minvalue,maxvalue,isoptimized=True):
        global impisd2
        import IMP.isd2 as impisd2
        
        nuisance=impisd2.Scale.setup_particle(IMP.Particle(m),initialvalue)
        nuisance.set_lower(minvalue)
        nuisance.set_upper(maxvalue)

        m.add_score_state(IMP.core.SingletonConstraint(impisd2.NuisanceRangeModifier(),None,nuisance))
        nuisance.set_is_optimized(nuisance.get_nuisance_key(),isoptimized)
        self.nuisance=nuisance

    def get_particle(self):
        return self.nuisance

class SetupWeight():

    def __init__(self,m,isoptimized=True):
        global impisd2
        import IMP.isd2 as impisd2    
        pw=IMP.Particle(m)
        self.weight=impisd2.Weight.setup_particle(pw)
        self.weight.set_weights_are_optimized(True)

    def get_particle(self):
        return self.weight



class ParticleToSampleList():


    def __init__(self,label="None"):

        self.dictionary_particle_type={}
        self.dictionary_particle_transformation={}
        self.dictionary_particle_name={}
        self.label=label

    def add_particle(self,particle,particle_type,particle_transformation,name):
        if not particle_type in ["Rigid_Bodies","Floppy_Bodies","Nuisances","X_coord","Weights"]:
            print "ParticleToSampleList: not the right particle type"
            exit()
        else:
            self.dictionary_particle_type[particle]=particle_type
            if particle_type=="Rigid_Bodies":
                if type(particle_transformation)==tuple and len(particle_transformation)==2 and type(particle_transformation[0])==float and type(particle_transformation[1])==float:
                    self.dictionary_particle_transformation[particle]=particle_transformation
                    self.dictionary_particle_name[particle]=name
                else:
                    print "ParticleToSampleList: not the right transformation format for Rigid_Bodies, should be a tuple a floats"
                    exit()
            else:
                if type(particle_transformation)==float:
                    self.dictionary_particle_transformation[particle]=particle_transformation
                    self.dictionary_particle_name[particle]=name
                else:
                    print "ParticleToSampleList: not the right transformation format sould be a float"
                    exit()

    def get_particles_to_sample(self):
        ps={}
        for particle in self.dictionary_particle_type:
            key=self.dictionary_particle_type[particle]+"ParticleToSampleList_"+self.dictionary_particle_name[particle]+"_"+self.label
            value=([particle],self.dictionary_particle_transformation[particle])
            ps[key]=value
        return ps

class Variance():
    def __init__(self, model, tau, niter, prot, th_profile, write_data=False):
        global sqrt,os,random
        from math import sqrt
        import os
        import random

        self.model = model
        self.write_data=write_data
        self.tau = tau
        self.niter = niter
        #! select particles from the model
        particles=IMP.atom.get_by_type(prot, IMP.atom.ATOM_TYPE)
        self.particles = particles
        #store reference coordinates and theoretical profile
        self.refpos = [ IMP.core.XYZ(p).get_coordinates() for p in particles ]
        self.model_profile = th_profile

    def perturb_particles(self, perturb=True):
        for i,p in enumerate(self.particles):
            newpos = array(self.refpos[i])
            if perturb:
                newpos += random.normal(0,self.tau,3)
            newpos = IMP.algebra.Vector3D(newpos)
            IMP.core.XYZ(p).set_coordinates(newpos)

    def get_profile(self):
        model_profile = self.model_profile
        p=model_profile.calculate_profile(self.particles, IMP.saxs.CA_ATOMS)
        return array( [ model_profile.get_intensity(i) for i in
                        xrange(model_profile.size()) ] )

    def init_variances(self):
        #create placeholders
        N = self.model_profile.size()
        a = self.profiles[0][:]
        self.m = matrix(a).T # Nx1
        self.V = self.m * self.m.T
        self.normm = linalg.norm(self.m)
        self.normV = linalg.norm(self.V)

    def update_variances(self):
        a = matrix(self.profiles[-1]) #1xN
        n = float(len(self.profiles))
        self.m = a.T/n + (n-1)/n * self.m
        self.V = a.T*a + self.V
        self.oldnormm = self.normm
        self.oldnormV = self.normV
        self.normm = linalg.norm(self.m)
        self.normV = linalg.norm(self.V)
        self.diffm = (self.oldnormm-self.normm)/self.oldnormm
        self.diffV = (self.oldnormV-self.normV)/self.oldnormV

    def get_direct_stats(self, a):
        nq = len(a[0])
        nprof = len(a)
        m = [0]*nq
        for prof in a:
            for q,I in enumerate(prof):
                m[q] += I
        m = array(m)/nprof
        V = matrix(a)
        V = V.T*V
        Sigma = (matrix(a-m))
        Sigma = Sigma.T*Sigma/(nprof-1)
        mi = matrix(diag(1./m))
        Sigmarel = mi.T*Sigma*mi
        return m,V,Sigma,Sigmarel

    def store_data(self):
        if not os.path.isdir('data'):
            os.mkdir('data')
        profiles = matrix(self.profiles)
        self.directm, self.directV, self.Sigma, self.Sigmarel = \
                self.get_direct_stats(array(profiles))
        directV = self.directV
        #print "V comparison",(linalg.norm(directV-self.V)/self.normV)
        save('data/profiles', profiles)
        #absolute profile differences
        fl=open('data/profiles.dat','w')
        for i,l in enumerate(array(profiles).T):
            self.model_profile.get_q(i)
            fl.write('%s ' % i)
            for k in l:
                fl.write('%s ' % (k-self.directm[i]))
            fl.write('\n')
        #relative profile differences
        fl=open('data/profiles_rel.dat','w')
        for i,l in enumerate(array(profiles).T):
            self.model_profile.get_q(i)
            fl.write('%s ' % i)
            for k in l:
                fl.write('%s ' % ((k-self.directm[i])/self.directm[i]))
            fl.write('\n')
        save('data/m', self.directm)
        save('data/V', self.directV)
        Sigma = self.Sigma
        save('data/Sigma', Sigma)
        #Sigma matrix
        fl=open('data/Sigma.dat', 'w')
        model_profile = self.model_profile
        for i in xrange(model_profile.size()):
            qi = model_profile.get_q(i)
            for j in xrange(model_profile.size()):
                qj = model_profile.get_q(j)
                vij = self.Sigma[i,j]
                fl.write('%s %s %s\n' % (qi,qj,vij))
            fl.write('\n')
        #Sigma eigenvalues
        fl=open('data/eigenvals','w')
        for i in linalg.eigvalsh(Sigma):
            fl.write('%s\n' % i)
        Sigmarel = self.Sigmarel
        save('data/Sigmarel', Sigmarel)
        #Sigmarel matrix
        fl=open('data/Sigmarel.dat', 'w')
        model_profile = self.model_profile
        for i in xrange(model_profile.size()):
            qi = model_profile.get_q(i)
            for j in xrange(model_profile.size()):
                qj = model_profile.get_q(j)
                vij = self.Sigmarel[i,j]
                fl.write('%s %s %s\n' % (qi,qj,vij))
            fl.write('\n')
        #Sigma eigenvalues
        fl=open('data/eigenvals_rel','w')
        for i in linalg.eigvalsh(Sigmarel):
            fl.write('%s\n' % i)
        #mean profile
        fl=open('data/mean.dat','w')
        for i in xrange(len(self.directm)):
            qi = self.model_profile.get_q(i)
            fl.write('%s ' % qi)
            fl.write('%s ' % self.directm[i])
            fl.write('%s ' % sqrt(self.Sigma[i,i]))
            fl.write('\n')

    def try_chol(self, jitter):
        Sigma=self.Sigma
        try:
            linalg.cholesky(Sigma+matrix(eye(len(Sigma)))*jitter)
        except linalg.LinAlgError:
            print "Decomposition failed with jitter =",jitter
            return
        print "Successful decomposition with jitter =",jitter

    def run(self):
        self.profiles = [self.get_profile()]
        #self.init_variances()
        for n in xrange(self.niter):
            self.perturb_particles()
            self.profiles.append(self.get_profile())
            #self.update_variances()
            #profiles = matrix(self.profiles)
            #print n,self.diffm,self.diffV
        #print
        #
        if self.write_data:
            self.store_data()
        #self.try_chol(0.)
        #for i in logspace(-7,0,num=8):
        #    self.try_chol(i)

    def get_cov(self, relative=True):
        if not relative:
            return self.Sigma
        else:
            return self.Sigmarel

    #-------------------------------

def get_cross_link_data(directory,filename,(distmin,distmax,ndist),
                                             (omegamin,omegamax,nomega),
                                            (sigmamin,sigmamax,nsigma),
                                            don=None,doff=None,prior=0,type_of_profile="gofr"):
    
    import IMP.isd2
    
    filen=IMP.isd2.get_data_path("CrossLinkPMFs.dict")
    xlpot=open(filen)

    for line in xlpot:
        dictionary=eval(line)
        break

    xpot=dictionary[directory][filename]["distance"]
    pot=dictionary[directory][filename][type_of_profile]

    dist_grid=get_grid(distmin, distmax, ndist, False)
    omega_grid=get_log_grid(omegamin, omegamax, nomega)
    sigma_grid=get_log_grid(sigmamin, sigmamax, nsigma)

    if don!=None and doff!=None:
        xlmsdata=IMP.isd2.CrossLinkData(dist_grid,omega_grid,sigma_grid,xpot,pot,don,doff,prior)
    else:
        xlmsdata=IMP.isd2.CrossLinkData(dist_grid,omega_grid,sigma_grid,xpot,pot)
    return xlmsdata

    #-------------------------------


def get_cross_link_data_from_length(length,(distmin,distmax,ndist),
                               (omegamin,omegamax,nomega),
                               (sigmamin,sigmamax,nsigma)):
    import IMP.isd2
    
    dist_grid=get_grid(distmin, distmax, ndist, False)
    omega_grid=get_log_grid(omegamin, omegamax, nomega)
    sigma_grid=get_log_grid(sigmamin, sigmamax, nsigma)

    xlmsdata=IMP.isd2.CrossLinkData(dist_grid,omega_grid,sigma_grid,length)
    return xlmsdata


def get_grid(gmin,gmax,ngrid,boundaries):
    grid=[]
    dx = ( gmax - gmin ) / float(ngrid)
    for i in range(0,ngrid+1):
        if(not boundaries and i==0): continue
        if(not boundaries and i==ngrid): continue
        grid.append( gmin + float(i) * dx )
    return grid

    #-------------------------------

def get_log_grid(gmin,gmax,ngrid):
    from math import exp, log
    grid=[]
    for i in range(0,ngrid+1):
        grid.append( gmin*exp(float(i)/ngrid*log(gmax/gmin)) )
    return grid

    #-------------------------------


def get_drmsd(prot0, prot1):
    drmsd=0.; npairs=0.;
    for i in range(0,len(prot0)-1):
        for j in range(i+1,len(prot0)):
            dist0=IMP.core.get_distance(prot0[i],prot0[j])
            dist1=IMP.core.get_distance(prot1[i],prot1[j])
            drmsd+=(dist0-dist1)**2
            npairs+=1.
    return math.sqrt(drmsd/npairs)

    #-------------------------------


def get_residue_index_and_chain_from_particle(p):
    rind=IMP.atom.Residue(IMP.atom.Atom(p).get_parent()).get_index()
    c=IMP.atom.Residue(IMP.atom.Atom(p).get_parent()).get_parent()
    cid=IMP.atom.Chain(c).get_id()
    return rind,cid

def get_particles_by_resolution(prot,resolution):
    #this function does not work with the root hierarchy, but 
    #individual proteins
    #for hier in prot.get_children():
    particles=[]
    residues=set() 

    for p in IMP.atom.get_leaves(prot):
        residues.update(IMP.atom.Fragment(p).get_residue_indexes())
       
    firstresn=min(residues)
    lastresn=max(residues)
    for nres in range(firstresn,lastresn+1):
        s=IMP.atom.Selection(prot,residue_index=nres)
        resolutions=[]

        # calculate the closest resolution for each set of particles that represent a residue 
        ps=s.get_selected_particles()

        if len(ps)>0:
            for p in ps:
                resolutions.append(IMP.pmi.Resolution(IMP.pmi.Resolution(p)).get_resolution())
            closestres=min(resolutions, key=lambda x:abs(float(x)-float(resolution)))
            
            # now we get the particle 
            for p in ps:
                if closestres==IMP.pmi.Resolution.get_resolution(IMP.pmi.Resolution(p)): 
                    if not p in particles:  
                        particles.append(p)     
        else: 
            print "get_particles_by_resolution: residue %d in molecule %s is missing" % (nres,prot.get_name())

    return list(particles)


    #-------------------------------

def select_calpha_or_residue(prot,chain,resid,ObjectName="None:",SelectResidue=False):
    #use calphas
    p=None
    s=IMP.atom.Selection(prot, chains=chain,
         residue_index=resid, atom_type=IMP.atom.AT_CA)

    ps=s.get_selected_particles()
    #check if the calpha selection is empty
    if ps:
        if len(ps)==1:
            p=ps[0]
        else:
            print ObjectName+" multiple residues selected for selection residue %s chain %s " % (resid,chain)
    else:
        #use the residue, in case of simplified representation
        s=IMP.atom.Selection(prot, chains=chain,
            residue_index=resid)
        ps=s.get_selected_particles()
        #check if the residue selection is empty
        if ps:
            if len(ps)==1:
                p=ps[0]
            else:
                print ObjectName+" multiple residues selected for selection residue %s chain %s " % (resid,chain)

        else:
            print ObjectName+" residue %s chain %s does not exist" % (resid,chain)
    return p


class map():
      def __init__(self):
          self.map={}
      
      def set_map_element(self,xvalue,yvalue):
          self.map[xvalue]=yvalue
      
      def get_map_element(self,invalue):
          n=0
          mindist=1
          for x in self.map:
              dist=(invalue-x)*(invalue-x)

              if n==0: 
                 mindist=dist
                 minx=x             
              if dist<mindist: 
                 mindist=dist
                 minx=x
              n+=1
          return self.map[minx]
          
########################
### Tools to simulate data
########################

def normal_density_function(expected_value,sigma,x):
    return 1/math.sqrt(2*math.pi)/sigma*math.exp(-(x-expected_value)**2/2/sigma/sigma)

def log_normal_density_function(expected_value,sigma,x):
    return 1/math.sqrt(2*math.pi)/sigma/x*math.exp(-(math.log(x/expected_value)**2/2/sigma/sigma))

def get_random_data_point(expected_value,ntrials,sensitivity,sigma,outlierprob,begin_end_nbins_tuple,log=False,loggrid=False):
    import random

    begin=begin_end_nbins_tuple[0]
    end=begin_end_nbins_tuple[1]
    nbins=begin_end_nbins_tuple[2]

    if not loggrid:
        fmod_grid=get_grid(begin,end,nbins,True)
    else:
        fmod_grid=get_log_grid(begin,end,nbins)


    norm=0
    cumul=[]
    cumul.append(0)

    a=[]
    for i in range(0,ntrials):
        a.append([random.random(),True])

    if sigma != 0.0:
        for j in range(1,len(fmod_grid)):
            fj=fmod_grid[j]
            fjm1=fmod_grid[j-1]
            df = fj - fjm1

            if not log:
                pj=normal_density_function(expected_value,sigma,fj)
                pjm1=normal_density_function(expected_value,sigma,fjm1)
            else:
                pj=log_normal_density_function(expected_value,sigma,fj)
                pjm1=log_normal_density_function(expected_value,sigma,fjm1)

            norm+= (pj+pjm1)/2.0*df;
            cumul.append(norm)
            #print fj, pj

        random_points=[]

        for i in range(len(cumul)):
            #print i,a, cumul[i], norm
            for aa in a:
                if (aa[0]<=cumul[i]/norm and aa[1]):
                    random_points.append(int(fmod_grid[i]/sensitivity)*sensitivity)
                    aa[1]=False

    else:
        random_points=[expected_value]*ntrials




    for i in range(len(random_points)):
        if random.random() < outlierprob:
            a=random.uniform(begin,end)
            random_points[i]=int(a/sensitivity)*sensitivity
    print random_points
    '''
    for i in range(ntrials):
      if random.random() > OUTLIERPROB_:
        r=truncnorm.rvs(0.0,1.0,expected_value,BETA_)
        if r>1.0: print r,expected_value,BETA_
      else:
        r=random.random()
      random_points.append(int(r/sensitivity)*sensitivity)
    '''

    rmean=0.; rmean2=0.
    for r in random_points:
        rmean+=r
        rmean2+=r*r

    rmean/=float(ntrials)
    rmean2/=float(ntrials)
    stddev=math.sqrt(max(rmean2-rmean*rmean,0.))
    return rmean,stddev
        
        
