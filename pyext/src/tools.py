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

class Output():

    def __init__(self,ascii=True):
        global os,RMF,imprmf,cPickle
        import cPickle as cPickle
        import os
        try:
            import RMF
            import IMP.rmf as imprmf
            self.rmf_library=True
        except ImportError:
            self.rmf_library=False


        self.dictionary_pdbs={}
        self.dictionary_rmfs={}
        self.dictionary_stats={}
        self.dictionary_stats2={}
        self.best_score_list=None
        self.nbestscoring=None
        self.suffix=None
        self.ascii=ascii
        self.initoutput={}

    def get_pdb_names(self):
        return self.dictionary_pdbs.keys()

    def get_rmf_names(self):
        return self.dictionary_rmfs.keys()

    def get_stat_names(self):
        return self.dictionary_stats.keys()

    def init_pdb(self,name,prot):
        flpdb=open(name,'w')
        flpdb.close()
        self.dictionary_pdbs[name]=prot

    def write_pdb(self,name,appendmode=True):
        if appendmode:
            flpdb=open(name,'a')
        else:
            flpdb=open(name,'w')
        IMP.atom.write_pdb(self.dictionary_pdbs[name],flpdb)
        flpdb.close()

    def write_pdbs(self,appendmode=True):
        for pdb in self.dictionary_pdbs.keys():
            self.write_pdb(pdb,appendmode)

    def init_pdb_best_scoring(self,suffix,prot,nbestscoring):
        # save only the nbestscoring conformations
        # create as many pdbs as needed
        self.best_score_list=[]
        self.suffix=suffix
        self.nbestscoring=nbestscoring
        for i in range(self.nbestscoring):
            name=suffix+"."+str(i)+".pdb"
            flpdb=open(name,'w')
            flpdb.close()
            self.dictionary_pdbs[name]=prot

    def write_pdb_best_scoring(self,score):
        if self.nbestscoring==None:
            print "Output.write_pdb_best_scoring: init_pdb_best_scoring not run"

        #update the score list
        if len(self.best_score_list)<self.nbestscoring:
            self.best_score_list.append(score)
            self.best_score_list.sort()
            index=self.best_score_list.index(score)
            for i in range(len(self.best_score_list)-2,index-1,-1):
                oldname=self.suffix+"."+str(i)+".pdb"
                newname=self.suffix+"."+str(i+1)+".pdb"
                os.rename(oldname, newname)
            filetoadd=self.suffix+"."+str(index)+".pdb"
            self.write_pdb(filetoadd,appendmode=False)

        else:
            if score<self.best_score_list[-1]:
                self.best_score_list.append(score)
                self.best_score_list.sort()
                self.best_score_list.pop(-1)
                index=self.best_score_list.index(score)
                for i in range(len(self.best_score_list)-1,index-1,-1):
                    oldname=self.suffix+"."+str(i)+".pdb"
                    newname=self.suffix+"."+str(i+1)+".pdb"
                    os.rename(oldname, newname)
                filenametoremove=self.suffix+"."+str(self.nbestscoring)+".pdb"
                os.remove(filenametoremove)
                filetoadd=self.suffix+"."+str(index)+".pdb"
                self.write_pdb(filetoadd,appendmode=False)

    def init_rmf(self,name,prot):
        if not self.rmf_library:
            print "Output error: neet rmf library to init rmf"
            exit()

        rh = RMF.create_rmf_file(name)
        imprmf.add_hierarchy(rh, prot)
        self.dictionary_rmfs[name]=rh

    def add_restraints_to_rmf(self,name,objectlist):
        for o in objectlist:
            rs=o.get_restraint()
            imprmf.add_restraints(self.dictionary_rmfs[name],rs.get_restraints())

    def add_geometries_to_rmf(self,name,objectlist):
        for o in objectlist:
            geos=o.get_geometries()
            imprmf.add_geometries(self.dictionary_rmfs[name],geos)

    
    def add_particle_pair_from_restraints_to_rmf(self,name,objectlist):
        for o in objectlist:
            print "here"
            
            pps=o.get_particle_pairs()
            for pp in pps:
              print type(IMP.core.EdgePairGeometry(pp))
              imprmf.add_geometry(self.dictionary_rmfs[name],IMP.core.EdgePairGeometry(pp))  

    def write_rmf(self,name,nframe):
        imprmf.save_frame(self.dictionary_rmfs[name],nframe)
        self.dictionary_rmfs[name].flush()

    def close_rmf(self,name):
        del self.dictionary_rmfs[name]

    def write_rmfs(self,nframe):
        for rmf in self.dictionary_rmfs.keys():
            self.write_rmf(rmf,nframe)

    def init_stat(self,name,listofobjects):
        if self.ascii:
            flstat=open(name,'w')
            flstat.close()
        else:
            flstat=open(name,'wb')
            flstat.close()

        #check that all objects in listofobjects have a  get_output method

        for l in listofobjects:
            if not "get_output" in dir(l):
                print "Output: object ", l, " doesn't have get_output() method"
                exit()
        self.dictionary_stats[name]=listofobjects

    def set_output_entry(self,key,value):
        self.initoutput.update({key:value})

    def write_stat(self,name,appendmode=True):
        output=self.initoutput
        for obj in self.dictionary_stats[name]:
            d=obj.get_output()
            #remove all entries that begin with _ (private entries)
            dfiltered=dict((k, v) for k, v in d.iteritems() if k[0]!="_")
            output.update(dfiltered)

        if appendmode:
            writeflag='a'
        else:
            writeflag='w'

        if self.ascii:
            flstat=open(name,writeflag)
            flstat.write("%s \n" % output)
            flstat.close()
        else:
            flstat=open(name,writeflag+'b')
            cPickle.dump(output,flstat,2)
            flstat.close()

    def write_stats(self):
        for stat in self.dictionary_stats.keys():
            self.write_stat(stat)

    def get_stat(self,name):
        output={}
        for obj in self.dictionary_stats[name]:
            output.update(obj.get_output())
        return output

    def write_test(self,name,listofobjects):
        '''
        write the test:
        output=tools.Output()
        output.write_test("test_modeling11_models.rmf_45492_11Sep13_veena_imp-020713.dat",outputobjects)
        run the test:
        output=tools.Output()        
        output.test("test_modeling11_models.rmf_45492_11Sep13_veena_imp-020713.dat",outputobjects)
        '''
        flstat=open(name,'w')    
        output=self.initoutput
        for l in listofobjects:
            if not "get_output" in dir(l):
                print "Output: object ", l, " doesn't have get_output() method"
                exit()
        self.dictionary_stats[name]=listofobjects
        for obj in self.dictionary_stats[name]:
            d=obj.get_output()
            #remove all entries that begin with _ (private entries)
            dfiltered=dict((k, v) for k, v in d.iteritems() if k[0]!="_")
            output.update(dfiltered)
        output.update({"ENVIRONMENT":str(self.get_environment_variables())})
        output.update({"IMP_VERSIONS":str(self.get_versions_of_relevant_modules())})             
        flstat.write("%s \n" % output)
        flstat.close()        
        

    def test(self,name,listofobjects):
        from numpy.testing import assert_approx_equal as aae
        output=self.initoutput
        for l in listofobjects:
            if not "get_output" in dir(l):
                print "Output: object ", l, " doesn't have get_output() method"
                exit()
        for obj in listofobjects:
            output.update(obj.get_output())
        output.update({"ENVIRONMENT":str(self.get_environment_variables())})
        output.update({"IMP_VERSIONS":str(self.get_versions_of_relevant_modules())})   

        flstat=open(name,'r')  
        for l in flstat:
            test_dict=eval(l)
        for k in test_dict:
            if k in output:
               if test_dict[k]!=output[k]: print str(k)+": test failed, old value: "+str(test_dict[k])+" new value "+str(output[k])
               #aae(float(test_dict[k]),
               #    float(output[k]),7,str(k)+": test failed, old value: "+str(test_dict[k])+" new value "+str(output[k]))
            else:
               print str(k)+" from old objects (file "+str(name)+") not in new objects"
               
    def get_environment_variables(self):
        import os 
        return str(os.environ)
    
    def get_versions_of_relevant_modules(self):
        import IMP
        versions={}        
        versions["IMP_VERSION"]=IMP.kernel.get_module_version()     
        try:
           import IMP.pmi
           versions["PMI_VERSION"]=IMP.pmi.get_module_version() 
        except (ImportError):
           pass                      
        try:
           import IMP.isd2
           versions["ISD2_VERSION"]=IMP.isd2.get_module_version()             
        except (ImportError):
           pass         
        return versions

       
    

#-------------------
      
    def init_stat2(self,name,listofobjects,extralabels=None,listofsummedobjects=None):
        #this is a new stat file that should be less 
        #space greedy!
        #listofsummedobjects must be in the form [([obj1,obj2,obj3,obj4...],label)]
        #extralabels
        
        if listofsummedobjects==None: listofsummedobjects=[]
        if extralabels==None: extralabels=[]
        flstat=open(name,'w')
        output={}
        stat2_keywords={"STAT2HEADER":"STAT2HEADER"}
        stat2_keywords.update({"STAT2HEADER_ENVIRON":str(self.get_environment_variables())})
        stat2_keywords.update({"STAT2HEADER_IMP_VERSIONS":str(self.get_versions_of_relevant_modules())})        
        stat2_inverse={}
        
        for l in listofobjects:
            if not "get_output" in dir(l):
                print "Output: object ", l, " doesn't have get_output() method"
                exit()
            else:
                d=l.get_output()
                #remove all entries that begin with _ (private entries)
                dfiltered=dict((k, v) for k, v in d.iteritems() if k[0]!="_")
                output.update(dfiltered)
        
        #check for customizable entries     
        for l in listofsummedobjects:
          for t in l[0]:
            if not "get_output" in dir(t):
                print "Output: object ", t, " doesn't have get_output() method"
                exit()
            else:
                if "_TotalScore" not in t.get_output():
                  print "Output: object ", t, " doesn't have _TotalScore entry to be summed"
                  exit()
                else:
                  output.update({l[1]:0.0})            
        
        for k in extralabels:
            output.update({k:0.0})  
        
        for n,k in enumerate(output):
            stat2_keywords.update({n:k})
            stat2_inverse.update({k:n})
        
        flstat.write("%s \n" % stat2_keywords)
        flstat.close()
        self.dictionary_stats2[name]=(listofobjects,stat2_inverse,listofsummedobjects,extralabels)

    def write_stat2(self,name,appendmode=True):
        output={}
        (listofobjects,stat2_inverse,listofsummedobjects,extralabels)=self.dictionary_stats2[name]

        #writing objects
        for obj in listofobjects:
            od=obj.get_output()
            dfiltered=dict((k, v) for k, v in od.iteritems() if k[0]!="_")
            for k in dfiltered: 
               output.update({stat2_inverse[k]:od[k]})
        
        #writing summedobjects
        for l in listofsummedobjects:
           partial_score=0.0
           for t in l[0]:
             d=t.get_output()
             partial_score+=float(d["_TotalScore"])
           output.update({stat2_inverse[l[1]]:str(partial_score)})
        
        #writing extralabels
        for k in extralabels:
           if k in self.initoutput:
              output.update({stat2_inverse[k]:self.initoutput[k]})
           else:
              output.update({stat2_inverse[k]:"None"})              
        
        if appendmode:
            writeflag='a'
        else:
            writeflag='w'

        flstat=open(name,writeflag)
        flstat.write("%s \n" % output)
        flstat.close()

    def write_stats2(self):
        for stat in self.dictionary_stats2.keys():
            self.write_stat2(stat)

class ProcessOutput():
    
    def __init__(self,filename):
        self.filename=filename
        self.isstat1=False
        self.isstat2=False
        
        #open the file
        if self.filename!=None:
           f=open(self.filename,"r")
        else:
           print "Error: No file name provided. Use -h for help"
           exit()
        
        #get the keys from the first line
        for line in f.readlines():
            d=eval(line)
            self.klist=d.keys()
            #check if it is a stat2 file
            if "STAT2HEADER" in self.klist: 
                import operator
                self.isstat2=True
                for k in self.klist:
                    if "STAT2HEADER" in str(k):
                       #if print_header: print k, d[k]
                       del d[k]
                stat2_dict=d
                #get the list of keys sorted by value
                kkeys=[k[0] for k in sorted(stat2_dict.iteritems(), key=operator.itemgetter(1))]
                self.klist=[k[1] for k in sorted(stat2_dict.iteritems(), key=operator.itemgetter(1))]
                self.invstat2_dict={}
                for k in kkeys:
                    self.invstat2_dict.update({stat2_dict[k]:k})
            else:
                self.isstat1=True
                self.klist.sort()
                
            break
        f.close()    
            
    def get_keys(self):      
        return self.klist
        
    def get_fields(self,fields):
           
           
           outdict={}
           for field in fields:
               outdict[field]=[]
           
           #print fields values
           f=open(self.filename,"r")
           line_number=0
           for line in f.readlines():
              line_number+=1
              try:
                 d=eval(line)
              except:
                 print "# Warning: skipped line number " + str(line_number) + " not a valid line"
                 continue
              if   self.isstat1: [outdict[field].append(d[field]) for field in fields]
              elif self.isstat2: 
                   if line_number==1: continue
                   [outdict[field].append(d[self.invstat2_dict[field]]) for field in fields]
           f.close()
           return outdict        

    def plot_fields(self,fields):
        import matplotlib.pyplot as plt
        
        plt.rc('lines', linewidth=4)
        fig, axs  = plt.subplots(nrows=len(fields))
        
        fig.set_size_inches(10.5,5.5*len(fields))
        
        plt.rc('axes', color_cycle=['r', 'g', 'b', 'y'])
        
        n=0
        for key in fields:
           x = range(len(fields[key]))
           y=[float(y) for y in fields[key]]
           if len(fields)>1:
              axs[n].plot(x,y)
              axs[n].set_title(key)
           else:
              axs.plot(x,y)
              axs.set_title(key)
           
           n+=1

        # Tweak spacing between subplots to prevent labels from overlapping
        plt.subplots_adjust(hspace=0.3)
        plt.show()     


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
   
    #for hier in prot.get_children():

    particles=[]
    resolutions=[]
    residues=set() 

    #calculate the closest resolution        
    for p in IMP.atom.get_leaves(prot):
        res=IMP.pmi.Resolution.get_resolution(IMP.pmi.Resolution(p))
        residues.update(IMP.atom.Fragment(p).get_residue_indexes())
        resolutions.append(res)
       
    closestres=min(resolutions, key=lambda x:abs(x-resolution))
    
    
    for p in IMP.atom.get_leaves(prot):
        if closestres==IMP.pmi.Resolution.get_resolution(IMP.pmi.Resolution(p)):  
          particles.append(p)
          for rindex in IMP.atom.Fragment(p).get_residue_indexes():
              residues.remove(rindex)
    
    #select the rest, residues which were not included because
    #they were not multi-res
    
    s=IMP.atom.Selection(prot, residue_indexes=list(residues))
    particles+=s.get_selected_particles()
    
    
    return particles

def set_floppy_body(p):
    if IMP.core.RigidMember.particle_is_instance(p):
        rb=IMP.core.RigidMember(p).get_rigid_body()
        rb.set_is_rigid_member(p.get_index(),False)

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

############################
#####Analysis tools
############################


#--------------------------------
#transform rmf-to-pdb

def rmf_to_pdb(rmffilename, frame):
   import RMF
   import IMP.rmf
   fh = RMF.open_rmf_file(rmffilename)
   traverse(fh.get_root_node())
   
def traverse(n):
   
   print "#####"
   
   c = n.get_children()
   
   
   if len(c) == 0:
       print "here"
       print IMP.atom.get_pdb_string(n)
   else:
       for ch in c:
           print "there", ch
           traverse(c)



# ----------------------------------
class GetModelDensity():
    def __init__(self, prot, dens_thresh=0.1, margin=20., voxel=5.):
        global impem,deepcopy,cdist
        import IMP.em as impem
        from copy import deepcopy
        from scipy.spatial.distance import cdist

        self.prot= prot
        self.dens_thresh= dens_thresh
        self.margin= margin
        self.voxel= voxel
        self.mgr= None
        self.densities= {}

    def get_grid_termini(self):
        minx,maxx,miny,maxy,minz,maxz = inf,-inf,inf,-inf,inf,-inf
        for p in IMP.atom.get_leaves(self.prot):
            a = IMP.core.XYZ(p)
            if a.get_x()<minx: minx=a.get_x()
            if a.get_x()>maxx: maxx=a.get_x()
            if a.get_y()<miny: miny=a.get_y()
            if a.get_y()>maxy: maxy=a.get_y()
            if a.get_z()<minz: minz=a.get_z()
            if a.get_z()>maxz: maxz=a.get_z()
        minx-=self.margin
        maxx+=self.margin
        miny-=self.margin
        maxy+=self.margin
        minz-=self.margin
        maxz+=self.margin
        mgr= mgrid[minx:maxx:self.voxel,\
                   miny:maxy:self.voxel,\
                   minz:maxz:self.voxel]
        mgr= reshape(mgr, (3,-1)).T
        self.mgr= mgr
        return self.mgr

    def set_mgr(self,mgr): self.mgr = mgr

    def get_subunit_density(self,name):
        coords= []
        radii= []
        for part in [IMP.atom.get_leaves(c) for c in self.prot.get_children()\
                     if c.get_name()==name][-1]:
            p= IMP.core.XYZR(part)
            coords.append(array([p.get_x(),p.get_y(),p.get_z()]))
            radii.append(p.get_radius())
        coords= array(coords)
        radii= array(radii)
        dists= cdist(self.mgr, coords)-radii
        dens= set(list(argwhere(dists<0)[:,0]))
        return dens

    def get_subunits_densities(self):
        for subunit in self.prot.get_children():
            subname= subunit.get_name()
            dens= self.get_subunit_density(subname)
            if subname not in self.densities:
                self.densities[subname]= array([1 if i in dens else 0 for i in xrange(len(self.mgr))])
            else:
                self.densities[subname]+= array([1 if i in dens else 0 for i in xrange(len(self.mgr))])
            #print self.densities[subname],self.densities[subname].max(),subname
        return self.densities

    def update_dict(self, dendict):
        self.densities= deepcopy(dendict)

    def write_mrc(self, outname):
        for subunit in self.densities:
            mdl= IMP.Model()
            apix=self.voxel
            resolution=6.
            bbox= IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(\
                          self.mgr[:,0].min(),self.mgr[:,1].min(),self.mgr[:,2].min()),\
                          IMP.algebra.Vector3D(\
                          self.mgr[:,0].max(),self.mgr[:,1].max(),self.mgr[:,2].max()))
            dheader = impem.create_density_header(bbox,apix)
            dheader.set_resolution(resolution)

            dmap = impem.SampledDensityMap(dheader)
            ps = []
            freqs= self.densities[subunit]
            for x,i in enumerate(self.mgr):
                if freqs[x]==0.: continue
                p=IMP.Particle(mdl)
                IMP.core.XYZR.setup_particle(p,\
                                     IMP.algebra.Sphere3D(i,\
                                     1.))#freqs[x]))
                s=IMP.atom.Mass.setup_particle(p,freqs[x])
                ps.append(p)
            dmap.set_particles(ps)
            dmap.resample()
            dmap.calcRMS() # computes statistic stuff about the map and insert it in the header
            #print subunit, len(ps)
            impem.write_map(dmap,outname+"_"+subunit+".mrc",impem.MRCReaderWriter())



# ----------------------------------
class Clustering():

    def __init__(self):
        global impem,deepcopy,cdist,array,argwhere,mgrid,shape,reshape,zeros,sqrt,diagonal,argsort
        import IMP.em as impem
        from numpy import array,argwhere,mgrid,shape,reshape,zeros,diagonal,argsort
        from copy import deepcopy
        from scipy.spatial.distance import cdist
        from math import sqrt
        self.all_coords = {}

    def set_prot(self, prot):
        self.prot = prot

    def get_subunit_coords(self,frame, align=0):
        coords= []
        for part in IMP.atom.get_leaves(self.prot):
            p= IMP.core.XYZR(part)
            #coords.append(array([p.get_x(),p.get_y(),p.get_z()]))
            coords.append(p.get_coordinates())
        #coords= array(coords)
        self.all_coords[frame]= coords

    def rmsd(self,mtr1,mtr2):
        return sqrt(sum(diagonal(cdist(mtr1,mtr2)**2)) / len(mtr1))
        #return IMP.atom.get_rmsd(mtr1,mtr2)

    def set_template(self, part_coords):
        self.tmpl_coords = part_coords

    def align_and_fill(self, frame, coords, assmb):
        transformation = IMP.algebra.get_transformation_aligning_first_to_second(coords,self.tmpl_coords)
        print IMP.atom.get_rmsd(coords, self.tmpl_coords),'###',
        coords = [transformation.get_transformed(n) for n in coords]
        assmb_coords = [transformation.get_transformed(IMP.core.XYZ(n).get_coordinates()) \
                        for n in IMP.atom.get_leaves(assmb)]
        print transformation,'###', IMP.atom.get_rmsd(coords, self.tmpl_coords)
        self.all_coords[frame]= assmb_coords

    def dist_matrix(self):
        K= self.all_coords.keys()
        M = zeros((len(K), len(K)))
        for f1 in xrange(len(K)-1):
            for f2 in xrange(f1,len(K)):
                r= self.rmsd(self.all_coords[K[f1]], self.all_coords[K[f2]])
                M[f1,f2]= r
                M[f2,f1]= r

        print M.max()
        from scipy.cluster import hierarchy as hrc
        import pylab as pl
        import pickle
        C = hrc.fclusterdata(M,0.5)
        outf = open('tmp_cluster_493.pkl','w')
        pickle.dump((K,M),outf)
        outf.close()
        #exit()
        C = list(argsort(C))
        M= M[C,:][:,C]
        M[0,0]=60.
        fig = pl.figure()
        ax = fig.add_subplot(111)
        cax = ax.imshow(M, interpolation='nearest')
        ax.set_yticks(range(len(K)))
        ax.set_yticklabels( [K[i] for i in C] )
        fig.colorbar(cax)
        pl.show()

# ----------------------------------
class GetModelDensity2():
    def __init__(self, align=0, dens_thresh=0.1, margin=20., voxel=5.):
        global impem,deepcopy,cdist,array,argwhere,mgrid,shape,reshape
        import IMP.em as impem
        from numpy import array,argwhere,mgrid,shape,reshape
        from copy import deepcopy
        from scipy.spatial.distance import cdist


        self.dens_thresh= dens_thresh
        self.margin= margin
        self.voxel= voxel
        self.mgr= None
        self.align= align
        self.densities= {}

    def set_align(ali):
        self.align = ali

    def set_grid(self, part_coords):
        coords = array([array(list(j)) for j in part_coords])
        minx,maxx,miny,maxy,minz,maxz = min(coords[:,0]),max(coords[:,0]),\
                                        min(coords[:,1]),max(coords[:,1]),\
                                        min(coords[:,2]),max(coords[:,2])
        minx-=self.margin
        maxx+=self.margin
        miny-=self.margin
        maxy+=self.margin
        minz-=self.margin
        maxz+=self.margin
        grid= mgrid[minx:maxx:self.voxel,\
                   miny:maxy:self.voxel,\
                   minz:maxz:self.voxel]
        grid= reshape(grid, (3,-1)).T
        self.grid= grid
        return self.grid

    def set_template(self, part_coords):
        self.tmpl_coords = part_coords

    def align_and_add(self, coords, prot):
        transformation = IMP.algebra.get_transformation_aligning_first_to_second(coords,self.tmpl_coords)
        print IMP.atom.get_rmsd(coords, self.tmpl_coords),'###',
        coords = [transformation.get_transformed(n) for n in coords]
        print transformation,'###', IMP.atom.get_rmsd(coords, self.tmpl_coords)
        self.get_subunits_densities(prot, transformation)

    def only_add(self, coords, prot):
        transformation = ''
        self.get_subunits_densities(prot, transformation)

    def get_subunit_density(self,name, prot, transformation):
        crds= []
        radii= []
        
        for part in [IMP.atom.get_leaves(c) for c in prot.get_children()\
                     if c.get_name()==name][-1]:
            p= IMP.core.XYZR(part)
            if transformation!='': crds.append(array(list(transformation.get_transformed((p.get_x(),p.get_y(),p.get_z())))))
            else:  crds.append(array([p.get_x(),p.get_y(),p.get_z()]))
            radii.append(p.get_radius())
        '''
        for subunit in prot.get_children():
            for sbu in subunit.get_children():
                subname= sbu.get_name()
                if subname==name:
                    for part in IMP.atom.get_leaves(sbu):
                        p= IMP.core.XYZR(part)
                        crds.append(array(list(transformation.get_transformed((p.get_x(),p.get_y(),p.get_z())))))
                        radii.append(p.get_radius())
        '''
        crds= array(crds)
        radii= array(radii)
        dists= cdist(self.grid, crds)-radii
        dens= set(list(argwhere(dists<0)[:,0]))
        return dens

    def get_subunits_densities(self, prot, transformation):
        for sbu in prot.get_children():
            if 1:#for sbu in subunit.get_children():
                subname= sbu.get_name()
                dens= self.get_subunit_density(subname, prot, transformation)
                if subname not in self.densities:
                    self.densities[subname]= array([1 if i in dens else 0 for i in xrange(len(self.grid))])
                else:
                    self.densities[subname]+= array([1 if i in dens else 0 for i in xrange(len(self.grid))])
        return self.densities

    def write_mrc(self, outname):
        for subunit in self.densities:
            mdl= IMP.Model()
            apix=self.voxel
            resolution=6.
            bbox= IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(\
                          self.grid[:,0].min(),self.grid[:,1].min(),self.grid[:,2].min()),\
                          IMP.algebra.Vector3D(\
                          self.grid[:,0].max(),self.grid[:,1].max(),self.grid[:,2].max()))
            dheader = impem.create_density_header(bbox,apix)
            dheader.set_resolution(resolution)

            dmap = impem.SampledDensityMap(dheader)
            ps = []
            freqs= self.densities[subunit]
            for x,i in enumerate(self.grid):
                if freqs[x]==0.: continue
                p=IMP.Particle(mdl)
                IMP.core.XYZR.setup_particle(p,\
                                     IMP.algebra.Sphere3D(i,\
                                     1.))#freqs[x]))
                s=IMP.atom.Mass.setup_particle(p,freqs[x])
                ps.append(p)
            dmap.set_particles(ps)
            dmap.resample()
            dmap.calcRMS() # computes statistic stuff about the map and insert it in the header
            print subunit, len(ps), subunit.rsplit('.',1)[0].split('/')[-1]
            impem.write_map(dmap,outname+"_"+subunit.rsplit('.',1)[0].split('/')[-1]+".mrc",impem.MRCReaderWriter())


# ----------------------------------

class GetContactMap():
    def __init__(self, distance=15.):
        global impem,deepcopy,cdist,array,argwhere,mgrid,shape,reshape,zeros,sqrt,diagonal,argsort,log
        import IMP.em as impem
        from numpy import array,argwhere,mgrid,shape,reshape,zeros,diagonal,argsort,log
        from copy import deepcopy
        from scipy.spatial.distance import cdist
        global itemgetter
        from operator import itemgetter
        
        self.distance = distance
        self.contactmap = ''
        self.namelist = []
        self.xlinks = 0
        self.XL = {}
        self.expanded = {}
        self.resmap = {}

    def set_prot(self, prot):
        self.prot = prot

    def get_subunit_coords(self,frame, align=0):
        coords= []
        radii= []
        namelist = []
        test,testr = [],[]
        for part in self.prot.get_children():
            SortedSegments = []
            for chl in part.get_children():
                start = IMP.atom.get_leaves(chl)[0]
                end   = IMP.atom.get_leaves(chl)[-1]

                startres = IMP.atom.Fragment(start).get_residue_indexes()[0]
                endres   = IMP.atom.Fragment(end).get_residue_indexes()[-1]
                SortedSegments.append((chl,startres))
            SortedSegments = sorted(SortedSegments, key=itemgetter(1))
            for sgmnt in SortedSegments:
                for leaf in IMP.atom.get_leaves(sgmnt[0]):
                    p= IMP.core.XYZR(leaf)
                    
                    if '.pdb' not in sgmnt[0].get_name(): 
                        d = IMP.core.NonRigidMember(leaf)
                        rf = d.get_rigid_body().get_reference_frame()
                        lc = rf.get_local_coordinates(d.get_coordinates())
                        #rf.get_transformation_to()
                        d.set_internal_coordinates(lc)

                        lc = array(list(lc))
                        crd = array([p.get_x(),p.get_y(),p.get_z()])
                    else: crd = array([p.get_x(),p.get_y(),p.get_z()])

                    coords.append(crd)
                    #coords.append(array([p.get_x(),p.get_y(),p.get_z()]))
                    radii.append(p.get_radius())
                    #if part.get_name()=='tfb1':
                        #print leaf,sgmnt[0].get_name(),crd,p.get_radius()
                        #test.append(crd)
                        #testr.append(p.get_radius())
                    new_name = part.get_name()+'_'+sgmnt[0].get_name()+\
                                    '_'+str(IMP.atom.Fragment(leaf).get_residue_indexes()[0])
                    namelist.append(new_name)
                    self.expanded[new_name] = len(IMP.atom.Fragment(leaf).get_residue_indexes())
                    if part.get_name() not in self.resmap: self.resmap[part.get_name()] = {}
                    for res in IMP.atom.Fragment(leaf).get_residue_indexes():
                        self.resmap[part.get_name()][res] = new_name
        '''
        test=array(test)
        testr=array(testr)
        tc= cdist(test,test)
        print tc,'\n'
        tc= (tc-testr).T - testr
        print tc
        import pylab as pl
        fig = pl.figure()
        ax = fig.add_subplot(111)
        cax = ax.imshow(tc, interpolation='nearest')
        fig.colorbar(cax)
        pl.show()

        #exit()'''

        coords = array(coords)
        radii = array(radii)
        if len(self.namelist)==0:
            self.namelist = namelist
            self.contactmap = zeros((len(coords), len(coords)))
        distances = cdist(coords, coords)
        distances = (distances-radii).T - radii
        distances = distances<=self.distance
        self.contactmap += distances


    def add_xlinks(self, filname):
        self.xlinks = 1
        data = open(filname)
        D = data.readlines()
        data.close()

        for d in D:
            d = d.strip().split()
            t1, t2 = (d[0],d[1]), (d[1],d[0])
            if t1 not in self.XL:
                self.XL[t1] = [(int(d[2])+1, int(d[3])+1)]
                self.XL[t2] = [(int(d[3])+1, int(d[2])+1)]
            else:
                self.XL[t1].append((int(d[2])+1, int(d[3])+1))
                self.XL[t2].append((int(d[3])+1, int(d[2])+1))

        


    def dist_matrix(self, skip_cmap=0, skip_xl=1):
        K= self.namelist
        M= self.contactmap
        C,R = [],[]
        L= sum(self.expanded.values())

        # exp new
        if skip_cmap==0:
            Matrices = {}
            proteins = [p.get_name() for p in self.prot.get_children()]
            missing = []
            for p1 in xrange(len(proteins)):
                for p2 in xrange(p1,len(proteins)):
                    pl1,pl2=max(self.resmap[proteins[p1]].keys()),max(self.resmap[proteins[p2]].keys())
                    pn1,pn2=proteins[p1],proteins[p2]
                    mtr=zeros((pl1+1,pl2+1))
                    print 'Creating matrix for: ',p1,p2,pn1,pn2,mtr.shape,pl1,pl2
                    for i1 in xrange(1,pl1+1):
                        for i2 in xrange(1,pl2+1):
                            try:
                                r1=K.index( self.resmap[pn1][i1] )
                                r2=K.index( self.resmap[pn2][i2] )
                                r=M[r1,r2]
                                mtr[i1-1,i2-1]=r
                            except KeyError: missing.append((pn1,pn2,i1,i2)); pass
                    Matrices[(pn1,pn2)]=mtr

        # add cross-links
        if skip_xl==0:
            if self.XL=={}: print "Error: cross-links were not provided, use add_xlinks function!"; exit()
            Matrices_xl = {}
            proteins = [p.get_name() for p in self.prot.get_children()]
            missing_xl = []
            for p1 in xrange(len(proteins)):
                for p2 in xrange(p1,len(proteins)):
                    pl1,pl2=max(self.resmap[proteins[p1]].keys()),max(self.resmap[proteins[p2]].keys())
                    pn1,pn2=proteins[p1],proteins[p2]
                    mtr=zeros((pl1+1,pl2+1))
                    flg=0
                    try: xls = self.XL[(pn1,pn2)]
                    except KeyError:
                        try: xls = self.XL[(pn2,pn1)]; flg=1
                        except KeyError: flg=2
                    if flg==0:
                        print 'Creating matrix for: ',p1,p2,pn1,pn2,mtr.shape,pl1,pl2
                        for xl1,xl2 in xls:
                            if xl1>pl1: print 'X'*10,xl1,xl2; xl1=pl1
                            if xl2>pl2: print 'X'*10,xl1,xl2; xl2=pl2
                            mtr[xl1-1,xl2-1]=100
                    elif flg==1:
                        print 'Creating matrix for: ',p1,p2,pn1,pn2,mtr.shape,pl1,pl2
                        for xl1,xl2 in xls:
                            if xl1>pl1: print 'X'*10,xl1,xl2; xl1=pl1
                            if xl2>pl2: print 'X'*10,xl1,xl2; xl2=pl2
                            mtr[xl2-1,xl1-1]=100
                    Matrices_xl[(pn1,pn2)]=mtr                

        # expand the matrix to individual residues
        NewM = []
        for x1 in xrange(len(K)):
            lst = []
            for x2 in xrange(len(K)):
                lst += [M[x1,x2]]*self.expanded[K[x2]]
            for i in xrange(self.expanded[K[x1]]): NewM.append(array(lst))
        NewM = array(NewM)

        # make list of protein names and create coordinate lists
        r=0
        for x,k in enumerate(K):
            if x==0:
                r += self.expanded[k]
                continue
            else:
                sb = k.split('_')[0]
                sbp= K[x-1].split('_')[0]
                if sb!=sbp:
                    C.append(sbp)
                    R.append(r)
                else: r += self.expanded[k]
        r += self.expanded[k]
        C.append(sbp)
        R.append(r)
        W = []
        for x,r in enumerate(R):
            if x==0: W.append(r)
            else: W.append(r-R[x-1])

        # start plotting
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        import scipy.sparse as sparse

        f = plt.figure()
        gs = gridspec.GridSpec(len(W), len(W),
                       width_ratios=W,
                       height_ratios=W
                       )

        cnt = 0
        for x1,r1 in enumerate(R):
            if x1==0: s1=0
            else: s1 = R[x1-1]
            for x2,r2 in enumerate(R):
                if x2==0: s2=0
                else: s2 = R[x2-1]

                ax = plt.subplot(gs[cnt])
                if skip_cmap==0:
                    try: mtr = Matrices[(C[x1],C[x2])]
                    except KeyError: mtr = Matrices[(C[x2],C[x1])].T
                    #cax = ax.imshow(log(NewM[s1:r1,s2:r2] / 1.), interpolation='nearest', vmin=0., vmax=log(NewM.max()))
                    cax = ax.imshow(log(mtr), interpolation='nearest', vmin=0., vmax=log(NewM.max()))
                    ax.set_xticks([])
                    ax.set_yticks([])
                if skip_xl==0:
                    try: mtr = Matrices_xl[(C[x1],C[x2])]
                    except KeyError: mtr = Matrices_xl[(C[x2],C[x1])].T
                    cax = ax.spy(sparse.csr_matrix(mtr), markersize=10, color='white', linewidth=100, alpha=0.5)
                    ax.set_xticks([])
                    ax.set_yticks([])
                                        

                '''
                if self.xlinks==1:
                    pa, pb = C[x1], C[x2]
                    try:
                        xls = self.XL[(pa,pb)]
                        xls = [i for i in xls if i[0] in range(1,r1-s1) and i[1] in range(1,r2-s2)]
                        xla = [ia[0] for ia in xls]
                        xlb = [ib[1] for ib in xls]
                        ax.scatter(xla, xlb, s=35, facecolor='none', linewidth=2)
                    except KeyError: 
                        try:
                            xls = self.XL[(pb,p1)]
                            xls = [(i[1],i[0]) for i in xls if i[1] in range(1,r1-s1) and i[0] in range(1,r2-s2)]
                            xla = [ia[0] for ia in xls]
                            xlb = [ib[1] for ib in xls]
                            ax.scatter(xla, xlb, s=35, facecolor='none', linewidth=2)
                        except KeyError: pass
                '''
                cnt+=1
                if x2==0: ax.set_ylabel(C[x1])
        plt.show()

        '''
        fig = pl.figure()
        ax = fig.add_subplot(111)
        cax = ax.imshow(log(M), interpolation='nearest')

        ax.set_yticks( R )
        ax.set_yticklabels( C )
        fig.colorbar(cax)
        pl.show()
        '''

def get_graph_from_hierarchy(hier):
    graph=[]
    depth_dict={}
    depth=0
    (graph,depth,depth_dict)=recursive_graph(hier,graph,depth,depth_dict)
    
    #filters node labels according to depth_dict
    node_labels_dict={}
    node_size_dict={}
    for key in depth_dict:
        node_size_dict=10/depth_dict[key]
        if depth_dict[key]<3: 
           node_labels_dict[key]=key
        else:
           node_labels_dict[key]=""
    draw_graph(graph,labels_dict=node_labels_dict)

def recursive_graph(hier,graph,depth,depth_dict):
    depth=depth+1
    nameh=IMP.atom.Hierarchy(hier).get_name()
    index=str(hier.get_particle().get_index())
    name1=nameh+"|#"+index
    depth_dict[name1]=depth
    
    children=IMP.atom.Hierarchy(hier).get_children()
    
    if len(children)==1 or children==None:
       depth=depth-1
       return (graph,depth,depth_dict)
    
    else:
      for c in children:
        (graph,depth,depth_dict)=recursive_graph(c,graph,depth,depth_dict)
        nameh=IMP.atom.Hierarchy(c).get_name()
        index=str(c.get_particle().get_index())
        namec=nameh+"|#"+index        
        graph.append((name1,namec))
      
      depth=depth-1
      return (graph,depth,depth_dict)
     
def draw_graph(graph, labels_dict=None, graph_layout='spring',
               node_size=5, node_color='blue', node_alpha=0.3,
               node_text_size=11,
               edge_color='blue', edge_alpha=0.3, edge_tickness=1,
               edge_text_pos=0.3,
               text_font='sans-serif'):

    import networkx as nx
    import matplotlib.pyplot as plt


    # create networkx graph
    G=nx.Graph()

    # add edges
    for edge in graph:
        G.add_edge(edge[0], edge[1])

    # these are different layouts for the network you may try
    # shell seems to work best
    if graph_layout == 'spring':
        graph_pos=nx.spring_layout(G)
    elif graph_layout == 'spectral':
        graph_pos=nx.spectral_layout(G)
    elif graph_layout == 'random':
        graph_pos=nx.random_layout(G)
    else:
        graph_pos=nx.shell_layout(G)

    # draw graph
    nx.draw_networkx_nodes(G,graph_pos,node_size=node_size, 
                           alpha=node_alpha, node_color=node_color)
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
                           alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(G, graph_pos,labels=labels_dict,font_size=node_text_size,
                            font_family=text_font)


    #if labels is None:
    #    labels = range(len(graph))

    #edge_labels = dict(zip(graph, labels))
    #nx.draw_networkx_edge_labels(G, graph_pos, edge_labels=edge_labels, 
    #                             label_pos=edge_text_pos)

    # show graph
    plt.show()

    
    
def draw_table():
    
    #still an example!
    
    from ipyD3 import d3object
    from IPython.display import display
    
    d3 = d3object(width=800,
              height=400,
              style='JFTable',
              number=1,
              d3=None,
              title='Example table with d3js',
              desc='An example table created created with d3js with data generated with Python.')
    data=[[1277.0, 654.0, 288.0, 1976.0, 3281.0, 3089.0, 10336.0, 4650.0, 4441.0, 4670.0, 944.0, 110.0],
    [1318.0, 664.0, 418.0, 1952.0, 3581.0, 4574.0, 11457.0, 6139.0, 7078.0, 6561.0, 2354.0, 710.0],
    [1783.0, 774.0, 564.0, 1470.0, 3571.0, 3103.0, 9392.0, 5532.0, 5661.0, 4991.0, 2032.0, 680.0],
    [1301.0, 604.0, 286.0, 2152.0, 3282.0, 3369.0, 10490.0, 5406.0, 4727.0, 3428.0, 1559.0, 620.0],
    [1537.0, 1714.0, 724.0, 4824.0, 5551.0, 8096.0, 16589.0, 13650.0, 9552.0, 13709.0, 2460.0, 720.0],
    [5691.0, 2995.0, 1680.0, 11741.0, 16232.0, 14731.0, 43522.0, 32794.0, 26634.0, 31400.0, 7350.0, 3010.0],
    [1650.0, 2096.0, 60.0, 50.0, 1180.0, 5602.0, 15728.0, 6874.0, 5115.0, 3510.0, 1390.0, 170.0],
    [72.0, 60.0, 60.0, 10.0, 120.0, 172.0, 1092.0, 675.0, 408.0, 360.0, 156.0, 100.0]]
    data=[list(i) for i in zip(*data)]
    sRows=[['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'Deecember']]
    sColumns=[['Prod {0}'.format(i) for i in xrange(1,9)],
          [None, '', None, None, 'Group 1', None, None, 'Group 2']]
    d3.addSimpleTable(   data, 
                     fontSizeCells=[12,],
                     sRows=sRows,
                     sColumns=sColumns,
                     sRowsMargins=[5,50,0],
                     sColsMargins=[5,20,10],
                     spacing=0,
                     addBorders=1,
                     addOutsideBorders=-1,
                     rectWidth=45,
                     rectHeight=0                   
                 )
    html=d3.render(mode=['html', 'show'])
    display(html)

        
        
