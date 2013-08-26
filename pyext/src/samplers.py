import IMP
import IMP.core

class MonteCarlo():
    def __init__(self,m,objects,temp):
        
        #check that isd2 is installed
        try:
            global impisd2
            import IMP.isd2 as impisd2 
            self.isd2_available = True
        except ImportError:
            self.isd2_available = False


        #check that the objects containts get_particles_to_sample methods
        #and the particle type is supported
        #list of particles to sample:
        self.losp=["Rigid_Bodies","Floppy_Bodies","Nuisances","X_coord","Weights"]
        self.simulated_annealing=False
        self.selfadaptive=False
        #that is -1 because mc has not yet run
        self.nframe=-1
        self.temp=temp
        self.mvs=[]
        self.mvslabels=[]
        self.label="None"
        self.m=m


        for ob in objects:
            try:
                ob.get_particles_to_sample()
            except:
                print "MonteCarlo: object ", ob, " doesn't have get_particles_to_sample() method"

            pts=ob.get_particles_to_sample()
            for k in pts.keys():

                if "Rigid_Bodies" in k:
                    mvs=self.get_rigid_body_movers(pts[k][0],pts[k][1],pts[k][2])
                    for mv in mvs: mv.set_name(k)
                    self.mvs+=mvs


                if "Floppy_Bodies" in k:
                    mvs=self.get_floppy_body_movers(pts[k][0],pts[k][1])
                    for mv in mvs: mv.set_name(k)
                    self.mvs+=mvs

                if "X_coord" in k:
                    mvs=self.get_X_movers(pts[k][0],pts[k][1])
                    for mv in mvs: mv.set_name(k)
                    self.mvs+=mvs

                if "Nuisances" in k:
                    if not self.isd2_available:
                        print "MonteCarlo: isd2 module needed to use nuisances"
                        exit()
                    mvs=self.get_nuisance_movers(pts[k][0],pts[k][1])
                    for mv in mvs: mv.set_name(k)
                    self.mvs+=mvs

                if "Weights" in k:
                    if not self.isd2_available:
                        print "MonteCarlo: isd2 module needed to use weights"
                        exit()
                    mvs=self.get_weight_movers(pts[k][0],pts[k][1])
                    for mv in mvs: mv.set_name(k)
                    self.mvs+=mvs

        # SerialMover
        self.smv=IMP.core.SerialMover(self.mvs)

        self.mc=IMP.core.MonteCarlo(self.m)
        self.mc.set_return_best(False)
        self.mc.set_kt(self.temp)
        self.mc.add_mover(self.smv)

    def get_mc(self):
        return self.mc

    def set_scoring_function(self,objectlist):
        rs = IMP.RestraintSet(self.m,1.0,'sfo')
        for ob in  objectlist:
            rs.add_restraint(ob.get_restraint())
        sf=IMP.core.RestraintsScoringFunction([rs])
        self.mc.set_scoring_function(sf)

    def set_simulated_annealing(self,tempmin,tempmax,timemin,timemax):
        self.simulated_annealing=True
        self.tempmin=tempmin
        self.tempmax=tempmax
        self.timemin=timemin
        self.timemax=timemax

    def set_self_adaptive(self,isselfadaptive=True):
        self.selfadaptive=isselfadaptive

    def get_nuisance_movers_parameters(self):
        '''returns a dictionary with the mover parameters
        for nuisance parameters'''
        output={}
        for i in range(self.smv.get_number_of_movers()):
            mv=self.smv.get_mover(i)
            name=mv.get_name()
            if "Nuisances" in name:
                stepsize = IMP.core.NormalMover.get_from(mv).get_sigma()
                output[name]=stepsize
        return output

    def get_number_of_movers(self):
        return self.smv.get_number_of_movers()

    def get_particle_types():
        return self.losp

    def run(self,nstep):
        self.nframe+=1
        self.mc.optimize(nstep*self.get_number_of_movers())

        #apply simulated annealing protocol
        if self.simulated_annealing:
            self.temp=self.temp_simulated_annealing()
            self.mc.set_kt(self.temp)

        #apply self adaptive protocol
        if self.selfadaptive:
            for i,mv in enumerate(self.smv.get_movers()):
                name=mv.get_name()


                if "Nuisances" in name:
                    mvacc=mv.get_number_of_accepted()
                    mvprp=mv.get_number_of_proposed()
                    accept=float(mvacc)/float(mvprp)
                    nmv=IMP.core.NormalMover.get_from(mv)
                    stepsize = nmv.get_sigma()

                    if 0.4 > accept or accept > 0.6:
                        nmv.set_sigma(stepsize*2*accept)
                    if accept < 0.05:
                        accept = 0.05
                        nmv.set_sigma(stepsize*2*accept)
                    if accept > 1.0:
                        accept = 1.0
                        nmv.set_sigma(stepsize*2*accept)

                if "Weights" in name:

                    mvacc=mv.get_number_of_accepted()
                    mvprp=mv.get_number_of_proposed()
                    accept=float(mvacc)/float(mvprp)
                    wmv=impisd2.WeightMover.get_from(mv)
                    stepsize = wmv.get_radius()

                    if 0.4 > accept or accept > 0.6:
                        wmv.set_radius(stepsize*2*accept)
                    if accept < 0.05:
                        accept = 0.05
                        wmv.set_radius(stepsize*2*accept)
                    if accept > 1.0:
                        accept = 1.0
                        wmv.set_radius(stepsize*2*accept)


    def get_nuisance_movers(self,nuisances,maxstep):
        mvs=[]
        for nuisance in nuisances:
            mvs.append(IMP.core.NormalMover([nuisance],IMP.FloatKeys([IMP.FloatKey("nuisance")]),maxstep))
        return mvs

    def get_rigid_body_movers(self,rbs,maxtrans,maxrot):
        mvs=[]
        for rb in rbs:
            mvs.append(IMP.core.RigidBodyMover(rb,maxtrans,maxrot))
        return mvs

    def get_floppy_body_movers(self,fbs,maxtrans):
        mvs=[]
        for fb in fbs:
            #check is that is a rigid body member:
            if IMP.core.NonRigidMember.particle_is_instance(fb):
            #if so force the particles to move anyway
                floatkeys=[IMP.FloatKey(4),IMP.FloatKey(5),IMP.FloatKey(6)]
                for fk in floatkeys:
                    fb.set_is_optimized(fk,True)
                mvs.append(IMP.core.BallMover([fb],IMP.FloatKeys(floatkeys),maxtrans))
            else:
                #otherwise use the normal ball mover
                mvs.append(IMP.core.BallMover([fb],maxtrans))
        return mvs

    def get_X_movers(self,fbs,maxtrans):
        mvs=[]
        Xfloatkey=IMP.core.XYZ.get_xyz_keys()[0]
        for fb in fbs:
            #check is that is a rigid body member:
            if IMP.core.NonRigidMember.particle_is_instance(fb):
                print "particle is part of a rigid body"
                exit()
            else:
                #otherwise use the normal ball mover
                mvs.append(IMP.core.NormalMover([fb],[Xfloatkey],maxtrans))
        return mvs



    def get_weight_movers(self,weights,maxstep):
        mvs=[]
        for weight in weights:
            if(weight.get_number_of_states()>1): mvs.append(impisd2.WeightMover(weight,maxstep))
        return mvs

    def temp_simulated_annealing(self):
        if self.nframe%(self.timemin+self.timemax)< self.timemin:
            value=0.0
        else:
            value=1.0
        temp=self.tempmin+(self.tempmax-self.tempmin)*value
        return temp

    def set_label(self,label):
        self.label=label

    def get_frame_number(self):
        return self.nframe

    def get_output(self):
        output={}
        acceptances=[]
        for i,mv in enumerate(self.smv.get_movers()):
            mvname=mv.get_name()
            mvacc=mv.get_number_of_accepted()
            mvprp=mv.get_number_of_proposed()
            try:
                mvacr=float(mvacc)/float(mvprp)
            except:
                mvacr=0.0
            output["MonteCarlo_Acceptance_"+mvname+"_"+str(i)]=str(mvacr)
            if "Nuisances" in mvname:
                output["MonteCarlo_StepSize_"+mvname+"_"+str(i)]=str(IMP.core.NormalMover.get_from(mv).get_sigma())
            if "Weights" in mvname:
                output["MonteCarlo_StepSize_"+mvname+"_"+str(i)]=str(impisd2.WeightMover.get_from(mv).get_radius())
        output["MonteCarlo_Temperature"]=str(self.mc.get_kt())
        output["MonteCarlo_Nframe"]=str(self.nframe)
        return output






class ConjugateGradients():
    def __init__(self,m,objects):
        self.m=m
        self.nframe=-1
        self.cg=IMP.core.ConjugateGradients(self.m)

    def set_label(self,label):
        self.label=label

    def get_frame_number(self):
        return self.nframe

    def run(self,nstep):
        self.nframe+=1
        self.cg.optimize(nstep)

    def set_scoring_function(self,objectlist):
        rs = IMP.RestraintSet(self.m,1.0,'sfo')
        for ob in  objectlist:
            rs.add_restraint(ob.get_restraint())
        sf=IMP.core.RestraintsScoringFunction([rs])
        self.cg.set_scoring_function(sf)

    def get_output(self):
        output={}
        acceptances=[]
        output["ConjugatedGradients_Nframe"]=str(self.nframe)
        return output







class PyMCMover():
    #only works if the sampled particles are rigid bodies
    def __init__(self, representation, mcchild, n_mc_steps):

        #mcchild must be pmi.samplers.MonteCarlo
        #representation must be pmi.representation

        self.rbs = representation.get_rigid_bodies()

        self.mc = mcchild
        self.n_mc_steps = n_mc_steps

    def store_move(self):
        #get all xyz coordinates of all rigid bodies of all copies
        self.oldcoords=[]
        for copy in self.rbs:
            crd=[]
            for rb in copy:
                crd.append(rb.get_reference_frame())
            self.oldcoords.append(crd)

    def propose_move(self, prob):
        self.mc.run(self.n_mc_steps)

    def reset_move(self):
        #reset each copy to crd
        for copy,crd in zip(self.rbs,self.oldcoords):
            for rb,ref in zip(copy,crd):
                rb.set_reference_frame(ref)

    def get_number_of_steps(self):
        return self.n_mc_steps

    def set_number_of_steps(self, nsteps):
        self.n_mc_steps = nsteps


class PyMC():

    def __init__(self,model):
        from math import exp
        import random

        self.m=model
        self.restraints=None
        self.first_call=True
        self.nframe=-1

    def add_mover(self,mv):
        self.mv = mv

    def set_kt(self,kT):
        self.kT=kT

    def set_return_best(self,thing):
        pass

    def set_move_probability(self,thing):
        pass

    def get_energy(self):
        if self.restraints:
            pot = sum([r.evaluate(False) for r in self.restraints])
        else:
            pot=self.m.evaluate(False)
        return pot

    def metropolis(self, old, new):
        deltaE=new-old
        print ": old %f new %f deltaE %f new_epot: %f" % (old,new,deltaE,
                self.m.evaluate(False)),
        kT=self.kT
        if deltaE < 0:
            return True
        else:
            return exp(-deltaE/kT) > random.uniform(0,1)

    def optimize(self,nsteps):
        self.naccept = 0
        self.nframe+=1
        print "=== new MC call"
        #store initial coordinates
        if self.first_call:
            self.mv.store_move()
            self.first_call=False
        for i in xrange(nsteps):
            print "MC step %d " % i,
            #get total energy
            old=self.get_energy()
            #make a MD move
            self.mv.propose_move(1)
            #get new total energy
            new=self.get_energy()
            if self.metropolis(old,new):
                #move was accepted: store new conformation
                self.mv.store_move()
                self.naccept += 1
                print "accepted "
            else:
                #move rejected: restore old conformation
                self.mv.reset_move()
                print " "

    def get_number_of_forward_steps(self):
        return self.naccept

    def set_restraints(self, restraints):
        self.restraints=restraints

    def set_scoring_function(self,objects):
        #objects should be pmi.restraints
        rs = IMP.RestraintSet(self.m,1.0,'sfo')
        for ob in  objects:
            rs.add_restraint(ob.get_restraint())
        self.set_restraints([rs])

    def get_output(self):
        output={}
        acceptances=[]
        output["PyMC_Temperature"]=str(self.kT)
        output["PyMC_Nframe"]=str(self.nframe)
        return output
