"""@namespace IMP.pmi.restraints.crosslinking_atomic
Restraints for handling crosslinking data at atomic resolution.
"""

from __future__ import print_function
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.isd
import IMP.pmi.sampling_tools as sampling_tools
from collections import defaultdict
import os.path

def setup_nuisance(m,rs,
                   init_val,
                   min_val_nuis,
                   max_val_nuis,
                   min_val_prior,
                   max_val_prior,
                   is_opt=True,
                   add_jeff=True):
    nuisance=IMP.isd.Scale.setup_particle(IMP.Particle(m),init_val)
    nuisance.set_lower(min_val_nuis)
    nuisance.set_upper(max_val_nuis)
    nuisance.set_is_optimized(nuisance.get_nuisance_key(),is_opt)
    rs.add_restraint(IMP.isd.UniformPrior(m,nuisance,1000000000.0,
                                          max_val_prior,min_val_prior))
    if add_jeff:
        rs.add_restraint(IMP.isd.JeffreysRestraint(m,nuisance.get_particle()))
    return nuisance


class RestraintSetupError(Exception):
    pass

class MyGetRestraint(object):
    def __init__(self,rs):
        self.rs=rs
    def get_restraint_for_rmf(self):
        return self.rs

class AtomicCrossLinkMSRestraint(object):
    def __init__(self,
                 root,
                 data,
                 extra_sel={'atom_type':IMP.atom.AtomType('NZ')},
                 length=10.0,
                 slope=0.01,
                 nstates=None,
                 label='',
                 max_dist=None,
                 nuisances_are_optimized=True,
                 sigma_init=5.0,
                 psi_init = 0.01,
                 one_psi=True):
        """Experimental ATOMIC XL restraint. Provide selections for the particles to restrain.
        Automatically creates one "sigma" per crosslinked residue and one "psis" per pair.
        Other nuisance options are available.
        \note Will return an error if the data+extra_sel don't specify two particles per XL pair.
        @param root      The root hierarchy on which you'll do selection
        @param data      CrossLinkData object
        @param extra_sel  Additional selections to add to each data point. Defaults to:
                          {'atom_type':IMP.atom.AtomType('NZ')}
        @param length     The XL linker length
        @param nstates    The number of states to model. Defaults to the number of states in root.
        @param label      The output label for the restraint
        @param nuisances_are_optimized Whether to optimize nuisances
        @param sigma_init The initial value for all the sigmas
        @param psi_init   The initial value for all the psis
        """

        self.mdl = root.get_model()
        self.root = root
        self.weight = 1.0
        self.label = label
        self.length = length
        self.nuis_opt = nuisances_are_optimized
        self.nstates = nstates
        if nstates is None:
            self.nstates = len(IMP.atom.get_by_type(root,IMP.atom.STATE_TYPE))
        elif nstates!=len(IMP.atom.get_by_type(root,IMP.atom.STATE_TYPE)):
            print("Warning: nstates is not the same as the number of states in root")

        self.rs = IMP.RestraintSet(self.mdl, 'xlrestr')
        self.rs_nuis = IMP.RestraintSet(self.mdl, 'prior_nuis')
        self.particles=defaultdict(list)
        self.one_psi = one_psi

        if one_psi:
            print('creating a single psi for all XLs')
        else:
            print('creating one psi for each XL')

        #### Setup two sigmas based on promiscuity of the residue ###
        psi_min_nuis = 1e-7
        psi_max_nuis = 0.4999999
        psi_min_prior = 0.01
        psi_max_prior = 0.49
        sigma_min_nuis = 1e-7
        sigma_max_nuis = 100.1
        sigma_min_prior = 1e-3
        sigma_max_prior = 100.0

        '''
        sig_threshold=4
        self.sig_low = setup_nuisance(self.mdl,self.rs_nuis,init_val=sigma_init,min_val=1.0,
                                      max_val=100.0,is_opt=self.nuis_opt)
        self.sig_high = setup_nuisance(self.mdl,self.rs_nuis,init_val=sigma_init,min_val=1.0,
                                       max_val=100.0,is_opt=self.nuis_opt)
        '''
        self.sigma = setup_nuisance(self.mdl,self.rs_nuis,
                                    init_val=sigma_init,
                                    min_val_nuis=sigma_min_nuis,
                                    max_val_nuis=sigma_max_nuis,
                                    min_val_prior=sigma_min_prior,
                                    max_val_prior=sigma_max_prior,
                                    is_opt=self.nuis_opt,
                                    add_jeff=False)
        if one_psi:
            self.psi = setup_nuisance(self.mdl,self.rs_nuis,
                                      init_val=psi_init,
                                      min_val_nuis=psi_min_nuis,
                                      max_val_nuis=psi_max_nuis,
                                      min_val_prior=psi_min_prior,
                                      max_val_prior=psi_max_prior,
                                      is_opt=self.nuis_opt,
                                      add_jeff=True)
        else:
            self.psis={}
            for unique_id in data:
                self.psis[unique_id]=setup_nuisance(self.mdl,self.rs_nuis,
                                                    init_val=psi_init,
                                                    min_val_nuis=psi_min_nuis,
                                                    max_val_nuis=psi_max_nuis,
                                                    min_val_prior=psi_min_prior,
                                                    max_val_prior=psi_max_prior,
                                                    is_opt=self.nuis_opt,
                                                    add_jeff=True)

        ### first read ahead to get the number of XL's per residue
        #num_xls_per_res=defaultdict(int)
        #for unique_id in data:
        #    for nstate in range(self.nstates):
        #        for xl in data[unique_id]:
        #            num_xls_per_res[str(xl.r1)]+=1
        #            num_xls_per_res[str(xl.r2)]+=1

        ### now create all the XL's, using the number of restraints to guide sigmas
        xlrs=[]
        for unique_id in data:
            # create restraint for this data point
            if one_psi:
                psip = self.psi.get_particle_index()
            else:
                psip = self.psis[unique_id].get_particle_index()

            r = IMP.isd.AtomicCrossLinkMSRestraint(self.mdl,
                                                   self.length,
                                                   psip,
                                                   slope,
                                                   True)
            xlrs.append(r)
            num_contributions=0

            # add a contribution for each XL ambiguity option within each state
            for nstate in range(self.nstates):
                for xl in data[unique_id]:
                    xl_pairs = xl.get_selection(root,state_index=nstate,
                                                 **extra_sel)

                    # figure out sig1 and sig2 based on num XLs
                    '''
                    num1=num_xls_per_res[str(xl.r1)]
                    num2=num_xls_per_res[str(xl.r2)]
                    if num1<sig_threshold:
                        sig1=self.sig_low
                    else:
                        sig1=self.sig_high
                    if num2<sig_threshold:
                        sig2=self.sig_low
                    else:
                        sig2=self.sig_high
                    '''
                    sig1 = self.sigma
                    sig2 = self.sigma

                    # add each copy contribution to restraint
                    for p1,p2 in xl_pairs:
                        self.particles[nstate]+=[p1,p2]
                        if max_dist is not None:
                            dist=IMP.core.get_distance(IMP.core.XYZ(p1),IMP.core.XYZ(p2))
                            if dist>max_dist:
                                continue
                        r.add_contribution([p1.get_index(),p2.get_index()],
                                           [sig1.get_particle_index(),sig2.get_particle_index()])
                        num_contributions+=1
                if num_contributions==0:
                    raise RestraintSetupError("No contributions!")

        print('created',len(xlrs),'XL restraints')
        self.rs=IMP.isd.LogWrapper(xlrs,self.weight)

    def set_weight(self,weight):
        self.weight = weight
        self.rs.set_weight(weight)

    def set_label(self, label):
        self.label = label

    def add_to_model(self):
        self.mdl.add_restraint(self.rs)
        self.mdl.add_restraint(self.rs_nuis)

    def get_hierarchy(self):
        return self.prot

    def get_restraint_set(self):
        return self.rs

    def create_restraints_for_rmf(self):
        """ create dummy harmonic restraints for each XL but don't add to model
        Makes it easy to see each contribution to each XL in RMF
        """
        dummy_mdl=IMP.Model()
        hps = IMP.core.HarmonicDistancePairScore(self.length,1.0)
        dummy_rs=[]
        for nxl in range(self.rs.get_number_of_restraints()):
            xl=IMP.isd.AtomicCrossLinkMSRestraint.get_from(self.rs.get_restraint(nxl))
            rs = IMP.RestraintSet(dummy_mdl, 'atomic_xl_'+str(nxl))
            for ncontr in range(xl.get_number_of_contributions()):
                ps=xl.get_contribution(ncontr)
                dr = IMP.core.PairRestraint(hps,[self.mdl.get_particle(p) for p in ps],
                                            'xl%i_contr%i'%(nxl,ncontr))
                rs.add_restraint(dr)
                dummy_rs.append(MyGetRestraint(rs))
        return dummy_rs


    def get_particles(self,state_num=0):
        """ Get particles involved in the restraint """
        return self.particles[state_num]

    def get_mc_sample_objects(self,max_step_sigma,max_step_psi):
        """ HACK! Make a SampleObjects class that can be used with PMI::samplers"""
        #ps=[[self.sig_low,self.sig_high,self.psi],max_step]
        psigma=[[self.sigma],max_step_sigma]
        if self.one_psi:
            ppsi=[[self.psi],max_step_psi]
        else:
            ppsi=[[self.psis[p] for p in self.psis],max_step_psi]
        ret = [sampling_tools.SampleObjects('Nuisances',psigma),
               sampling_tools.SampleObjects('Nuisances',ppsi)]
        return ret

    def __repr__(self):
        return 'XL restraint with '+str(len(self.rs.get_restraint(0).get_number_of_restraints())) \
            + ' data points'

    def load_nuisances_from_stat_file(self,in_fn,nframe):
        """Read a stat file and load all the sigmas.
        This is potentially quite stupid.
        It's also a hack since the sigmas should be stored in the RMF file.
        Also, requires same sigma for each contribution.
        """
        import subprocess
        for nxl in range(self.rs.get_number_of_restraints()):
            xl=IMP.isd.AtomicCrossLinkMSRestraint.get_from(self.rs.get_restraint(nxl))
            sig1_val = float(subprocess.check_output(["process_output.py","-f",in_fn,
                                    "-s","AtomicXLRestraint_%i_Sig1"%nxl]).split('\n>')[1+nframe])
            sig2_val = float(subprocess.check_output(["process_output.py","-f",in_fn,
                                    "-s","AtomicXLRestraint_%i_Sig2"%nxl]).split('\n>')[1+nframe])

            for contr in range(xl.get_number_of_contributions()):
                sig1,sig2=xl.get_contribution_sigmas(contr)
                IMP.isd.Scale(self.mdl,sig1).set_scale(sig1_val)
                IMP.isd.Scale(self.mdl,sig2).set_scale(sig2_val)
        print('loaded nuisances from file')

    def plot_violations(self,out_prefix,
                        max_prob_for_violation=0.1,
                        min_dist_for_violation=1e9,
                        coarsen=False):
        """Create CMM files, one for each state, of all xinks.
        will draw in GREEN if non-violated in all states (or if only one state)
        will draw in PURPLE if non-violated only in a subset of states (draws nothing elsewhere)
        will draw in RED in ALL states if all violated
        (if only one state, you'll only see green and red)

        @param out_prefix             Output xlink files prefix
        @param max_prob_for_violation It's a violation if the probability is below this
        @param min_dist_for_violation It's a violation if the min dist is above this
        @param coarsen                Use CA positions
        """
        print('going to calculate violations and plot CMM files')
        all_dists=[]

        # prepare one output file per state
        out_fns=[]
        out_nvs=[]
        for nstate in range(self.nstates):
            outf=open(out_prefix+str(nstate)+'.cmm','w')
            outf.write('<marker_set name="xlinks_state%i"> \n' % nstate)
            out_fns.append(outf)
            out_nvs.append(0)

        # for each crosslink, evaluate probability and lowest distance for each state
        for nxl in range(self.rs.get_number_of_restraints()):
            print(nxl)
            xl=IMP.isd.AtomicCrossLinkMSRestraint.get_from(self.rs.get_restraint(nxl))
            best_contr_per_state=[[1e6,None] for i in range(self.nstates)] # [low dist, idx low contr]
            ncontr_per_state=[[] for i in range(self.nstates)]     # [list of the ncontr in each state]

            # get best contributions and probabilities for each state
            for ncontr in range(xl.get_number_of_contributions()):
                idxs = xl.get_contribution(ncontr)
                nstate = IMP.atom.get_state_index(IMP.atom.Atom(self.mdl,idxs[0]))
                ncontr_per_state[nstate].append(ncontr)
                idx1,idx2,dist=self._get_contribution_info(xl,ncontr,use_CA=coarsen)
                if dist<best_contr_per_state[nstate][0]*0.9:
                    if best_contr_per_state[nstate][1] is not None:
                        c1=IMP.atom.get_chain_id(IMP.atom.Atom(self.mdl,idx1))
                        c2=IMP.atom.get_chain_id(IMP.atom.Atom(self.mdl,idx2))
                        if c1==c2:
                            continue
                    best_contr_per_state[nstate] = [dist,ncontr]
            prob_per_state = [xl.evaluate_for_contributions(cr,None) for cr in ncontr_per_state]
            prob_global = xl.unprotected_evaluate(None)
            all_dists.append(min(best_contr_per_state)[0])

            # now check each state and see how many pass
            npass=[]
            nviol=[]
            for nstate in range(self.nstates):
                prob = prob_per_state[nstate]
                idx1,idx2,low_dist=self._get_contribution_info(xl,best_contr_per_state[nstate][1],
                                                               use_CA=coarsen)
                if prob<max_prob_for_violation or low_dist>min_dist_for_violation:
                    nviol.append(nstate)
                else:
                    npass.append(nstate)

            # special case when all pass or all fail
            all_pass=False
            all_viol=False

            if len(npass)==self.nstates:
                all_pass=True
            elif len(nviol)==self.nstates:
                all_viol=True
            print(xl,'prob:',prob_global,'prob per state: ',prob_per_state,'best dists',best_contr_per_state,all_viol)
            # finally, color based on above info
            for nstate in range(self.nstates):
                if all_pass:
                    r=0.365; g=0.933; b=0.365;
                elif all_viol:
                    r=0.980; g=0.302; b=0.247;
                else:
                    if nstate in nviol:
                        continue
                    else:
                        r=0.9; g=0.34; b=0.9;
                idx1,idx2,low_dist=self._get_contribution_info(xl,best_contr_per_state[nstate][1],
                                                               use_CA=coarsen)
                c1=IMP.core.XYZ(self.mdl,idx1).get_coordinates()
                c2=IMP.core.XYZ(self.mdl,idx2).get_coordinates()
                a1=IMP.atom.Atom(self.mdl,idx1)
                a2=IMP.atom.Atom(self.mdl,idx2)

                outf = out_fns[nstate]
                nv = out_nvs[nstate]
                outf.write('<marker id= "%d" x="%.3f" y="%.3f" z="%.3f" radius="0.8" '
                           'r="%.2f" g="%.2f" b="%.2f"/> \n' % (nv,c1[0],c1[1],c1[2],r,g,b))
                outf.write('<marker id= "%d" x="%.3f" y="%.3f" z="%.3f" radius="0.8"  '
                           'r="%.2f" g="%.2f" b="%.2f"/> \n' % (nv+1,c2[0],c2[1],c2[2],r,g,b))
                outf.write('<link id1= "%d" id2="%d" radius="0.8" '
                           'r="%.2f" g="%.2f" b="%.2f"/> \n' % (nv,nv+1,r,g,b))
                out_nvs[nstate]+=2

        for nstate in range(self.nstates):
            out_fns[nstate].write('</marker_set>\n')
            out_fns[nstate].close()

        return all_dists
    def _get_contribution_info(self,xl,ncontr,use_CA=False):
        """Return the particles at that contribution. If requested will return CA's instead"""
        idx1=xl.get_contribution(ncontr)[0]
        idx2=xl.get_contribution(ncontr)[1]
        if use_CA:
            idx1 = IMP.atom.Selection(IMP.atom.get_residue(IMP.atom.Atom(self.mdl,idx1)),
                                      atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes()[0]
            idx2 = IMP.atom.Selection(IMP.atom.get_residue(IMP.atom.Atom(self.mdl,idx2)),
                                      atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes()[0]
        dist = IMP.algebra.get_distance(IMP.core.XYZ(self.mdl,idx1).get_coordinates(),
                                        IMP.core.XYZ(self.mdl,idx2).get_coordinates())
        return idx1,idx2,dist

    def print_stats(self):
        print("XL restraint statistics\n<num> <prob> <bestdist> <sig1> <sig2> <is_viol>")
        for nxl in range(self.rs.get_number_of_restraints()):
            xl=IMP.isd.AtomicCrossLinkMSRestraint.get_from(self.rs.get_restraint(nxl))
            prob = xl.unprotected_evaluate(None)
            if prob<0.05:
                is_viol=1
            else:
                is_viol=0
            low_dist=1e6
            low_contr=None
            for contr in range(xl.get_number_of_contributions()):
                pp = xl.get_contribution(contr)
                dist = IMP.core.get_distance(IMP.core.XYZ(self.mdl,pp[0]),
                                             IMP.core.XYZ(self.mdl,pp[1]))
                if dist<low_dist:
                    low_dist=dist
                    low_contr=contr
            #dist,sig1,sig2,psi = xl.get_contribution_scores(low_contr)
            #print('%i %.4f %.4f %.4f %.4f %i'%(nxl,prob,low_dist,sig1,sig2,is_viol))
            print(low_dist)


    def get_output(self):
        self.mdl.update()
        output = {}
        score = self.weight * self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["AtomicXLRestraint" + self.label] = str(score)

        ### HACK to make it easier to see the few sigmas
        #output["AtomicXLRestraint_sig_low"] = self.sig_low.get_scale()
        #output["AtomicXLRestraint_sig_high"] = self.sig_high.get_scale()
        output["AtomicXLRestraint_sigma"] = self.sigma.get_scale()
        output["AtomicXLRestraint_priors"] = self.rs_nuis.unprotected_evaluate(None)
        if self.one_psi:
            output["AtomicXLRestraint_psi"] = self.psi.get_scale()
        ######

        # count distances above length
        bad_count=0
        for nxl in range(self.rs.get_number_of_restraints()):
            xl=IMP.isd.AtomicCrossLinkMSRestraint.get_from(self.rs.get_restraint(nxl))
            prob = xl.unprotected_evaluate(None)
            if prob<0.1:
                bad_count+=1
            low_dist=1e6
            low_contr=None
            for contr in range(xl.get_number_of_contributions()):
                pp = xl.get_contribution(contr)
                dist = IMP.core.get_distance(IMP.core.XYZ(self.mdl,pp[0]),
                                             IMP.core.XYZ(self.mdl,pp[1]))
                if dist<low_dist:
                    low_dist=dist
                    low_contr=contr

            output["AtomicXLRestraint_%i_%s"%(nxl,"Prob")]=str(prob)
            output["AtomicXLRestraint_%i_%s"%(nxl,"BestDist")]=str(low_dist)
            if not self.one_psi:
                pval = IMP.isd.Scale(self.mdl,xl.get_psi()).get_scale()
                output["AtomicXLRestraint_%i_%s"%(nxl,"psi")]=str(pval)
            #dist,sig1,sig2,psi = xl.get_contribution_scores(low_contr)
            #output["AtomicXLRestraint_%i_%s"%(nxl,"Sig1")]=str(sig1)
            #output["AtomicXLRestraint_%i_%s"%(nxl,"Sig2")]=str(sig2)

        output["AtomicXLRestraint_NumViol"] = str(bad_count)
        return output
