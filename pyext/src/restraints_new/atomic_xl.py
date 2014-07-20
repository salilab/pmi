#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.isd_emxl
import IMP.pmi.hierarchy_tools as hierarchy_tools
import itertools
from collections import defaultdict

def setup_nuisance(m,initialvalue,minvalue,maxvalue,isoptimized=True):
    nuisance=IMP.isd.Scale.setup_particle(IMP.Particle(m),initialvalue)
    if minvalue:
        nuisance.set_lower(minvalue)
    if maxvalue:
        nuisance.set_upper(maxvalue)
    nuisance.set_is_optimized(nuisance.get_nuisance_key(),isoptimized)
    return nuisance

class SampleObjects(object):
    """ hack class to provide things to sample for PMI::samplers """
    def __init__(self,dict_name,pack_in_dict):
        self.d={dict_name:pack_in_dict}
    def get_particles_to_sample(self):
        return self.d

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
                 nuisances_are_optimized=True):
        """Create XL restraint. Provide selections for the particles to restrain.
        Automatically creates one "sigma" per crosslinked residue and one "psi" per pair.
        Other nuisance options are available.
        \note Will return an error if the data+extra_sel don't specify two particles per XL pair.
        @param root      The root hierarchy on which you'll do selection
        @param data      Data dict. each key is a unique XL ID,
                         entries are lists of residue pairs and score. Example:

           data[1030] =
              [ { 'r1': {'molecule':'A','residue_index':5},
                  'r2': {'molecule':'B','residue_index':100},
                   'Score': 123 },
                { 'r1': {'molecule':'C','residue_index':63},
                  'r2': {'molecule':'D','residue_index':94},
                  'Score': 600 }
              ]
        since the molecule may have multiple copies, these selection arguments could return
          multiple possible contributions. will have to check that there is only ONE per copy.
        each contribution will be assigned a PSI based on the score

        @param extra_sel Additional selections to add to each data point. Defaults to:
                         {'atom_type':IMP.atom.AtomType('NZ')}
        @param length    The XL linker length
        @param nstates   The number of states to model. Defaults to the number of states in root.
        @param one_sigma Only create one sigma particle and use it for all XL
        @param one_psi   Only create one psi particle and use it for all XL
        @param label     The output label for the restraint
        @param nuisances_are_optimized
        """

        self.mdl = root.get_model()
        self.root = root
        self.weight = 1.0
        self.label = label
        self.length = length
        self.nuis_opt = nuisances_are_optimized
        if nstates is None:
            nstates = len(IMP.atom.get_by_type(root,IMP.atom.STATE_TYPE))
        elif nstates!=len(IMP.atom.get_by_type(root,IMP.atom.STATE_TYPE)):
            print "Warning: nstates is not the same as the number of states in root"

        self.rs = IMP.RestraintSet(self.mdl, 'xlrestr')
        self.particles=[]

        #### FIX THIS NUISANCE STUFF ###
        one_sigma=True
        one_psi=True
        psi_init=0.01
        psi_min=0.0
        psi_max=0.5
        sig_threshold=4
        self.sig_low = setup_nuisance(self.mdl,initialvalue=5.0,minvalue=1.0,
                                      maxvalue=100.0,isoptimized=self.nuis_opt)
        self.sig_high = setup_nuisance(self.mdl,initialvalue=10.0,minvalue=1.0,
                                       maxvalue=100.0,isoptimized=self.nuis_opt)
        psi = setup_nuisance(self.mdl,psi_init,psi_min,psi_max,self.nuis_opt)


        ### first read ahead to get the number of XL's per residue
        num_xls_per_res=defaultdict(int)
        for unique_id in data:
            for nstate in range(nstates):
                for xl in data[unique_id]:
                    num_xls_per_res[str(xl['r1'])]+=1
                    num_xls_per_res[str(xl['r2'])]+=1

        ### now create all the XL's, using the number of restraints to guide sigmas
        xlrs=[]
        for unique_id in data:

            # create restraint for this data point
            print 'creating xl with unique id',unique_id
            r = IMP.isd_emxl.AtomicCrossLinkMSRestraint(self.mdl,self.length,slope,True)
            xlrs.append(r)
            num_contributions=0

            # add a contribution for each XL ambiguity option within each state
            for nstate in range(nstates):
                print '\tstate',nstate
                # select the state
                esel = extra_sel.copy()
                esel['state_index'] = nstate

                for xl in data[unique_id]:
                    # select the particles (should grab all copies)
                    print '\t',xl
                    sel1 = IMP.atom.Selection(root,**hierarchy_tools.combine_dicts(
                        xl['r1'],esel)).get_selected_particles()
                    print '\tsel1:',sel1

                    sel2 = IMP.atom.Selection(root,**hierarchy_tools.combine_dicts(
                        xl['r2'],esel)).get_selected_particles()
                    print '\tsel2:',sel2
                    self.particles+=sel1+sel2

                    # check to make sure all particles in selections are from different copies
                    if len(sel1)==0 or len(sel2)==0:
                        raise RestraintSetupError("at least one selection is empty")

                    for s in (sel1,sel2):
                        idxs=[IMP.atom.get_copy_index(p) for p in s]
                        if len(idxs)!=len(set(idxs)):
                            raise RestraintSetupError("this XL is selecting more than one particle per copy")

                    # add each copy contribution to restraint
                    for p1,p2 in itertools.product(sel1,sel2):
                        if max_dist is not None:
                            dist=IMP.core.get_distance(IMP.core.XYZ(p1),IMP.core.XYZ(p2))
                            if dist>max_dist:
                                continue
                        num1=num_xls_per_res[str(xl['r1'])]
                        num2=num_xls_per_res[str(xl['r2'])]
                        if num1<sig_threshold:
                            sig1=self.sig_low
                            print "\tsig1 is low"
                        else:
                            sig1=self.sig_high
                            print "\tsig1 is high"
                        if num2<sig_threshold:
                            sig2=self.sig_low
                            print "\tsig2 is low"
                        else:
                            sig2=self.sig_high
                            print "\tsig2 is high"
                        r.add_contribution([p1.get_index(),p2.get_index()],
                                           [sig1.get_particle_index(),sig2.get_particle_index()],
                                           psi.get_particle_index())
                        dist,s1,s2,psv=r.get_contribution_scores(num_contributions)
                        print '\tadding contribution, init dist',dist
                        num_contributions+=1
                if num_contributions==0:
                    raise RestraintSetupError("No contributions!")
                print '\tCur XL score:',r.evaluate(False)

        #print '>>> Number of XLs per residue'
        #for key in num_xls_per_res:
        #    print key,num_xls_per_res[key]

        self.rs=IMP.isd_emxl.LogWrapper(xlrs,self.weight)

    def set_weight(self,weight):
        self.weight = weight
        self.rs.set_weight(weight)

    def set_label(self, label):
        self.label = label

    def add_to_model(self):
        self.mdl.add_restraint(self.rs)

    def get_hierarchy(self):
        return self.prot

    def get_restraint_set(self):
        return self.rs

    def enable_sampling(self):
        """ HACK! Adds necessary attributes to the selected residues for MD sampling"""
        vxkey = IMP.FloatKey('vx')
        vykey = IMP.FloatKey('vy')
        vzkey = IMP.FloatKey('vz')
        for p in self.particles:
            ps=IMP.atom.Hierarchy(p).get_parent().get_children()
            for pp in ps:
                IMP.core.XYZ(pp).set_coordinates_are_optimized(True)
                pp.add_attribute(vxkey, 0.0)
                pp.add_attribute(vykey, 0.0)
                pp.add_attribute(vzkey, 0.0)

    def get_md_sample_objects(self):
        """ HACK! Make a SampleObjects class that can be used with PMI::samplers"""
        ps=[]
        for p in self.particles:
            ps+=IMP.atom.Hierarchy(p).get_parent().get_children()
        return [SampleObjects('Floppy_Bodies_SimplifiedModel',[ps])]

    def get_mc_sample_objects(self,max_step):
        """ HACK! Make a SampleObjects class that can be used with PMI::samplers"""
        ps=[[self.sig_low,self.sig_high],max_step]
        return [SampleObjects('Nuisances',ps)]

    def __repr__(self):
        return 'XL restraint with '+str(len(self.rs.get_restraint(0).get_number_of_restraints())) \
            + ' data points'

    def get_output(self):
        self.mdl.update()
        output = {}
        score = self.weight * self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["AtomicXLRestraint" + self.label] = str(score)

        ### HACK to make it easier to see the few sigmas
        output["AtomicXLRestraint_sig_low"] = self.sig_low.get_scale()
        output["AtomicXLRestraint_sig_high"] = self.sig_high.get_scale()
        ######

        # count distances above length
        bad_count=0
        bad_av=0.0
        for nxl in range(self.rs.get_number_of_restraints()):
            xl=IMP.isd_emxl.AtomicCrossLinkMSRestraint.cast(self.rs.get_restraint(nxl))
            low_dist=1e6
            for contr in range(xl.get_number_of_contributions()):
                dist,sig1,sig2,psi = xl.get_contribution_scores(contr)
                if self.label!='':
                    output["AtomicXLRestraint_xl%i_cont%i_%s_%s"%(nxl,contr,"Dist",self.label)]=str(dist)
                    #output["AtomicXLRestraint_xl%i_cont%i_%s_%s"%(nxl,contr,"Psi",self.label)]=str(psi)
                    #output["AtomicXLRestraint_xl%i_cont%i_%s_%s"%(nxl,contr,"Sig1",self.label)]=str(sig1)
                    #output["AtomicXLRestraint_xl%i_cont%i_%s_%s"%(nxl,contr,"Sig2",self.label)]=str(sig2)
                else:
                    output["AtomicXLRestraint_xl%i_cont%i_%s"%(nxl,contr,"Dist")]=str(dist)
                    #output["AtomicXLRestraint_xl%i_cont%i_%s"%(nxl,contr,"Psi")]=str(psi)
                    #output["AtomicXLRestraint_xl%i_cont%i_%s"%(nxl,contr,"Sig1")]=str(sig1)
                    #output["AtomicXLRestraint_xl%i_cont%i_%s"%(nxl,contr,"Sig2")]=str(sig2)
                if dist<low_dist:
                    low_dist=dist
            if low_dist>self.length:
                bad_count+=1
                bad_av+=low_dist
        bad_av/=bad_count
        output["AtomicXLRestraint_NumViol"] = str(bad_count)
        output["AtomicXLRestraint_AvViolDist"] = str(bad_av)
        return output
