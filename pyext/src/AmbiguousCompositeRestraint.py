#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi

class AmbiguousCompositeRestraint():

    '''
    this restraint allows ambiguous crosslinking between multiple copies
    it is a variant of the SimplifiedCrossLinkMS
    '''

    def __init__(
        self,
        prot,
        hiers,
        restraints_file,
        cut_off=5.0,
        lam=1.0,
        plateau=0.0,
        label="None"):

        self.weight = 1.0
        self.prot = prot
        self.m = self.prot.get_model()     
        self.rs = IMP.RestraintSet(self.m, 'data')           
        self.label = "None"
        self.pairs = []

        self.outputlevel = "low"
        
        
        particles=[]
        for hier in hiers:
            particles += IMP.atom.get_leaves(hier)

        # fill the cross-linker pmfs
        # to accelerate the init the list listofxlinkertypes might contain only
        # yht needed crosslinks

        for line in open(restraints_file):

            tokens = line.split()
            # skip character
            if (tokens[0] == "#"):
                continue
            r1 = int(tokens[2])
            c1 = tokens[0]
            r2 = int(tokens[3])
            c2 = tokens[1]

            hrc1 = []
            for h in self.prot.get_children():
                if c1 in h.get_name():
                    hrc1.append(h)

            #hrc1 = [h for h in self.prot.get_children() if c1 in h.get_name()][0]
            # print line

            s1 = IMP.atom.Selection(hierarchies=hrc1, residue_index=r1)
            ps1 = s1.get_selected_particles()
            ps1nosym=[p for p in ps1 if IMP.pmi.Symmetric(p).get_symmetric()==0]
            
            if len(ps1) == 0:
                print "ConnectivityDistanceRestraint: WARNING> residue %d of chain %s is not there" % (r1, c1)
                continue
            #hrc2 = [h for h in self.prot.get_children() if c2 in h.get_name()][0]

            hrc2 = []
            for h in self.prot.get_children():
                if c2 in h.get_name():
                    hrc2.append(h)

            s2 = IMP.atom.Selection(hierarchies=hrc2, residue_index=r2)
            ps2 = s2.get_selected_particles()
            ps2nosym=[p for p in ps2 if IMP.pmi.Symmetric(p).get_symmetric()==0]
            
            if len(ps2) == 0:
                print "ConnectivityDistanceRestraint: WARNING> residue %d of chain %s is not there" % (r2, c2)
                continue

            ps1 = (list(set(ps1) & set(particles)))
            ps2 = (list(set(ps2) & set(particles)))
            ps1nosym = (list(set(ps1nosym) & set(particles)))
            ps2 = (list(set(ps2) & set(particles)))
            ps2nosym = (list(set(ps2nosym) & set(particles)))
 
            cr = IMP.pmi.CompositeRestraint(
                 self.m,
                 ps1nosym,
                 cut_off,
                 lam,
                 False,
                 plateau)
            cr.add_composite_particle(ps2)

            self.rs.add_restraint(cr)
            self.pairs.append((ps1nosym, hrc1, c1, r1, ps2, hrc2, c2, r2, cr))

            cr = IMP.pmi.CompositeRestraint(
                 self.m,
                 ps1,
                 cut_off,
                 lam,
                 False,
                 plateau)
            cr.add_composite_particle(ps2nosym)

            self.rs.add_restraint(cr)
            self.pairs.append((ps1, hrc1, c1, r1, ps2nosym, hrc2, c2, r2, cr))


    def plot_restraint(
        self,
        uncertainty1,
        uncertainty2,
        maxdist=50,
            npoints=10):
        import IMP.pmi.output as output

        p1 = IMP.Particle(self.m)
        p2 = IMP.Particle(self.m)
        d1 = IMP.core.XYZR.setup_particle(p1)
        d2 = IMP.core.XYZR.setup_particle(p2)
        d1.set_radius(uncertainty1)
        d2.set_radius(uncertainty2)
        s1 = IMP.atom.Selection(p1)
        s2 = IMP.atom.Selection(p2)
        sels = [s1, s2]
        strength = 1 / (uncertainty1 ** 2 + uncertainty2 ** 2)
        cr = IMP.atom.create_connectivity_restraint(
            sels,
            self.expdistance,
            strength)
        dists = []
        scores = []
        for i in range(npoints):
            d2.set_coordinates(
                IMP.algebra.Vector3D(maxdist / npoints * float(i), 0, 0))
            dists.append(IMP.core.get_distance(d1, d2))
            scores.append(cr.unprotected_evaluate(None))
        output.plot_xy_data(dists, scores)

    def set_label(self, label):
        self.label = label
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

    def set_output_level(self, level="low"):
            # this might be "low" or "high"
        self.outputlevel = level

    def set_weight(self, weight):
        self.weight = weight
        self.rs.set_weight(weight)

    def get_output(self):
        # content of the crosslink database pairs
        # self.pairs.append((p1,p2,dr,r1,c1,r2,c2))
        self.m.update()

        output = {}
        score = self.weight * self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["ConnectivityDistanceRestraint_Score_" + self.label] = str(score)
        for n, p in enumerate(self.pairs):

            ps1 = p[0]
            hrc1 = p[1]
            c1 = p[2]
            r1 = p[3]
            ps2 = p[4]
            hrc2 = p[5]
            c2 = p[6]
            r2 = p[7]
            cr = p[8]
            for n1, p1 in enumerate(ps1):
                name1 = hrc1[n1].get_name()

                for n2, p2 in enumerate(ps2):
                    name2 = hrc2[n2].get_name()
                    d1 = IMP.core.XYZR(p1)
                    d2 = IMP.core.XYZR(p2)
                    label = str(r1) + ":" + name1 + "_" + str(r2) + ":" + name2
                    output["ConnectivityDistanceRestraint_Distance_" +
                           label] = str(IMP.core.get_distance(d1, d2))

            label = str(r1) + ":" + c1 + "_" + str(r2) + ":" + c2
            output["ConnectivityDistanceRestraint_Score_" +
                   label] = str(self.weight * cr.unprotected_evaluate(None))

        return output

