#!/usr/bin/env python
import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

class ExcludedVolumeSphere():
    '''
    all leaves of the input hierarchies will be input in the 
    restraint. If other_hierarchies is defined, then a Bipartite container
    between "hierarchies" and "other_hierarchies" leaves is initialized
    '''

    def __init__(self, representation,
                       hierarchies=None,
                       other_hierarchies=None, 
                       resolution=None, kappa=1.0):

        self.m = representation.get_model()
        self.rs = IMP.RestraintSet(self.m, 'excluded_volume')
        self.weight = 1.0
        self.kappa = kappa

        self.label = "None"
        self.cpc=None
        
        ssps=IMP.core.SoftSpherePairScore(self.kappa)
        lsa = IMP.container.ListSingletonContainer(self.m)
        particles=IMP.pmi.tools.select(representation,resolution=resolution,hierarchies=hierarchies)
        lsa.add_particles(particles)          
           
        if other_hierarchies==None:
           self.cpc=IMP.container.ClosePairContainer(lsa,0.0,10.0)
           evr=IMP.container.PairsRestraint(ssps,self.cpc)

        else:
           other_lsa = IMP.container.ListSingletonContainer(self.m)
           other_particles=IMP.pmi.tools.select(representation,resolution=resolution,hierarchies=other_hierarchies)
           other_lsa.add_particles(particles)
           self.cpc=IMP.container.CloseBipartitePairContainer(lsa,other_lsa,0.0,10.0)
           evr=IMP.container.PairsRestraint(ssps,self.cpc)

        self.rs.add_restraint(evr)

    def add_excluded_particle_pairs(self, excluded_particle_pairs):
        # add pairs to be filtered when calculating  the score
        lpc = IMP.container.ListPairContainer(self.m)
        lpc.add_particle_pairs(excluded_particle_pairs)
        icpf = IMP.container.InContainerPairFilter(lpc)
        self.cpc.add_pair_filter(icpf)

    def set_label(self, label):
        self.label = label

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_restraint(self):
        return self.rs

    def set_weight(self, weight):
        self.weight = weight
        self.rs.set_weight(weight)

    def get_output(self):
        self.m.update()
        output = {}
        score = self.weight * self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["ExcludedVolumeSphere_" + self.label] = str(score)
        return output




class ResidueBondRestraint():
    '''
    add bond restraint between pair of consecutive 
    residues/beads to enforce the stereochemistry.
    '''
    import IMP.pmi.tools
    from math import pi as pi
    
    def __init__(self,representation,selection_tuple,distance=3.78,strength=10.0):
        self.m=representation.get_model()
        self.rs = IMP.RestraintSet(self.m, "Bonds")
        self.weight=1
        self.label="None"
        self.pairslist=[]
        
        particles=IMP.pmi.tools.select_by_tuple(representation,selection_tuple,resolution=1)
        
        ts=IMP.core.Harmonic(distance,strength)
        

        for ps in IMP.pmi.tools.sublist_iterator(particles,2,2):
            pair=[]
            if len(ps)!=2: print "ResidueBondRestraint: wrong length of pair"; exit()
            for p in ps:
                if not IMP.atom.Residue.particle_is_instance(p):
                   print "ResidueBondRestraint: not a residue"; exit()
                else:
                   pair.append(p)
            print "ResidueBondRestraint: adding a restraint between %s %s" % (pair[0].get_name(),pair[1].get_name())       
            self.rs.add_restraint(IMP.core.DistanceRestraint(ts,pair[0],pair[1]))
            self.pairslist.append(IMP.ParticlePair(pair[0], pair[1]))
            self.pairslist.append(IMP.ParticlePair(pair[1], pair[0]))

    def set_label(self, label):
        self.label = label
        self.rs.set_name(label)
        for r in self.rs.get_restraints():
            r.set_name(label)

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_restraint(self):
        return self.rs

    def set_weight(self, weight):
        self.weight = weight
        self.rs.set_weight(weight)

    def get_excluded_pairs(self):
        return self.pairslist

    def get_output(self):
        self.m.update()
        output = {}
        score = self.weight * self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["ResidueBondRestraint_" + self.label] = str(score)
        return output
    


class ResidueAngleRestraint():
    '''
    add angular restraint between triplets of consecutive 
    residues/beads to enforce the stereochemistry.
    '''
    import IMP.pmi.tools
    from math import pi as pi
    
    def __init__(self,representation,selection_tuple,anglemin=100.0,anglemax=140.0,strength=10.0):
        self.m=representation.get_model()
        self.rs = IMP.RestraintSet(self.m, "Angles")
        self.weight=1
        self.label="None"
        self.pairslist=[]

        particles=IMP.pmi.tools.select_by_tuple(representation,selection_tuple,resolution=1)
        
        ts=IMP.core.HarmonicWell((self.pi*anglemin/180.0,self.pi*anglemax/180.0),strength)
        
        for ps in IMP.pmi.tools.sublist_iterator(particles,3,3):
            triplet=[]
            if len(ps)!=3: print "ResidueAngleRestraint: wrong length of triplet"; exit()
            for p in ps:
                if not IMP.atom.Residue.particle_is_instance(p):
                   print "ResidueAngleRestraint: not a residue"; exit()
                else:
                   triplet.append(p)
            print "ResidueAngleRestraint: adding a restraint between %s %s %s" % (triplet[0].get_name(),triplet[1].get_name(),triplet[2].get_name())       
            self.rs.add_restraint(IMP.core.AngleRestraint(ts,triplet[0],triplet[1],triplet[2]))
            self.pairslist.append(IMP.ParticlePair(triplet[0], triplet[2]))
            self.pairslist.append(IMP.ParticlePair(triplet[2], triplet[0]))

    def set_label(self, label):
        self.label = label
        self.rs.set_name(label)
        for r in self.rs.get_restraints():
            r.set_name(label)

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_restraint(self):
        return self.rs

    def set_weight(self, weight):
        self.weight = weight
        self.rs.set_weight(weight)

    def get_excluded_pairs(self):
        return self.pairslist

    def get_output(self):
        self.m.update()
        output = {}
        score = self.weight * self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["ResidueAngleRestraint_" + self.label] = str(score)
        return output

class ResidueDihedralRestraint():
    '''
    add dihedral restraints between quatruplet of consecutive 
    residues/beads to enforce the stereochemistry.
    Give as input a string of "C" and "T", meaning cys (0+-40) or trans (180+-40)
    dihedral. The length of the string mush be #residue-3.
    Without the string, the dihedral will be assumed trans.
    '''
    import IMP.pmi.tools
    from math import pi as pi
    
    def __init__(self,representation,selection_tuple,stringsequence=None,strength=10.0):
        self.m=representation.get_model()
        self.rs = IMP.RestraintSet(self.m, "Angles")
        self.weight=1
        self.label="None"
        self.pairslist=[]
       
        particles=IMP.pmi.tools.select_by_tuple(representation,selection_tuple,resolution=1)
        
        if stringsequence==None:
           stringsequence="T"*(len(particles)-3)        

        for n,ps in enumerate(IMP.pmi.tools.sublist_iterator(particles,4,4)):
            quadruplet=[]
            if len(ps)!=4: print "ResidueDihedralRestraint: wrong length of quadruplet"; exit()
            for p in ps:
                if not IMP.atom.Residue.particle_is_instance(p):
                   print "ResidueDihedralRestraint: not a residue"; exit()
                else:
                   quadruplet.append(p)
            dihedraltype=stringsequence[n]
            if dihedraltype=="C":
               anglemin=-70.0
               anglemax=70.0
               ts=IMP.core.HarmonicWell((self.pi*anglemin/180.0,self.pi*anglemax/180.0),strength)
               print "ResidueDihedralRestraint: adding a CYS restraint between %s %s %s %s" % (quadruplet[0].get_name(),quadruplet[1].get_name(),
               quadruplet[2].get_name(),quadruplet[3].get_name())       
            if dihedraltype=="T":
               anglemin=180-70.0
               anglemax=180+70.0
               ts=IMP.core.HarmonicWell((self.pi*anglemin/180.0,self.pi*anglemax/180.0),strength)
               print "ResidueDihedralRestraint: adding a TRANS restraint between %s %s %s %s" % (quadruplet[0].get_name(),quadruplet[1].get_name(),
               quadruplet[2].get_name(),quadruplet[3].get_name()) 
            self.rs.add_restraint(IMP.core.DihedralRestraint(ts,quadruplet[0],quadruplet[1],quadruplet[2],quadruplet[3]))
            self.pairslist.append(IMP.ParticlePair(quadruplet[0], quadruplet[3]))
            self.pairslist.append(IMP.ParticlePair(quadruplet[3], quadruplet[0]))

    def set_label(self, label):
        self.label = label
        self.rs.set_name(label)
        for r in self.rs.get_restraints():
            r.set_name(label)

    def add_to_model(self):
        self.m.add_restraint(self.rs)

    def get_restraint(self):
        return self.rs

    def set_weight(self, weight):
        self.weight = weight
        self.rs.set_weight(weight)

    def get_excluded_pairs(self):
        return self.pairslist

    def get_output(self):
        self.m.update()
        output = {}
        score = self.weight * self.rs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score)
        output["ResidueDihedralRestraint_" + self.label] = str(score)
        return output
            

#
class SecondaryStructure():
    import IMP.isd2
    from math import pi
    from math import log    
    
    def __init__(
        self,
        representation,
        selection_tuple,
        ssstring,
        mixture=False,
        nativeness=1.0,
            kt_caff=0.1):
        # check that the secondary structure string
        # is compatible with the ssstring


        self.particles=IMP.pmi.tools.select_by_tuple(representation,selection_tuple,resolution=1)
        self.m = representation.get_model()
        self.dihe_dict = {}
        self.ang_dict = {}
        self.do_mix = {}
        self.anglfilename = IMP.isd2.get_data_path("CAAngleRestraint.dat")
        self.dihefilename = IMP.isd2.get_data_path("CADihedralRestraint.dat")
        self.nativeness = nativeness
        self.kt_caff = kt_caff
        self.anglrs = IMP.RestraintSet(self.m, "Angles")
        self.dihers = IMP.RestraintSet(self.m, "Dihedrals")
        self.bondrs = IMP.RestraintSet(self.m, "Bonds")
        self.label = "None"

        if len(self.particles) != len(ssstring):
            print len(self.particles), len(ssstring)
            print "SecondaryStructure: residue range and SS string incompatible"
        self.ssstring = ssstring

        (bondrslist, anglrslist, diherslist,
         pairslist) = self.get_CA_force_field()
        self.pairslist = pairslist

        #print anglrslist, diherslist, bondrslist, self.particles
        self.anglrs.add_restraints(anglrslist)
        self.dihers.add_restraints(diherslist)
        self.bondrs.add_restraints(bondrslist)

    def set_label(self, label):
        self.label = label

    def add_to_model(self):
        self.m.add_restraint(self.anglrs)
        self.m.add_restraint(self.dihers)
        self.m.add_restraint(self.bondrs)

    def get_CA_force_field(self):
        bondrslist = []
        anglrslist = []
        diherslist = []
        pairslist = []
        # add bonds
        for res in range(0, len(self.particles) - 1):

            ps = self.particles[res:res + 2]
            pairslist.append(IMP.ParticlePair(ps[0], ps[1]))
            pairslist.append(IMP.ParticlePair(ps[1], ps[0]))
            br = self.get_distance_restraint(ps[0], ps[1], 3.78, 416.0)
            br.set_name('Bond_restraint')
            bondrslist.append(br)
        # add dihedrals
        for res in range(0, len(self.particles) - 4):

            # if res not in dihe_dict: continue
            # get the appropriate parameters
            # get the particles
            ps = self.particles[res:res + 5]
            [phi0,
             phi1,
             score_dih] = self.read_potential_dihedral(
                self.ssstring[res:res + 4],
                True)
            pairslist.append(IMP.ParticlePair(ps[0], ps[3]))
            pairslist.append(IMP.ParticlePair(ps[3], ps[0]))
            pairslist.append(IMP.ParticlePair(ps[1], ps[4]))
            pairslist.append(IMP.ParticlePair(ps[4], ps[1]))
            dr = IMP.isd2.CADihedralRestraint(
                ps[0],
                ps[1],
                ps[2],
                ps[3],
                ps[4],
                phi0,
                phi1,
                score_dih)
            dr.set_name('Dihedral restraint')
            diherslist.append(dr)
        # add angles
        for res in range(0, len(self.particles) - 2):
            ps = self.particles[res:res + 3]
            [psi, score_ang] = self.read_potential_angle(
                self.ssstring[res:res + 2], True)
            pairslist.append(IMP.ParticlePair(ps[0], ps[2]))
            pairslist.append(IMP.ParticlePair(ps[2], ps[0]))
            dr = IMP.isd2.CAAngleRestraint(ps[0], ps[1], ps[2], psi, score_ang)
            dr.set_name('Angle restraint')
            anglrslist.append(dr)
        return (bondrslist, anglrslist, diherslist, pairslist)

    def read_potential_dihedral(self, string, mix=False):
    # read potentials for dihedral
        score_dih = []
        phi0 = []
        phi1 = []
        for i in range(0, 36):
            phi0.append(i * 10.0 / 180.0 * self.pi)
            phi1.append(i * 10.0 / 180.0 * self.pi)
            for j in range(0, 36):
                score_dih.append(0.0)
        # open file
        if not mix:
            f = open(self.dihefilename, 'r')
            for line in f.readlines():
                riga = (line.strip()).split()
                if (len(riga) == 4 and riga[0] == string):
                    ii = int(float(riga[1]) / 10.0)
                    jj = int(float(riga[2]) / 10.0)
                    score_dih[ii * len(phi0) + jj] = - \
                        self.kt_caff * self.log(float(riga[3]))
            f.close()
        if mix:
            # mix random coil and native secondary structure
            counts = []
            for i in range(0, 36):
                for j in range(0, 36):
                    counts.append(1.0)
            f = open(self.dihefilename, 'r')
            for line in f.readlines():
                riga = (line.strip()).split()
                if (len(riga) == 4 and riga[0] == string):
                    ii = int(float(riga[1]) / 10.0)
                    jj = int(float(riga[2]) / 10.0)
                    counts[ii * len(phi0) + jj] += self.nativeness * \
                        float(riga[3])
                if (len(riga) == 4 and riga[0] == "-----"):
                    ii = int(float(riga[1]) / 10.0)
                    jj = int(float(riga[2]) / 10.0)
                    counts[ii * len(phi0) + jj] += (1.0 - self.nativeness) * \
                        float(riga[3])
            f.close()
            for i in range(len(counts)):
                score_dih[i] = -self.kt_caff * self.log(counts[i])
        return [phi0, phi1, score_dih]

    def read_potential_angle(self, string, mix=False):
    # read potentials for angles
        score_ang = []
        psi = []
        for i in range(0, 180):
            psi.append(i / 180.0 * self.pi)
            score_ang.append(0.0)
        # read file
        if not mix:
            f = open(self.anglfilename, 'r')
            for line in f.readlines():
                riga = (line.strip()).split()
                if (len(riga) == 3 and riga[0] == string):
                    ii = int(riga[1])
                    score_ang[ii] = -self.kt_caff * self.log(float(riga[2]))
            f.close()
        if mix:
            # mix random coil and native secondary structure
            counts = []
            for i in range(0, 180):
                counts.append(1.0)

            f = open(self.anglfilename, 'r')
            for line in f.readlines():
                riga = (line.strip()).split()
                if (len(riga) == 3 and riga[0] == string):
                    ii = int(riga[1])
                    counts[ii] += self.nativeness * float(riga[2])
                if (len(riga) == 3 and riga[0] == "---"):
                    ii = int(riga[1])
                    counts[ii] += (1.0 - self.nativeness) * float(riga[2])
            f.close()
            for i in range(0, 180):
                score_ang[i] = -self.kt_caff * self.log(counts[i])
        return [psi, score_ang]

    def get_excluded_pairs(self):
        return self.pairslist

    def get_restraint(self):
        tmprs = IMP.RestraintSet(self.m, 'tmp')
        tmprs.add_restraint(self.anglrs)
        tmprs.add_restraint(self.dihers)
        tmprs.add_restraint(self.bondrs)
        return tmprs

    def get_distance_restraint(self, p0, p1, d0, kappa):
        h = IMP.core.Harmonic(d0, kappa)
        dps = IMP.core.DistancePairScore(h)
        pr = IMP.core.PairRestraint(dps, IMP.ParticlePair(p0, p1))
        return pr

    def get_output(self):
        output = {}
        self.m.update()
        score_angle = self.anglrs.unprotected_evaluate(None)
        score_dihers = self.dihers.unprotected_evaluate(None)
        score_bondrs = self.bondrs.unprotected_evaluate(None)
        output["_TotalScore"] = str(score_angle + score_dihers + score_bondrs)

        output["SecondaryStructure_Angles_" + self.label] = str(score_angle)
        output["SecondaryStructure_Dihedrals_" +
               self.label] = str(score_dihers)
        output["SecondaryStructure_Bonds_" + self.label] = str(score_bondrs)
        return output
