from __future__ import print_function
import IMP
import IMP.test
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container
import IMP.pmi1.restraints.stereochemistry
import IMP.pmi1.restraints.basic
import IMP.pmi1.restraints.proteomics
import IMP.pmi1.restraints.crosslinking
import IMP.pmi1.restraints.em
import IMP.pmi1.representation
import IMP.pmi1.tools
import IMP.pmi1.macros

class TestPMI(IMP.test.TestCase):
    def test_em_pmi(self):
        """Test Gaussian setup and restraint in PMI1"""
        try:
            import sklearn
        except ImportError:
            self.skipTest("no sklearn package")

        outputobjects = []
        sampleobjects = []

        # setting up topology

        m = IMP.Model()
        simo = IMP.pmi1.representation.Representation(m,upperharmonic=True,disorderedlength=False)

        fastadirectory = self.get_input_file_name("mediator/")
        pdbdirectory=self.get_input_file_name("mediator/")
        gmmdirectory=self.get_input_file_name("mediator/")
        midpdb="cr_mid_fullmed10.pdb"

        # compname  hier_name    color         fastafile              fastaid          pdbname      chain    resrange      read    "BEADS"ize rigid_body super_rigid_body emnum_components emtxtfilename  emmrcfilename chain of super rigid bodies

        domains_middle= [("med4",  "med4_1",    0.10,  fastadirectory+"med4.fasta",  "med4",   pdbdirectory+midpdb,   "D",    (1,131,0),    True,       20,      1,         [19,1,2],     2,   gmmdirectory+"med4_1.txt",  gmmdirectory+"med4_1.mrc",   [0]),
                         ("med4",  "med4_2",    0.10,  fastadirectory+"med4.fasta",  "med4",   "BEADS",               None,   (132,284,0),  True,       20,      2,         [19,1,2],     0,   None,  None,   [0])]

        domains=domains_middle

        with IMP.allow_deprecated():
            bm=IMP.pmi1.macros.BuildModel1(simo)
        bm.build_model(domains)
        bm.scale_bead_radii(40,0.8)

        resdensities_middle=bm.get_density_hierarchies([t[1] for t in domains_middle])

        # randomize the initial configuration
        simo.shuffle_configuration(100)

        # defines the movers
        simo.set_rigid_bodies_max_rot(1.0)
        simo.set_floppy_bodies_max_trans(0.1)
        simo.set_rigid_bodies_max_trans(0.1)
        outputobjects.append(simo)
        sampleobjects.append(simo)

        # scoring function
        #simo.optimize_floppy_bodies(200)

        # middle module em density
        middle_mass=sum((IMP.atom.Mass(p).get_mass() for h in resdensities_middle for p in IMP.atom.get_leaves(h)))
        gemh = IMP.pmi1.restraints.em.GaussianEMRestraint(
            resdensities_middle,
            gmmdirectory+'target_gmm.txt',
            target_mass_scale=middle_mass,
            slope=0.000001,
            target_radii_scale=3.0)
        gemh.get_restraint_set().set_was_used(True)

if __name__ == '__main__':
    IMP.test.main()
