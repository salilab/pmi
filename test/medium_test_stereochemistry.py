import IMP
import IMP.atom
import IMP.test
import IMP.pmi1
import IMP.pmi1.io
import IMP.pmi1.topology
import IMP.pmi1.restraints.stereochemistry


class Tests(IMP.test.TestCase):

    def test_parse_dssp(self):
        """Test reading DSSP files"""
        sses = IMP.pmi1.io.parse_dssp(self.get_input_file_name('chainA.dssp'),'A')
        self.assertEqual(sorted(sses.keys()),sorted(['helix','beta','loop']))
        self.assertEqual(sses['helix'][1][0],[100,126,'A'])
        self.assertEqual(sses['beta'][0],[[76,78,'A'],[91,93,'A']])
        self.assertEqual(len(sses['helix']),20)
        self.assertEqual(len(sses['beta']),3)
        self.assertEqual(len(sses['loop']),32)

    def test_excluded_volume_sphere(self):
        """Test excluded volume in PMI1"""
        pdbfile = self.get_input_file_name("mini.pdb")
        fastafile = self.get_input_file_name("mini.fasta")

        comps = ["P1", "P2", "P3"]
        chains = "ABB"
        colors = [0., 0.5, 1.0]
        beadsize = 20
        fastids = IMP.pmi1.tools.get_ids_from_fasta_file(fastafile)

        m = IMP.Model()
        simo = IMP.pmi1.representation.Representation(m)

        simo.create_component(comps[0], color=colors[0])
        simo.add_component_sequence(comps[0], fastafile, id=fastids[0])
        simo.autobuild_model(comps[0], pdbfile, chains[0],
                             resolutions=[1], missingbeadsize=beadsize)

        simo.create_component(comps[1], color=colors[1])
        simo.add_component_sequence(comps[1], fastafile, id=fastids[1])
        simo.autobuild_model(comps[1], pdbfile, chains[1],
                             resolutions=[10], missingbeadsize=beadsize)

        simo.create_component(comps[2], color=colors[2])
        simo.add_component_sequence(comps[2], fastafile, id=fastids[1])
        simo.autobuild_model(comps[2], pdbfile, chains[2],
                             resolutions=[1,10], missingbeadsize=2)

        ev1 = IMP.pmi1.restraints.stereochemistry.ExcludedVolumeSphere(simo)
        ev1.add_to_model()
        print(ev1.get_output())
        sf = IMP.core.RestraintsScoringFunction(
                                   IMP.pmi1.tools.get_restraint_set(m))
        print(sf.evaluate(False))
        self.assertEqual(len(ev1.cpc.get_all_possible_indexes()),16)

        ev2 = IMP.pmi1.restraints.stereochemistry.ExcludedVolumeSphere(
            simo,
            [simo.hier_dict["P2"]],
            [simo.hier_dict["P3"]])
        ev2.add_to_model()
        print(ev2.get_output())
        sf = IMP.core.RestraintsScoringFunction(
            IMP.pmi1.tools.get_restraint_set(m))
        print(sf.evaluate(False))
        self.assertEqual(len(ev2.cpc.get_all_possible_indexes()),8)


if __name__ == '__main__':
    IMP.test.main()
