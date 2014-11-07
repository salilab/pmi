import os
import IMP
import IMP.test

import IMP.pmi.restraints.stereochemistry
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.output

class Tests(IMP.test.TestCase):
    def test_pdb_writing(self):
        """Test writing of PDB files"""

        # input parameter
        pdbfile = self.get_input_file_name("mini.pdb")
        fastafile = self.get_input_file_name("mini.fasta")

        components = ["Rpb1", "Rpb2" ]

        chains = "AB"

        colors = [0., 1.0]

        beadsize = 20

        fastids = IMP.pmi.tools.get_ids_from_fasta_file(fastafile)


        m = IMP.Model()
        simo = IMP.pmi.representation.SimplifiedModel(m)


        simo.add_component_name("Rpb1", color=colors[0])
        simo.add_component_sequence("Rpb1", fastafile, id=fastids[0])
        simo.autobuild_model("Rpb1", pdbfile, "A",
                             resolutions=[1], beadsize=beadsize)
        simo.setup_component_sequence_connectivity("Rpb1", 1)

        simo.add_component_name("Rpb2", color=colors[1])
        simo.add_component_sequence("Rpb2", fastafile, id=fastids[1])
        simo.autobuild_model("Rpb2", pdbfile, "B",
                             resolutions=[10], beadsize=beadsize)
        simo.setup_component_sequence_connectivity("Rpb2", 1)

        output = IMP.pmi.output.Output()
        output.init_pdb("test_pdb_writing.pdb", simo.prot)
        output.write_pdbs()
        output.init_pdb_best_scoring("test_pdb_writing", simo.prot, 10)
        for i in range(20):
            score = -float(i)
            output.write_pdb_best_scoring(score)
        os.unlink('test_pdb_writing.pdb')
        for i in range(10):
            os.unlink('test_pdb_writing.%d.pdb' % i)
        self.assertFalse(os.path.exists('test_pdb_writing.10.pdb'))


if __name__ == '__main__':
    IMP.test.main()
