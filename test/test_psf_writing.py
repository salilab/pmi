import os
import IMP
import IMP.test

import IMP.pmi.restraints.stereochemistry
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.output

class Tests(IMP.test.TestCase):
    def test_pdb_writing(self):
        """Test writing of PSF files"""

        # input parameter
        pdbfile = self.get_input_file_name("mini.pdb")
        fastafile = self.get_input_file_name("mini.fasta")

        components = ["Rpb1", "Rpb2" ]

        chains = "AB"

        colors = [0., 1.0]

        beadsize = 1

        fastids = IMP.pmi.tools.get_ids_from_fasta_file(fastafile)


        m = IMP.Model()
        simo = IMP.pmi.representation.SimplifiedModel(m)


        simo.add_component_name("Rpb1", color=colors[0])
        simo.add_component_sequence("Rpb1", fastafile, id=fastids[0])
        simo.autobuild_model("Rpb1", pdbfile, "A",
                             resolutions=[1], missingbeadsize=beadsize)
        simo.setup_component_sequence_connectivity("Rpb1", 1)

        simo.add_component_name("Rpb2", color=colors[1])
        simo.add_component_sequence("Rpb2", fastafile, id=fastids[1])
        simo.autobuild_model("Rpb2", pdbfile, "B",
                             resolutions=[1], missingbeadsize=beadsize)
        simo.setup_component_sequence_connectivity("Rpb2", 1)

        output = IMP.pmi.output.Output()
        output.init_pdb("test_psf_writing.pdb", simo.prot)
        output.write_pdb("test_psf_writing.pdb")
        output.write_psf("test_psf_writing.psf","test_psf_writing.pdb")
        self.assertTrue(os.path.exists('test_psf_writing.psf'))


if __name__ == '__main__':
    IMP.test.main()
