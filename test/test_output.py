import os
import IMP
import IMP.test
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.output

class Tests(IMP.test.TestCase):
    def test_get_particle_infos(self):
        """Test get_particle_infos_for_pdb_writing with no particles"""
        m = IMP.Model()
        simo = IMP.pmi.representation.Representation(m)
        output = IMP.pmi.output.Output()
        output.init_pdb("test_output.pdb", simo.prot)
        info, center = output.get_particle_infos_for_pdb_writing(
                                              "test_output.pdb")
        self.assertEqual(len(info), 0)
        self.assertAlmostEqual(center[0], 0., delta=1e-5)
        self.assertAlmostEqual(center[1], 0., delta=1e-5)
        self.assertAlmostEqual(center[2], 0., delta=1e-5)
        os.unlink('test_output.pdb')
        
    def test_process_output(self):
        import numpy
        """test reading stat files"""
        instat = self.get_input_file_name("./output1/stat.0.out")

        po = IMP.pmi.output.ProcessOutput(instat)

        categories = po.get_keys()
        self.assertEqual(len(categories), 25)

        criteria = [("rmf_frame_index", 5, "<")]
        self.assertEqual(len(po.return_models_satisfying_criteria(criteria)), 11)

        criteria = [("rmf_frame_index", 5, "<"), ('AtomicXLRestraint', 10.0, ">")]
        self.assertEqual(len(po.return_models_satisfying_criteria(criteria)), 4)

        vals = po.get_fields(["AtomicXLRestraint"])["AtomicXLRestraint"]
        self.assertAlmostEqual(numpy.average(numpy.array(vals).astype(numpy.float)), 10.1270600392)

if __name__ == '__main__':
    IMP.test.main()
