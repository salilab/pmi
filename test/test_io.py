import IMP
import IMP.test
import IMP.pmi
import IMP.pmi.io

class Tests(IMP.test.TestCase):

    def test_parse_dssp_single_chain(self):
        """Test reading DSSP files, single chain"""
        sses = IMP.pmi.io.parse_dssp(self.get_input_file_name('chainA.dssp'),
                                     'A')
        self.assertEqual(sorted(sses.keys()),sorted(['helix','beta','loop']))
        self.assertEqual(sses['helix'][1][0],[100,126,'A'])
        self.assertEqual(sses['beta'][0],[[76,78,'A'],[91,93,'A']])
        self.assertEqual(len(sses['helix']),20)
        self.assertEqual(len(sses['beta']),3)
        self.assertEqual(len(sses['loop']),32)

    def test_parse_dssp_multiple_chain(self):
        """Test reading DSSP files, single chain"""
        sses = IMP.pmi.io.parse_dssp(self.get_input_file_name('chainA.dssp'),
                                     name_map={'A':'prot1'})
        self.assertEqual(sorted(sses.keys()),sorted(['helix','beta','loop']))
        self.assertEqual(sses['helix'][1][0],[100,126,'prot1'])
        self.assertEqual(sses['helix'][-1][0], [443, 444, 'I'])
        self.assertEqual(sses['beta'][0],[[76,78,'prot1'],[91,93,'prot1']])
        self.assertEqual(len(sses['helix']),121)
        self.assertEqual(len(sses['beta']),18)
        self.assertEqual(len(sses['loop']),183)

if __name__ == '__main__':
    IMP.test.main()
