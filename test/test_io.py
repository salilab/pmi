import IMP
import IMP.test
import IMP.pmi
import IMP.pmi.io

class Tests(IMP.test.TestCase):

    def test_parse_dssp(self):
        """Test reading DSSP files"""
        sses = IMP.pmi.io.parse_dssp(self.get_input_file_name('chainA.dssp'),
                                     'A')
        self.assertEqual(sorted(sses.keys()),sorted(['helix','beta','loop']))
        self.assertEqual(sses['helix'][1][0],[100,126,'A'])
        self.assertEqual(sses['beta'][0],[[76,78,'A'],[91,93,'A']])
        self.assertEqual(len(sses['helix']),20)
        self.assertEqual(len(sses['beta']),3)
        self.assertEqual(len(sses['loop']),32)

if __name__ == '__main__':
    IMP.test.main()
