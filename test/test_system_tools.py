from __future__ import print_function
import IMP.test
from IMP.pmi.topology.system_tools import _get_color_for_representation


class Tests(IMP.test.TestCase):

    def test_get_color_for_representation(self):
        """Test _get_color_for_representation()"""
        class MockRepresentation(object):
            def __init__(self, c):
                self.color = c

        c = _get_color_for_representation(MockRepresentation(0.1))
        self.assertAlmostEqual(c.get_red(), 0.0, delta=0.01)
        self.assertAlmostEqual(c.get_green(), 0.2, delta=0.01)
        self.assertAlmostEqual(c.get_blue(), 0.8, delta=0.01)

        c = _get_color_for_representation(MockRepresentation("red"))
        self.assertAlmostEqual(c.get_red(), 1.0, delta=0.01)
        self.assertAlmostEqual(c.get_green(), 0.0, delta=0.01)
        self.assertAlmostEqual(c.get_blue(), 0.0, delta=0.01)

        c = _get_color_for_representation(MockRepresentation("#fa8072"))
        self.assertAlmostEqual(c.get_red(), 0.980, delta=0.01)
        self.assertAlmostEqual(c.get_green(), 0.502, delta=0.01)
        self.assertAlmostEqual(c.get_blue(), 0.447, delta=0.01)

        c = _get_color_for_representation(MockRepresentation((0.1, 0.2, 0.3)))
        self.assertAlmostEqual(c.get_red(), 0.1, delta=0.01)
        self.assertAlmostEqual(c.get_green(), 0.2, delta=0.01)
        self.assertAlmostEqual(c.get_blue(), 0.3, delta=0.01)

        d = IMP.display.Color(0.2, 0.3, 0.4)
        c = _get_color_for_representation(MockRepresentation(d))
        self.assertAlmostEqual(c.get_red(), 0.2, delta=0.01)
        self.assertAlmostEqual(c.get_green(), 0.3, delta=0.01)
        self.assertAlmostEqual(c.get_blue(), 0.4, delta=0.01)

        self.assertRaises(TypeError,
                          _get_color_for_representation, MockRepresentation(42))

        r = MockRepresentation(None)
        self.assertIsNone(_get_color_for_representation(r))

    def test_get_structure_pdb(self):
        """Test get_structure given a single-model PDB"""
        m = IMP.Model()
        rs = IMP.pmi.topology.system_tools.get_structure(
                m, self.get_input_file_name('mini.pdb'), 'A', [4, 5])
        self.assertEqual(len(rs), 2)

        # Empty residue range
        rs = IMP.pmi.topology.system_tools.get_structure(
                m, self.get_input_file_name('mini.pdb'), 'A')
        self.assertEqual(len(rs), 7)

        rs = IMP.pmi.topology.system_tools.get_structure(
                m, self.get_input_file_name('mini.pdb'), 'A', [6, 'END'])
        self.assertEqual(len(rs), 3)

        rs = IMP.pmi.topology.system_tools.get_structure(
                m, self.get_input_file_name('mini.pdb'), 'A', [4, 5],
                ca_only=True)
        self.assertEqual(len(rs), 2)

        # Invalid range
        with self.assertWarns(IMP.pmi.StructureWarning):
            rs = IMP.pmi.topology.system_tools.get_structure(
                   m, self.get_input_file_name('mini.pdb'), 'A', [40, 50])
            self.assertEqual(len(rs), 0)

    def test_get_structure_mmcif(self):
        """Test get_structure given a single-model mmCIF"""
        m = IMP.Model()
        rs = IMP.pmi.topology.system_tools.get_structure(
                m, self.get_input_file_name('mini.cif'), 'A', [4, 5])
        self.assertEqual(len(rs), 2)

        rs = IMP.pmi.topology.system_tools.get_structure(
                m, self.get_input_file_name('mini.cif'), 'A', [4, 5],
                ca_only=True)
        self.assertEqual(len(rs), 2)

        with self.assertWarns(IMP.pmi.StructureWarning):
            rs = IMP.pmi.topology.system_tools.get_structure(
                   m, self.get_input_file_name('mini.cif'), 'A', [40, 50])
            self.assertEqual(len(rs), 0)

    def test_get_structure_multi_pdb(self):
        """Test get_structure given a multi-model PDB"""
        m = IMP.Model()
        rs = IMP.pmi.topology.system_tools.get_structure(
                m, self.get_input_file_name('multi.pdb'), 'A', [56, 57],
                model_num=1)
        self.assertEqual(len(rs), 2)

        self.assertRaises(IndexError,
                IMP.pmi.topology.system_tools.get_structure,
                m, self.get_input_file_name('multi.pdb'), 'A', [56, 57],
                model_num=2)

if __name__ == '__main__':
    IMP.test.main()
