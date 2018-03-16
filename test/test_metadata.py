from __future__ import print_function
import IMP.test
import IMP.pmi.metadata
import IMP.pmi.representation
import ihm.location
import os

class Tests(IMP.test.TestCase):

    def test_software(self):
        """Test metadata.Software"""
        s = IMP.pmi.metadata.Software(name='test', classification='test code',
                                      description='Some test program',
                                      url='http://salilab.org')
        self.assertEqual(s.name, 'test')

    def test_citation(self):
        """Test metadata.Citation"""
        s = IMP.pmi.metadata.Citation(title='Test paper', journal='J Mol Biol',
                                      volume=45, page_range=(1,20), year=2016,
                                      authors=['Smith, A.', 'Jones, B.'],
                                      doi='10.2345/S1384107697000225',
                                      pmid='1234')
        self.assertEqual(s.title, 'Test paper')

    def test_python_script(self):
        """Test metadata.PythonScript"""
        r = ihm.location.Repository(doi='10.5281/zenodo.46266')
        f = ihm.location.WorkflowFileLocation(repo=r, path='foo')
        s = IMP.pmi.metadata.PythonScript(location=f)

    def test_chimerax_command_script(self):
        """Test metadata.ChimeraXCommandScript"""
        r = ihm.location.Repository(doi='10.5281/zenodo.46266')
        f = ihm.location.WorkflowFileLocation(repo=r, path='foo')
        s = IMP.pmi.metadata.ChimeraXCommandScript(location=f)

    def test_repr_add(self):
        """Test Representation.add_metadata()"""
        m = IMP.Model()
        r = IMP.pmi.representation.Representation(m)
        r.add_metadata(ihm.location.Repository(
                                   doi='10.5281/zenodo.46266', root='..'))
        self.assertEqual(r._metadata[0]._root, os.path.abspath('..'))

    def test_dataset_add_parent(self):
        """Test Dataset.add_parent()"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d1 = IMP.pmi.metadata.CXMSDataset(loc)
        loc = ihm.location.InputFileLocation(repo='mydoi', path='b')
        d2 = IMP.pmi.metadata.MassSpecDataset(loc)
        d1.add_parent(d2)
        self.assertEqual(d1._parents, {d2:None})
        # Ignore duplicates
        d1.add_parent(d2)
        self.assertEqual(d1._parents, {d2:None})

    def test_dataset_add_primary_no_parents(self):
        """Test Dataset.add_primary() with no parents"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d1 = IMP.pmi.metadata.CXMSDataset(loc)
        loc = ihm.location.InputFileLocation(repo='mydoi', path='b')
        d2 = IMP.pmi.metadata.MassSpecDataset(loc)
        d1.add_primary(d2)
        self.assertEqual(d1._parents, {d2:None})

    def test_dataset_add_primary_one_parent(self):
        """Test Dataset.add_primary() with one parent"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d1 = IMP.pmi.metadata.CXMSDataset(loc)
        loc = ihm.location.InputFileLocation(repo='mydoi', path='b')
        d2 = IMP.pmi.metadata.MassSpecDataset(loc)
        d1.add_parent(d2)
        loc = ihm.location.InputFileLocation(repo='mydoi', path='c')
        d3 = IMP.pmi.metadata.MassSpecDataset(loc)
        d1.add_primary(d3)
        self.assertEqual(d1._parents, {d2:None})
        self.assertEqual(d2._parents, {d3:None})

    def test_dataset_add_primary_two_parents(self):
        """Test Dataset.add_primary() with two parents"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d1 = IMP.pmi.metadata.CXMSDataset(loc)
        loc = ihm.location.InputFileLocation(repo='mydoi', path='b')
        d2 = IMP.pmi.metadata.MassSpecDataset(loc)
        d1.add_parent(d2)
        loc = ihm.location.InputFileLocation(repo='mydoi', path='c')
        d3 = IMP.pmi.metadata.MassSpecDataset(loc)
        d1.add_parent(d3)
        loc = ihm.location.InputFileLocation(repo='mydoi', path='d')
        d4 = IMP.pmi.metadata.MassSpecDataset(loc)
        self.assertRaises(ValueError, d1.add_primary, d4)

    def test_cxms_dataset(self):
        """Test CXMSDataset"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d = IMP.pmi.metadata.CXMSDataset(loc)
        self.assertEqual(d._data_type, 'CX-MS data')

    def test_mass_spec_dataset(self):
        """Test MassSpecDataset"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d = IMP.pmi.metadata.MassSpecDataset(loc)
        self.assertEqual(d._data_type, 'Mass Spectrometry data')

    def test_em_density_dataset(self):
        """Test EMDensityDataset"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d = IMP.pmi.metadata.EMDensityDataset(loc)
        self.assertEqual(d._data_type, '3DEM volume')

    def test_pdb_dataset(self):
        """Test PDBDataset"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d = IMP.pmi.metadata.PDBDataset(loc)
        self.assertEqual(d._data_type, 'Experimental model')

    def test_comp_model_dataset(self):
        """Test ComparativeModelDataset"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d = IMP.pmi.metadata.ComparativeModelDataset(loc)
        self.assertEqual(d._data_type, 'Comparative model')

    def test_int_model_dataset(self):
        """Test IntegrativeModelDataset"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d = IMP.pmi.metadata.IntegrativeModelDataset(loc)
        self.assertEqual(d._data_type, 'Integrative model')

    def test_em_micrographs_dataset(self):
        """Test EMMicrographsDataset"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d = IMP.pmi.metadata.EMMicrographsDataset(loc, 400)
        self.assertEqual(d._data_type, 'EM raw micrographs')
        self.assertEqual(d.number, 400)
        d2 = IMP.pmi.metadata.EMMicrographsDataset(loc, 400)
        self.assertEqual(d, d2)
        # Not equal if number differs
        d3 = IMP.pmi.metadata.EMMicrographsDataset(loc, 600)
        self.assertNotEqual(d, d3)

    def test_em2d_class_dataset(self):
        """Test EM2DClassDataset"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d = IMP.pmi.metadata.EM2DClassDataset(loc)
        self.assertEqual(d._data_type, '2DEM class average')

    def test_sas_dataset(self):
        """Test SASDataset"""
        loc = ihm.location.InputFileLocation(repo='mydoi', path='a')
        d = IMP.pmi.metadata.SASDataset(loc)
        self.assertEqual(d._data_type, 'SAS data')

if __name__ == '__main__':
    IMP.test.main()
