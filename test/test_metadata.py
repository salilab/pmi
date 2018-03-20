from __future__ import print_function
import IMP.test
import IMP.pmi.metadata
import IMP.pmi.representation
import ihm.location
import os

class Tests(IMP.test.TestCase):

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

if __name__ == '__main__':
    IMP.test.main()
