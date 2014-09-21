import IMP.pmi
from IMP.pmi.io.data_storage import SelectionDict,SubsequenceData,CrossLinkData
import IMP.pmi.io.data_parsers as data_parsers
import IMP.test
import os,sys

class DataStorage(IMP.test.TestCase):
    def setUp(self):
        IMP.test.TestCase.setUp(self)
        self.mdl = IMP.Model()

    def test_selection_dictionary(self):
        sd = SelectionDict(self.mdl,{'chain':'A'},residue_index=5)
        sd.add_fields(molecule='himom')
        self.assertEqual(cmp(sd.get_data(),{
            'residue_index': 5, 'molecule': 'himom', 'chain': 'A'}),0)

    def test_subsequence_data(self):
        global_sdict = SelectionDict(self.mdl,
                                     atom_type=IMP.atom.AtomType("CA"))
        ssd = SubsequenceData(self.mdl,global_sdict)
        ssd.add_subsequence('helix',[SelectionDict(self.mdl,
                                                   molecule='himom',
                                                   residue_indexes=range(55,58))])
        ssd.add_subsequence('beta',[SelectionDict(self.mdl,
                                                  molecule='himom',
                                                  residue_indexes=range(5,10)),
                                    SelectionDict(self.mdl,
                                                  molecule='hidad',
                                                  residue_indexes=range(95,100))])

        self.assertEqual(cmp(ssd['helix'][0],
                             [SelectionDict(self.mdl,
                                            residue_indexes=range(55,58),
                                            molecule='himom',
                                            atom_type=IMP.atom.AtomType("CA"))]),0)
        self.assertEqual(cmp(ssd['beta'][0],
                             [SelectionDict(self.mdl,
                                            residue_indexes=range(5,10),
                                            molecule='himom',
                                            atom_type=IMP.atom.AtomType("CA")),
                              SelectionDict(self.mdl,
                                            residue_indexes=range(95,100),
                                            molecule='hidad',
                                            atom_type=IMP.atom.AtomType("CA"))]),0)
    def test_dssp_parsing(self):
        """Test reading DSSP files"""
        sses = data_parsers.parse_dssp(self.mdl,self.get_input_file_name('chainA.dssp'),'A')
        self.assertEqual(sorted(sses.keys()),sorted(['helix','beta','loop']))
        self.assertEqual(len(sses['helix']),20)
        self.assertEqual(len(sses['beta']),3)
        self.assertEqual(len(sses['loop']),32)

    def test_crosslink_data(self):
        """Test the CrossLinkData storage class"""
        global_sdict = SelectionDict(self.mdl,
                                     atom_type=IMP.atom.AtomType("CA"))
        xld = CrossLinkData(self.mdl,global_sdict)
        xld.add_cross_link(0,
                           SelectionDict(self.mdl,
                                         residue_index=5,
                                         molecule='himom'),
                           SelectionDict(self.mdl,
                                         residue_index=10,
                                         molecule='himom'),
                           score=5.0)
        xld.add_cross_link(0,
                           SelectionDict(self.mdl,
                                         residue_index=25,
                                         molecule='himom'),
                           SelectionDict(self.mdl,
                                         residue_index=35,
                                         molecule='himom'),
                           score=5.0)
        xld.add_cross_link(1,
                           SelectionDict(self.mdl,
                                         residue_index=225,
                                         molecule='hidad'),
                           SelectionDict(self.mdl,
                                         residue_index=335,
                                         molecule='hidad'),
                           score=5.0)

    def test_xl_parsing(self):
        """Test the XL davis parser"""
        data = data_parsers.parse_xlinks_davis(self.mdl,
                                               self.get_input_file_name('xls_davis.txt'),
                                               name_map={'His-TEV-Tub4':'ytub'},
                                               named_offsets={'ytub':-33})
        self.assertEqual(len(data),41)
        self.assertEqual(cmp(data[0],
                             [{'r1':SelectionDict(self.mdl,
                                                 residue_index=337,
                                                 molecule='ytub'),
                              'r2':SelectionDict(self.mdl,
                                                 residue_index=831,
                                                 molecule='Spc98'),
                              'score':10.88819}]),0)

"""
class InputTest(IMP.test.TestCase):
    def test_get_best_models(self):
        rank = 0
        number_of_processes = 1
        comm = None

        root_dir=os.path.dirname(os.path.abspath(sys.argv[0]))
        dir1 = os.path.join(root_dir,"input/output1")
        dir2 = os.path.join(root_dir,"input/output2")
        gbm = IMP.pmi.io.input.GetBestModels(data_dirs=[dir1,dir2],
                                             score_key='AtomicXLRestraint',
                                             number_of_processes=number_of_processes,
                                             number_of_best_scoring_models=3)
        result = gbm.get_best_models(rank)
        test_set=set()
        for r in result:
            test_set.add(tuple(r[1:3]))
        self.assertEqual(test_set,set(((os.path.join(dir1,'rmfs/1.rmf3'),8),
                                       (os.path.join(dir1,'rmfs/1.rmf3'),15),
                                       (os.path.join(dir1,'rmfs/1.rmf3'),14))))
"""

if __name__ == '__main__':
    IMP.test.main()
