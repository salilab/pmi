import IMP.pmi
import IMP.pmi.io.input
import IMP.test
import os,sys


class InputTest(IMP.test.TestCase):
    def test_get_best_models(self):
        input_dir=os.path.dirname(self.get_input_file_name('chainA.pdb'))
        stat_files = [os.path.join(input_dir,'output1','stat.0.out'),
                      os.path.join(input_dir,'output1','stat.1.out'),
                      os.path.join(input_dir,'output2','stat.0.out'),
                      os.path.join(input_dir,'output2','stat.1.out')]
        feature_keys = ['CHARMM']
        results = IMP.pmi.io.input.get_best_models(stat_files,
                                                   "SimplifiedModel_Total_Score",
                                                   feature_keys,
                                                   prefiltervalue=240.0)

        rmf_file_list,rmf_file_frame_list,score_list,feature_keyword_list_dict=results
        self.assertEqual(len(rmf_file_list),9)
        self.assertEqual(len(rmf_file_frame_list),9)
        self.assertEqual(len(score_list),9)
        self.assertEqual('CHARMM_BONDS' in feature_keyword_list_dict,True)



        #self.assertEqual(rmf_file_,set(((os.path.join(dir1,'rmfs/1.rmf3'),8),
        #                               (os.path.join(dir1,'rmfs/1.rmf3'),15),
        #                               (os.path.join(dir1,'rmfs/1.rmf3'),14))))


if __name__ == '__main__':
    IMP.test.main()
