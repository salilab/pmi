import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.io.input
import IMP.test
import RMF
import IMP.rmf
import os,sys


class AnalysisTest(IMP.test.TestCase):
    def setUp(self):
        IMP.test.TestCase.setUp(self)
        self.model = IMP.Model()

    def test_clustering(self):
        """Test clustering can calculate distance matrix, align, and cluster correctly"""
        pass

    def test_get_model_density(self):
        """Test GetModelDensity correctly creates and adds density maps"""
        custom_ranges={'med2':[(1,100,'med2')],
                       'med16':['med16']}
        mdens = IMP.pmi.analysis.GetModelDensity(custom_ranges)
        rmf_file=self.get_input_file_name('output/rmfs/2.rmf3')
        rh = RMF.open_rmf_file_read_only(rmf_file)
        prots = IMP.rmf.create_hierarchies(rh,self.model)
        IMP.rmf.load_frame(rh,0)
        mdens.add_subunits_density(prots[0])
        self.assertEqual(mdens.get_density_keys(),['med2','med16'])
        med2_coords=[]
        med16_coords=[]
        for i in range(4):
            IMP.rmf.load_frame(rh,i)
            s2=[child for child in prots[0].get_children()
                      if child.get_name() == 'med2']
            med2_coords+=[IMP.core.XYZ(p).get_coordinates() for p in
                         IMP.atom.Selection(s2,residue_indexes=range(1,100+1)).get_selected_particles()]
            s16=[child for child in prots[0].get_children()
                      if child.get_name() == 'med16']
            med16_coords+=[IMP.core.XYZ(p).get_coordinates() for p in
                           IMP.atom.Selection(s16).get_selected_particles()]
            mdens.add_subunits_density(prots[0])

        bbox2=IMP.algebra.BoundingBox3D(med2_coords)
        bbox16=IMP.algebra.BoundingBox3D(med16_coords)
        self.assertTrue(IMP.em.get_bounding_box(mdens.get_density('med2')).get_contains(bbox2))
        self.assertTrue(IMP.em.get_bounding_box(mdens.get_density('med16')).get_contains(bbox16))

    def test_analysis_macro(self):
        """Test the analysis macro does everything correctly"""
        pass

if __name__ == '__main__':
    IMP.test.main()
