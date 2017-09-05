import IMP
import IMP.pmi
import IMP.pmi.macros
import IMP.test
import glob

class Tests(IMP.test.TestCase):

    def test_analysis_replica_exchange(self):
        model=IMP.Model()
        sts=sorted(glob.glob(self.get_input_file_name("output_test/stat.0.out").replace(".0.",".*.")))
        are=IMP.pmi.macros.AnalysisReplicaExchange(model,sts,10)

        are.set_alignment_selection(molecule="Rpb4")
        are.cluster(20)

        print(are)

        dcr={"Rpb4":["Rpb4"],"Rpb7":["Rpb7"],"All":["Rpb4","Rpb7"]}

        self.assertEqual(len(are),4)

        for cluster in are:
            print(cluster)
            are.save_coordinates(cluster)
            are.save_densities(cluster,dcr)
            for member in cluster:
                print(member)



if __name__ == '__main__':
    IMP.test.main()
