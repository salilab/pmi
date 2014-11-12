import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.macros
import os
import IMP.test
import tempfile

class TopologyReaderTests(IMP.test.TestCase):
    def test_reading(self):
        '''Test basic reading'''
        topology_file=self.get_input_file_name("topology.txt")
        t=IMP.pmi.topology.TopologyReader(topology_file)
        self.assertEqual(len(t.component_list),15)
        self.assertEqual(t.component_list[4].domain_name,"Rpb2_2")
        self.assertEqual(t.component_list[5].name,"Rpb3")
    def test_round_trip(self):
        '''Test reading and writing'''
        topology_file=self.get_input_file_name("topology.txt")
        outfile = tempfile.NamedTemporaryFile('w')
        t=IMP.pmi.topology.TopologyReader(topology_file)
        t.write_topology_file(outfile.name)

        t=IMP.pmi.topology.TopologyReader(outfile.name)
        self.assertEqual(len(t.component_list),15)
        self.assertEqual(t.component_list[4].domain_name,"Rpb2_2")
        self.assertEqual(t.component_list[5].name,"Rpb3")
    def test_build(self):
        '''Test building with macro BuildModel1 using a topology file'''
        pass
        #m = IMP.Model()
        #simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)
        #bm=IMP.pmi.macros.BuildModel1(simo)
        #bm.build_model(component_topologies=t.component_list[6:8])

if __name__=="__main__":
    IMP.test.main()
