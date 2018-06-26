import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.macros
import os
import IMP.test
import IMP.rmf
import IMP.pmi.plotting
import IMP.pmi.plotting.topology

def children_as_dict(h):
    cdict={}
    for c in h.get_children():
        cdict[c.get_name()]=c
    return cdict

class Tests(IMP.test.TestCase):
    def test_old(self):
        """Test reading of old-style topology file"""
        topology_file = self.get_input_file_name("topology.txt")
        with IMP.allow_deprecated():
            t = IMP.pmi.topology.TopologyReader(topology_file)
        c = t.get_components()
        self.assertEqual(len(c), 3)
        self.assertEqual(c[0].molname, "Prot1")
        self.assertEqual(os.path.abspath(c[0].fasta_file),
                         self.get_input_file_name("seqs.fasta"))
        self.assertEqual(c[1].molname, "Prot2")
        with IMP.allow_deprecated():
            # Test deprecated interface
            self.assertEqual(c[1].name, "Prot2")
            self.assertEqual(c[1].domain_name, "Prot2A")
        self.assertEqual(c[1].get_unique_name(), "Prot2..0")
        self.assertEqual(c[2].get_unique_name(), "Prot2..1")

    def test_reading(self):
        """Test basic reading"""
        topology_file = self.get_input_file_name("topology_new.txt")
        t = IMP.pmi.topology.TopologyReader(topology_file)
        self.assertEqual(list(t.molecules.keys()),
                         ['Prot1', 'Prot2', 'Prot3', 'Prot4', 'Prot5'])
        c1 = t.get_components()
        with IMP.allow_deprecated():
            # Test deprecated interface
            c2 = t.component_list
        for c in (c1, c2):
            self.assertEqual(len(c),9)
            self.assertEqual(c[0].molname,"Prot1")
            self.assertEqual(c[1].molname,"Prot1")
            self.assertEqual(c[1].copyname,"1")
            self.assertEqual(c[2].get_unique_name(),"Prot2..0")
            self.assertEqual(c[3].get_unique_name(),"Prot2..1")
            self.assertEqual(c[5].get_unique_name(),"Prot2.1.1")

    def test_round_trip(self):
        """Test reading and writing"""
        topology_file=self.get_input_file_name("topology_new.txt")
        outfile = self.get_tmp_file_name("ttest.txt")
        t=IMP.pmi.topology.TopologyReader(topology_file)
        t.write_topology_file(outfile)

        tnew = IMP.pmi.topology.TopologyReader(outfile)
        c = tnew.get_components()
        self.assertEqual(len(c),9)
        self.assertEqual(c[0].molname,"Prot1")
        self.assertEqual(c[1].molname,"Prot1")
        self.assertEqual(c[1].copyname,"1")
        self.assertEqual(c[5].get_unique_name(),"Prot2.1.1")

    def test_set_movers(self):
        """Check if rigid bodies etc are set up as requested"""
        try:
            import sklearn
        except ImportError:
            self.skipTest("no sklearn package")
        mdl = IMP.Model()
        tfile = self.get_input_file_name('topology_new.txt')
        input_dir = os.path.dirname(tfile)
        t = IMP.pmi.topology.TopologyReader(tfile,
                                            pdb_dir=input_dir,
                                            fasta_dir=input_dir,
                                            gmm_dir=input_dir)
        comps = t.get_components()
        self.assertEqual(comps[0].pdb_file, os.path.join(input_dir,'prot.pdb'))
        rbs = t.get_rigid_bodies()
        srbs = t.get_super_rigid_bodies()
        csrbs = t.get_chains_of_super_rigid_bodies()

        expected_rbs = [['Prot1.1.0','Prot1..0'],
                        ['Prot2..0','Prot2..1','Prot2.1.0','Prot2.1.1'],
                        ['Prot4..0', 'Prot5..0']]
        expected_srbs = [['Prot1.1.0','Prot1..0','Prot2..0','Prot2..1',
                          'Prot2.1.0','Prot2.1.1','Prot4..0','Prot3..0', 'Prot5..0'],
                         ['Prot1.1.0','Prot1..0','Prot3..0']]

        found1 = set(tuple(sorted(i)) for i in rbs)
        found2 = set(tuple(sorted(i)) for i in expected_rbs)
        self.assertEqual(found1,found2)

        found1 = set(tuple(sorted(i)) for i in srbs)
        found2 = set(tuple(sorted(i)) for i in expected_srbs)
        self.assertEqual(found1,found2)


if __name__=="__main__":
    IMP.test.main()
