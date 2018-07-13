from __future__ import print_function
import IMP.test
import IMP.pmi1.output
import IMP.pmi1.representation

class TestPO(IMP.pmi1.output.ProtocolOutput):
    _file_datasets = []
    _each_metadata = []
    def _add_state(self, obj):
        return self

class Tests(IMP.test.TestCase):

    def test_prot_add(self):
        """Test Representation.add_protocol_output()"""
        m = IMP.Model()
        r = IMP.pmi1.representation.Representation(m)
        po = TestPO()
        r.add_protocol_output(po)

if __name__ == '__main__':
    IMP.test.main()
