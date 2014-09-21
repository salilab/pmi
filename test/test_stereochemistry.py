import IMP
import IMP.atom
import IMP.pmi
import IMP.pmi.representation_new as r
import IMP.pmi.restraints_new.stereochemistry
import IMP.pmi.data_tools
import IMP.test
import IMP.pmi.sequence_tools


class StereochemistryTests(IMP.test.TestCase):
    def test_elastic_network(self):
        """ test PMI setup of elastic nets """

        # create two states, each with two copies of the protein
        s=r.System()
        seqs=r.Sequences(self.get_input_file_name('chainA.fasta'),
                         name_map={'GCP2_YEAST':'Prot1'})
        # build state
        st1=s.create_state()
        m1=st1.create_molecule("Prot1",sequence=seqs["Prot1"],chain_id='A')
        atomic_res=m1.add_structure(self.get_input_file_name('chainA.pdb'),
                                    chain_id='A',
                                    model_num=0)
        m1.add_representation(atomic_res,resolutions=[0])

        hier = s.build(merge_type="backbone")

        # get secondary structure
        sse_selections=IMP.pmi.data_tools.parse_dssp(
            self.get_input_file_name('chainA.dssp'),'A')

        # create elastic network from some SSEs
        data=sse_selections['helix'][0][0]
        er = IMP.pmi.restraints_new.stereochemistry.ElasticNetworkRestraint(hier,
                                                selection_dict=data,
                                                extra_sel={'atom_type':IMP.atom.AtomType("CA")},
                                                strength=10.0,
                                                dist_cutoff=5.0)
        print er.get_restraint()
        self.assertEqual(er.get_restraint().get_number_of_restraints(),12)

if __name__ == '__main__':
    IMP.test.main()
