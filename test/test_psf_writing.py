import os
import IMP
import IMP.test

import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.em
import IMP.pmi.representation
import IMP.pmi.tools
import IMP.pmi.output

class Tests(IMP.test.TestCase):
    def test_psf_writing(self):
        """Test writing of PSF files"""

        # input parameter
        pdbfile = self.get_input_file_name("mini.pdb")
        fastafile = self.get_input_file_name("mini.fasta")

        components = ["Rpb1", "Rpb2" ]

        chains = "AB"

        colors = [0., 1.0]

        beadsize = 1

        fastids = IMP.pmi.tools.get_ids_from_fasta_file(fastafile)


        m = IMP.Model()
        simo = IMP.pmi.representation.Representation(m)


        simo.create_component("Rpb1", color=colors[0])
        simo.add_component_sequence("Rpb1", fastafile, id=fastids[0])
        simo.autobuild_model("Rpb1", pdbfile, "A",
                             resolutions=[1], missingbeadsize=beadsize)
        simo.setup_component_sequence_connectivity("Rpb1", 1)

        simo.create_component("Rpb2", color=colors[1])
        simo.add_component_sequence("Rpb2", fastafile, id=fastids[1])
        simo.autobuild_model("Rpb2", pdbfile, "B",
                             resolutions=[1], missingbeadsize=beadsize)
        simo.setup_component_sequence_connectivity("Rpb2", 1)

        output = IMP.pmi.output.Output()
        output.init_pdb("test_psf_writing.pdb", simo.prot)
        output.write_pdb("test_psf_writing.pdb")
        output.write_psf("test_psf_writing.psf","test_psf_writing.pdb")
        psf_content='''PSF CMAP CHEQ
17 !NATOM
       1 A    1    "MET" C    C         1.000000      0.000000       0      0.000000      0.000000
       2 A    2    "VAL" C    C         1.000000      0.000000       0      0.000000      0.000000
       3 A    3    "GLY" C    C         1.000000      0.000000       0      0.000000      0.000000
       4 A    4    "GLN" C    C         1.000000      0.000000       0      0.000000      0.000000
       5 A    5    "GLN" C    C         1.000000      0.000000       0      0.000000      0.000000
       6 A    6    "TYR" C    C         1.000000      0.000000       0      0.000000      0.000000
       7 A    7    "SER" C    C         1.000000      0.000000       0      0.000000      0.000000
       8 A    8    "SER" C    C         1.000000      0.000000       0      0.000000      0.000000
       9 B    1    "ALA" C    C         1.000000      0.000000       0      0.000000      0.000000
      10 B    2    "ALA" C    C         1.000000      0.000000       0      0.000000      0.000000
      11 B    3    "ASP" C    C         1.000000      0.000000       0      0.000000      0.000000
      12 B    4    "GLU" C    C         1.000000      0.000000       0      0.000000      0.000000
      13 B    5    "SER" C    C         1.000000      0.000000       0      0.000000      0.000000
      14 B    6    "ALA" C    C         1.000000      0.000000       0      0.000000      0.000000
      15 B    7    "PRO" C    C         1.000000      0.000000       0      0.000000      0.000000
      16 B    8    "ILE" C    C         1.000000      0.000000       0      0.000000      0.000000
      17 B    9    "THR" C    C         1.000000      0.000000       0      0.000000      0.000000
15 !NBOND: bonds
       1       2       2       3       3       4       4       5
       5       6       6       7       7       8       9      10
      10      11      11      12      12      13      13      14
      14      15      15      16      16      17'''.split("\n")

        with open("test_psf_writing.psf") as f:
            for nl, l in enumerate(f):
                self.assertEqual(psf_content[nl],l.replace('\n',''))
        os.unlink('test_psf_writing.psf')
        os.unlink('test_psf_writing.pdb')

    def test_pmi2_pdb_writing_with_gaussians(self):
        """Test writing of PDB/PSF files"""

        # input parameter
        beadsize = 1
        mdl = IMP.Model()
        s = IMP.pmi.topology.System(mdl)
        st1 = s.create_state()
        seqs = IMP.pmi.topology.Sequences(self.get_input_file_name('seqs.fasta'))

        m1 = st1.create_molecule("Prot1",sequence=seqs["Protein_1"])
        atomic_res = m1.add_structure(self.get_input_file_name('prot.pdb'),
                                      chain_id='A',res_range=(55,63),offset=-54)

        fname = self.get_tmp_file_name('test_gmm')
        m1.add_representation(atomic_res,resolutions=[1,10],
                              density_residues_per_component=2,
                              density_voxel_size=3.0,
                              density_prefix=fname)

        m1.add_representation(m1.get_non_atomic_residues(),resolutions=[5])
        hier = m1.build()
        densities = IMP.atom.Selection(hier,representation_type=IMP.atom.DENSITIES).get_selected_particles() +\
                    [r.get_hierarchy() for r in m1.get_non_atomic_residues()]
        gem = IMP.pmi.restraints.em.GaussianEMRestraint(densities,
                                                        target_fn=self.get_input_file_name('prot_gmm.txt'),
                                                        target_is_rigid_body=True)

        gem.add_to_model()
        gem.add_target_density_to_hierarchy(st1.get_hierarchy())

        output = IMP.pmi.output.Output()
        output.init_pdb("test_psf_writing.pdb", s.get_hierarchy())
        output.write_pdb("test_psf_writing.pdb")
        output.write_psf("test_psf_writing.psf","test_psf_writing.pdb")

        os.unlink('test_psf_writing.psf')
        os.unlink('test_psf_writing.pdb')

if __name__ == '__main__':
    IMP.test.main()
