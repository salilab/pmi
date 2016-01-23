## \example pmi/multiscale.py

"""This script shows how to represent a system
At multiple scales, including electron densities.
"""

import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.test
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry

IMP.set_log_level(IMP.SILENT)

class MultiscaleTopologyTest(IMP.test.TestCase):

      def initialize_system(self, mdl):

            s = IMP.pmi.topology.System(mdl)
            st1 = s.create_state()

            # Read sequences and create Molecules
            seqs = IMP.pmi.topology.Sequences('../examples/data/gcp2.fasta')
            mol = st1.create_molecule("GCP2",sequence=seqs["GCP2_YEAST"],chain_id='A')

            # Add structure. This function returns a list of the residues that now have structure
            a1 = mol.add_structure('../examples/data/gcp2.pdb',
                                    chain_id='A')

            # Add representations. For structured regions, created a few beads as well as densities
            #  For unstructured regions, create a single bead level and set those up as densities
            mol.add_representation(a1,
                                    resolutions=[10,100],
                                    density_prefix='../examples/data/gcp2_gmm',
                                    density_residues_per_component=20,
                                    density_voxel_size=3.0)
            mol.add_representation(mol.get_non_atomic_residues(),
                                    resolutions=[10],
                                    setup_particles_as_densities=True)

            # When you call build, this actually makes the beads and fits the GMMs
            #  This returns a canonical IMP hierarchy
            hier = s.build()
            #IMP.atom.show_molecular_hierarchy(hier)

            return a1, hier, mol


      def test_num_residues(self):
            mdl = IMP.Model()
            (a1, hier, mol)=self.initialize_system(mdl)

            num_res_get_res = len(mol.get_residues())
            num_res = len(mol.residues)
            len_struct_list = len(a1)
            len_struct_list_mol = len(mol.get_atomic_residues())

            self.assertEqual(len(a1), len(mol.get_atomic_residues()))
            self.assertEqual(len(mol.get_residues()), len(mol.residues))
            self.assertEqual(823, len(mol.residues))
            self.assertEqual(581, len(a1))
            self.assertEqual(581, len(mol.get_atomic_residues()))

      def test_num_unstruct_res(self):
            mdl = IMP.Model()
            (a1, hier, mol)=self.initialize_system(mdl)

            struct_res = 0
            unstruct_res = 0
            for res in mol.residues:
                  if res.get_has_structure():
                        struct_res+=1
                  else:
                        unstruct_res+=1

            self.assertEqual(242, unstruct_res) # Should be 843 - 581 residues
            self.assertEqual(581, struct_res)
            self.assertEqual(len(mol.residues), unstruct_res + struct_res)

      def test_num_struct_res(self):
            mdl = IMP.Model()
            (a1, hier, mol)=self.initialize_system(mdl)

            struct_res = 0
            unstruct_res = 0
            for res in mol.residues:
                  if res.get_has_structure():
                        struct_res+=1
                  else:
                        unstruct_res+=1

            self.assertEqual(581, struct_res)
            self.assertEqual(len(mol.residues), unstruct_res + struct_res)

      def test_residue_type(self):
            # The self.hier object in TempResidue is an IMP.atom.Hierarchy
            mdl = IMP.Model()
            (a1, hier, mol)=self.initialize_system(mdl)

            res = mol.residues[0]

            self.assertEqual(IMP.atom.Residue, type(res.hier))


  
      def test_residue_print(self):
            # PMI Residues cannot print their name
            # The self.hier object in TempResidue is an IMP.atom.Hierarchy
            # IMP.atom.Hierarchy does not have a function get_residue_type()

            mdl = IMP.Model()
            (a1, hier, mol)=self.initialize_system(mdl)

            res = mol.residues[0]

            try:
                  print(res)
            except:
                  self.fail("Cannot print Residue")

            try:
                  print(res.get_code())
            except:
                  self.fail("Cannot print residue code")


      def test_molecule_rigid_members(self):
            # None of the leaves of the molecule are RigidMembers

            mdl = IMP.Model()
            (a1, hier, mol)=self.initialize_system(mdl)

            dof = IMP.pmi.dof.DegreesOfFreedom(mdl)
            dof.create_rigid_body(mol,
                      nonrigid_parts=mol.get_non_atomic_residues())

            rb = IMP.core.RigidBodyMember(IMP.atom.get_leaves(mol.get_hierarchy())[0]).get_rigid_body()
            nrms = 0
            rms = 0
            all_things = IMP.atom.get_leaves(mol.get_hierarchy())
            
            for part in all_things:
                  if IMP.core.NonRigidMember.get_is_setup(part):
                        nrms += 1
                  elif IMP.core.RigidMember.get_is_setup(part):
                        rms += 1
                  else:
                        self.fail("Particle not a RigidMember or a NonRigidMember")

            self.assertNotEqual(0,rms)

      def test_molecule_rigid_members(self):
            # None of the leaves of the selection (Resolution=10) are RigidMembers

            mdl = IMP.Model()
            (a1, hier, mol)=self.initialize_system(mdl)

            dof = IMP.pmi.dof.DegreesOfFreedom(mdl)

            # Create rigid body
            dof.create_rigid_body(mol,
                      nonrigid_parts=mol.get_non_atomic_residues())

            rb = IMP.core.RigidBodyMember(IMP.atom.get_leaves(mol.get_hierarchy())[0]).get_rigid_body()
            nrms = 0
            rms = 0

            selection = IMP.atom.Selection(hierarchy=mol.get_hierarchy(),resolution=10).get_hierarchies()

            all_things = IMP.atom.get_leaves(selection[0])
            
            for part in all_things:
                  if IMP.core.NonRigidMember.get_is_setup(part):
                        nrms += 1
                  elif IMP.core.RigidMember.get_is_setup(part):
                        rms += 1
                  else:
                        self.fail("Particle not a RigidMember or a NonRigidMember")

            self.assertNotEqual(0,rms)


      def test_molecule_rigid_members2(self):
            # When the rigid body is created with no assigned nonrigid members
            # all leaves of the molecule are of type RigidMember

            mdl = IMP.Model()
            (a1, hier, mol)=self.initialize_system(mdl)

            dof = IMP.pmi.dof.DegreesOfFreedom(mdl)

            # Create rigid body
            dof.create_rigid_body(mol)

            rb = IMP.core.RigidBodyMember(IMP.atom.get_leaves(mol.get_hierarchy())[0]).get_rigid_body()
            nrms = 0
            rms = 0

            selection = IMP.atom.Selection(hierarchy=mol.get_hierarchy(),resolution=10).get_hierarchies()

            all_things = IMP.atom.get_leaves(selection[0])
            
            for part in all_things:
                  if IMP.core.NonRigidMember.get_is_setup(part):
                        nrms += 1
                  elif IMP.core.RigidMember.get_is_setup(part):
                        rms += 1
                  else:
                        self.fail("Particle not a RigidMember or a NonRigidMember")

            self.assertEqual(0,nrms)


if __name__ == '__main__':
    IMP.test.main()   



