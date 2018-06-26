from __future__ import print_function
import IMP
import os
import IMP.test
import IMP.core
import IMP.container
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.io
import IMP.pmi.io.crosslink
import IMP.pmi.representation
import IMP.pmi.restraints
import IMP.pmi.restraints.em
import IMP.pmi.restraints.crosslinking
import IMP.pmi.macros
import RMF
import IMP.rmf
from math import *

class Tests(IMP.test.TestCase):
    def test_particle_to_sample_filter(self):
        """Test ParticleToSampleFilter"""
        class MockRestraint(object):
            def __init__(self, sos):
                self.sos = sos
            def get_particles_to_sample(self):
                return self.sos
        r1 = MockRestraint({'Nuisances_Sigma': ('p0', 'p1'),
                            'Nuisances_Psi': ('p2', 'p3')})
        r2 = MockRestraint({'Nuisances_Sigma': ('p0', 'p4')})
        with IMP.allow_deprecated():
            p = IMP.pmi.tools.ParticleToSampleFilter([r1, r2])
        p.add_filter('Sigma')
        ps = p.get_particles_to_sample()
        self.assertEqual(list(ps.keys()), ['Nuisances_Sigma'])
        val = ps['Nuisances_Sigma']
        self.assertEqual(sorted(val), ['p0', 'p0', 'p1', 'p4'])

    def test_particle_to_sample_list(self):
        """Test ParticleToSampleList"""
        p = IMP.pmi.tools.ParticleToSampleList()
        self.assertEqual(p.label, 'None')
        self.assertRaises(TypeError, p.add_particle, 'P0', 'bad_type', 1, 'foo')

        p.add_particle('RB1', 'Rigid_Bodies', (1., 2.), 'myRB1')
        # Test bad rigid body transformation
        self.assertRaises(TypeError, p.add_particle,
                          'RB1', 'Rigid_Bodies', [1., 2.], 'myRB1')

        p.add_particle('S1', 'Surfaces', (1., 2., 3.), 'myS1')
        self.assertRaises(TypeError, p.add_particle,
                          'S1', 'Surfaces', [1., 2.], 'myS1')

        p.add_particle('F1', 'Floppy_Bodies', 1., 'myF1')
        self.assertRaises(TypeError, p.add_particle,
                          'F1', 'Floppy_Bodies', 'badtransform', 'myF1')

        self.assertEqual(p.get_particles_to_sample(),
                {'SurfacesParticleToSampleList_myS1_None':
                        (['S1'], (1.0, 2.0, 3.0)),
                 'Rigid_BodiesParticleToSampleList_myRB1_None':
                        (['RB1'], (1.0, 2.0)),
                 'Floppy_BodiesParticleToSampleList_myF1_None': (['F1'], 1.0)})

    def test_input_adaptor_non_pmi(self):
        mdl = IMP.Model()
        root=IMP.atom.Hierarchy(IMP.Particle(mdl))
        for i in range(10):
            p=IMP.Particle(mdl)
            h=IMP.atom.Hierarchy.setup_particle(p)
            IMP.atom.Mass.setup_particle(p,1.0)
            xyzr=IMP.core.XYZR.setup_particle(p)
            xyzr.set_coordinates((0,0,0))
            xyzr.set_radius(1.0)
            root.add_child(h)
        hs=IMP.pmi.tools.input_adaptor(root)
        self.assertEqual([IMP.atom.get_leaves(root)],hs)
        hs=IMP.pmi.tools.input_adaptor(root,pmi_resolution=1)
        self.assertEqual([IMP.atom.get_leaves(root)],hs)

    def test_Segments(self):
        s=IMP.pmi.tools.Segments(1)
        self.assertEqual(s.segs,[[1]])
        s=IMP.pmi.tools.Segments([1])
        self.assertEqual(s.segs,[[1]])
        s=IMP.pmi.tools.Segments([1,1])
        self.assertEqual(s.segs,[[1]])
        s=IMP.pmi.tools.Segments([1,2])
        self.assertEqual(s.segs,[[1,2]])
        s=IMP.pmi.tools.Segments([1,2,3])
        self.assertEqual(s.segs,[[1,2,3]])
        s=IMP.pmi.tools.Segments([1,2,3,5])
        self.assertEqual(s.segs,[[1,2,3],[5]])
        s.add(6)
        self.assertEqual(s.segs,[[1,2,3],[5,6]])
        s.add(0)
        self.assertEqual(s.segs,[[0,1,2,3],[5,6]])
        s.add(3)
        self.assertEqual(s.segs,[[0,1,2,3],[5,6]])
        s.add(4)
        self.assertEqual(s.segs,[[0,1,2,3,4,5,6]])
        s.add([-3,-4])
        self.assertEqual(s.segs,[[-4,-3],[0,1,2,3,4,5,6]])
        s.remove(2)
        self.assertEqual(s.segs,[[-4,-3],[0,1],[3,4,5,6]])
        s.remove(5)
        self.assertEqual(s.segs,[[-4,-3],[0,1],[3,4],[6]])
        s.remove(5)
        self.assertEqual(s.segs,[[-4,-3],[0,1],[3,4],[6]])
        s.add(-1)
        self.assertEqual(s.segs,[[-4,-3],[-1,0,1],[3,4],[6]])

    def assertEqualUnordered(self, a, b):
        """Compare two unordered lists; i.e. each list must have the
           same elements, but possibly in a different order"""
        self.assertEqual(len(a), len(b))
        for i in a + b:
            self.assertIn(i, a)
            self.assertIn(i, b)

    def test_threetoone(self):

        import string
        import random
        def id_generator(size=3, chars=string.ascii_uppercase + string.digits):
            return ''.join(random.choice(chars) for _ in range(size))

        threetoone = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 'UNK': 'X'}

        tto=IMP.pmi.tools.ThreeToOneConverter(is_nucleic=False)

        for key in threetoone:
            self.assertEqual(threetoone[key],tto[key])

        for s in range(10):
            id=id_generator()
            if id in threetoone:
                self.assertEqual(threetoone[id], tto[id])
            else:
                self.assertEqual("X",tto[id])

        threetoone = {'ADE': 'A', 'URA': 'U', 'CYT': 'C', 'GUA': 'G',
                      'THY': 'T', 'UNK': 'X'}

        tto = IMP.pmi.tools.ThreeToOneConverter(is_nucleic=True)

        for key in threetoone:
            self.assertEqual(threetoone[key], tto[key])


        for s in range(10):
            id = id_generator()
            if id in threetoone:
                self.assertEqual(threetoone[id], tto[id])
            else:
                self.assertEqual("X", tto[id])

    def test_get_restraint_set(self):
        """Test get_restraint_set()"""
        m = IMP.Model()
        # Should make an empty set
        rs = IMP.pmi.tools.get_restraint_set(m)
        self.assertEqual(rs.get_number_of_restraints(), 0)

        for rmf in (True, False):
            rs = IMP.pmi.tools.get_restraint_set(m, rmf)
            self.assertEqual(rs.get_number_of_restraints(), 0)

    def test_add_restraint(self):
        """Test add_restraint_to_model()"""
        m = IMP.Model()

        r1 = IMP._ConstRestraint(m, [], 1)
        IMP.pmi.tools.add_restraint_to_model(m, r1, add_to_rmf=False)
        r2 = IMP._ConstRestraint(m, [], 1)
        IMP.pmi.tools.add_restraint_to_model(m, r2, add_to_rmf=True)

        rs = IMP.pmi.tools.get_restraint_set(m, rmf=False)
        self.assertEqual(rs.get_number_of_restraints(), 2)

        rs = IMP.pmi.tools.get_restraint_set(m, rmf=True)
        self.assertEqual(rs.get_number_of_restraints(), 1)


if __name__ == '__main__':
    IMP.test.main()
