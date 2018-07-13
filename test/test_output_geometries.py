import os

import IMP
import IMP.test
import IMP.display
import IMP.algebra
import IMP.atom
import IMP.core
import IMP.pmi1.macros
import IMP.pmi1.output
import IMP.rmf
import RMF


class Tests(IMP.test.TestCase):

    def test_round_trip(self):
        mdl = IMP.Model()
        p = IMP.Particle(mdl)
        IMP.core.XYZR.setup_particle(p)
        IMP.atom.Mass.setup_particle(p, 1)
        h = IMP.atom.Hierarchy.setup_particle(p)
        geo = IMP.display.SphereGeometry(
            IMP.algebra.Sphere3D(IMP.algebra.Vector3D(1, 2, 3), 4.))
        rmf_fn = self.get_tmp_file_name("geometrytest.rmf")
        o = IMP.pmi1.output.Output()
        o.init_rmf(rmf_fn, [h], geometries=[geo])
        o.write_rmf(rmf_fn)
        o.close_rmf(rmf_fn)
        f = RMF.open_rmf_file_read_only(rmf_fn)
        rmf_geos = IMP.rmf.create_geometries(f)
        IMP.rmf.load_frame(f, RMF.FrameID(0))
        self.assertEqual(len(rmf_geos), 1)
        self.assertEqual(len(rmf_geos[0].get_components()), 1)


if __name__ == "__main__":
    IMP.test.main()
