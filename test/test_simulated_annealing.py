import IMP
import IMP.test
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.samplers
import IMP.pmi.tools
import IMP.pmi.restraints.basic
import ast
import os
import shutil


class TwoReplicaExchange(IMP.pmi.samplers._SerialReplicaExchange):
    """Act as the first of two replicas (temperatures tmin and tmax), and
       record the arguments of every attempted exchange"""
    def __init__(self, accept):
        super().__init__()
        self.accept = accept
        self.exchanges = []

    def get_number_of_replicas(self):
        return 2

    def create_temperatures(self, tmin, tmax, nrep):
        self.temperatures = [tmin, tmax]
        return self.temperatures

    def get_friend_index(self, step):
        return 1

    def get_friend_parameter(self, key, findex):
        return [self.temperatures[findex]]

    def do_exchange(self, myscore, fscore, findex):
        self.exchanges.append((myscore, fscore))
        return self.accept


def make_particle(m):
    p = IMP.Particle(m)
    d = IMP.core.XYZ.setup_particle(p, IMP.algebra.Vector3D(0, 0, 0))
    d.set_coordinates_are_optimized(True)
    IMP.atom.Mass.setup_particle(p, 1.0)
    r = IMP.core.SingletonRestraint(
        m, IMP.core.DistanceToSingletonScore(
            IMP.core.Harmonic(0., 1.), IMP.algebra.Vector3D(0, 0, 0)), p)
    IMP.pmi.tools.add_restraint_to_model(m, r)
    return p


def make_mc(m):
    p = make_particle(m)
    return IMP.pmi.samplers.MonteCarlo(m, [IMP.core.BallMover(m, p, 1.0)],
                                       2.0)


class Tests(IMP.test.TestCase):

    def test_mc_temperatures(self):
        """Test that simulated annealing scales the MC temperature"""
        m = IMP.Model()
        mc = make_mc(m)
        mc.set_simulated_annealing(1.0, 2.5, 2, 1)
        kts = []
        for i in range(6):
            mc.optimize(1)
            kts.append(mc.get_mc().get_kt())
            # Stat file should report the temperature just used
            self.assertAlmostEqual(
                float(mc.get_output()["MonteCarlo_Temperature"]), kts[-1],
                delta=1e-6)
        for kt, expected in zip(kts, [2.0, 2.0, 5.0, 2.0, 2.0, 5.0]):
            self.assertAlmostEqual(kt, expected, delta=1e-6)
        # The factor applies to any new temperature (e.g. from a replica
        # exchange swap)
        mc.set_kt(3.0)
        self.assertAlmostEqual(mc.get_mc().get_kt(), 7.5, delta=1e-6)
        mc.optimize(1)
        self.assertAlmostEqual(mc.get_mc().get_kt(), 3.0, delta=1e-6)
        self.assertRaises(ValueError, mc.set_simulated_annealing,
                          0.0, 2.5, 2, 1)

    def test_md_temperatures(self):
        """Test that simulated annealing scales the MD temperature"""
        m = IMP.Model()
        p = make_particle(m)
        md = IMP.pmi.samplers.MolecularDynamics(m, [p], 2.0)
        md.set_simulated_annealing(1.0, 2.5, 2, 1)
        kts = []
        for i in range(3):
            md.optimize(1)
            kts.append(md.ltstate.get_temperature() * 0.0019872041)
        for kt, expected in zip(kts, [2.0, 2.0, 5.0]):
            self.assertAlmostEqual(kt, expected, delta=1e-6)

    def test_swap(self):
        """Test replica exchange swaps with simulated annealing"""
        for accept in (False, True):
            m = IMP.Model()
            mc = make_mc(m)
            mc.set_simulated_annealing(1.0, 2.5, 1, 1)
            rem = TwoReplicaExchange(accept)
            rex = IMP.pmi.samplers.ReplicaExchange(
                m, 1.0, 2.0, [mc], replica_exchange_object=rem)
            self.assertAlmostEqual(mc.get_mc().get_kt(), 1.0, delta=1e-6)
            # Frame 0 (factor 1) and frame 1 (factor 2.5)
            for frame, factor in ((0, 1.0), (1, 2.5)):
                mc.optimize(1)
                self.assertAlmostEqual(rex.get_temperature_factor(), factor,
                                       delta=1e-6)
                rex.swap_temp(frame, score=10.0)
                # Exchange is attempted with the temperatures actually used
                myscore, fscore = rem.exchanges[-1]
                self.assertAlmostEqual(myscore, 10.0 / (1.0 * factor),
                                       delta=1e-6)
                self.assertAlmostEqual(fscore, 10.0 / (2.0 * factor),
                                       delta=1e-6)
                # An accepted exchange keeps the annealing factor
                expected_kt = 2.0 * factor if accept else 1.0 * factor
                self.assertAlmostEqual(mc.get_mc().get_kt(), expected_kt,
                                       delta=1e-6)

    def test_macro_output(self):
        """Test that frames at high annealing temperature are not saved"""
        m = IMP.Model()
        s = IMP.pmi.topology.System(m)
        st1 = s.create_state()
        nup84 = st1.create_molecule("Nup84", "MELS", "X")
        nup84.add_representation(resolutions=[1])
        hier = s.build()
        dof = IMP.pmi.dof.DegreesOfFreedom(m)
        dof.create_flexible_beads(nup84, max_trans=1.0)
        dr = IMP.pmi.restraints.basic.DistanceRestraint(
            root_hier=hier, tuple_selection1=(1, 1, "Nup84"),
            tuple_selection2=(4, 4, "Nup84"), distancemin=5, distancemax=5)
        dr.add_to_model()
        outdir = 'test_sa_output/'
        rex = IMP.pmi.macros.ReplicaExchange(
            m, root_hier=hier, monte_carlo_sample_objects=dof.get_movers(),
            output_objects=[dr], monte_carlo_steps=5, number_of_frames=6,
            number_of_best_scoring_models=0,
            simulated_annealing=True,
            simulated_annealing_minimum_temperature=1.0,
            simulated_annealing_maximum_temperature=2.5,
            simulated_annealing_minimum_temperature_nframes=2,
            simulated_annealing_maximum_temperature_nframes=1,
            global_output_directory=outdir)
        rex.execute_macro()
        with open(os.path.join(outdir, 'stat.0.out')) as fh:
            header = ast.literal_eval(fh.readline())
            key = [k for k, v in header.items()
                   if v == 'MonteCarlo_Temperature'][0]
            temps = [float(ast.literal_eval(line)[key]) for line in fh]
        shutil.rmtree(outdir)
        # Frames 2 and 5 are sampled at the high temperature
        self.assertEqual(len(temps), 4)
        for t in temps:
            self.assertAlmostEqual(t, 1.0, delta=1e-6)


if __name__ == '__main__':
    IMP.test.main()
