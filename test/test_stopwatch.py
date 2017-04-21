from __future__ import print_function
import IMP.test
import IMP.pmi.tools

class Tests(IMP.test.TestCase):
    def test_stopwatch_with_delta(self):
        """Test Stopwatch with delta"""
        s = IMP.pmi.tools.Stopwatch()
        self.assertEqual(s.label, "None")
        s.set_label("foo")
        self.assertEqual(s.label, "foo")
        output = s.get_output()
        self.assertEqual(list(output.keys()), ['Stopwatch_foo_delta_seconds'])
        self.assertIsInstance(list(output.values())[0], str)

    def test_stopwatch_elapsed(self):
        """Test Stopwatch without delta"""
        s = IMP.pmi.tools.Stopwatch(isdelta=False)

        times = [float(s.get_output()['Stopwatch_None_elapsed_seconds'])
                 for _ in range(50)]
        # times should be cumulative
        self.assertGreater(times[-1], times[0])

if __name__ == '__main__':
    IMP.test.main()
