from __future__ import print_function
import time
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

        times = []
        for i in range(5):
            o = s.get_output()
            times.append(float(o['Stopwatch_None_elapsed_seconds']))
            time.sleep(0.01)
        # times should be cumulative
        self.assertGreater(times[-1], times[0])

if __name__ == '__main__':
    IMP.test.main()
