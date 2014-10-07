import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.io.input
import IMP.pmi.plotting.crosslinking
import IMP.test
import RMF
import IMP.rmf
import os,sys

class PlottingTests(IMP.test.TestCase):
    def test_graphxl(self):
        dd={"med14-NTD":[(1,711,"med14")],
            "med14-CTD":[(712,1082,"med14")],
            "med16":["med16"],
            "med2":["med2"],
            "med3":["med3"],
            "med5":["med5"],
            "med15":["med15"]}

        """
        dd={"med6":["med6"],
            "med8":["med8"],
            "med11":["med11"],
            "med17":["med17"],
            "med18":["med18"],
            "med20":["med20"],
            "med22":["med22"],
            "med4":["med4"],
            "med7":["med7"],
            "med9":["med9"],
            "med31":["med31"],
            "med21":["med21"],
            "med10":["med10"],
            "med1":["med1"],
            "med14-NTD":[(1,711,"med14")],
            "med14-CTD":[(712,1000,"med14")],
            "med19":["med19"],
            "med2":["med2"],
            "med3":["med3"],
            "med5":["med5"],
            "med15":["med15"],
            "med16":["med16"]}
        """
        g = IMP.pmi.plotting.crosslinking.GraphXL(IMP.Model(),dd,50.0)
        g.add_rmf('0.rmf3',0)
        g.make_plot('out.png')

if __name__ == '__main__':
    IMP.test.main()
