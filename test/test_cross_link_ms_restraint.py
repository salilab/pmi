import IMP
import IMP.test
import IMP.core
import IMP.container
import IMP.pmi
import IMP.pmi.representation
import IMP.pmi.restraints
import IMP.pmi.restraints.crosslinking
from math import *

class ISDCrossMSTest(IMP.test.TestCase,):

    def setUp(self):
        self.m = IMP.Model()
        IMP.test.TestCase.setUp(self)       
        self.rcomplex=self.init_representation_complex()  
        self.rbeads=self.init_representation_beads()  

    def sphere_cap(self,r1, r2, d):
        sc = 0.0
        if d <= max(r1, r2) - min(r1, r2):
            sc = min(4.0 / 3 * pi * r1 * r1 * r1,
                          4.0 / 3 * pi * r2 * r2 * r2)
        elif d >= r1 + r2 :
            sc = 0
        else:
            sc = (pi / 12 / d * (r1 + r2 - d) * (r1 + r2 - d)) * \
                 (d * d + 2 * d * r1 - 3 * r1 * r1 + 2 * d * r2 + 6 * r1 * r2 -
                  3 * r2 * r2)
        return sc

    def get_probability(self,xyz1s,xyz2s,sigma1s,sigma2s,psi,length,slope):
        onemprob = 1.0
        psi = psi.get_scale()
        for n in range(len(xyz1s)): 
          xyz1=xyz1s[n]
          xyz2=xyz2s[n] 
          sigma1=sigma1s[n]
          sigma2=sigma2s[n]
          dist=IMP.core.get_distance(xyz1, xyz2)

          sigmai = sigma1.get_scale()
          sigmaj = sigma2.get_scale()
          voli = 4.0 / 3.0 * pi * sigmai * sigmai * sigmai
          volj = 4.0 / 3.0 * pi * sigmaj * sigmaj * sigmaj
          fi = 0
          fj = 0
          if dist < sigmai + sigmaj :
              xlvol = 4.0 / 3.0 * pi * (length / 2) * (length / 2) * \
                             (length / 2)
              fi = min(voli, xlvol)
              fj = min(volj, xlvol)
          else:
              di = dist - sigmaj - length / 2
              dj = dist - sigmai - length / 2
              fi = self.sphere_cap(sigmai, length / 2, abs(di))
              fj = self.sphere_cap(sigmaj, length / 2, abs(dj))
          pofr = fi * fj / voli / volj 
          onemprob = onemprob * (1.0 - (psi * (1.0 - pofr) + pofr * (1 - psi))*exp(-slope*dist))
        prob = 1.0 - onemprob
        return prob


    def init_representation_complex(self):
        pdbfile = IMP.pmi.get_data_path("1WCM.pdb")
        fastafile = IMP.pmi.get_data_path("1WCM.fasta.txt")
        components = ["Rpb1","Rpb2","Rpb3","Rpb4"]
        chains = "ABCD"
        colors = [0.,0.1,0.5,1.0]
        beadsize = 20
        fastids = IMP.pmi.tools.get_ids_from_fasta_file(fastafile)
        
        r = IMP.pmi.representation.Representation(self.m)
        hierarchies = {}
        for n in range(len(components)):
            r.create_component(components[n], color=colors[n])
            r.add_component_sequence(components[n], fastafile, id="1WCM:"+chains[n]+"|PDBID|CHAIN|SEQUENCE")
            hierarchies[components[n]] = r.autobuild_model(
                components[n], pdbfile, chains[n],
                resolutions=[1, 10, 100], missingbeadsize=beadsize)
            r.setup_component_sequence_connectivity(components[n], 1)
        return r
    
    def init_representation_beads(self):
        r = IMP.pmi.representation.Representation(self.m)
        r.add_component_name("ProtA",color=1.0)
        r.add_component_beads("ProtA", [(1,10)],incoord=(0,0,0))
        r.add_component_beads("ProtA", [(11,20)],incoord=(10,0,0))    
        r.add_component_beads("ProtA", [(21,30)],incoord=(20,0,0)) 
        r.add_component_name("ProtB",color=1.0)
        r.add_component_beads("ProtB", [(1,10)],incoord=(0,10,0))
        r.add_component_beads("ProtB", [(11,20)],incoord=(10,10,0))    
        r.add_component_beads("ProtB", [(21,30)],incoord=(20,10,0))  
        return r       

    def setup_crosslinks_beads(self,representation,mode):
        
        restraints=IMP.pmi.tools.get_random_cross_link_dataset(representation,
                                                    number_of_cross_links=100,
                                                    resolution=1.0,
                                                    ambiguity_probability=0.2,
                                                    confidence_score_range=[0,100])
        
        ids_map=IMP.pmi.tools.map()
        ids_map.set_map_element(25.0,0.1)        
        ids_map.set_map_element(75,0.01)   

        xl = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(
            representation,
            restraints,
            21,
            label="XL",
            ids_map=ids_map,
            resolution=1,
            inner_slope=0.01)

        sig = xl.get_sigma(1.0)[0]
        psi1 = xl.get_psi(25.0)[0]        
        psi2 = xl.get_psi(75.0)[0]    
        sig.set_scale(10.0)
        psi1.set_scale(0.1)
        psi2.set_scale(0.01)
        
        return xl
        

    def setup_crosslinks_complex(self,representation,mode):
        
        if mode=="single_category":
          columnmap={}
          columnmap["Protein1"]="pep1.accession"
          columnmap["Protein2"]="pep2.accession"
          columnmap["Residue1"]="pep1.xlinked_aa"
          columnmap["Residue2"]="pep2.xlinked_aa"
          columnmap["IDScore"]=None
          columnmap["XLUniqueID"]=None

          ids_map=IMP.pmi.tools.map()
          ids_map.set_map_element(1.0,1.0)

        xl = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representation,
                                   IMP.pmi.get_data_path("polii_xlinks.csv"),
                                   length=21.0,
                                   slope=0.01,
                                   columnmapping=columnmap,
                                   ids_map=ids_map,
                                   resolution=1.0,
                                   label="XL",
                                   csvfile=True)
        xl.add_to_model()
        xl.set_label("XL")
        psi=xl.get_psi(1.0)[0]
        psi.set_scale(0.05)
        sigma=xl.get_sigma(1)[0]
        sigma.set_scale(10.0)
        return xl

   


    def test_score_beads(self):
        self.xlb=self.setup_crosslinks_beads(self.rbeads,"single_category")
      
        cross_link_dict={}
        
        # randomize coordinates
        
        for p in self.xlb.pairs:
            p0 = p[0]
            p1 = p[1]
            prob = p[2].get_probability() 
            resid1 = p[3]
            chain1 = p[4]
            resid2 = p[5]
            chain2 = p[6]
            attribute = p[7]
            xlid=p[11]
            
            d0 = IMP.core.XYZ(p0)
            d1 = IMP.core.XYZ(p1)
            
            sig1 = self.xlb.get_sigma(p[8])[0]
            sig2 = self.xlb.get_sigma(p[9])[0]
            psi = self.xlb.get_psi(p[10])[0]

            if xlid not in cross_link_dict:
               cross_link_dict[xlid]=([d0],[d1],[sig1],[sig2],psi,prob)
            else:
               cross_link_dict[xlid][0].append(d0)
               cross_link_dict[xlid][1].append(d1)               
               cross_link_dict[xlid][2].append(sig1)
               cross_link_dict[xlid][3].append(sig2)
                              
        for xlid in cross_link_dict:

            test_prob=self.get_probability(cross_link_dict[xlid][0],
                                           cross_link_dict[xlid][1],
                                           cross_link_dict[xlid][2],
                                           cross_link_dict[xlid][3],
                                           cross_link_dict[xlid][4],21.0,0.01)
            
            prob=cross_link_dict[xlid][5]
            
            # decrease the tolerance
            
            self.assertAlmostEqual(prob/test_prob,1.0, delta=0.04)

    
    def test_score_complex(self):
        self.xlc=self.setup_crosslinks_complex(self.rcomplex,"single_category")
        
        for p in self.xlc.pairs:
            p0 = p[0]
            p1 = p[1]
            prob = p[2].get_probability() 
            resid1 = p[3]
            chain1 = p[4]
            resid2 = p[5]
            chain2 = p[6]
            attribute = p[7]
            d0 = IMP.core.XYZ(p0)
            d1 = IMP.core.XYZ(p1)
            dist=IMP.core.get_distance(d0, d1)
            
            sig1 = self.xlc.get_sigma(p[8])[0]
            sig2 = self.xlc.get_sigma(p[9])[0]
            psi = self.xlc.get_psi(p[10])[0]
            test_prob=self.get_probability([d0],[d1],[sig1],[sig2],psi,21.0,0.0)
            
            print resid1,resid2,chain1,chain2,attribute,sig1.get_scale(),sig2.get_scale(),psi.get_scale(),prob,test_prob
            self.assertAlmostEqual(prob, test_prob, delta=0.00001)

if __name__ == '__main__':
    IMP.test.main()
