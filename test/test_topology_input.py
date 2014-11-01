import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.topology.topology_io
import IMP.pmi.macros

topology_file="/flute1/home/saltzberg/swr/imp/modules/pmi/data/topology.txt"

t=IMP.pmi.topology.topology_io.TopologyReader(topology_file)

for c in t.component_list:
    print c.name, c.pdb_file, c.residue_range, c.gmm_file


m = IMP.Model()
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)
bm=IMP.pmi.macros.BuildModel1(simo)
bm.build_model(component_topologies=t.component_list[6:8])
