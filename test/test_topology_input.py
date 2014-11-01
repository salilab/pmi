import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.topology.topology_io
import IMP.pmi.macros

topology_file="/flute1/home/saltzberg/swr/imp/modules/pmi/data/topology.txt"

t=IMP.pmi.topology.topology_io.TopologyReader(topology_file)

for c in t.component_list:
    print c.name, c.pdb_file, c.residue_range, c.gmm_file
