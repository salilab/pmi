import IMP
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.macros
import os

topology_file=IMP.pmi.get_data_path("topology.txt")

t=IMP.pmi.topology.TopologyReader(topology_file)

#for c in t.component_list:
#    print c.name, c.pdb_file, c.residue_range, c.gmm_file

topology_out=os.getcwd()+"/topology_out.txt"

t.write_topology_file(topology_out)


t=IMP.pmi.topology.topology_io.Topology(topology_out)


m = IMP.Model()
simo = IMP.pmi.representation.Representation(m,upperharmonic=True,disorderedlength=False)
bm=IMP.pmi.macros.BuildModel1(simo)
bm.build_model(component_topologies=t.component_list[6:8])


os.remove(topology_out)
