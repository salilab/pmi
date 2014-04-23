## \example pmi/rnapolii/crosslink/model_crosslinks.py

import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi.restraints.stereochemistry as stereochemistry
import IMP.pmi.restraints.crosslinking as crosslinking
import IMP.pmi.representation as representation
import IMP.pmi.tools as tools
import IMP.pmi.samplers as samplers
import IMP.pmi.output as output

IMP.base.set_log_level(IMP.base.SILENT)

# input parameter

pdbfile = IMP.pmi.get_data_path("1WCM.pdb")
fastafile = IMP.pmi.get_data_path("1WCM.fasta.txt")
fastids = tools.get_ids_from_fasta_file(fastafile)
missing_bead_size = 20


#         Component  pdbfile    chainid  rgb color     fastafile     sequence id
# in fastafile
data = [("Rpb1",     pdbfile,   "A",     0.00000000,  (fastafile,    0)),
        ("Rpb2",     pdbfile,   "B",     0.09090909,  (fastafile,    1)),
        ("Rpb3",     pdbfile,   "C",     0.18181818,  (fastafile,    2)),
        ("Rpb4",     pdbfile,   "D",     0.27272727,  (fastafile,    3)),
        ("Rpb5",     pdbfile,   "E",     0.36363636,  (fastafile,    4)),
        ("Rpb6",     pdbfile,   "F",     0.45454545,  (fastafile,    5)),
        ("Rpb7",     pdbfile,   "G",     0.54545455,  (fastafile,    6)),
        ("Rpb8",     pdbfile,   "H",     0.63636364,  (fastafile,    7)),
        ("Rpb9",     pdbfile,   "I",     0.72727273,  (fastafile,    8)),
        ("Rpb10",    pdbfile,   "L",     0.81818182,  (fastafile,    9)),
        ("Rpb11",    pdbfile,   "J",     0.90909091,  (fastafile,   10)),
        ("Rpb12",    pdbfile,   "K",     1.00000000,  (fastafile,   11))]

# mapping the file columns for crosslinks
crosslink_file = 'polii_xlinks.csv'

columnmap = {}
columnmap["Protein1"] = "pep1.accession"
columnmap["Protein2"] = "pep2.accession"
columnmap["Residue1"] = "pep1.xlinked_aa"
columnmap["Residue2"] = "pep2.xlinked_aa"


# create the representation
log_objects = []
optimizable_objects = []

m = IMP.Model()
r = representation.Representation(m)

hierarchies = {}

for d in data:
    component_name = d[0]
    pdb_file = d[1]
    chain_id = d[2]
    color_id = d[3]
    fasta_file = d[4][0]
    fasta_file_id = d[4][1]
    # avoid to add a component with the same name
    r.create_component(component_name,
                       color=color_id)

    r.add_component_sequence(component_name,
                             fasta_file,
                             id=fastids[fasta_file_id])

    hierarchies = r.autobuild_model(component_name,
                                    pdb_file,
                                    chain_id,
                                    resolutions=[1, 10],
                                    missingbeadsize=missing_bead_size)

    r.show_component_table(component_name)

    r.set_rigid_bodies([component_name])

    # r.set_chain_of_super_rigid_bodies(hierarchies,min_length=2,max_length=2)

    r.setup_component_sequence_connectivity(component_name, resolution=1)


# put it at the end of rigid bodies
r.set_floppy_bodies()


r.shuffle_configuration(100, avoidcollision=True)

log_objects.append(r)

ev = stereochemistry.ExcludedVolumeSphere(r, resolution=10)
ev.add_to_model()
log_objects.append(ev)

xl = crosslinking.SigmoidalCrossLinkMS(
    r,
    crosslink_file,
    inflection=25,
    slope=2.0,
    amplitude=5.0,
    linear_slope=0.03,
    resolution=1,
    columnmapping=columnmap,
    csvfile=True)
xl.add_to_model()
log_objects.append(xl)

# move the movers outside ang get the movers from representation and restraints
# helper function to generate samplers samplers.get_movers([r,xl])
# something like: mc = samplers.MonteCarlo(m,samplers.get_movers([r,xl]), 1.0)
mc = samplers.MonteCarlo(m, [r], 1.0)
mc.set_simulated_annealing(min_temp=1.0,
                           max_temp=5.0,
                           min_temp_time=200,
                           max_temp_time=50)
log_objects.append(mc)


o = output.Output()
rmf = o.init_rmf("conformations.rmf3", [r.prot])
# there is a function in IMP: redundant!
# for the pdb use write_pdb_of_c_alphas
# write_pdb_of_c_alphas(IMP.pmi.tools.select(r,resolution=1))
o.add_restraints_to_rmf("conformations.rmf3", [xl])
o.init_stat2("modeling.stat", log_objects)
# remove the frame index
o.write_rmf("conformations.rmf3")

o.write_rmf("conformations.rmf3")


for i in range(2, 100):
    print i
    mc.optimize(10)
    o.write_rmf("conformations.rmf3")
    o.write_stats2()
o.close_rmf("conformations.rmf3")
