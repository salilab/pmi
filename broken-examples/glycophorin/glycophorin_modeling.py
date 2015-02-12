## \example pmi/glycophorin/glycophorin.modeling.py

import IMP
import IMP.core
import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi
import IMP.pmi.restraints.basic
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.representation
import IMP.pmi.macros

import os

def get_data(fname):
    return IMP.pmi.get_example_path("glycophorin/data/" + fname)

log_objects = []
sample_objects = []

m = IMP.Model()

r = IMP.pmi.representation.Representation(m)

r.create_component("chainA",color=0.25)
r.add_component_sequence("chainA", "./data/seq.fasta", id="chainA",offs=61)
a1=[]
for i in range(62,72+1):
    a1+=r.add_component_beads("chainA",[(i,i)])
a2 =r.add_component_ideal_helix("chainA",resolutions=[1],resrange=(73,99))
a3=[]
for i in range(100,101+1):
    a3+=r.add_component_beads("chainA",[(i,i)])

r.setup_component_geometry("chainA")
r.show_component_table("chainA")

r.create_component("chainB",color=0.5)
r.add_component_sequence("chainB","./data/seq.fasta", id="chainB",offs=61)

b1=[]
for i in range(62,72+1):
    b1+=r.add_component_beads("chainB",[(i,i)])
b2=r.add_component_ideal_helix("chainB",resolutions=[1],resrange=(73,99))
b3=[]
for i in range(100,101+1):
    b3+=r.add_component_beads("chainB",[(i,i)])

r.setup_component_geometry("chainB")
r.show_component_table("chainB")

r.set_rigid_body_from_hierarchies(a2)
r.set_rigid_body_from_hierarchies(b2)
r.set_chain_of_super_rigid_bodies([a1,a2,a3],2,5)
r.set_chain_of_super_rigid_bodies([b1,b2,b3],2,5)
r.set_super_rigid_bodies(["chainA"])
r.set_super_rigid_bodies(["chainB"])
r.setup_bonds()
r.set_floppy_bodies()

r.shuffle_configuration(10)
r.optimize_floppy_bodies(1000)

log_objects.append(r)
sample_objects.append(r)

listofexcludedpairs = []

lof = [(62, 75, "chainA"), (96, 101, "chainA"),
       (62, 75, "chainB"), (96, 101, "chainB")]

# add bonds and angles
for l in lof:

    rbr = IMP.pmi.restraints.stereochemistry.ResidueBondRestraint(r, l)
    rbr.add_to_model()
    listofexcludedpairs += rbr.get_excluded_pairs()
    log_objects.append(rbr)

    rar = IMP.pmi.restraints.stereochemistry.ResidueAngleRestraint(r, l)
    rar.add_to_model()
    listofexcludedpairs += rar.get_excluded_pairs()
    log_objects.append(rar)

ev = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(r,resolution=1.0)
ev.add_excluded_particle_pairs(listofexcludedpairs)
ev.add_to_model()
log_objects.append(ev)

eb = IMP.pmi.restraints.basic.ExternalBarrier(r,50)
eb.add_to_model()
log_objects.append(eb)

xl = IMP.pmi.restraints.crosslinking.CysteineCrossLinkRestraint([r],
                    filename="./data/expdata.txt", cbeta=True)
xl.add_to_model()
log_objects.append(xl)
sample_objects.append(xl)

mc=IMP.pmi.macros.ReplicaExchange0(m,r,
                      monte_carlo_sample_objects=sample_objects,
                      output_objects=log_objects,
                      crosslink_restraints=None,
                      monte_carlo_temperature=1.0,
                      replica_exchange_minimum_temperature=1.0,
                      replica_exchange_maximum_temperature=2.5,
                      number_of_best_scoring_models=500,
                      monte_carlo_steps=10,
                      number_of_frames=10000,
                      write_initial_rmf=True,
                      initial_rmf_name_suffix="initial",
                      stat_file_name_suffix="stat",
                      best_pdb_name_suffix="model",
                      do_clean_first=True,
                      do_create_directories=True,
                      global_output_directory="./output/",
                      rmf_dir="/rmfs/",
                      best_pdb_dir="/pdbs/",
                      replica_stat_file_suffix="stat_replica" )

mc.execute_macro()
