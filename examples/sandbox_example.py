import IMP
import IMP.pmi
import IMP.pmi.sandbox

# the job of PMI is to create and maintain this hierarchy:
# System: stores all the states
# State: stores components (and their copies)
# component: hierarchy decorated as Copy, stores a set of:
# Fragment/Representation an arbitrary subset of a component, has its own
# set of resolutions


# Notes:
# MUST be able to create this system using a table (or a set of them)
# should be sequential - you can't recreate things with the same name or sequence
# also get an error if bodies intersect
# super rigid bodies
# you can make copy index
# to allow symmetry: all PMI hierarchies use child order attribute
# allow guess of number of Gaussians, also an option to do per-residue (R) or per-atom (A)
# always: N gaussians per rigid body
# automatically add beads, flag to NOT do it (leave adding back in for later)

system = IMP.pmi.sandbox.System(model_gaussian_dir='model_gaussians/',
                                target_gaussian_dir='target_gaussians/')
system.add_states(num_states=1)

# BIG QUESTION: can we edit a hierarchy? create a single fragment when you load it, then create new ones? weird.
# alternative: temporary handles for "registering fragments"
# in that case you want to be careful that PMI doesn't screw it up. that's why the graph is nice.
# so what if we create a little object with a clear relationship to the
# hierarchy, and create it only ONCE (also reversible).


# registers a contiguous sequence
# returns a "handle" to the whole sequence (a kind of symbolic selection, can access with residue range etc, see modeller)
# actually this should return a "sequence" object amenable to selection (not decorator, just a little class w/ selection)
# just check how you can select with modeller.sequence
cA = system.add_component(state_num=0,
                          component_name="1oel_A",
                          fasta_fn='data/1oel_A.fasta')

cB = system.add_component(state_num=0,
                          component_name="1oel_B",
                          fasta_fn='data/1oel_B.fasta')

# parses PDB files just to get the init coordinates - store 'em in the sequence!
# should also be able to call component.set_coordinates()
system.set_coordinates(component_name='1oel_A',
                       state_num=0,
                       copy_num=0,
                       pdb_fn='data/1oel.pdb',
                       chain='A')

system.set_coordinates(component_name='1oel_B',
                       state_num=0,
                       copy_num=0,
                       pdb_fn='data/1oel.pdb',
                       chain='B')

# make some "selections" (these should be like python sets) (can also just
# use python indexing)
cA_s1 = cA.residue_range(1, 95)
cA_s2 = cA.residue_range(96, 100)
cA_s3 = cA.residue_range(101, 200)

# create a "body" from parts of A and B. not necessarily rigid, but it's a logical unit.
# key: allows you to set up parts without worrying about what's in your PDB file, what your resolutions are, etc
# residues within the bodies will be modeled at the same set of resolutions
bodyAB = system.create_body(rigid_components=[cA_s1, cB],
                            flexible_components=cA_s2,
                            resolutions=[1, 10],
                            unstructured_resolution=10,
                            representation_type='Gaussian')
bodyA = system.create_body(rigid_components=cA_s3,
                           resolutions=[1, 10],
                           unstructured_resolution=10,
                           representation_type='Gaussian')


# create() loops through the bodies and creates all resolutions (fitting Gaussians if needed)
# question: what about sequences or coordinates that aren't in bodies?
system.create()
