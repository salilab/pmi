import IMP
import IMP.atom

def get_structure(mdl,pdb_fn,chain,res_range=None,offset=0):
    """return the residues from a pdb file with particular chain and range"""
    mh=IMP.atom.read_pdb(pdb_fn,mdl,IMP.atom.NonWaterNonHydrogenPDBSelector())

    # first update using offset:
    for rr in IMP.atom.get_by_type(mh,IMP.atom.RESIDUE_TYPE):
        IMP.atom.Residue(rr).set_index(IMP.atom.Residue(rr).get_index()+offset)

    # get requested chain and residue range
    if res_range is not None:
        res_range=range(res_range[0],res_range[1]+1)
    sel=IMP.atom.Selection(mh,chain=chain,residue_indexes=res_range,
                               atom_type=IMP.atom.AT_CA)
    ret=[]
    for p in sel.get_selected_particles():
        res=IMP.atom.Residue(IMP.atom.Atom(p).get_parent())
        ret.append(res)
    return ret

def fill_in_missing_backbone(mdl,residues):
    """Guess CA position based on surroundings.
    Only does this if the Residue has associated Representations
    for each residue:
      if residue.representations != defaultdict(set)
        create a new Residue and CAlpha and guess position
        store as residue.res
    """
    pass


def build_along_backbone(mdl,root,residues,rep_type,ca_centers=True):
    """Group residues along the backbone, adding beads as needed.
    Currently this first groups into contiguous fragment regions ('folders')
    with identical sets of resolutions. However this behavior may
    The current resolutions used are 0 for atomic, and N for N residues per ball.
    @param mdl        the model
    @param root       the hierarchy to which all the fragments and resolutions will be added
    @param residues   list of PMI Residues, labeled with resolution
    @param rep_type   Representation type (currently supports IMP.atom.BALLS)
    @param ca_centers If true, when making residue beads, will place the bead at the CA position
    """
    allowed_reps=[IMP.atom.BALLS]
    if rep_type not in allowed_reps:
        print "Only supported representation types are",allowed_types
    prev_rep = None
    cur_fragment=[]
    fragments=[]

    # group into fragments with identical resolutions
    for res in residues:
        if prev_rep is None:
            prev_rep = res.representations
        rep = res.representations
        if rep==prev_rep:
            cur_fragment.append(res)
        else:
            fragments.append(cur_fragment)
            cur_fragment=[res]
        prev_rep=rep
    fragments.append(cur_fragment)

    # build the representations within each fragment
    for frag_res in fragments:
        res_nums=[r.num for r in frag_res]
        this_rep=frag_res[0].representations
        if len(this_rep)==0:
            continue
        print 'building frag',frag_res,'pdb nums',res_nums
        frag = IMP.atom.Fragment.setup_particle(mdl,mdl.add_particle("fragment root"),res_nums)
        root.add_child(frag)
        frep = IMP.atom.Representation.setup_particle(frag,0)

        if 'balls' in this_rep:
            # if atomic, add the residues as child (inside another fragment - check RMF)
            if 0 in this_rep['balls']:
                f=IMP.atom.Fragment.setup_particle(mdl,mdl.add_particle("resolution 0"),
                                                   res_nums)
                for residue in frag_res:
                    #frag.add_child(residue.res)
                    f.add_child(residue.res)
                frag.add_child(f)

            # add one-residue-per-bead
            if 1 in this_rep['balls']:
                res1=IMP.atom.Fragment.setup_particle(mdl,mdl.add_particle("resolution 1"),res_nums)
                for residue in frag_res:
                    if ca_centers==True:
                        rp1 = IMP.Particle(mdl)
                        rp1.set_name("res1_idx%i"%residue.num)
                        rt=residue.res.get_residue_type()
                        r1 = IMP.atom.Residue.setup_particle(rp1,rt,residue.num)
                        res1.add_child(r1)
                        try:
                            vol = IMP.atom.get_volume_from_residue_type(rt)
                        except IMP.base.ValueException:
                            vol = IMP.atom.get_volume_from_residue_type(
                                IMP.atom.ResidueType("ALA"))
                        try:
                            mass = IMP.atom.get_mass_from_residue_type(rt)
                        except:
                            mass = IMP.atom.get_mass_from_residue_type(IMP.atom.ResidueType("ALA"))
                        calpha = IMP.atom.Selection(residue.res,atom_type=IMP.atom.AT_CA). \
                                   get_selected_particles()[0]
                        radius = IMP.algebra.get_ball_radius_from_volume_3d(vol)
                        shape = IMP.algebra.Sphere3D(IMP.core.XYZ(calpha).get_coordinates(),radius)
                        IMP.core.XYZR.setup_particle(rp1,shape)
                        IMP.atom.Mass.setup_particle(rp1,mass)
                frep.add_representation(res1,IMP.atom.BALLS,1)


            # add all other resolutions
            for resolution in set(this_rep['balls']) - set([0,1]):
                pass

def show_representation(node):
    if IMP.atom.Representation.get_is_setup(node):
        repr=IMP.atom.Representation(node)
        resolutions=repr.get_resolutions()
        for r in resolutions:
            print '---- resolution %i ----' %r
            IMP.atom.show_molecular_hierarchy(repr.get_representation(r))
        return True
    else:
        return False

def recursive_show_representations(root):
    pass
