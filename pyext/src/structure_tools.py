import IMP
import IMP.atom

def get_structure(mdl,pdb_fn,chain,res_range=None,offset=0):
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
        #res.set_index(res.get_index()+offset)
        ret.append(res)
    return ret
