import RMF
import IMP
import IMP.rmf
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.analysis
from optparse import OptionParser
import sys,os

def parse_args():
    usage = """%prog [options] <rmf_fn> <frame number> <out_fn.cmm>
    """
    parser = OptionParser(usage)
    parser.add_option("-t","--threshold",dest="threshold",
                      default=35.0,
                      type=float,
                      help="Above this threshold, use a different color. "
                      "Default is 35A")
    parser.add_option("-r","--radius",dest="radius",
                      default=2.0,
                      type=float,
                      help="Radius for the XL. Default is 2.0")
    parser.add_option("-c","--color",dest="color",
                      default='93,238,93',
                      help="RGB colors for non-violated XL"
                      "format is R,G,B where each value is out of 255"
                      "Default is green")
    parser.add_option("-v","--color_viol",dest="color_viol",
                      default='250,77,63',
                      help="RGB colors for violated XL (above options.threshold)"
                      "format is R,G,B where each value is out of 255"
                      "Default is red")

    (options, args) = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")
    return options,args

def run():
    options,args = parse_args()
    rmf_fn,frame_num,out_fn = args
    try:
        color = map(lambda x: float(x)/255,options.color.split(','))
        vcolor = map(lambda x: float(x)/255,options.color_viol.split(','))
        t = color[2]
        t = vcolor[2]
        threshold = float(options.threshold)
        radius = float(options.radius)
    except:
        raise InputError("Wrong color format")

    outf=open(out_fn,'w')
    outf.write('<marker_set name="xlinks">\n')

    mdl = IMP.Model()
    rh = RMF.open_rmf_file_read_only(rmf_fn)
    prots = IMP.rmf.create_hierarchies(rh, mdl)
    prot=prots[0]
    rs = IMP.rmf.create_restraints(rh, mdl)
    IMP.rmf.load_frame(rh,int(frame_num))
    nv=0
    for r in rs:
        ps = r.get_inputs()
        c1,c2 = [IMP.core.XYZ(IMP.kernel.Particle.get_from(p)).get_coordinates() for p in ps]
        dist = IMP.algebra.get_distance(c1,c2)
        if dist<threshold:
            r,g,b = color
        else:
            r,g,b = vcolor
        outf.write('<marker id= "%d" x="%.3f" y="%.3f" z="%.3f" radius="%.2f" '
                   'r="%.2f" g="%.2f" b="%.2f"/> \n' % (nv,c1[0],c1[1],c1[2],radius,r,g,b))
        outf.write('<marker id= "%d" x="%.3f" y="%.3f" z="%.3f" radius="%.2f" '
                   'r="%.2f" g="%.2f" b="%.2f"/> \n' % (nv+1,c2[0],c2[1],c2[2],radius,r,g,b))
        outf.write('<link id1= "%d" id2="%d" radius="%.2f" '
                   'r="%.2f" g="%.2f" b="%.2f"/> \n' % (nv,nv+1,radius,r,g,b))
        nv+=2
    outf.write('</marker_set>\n')
    outf.close()
    print 'wrote',len(rs),'XLs to',out_fn

if __name__=="__main__":
    run()
