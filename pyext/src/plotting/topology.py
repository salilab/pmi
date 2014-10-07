import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
from collections import defaultdict
import numpy as np
from scipy.spatial.distance import cdist

class TopologyPlot(object):
    """A class to read RMF files and make a network contact map"""
    def __init__(self,model,selections,cutoff):
        """Set up a new graphXL object
        @param model          The IMP model
        @param selection_dict A dictionary containing component names.
                              Keys are labels
                              values are either moleculename or start,stop,moleculename
        @param cutoff        The
        """
        self.mdl = model
        self.selections = selections
        self.contact_counts={}
        self.edges=defaultdict(int)
        self.cutoff = cutoff
        self.gcpf = IMP.core.GridClosePairsFinder()
        self.gcpf.set_distance(self.cutoff)
        self.names = self.selections.keys()
        self.num_rmf=0

    def add_rmf(self,rmf_fn,nframe):
        """Add selections from an RMF file"""
        print 'reading from RMF file',rmf_fn
        rh = RMF.open_rmf_file_read_only(rmf_fn)
        prots = IMP.rmf.create_hierarchies(rh, self.mdl)
        hier = prots[0]
        IMP.rmf.load_frame(rh,0)
        ps_per_component=defaultdict(list)
        self.size_per_component=defaultdict(int)
        self.mdl.update()

        #gathers particles for all components
        part_dict = IMP.pmi.analysis.get_particles_at_resolution_one(hier)
        all_particles_by_resolution = []
        for name in part_dict:
            all_particles_by_resolution += part_dict[name]

        for component_name in self.selections:
            for seg in self.selections[component_name]:
                if type(seg) == str:
                    s = IMP.atom.Selection(hier,molecule=seg)
                elif type(seg) == tuple:
                    s = IMP.atom.Selection(hier,molecule=seg[2],
                        residue_indexes=range(seg[0], seg[1] + 1))
                else:
                    raise Exception('could not understand selection tuple '+str(seg))
                parts = list(set(s.get_selected_particles()) & set(all_particles_by_resolution))
                ps_per_component[component_name] += IMP.get_indexes(parts)
                if self.num_rmf==0:
                    self.size_per_component[component_name] += sum(len(IMP.pmi.tools.get_residue_indexes(p)) for p in parts)

        for n1,name1 in enumerate(self.names):
            for name2 in self.names[n1+1:]:
                ncontacts = len(self.gcpf.get_close_pairs(self.mdl,
                                                     ps_per_component[name1],
                                                     ps_per_component[name2]))
                if ncontacts>0:
                    self.edges[tuple(sorted((name1,name2)))]+=ncontacts
        self.num_rmf+=1

    def make_plot(self,out_fn):
        edges=[]
        weights=[]
        print 'num edges',len(self.edges)
        for edge,count in self.edges.iteritems():
            edges.append(edge)
            weights.append(count)
        for nw,w in enumerate(weights):
            weights[nw]=float(weights[nw])/max(weights)
        IMP.pmi.output.draw_graph(edges,#node_size=1000,
                                  node_size=dict(self.size_per_component),
                                  edge_thickness=1, #weights,
                                  edge_alpha=0.3,
                                  out_filename=out_fn)
