"""@namespace IMP.pmi.input
   Utility classes functions for reading PMI files
"""

import IMP
import IMP.pmi
import IMP.rmf
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.pmi.tools
import RMF
import sys,os
import numpy as np
from collections import defaultdict

def _append_to_dict(d0,d1):
    for k in d0:
        d0[k]+=d1[k]

def save_best_models(mdl,
                     out_dir,
                     stat_files,
                     number_of_best_scoring_models=10,
                     score_key="SimplifiedModel_Total_Score_None",
                     feature_keys=None,
                     rmf_file_key="rmf_file",
                     rmf_file_frame_key="rmf_frame_index",
                     prefiltervalue=None,
                     get_every=1,
                     override_rmf_dir=None):
    """Given a list of stat files, read them all and find the best models.
    Save to a single RMF along with a stat file and a list of the extracted frames.
    """

    # start by splitting into jobs
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        number_of_processes = comm.size
    except ImportError:
        rank = 0
        number_of_processes = 1
    my_stat_files=IMP.pmi.tools.chunk_list_into_segments(
        stat_files,number_of_processes)[rank]

    # filenames
    out_stat_fn = os.path.join(out_dir,"top_"+str(number_of_best_scoring_models)+".out")
    out_rmf_fn = os.path.join(out_dir,"top_"+str(number_of_best_scoring_models)+".rmf3")

    # extract all the models
    all_fields=[]
    for nsf,sf in enumerate(my_stat_files):

        # get list of keywords
        root_directory_of_stat_file = os.path.dirname(os.path.dirname(sf))
        print "getting data from file %s" % sf
        po = IMP.pmi.output.ProcessOutput(sf)
        all_keys = [score_key,
                    rmf_file_key,
                    rmf_file_frame_key]
        for k in po.get_keys():
            for fk in feature_keys:
                if fk in k:
                    all_keys.append(k)

        # optionally filter based on score
        if prefiltervalue is None:
            fields = po.get_fields(all_keys,
                                   get_every=get_every)
        else:
            fields = po.get_fields(all_keys,
                                   filtertuple=(score_key,"<",prefiltervalue),
                                   get_every=get_every)

        # check that all lengths are all equal
        length_set = set([len(fields[f]) for f in fields])
        minlen = min(length_set)

        # if some of the fields are missing, truncate
        # the feature files to the shortest one
        if len(length_set) > 1:
            print "get_best_models: the statfile is not synchronous"
            minlen = min(length_set)
            for f in fields:
                fields[f] = fields[f][0:minlen]
        if nsf==0:
            all_fields=fields
        else:
            for k in fields:
                all_fields[k]+=fields[k]

        if override_rmf_dir is not None:
            for i in range(minlen):
                all_fields[rmf_file_key][i]=os.path.join(
                    override_rmf_dir,os.path.basename(all_fields[rmf_file_key][i]))

    # gather info, sort, write
    comm.Barrier()
    if rank!=0:
        comm.send(all_fields, dest=0, tag=11)
    else:
        for i in range(1,number_of_processes):
            data_tmp = comm.recv(source=i, tag=11)
            for k in all_fields:
                all_fields[k]+=data_tmp[k]

        # sort by total score
        order = sorted(range(len(all_fields[score_key])),
                       key=lambda i: float(all_fields[score_key][i]))

        # write the stat and RMF files
        stat = open(out_stat_fn,'w')
        rh0 = RMF.open_rmf_file_read_only(
            os.path.join(root_directory_of_stat_file,all_fields[rmf_file_key][0]))
        prots = IMP.rmf.create_hierarchies(rh0,mdl)
        del rh0
        outf = RMF.create_rmf_file(out_rmf_fn)
        IMP.rmf.add_hierarchies(outf,prots)
        for i in order[:number_of_best_scoring_models]:
            dline=dict((k,all_fields[k][i]) for k in all_fields)
            rh = RMF.open_rmf_file_read_only(
                os.path.join(root_directory_of_stat_file,all_fields[rmf_file_key][i]))
            IMP.rmf.link_hierarchies(rh,prots)
            IMP.rmf.load_frame(rh,all_fields[rmf_file_frame_key][i])
            IMP.rmf.save_frame(outf)
            del rh
            stat.write(str(dline)+'\n')
        stat.close()
        print 'wrote stats to',out_stat_fn
        print 'wrote rmfs to',out_rmf_fn




def get_best_models(stat_files,
                    score_key="SimplifiedModel_Total_Score_None",
                    feature_keys=None,
                    rmf_file_key="rmf_file",
                    rmf_file_frame_key="rmf_frame_index",
                    prefiltervalue=None,
                    get_every=1):
    """ Given a list of stat files, read them all and find the best models.
    Returns the best rmf filenames, frame numbers, scores, and values for feature keywords
    """
    rmf_file_list=[]              # best RMF files
    rmf_file_frame_list=[]        # best RMF frames
    score_list=[]                 # best scores
    feature_keyword_list_dict=defaultdict(list)  # best values of the feature keys
    for sf in stat_files:
        root_directory_of_stat_file = os.path.dirname(os.path.dirname(sf))
        print "getting data from file %s" % sf
        po = IMP.pmi.output.ProcessOutput(sf)
        keywords = po.get_keys()

        feature_keywords = [score_key,
                            rmf_file_key,
                            rmf_file_frame_key]

        for k in keywords:
            for fk in feature_keys:
                if fk in k:
                    feature_keywords.append(k)

        if prefiltervalue is None:
            fields = po.get_fields(feature_keywords,
                                   get_every=get_every)
        else:
            fields = po.get_fields(feature_keywords,
                                   filtertuple=(score_key,"<",prefiltervalue),
                                   get_every=get_every)

        # check that all lengths are all equal
        length_set = set()
        for f in fields:
            length_set.add(len(fields[f]))

        # if some of the fields are missing, truncate
        # the feature files to the shortest one
        if len(length_set) > 1:
            print "get_best_models: the statfile is not synchronous"
            minlen = min(length_set)
            for f in fields:
                fields[f] = fields[f][0:minlen]

        # append to the lists
        score_list += fields[score_key]
        for rmf in fields[rmf_file_key]:
            rmf_file_list.append(os.path.join(root_directory_of_stat_file,rmf))

        rmf_file_frame_list += fields[rmf_file_frame_key]

        if feature_keywords is not None:
            for k in feature_keywords:
                feature_keyword_list_dict[k] += fields[k]

    return rmf_file_list,rmf_file_frame_list,score_list,feature_keyword_list_dict


def read_coordinates_of_rmfs(model,
                             rmf_tuples,
                             alignment_components=None,
                             rmsd_calculation_components=None):
    """ Read in coordinates of a set of RMF tuples.
    Returns the coordinates split as requested (all, alignment only, rmsd only) as well as
    RMF file names (as keys in a dictionary, with values being the rank number) and just a plain list
    @param model      The IMP model
    @param rmf_tuples [score,filename,frame number,original order number, rank]
    @param alignment_components Tuples to specify what you're aligning on
    @param rmsd_calculation_components Tuples to specify what components are used for RMSD calc
    """
    all_coordinates = []
    rmsd_coordinates = []
    alignment_coordinates = []
    all_rmf_file_names = []
    rmf_file_name_index_dict = {} # storing the features

    for cnt, tpl in enumerate(rmf_tuples):
        rmf_file = tpl[1]
        frame_number = tpl[2]

        prot = IMP.pmi.analysis.get_hier_from_rmf(model,
                                                  frame_number,
                                                  rmf_file)
        if not prot:
            continue

        # getting the particles
        part_dict = IMP.pmi.analysis.get_particles_at_resolution_one(prot)
        all_particles=[pp for key in part_dict for pp in part_dict[key]]
        all_ps_set=set(all_particles)
        # getting the coordinates
        model_coordinate_dict = {}
        template_coordinate_dict={}
        rmsd_coordinate_dict={}
        for pr in part_dict:
            model_coordinate_dict[pr] = np.array(
               [np.array(IMP.core.XYZ(i).get_coordinates()) for i in part_dict[pr]])

        if alignment_components is not None:
            for pr in alignment_components:
                if type(alignment_components[pr]) is str:
                    name=alignment_components[pr]
                    s=IMP.atom.Selection(prot,molecule=name)
                elif type(alignment_components[pr]) is tuple:
                    name=alignment_components[pr][2]
                    rend=alignment_components[pr][1]
                    rbegin=alignment_components[pr][0]
                    s=IMP.atom.Selection(prot,molecule=name,residue_indexes=range(rbegin,rend+1))
                ps=s.get_selected_particles()
                filtered_particles=[p for p in ps if p in all_ps_set]
                template_coordinate_dict[pr] = \
                    [map(float,IMP.core.XYZ(i).get_coordinates()) for i in filtered_particles]

        if rmsd_calculation_components is not None:
            for pr in rmsd_calculation_components:
                if type(rmsd_calculation_components[pr]) is str:
                    name=rmsd_calculation_components[pr]
                    s=IMP.atom.Selection(prot,molecule=name)
                elif type(rmsd_calculation_components[pr]) is tuple:
                    name=rmsd_calculation_components[pr][2]
                    rend=rmsd_calculation_components[pr][1]
                    rbegin=rmsd_calculation_components[pr][0]
                    s=IMP.atom.Selection(prot,molecule=name,residue_indexes=range(rbegin,rend+1))
                ps=s.get_selected_particles()
                filtered_particles=[p for p in ps if p in all_ps_set]
                rmsd_coordinate_dict[pr] = \
                    [map(float,IMP.core.XYZ(i).get_coordinates()) for i in filtered_particles]

        all_coordinates.append(model_coordinate_dict)
        alignment_coordinates.append(template_coordinate_dict)
        rmsd_coordinates.append(rmsd_coordinate_dict)
        frame_name = rmf_file + '|' + str(frame_number)
        all_rmf_file_names.append(frame_name)
        rmf_file_name_index_dict[frame_name] = tpl[4]
    return all_coordinates,alignment_coordinates,rmsd_coordinates,rmf_file_name_index_dict,all_rmf_file_names

def get_bead_sizes(model,rmf_tuple,rmsd_calculation_components=None):
    '''
    @param model      The IMP model
    @param rmf_tuple  score,filename,frame number,original order number, rank
    @param rmsd_calculation_components Tuples to specify what components are used for RMSD calc
    '''
    rmf_file = rmf_tuple[1]
    frame_number = rmf_tuple[2]

    prot = IMP.pmi.analysis.get_hier_from_rmf(model,
                                              frame_number,
                                              rmf_file)

    # getting the particles
    part_dict = IMP.pmi.analysis.get_particles_at_resolution_one(prot)
    all_particles=[pp for key in part_dict for pp in part_dict[key]]
    all_ps_set=set(all_particles)
    # getting the coordinates
    rmsd_bead_size_dict={}

    if rmsd_calculation_components is not None:
        for pr in rmsd_calculation_components:
            if type(rmsd_calculation_components[pr]) is str:
                name=rmsd_calculation_components[pr]
                s=IMP.atom.Selection(prot,molecule=name)
            elif type(rmsd_calculation_components[pr]) is tuple:
                name=rmsd_calculation_components[pr][2]
                rend=rmsd_calculation_components[pr][1]
                rbegin=rmsd_calculation_components[pr][0]
                s=IMP.atom.Selection(prot,molecule=name,residue_indexes=range(rbegin,rend+1))
            ps=s.get_selected_particles()
            filtered_particles=[p for p in ps if p in all_ps_set]
            rmsd_bead_size_dict[pr] = \
                [len(IMP.pmi.tools.get_residue_indexes(p)) for p in filtered_particles]


    return rmsd_bead_size_dict
