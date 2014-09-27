"""@namespace IMP.pmi.input
   Utility classes functions for reading PMI files
"""

import IMP
import IMP.pmi
import IMP.rmf
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.pmi.tools
import glob
import sys,os
import numpy as np
from collections import defaultdict

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
    """ Read in coordinates of a set of RMF tuples
    @param model      The IMP model
    @param rmf_tuples [score,filename,frame number,original order number, rank]
    @param alignment_components Tuples to specify what you're aligning on
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
            filtered_particles=list(set(ps)&set(all_particles))
            template_coordinate_dict[pr] = np.array(
                [np.array(IMP.core.XYZ(i).get_coordinates()) for i in filtered_particles])

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
            filtered_particles=list(set(ps)&set(all_particles))
            rmsd_coordinate_dict[pr] = np.array(
                [np.array(IMP.core.XYZ(i).get_coordinates()) for i in filtered_particles])

        all_coordinates.append(model_coordinate_dict)
        alignment_coordinates.append(template_coordinate_dict)
        rmsd_coordinates.append(rmsd_coordinate_dict)
        frame_name = rmf_file + '|' + str(frame_number)
        all_rmf_file_names.append(frame_name)
        rmf_file_name_index_dict[frame_name] = tpl[4]
    return all_coordinates,alignment_coordinates,rmsd_coordinates,rmf_file_name_index_dict,all_rmf_file_names
