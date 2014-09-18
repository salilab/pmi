import IMP
import IMP.em
import IMP.pmi
import IMP.rmf
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.pmi.macros
import RMF
import numpy as np
from operator import itemgetter
from optparse import OptionParser
import sys,os
from glob import glob
import re

def parse_args():
    usage = """%prog [options] <analysis_option> <input file>
    Currently implemented analysis options: cluster
    The input file contains all the setup options:

    (general)
    subunits
    global_output_dir = "./"

    (cluster)
    feature_keys = None
    outputdir = "kmeans_2_1/"
    number_of_best_scoring_models = 100
    alignment_components = None
    rmsd_calculation_components = components_names
    distance_matrix_file = "distance.rawmatrix.pkl"
    feature_keys = []
    load_distance_matrix_file = False,
    skip_clustering = 0
    display_plot = 0
    exit_after_display = 0
    get_every = 1
    is_mpi = 1
    number_of_clusters = 1
    voxel_size = 3.0
    """
    parser = OptionParser(usage)
    (options, args) = parser.parse_args()
    if len(args) != 2:
        parser.error("incorrect number of arguments")
    return options,args

def str2bool(s):
    if s=='1':
        return True
    else:
        return False

def run():
    # handle input
    options,args = parse_args()
    inf = open(args[1],'r')
    mdl = IMP.Model()
    if args[0]=="cluster":
        cluster_dict={'merge_directories' : './',
                      'global_output_dir' : './',
                      'feature_keys':None,
                      'number_of_best_scoring_models' : 100,
                      'distance_matrix_file' : "distance.rawmatrix.pkl",
                      'outputdir' : "kmeans_2_1/" ,
                      'load_distance_matrix_file' : 0,
                      'skip_clustering' : 0,
                      'display_plot' : 0,
                      'exit_after_display' : 0,
                      'get_every' : 0,
                      'is_mpi' : 1,
                      'number_of_clusters' : 1,
                      'voxel_size' : 3.0}
        new_dict={}
        for l in inf:
            fields = l.split()
            key = fields[0]
            if len(fields)>2:
                new_dict[key]=fields[1:]
            else:
                new_dict[key]=fields[1]
        cluster_dict.update(new_dict)
        cluster_dict['number_of_best_scoring_models']=int(cluster_dict['number_of_best_scoring_models'])
        cluster_dict['load_distance_matrix_file']=str2bool(cluster_dict['load_distance_matrix_file'])
        cluster_dict['skip_clustering']=str2bool(cluster_dict['skip_clustering'])
        cluster_dict['display_plot']=str2bool(cluster_dict['display_plot'])
        cluster_dict['exit_after_display']=str2bool(cluster_dict['exit_after_display'])
        cluster_dict['get_every']=int(cluster_dict['get_every'])
        cluster_dict['is_mpi']=str2bool(cluster_dict['is_mpi'])
        cluster_dict['number_of_clusters']=int(cluster_dict['number_of_clusters'])
        cluster_dict['voxel_size']=float(cluster_dict['voxel_size'])
        cluster_dict['prefilter_value']=float(cluster_dict['prefilter_value'])

        sels={}
        for s in cluster_dict['subunits']:
            sels[s]=[s]


        print '\nRUNNING CLUSTERING WITH THESE OPTIONS'
        for k in cluster_dict:
            print k,':',cluster_dict[k]
        print 'density custom dict',sels

        mc=IMP.pmi.macros.AnalysisReplicaExchange0(mdl,
                                                   stat_file_name_suffix="stat",
                                                   merge_directories=cluster_dict['merge_directories'],
                                                   global_output_directory=cluster_dict['global_output_dir'],
                                                   rmf_dir="rmfs/")

        mc.clustering("SimplifiedModel_Total_Score_None",
                      "rmf_file",
                      "rmf_frame_index",
                      prefiltervalue=cluster_dict['prefilter_value'],
                      number_of_best_scoring_models=cluster_dict['number_of_best_scoring_models'],
                      alignment_components=None,
                      rmsd_calculation_components=cluster_dict['subunits'],
                      distance_matrix_file="distance.rawmatrix.pkl",
                      outputdir=cluster_dict['output_dir'],
                      feature_keys=cluster_dict['feature_keys'],
                      load_distance_matrix_file=cluster_dict['load_distance_matrix_file'],
                      skip_clustering=cluster_dict['skip_clustering'],
                      display_plot=cluster_dict['display_plot'],
                      exit_after_display=cluster_dict['exit_after_display'],
                      get_every=cluster_dict['get_every'],
                      is_mpi=cluster_dict['is_mpi'],
                      number_of_clusters=cluster_dict['number_of_clusters'],
                      voxel_size=cluster_dict['voxel_size'],
                      density_custom_ranges=sels)


if __name__=="__main__":
    run()
