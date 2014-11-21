#!/usr/bin/env python

"""@namespace IMP.pmi.analysis
   Clustering tools
"""
import IMP
import IMP.algebra
import IMP.em
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.output
import IMP.rmf
import RMF
import IMP.pmi.analysis
from operator import itemgetter
from copy import deepcopy
from math import log,sqrt
import itertools
import numpy as np


class Alignment(object):
    """Performs alignment and RMSD calculation for two sets of coordinates
    Inputs:

      - query = {'p1':coords(L,3), 'p2':coords(L,3)}
      - template = {'p1':coords(L,3), 'p2':coords(L,3)}

    The class also takes into accout non-equal stoichiometry of the proteins. If this
    is the case, the protein names of protein in multiple copies should be delivered
    in the following form: nameA..1, nameA..2 (note two dots).
    """

    def __init__(self, template, query, weights=None):

        self.query = query
        self.template = template
        self.weights=weights

        if len(self.query.keys()) != len(self.template.keys()):
            raise ValueError('''the number of proteins
                               in template and query does not match!''')

    def permute(self):

        self.proteins = sorted(self.query.keys())
        prots_uniq = [i.split('..')[0] for i in self.proteins]
        P = {}
        for p in prots_uniq:
            np = prots_uniq.count(p)
            copies = [i for i in self.proteins if i.split('..')[0] == p]
            prmts = list(itertools.permutations(copies, len(copies)))
            P[p] = prmts
        self.P = P
        self.Product = list(itertools.product(*P.values()))

    def get_rmsd(self):

        self.permute()

        template_xyz = []
        weights = []
        torder = sum([list(i) for i in self.Product[0]], [])
        for t in torder:
            template_xyz += [IMP.algebra.Vector3D(i) for i in self.template[t]]
            if self.weights is not None:
                weights += [i for i in self.weights[t]]
        #template_xyz = np.array(template_xyz)

        self.rmsd = 10000000000.
        for comb in self.Product:



            order = sum([list(i) for i in comb], [])
            query_xyz = []
            for p in order:
                query_xyz += [IMP.algebra.Vector3D(i) for i in self.query[p]]
            #query_xyz = np.array(query_xyz)
            #if len(template_xyz) != len(query_xyz):
            #    print '''Alignment.get_rmsd: ERROR: the number of coordinates
            #                   in template and query does not match!'''
            #    exit()

            if self.weights is not None:
                dist=IMP.algebra.get_weighted_rmsd(template_xyz, query_xyz, weights)
            else:
                dist=IMP.algebra.get_rmsd(template_xyz, query_xyz)
            #dist = sqrt(
            #    sum(np.diagonal(cdist(template_xyz, query_xyz) ** 2)) / len(template_xyz))
            if dist < self.rmsd:
                self.rmsd = dist
        return self.rmsd

    def align(self):
        from scipy.spatial.distance import cdist

        self.permute()

        template_xyz = []
        torder = sum([list(i) for i in self.Product[0]], [])
        for t in torder:
            template_xyz += [IMP.algebra.Vector3D(i) for i in self.template[t]]
        #template_xyz = np.array(template_xyz)

        self.rmsd, Transformation = 10000000000., ''
        for comb in self.Product:
            order = sum([list(i) for i in comb], [])
            query_xyz = []
            for p in order:
                query_xyz += [IMP.algebra.Vector3D(i) for i in self.query[p]]
            #query_xyz = np.array(query_xyz)

            if len(template_xyz) != len(query_xyz):
                raise ValueError('''the number of coordinates
                               in template and query does not match!''')

            transformation = IMP.algebra.get_transformation_aligning_first_to_second(
                query_xyz,
                template_xyz)
            query_xyz_tr = [transformation.get_transformed(n)
                            for n in query_xyz]

            dist = sqrt(
                sum(np.diagonal(cdist(template_xyz, query_xyz_tr) ** 2)) / len(template_xyz))
            if dist < self.rmsd:
                self.rmsd = dist
                Transformation = transformation

        return (self.rmsd, Transformation)


# TEST for the alignment ###
"""
Proteins = {'a..1':np.array([np.array([-1.,1.])]),
            'a..2':np.array([np.array([1.,1.,])]),
            'a..3':np.array([np.array([-2.,1.])]),
            'b':np.array([np.array([0.,-1.])]),
            'c..1':np.array([np.array([-1.,-1.])]),
            'c..2':np.array([np.array([1.,-1.])]),
            'd':np.array([np.array([0.,0.])]),
            'e':np.array([np.array([0.,1.])])}

Ali = Alignment(Proteins, Proteins)
Ali.permute()
if Ali.get_rmsd() == 0.0: print 'successful test!'
else: print 'ERROR!'; exit()
"""


# ----------------------------------
class Violations(object):

    def __init__(self, filename):

        self.violation_thresholds = {}
        self.violation_counts = {}

        data = open(filename)
        D = data.readlines()
        data.close()

        for d in D:
            d = d.strip().split()
            self.violation_thresholds[d[0]] = float(d[1])

    def get_number_violated_restraints(self, rsts_dict):
        num_violated = 0
        for rst in self.violation_thresholds:
            if rst not in rsts_dict:
                continue  # print rst;
            if float(rsts_dict[rst]) > self.violation_thresholds[rst]:
                num_violated += 1
                if rst not in self.violation_counts:
                    self.violation_counts[rst] = 1
                else:
                    self.violation_counts[rst] += 1
        return num_violated


# ----------------------------------
class Clustering(object):
    """A class to cluster structures.
    Uses scipy's cdist function to compute distance matrices
    And sklearn's kmeans clustering module.
    """
    def __init__(self,rmsd_weights=None):
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            self.rank = comm.Get_rank()
            self.number_of_processes = comm.size
        except ImportError:
            self.number_of_processes = 1
            self.rank = 0
        self.all_coords = {}
        self.structure_cluster_ids = None
        self.tmpl_coords = None
        self.rmsd_weights=rmsd_weights

    def set_template(self, part_coords):

        self.tmpl_coords = part_coords

    def fill(self, frame, Coords):
        """
        fill stores coordinates of a model into a dictionary all_coords,
        containing coordinates for all models.
        """

        self.all_coords[frame] = Coords

    def dist_matrix(self):

        self.model_list_names = self.all_coords.keys()
        self.model_indexes = range(len(self.model_list_names))
        self.model_indexes_dict = dict(
            zip(self.model_list_names, self.model_indexes))
        model_indexes_unique_pairs = list(itertools.combinations(self.model_indexes, 2))

        my_model_indexes_unique_pairs = IMP.pmi.tools.chunk_list_into_segments(
            model_indexes_unique_pairs,
            self.number_of_processes)[self.rank]

        print "process %s assigned with %s pairs" % (str(self.rank), str(len(my_model_indexes_unique_pairs)))

        (raw_distance_dict, self.transformation_distance_dict) = self.matrix_calculation(self.all_coords,
                                                                                         self.tmpl_coords,
                                                                                         my_model_indexes_unique_pairs)

        if self.number_of_processes > 1:
            raw_distance_dict = IMP.pmi.tools.scatter_and_gather(
                raw_distance_dict)
            pickable_transformations = self.get_pickable_transformation_distance_dict(
            )
            pickable_transformations = IMP.pmi.tools.scatter_and_gather(
                pickable_transformations)
            self.set_transformation_distance_dict_from_pickable(
                pickable_transformations)

        self.raw_distance_matrix = np.zeros(
            (len(self.model_list_names), len(self.model_list_names)))
        for item in raw_distance_dict:
            (f1, f2) = item
            self.raw_distance_matrix[f1, f2] = raw_distance_dict[item]
            self.raw_distance_matrix[f2, f1] = raw_distance_dict[item]

    def get_dist_matrix(self):
        return self.raw_distance_matrix

    def do_cluster(self, number_of_clusters,seed=None):
        """Run K-means clustering
        @param number_of_clusters Num means
        @param seed the random seed
        """
        from sklearn.cluster import KMeans
        if seed is not None:
            np.random.seed(seed)
        try:
            # check whether we have the right version of sklearn
            kmeans = KMeans(n_clusters=number_of_clusters)
        except TypeError:
            # sklearn older than 0.12
            kmeans = KMeans(k=number_of_clusters)
        kmeans.fit_predict(self.raw_distance_matrix)

        self.structure_cluster_ids = kmeans.labels_

    def get_pickable_transformation_distance_dict(self):
        pickable_transformations = {}
        for label in self.transformation_distance_dict:
            tr = self.transformation_distance_dict[label]
            trans = tuple(tr.get_translation())
            rot = tuple(tr.get_rotation().get_quaternion())
            pickable_transformations[label] = (rot, trans)
        return pickable_transformations

    def set_transformation_distance_dict_from_pickable(
        self,
            pickable_transformations):
        self.transformation_distance_dict = {}
        for label in pickable_transformations:
            tr = pickable_transformations[label]
            trans = IMP.algebra.Vector3D(tr[1])
            rot = IMP.algebra.Rotation3D(tr[0])
            self.transformation_distance_dict[
                label] = IMP.algebra.Transformation3D(rot, trans)

    def save_distance_matrix_file(self, file_name='cluster.rawmatrix.pkl'):
        import pickle
        outf = open(file_name + ".data", 'w')

        # to pickle the transformation dictionary
        # you have to save the arrays correposnding to
        # the transformations

        pickable_transformations = self.get_pickable_transformation_distance_dict(
        )
        pickle.dump(
            (self.structure_cluster_ids,
             self.model_list_names,
             pickable_transformations),
            outf)

        np.save(file_name + ".npy", self.raw_distance_matrix)

    def load_distance_matrix_file(self, file_name='cluster.rawmatrix.pkl'):
        import pickle

        inputf = open(file_name + ".data", 'r')
        (self.structure_cluster_ids, self.model_list_names,
         pickable_transformations) = pickle.load(inputf)
        inputf.close()

        self.raw_distance_matrix = np.load(file_name + ".npy")

        self.set_transformation_distance_dict_from_pickable(
            pickable_transformations)
        self.model_indexes = range(len(self.model_list_names))
        self.model_indexes_dict = dict(
            zip(self.model_list_names, self.model_indexes))

    def plot_matrix(self, figurename="clustermatrix.pdf"):
        import pylab as pl
        from scipy.cluster import hierarchy as hrc

        fig = pl.figure()
        ax = fig.add_subplot(211)
        dendrogram = hrc.dendrogram(
            hrc.linkage(self.raw_distance_matrix),
            color_threshold=7,
            no_labels=True)
        leaves_order = dendrogram['leaves']

        ax = fig.add_subplot(212)
        cax = ax.imshow(
            self.raw_distance_matrix[leaves_order,
                                     :][:,
                                        leaves_order],
            interpolation='nearest')
        # ax.set_yticks(range(len(self.model_list_names)))
        #ax.set_yticklabels( [self.model_list_names[i] for i in leaves_order] )
        fig.colorbar(cax)
        pl.savefig(figurename, dpi=300)
        #pl.show()

    def get_model_index_from_name(self, name):
        return self.model_indexes_dict[name]

    def get_cluster_labels(self):
        # this list
        return list(set(self.structure_cluster_ids))

    def get_number_of_clusters(self):
        return len(self.get_cluster_labels())

    def get_cluster_label_indexes(self, label):
        return (
            [i for i, l in enumerate(self.structure_cluster_ids) if l == label]
        )

    def get_cluster_label_names(self, label):
        return (
            [self.model_list_names[i]
                for i in self.get_cluster_label_indexes(label)]
        )

    def get_cluster_label_average_rmsd(self, label):

        indexes = self.get_cluster_label_indexes(label)

        if len(indexes) > 1:
            sub_distance_matrix = self.raw_distance_matrix[
                indexes, :][:, indexes]
            average_rmsd = np.sum(sub_distance_matrix) / \
                (len(sub_distance_matrix)
                 ** 2 - len(sub_distance_matrix))
        else:
            average_rmsd = 0.0
        return average_rmsd

    def get_cluster_label_size(self, label):
        return len(self.get_cluster_label_indexes(label))

    def get_transformation_to_first_member(
        self,
        cluster_label,
            structure_index):
        reference = self.get_cluster_label_indexes(cluster_label)[0]
        return self.transformation_distance_dict[(reference, structure_index)]

    def matrix_calculation(self, all_coords, template_coords, list_of_pairs):

        model_list_names = all_coords.keys()
        rmsd_protein_names = all_coords[model_list_names[0]].keys()
        raw_distance_dict = {}
        transformation_distance_dict = {}
        if template_coords is None:
            do_alignment = False
        else:
            do_alignment = True
            alignment_template_protein_names = template_coords.keys()

        for (f1, f2) in list_of_pairs:

            if not do_alignment:
                # here we only get the rmsd,
                # we need that for instance when you want to cluster conformations
                # globally, eg the EM map is a reference
                transformation = IMP.algebra.get_identity_transformation_3d()

                coords_f1 = dict([(pr, all_coords[model_list_names[f1]][pr])
                                 for pr in rmsd_protein_names])
                coords_f2 = {}
                for pr in rmsd_protein_names:
                    coords_f2[pr] = all_coords[model_list_names[f2]][pr]

                Ali = IMP.pmi.analysis.Alignment(coords_f1, coords_f2, self.rmsd_weights)
                rmsd = Ali.get_rmsd()

            elif do_alignment:
                # here we actually align the conformations first
                # and than calculate the rmsd. We need that when the
                # protein(s) is the reference
                coords_f1 = dict([(pr, all_coords[model_list_names[f1]][pr])
                                 for pr in alignment_template_protein_names])
                coords_f2 = dict([(pr, all_coords[model_list_names[f2]][pr])
                                 for pr in alignment_template_protein_names])

                Ali = IMP.pmi.analysis.Alignment(coords_f1, coords_f2)
                template_rmsd, transformation = Ali.align()

                # here we calculate the rmsd
                # we will align two models based n the nuber of subunits provided
                # and transform coordinates of model 2 to model 1
                coords_f1 = dict([(pr, all_coords[model_list_names[f1]][pr])
                                 for pr in rmsd_protein_names])
                coords_f2 = {}
                for pr in rmsd_protein_names:
                    coords_f2[pr] = [transformation.get_transformed(
                        i) for i in all_coords[model_list_names[f2]][pr]]

                Ali = IMP.pmi.analysis.Alignment(coords_f1, coords_f2, self.rmsd_weights)
                rmsd = Ali.get_rmsd()

            raw_distance_dict[(f1, f2)] = rmsd
            raw_distance_dict[(f2, f1)] = rmsd
            transformation_distance_dict[(f1, f2)] = transformation
            transformation_distance_dict[(f2, f1)] = transformation

        return raw_distance_dict, transformation_distance_dict
