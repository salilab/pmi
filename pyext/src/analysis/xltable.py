#!/usr/bin/env python

"""@namespace IMP.pmi.analysis
   A class for displaying crosslink data
"""

import analysis
import IMP
import IMP.algebra
import IMP.em
import IMP.pmi
import IMP.pmi.tools
import IMP.pmi.output
import IMP.rmf
import RMF
from operator import itemgetter
from copy import deepcopy
from math import log,sqrt
import itertools
import numpy as np


class CrossLinkTable(object):

    def __init__(self):
        self.crosslinks = []
        self.external_csv_data = None
        self.crosslinkedprots = set()
        self.mindist = +10000000.0
        self.maxdist = -10000000.0
        self.contactmap = None

    def set_hierarchy(self, prot):
        self.prot_length_dict = {}
        self.model=prot.get_model()

        for i in prot.get_children():
            name = i.get_name()
            residue_indexes = []
            for p in IMP.atom.get_leaves(i):
                residue_indexes += IMP.pmi.tools.get_residue_indexes(p)

            if len(residue_indexes) != 0:
                self.prot_length_dict[name] = max(residue_indexes)

    def set_coordinates_for_contact_map(self, rmf_name,rmf_frame_index):
        from scipy.spatial.distance import cdist

        rh= RMF.open_rmf_file_read_only(rmf_name)
        prots=IMP.rmf.create_hierarchies(rh, self.model)
        IMP.rmf.load_frame(rh, rmf_frame_index)
        print "getting coordinates for frame %i rmf file %s" % (rmf_frame_index, rmf_name)
        del rh


        coords = []
        radii = []
        namelist = []

        particles_dictionary = analysis.get_particles_at_resolution_one(prots[0])

        resindex = 0
        self.index_dictionary = {}

        for name in particles_dictionary:
            residue_indexes = []
            for p in particles_dictionary[name]:
                print p.get_name()
                residue_indexes = IMP.pmi.tools.get_residue_indexes(p)
                #residue_indexes.add( )

                if len(residue_indexes) != 0:

                    for res in range(min(residue_indexes), max(residue_indexes) + 1):
                        d = IMP.core.XYZR(p)

                        crd = np.array([d.get_x(), d.get_y(), d.get_z()])
                        coords.append(crd)
                        radii.append(d.get_radius())
                        if name not in self.index_dictionary:
                            self.index_dictionary[name] = [resindex]
                        else:
                            self.index_dictionary[name].append(resindex)
                        resindex += 1

        coords = np.array(coords)
        radii = np.array(radii)

        distances = cdist(coords, coords)
        distances = (distances - radii).T - radii

        distances = np.where(distances <= 20.0, 1.0, 0)
        if self.contactmap is None:
            self.contactmap = np.zeros((len(coords), len(coords)))
        self.contactmap += distances

        for prot in prots: IMP.atom.destroy(prot)

    def set_crosslinks(
        self, data_file, search_label='ISDCrossLinkMS_Distance_',
        mapping=None,
        filter_label=None,
        filter_rmf_file_names=None, #provide a list of rmf base names to filter the stat file
        external_csv_data_file=None,
        external_csv_data_file_unique_id_key="Unique ID"):

        # example key ISDCrossLinkMS_Distance_intrarb_937-State:0-108:RPS3_55:RPS30-1-1-0.1_None
        # mapping is a dictionary that maps standard keywords to entry positions in the key string
        # confidence class is a filter that
        # external datafile is a datafile that contains further information on the crosslinks
        # it will use the unique id to create the dictionary keys

        po = IMP.pmi.output.ProcessOutput(data_file)
        keys = po.get_keys()

        xl_keys = [k for k in keys if search_label in k]

        if filter_rmf_file_names is not None:
            rmf_file_key="local_rmf_file_name"
            fs = po.get_fields(xl_keys+[rmf_file_key])
        else:
            fs = po.get_fields(xl_keys)

        # this dictionary stores the occurency of given crosslinks
        self.cross_link_frequency = {}

        # this dictionary stores the series of distances for given crosslinked
        # residues
        self.cross_link_distances = {}

        # this dictionary stores the series of distances for given crosslinked
        # residues
        self.cross_link_distances_unique = {}

        if not external_csv_data_file is None:
            # this dictionary stores the further information on crosslinks
            # labeled by unique ID
            self.external_csv_data = {}
            xldb = IMP.pmi.tools.get_db_from_csv(external_csv_data_file)

            for xl in xldb:
                self.external_csv_data[
                    xl[external_csv_data_file_unique_id_key]] = xl

        # this list keeps track the tuple of cross-links and sample
        # so that we don't count twice the same crosslinked residues in the
        # same sample
        cross_link_frequency_list = []

        self.unique_cross_link_list = []

        for key in xl_keys:
            print key
            keysplit = key.replace(
                "_",
                " ").replace(
                "-",
                " ").replace(
                ":",
                " ").split(
            )

            if filter_label!=None:
                if filter_label not in keysplit: continue

            if mapping is None:
                r1 = int(keysplit[5])
                c1 = keysplit[6]
                r2 = int(keysplit[7])
                c2 = keysplit[8]
                try:
                    confidence = keysplit[12]
                except:
                    confidence = '0.0'
                try:
                    unique_identifier = keysplit[3]
                except:
                    unique_identifier = '0'
            else:
                r1 = int(keysplit[mapping["Residue1"]])
                c1 = keysplit[mapping["Protein1"]]
                r2 = int(keysplit[mapping["Residue2"]])
                c2 = keysplit[mapping["Protein2"]]
                try:
                    confidence = keysplit[mapping["Confidence"]]
                except:
                    confidence = '0.0'
                try:
                    unique_identifier = keysplit[mapping["Unique Identifier"]]
                except:
                    unique_identifier = '0'

            self.crosslinkedprots.add(c1)
            self.crosslinkedprots.add(c2)

            # construct the list of distances corresponding to the input rmf
            # files

            dists=[]
            if filter_rmf_file_names is not None:
                for n,d in enumerate(fs[key]):
                    if fs[rmf_file_key][n] in filter_rmf_file_names:
                        dists.append(float(d))
            else:
                dists=[float(f) for f in fs[key]]

            # check if the input confidence class corresponds to the
            # one of the cross-link

            mdist = self.median(dists)

            stdv = np.std(np.array(dists))
            if self.mindist > mdist:
                self.mindist = mdist
            if self.maxdist < mdist:
                self.maxdist = mdist

            # calculate the frequency of unique crosslinks within the same
            # sample
            if not self.external_csv_data is None:
                sample = self.external_csv_data[unique_identifier]["Sample"]
            else:
                sample = "None"

            if (r1, c1, r2, c2,mdist) not in cross_link_frequency_list:
                if (r1, c1, r2, c2) not in self.cross_link_frequency:
                    self.cross_link_frequency[(r1, c1, r2, c2)] = 1
                    self.cross_link_frequency[(r2, c2, r1, c1)] = 1
                else:
                    self.cross_link_frequency[(r2, c2, r1, c1)] += 1
                    self.cross_link_frequency[(r1, c1, r2, c2)] += 1
                cross_link_frequency_list.append((r1, c1, r2, c2))
                cross_link_frequency_list.append((r2, c2, r1, c1))
                self.unique_cross_link_list.append(
                    (r1, c1, r2, c2,mdist))

            if (r1, c1, r2, c2) not in self.cross_link_distances:
                self.cross_link_distances[(
                    r1,
                    c1,
                    r2,
                    c2,
                    mdist,
                    confidence)] = dists
                self.cross_link_distances[(
                    r2,
                    c2,
                    r1,
                    c1,
                    mdist,
                    confidence)] = dists
                self.cross_link_distances_unique[(r1, c1, r2, c2)] = dists
            else:
                self.cross_link_distances[(
                    r2,
                    c2,
                    r1,
                    c1,
                    mdist,
                    confidence)] += dists
                self.cross_link_distances[(
                    r1,
                    c1,
                    r2,
                    c2,
                    mdist,
                    confidence)] += dists

            self.crosslinks.append(
                (r1,
                 c1,
                 r2,
                 c2,
                 mdist,
                 stdv,
                 confidence,
                 unique_identifier,
                 'original'))
            self.crosslinks.append(
                (r2,
                 c2,
                 r1,
                 c1,
                 mdist,
                 stdv,
                 confidence,
                 unique_identifier,
                 'reversed'))

        self.cross_link_frequency_inverted = {}
        for xl in self.unique_cross_link_list:
            (r1, c1, r2, c2, mdist) = xl
            frequency = self.cross_link_frequency[(r1, c1, r2, c2)]
            if frequency not in self.cross_link_frequency_inverted:
                self.cross_link_frequency_inverted[
                    frequency] = [(r1, c1, r2, c2)]
            else:
                self.cross_link_frequency_inverted[
                    frequency].append((r1, c1, r2, c2))

        # -------------

    def median(self, mylist):
        sorts = sorted(mylist)
        length = len(sorts)
        print length
        if length == 1:
            return mylist[0]
        if not length % 2:
            return (sorts[length / 2] + sorts[length / 2 - 1]) / 2.0
        return sorts[length / 2]

    def set_threshold(self,threshold):
        self.threshold=threshold

    def set_tolerance(self,tolerance):
        self.tolerance=tolerance

    def colormap(self, dist):
        if dist < self.threshold - self.tolerance:
            return "Green"
        elif dist >= self.threshold + self.tolerance:
            return "Orange"
        else:
            return "Red"

    def write_cross_link_database(self, filename, format='csv'):
        import csv

        fieldnames = [
            "Unique ID", "Protein1", "Residue1", "Protein2", "Residue2",
            "Median Distance", "Standard Deviation", "Confidence", "Frequency", "Arrangement"]

        if not self.external_csv_data is None:
            keys = self.external_csv_data.keys()
            innerkeys = self.external_csv_data[keys[0]].keys()
            innerkeys.sort()
            fieldnames += innerkeys

        dw = csv.DictWriter(
            open(filename,
                 "w"),
            delimiter=',',
            fieldnames=fieldnames)
        dw.writeheader()
        for xl in self.crosslinks:
            (r1, c1, r2, c2, mdist, stdv, confidence,
             unique_identifier, descriptor) = xl
            if descriptor == 'original':
                outdict = {}
                outdict["Unique ID"] = unique_identifier
                outdict["Protein1"] = c1
                outdict["Protein2"] = c2
                outdict["Residue1"] = r1
                outdict["Residue2"] = r2
                outdict["Median Distance"] = mdist
                outdict["Standard Deviation"] = stdv
                outdict["Confidence"] = confidence
                outdict["Frequency"] = self.cross_link_frequency[
                    (r1, c1, r2, c2)]
                if c1 == c2:
                    arrangement = "Intra"
                else:
                    arrangement = "Inter"
                outdict["Arrangement"] = arrangement
                if not self.external_csv_data is None:
                    outdict.update(self.external_csv_data[unique_identifier])

                dw.writerow(outdict)

    def plot(self, prot_listx=None, prot_listy=None, no_dist_info=False,
             no_confidence_info=False, filter=None, layout="whole", crosslinkedonly=False,
             filename=None, confidence_classes=None, alphablend=0.1, scale_symbol_size=1.0,
             gap_between_components=0,
             rename_protein_map=None):
        # layout can be:
        #                "lowerdiagonal"  print only the lower diagonal plot
        #                "upperdiagonal"  print only the upper diagonal plot
        #                "whole"  print all
        # crosslinkedonly: plot only components that have crosslinks
        # no_dist_info: if True will plot only the cross-links as grey spots
        # filter = tuple the tuple contains a keyword to be search in the database
        #                a relationship ">","==","<"
        #                and a value
        #                example ("ID_Score",">",40)
        # scale_symbol_size rescale the symbol for the crosslink
        # rename_protein_map is a dictionary to rename proteins

        import matplotlib.pyplot as plt
        import matplotlib.cm as cm

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)

        ax.set_xticks([])
        ax.set_yticks([])

        # set the list of proteins on the x axis
        if prot_listx is None:
            if crosslinkedonly:
                prot_listx = list(self.crosslinkedprots)
            else:
                prot_listx = self.prot_length_dict.keys()
            prot_listx.sort()

        nresx = gap_between_components + \
            sum([self.prot_length_dict[name]
                + gap_between_components for name in prot_listx])

        # set the list of proteins on the y axis

        if prot_listy is None:
            if crosslinkedonly:
                prot_listy = list(self.crosslinkedprots)
            else:
                prot_listy = self.prot_length_dict.keys()
            prot_listy.sort()

        nresy = gap_between_components + \
            sum([self.prot_length_dict[name]
                + gap_between_components for name in prot_listy])

        # this is the residue offset for each protein
        resoffsetx = {}
        resendx = {}
        res = gap_between_components
        for prot in prot_listx:
            resoffsetx[prot] = res
            res += self.prot_length_dict[prot]
            resendx[prot] = res
            res += gap_between_components

        resoffsety = {}
        resendy = {}
        res = gap_between_components
        for prot in prot_listy:
            resoffsety[prot] = res
            res += self.prot_length_dict[prot]
            resendy[prot] = res
            res += gap_between_components

        resoffsetdiagonal = {}
        res = gap_between_components
        for prot in IMP.pmi.tools.OrderedSet(prot_listx + prot_listy):
            resoffsetdiagonal[prot] = res
            res += self.prot_length_dict[prot]
            res += gap_between_components

        # plot protein boundaries

        xticks = []
        xlabels = []
        for n, prot in enumerate(prot_listx):
            res = resoffsetx[prot]
            end = resendx[prot]
            for proty in prot_listy:
                resy = resoffsety[proty]
                endy = resendy[proty]
                ax.plot([res, res], [resy, endy], 'k-', lw=0.4)
                ax.plot([end, end], [resy, endy], 'k-', lw=0.4)
            xticks.append((float(res) + float(end)) / 2)
            if rename_protein_map is not None:
                if prot in rename_protein_map:
                    prot=rename_protein_map[prot]
            xlabels.append(prot)

        yticks = []
        ylabels = []
        for n, prot in enumerate(prot_listy):
            res = resoffsety[prot]
            end = resendy[prot]
            for protx in prot_listx:
                resx = resoffsetx[protx]
                endx = resendx[protx]
                ax.plot([resx, endx], [res, res], 'k-', lw=0.4)
                ax.plot([resx, endx], [end, end], 'k-', lw=0.4)
            yticks.append((float(res) + float(end)) / 2)
            if rename_protein_map is not None:
                if prot in rename_protein_map:
                    prot=rename_protein_map[prot]
            ylabels.append(prot)

        # plot the contact map
        print prot_listx, prot_listy

        if not self.contactmap is None:
            import matplotlib.cm as cm
            tmp_array = np.zeros((nresx, nresy))

            for px in prot_listx:
                print px
                for py in prot_listy:
                    print py
                    resx = resoffsety[px]
                    lengx = resendx[px] - 1
                    resy = resoffsety[py]
                    lengy = resendy[py] - 1
                    indexes_x = self.index_dictionary[px]
                    minx = min(indexes_x)
                    maxx = max(indexes_x)
                    indexes_y = self.index_dictionary[py]
                    miny = min(indexes_y)
                    maxy = max(indexes_y)

                    print px, py, minx, maxx, miny, maxy

                    try:
                        tmp_array[
                            resx:lengx,
                            resy:lengy] = self.contactmap[
                            minx:maxx,
                            miny:maxy]
                    except:
                        continue


            ax.imshow(tmp_array,
                      cmap=cm.binary,
                      origin='lower',
                      interpolation='nearest')

        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels, rotation=90)
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels)
        ax.set_xlim(0,nresx)
        ax.set_ylim(0,nresy)


        # set the crosslinks

        already_added_xls = []

        for xl in self.crosslinks:

            (r1, c1, r2, c2, mdist, stdv, confidence,
             unique_identifier, descriptor) = xl

            if confidence_classes is not None:
                if confidence not in confidence_classes:
                    continue

            try:
                pos1 = r1 + resoffsetx[c1]
            except:
                continue
            try:
                pos2 = r2 + resoffsety[c2]
            except:
                continue

            if not filter is None:
                xldb = self.external_csv_data[unique_identifier]
                xldb.update({"Distance": mdist})
                xldb.update({"Distance_stdv": stdv})

                if filter[1] == ">":
                    if float(xldb[filter[0]]) <= float(filter[2]):
                        continue

                if filter[1] == "<":
                    if float(xldb[filter[0]]) >= float(filter[2]):
                        continue

                if filter[1] == "==":
                    if float(xldb[filter[0]]) != float(filter[2]):
                        continue

            # all that below is used for plotting the diagonal
            # when you have a rectangolar plots

            pos_for_diagonal1 = r1 + resoffsetdiagonal[c1]
            pos_for_diagonal2 = r2 + resoffsetdiagonal[c2]

            if layout == 'lowerdiagonal':
                if pos_for_diagonal1 <= pos_for_diagonal2:
                    continue
            if layout == 'upperdiagonal':
                if pos_for_diagonal1 >= pos_for_diagonal2:
                    continue

            already_added_xls.append((r1, c1, r2, c2))

            if not no_confidence_info:
                if confidence == '0.01':
                    markersize = 14 * scale_symbol_size
                elif confidence == '0.05':
                    markersize = 9 * scale_symbol_size
                elif confidence == '0.1':
                    markersize = 6 * scale_symbol_size
                else:
                    markersize = 15 * scale_symbol_size
            else:
                markersize = 5 * scale_symbol_size

            if not no_dist_info:
                color = self.colormap(mdist)
            else:
                color = "gray"

            ax.plot(
                [pos1],
                [pos2],
                'o',
                c=color,
                alpha=alphablend,
                markersize=markersize)



        fig.set_size_inches(0.004 * nresx, 0.004 * nresy)

        [i.set_linewidth(2.0) for i in ax.spines.itervalues()]

        #plt.tight_layout()

        if filename:
            plt.savefig(filename + ".pdf", dpi=300, transparent="False")


        plt.show()

    def get_frequency_statistics(self, prot_list,
                                 prot_list2=None):

        violated_histogram = {}
        satisfied_histogram = {}
        unique_cross_links = []

        for xl in self.unique_cross_link_list:
            (r1, c1, r2, c2, mdist) = xl

            # here we filter by the protein
            if prot_list2 is None:
                if not c1 in prot_list:
                    continue
                if not c2 in prot_list:
                    continue
            else:
                if c1 in prot_list and c2 in prot_list2:
                    pass
                elif c1 in prot_list2 and c2 in prot_list:
                    pass
                else:
                    continue

            frequency = self.cross_link_frequency[(r1, c1, r2, c2)]

            if (r1, c1, r2, c2) not in unique_cross_links:
                if mdist > 35.0:
                    if frequency not in violated_histogram:
                        violated_histogram[frequency] = 1
                    else:
                        violated_histogram[frequency] += 1
                else:
                    if frequency not in satisfied_histogram:
                        satisfied_histogram[frequency] = 1
                    else:
                        satisfied_histogram[frequency] += 1
                unique_cross_links.append((r1, c1, r2, c2))
                unique_cross_links.append((r2, c2, r1, c1))

        print "# satisfied"

        total_number_of_crosslinks=0

        for i in satisfied_histogram:
            # if i in violated_histogram:
            #   print i, satisfied_histogram[i]+violated_histogram[i]
            # else:
            if i in violated_histogram:
                print i, violated_histogram[i]+satisfied_histogram[i]
            else:
                print i, satisfied_histogram[i]
            total_number_of_crosslinks+=i*satisfied_histogram[i]

        print "# violated"

        for i in violated_histogram:
            print i, violated_histogram[i]
            total_number_of_crosslinks+=i*violated_histogram[i]

        print total_number_of_crosslinks


# ------------
    def print_cross_link_binary_symbols(self, prot_list,
                                        prot_list2=None):
        tmp_matrix = []
        confidence_list = []
        for xl in self.crosslinks:
            (r1, c1, r2, c2, mdist, stdv, confidence,
             unique_identifier, descriptor) = xl

            if prot_list2 is None:
                if not c1 in prot_list:
                    continue
                if not c2 in prot_list:
                    continue
            else:
                if c1 in prot_list and c2 in prot_list2:
                    pass
                elif c1 in prot_list2 and c2 in prot_list:
                    pass
                else:
                    continue

            if descriptor != "original":
                continue

            confidence_list.append(confidence)

            dists = self.cross_link_distances_unique[(r1, c1, r2, c2)]
            tmp_dist_binary = []
            for d in dists:
                if d < 35:
                    tmp_dist_binary.append(1)
                else:
                    tmp_dist_binary.append(0)
            tmp_matrix.append(tmp_dist_binary)

        matrix = zip(*tmp_matrix)

        satisfied_high_sum = 0
        satisfied_mid_sum = 0
        satisfied_low_sum = 0
        total_satisfied_sum = 0
        for k, m in enumerate(matrix):
            satisfied_high = 0
            total_high = 0
            satisfied_mid = 0
            total_mid = 0
            satisfied_low = 0
            total_low = 0
            total_satisfied = 0
            total = 0
            for n, b in enumerate(m):
                if confidence_list[n] == "0.01":
                    total_high += 1
                    if b == 1:
                        satisfied_high += 1
                        satisfied_high_sum += 1
                elif confidence_list[n] == "0.05":
                    total_mid += 1
                    if b == 1:
                        satisfied_mid += 1
                        satisfied_mid_sum += 1
                elif confidence_list[n] == "0.1":
                    total_low += 1
                    if b == 1:
                        satisfied_low += 1
                        satisfied_low_sum += 1
                if b == 1:
                    total_satisfied += 1
                    total_satisfied_sum += 1
                total += 1
            print k, satisfied_high, total_high
            print k, satisfied_mid, total_mid
            print k, satisfied_low, total_low
            print k, total_satisfied, total
        print float(satisfied_high_sum) / len(matrix)
        print float(satisfied_mid_sum) / len(matrix)
        print float(satisfied_low_sum) / len(matrix)
# ------------

    def get_unique_crosslinks_statistics(self, prot_list,
                                         prot_list2=None):

        print prot_list
        print prot_list2
        satisfied_high = 0
        total_high = 0
        satisfied_mid = 0
        total_mid = 0
        satisfied_low = 0
        total_low = 0
        total = 0
        tmp_matrix = []
        satisfied_string = []
        for xl in self.crosslinks:
            (r1, c1, r2, c2, mdist, stdv, confidence,
             unique_identifier, descriptor) = xl

            if prot_list2 is None:
                if not c1 in prot_list:
                    continue
                if not c2 in prot_list:
                    continue
            else:
                if c1 in prot_list and c2 in prot_list2:
                    pass
                elif c1 in prot_list2 and c2 in prot_list:
                    pass
                else:
                    continue

            if descriptor != "original":
                continue

            total += 1
            if confidence == "0.01":
                total_high += 1
                if mdist <= 35:
                    satisfied_high += 1
            if confidence == "0.05":
                total_mid += 1
                if mdist <= 35:
                    satisfied_mid += 1
            if confidence == "0.1":
                total_low += 1
                if mdist <= 35:
                    satisfied_low += 1
            if mdist <= 35:
                satisfied_string.append(1)
            else:
                satisfied_string.append(0)

            dists = self.cross_link_distances_unique[(r1, c1, r2, c2)]
            tmp_dist_binary = []
            for d in dists:
                if d < 35:
                    tmp_dist_binary.append(1)
                else:
                    tmp_dist_binary.append(0)
            tmp_matrix.append(tmp_dist_binary)

        print "unique satisfied_high/total_high", satisfied_high, "/", total_high
        print "unique satisfied_mid/total_mid", satisfied_mid, "/", total_mid
        print "unique satisfied_low/total_low", satisfied_low, "/", total_low
        print "total", total

        matrix = zip(*tmp_matrix)
        satisfied_models = 0
        satstr = ""
        for b in satisfied_string:
            if b == 0:
                satstr += "-"
            if b == 1:
                satstr += "*"

        for m in matrix:
            all_satisfied = True
            string = ""
            for n, b in enumerate(m):
                if b == 0:
                    string += "0"
                if b == 1:
                    string += "1"
                if b == 1 and satisfied_string[n] == 1:
                    pass
                elif b == 1 and satisfied_string[n] == 0:
                    pass
                elif b == 0 and satisfied_string[n] == 0:
                    pass
                elif b == 0 and satisfied_string[n] == 1:
                    all_satisfied = False
            if all_satisfied:
                satisfied_models += 1
            print string
            print satstr, all_satisfied
        print "models that satisfies the median satisfied crosslinks/total models", satisfied_models, len(matrix)

    def plot_matrix_cross_link_distances_unique(self, figurename, prot_list,
                                                prot_list2=None):

        import pylab as pl

        tmp_matrix = []
        for kw in self.cross_link_distances_unique:
            (r1, c1, r2, c2) = kw
            dists = self.cross_link_distances_unique[kw]

            if prot_list2 is None:
                if not c1 in prot_list:
                    continue
                if not c2 in prot_list:
                    continue
            else:
                if c1 in prot_list and c2 in prot_list2:
                    pass
                elif c1 in prot_list2 and c2 in prot_list:
                    pass
                else:
                    continue
            # append the sum of dists to order by that in the matrix plot
            dists.append(sum(dists))
            tmp_matrix.append(dists)

        tmp_matrix.sort(key=itemgetter(len(tmp_matrix[0]) - 1))

        # print len(tmp_matrix),  len(tmp_matrix[0])-1
        matrix = np.zeros((len(tmp_matrix), len(tmp_matrix[0]) - 1))

        for i in range(len(tmp_matrix)):
            for k in range(len(tmp_matrix[i]) - 1):
                matrix[i][k] = tmp_matrix[i][k]

        print matrix

        fig = pl.figure()
        ax = fig.add_subplot(211)

        cax = ax.imshow(matrix, interpolation='nearest')
        # ax.set_yticks(range(len(self.model_list_names)))
        #ax.set_yticklabels( [self.model_list_names[i] for i in leaves_order] )
        fig.colorbar(cax)
        pl.savefig(figurename, dpi=300)
        pl.show()

    def plot_bars(
        self,
        filename,
        prots1,
        prots2,
        nxl_per_row=20,
        arrangement="inter",
            confidence_input="None"):

        data = []
        for xl in self.cross_link_distances:
            (r1, c1, r2, c2, mdist, confidence) = xl
            if c1 in prots1 and c2 in prots2:
                if arrangement == "inter" and c1 == c2:
                    continue
                if arrangement == "intra" and c1 != c2:
                    continue
                if confidence_input == confidence:
                    label = str(c1) + ":" + str(r1) + \
                        "-" + str(c2) + ":" + str(r2)
                    values = self.cross_link_distances[xl]
                    frequency = self.cross_link_frequency[(r1, c1, r2, c2)]
                    data.append((label, values, mdist, frequency))

        sort_by_dist = sorted(data, key=lambda tup: tup[2])
        sort_by_dist = zip(*sort_by_dist)
        values = sort_by_dist[1]
        positions = range(len(values))
        labels = sort_by_dist[0]
        frequencies = map(float, sort_by_dist[3])
        frequencies = [f * 10.0 for f in frequencies]

        nchunks = int(float(len(values)) / nxl_per_row)
        values_chunks = IMP.pmi.tools.chunk_list_into_segments(values, nchunks)
        positions_chunks = IMP.pmi.tools.chunk_list_into_segments(
            positions,
            nchunks)
        frequencies_chunks = IMP.pmi.tools.chunk_list_into_segments(
            frequencies,
            nchunks)
        labels_chunks = IMP.pmi.tools.chunk_list_into_segments(labels, nchunks)

        for n, v in enumerate(values_chunks):
            p = positions_chunks[n]
            f = frequencies_chunks[n]
            l = labels_chunks[n]
            IMP.pmi.output.plot_fields_box_plots(
                filename + "." + str(n), v, p, f,
                valuename="Distance (Ang)", positionname="Unique " + arrangement + " Crosslinks", xlabels=l)

    def crosslink_distance_histogram(self, filename,
                                     prot_list=None,
                                     prot_list2=None,
                                     confidence_classes=None,
                                     bins=40,
                                     color='#66CCCC',
                                     yplotrange=[0, 1],
                                     format="png",
                                     normalized=False):
        if prot_list is None:
            prot_list = self.prot_length_dict.keys()

        distances = []
        for xl in self.crosslinks:
            (r1, c1, r2, c2, mdist, stdv, confidence,
             unique_identifier, descriptor) = xl

            if not confidence_classes is None:
                if confidence not in confidence_classes:
                    continue

            if prot_list2 is None:
                if not c1 in prot_list:
                    continue
                if not c2 in prot_list:
                    continue
            else:
                if c1 in prot_list and c2 in prot_list2:
                    pass
                elif c1 in prot_list2 and c2 in prot_list:
                    pass
                else:
                    continue

            distances.append(mdist)

        IMP.pmi.output.plot_field_histogram(
            filename, distances, valuename="C-alpha C-alpha distance [Ang]",
            bins=bins, color=color,
            format=format,
            reference_xline=35.0,
            yplotrange=yplotrange, normalized=normalized)

    def scatter_plot_xl_features(self, filename,
                                 feature1=None,
                                 feature2=None,
                                 prot_list=None,
                                 prot_list2=None,
                                 yplotrange=None,
                                 reference_ylines=None,
                                 distance_color=True,
                                 format="png"):
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)

        for xl in self.crosslinks:
            (r1, c1, r2, c2, mdist, stdv, confidence,
             unique_identifier, arrangement) = xl

            if prot_list2 is None:
                if not c1 in prot_list:
                    continue
                if not c2 in prot_list:
                    continue
            else:
                if c1 in prot_list and c2 in prot_list2:
                    pass
                elif c1 in prot_list2 and c2 in prot_list:
                    pass
                else:
                    continue

            xldb = self.external_csv_data[unique_identifier]
            xldb.update({"Distance": mdist})
            xldb.update({"Distance_stdv": stdv})

            xvalue = float(xldb[feature1])
            yvalue = float(xldb[feature2])

            if distance_color:
                color = self.colormap(mdist)
            else:
                color = "gray"

            ax.plot([xvalue], [yvalue], 'o', c=color, alpha=0.1, markersize=7)

        if not yplotrange is None:
            ax.set_ylim(yplotrange)
        if not reference_ylines is None:
            for rl in reference_ylines:
                ax.axhline(rl, color='red', linestyle='dashed', linewidth=1)

        if filename:
            plt.savefig(filename + "." + format, dpi=150, transparent="False")

        plt.show()
