#!/usr/bin/env python3
import os

import h5py

from functools import partial
import pybedtools
import json
import numpy as np
import torch
from torch.multiprocessing import Pool, Process
import gc
import time
import fire
from bs4 import BeautifulSoup
from rankagg import FullListRankAggregator


# def decide_cofactor_motif_set(outfilename, n=10, *motif_scan_json_files):
#     list_cofactor = []
#     for motif_scan_json_file in motif_scan_json_files:
#         dic_tf_zscore = {}
#         with open(motif_scan_json_file.replace("mdseqpos_index.html", "motif_list.json"), "r") as infile:
#             for line in infile:
#                 motif_info = json.loads(line.strip())
#                 dic_tf_zscore[motif_info['id']] = motif_info['seqpos_results']['zscore']
#         sorted_tf_zscore = sorted(dic_tf_zscore.items(), key=lambda x: x[1], reverse=False)
#         list_cofactor += list(zip(*sorted_tf_zscore))[0][:n]
#     list_cofactor = set([cofactor for cofactor in list_cofactor if not cofactor.startswith('denovo')])
#     list_cofactor = list(list_cofactor)
#     with open(outfilename, "w") as outfile:
#         json.dump(list_cofactor, outfile)
#
# def decide_cofactor_motif_set(outfilename, n=10, *motif_scan_json_files):
#     list_cofactor = []
#     while 1:
#         for motif_scan_json_file in motif_scan_json_files:
#             with open(motif_scan_json_file.replace("mdseqpos_index.html", "table.html"), "r") as infile:
#                 tree = BeautifulSoup(infile.read(), 'html')
#             motifs = [i('td')[1]('a')[0].string.strip() for i in tree("tbody")[0]('tr') if 'rowspan' in i('td')[0].attrs]
#             motifs = [i for i in motifs if not i.startswith('denovo')]
#             list_cofactor += motifs[:n]
#         list_cofactor = set([cofactor for cofactor in list_cofactor if not cofactor.startswith('denovo')])
#         list_cofactor = list(list_cofactor)
#         if len(list_cofactor) <= 40:
#             break
#         else:
#             n -= 1
#     with open(outfilename, "w") as outfile:
#         json.dump(list_cofactor, outfile)


def decide_cofactor_motif_set(outfilename, n=10, *motif_scan_json_files):
    list_cofactor = []

    for motif_scan_json_file in motif_scan_json_files:
        with open(motif_scan_json_file.replace("mdseqpos_index.html", "table.html"), "r") as infile:
            tree = BeautifulSoup(infile.read(), 'html')
        print(motif_scan_json_file)
        list_cluster_top_index = [ind for ind, i in enumerate(tree("tbody")[0]('tr')) if 'rowspan' in i('td')[0].attrs]
        list_motif_name = []
        for ind, cluster_top_index in enumerate(list_cluster_top_index):
            if len(list_motif_name) >= 10:
                break
            index = cluster_top_index
            while index < list_cluster_top_index[ind + 1]:
                info = tree("tbody")[0]('tr')[index]
                if index == cluster_top_index:
                    motif_name = info('td')[1]('a')[0].string.strip()
                else:
                    motif_name = info('td')[0]('a')[0].string.strip()
                if not motif_name.startswith('denovo'):
                    list_motif_name.append(motif_name)
                    break
                else:
                    index += 1
        list_cofactor += list_motif_name
    list_cofactor = set([cofactor for cofactor in list_cofactor if not cofactor.startswith('denovo')])
    list_cofactor = list(list_cofactor)

    list_dic_motif_name_score = []
    for motif_scan_json_file in motif_scan_json_files:
        # with open(motif_scan_json_file.replace("mdseqpos_index.html", "table.html"), "r") as infile:
        #     tree = BeautifulSoup(infile.read(), 'html')
        # list_cluster_top_index = [ind for ind, i in enumerate(tree("tbody")[0]('tr')) if 'rowspan' in i('td')[0].attrs]
        dic_motif_name_score = {}
        # for ind, info in enumerate(tree("tbody")[0]('tr')):
        #     if ind in list_cluster_top_index:
        #         motif_name = info('td')[1]('a')[0].string.strip()
        #     else:
        #         motif_name = info('td')[0]('a')[0].string.strip()
        #     if motif_name in list_cofactor:
        #         if ind in list_cluster_top_index:
        #             dic_motif_name_score[motif_name] = -float(info('td')[6].string.strip())
        #         else:
        #             dic_motif_name_score[motif_name] = -float(info('td')[5].string.strip())
        dic_tf_zscore = {}
        with open(motif_scan_json_file.replace("mdseqpos_index.html", "motif_list.json"), "r") as infile:
            for line in infile:
                motif_info = json.loads(line.strip())
                dic_tf_zscore[motif_info['id']] = motif_info['seqpos_results']['zscore']
        for motif_name in list_cofactor:
            if motif_name in dic_tf_zscore:
                dic_motif_name_score[motif_name] = -dic_tf_zscore[motif_name]
            else:
                dic_motif_name_score[motif_name] = 0
        list_dic_motif_name_score.append(dic_motif_name_score)

    FLRA = FullListRankAggregator()

    aggRanks = FLRA.aggregate_ranks(list_dic_motif_name_score, method='robust')
    selected_cofactor = [i for i in aggRanks if aggRanks[i] <= n]
    # print(aggRanks)
    with open(outfilename.replace(".json",".motif_agg_rank.json"), "w") as outfile:
        json.dump(aggRanks, outfile)
    with open(outfilename, "w") as outfile:
        json.dump(selected_cofactor, outfile)

class MotifFeatures(object):
    def __init__(self, chrom_all, chrom_size_file, genome_sequence_fa, motif_pwm_path, motif_feature_path=None,
                 batch=10000000):
        self.chrom_all = chrom_all
        self.dic_chrom_length = {}
        with open(chrom_size_file, "r") as infile:
            for line in infile:
                line = line.strip().split("\t")
                if line[0] in chrom_all:
                    self.dic_chrom_length[line[0]] = int(line[1])

        self.genome_sequence_fa = genome_sequence_fa
        self.motif_pwm_path = motif_pwm_path
        if motif_feature_path is None:
            self.motif_feature_path = "./hdf5s"
        else:
            self.motif_feature_path = motif_feature_path

        self.batch = batch

    def read_motif_pwm(self):
        tf_length_max = 0
        motif_names = []
        motif_pwms = []
        # motif_lengths = []
        for motif_file in os.listdir(self.motif_pwm_path):
            tf_matrix = np.loadtxt(self.motif_pwm_path + motif_file, comments=["#", ">"], delimiter='\t',
                                   dtype='float32')
            tf_length = tf_matrix.shape[0]
            # motif_lengths.append(tf_length)
            if tf_length_max < tf_length:
                tf_length_max = tf_length
            if motif_file.endswith('.pwm'):
                tf_name = motif_file.replace('.pwm', '')
            else:
                tf_name = motif_file
            motif_names.append(tf_name)
            motif_pwms.append(tf_matrix)

        for ind, pwm in enumerate(motif_pwms):
            zero_padding_length = tf_length_max - pwm.shape[0]
            motif_pwms[ind] = np.vstack((
                                         np.zeros((int(zero_padding_length/2), 4)),
                                         pwm,
                                         np.zeros((zero_padding_length - int(zero_padding_length/2), 4))
                                ))

        motif_pwms = np.array(motif_pwms)
        return motif_pwms, motif_names

    @staticmethod
    def chrom_loc2motif_scores(fasta_instance, chrom, start, end, strand, motif_pwms, num_threads=1):
        torch.set_num_threads(num_threads)
        sequence = pybedtools.BedTool.seq((chrom,
                                           max(start - int(motif_pwms.shape[1]/2), 0),
                                           end + motif_pwms.shape[1] + 200 - 50 - 1 - int(motif_pwms.shape[1]/2)
                                           ), fasta_instance)
        sequence = sequence.replace("a", "A")
        sequence = sequence.replace("c", "C")
        sequence = sequence.replace("g", "G")
        sequence = sequence.replace("t", "T")
        sequence = sequence.replace("N", "A")
        sequence = sequence.replace("n", "A")
        if strand == "+":
            sequence = sequence.replace("A", "0")
            sequence = sequence.replace("C", "1")
            sequence = sequence.replace("G", "2")
            sequence = sequence.replace("T", "3")
        else:
            motif_pwms = np.flip(motif_pwms, axis=1).copy()
            sequence = sequence.replace("A", "3")
            sequence = sequence.replace("C", "2")
            sequence = sequence.replace("G", "1")
            sequence = sequence.replace("T", "0")

        sequence = np.asarray([0] * max(int(motif_pwms.shape[1]/2) - start, 0) + list(sequence), dtype='int')
        sequence_one_hot = torch.nn.functional.one_hot(torch.from_numpy(sequence), num_classes=4)

        sequence_one_hot = sequence_one_hot.reshape((1, 1, -1, 4)).float()
        motif_filters = torch.from_numpy(motif_pwms.reshape((motif_pwms.shape[0], 1, motif_pwms.shape[1], 4))).float()
        scores = torch.nn.functional.conv2d(sequence_one_hot, motif_filters, padding=0, stride=1)

        #     scores=scores.reshape((motif_pwms.shape[0], -1)).numpy()
        scores = scores.reshape((motif_pwms.shape[0], -1))

        return scores

    @staticmethod
    def torch_topk(tensor, k, dim, largest=None, sorted=None):
        if largest is None:
            largest = True
        if sorted is None:
            sorted = True
        return torch.topk(tensor, k, dim, largest, sorted)[0]

    def prepare_motif_top4_feature(self, chrom, num_threads=1):
        torch.set_num_threads(num_threads)

        batch = self.batch

        chrom_length = self.dic_chrom_length[chrom]
        motif_pwms, motif_names = self.read_motif_pwm()
        fasta_instance = pybedtools.example_filename(self.genome_sequence_fa)

        time_start = time.time()
        # print("Start %s" % chrom)
        with h5py.File("%s/%s_motifs_top4_scores.h5" % (self.motif_feature_path, chrom), "w") as outfile:
            outfile.create_dataset("motif_names", data=np.array(motif_names, dtype='S'), shape=(len(motif_names),),
                                   dtype='S200', compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            h5_scores = outfile.create_dataset("scores", dtype=np.float32,
                                               shape=(len(motif_names), int((chrom_length - 200) / 50) + 1 , 4),
                                               maxshape=(len(motif_names), None, 4),
                                               compression='gzip', shuffle=True, fletcher32=True, compression_opts=4)
            time_scan = time.time()
            for i in range(0, chrom_length, batch):
                time_one = time.time()
                # forward strand
                motif_scores = self.chrom_loc2motif_scores(fasta_instance, chrom, i,
                                                           min((i + batch + motif_pwms.shape[1]),
                                                               chrom_length - motif_pwms.shape[1] - 200 + 50 + 1), '+',
                                                           motif_pwms, num_threads)
                motif_scores_unfold = motif_scores.unfold(dimension=1, size=200, step=50)
                motif_scores_unfold.share_memory_()
                pool = Pool(processes=num_threads)
                torch_topk_new = partial(self.torch_topk, k=3, dim=2, largest=True, sorted=True)
                results = pool.map_async(torch_topk_new, torch.split(motif_scores_unfold, split_size_or_sections=(
                        int(motif_scores_unfold.shape[1] / num_threads) + 1), dim=1))
                motif_scores_top_3_forward = torch.cat(results.get(), dim=1)
                pool.close()
                pool.join()

                del motif_scores, motif_scores_unfold
                gc.collect()

                # reverse strand
                motif_scores = self.chrom_loc2motif_scores(fasta_instance, chrom, i,
                                                           min((i + batch + motif_pwms.shape[1]),
                                                               chrom_length - motif_pwms.shape[1] - 200 + 50 + 1), '-',
                                                           motif_pwms, num_threads)
                motif_scores_unfold = motif_scores.unfold(dimension=1, size=200, step=50)
                motif_scores_unfold.share_memory_()
                pool = Pool(processes=num_threads)
                torch_topk_new = partial(self.torch_topk, k=3, dim=2, largest=True, sorted=True)
                results = pool.map_async(torch_topk_new, torch.split(motif_scores_unfold, split_size_or_sections=(
                        int(motif_scores_unfold.shape[1] / num_threads) + 1), dim=1))
                motif_scores_top_3_reverse = torch.cat(results.get(), dim=1)
                pool.close()
                pool.join()

                del motif_scores, motif_scores_unfold
                gc.collect()

                # top 4
                motif_scores_top_6 = torch.cat((motif_scores_top_3_forward, motif_scores_top_3_reverse), dim=2)
                motif_scores_top_6.share_memory_()
                pool = Pool(processes=num_threads)
                torch_topk_new = partial(self.torch_topk, k=4, dim=2, largest=True, sorted=True)
                results = pool.map_async(torch_topk_new, torch.split(motif_scores_top_6, split_size_or_sections=(
                        int(motif_scores_top_6.shape[1] / num_threads) + 1), dim=1))
                motif_scores_top_4 = torch.cat(results.get(), dim=1)

                pool.close()
                pool.join()

                time_end = time.time()

                print("%.2f percents of chrom have been scanned...(cost %ds)" % (
                    float(min((i + batch + motif_pwms.shape[1]), chrom_length)) / chrom_length * 100,
                    (time_end - time_one)))

                #         h5_scores[:,i:min(i+batch,chrom_length-tf_length_max+1)]=scores
                #             print i,i+motif_scores_top_4.shape[1]
                #             print torch.max(motif_scores_top_4[0,:,0])
                h5_scores[:, int(i / 50):(int(i / 50) + motif_scores_top_4.shape[1]), :] = motif_scores_top_4.numpy()

                outfile.flush()

                time_end = time.time()
                print("%.2f percents of chrom have been scanned and written...(%.2fmins to finish, cost %.2fmins)" %
                      (float(min((i + batch + motif_pwms.shape[1]), chrom_length)) / chrom_length * 100,
                       ((time_end - time_scan) * chrom_length / float(
                           min((i + batch + motif_pwms.shape[1]), chrom_length)) - (time_end - time_scan)) / 60,
                       float(time_end - time_one) / 60))

        print("Done! %s" % chrom)

        time_end = time.time()

        print("total time:%.2fmins" % ((time_end - time_start) / 60.0))


if __name__ == '__main__':
    # fire.Fire(MotifFeatures)
    fire.Fire()

