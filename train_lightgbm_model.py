#!/usr/bin/env python3
import gc
import os
import pickle

import fire
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
from hyperopt.fmin import generate_trials_to_calculate
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import precision_recall_curve
from numpy import linalg as LA
import sklearn.metrics as metrics
import json
import lightgbm as lgb

import numpy as np

import pandas as pd
import glob
from sklearn.preprocessing import QuantileTransformer
import yaml

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from sklearn.metrics import average_precision_score

from early_stopping_avg import early_stopping
from hyperopt import STATUS_OK
from hyperopt import hp
from timeit import default_timer as timer
import numpy as np
from hyperopt import tpe
from hyperopt import Trials
from hyperopt import fmin


def focal_isoform_binary_object(pred, dtrain, alpha=0.5, beta=0.0, gamma=2.0):
    #     alpha controls weight of positives
    #     (0,1) less
    #     >1 more or(0-0.5 less, 0.5-1 more)
    #     beta controls the shift of loss function
    #     >0 to left(less weight to well-trained samples)
    #     gamma controls the steepness of loss function
    #     >0
    label = dtrain.get_label()
    x = beta + (2.0 * label - 1) * gamma * pred
    p = 1. / (1. + np.exp(-x))
    # grad = (1 + (alpha - 1) * label) * (2 * label - 1) * (p - 1)
    grad = (1 - label + (label * 2 - 1) * alpha) * (2 * label - 1) * (p - 1)
    # hess = (1 + (alpha - 1) * label) * gamma * (1 - p) * p
    hess = (1 - label + (label * 2 - 1) * alpha) * gamma * (1 - p) * p
    return grad, hess


def lgb_auprc_score(y_hat, data):
    y_true = data.get_label()
    # TODO try not to round yhat
    # y_hat = np.round(y_hat)  # scikits f1 doesn't like probabilities
    return 'auprc', average_precision_score(y_true, y_hat), True


class LightGBMModel(object):
    def __init__(self, config_file, training_tf_name=None,
                 cofactor_motif_set_file=None, quantile_transformer_path=None, dnase_feature_path=None,
                 motif_feature_path=None, selected_motif_feature_path=None, step=120):
        with open(config_file, "r") as infile:
            config = yaml.load(infile, Loader=Loader)
        self.config = config
        self.chrom_all = config['chrom_all']
        # self.region_topic_model_h5 = config['region_topic_model_h5']
        self.dic_chrom_length = {}
        self.chrom_sets = config['chrom_sets']
        self.training_tf_name = training_tf_name
        self.dic_chrom_length = {}
        with open(config['chrom_size_file'], "r") as infile:
            for line in infile:
                line = line.strip().split("\t")
                if line[0] in self.chrom_all:
                    self.dic_chrom_length[line[0]] = int(line[1])

        # if regions_all_file is not None:
        #     self.df_all_regions = pd.read_csv(regions_all_file, sep="\t", header=None)
        #     self.df_all_regions.columns = ['chr', 'start', 'stop']
        # else:
        #     self.df_all_regions = None

        if training_tf_name is not None:
            print('ok')
            self.df_all_regions_label = pd.read_csv('/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human6/label/train/FOXA1.train.labels.tsv', sep="\t", header=0)
            #self.df_all_regions_label = pd.read_csv(
            #    "%s/%s.%s" % (
            #        config['training_cell_types_regions_label_path'], training_tf_name,
            #        config['training_cell_types_regions_label_name']),
            #    sep="\t", header=0)
        else:
            self.df_all_regions_label = None

        if cofactor_motif_set_file is not None:
            print('ok')
            with open(cofactor_motif_set_file, "r") as infile:
                self.cofactor_motif_set = json.load(infile)
        else:
            self.cofactor_motif_set = None

        if quantile_transformer_path is None:
            self.quantile_transformer_path = "./train/quantile_transformer"
        else:
            self.quantile_transformer_path = quantile_transformer_path

        if dnase_feature_path is None:
            self.dnase_feature_path = "./hdf5s/DNase"
        else:
            self.dnase_feature_path = dnase_feature_path

        if motif_feature_path is None:
            self.motif_feature_path = "./hdf5s/motif"
        else:
            self.motif_feature_path = motif_feature_path

        if selected_motif_feature_path is None:
            self.selected_motif_feature_path = "./hdf5s/motif"
        else:
            self.selected_motif_feature_path = selected_motif_feature_path

        self.step = step

    def prepare_motif_h5_data(self, chrom):

        df_temp = self.df_all_regions_label[self.df_all_regions_label['chr'] == chrom].copy()
        df_temp = df_temp.iloc[:, :3]
        print('check1')
        with h5py.File("%s/%s_motifs_top4_scores.h5" % (self.motif_feature_path, chrom), "r") as infile:
            motif_names = infile['motif_names'][...]
            motif_names = list(map(lambda x: x.decode('UTF-8'), motif_names))
            row_index = [i for i, v in enumerate(motif_names) if v in self.cofactor_motif_set]
            selected_motifs = [motif_names[i] for i in row_index]
            scores = infile["scores"][row_index, :, :]
            print(scores.shape)
            for i in [-7, -5, -3, -1, 0, 1, 3, 5, 7]:
                print("%s %d" % (chrom, i))
                region_index = np.array(list(map(lambda x: x / 50 + i, df_temp["start"])))
                region_index = np.clip(region_index, 0, scores.shape[1] - 1)
                scores_region = scores[:, region_index.astype(int), :]
                for ind, j in enumerate(selected_motifs):
                    for k in range(4):
                        df_temp["%s_offset_%d_top_%d" % (j, i, k)] = scores_region[ind, :, k]
        print(np.array(df_temp).shape)
        print('%s/%s_motif_features_lightGBM.h5' % (self.selected_motif_feature_path, chrom))
        print(1)
        with h5py.File('%s%s_motif_features_lightGBM.h5' % (self.selected_motif_feature_path, chrom), "w") as outfile:
            outfile.create_dataset("feature_names", data=np.array(df_temp.iloc[:, 3:].columns, dtype='S'),
                                   shape=(df_temp.shape[1] - 3,),
                                   dtype='S200', compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            print(2)
            outfile.create_dataset("starts", data=df_temp['start'].tolist(), shape=(df_temp.shape[0],),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            print(3)
            outfile.create_dataset("scores", data=df_temp.iloc[:, 3:].values, dtype=np.float32,
                                   shape=(df_temp.shape[0], df_temp.shape[1] - 3),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=4)

    def prepare_dnase_autoencoder_h5_data(self, cell_line, chrom, outfile_path):

        df_temp = self.df_all_regions_label[self.df_all_regions_label['chr'] == chrom].copy()
        df_temp = df_temp.iloc[:, :3]
        # for cell_line in self.config['cell_types']:
        with h5py.File(
                "%s/DNASE.%s.merge.binSize.1.corrected_sorted_hg19_25bpbin_bwaverage_transformed_%s_scanned_with_autoencoder_v4.hdf5" % (
                        '/n/scratchlfs/xiaoleliu_lab/cchen/Cistrome_imputation/encode/data/DNase_scanning/scan_result',
                        cell_line, chrom), "r") as infile:
            # print(row_index)
            scores = infile["DNase_feature_scanning"][:, :]
            #         for i in [-13,-11,-9,-7,-5,-3,-1,0,1,3,5,7,9,11,13]:
            for i in [-12, -8, -4, 0, 4, 8, 12]:
                # print("%s %d" % (chrom, i))
                region_index = np.array(list(map(lambda x: x / 50 + i, df_temp["start"])))
                region_index = np.clip(region_index, 0, scores.shape[0] - 1)
                scores_region = scores[region_index.astype(int), :]
                for k in range(32):
                    df_temp["DNase_autoencoder_offset_%d_%d" % (i, k)] = scores_region[:, k]
        print('%s/DNASE_autoencoder_lightGBM.%s.%s.h5' % (outfile_path, chrom, cell_line))
        with h5py.File('%s/DNASE_autoencoder_lightGBM.%s.%s.h5' % (outfile_path, chrom, cell_line), "w") as outfile:
            outfile.create_dataset("feature_names", data=np.array(df_temp.iloc[:, 3:].columns, dtype='S'),
                                   shape=(df_temp.shape[1] - 3,),
                                   dtype='S200', compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            outfile.create_dataset("starts", data=df_temp['start'].tolist(), shape=(df_temp.shape[0],),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            outfile.create_dataset("scores", data=df_temp.iloc[:, 3:].values, dtype=np.float32,
                                   shape=(df_temp.shape[0], df_temp.shape[1] - 3),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=4)

    def get_dnase_features(self, cell_line, chrom, dir_dnase_feature_median, selected_bin_index_file=None):
        with h5py.File("%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_all_cell_types.h5" % (
                self.dnase_feature_path, chrom), "r") as infile:
            samples = list(infile['samples'][...])
            cell_line = str(cell_line)
            samples = list(map(lambda x: x.decode('UTF-8'), samples))
            cell_line_index = np.where(np.array(samples) == cell_line)[0][0]
            if selected_bin_index_file is None:
                cell_line_scores = infile[chrom][cell_line_index, :, :]
            else:
                selected_bin_index = np.load(selected_bin_index_file)
                cell_line_scores = infile[chrom][cell_line_index, selected_bin_index, :]
        # with h5py.File("%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_median.h5" % (
        #         dir_dnase_feature_median, chrom), "r") as infile:
        #     if selected_bin_index_file is None:
        #         median_scores = infile[chrom][:, :]
        #     else:
        #         selected_bin_index = np.load(selected_bin_index_file)
        #         median_scores = infile[chrom][selected_bin_index, :]
        # scores = np.hstack((cell_line_scores, median_scores))
        # return scores
        return cell_line_scores

    def get_dnase_features_autoencoder(self, cell_line, chrom, feature_path, selected_bin_index_file=None):
        with h5py.File("%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_all_cell_types.h5" % (
                self.dnase_feature_path, chrom), "r") as infile:
            size = infile[chrom].shape[1]
        with h5py.File("%s/DNASE_autoencoder_lightGBM.%s.%s.h5" % (feature_path, chrom, cell_line), "r") as infile:
            if selected_bin_index_file is None:
                cell_line_scores = infile['scores'][:size, :]
            else:
                selected_bin_index = np.load(selected_bin_index_file)
                cell_line_scores = infile['scores'][:size, :]
        # with h5py.File("%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_median.h5" % (
        #         dir_dnase_feature_median, chrom), "r") as infile:
        #     if selected_bin_index_file is None:
        #         median_scores = infile[chrom][:, :]
        #     else:
        #         selected_bin_index = np.load(selected_bin_index_file)
        #         median_scores = infile[chrom][selected_bin_index, :]
        # scores = np.hstack((cell_line_scores, median_scores))
        # return scores
        return cell_line_scores

    def prepare_lightgbm_binary_dnase_feature(self, cell_line, chrom_set_name, dir_dnase_feature_median, dir_out,
                                              reference=None, selected_bin_index_file=None, ATAC_long_short=False):
        cell_line = str(cell_line)
        chrom = "chr19"
        # TODO change to 50bp or 100bp
        with h5py.File("%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_all_cell_types.h5" % (
                self.dnase_feature_path, chrom), "r") as infile:
            needed_feature_names = list(infile['feature_names'][...])
            needed_feature_names = list(map(lambda x: x.decode('UTF-8'), needed_feature_names))
        # with h5py.File("%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_median.h5" % (
        #         #         dir_dnase_feature_median, chrom), "r") as infile:
        #         #     median_feature_names = list(infile['feature_names'][...])
        #         #     median_feature_names = list(map(lambda x: x.decode('UTF-8'), median_feature_names))
        #         # needed_feature_names += median_feature_names
        list_scores = []
        chrom_set = self.chrom_sets[chrom_set_name]
        if not ATAC_long_short:
            for chrom in chrom_set:
                scores = self.get_dnase_features(cell_line, chrom, dir_dnase_feature_median, selected_bin_index_file)
                list_scores.append(scores)
                print(np.array(scores).shape)
                # print(cell_line, chrom, subset_index)
            all_score = np.vstack(list_scores)
            # TODO change to 50bp or 100bp
            with open("%s/%s_variable_bp_quantile_map.pkl" % (self.quantile_transformer_path, cell_line),
                      'rb') as fin:
                qt = pickle.load(fin, encoding='latin1')
            _ = qt.transform(all_score)
            # _ = qt.transform(all_score[:, :int(all_score.shape[1] / 2)])
            # _ = qt.transform(all_score[:, :cell_line_scores.shape[1]])
            if reference is not None:
                # reference = lgb.Dataset(glob.glob("%s/lightGBM.dnase.*.*.bin" % reference)[0])
                reference = lgb.Dataset(reference)
            print('all_score.shape')
            print(all_score.shape)
            train_data = lgb.Dataset(all_score, feature_name=list(needed_feature_names), reference=reference)
            train_data.save_binary('a.bin', max_bin=255)
        else:
            list_scores_short_long = []
            for frag_size in ['short','long']:
                list_scores = []
                for chrom in chrom_set:
                    scores = self.get_dnase_features("%s_%s" % (cell_line, frag_size), chrom, dir_dnase_feature_median,
                                                     selected_bin_index_file)
                    list_scores.append(scores)
                    # print(cell_line, chrom, subset_index)
                all_score = np.vstack(list_scores)
                # TODO change to 50bp or 100bp
                with open("%s/%s_variable_bp_quantile_map.pkl" % (self.quantile_transformer_path, "%s_%s" % (cell_line, frag_size)),
                          'rb') as fin:
                    qt = pickle.load(fin, encoding='latin1')
                _ = qt.transform(all_score)
                # _ = qt.transform(all_score[:, :int(all_score.shape[1] / 2)])
                # _ = qt.transform(all_score[:, :cell_line_scores.shape[1]])
                list_scores_short_long.append(all_score)
            all_score_short_long = np.hstack(list_scores_short_long)
            if reference is not None:
                # reference = lgb.Dataset(glob.glob("%s/lightGBM.dnase.*.*.bin" % reference)[0])
                reference = lgb.Dataset(reference)
            needed_feature_names_short_long = ['%s_%s' % (feature_name, frag_size) for frag_size in ['short','long']
                                               for feature_name in needed_feature_names
                                               ]
            train_data = lgb.Dataset(all_score_short_long, feature_name=list(needed_feature_names_short_long), reference=reference)
        #train_data.save_binary("%s/lightGBM.dnase.%s.%s.bin" % (dir_out, cell_line, chrom_set_name))

    def prepare_lightgbm_binary_dnase_feature_autoencoder(self, cell_line, chrom_set_name, feature_path,
                                                          dir_out,
                                                          reference=None, selected_bin_index_file=None):
        list_scores = []
        chrom_set = self.chrom_sets[chrom_set_name]
        for chrom in chrom_set:
            scores = self.get_dnase_features_autoencoder(cell_line, chrom, feature_path,
                                                         selected_bin_index_file)
            list_scores.append(scores)
            # print(cell_line, chrom, subset_index)
        all_score = np.vstack(list_scores)
        if reference is not None:
            reference = lgb.Dataset(glob.glob("%s/lightGBM.autoencoder.dnase.*.*.bin" % reference)[0])
        needed_feature_names = []
        for i in [-12, -8, -4, 0, 4, 8, 12]:
            for k in range(32):
                needed_feature_names.append("DNase_autoencoder_offset_%d_%d" % (i, k))
        train_data = lgb.Dataset(all_score, feature_name=list(needed_feature_names), reference=reference)
        train_data.save_binary("%s/lightGBM.autoencoder.dnase.%s.%s.bin" % (dir_out, cell_line, chrom_set_name))

    def prepare_lightgbm_binary_data_motif_feature_subset(self, chrom_set_name, subset_index, dir_out,
                                                          selected_bin_index_file=None, reference=None):
        chrom = "chr19"
        with h5py.File("%s/%s_motif_features_lightGBM.h5" % (self.selected_motif_feature_path, chrom),
                       "r") as infile:
            all_feature_names = list(infile['feature_names'][...])
            all_feature_names = list(map(lambda x: x.decode('UTF-8'), all_feature_names))
            needed_feature_names = [all_feature_names[i:i + self.step]
                                    for i in range(0, len(all_feature_names), self.step)][subset_index - 1]
            feature_index = [list(range(i, min(i + self.step, len(all_feature_names)))) for i in
                             range(0, len(all_feature_names), self.step)][subset_index - 1]
            # needed_feature_names = list(map(lambda x: x.decode('UTF-8'), needed_feature_names))
        list_scores = []
        chrom_set = self.chrom_sets[chrom_set_name]
        # with h5py.File(self.region_topic_model_h5, "r") as region_topic_infile:
        for chrom in chrom_set:
            print(chrom)
            with h5py.File("%s/%s_motif_features_lightGBM.h5" % (self.selected_motif_feature_path, chrom),
                           "r") as infile:
                # feature_names = list(infile['feature_names'][...])
                # feature_names = list(map(lambda x: x.decode('UTF-8'), feature_names))
                # feature_index = [i for i, v in enumerate(feature_names) if (v in needed_feature_names)]
                if selected_bin_index_file is None:
                    scores = infile['scores'][:, feature_index]
                    # if subset_index == 1:
                    #     scores = np.hstack([region_topic_infile[chrom][:, :], scores])
                else:
                    selected_bin_index = np.load(selected_bin_index_file)
                    scores = infile['scores'][selected_bin_index, feature_index]
                    # if subset_index == 1:
                    #     scores = np.hstack([region_topic_infile[chrom][selected_bin_index, :], scores])

            list_scores.append(scores)
            # print(cell_line, chrom, subset_index)
        all_score = np.vstack(list_scores)
        print(all_score.shape)
        if reference is not None:
            reference = lgb.Dataset(glob.glob("%s/lightGBM.motif.*.%d.bin" % (reference, subset_index))[0])
        # if subset_index == 1:
        #     # needed_feature_names = ["topic_%d" % topic_id for topic_id in range(9)] \
        #     #                        + needed_feature_names
        #     # train_data = lgb.Dataset(all_score, categorical_feature=[8],
        #     #                          feature_name=list(needed_feature_names), reference=reference)
        #     needed_feature_names = ["topic_%d" % topic_id for topic_id in range(1)] \
        #                            + needed_feature_names
        #     train_data = lgb.Dataset(all_score[:, 8:],
        #                              categorical_feature=[0],
        #                              feature_name=list(needed_feature_names), reference=reference)
        # else:
        #     train_data = lgb.Dataset(all_score, feature_name=list(needed_feature_names), reference=reference)
        train_data = lgb.Dataset(all_score, feature_name=list(needed_feature_names), reference=reference)
        train_data.save_binary("%s/lightGBM.motif.%s.%d.bin" % (dir_out, chrom_set_name, subset_index))

    # def merge_lightgbm_binary_data(self, cell_line, chrom_set_name, dir_out):
    #     all_feature_names = []
    #     chrom = "chr22"
    #     # TODO change to 50bp or 100bp
    #     # with h5py.File("%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_all_cell_types.h5" % (
    #     #         self.dnase_feature_path, chrom), "r") as infile:
    #     #     all_feature_names += list(infile['feature_names'][...])
    #     # chrom = "chr22"
    #     with h5py.File("%s/%s_motif_features_lightGBM.h5" % (self.selected_motif_feature_path, chrom),
    #                    "r") as infile:
    #         all_feature_names += list(infile['feature_names'][...])
    #     all_feature_names = list(map(lambda x: x.decode('UTF-8'), all_feature_names))
    #     # for cell_line in self.df_all_regions_label.columns.tolist()[3:]:
    #     for cell_line in [cell_line]:
    #         train_data_all = None
    #         for subset_index in range(int(np.ceil(len(all_feature_names) / self.step) + 1)):
    #             train_data = lgb.Dataset("%s/lightGBM.%s.%s.%d.bin" %
    #                                      (dir_out, cell_line, chrom_set_name, subset_index)).construct()
    #             if train_data_all is None:
    #                 train_data_all = train_data
    #             else:
    #                 # train_data_all=train_data_all.add_features_from(train_data)
    #                 train_data_all.add_features_from(train_data)
    #             # print(subset_index)
    #         train_data_all.save_binary("%s/lightGBM_all.%s.%s.bin" % (dir_out, cell_line, chrom_set_name))
    #     print(cell_line, chrom_set_name)

    def merge_lightgbm_binary_data(self, cell_line, chrom_set_name, dir_out=None, lightgbm_dnase_binary_files_path=None,
                                   lightgbm_motif_binary_files_path=None):
        if dir_out is None:
            dir_out = "./train/%s/binary_files" % self.training_tf_name
        if lightgbm_motif_binary_files_path is None:
            lightgbm_motif_binary_files_path = "./train/%s/binary_files" % self.training_tf_name
        if lightgbm_dnase_binary_files_path is None:
            lightgbm_dnase_binary_files_path = "./train/data/dnase_feature_binary_files"
        cell_line = str(cell_line)
        all_feature_names = []
        chrom = "chr19"
        # TODO change to 50bp or 100bp
        with h5py.File("%s/%s_motif_features_lightGBM.h5" % (self.selected_motif_feature_path, chrom),
                       "r") as infile:
            all_feature_names += list(infile['feature_names'][...])
        all_feature_names = list(map(lambda x: x.decode('UTF-8'), all_feature_names))
        train_data_all = lgb.Dataset("%s/lightGBM.dnase.%s.%s.bin" %
                                     (lightgbm_dnase_binary_files_path, cell_line, chrom_set_name)).construct()
        for subset_index in range(int(np.ceil(len(all_feature_names) / self.step))):
            train_data = lgb.Dataset("%s/lightGBM.motif.%s.%d.bin" %
                                     (lightgbm_motif_binary_files_path, chrom_set_name, subset_index + 1)).construct()
            train_data_all.add_features_from(train_data)
        temp = []
        chrom_set = self.chrom_sets[chrom_set_name]
        for chrom in chrom_set:
            df_temp = self.df_all_regions_label.loc[self.df_all_regions_label['chr'] == chrom, :]
            temp.append(df_temp)
        df_all_temp = pd.concat(temp, ignore_index=True)
        # selected_index = np.where(df_all_temp[cell_line] != "A")[0]
        # ignore_index = np.where(df_all_temp[cell_line] == "A")[0]
        # label_b_u = np.delete(np.array(df_all_temp[cell_line]), ignore_index, axis=0)
        # labels = list(map(lambda x: 1 if x == "B" else 0, label_b_u))
        # train_data_all_subset = train_data_all.subset(selected_index)
        # train_data_all_subset.set_label(labels)
        # return train_data_all_subset
        weight = (np.array(df_all_temp[cell_line]) != "A").astype(int)
        train_data_all.set_weight(weight)
        labels = (np.array(df_all_temp[cell_line]) == "B").astype(int)
        train_data_all.set_label(labels)
        # return train_data_all
        train_data_all.save_binary("%s/lightGBM.all.%s.%s.bin" % (dir_out, cell_line, chrom_set_name))

    def merge_lightgbm_binary_data_autoencoder(self, cell_line, chrom_set_name, dir_out=None,
                                               lightgbm_dnase_binary_files_path=None,
                                               lightgbm_motif_binary_files_path=None):
        if dir_out is None:
            dir_out = "./train/%s/binary_files" % self.training_tf_name
        if lightgbm_motif_binary_files_path is None:
            lightgbm_motif_binary_files_path = "./train/%s/binary_files" % self.training_tf_name
        if lightgbm_dnase_binary_files_path is None:
            lightgbm_dnase_binary_files_path = "./train/data/dnase_feature_binary_files"
        all_feature_names = []
        chrom = "chr19"
        # TODO change to 50bp or 100bp
        with h5py.File("%s/%s_motif_features_lightGBM.h5" % (self.selected_motif_feature_path, chrom),
                       "r") as infile:
            all_feature_names += list(infile['feature_names'][...])
        all_feature_names = list(map(lambda x: x.decode('UTF-8'), all_feature_names))
        train_data_all = lgb.Dataset("%s/lightGBM.autoencoder.dnase.%s.%s.bin" %
                                     (lightgbm_dnase_binary_files_path, cell_line, chrom_set_name)).construct()
        # train_data = lgb.Dataset("%s/lightGBM.dnase.%s.%s.bin" %
        #                          (lightgbm_dnase_binary_files_path, cell_line, chrom_set_name)).construct()
        # train_data_all.add_features_from(train_data)
        for subset_index in range(int(np.ceil(len(all_feature_names) / self.step))):
            train_data = lgb.Dataset("%s/lightGBM.motif.%s.%d.bin" %
                                     (lightgbm_motif_binary_files_path, chrom_set_name, subset_index + 1)).construct()
            train_data_all.add_features_from(train_data)
        temp = []
        chrom_set = self.chrom_sets[chrom_set_name]
        for chrom in chrom_set:
            df_temp = self.df_all_regions_label.loc[self.df_all_regions_label['chr'] == chrom, :]
            temp.append(df_temp)
        df_all_temp = pd.concat(temp, ignore_index=True)
        # selected_index = np.where(df_all_temp[cell_line] != "A")[0]
        # ignore_index = np.where(df_all_temp[cell_line] == "A")[0]
        # label_b_u = np.delete(np.array(df_all_temp[cell_line]), ignore_index, axis=0)
        # labels = list(map(lambda x: 1 if x == "B" else 0, label_b_u))
        # train_data_all_subset = train_data_all.subset(selected_index)
        # train_data_all_subset.set_label(labels)
        # return train_data_all_subset
        weight = (np.array(df_all_temp[cell_line]) != "A").astype(int)
        train_data_all.set_weight(weight)
        labels = (np.array(df_all_temp[cell_line]) == "B").astype(int)
        train_data_all.set_label(labels)
        # return train_data_all
        train_data_all.save_binary("%s/lightGBM.autoencoder.all.%s.%s.bin" % (dir_out, cell_line, chrom_set_name))

    def train_models(self, cell_line, chrom_set_name, lightgbm_dnase_binary_files_path=None,
                     lightgbm_motif_binary_files_path=None, dir_out=None, num_threads=16):
        if lightgbm_motif_binary_files_path is None:
            lightgbm_motif_binary_files_path = "./train/%s/binary_files" % self.training_tf_name
        if lightgbm_dnase_binary_files_path is None:
            lightgbm_dnase_binary_files_path = "./train/data/dnase_feature_binary_files"
        if dir_out is None:
            dir_out = "./train/%s/models/" % self.training_tf_name
        cell_line = str(cell_line)
        params = {
            'boosting_type': 'gbdt',
            # 'boosting_type': 'dart',
            # 'drop_rate': 0.3,
            # 'max_drop': 50,
            # 'skip_drop': 0.5,
            # 'drop_seed': 6,
            # 'pos_bagging_fraction': 1,
            # 'neg_bagging_fraction': 0.01,
            # 'bagging_freq': 10000,
            # 'bagging_seed': 6,
            'objective': 'binary',
            # 'objective': focal_binary_object,
            # 'metric': ['binary_error', 'binary_logloss', "auc"],
            'metric': ["auc"],
            # 'is_unbalance': True,
            # "scale_pos_weight": 100,
            'metric_freq': 10,
            'num_leaves': 63,
            # 'num_leaves': 7,
            # 'max_bin': 255,
            'num_threads': num_threads,
            'learning_rate': 0.05,
            'feature_fraction': 1,
            'boost_from_average': False,
            'verbose': 1
        }
        # other_set_type = list(set(self.chrom_sets.keys()) - {chrom_set_name})[0]
        if len(self.df_all_regions_label.columns[3:]) > 1:
            other_cell_lines = list(set(self.df_all_regions_label.columns[3:]) - {cell_line})
        else:
            other_cell_lines = [self.df_all_regions_label.columns[3]]
        # train_data = self.merge_lightgbm_binary_data(cell_line, chrom_set_name, lightgbm_dnase_binary_files_path,
        #                                              lightgbm_motif_binary_files_path)
        train_data = lgb.Dataset(
            "%s/lightGBM.all.%s.%s.bin" % (lightgbm_motif_binary_files_path, cell_line, chrom_set_name))
        list_validation_data = []
        for other_cell_line in other_cell_lines:
            # validation_data = self.merge_lightgbm_binary_data(other_cell_line, "chrom_set_test",
            #                                                   lightgbm_dnase_binary_files_path,
            #                                                   lightgbm_motif_binary_files_path)
            validation_data = lgb.Dataset("%s/lightGBM.all.%s.%s.bin" % (
                lightgbm_motif_binary_files_path, other_cell_line, "chrom_set_test"), reference=train_data)
            list_validation_data.append(validation_data)
        evals_result = {}
        train_data = train_data.construct()
        # see: https://arxiv.org/pdf/1909.04868.pdf
        beta = - np.log10(2 * train_data.num_data()/np.where(train_data.get_label() > 0)[0].shape[0] - 1)
        gbm = lgb.train(params=params,
                        train_set=train_data,
                        fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=0.5, beta=beta, gamma=1),
                        #                 fobj=lambda x,y:logistic_obj(x,y,imbalance_alpha=1.0),
                        valid_sets=[train_data] + list_validation_data,
                        valid_names=['train'] + ["%s_%s" % (other_cell_line, "set_test") \
                                                 for other_cell_line in other_cell_lines],
                        feval=lgb_auprc_score,
                        # early_stopping_rounds=20,
                        evals_result=evals_result,
                        num_boost_round=200,
                        keep_training_booster=False,
                        callbacks=[early_stopping(20, first_metric_only=False, verbose=True)]
                        )
        with open("%s/%s.%s.%s_model.pkl" % (
                dir_out, self.training_tf_name, cell_line, chrom_set_name), 'wb') as fout:
            pickle.dump(gbm, fout)
        with open("%s/%s.%s.%s_evals_result.pkl" % (
                dir_out, self.training_tf_name, cell_line, chrom_set_name),
                  'wb') as outfile_evals_result:
            pickle.dump(evals_result, outfile_evals_result, pickle.HIGHEST_PROTOCOL)

    def train_models_hyperopt(self, cell_line, chrom_set_name, lightgbm_dnase_binary_files_path=None,
                              lightgbm_motif_binary_files_path=None, dir_out=None, num_threads=16):
        if lightgbm_motif_binary_files_path is None:
            lightgbm_motif_binary_files_path = "./train/%s/binary_files" % self.training_tf_name
        if lightgbm_dnase_binary_files_path is None:
            lightgbm_dnase_binary_files_path = "./train/data/dnase_feature_binary_files"
        if dir_out is None:
            dir_out = "./train/%s/models/" % self.training_tf_name

        # other_set_type = list(set(self.chrom_sets.keys()) - {chrom_set_name})[0]
        if len(self.df_all_regions_label.columns[3:]) > 1:
            other_cell_lines = list(set(self.df_all_regions_label.columns[3:]) - {cell_line})
        else:
            other_cell_lines = [self.df_all_regions_label.columns[3]]
        # train_data = self.merge_lightgbm_binary_data(cell_line, chrom_set_name, lightgbm_dnase_binary_files_path,
        #                                              lightgbm_motif_binary_files_path)
        train_data = lgb.Dataset(
            "%s/lightGBM.all.%s.%s.bin" % (lightgbm_motif_binary_files_path, cell_line, chrom_set_name))
        list_validation_data = []
        for other_cell_line in other_cell_lines:
            # validation_data = self.merge_lightgbm_binary_data(other_cell_line, "chrom_set_test",
            #                                                   lightgbm_dnase_binary_files_path,
            #                                                   lightgbm_motif_binary_files_path)
            validation_data = lgb.Dataset("%s/lightGBM.all.%s.%s.bin" % (
                lightgbm_motif_binary_files_path, other_cell_line, "chrom_set_test"), reference=train_data)
            list_validation_data.append(validation_data)

        def hyperopt_objective(argsDict):
            """Objective function for Gradient Boosting Machine Hyperparameter Optimization"""
            # Keep track of evals
            global ITERATION
            ITERATION += 1

            global bayes_trials
            with open("%s/%s.%s.%s_bayes_trials.pkl" % (
                    dir_out, self.training_tf_name, cell_line, chrom_set_name), 'wb') as fout:
                pickle.dump(bayes_trials, fout)
            # Make sure parameters that need to be integers are integers
            for parameter_name in ['num_leaves',
                                   'min_data_in_leaf',
                                   # 'max_depth'
                                   ]:
                argsDict[parameter_name] = int(argsDict[parameter_name])
            start = timer()
            params = {
                'boosting_type': 'gbdt',
                # 'ignore_column': list(range(500)),
                # 'boosting_type': 'dart',
                # 'drop_rate': 0.3,
                # 'max_drop': 50,
                # 'skip_drop': 0.5,
                # 'drop_seed': 6,
                # 'pos_bagging_fraction': 0.001,
                # 'neg_bagging_fraction': 0.001,
                # 'bagging_freq': 10000,
                # 'bagging_seed': 6,
                'objective': 'binary',
                # 'objective': focal_binary_object,
                # 'metric': ['binary_error', 'binary_logloss', "auc"],
                'metric': ["auc"],
                # 'first_metric_only': True,
                # 'is_unbalance': True,
                # "scale_pos_weight": 100,
                # 'feature_fraction_bynode': True,
                'metric_freq': 10,
                'num_leaves': argsDict['num_leaves'],
                'min_data_in_leaf': argsDict['min_data_in_leaf'],
                # 'min_data_in_leaf': 20,
                # 'max_depth': argsDict['max_depth'],
                # 'min_sum_hessian_in_leaf': argsDict['min_sum_hessian_in_leaf'],
                # 'bagging_fraction': argsDict['bagging_fraction'],
                # 'feature_fraction': argsDict['feature_fraction'],
                # 'lambda_l1': argsDict['lambda_l1'],
                # 'lambda_l2': argsDict['lambda_l2'],
                # 'max_bin': 255,
                'num_threads': num_threads,
                # 'learning_rate': argsDict['learning_rate'],
                'learning_rate': 0.1,
                'bagging_freq': 1,
                'boost_from_average': False,
                'verbose': 1
            }
            evals_result = {}
            valid_names = ['train'] + ["%s_%s" % (other_cell_line, "set_test") \
                                       for other_cell_line in other_cell_lines]
            gbm = lgb.train(params,
                            train_set=train_data,
                            fobj=lambda x, y: focal_isoform_binary_object(x, y,
                                                                          # alpha=float(
                                                                          #     np.clip(argsDict['alpha'], 0.001, 0.999)),
                                                                          alpha=1. / (1. + np.exp(
                                                                              -argsDict['alpha_isoform'])),
                                                                          beta=argsDict['beta'],
                                                                          gamma=argsDict['gamma']),
                            valid_sets=[train_data] + list_validation_data,
                            valid_names=valid_names,
                            feval=lgb_auprc_score,
                            num_boost_round=300,
                            # early_stopping_rounds=20,
                            evals_result=evals_result,
                            keep_training_booster=False,
                            callbacks=[early_stopping(20, first_metric_only=False, verbose=True)],
                            )
            run_time = timer() - start
            # Extract the best score
            # best_score = np.max(cv_results['auprc-mean'])
            auprc_sum = None
            n = 0
            for valid_name in valid_names:
                if valid_name != "train":
                    if auprc_sum is None:
                        auprc_sum = np.array(evals_result[valid_name]['auprc'])
                    else:
                        auprc_sum += np.array(evals_result[valid_name]['auprc'])
                    n += 1
            best_score = np.max(auprc_sum / n)

            # Loss must be minimized
            loss = 1 - best_score
            # Boosting rounds that returned the highest cv score
            # n_estimators = int(np.argmax(cv_results['auprc-mean']) + 1)
            n_estimators = int(np.argmax(auprc_sum) + 1)
            print('auprc:{} ITERATION:{} n_estimators:{} run_time:{}'.format(best_score, ITERATION, n_estimators,
                                                                             run_time),
                  end="\n")
            # Dictionary with information for evaluation
            return {'loss': loss,
                    'params': argsDict,
                    'iteration': ITERATION,
                    'estimators': n_estimators,
                    'gbm': gbm,
                    'evals_result': evals_result,
                    'train_time': run_time,
                    'status': STATUS_OK}
            # return loss

        # Define the search space
        space = {
            # 'class_weight': hp.choice('class_weight', [None, 'balanced']),
            'num_leaves': hp.qloguniform('num_leaves', np.log(15), np.log(1023), 5),
            # 'max_depth': hp.quniform('max_depth', 3, 63, 1),
            # 'min_sum_hessian_in_leaf': hp.loguniform('min_sum_hessian_in_leaf', np.log(0.001), np.log(1)),
            # 'learning_rate': hp.loguniform('learning_rate', np.log(0.001), np.log(0.2)),
            'min_data_in_leaf': hp.quniform('min_data_in_leaf', 20, 1000, 5),
            # 'lambda_l1': hp.uniform('lambda_l1', 0.0, 1.0),
            # 'lambda_l2': hp.uniform('lambda_l2', 0.0, 1.0),
            # 'bagging_fraction': hp.uniform('bagging_fraction', 0.4, 1.0),
            # 'feature_fraction': hp.uniform('feature_fraction', 0.4, 1.0),
            # 'alpha': hp.loguniform('alpha', np.log(1), np.log(100)),
            # 'alpha': hp.normal('alpha', 0.5, 0.15),
            'alpha_isoform': hp.normal('alpha_isoform', 0, 3),
            'beta': hp.uniform('beta', -10, 10),
            'gamma': hp.loguniform('gamma', np.log(1), np.log(20)),
        }


        # Keep track of results
        global bayes_trials
        # bayes_trials = Trials()
        bayes_trials_file_path = "%s/%s.%s.%s_bayes_trials.pkl" % (
                dir_out, self.training_tf_name, cell_line, chrom_set_name)
        if os.path.exists(bayes_trials_file_path):
            with open(bayes_trials_file_path, 'rb') as fin:
                bayes_trials = pickle.load(fin, encoding='latin1')
        else:
            bayes_trials = generate_trials_to_calculate(
                [{'num_leaves': 63, 'min_data_in_leaf': 20, 'alpha_isoform': 0, 'beta': -1.5, 'gamma': 1.01}])

        # Global variable
        global ITERATION

        ITERATION = 0

        # Run optimization
        best = fmin(fn=hyperopt_objective, space=space, algo=tpe.suggest,
                    max_evals=len(bayes_trials.tids)+30, trials=bayes_trials, rstate=np.random.RandomState(6))

        # Sort the trials with lowest loss (highest AUC) first

        bayes_trials_results = sorted(bayes_trials.results, key=lambda x: x['loss'])
        # bayes_trials_results[:10]

        with open("%s/%s.%s.%s_model.hyperopt.pkl" % (
                dir_out, self.training_tf_name, cell_line, chrom_set_name), 'wb') as fout:
            pickle.dump(bayes_trials_results[0]['gbm'], fout)
        with open("%s/%s.%s.%s_evals_result.hyperopt.pkl" % (
                dir_out, self.training_tf_name, cell_line, chrom_set_name),
                  'wb') as outfile_evals_result:
            pickle.dump(bayes_trials_results[0]['evals_result'], outfile_evals_result, pickle.HIGHEST_PROTOCOL)
        with open("%s/%s.%s.%s_bayes_trials.pkl" % (
                dir_out, self.training_tf_name, cell_line, chrom_set_name), 'wb') as fout:
            pickle.dump(bayes_trials, fout)

    def train_models_autoencoder(self, cell_line, chrom_set_name, lightgbm_dnase_binary_files_path=None,
                                 lightgbm_motif_binary_files_path=None, dir_out=None, num_threads=16):
        if lightgbm_motif_binary_files_path is None:
            lightgbm_motif_binary_files_path = "./train/%s/binary_files" % self.training_tf_name
        if lightgbm_dnase_binary_files_path is None:
            lightgbm_dnase_binary_files_path = "./train/data/dnase_feature_binary_files"
        if dir_out is None:
            dir_out = "./train/%s/models/" % self.training_tf_name
        params = {
            'boosting_type': 'gbdt',
            # 'boosting_type': 'dart',
            # 'drop_rate': 0.3,
            # 'max_drop': 50,
            # 'skip_drop': 0.5,
            # 'drop_seed': 6,
            # 'pos_bagging_fraction': 1,
            # 'neg_bagging_fraction': 0.01,
            # 'bagging_freq': 10000,
            # 'bagging_seed': 6,
            'objective': 'binary',
            # 'objective': focal_binary_object,
            # 'metric': ['binary_error', 'binary_logloss', "auc"],
            'metric': ["auc"],
            # 'is_unbalance': True,
            # "scale_pos_weight": 100,
            'metric_freq': 10,
            'num_leaves': 63,
            # 'max_bin': 255,
            'num_threads': num_threads,
            'learning_rate': 0.1,
            'feature_fraction': 1,
            'boost_from_average': False,
            'verbose': 1
        }
        # other_set_type = list(set(self.chrom_sets.keys()) - {chrom_set_name})[0]
        if len(self.df_all_regions_label.columns[3:]) > 1:
            other_cell_lines = list(set(self.df_all_regions_label.columns[3:]) - {cell_line})
        else:
            other_cell_lines = [self.df_all_regions_label.columns[3]]
        # train_data = self.merge_lightgbm_binary_data(cell_line, chrom_set_name, lightgbm_dnase_binary_files_path,
        #                                              lightgbm_motif_binary_files_path)
        train_data = lgb.Dataset(
            "%s/lightGBM.autoencoder.all.%s.%s.bin" % (lightgbm_motif_binary_files_path, cell_line, chrom_set_name))
        list_validation_data = []
        for other_cell_line in other_cell_lines:
            # validation_data = self.merge_lightgbm_binary_data(other_cell_line, "chrom_set_test",
            #                                                   lightgbm_dnase_binary_files_path,
            #                                                   lightgbm_motif_binary_files_path)
            validation_data = lgb.Dataset("%s/lightGBM.autoencoder.all.%s.%s.bin" % (
                lightgbm_motif_binary_files_path, other_cell_line, "chrom_set_test"), reference=train_data)
            list_validation_data.append(validation_data)
        evals_result = {}
        gbm = lgb.train(params=params,
                        train_set=train_data,
                        fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=0.5, beta=-1.5, gamma=1.01),
                        #                 fobj=lambda x,y:logistic_obj(x,y,imbalance_alpha=1.0),
                        valid_sets=[train_data] + list_validation_data,
                        valid_names=['train'] + ["%s_%s" % (other_cell_line, "set_test") \
                                                 for other_cell_line in other_cell_lines],
                        feval=lgb_auprc_score,
                        early_stopping_rounds=20,
                        evals_result=evals_result,
                        keep_training_booster=False)
        with open("%s/%s.%s.%s_model.autoencoder.pkl" % (
                dir_out, self.training_tf_name, cell_line, chrom_set_name), 'wb') as fout:
            pickle.dump(gbm, fout)
        with open("%s/%s.%s.%s_evals_result.autoencoder.pkl" % (
                dir_out, self.training_tf_name, cell_line, chrom_set_name),
                  'wb') as outfile_evals_result:
            pickle.dump(evals_result, outfile_evals_result, pickle.HIGHEST_PROTOCOL)

    def make_prediction(self, cell_line, chrom, dir_dnase_feature_median=None, lightgbm_model_files_path=None,
                        dir_out=None):
        if dir_dnase_feature_median is None:
            dir_dnase_feature_median = "./hdf5s/DNase/median"
        if lightgbm_model_files_path is None:
            lightgbm_model_files_path = "./train/%s/models/" % self.training_tf_name
        if dir_out is None:
            dir_out = "./train/%s/predictions/" % self.training_tf_name
        cell_line = str(cell_line)
        scores = self.get_dnase_features(cell_line, chrom, dir_dnase_feature_median)
        with open("%s/%s_variable_bp_quantile_map.pkl" % (self.quantile_transformer_path, cell_line),
                  'rb') as fin:
            qt = pickle.load(fin, encoding='latin1')
        # _ = qt.transform(scores[:, :int(scores.shape[1] / 2)])
        _ = qt.transform(scores)
        # with h5py.File(self.region_topic_model_h5, "r") as region_topic_infile:
        #     scores_topic = region_topic_infile[chrom][...]
        with h5py.File("%s/%s_motif_features_lightGBM.h5" % (self.selected_motif_feature_path, chrom),
                       "r") as infile:
            scores_motif = infile['scores'][...]
        # scores_all = np.hstack((scores, scores_topic[:, 8:], scores_motif))
        scores_all = np.hstack((scores, scores_motif))
        # model_files = glob.glob("%s/%s_*_model.pkl" % (lightgbm_model_files_path, self.training_tf_name))
        model_files = ["%s/%s.%s.%s_model.pkl" % (
            lightgbm_model_files_path, self.training_tf_name, training_cell_line, training_chrom_set_name)
                       for training_cell_line in list(self.df_all_regions_label.columns[3:])
                       for training_chrom_set_name in sorted(list(set(self.chrom_sets.keys()) - {'chrom_set_test'}))]
        model_files = np.array(model_files, dtype='S')
        preds = np.zeros((scores_all.shape[0], len(model_files)))
        for ind_model_file, model_file in enumerate(model_files):
            with open(model_file, 'rb') as fin:
                gbm = pickle.load(fin, encoding='latin1')
            ypred = gbm.predict(scores_all)
            preds[:, ind_model_file] = ypred
        with h5py.File('%s/%s.%s.%s_preds.h5' % (dir_out, self.training_tf_name, cell_line, chrom),
                       "w") as outfile:
            outfile.create_dataset("model_files", data=model_files,
                                   shape=(len(model_files),),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            outfile.create_dataset("preds", data=preds,
                                   shape=preds.shape,
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)

    def make_prediction_hyperopt(self, cell_line, chrom, dir_dnase_feature_median=None, lightgbm_model_files_path=None,
                        dir_out=None):
        if dir_dnase_feature_median is None:
            dir_dnase_feature_median = "./hdf5s/DNase/median"
        if lightgbm_model_files_path is None:
            lightgbm_model_files_path = "./train/%s/models/" % self.training_tf_name
        if dir_out is None:
            dir_out = "./train/%s/predictions/" % self.training_tf_name
        scores = self.get_dnase_features(cell_line, chrom, dir_dnase_feature_median)
        with open("%s/%s_variable_bp_quantile_map.pkl" % (self.quantile_transformer_path, cell_line),
                  'rb') as fin:
            qt = pickle.load(fin, encoding='latin1')
        # _ = qt.transform(scores[:, :int(scores.shape[1] / 2)])
        _ = qt.transform(scores)
        with h5py.File("%s/%s_motif_features_lightGBM.h5" % (self.selected_motif_feature_path, chrom),
                       "r") as infile:
            scores_motif = infile['scores'][...]
        scores_all = np.hstack((scores, scores_motif))
        # model_files = glob.glob("%s/%s_*_model.pkl" % (lightgbm_model_files_path, self.training_tf_name))
        model_files = ["%s/%s.%s.%s_model.hyperopt.pkl" % (
            lightgbm_model_files_path, self.training_tf_name, training_cell_line, training_chrom_set_name)
                       for training_cell_line in list(self.df_all_regions_label.columns[3:])
                       for training_chrom_set_name in sorted(list(set(self.chrom_sets.keys()) - {'chrom_set_test'}))]
        model_files = np.array(model_files, dtype='S')
        preds = np.zeros((scores_all.shape[0], len(model_files)))
        for ind_model_file, model_file in enumerate(model_files):
            with open(model_file, 'rb') as fin:
                gbm = pickle.load(fin, encoding='latin1')
            ypred = gbm.predict(scores_all)
            preds[:, ind_model_file] = ypred
        with h5py.File('%s/%s.%s.%s_preds.hyperopt.h5' % (dir_out, self.training_tf_name, cell_line, chrom),
                       "w") as outfile:
            outfile.create_dataset("model_files", data=model_files,
                                   shape=(len(model_files),),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            outfile.create_dataset("preds", data=preds,
                                   shape=preds.shape,
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)

    def make_prediction_leaf(self, cell_line, chrom, dir_dnase_feature_median=None, lightgbm_model_files_path=None,
                             dir_out=None):
        if dir_dnase_feature_median is None:
            dir_dnase_feature_median = "./hdf5s/DNase/median"
        if lightgbm_model_files_path is None:
            lightgbm_model_files_path = "./train/%s/models/" % self.training_tf_name
        if dir_out is None:
            dir_out = "./train/%s/predictions/" % self.training_tf_name
        scores = self.get_dnase_features(cell_line, chrom, dir_dnase_feature_median)
        with open("%s/%s_variable_bp_quantile_map.pkl" % (self.quantile_transformer_path, cell_line),
                  'rb') as fin:
            qt = pickle.load(fin, encoding='latin1')
        # _ = qt.transform(scores[:, :int(scores.shape[1] / 2)])
        _ = qt.transform(scores)
        with h5py.File("%s/%s_motif_features_lightGBM.h5" % (self.selected_motif_feature_path, chrom),
                       "r") as infile:
            scores_motif = infile['scores'][...]
        scores_all = np.hstack((scores, scores_motif))
        # model_files = glob.glob("%s/%s_*_model.pkl" % (lightgbm_model_files_path, self.training_tf_name))
        model_files = ["%s/%s.%s.%s_model.pkl" % (
            lightgbm_model_files_path, self.training_tf_name, training_cell_line, training_chrom_set_name)
                       for training_cell_line in list(self.df_all_regions_label.columns[3:])
                       for training_chrom_set_name in sorted(list(set(self.chrom_sets.keys()) - {'chrom_set_test'}))]
        with h5py.File('%s/%s.%s.%s_pred_leafs.h5' % (dir_out, self.training_tf_name, cell_line, chrom),
                       "w") as outfile:
            for ind_model_file, model_file in enumerate(model_files):
                with open(model_file, 'rb') as fin:
                    gbm = pickle.load(fin, encoding='latin1')
                leafs = gbm.predict(scores_all, raw_score=False, pred_leaf=True, pred_contrib=False)
                leaf_outputs = np.zeros(leafs.shape)
                for i in range(leafs.shape[0]):
                    for j in range(leafs.shape[1]):
                        leaf_outputs[i, j] = gbm.get_leaf_output(j, leafs[i, j])
                gc.collect()
                outfile.create_dataset("%s/%s" % (model_file.split("/")[-1], cell_line),
                                       data=leaf_outputs,
                                       shape=leaf_outputs.shape,
                                       compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
                outfile.flush()

                training_cell_line = model_file.split("/")[-1].split('.')[1]
                source_scores = self.get_dnase_features(training_cell_line, chrom, dir_dnase_feature_median)
                with open("%s/%s_variable_bp_quantile_map.pkl" % (self.quantile_transformer_path, training_cell_line),
                          'rb') as fin:
                    qt = pickle.load(fin, encoding='latin1')
                # _ = qt.transform(scores[:, :int(scores.shape[1] / 2)])
                _ = qt.transform(source_scores)
                # source_scores_all = np.hstack((source_scores, scores_motif))
                gc.collect()
                leafs = gbm.predict(np.hstack((source_scores, scores_motif)), raw_score=False, pred_leaf=True,
                                    pred_contrib=False)
                leaf_outputs = np.zeros(leafs.shape)
                for i in range(leafs.shape[0]):
                    for j in range(leafs.shape[1]):
                        leaf_outputs[i, j] = gbm.get_leaf_output(j, leafs[i, j])
                outfile.create_dataset("%s/%s" % (model_file.split("/")[-1], training_cell_line),
                                       data=leaf_outputs,
                                       shape=leaf_outputs.shape,
                                       compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
                outfile.flush()
            model_files = np.array(model_files, dtype='S')
            outfile.create_dataset("model_files", data=model_files,
                                   shape=(len(model_files),),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)

    def make_prediction_autoencoder(self, cell_line, chrom, dir_dnase_feature_median=None,
                                    lightgbm_model_files_path=None,
                                    dir_out=None):
        if dir_dnase_feature_median is None:
            dir_dnase_feature_median = "./hdf5s/DNase/median"
        if lightgbm_model_files_path is None:
            lightgbm_model_files_path = "./train/%s/models/" % self.training_tf_name
        if dir_out is None:
            dir_out = "./train/%s/predictions/" % self.training_tf_name
        scores_autoencoder = self.get_dnase_features_autoencoder(cell_line, chrom, './hdf5s/DNase')
        # scores = self.get_dnase_features(cell_line, chrom, dir_dnase_feature_median)
        # with open("%s/%s_variable_bp_quantile_map.pkl" % (self.quantile_transformer_path, cell_line),
        #           'rb') as fin:
        #     qt = pickle.load(fin, encoding='latin1')
        # # _ = qt.transform(scores[:, :int(scores.shape[1] / 2)])
        # _ = qt.transform(scores)
        with h5py.File("%s/%s_motif_features_lightGBM.h5" % (self.selected_motif_feature_path, chrom),
                       "r") as infile:
            scores_motif = infile['scores'][...]
        # scores_all = np.hstack((scores_autoencoder, scores, scores_motif))
        scores_all = np.hstack((scores_autoencoder, scores_motif))
        # model_files = glob.glob("%s/%s_*_model.pkl" % (lightgbm_model_files_path, self.training_tf_name))
        model_files = ["%s/%s.%s.%s_model.autoencoder.pkl" % (
            lightgbm_model_files_path, self.training_tf_name, training_cell_line, training_chrom_set_name)
                       for training_cell_line in list(self.df_all_regions_label.columns[3:])
                       for training_chrom_set_name in sorted(list(set(self.chrom_sets.keys()) - {'chrom_set_test'}))]
        model_files = np.array(model_files, dtype='S')
        preds = np.zeros((scores_all.shape[0], len(model_files)))
        for ind_model_file, model_file in enumerate(model_files):
            with open(model_file, 'rb') as fin:
                gbm = pickle.load(fin, encoding='latin1')
            ypred = gbm.predict(scores_all)
            preds[:, ind_model_file] = ypred
        with h5py.File('%s/%s.%s.%s_preds.autoencoder.h5' % (dir_out, self.training_tf_name, cell_line, chrom),
                       "w") as outfile:
            outfile.create_dataset("model_files", data=model_files,
                                   shape=(len(model_files),),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            outfile.create_dataset("preds", data=preds,
                                   shape=preds.shape,
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)

    def evaluation(self, cell_line, lightgbm_preds_files_path=None, dir_out=None):
        if dir_out is None:
            dir_out = "./train/%s/evaluations/" % self.training_tf_name
        cell_line = str(cell_line)
        df_test_regions_label = pd.read_csv(
            "%s/%s.%s" % (
                self.config['test_cell_types_regions_label_path'], self.training_tf_name,
                self.config['test_cell_types_regions_label_name']), sep="\t", header=0)
        list_preds_binary = []
        # list_preds_binary_2 = []
        list_labels = []
        list_preds_matrix = []
        list_chroms = []
        list_starts = []
        for chrom in self.chrom_all:
            with h5py.File(
                    '%s/%s.%s.%s_preds.h5' % (lightgbm_preds_files_path, self.training_tf_name, cell_line, chrom),
                    "r") as infile:
                model_files = infile['model_files'][...]
                preds = infile['preds'][...]
                labels = np.array(df_test_regions_label.loc[df_test_regions_label['chr'] == chrom, :][cell_line])
                list_preds_matrix.append(preds)
                preds_binary = np.mean(1. / (1. + np.exp(-preds)), axis=1)
                # preds_binary_2 = 1. / (1. + np.exp(-np.mean(preds, axis=1)))
                list_preds_binary.append(preds_binary)
                # list_preds_binary_2.append(preds_binary_2)
                list_labels.append(labels)
            with h5py.File("%s/%s_motif_features_lightGBM.h5" % (self.selected_motif_feature_path, chrom),
                           "r") as infile:
                list_starts.append(infile['starts'][...])
                list_chroms.append(np.array([chrom] * infile['starts'].shape[0]))
        labels = np.hstack(list_labels)
        preds = np.hstack(list_preds_binary)
        # preds_2 = np.hstack(list_preds_binary_2)
        preds_matrix = np.vstack(list_preds_matrix)
        starts = np.hstack(list_starts)
        chroms = np.hstack(list_chroms)
        ignore_index = np.where(labels == "A")[0]
        preds_matrix = np.delete(preds_matrix, ignore_index, axis=0)
        preds = np.delete(preds, ignore_index, axis=0)
        label_b_u = np.delete(labels, ignore_index, axis=0)
        starts = np.delete(starts, ignore_index, axis=0)
        chroms = np.delete(chroms, ignore_index, axis=0)
        label_b_u = np.array(list(map(lambda x: 1 if x == "B" else 0, label_b_u)))
        with open("%s/%s.%s_performance.txt" % (dir_out, self.training_tf_name, cell_line), "w") as outfile:
            fpr, tpr, thresholds = metrics.roc_curve(label_b_u, preds, pos_label=1)
            auc = metrics.auc(fpr, tpr)
            auprc = average_precision_score(label_b_u, preds)
            outfile.write("average model: auc:%.6f auprc:%.6f\n" % (auc, auprc))
            temp = []
            for i in range(preds_matrix.shape[1]):
                fpr, tpr, thresholds = metrics.roc_curve(label_b_u, preds_matrix[:, i], pos_label=1)
                auc = metrics.auc(fpr, tpr)
                # auprc = average_precision_score(label_b_u, preds_matrix[:, i])
                auprc = average_precision_score(label_b_u, 1. / (1. + np.exp(-preds_matrix[:, i])))
                outfile.write("%s model: auc:%.6f auprc:%.6f\n" % (
                    model_files[i].decode().split("/")[-1].replace('_model.pkl', ''), auc, auprc))
                precision, recall, thresholds = precision_recall_curve(label_b_u,
                                                                       1. / (1. + np.exp(-preds_matrix[:, i])),
                                                                       pos_label=1)
                df_temp = pd.DataFrame(None)
                df_temp["precision"] = precision
                df_temp["recall"] = recall
                df_temp["model"] = model_files[i].decode().split("/")[-1]
                temp.append(df_temp.sample(n=min(100000, df_temp.shape[0])))
            df_plot = pd.concat(temp, ignore_index=True)
            plt.figure(figsize=(8, 6))
            ax = sns.lineplot(x="recall", y="precision", data=df_plot, hue='model', palette="tab10")
            ax.set_title("%s in %s" % (self.training_tf_name, cell_line))
            ax.set_xlabel("Recall")
            ax.set_ylabel("Precision")
            ax.get_figure().savefig("%s/%s_%s_PRC.pdf" % (dir_out, self.training_tf_name, cell_line))
            df_plot.to_csv(
                    "%s/df_plot.PRC.%s.xls" % (
                        dir_out, cell_line),
                    sep="\t", header=True, index=False)

        with open("%s/%s.%s_confusion_matrix.txt" % (dir_out, self.training_tf_name, cell_line), "w") as outfile:
            for i in range(preds_matrix.shape[1]):
                one_preds = 1. / (1. + np.exp(-preds_matrix[:, i]))
                cutoff = 0.5
                true_positive = np.where((one_preds >= cutoff) & (label_b_u == 1))[0]
                false_positive = np.where((one_preds >= cutoff) & (label_b_u == 0))[0]
                false_negative = np.where((one_preds < cutoff) & (label_b_u == 1))[0]
                outfile.write("%s model: all_regions:%d true_positive:%d false_positive:%d false_negative:%d\n" % (
                    model_files[i].decode().split("/")[-1].replace('_model.pkl', ''), len(one_preds),
                    len(true_positive), len(false_positive), len(false_negative)))
                df = pd.DataFrame(None)
                df["chrom"] = np.hstack((chroms[true_positive], chroms[false_positive], chroms[false_negative]))
                df["start"] = np.hstack((starts[true_positive], starts[false_positive], starts[false_negative]))
                df["preds"] = np.hstack(
                    (one_preds[true_positive], one_preds[false_positive], one_preds[false_negative]))
                df["label"] = np.hstack(
                    (label_b_u[true_positive], label_b_u[false_positive], label_b_u[false_negative]))
                df["class"] = ['true_positive'] * len(true_positive) + ['false_positive'] * len(false_positive) + [
                    'false_negative'] * len(false_negative)
                df.to_csv(
                    "%s/df.%s.regions.%s.xls" % (
                        dir_out, model_files[i].decode().split("/")[-1].replace('_model.pkl', ''), cell_line),
                    sep="\t", header=True, index=False)

    def evaluation_hyperopt(self, cell_line, lightgbm_preds_files_path=None, dir_out=None):
        if dir_out is None:
            dir_out = "./train/%s/evaluations/" % self.training_tf_name
        df_test_regions_label = pd.read_csv(
            "%s/%s.%s" % (
                self.config['test_cell_types_regions_label_path'], self.training_tf_name,
                self.config['test_cell_types_regions_label_name']), sep="\t", header=0)
        list_preds_binary = []
        # list_preds_binary_2 = []
        list_labels = []
        list_preds_matrix = []
        list_chroms = []
        list_starts = []
        for chrom in self.chrom_all:
            with h5py.File(
                    '%s/%s.%s.%s_preds.hyperopt.h5' % (lightgbm_preds_files_path, self.training_tf_name, cell_line, chrom),
                    "r") as infile:
                model_files = infile['model_files'][...]
                preds = infile['preds'][...]
                labels = np.array(df_test_regions_label.loc[df_test_regions_label['chr'] == chrom, :][cell_line])
                list_preds_matrix.append(preds)
                preds_binary = np.mean(1. / (1. + np.exp(-preds)), axis=1)
                # preds_binary_2 = 1. / (1. + np.exp(-np.mean(preds, axis=1)))
                list_preds_binary.append(preds_binary)
                # list_preds_binary_2.append(preds_binary_2)
                list_labels.append(labels)
            with h5py.File("%s/%s_motif_features_lightGBM.h5" % (self.selected_motif_feature_path, chrom),
                           "r") as infile:
                list_starts.append(infile['starts'][...])
                list_chroms.append(np.array([chrom] * infile['starts'].shape[0]))
        labels = np.hstack(list_labels)
        preds = np.hstack(list_preds_binary)
        # preds_2 = np.hstack(list_preds_binary_2)
        preds_matrix = np.vstack(list_preds_matrix)
        starts = np.hstack(list_starts)
        chroms = np.hstack(list_chroms)
        ignore_index = np.where(labels == "A")[0]
        preds_matrix = np.delete(preds_matrix, ignore_index, axis=0)
        preds = np.delete(preds, ignore_index, axis=0)
        label_b_u = np.delete(labels, ignore_index, axis=0)
        starts = np.delete(starts, ignore_index, axis=0)
        chroms = np.delete(chroms, ignore_index, axis=0)
        label_b_u = np.array(list(map(lambda x: 1 if x == "B" else 0, label_b_u)))
        with open("%s/%s.%s_performance.hyperopt.txt" % (dir_out, self.training_tf_name, cell_line), "w") as outfile:
            fpr, tpr, thresholds = metrics.roc_curve(label_b_u, preds, pos_label=1)
            auc = metrics.auc(fpr, tpr)
            auprc = average_precision_score(label_b_u, preds)
            outfile.write("average model: auc:%.6f auprc:%.6f\n" % (auc, auprc))
            temp = []
            for i in range(preds_matrix.shape[1]):
                fpr, tpr, thresholds = metrics.roc_curve(label_b_u, preds_matrix[:, i], pos_label=1)
                auc = metrics.auc(fpr, tpr)
                # auprc = average_precision_score(label_b_u, preds_matrix[:, i])
                auprc = average_precision_score(label_b_u, 1. / (1. + np.exp(-preds_matrix[:, i])))
                outfile.write("%s model: auc:%.6f auprc:%.6f\n" % (
                    model_files[i].decode().split("/")[-1].replace('_model.pkl', ''), auc, auprc))
                precision, recall, thresholds = precision_recall_curve(label_b_u, preds, pos_label=1)
                df_temp = pd.DataFrame(None)
                df_temp["precision"] = precision
                df_temp["recall"] = recall
                df_temp["model"] = model_files[i].decode().split("/")[-1]
                temp.append(df_temp.sample(n=100000))
            df_plot = pd.concat(temp, ignore_index=True)
            plt.figure(figsize=(8, 6))
            ax = sns.lineplot(x="recall", y="precision", data=df_plot, hue='model', palette="tab10")
            ax.set_title("%s in %s" % (self.training_tf_name, cell_line))
            ax.set_xlabel("Recall")
            ax.set_ylabel("Precision")
            ax.get_figure().savefig("%s/%s_%s_PRC.pdf" % (dir_out, self.training_tf_name, cell_line))

        with open("%s/%s.%s_confusion_matrix.hyperopt.txt" % (dir_out, self.training_tf_name, cell_line), "w") as outfile:
            for i in range(preds_matrix.shape[1]):
                one_preds = 1. / (1. + np.exp(-preds_matrix[:, i]))
                cutoff = 0.5
                true_positive = np.where((one_preds >= cutoff) & (label_b_u == 1))[0]
                false_positive = np.where((one_preds >= cutoff) & (label_b_u == 0))[0]
                false_negative = np.where((one_preds < cutoff) & (label_b_u == 1))[0]
                outfile.write("%s model: all_regions:%d true_positive:%d false_positive:%d false_negative:%d\n" % (
                    model_files[i].decode().split("/")[-1].replace('_model.pkl', ''), len(one_preds),
                    len(true_positive), len(false_positive), len(false_negative)))
                df = pd.DataFrame(None)
                df["chrom"] = np.hstack((chroms[true_positive], chroms[false_positive], chroms[false_negative]))
                df["start"] = np.hstack((starts[true_positive], starts[false_positive], starts[false_negative]))
                df["preds"] = np.hstack(
                    (one_preds[true_positive], one_preds[false_positive], one_preds[false_negative]))
                df["label"] = np.hstack(
                    (label_b_u[true_positive], label_b_u[false_positive], label_b_u[false_negative]))
                df["class"] = ['true_positive'] * len(true_positive) + ['false_positive'] * len(false_positive) + [
                    'false_negative'] * len(false_negative)
                df.to_csv(
                    "%s/df.%s.regions.xls" % (
                        dir_out, model_files[i].decode().split("/")[-1].replace('_model.pkl', '')),
                    sep="\t", header=True, index=False)

    @staticmethod
    def bins_with_dnase_peak(cell_line):
        f = os.popen(
            'zcat /n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/DNASE/'
            'peaks/conservative/DNASE.%s.conservative.narrowPeak.gz | bedtools sort -i stdin | '
            'bedtools intersect -c -b stdin -a /n/scratchlfs/xiaoleliu_lab/Jingyu/'
            'impute_cistrome/ENCODE_DREAM/annotations/test_regions.blacklistfiltered.sorted.bed' % cell_line)
        df = pd.read_csv(f, sep="\t", header=None)
        df.columns = ['chr', 'start', 'end', 'peak_num']
        return df

    def evaluation_leafs(self, cell_line, lightgbm_preds_files_path=None, dir_out=None):
        if dir_out is None:
            dir_out = "./train/%s/evaluations/" % self.training_tf_name
        dic_model_cell_type_cov = {}
        chrom = "chr19"
        with h5py.File(
                '%s/%s.%s.%s_pred_leafs.h5' % (lightgbm_preds_files_path, self.training_tf_name, cell_line, chrom),
                "r") as infile:
            model_files = infile['model_files'][...]
        model_files = list(map(lambda x: x.decode(), model_files))
        df_temp = self.bins_with_dnase_peak(cell_line)
        for model_file in model_files:
            list_pred_leafs = []
            list_labels = []
            for chrom in self.chrom_all:
                with h5py.File(
                        '%s/%s.%s.%s_pred_leafs.h5' % (
                                lightgbm_preds_files_path, self.training_tf_name, cell_line, chrom),
                        "r") as infile:
                    pred_leafs = infile["%s/%s" % (model_file.split("/")[-1], cell_line)][...]
                    labels = np.array(df_temp.loc[df_temp['chr'] == chrom, :]['peak_num'])
                    list_pred_leafs.append(pred_leafs)
                    list_labels.append(labels)
            labels = np.hstack(list_labels)
            pred_leafs = np.vstack(list_pred_leafs)
            ignore_index = np.where(labels == 0)[0]
            pred_leafs = np.delete(pred_leafs, ignore_index, axis=0)
            scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
            scaler.fit(pred_leafs)
            pred_leafs_scaled = scaler.transform(pred_leafs)
            dic_model_cell_type_cov[model_file.split("/")[-1]] = {}
            dic_model_cell_type_cov[model_file.split("/")[-1]][cell_line] = np.cov(pred_leafs_scaled.T)

        for model_file in model_files:
            training_cell_line = model_file.split("/")[-1].split('.')[1]
            df_temp = self.bins_with_dnase_peak(training_cell_line)
            list_pred_leafs = []
            list_labels = []
            for chrom in self.chrom_all:
                with h5py.File(
                        '%s/%s.%s.%s_pred_leafs.h5' % (
                                lightgbm_preds_files_path, self.training_tf_name, cell_line, chrom),
                        "r") as infile:
                    pred_leafs = infile["%s/%s" % (model_file.split("/")[-1], training_cell_line)][...]
                    labels = np.array(df_temp.loc[df_temp['chr'] == chrom, :]['peak_num'])
                    list_pred_leafs.append(pred_leafs)
                    list_labels.append(labels)
            labels = np.hstack(list_labels)
            pred_leafs = np.vstack(list_pred_leafs)
            ignore_index = np.where(labels == 0)[0]
            pred_leafs = np.delete(pred_leafs, ignore_index, axis=0)
            scaler = StandardScaler(copy=True, with_mean=True, with_std=True)
            scaler.fit(pred_leafs)
            pred_leafs_scaled = scaler.transform(pred_leafs)
            dic_model_cell_type_cov[model_file.split("/")[-1]][training_cell_line] = np.cov(pred_leafs_scaled.T)

        with open("%s/%s.%s_models_fro_norm.txt" % (dir_out, self.training_tf_name, cell_line), "w") as outfile:
            for model_file in model_files:
                training_cell_line = model_file.split("/")[-1].split('.')[1]
                outfile.write("model_name:%s fro_norm:%.6f\n" % (
                    model_file.split("/")[-1], LA.norm(
                        dic_model_cell_type_cov[model_file.split("/")[-1]][training_cell_line] -
                        dic_model_cell_type_cov[model_file.split("/")[-1]][cell_line], 'fro')))

    def evaluation_autoencoder(self, cell_line, lightgbm_preds_files_path=None, dir_out=None):
        if dir_out is None:
            dir_out = "./train/%s/evaluations/" % self.training_tf_name
        df_test_regions_label = pd.read_csv(
            "%s/%s.%s" % (
                self.config['test_cell_types_regions_label_path'], self.training_tf_name,
                self.config['test_cell_types_regions_label_name']), sep="\t", header=0)
        list_preds_binary = []
        # list_preds_binary_2 = []
        list_labels = []
        list_preds_matrix = []
        for chrom in self.chrom_all:
            with h5py.File(
                    '%s/%s.%s.%s_preds.autoencoder.h5' % (
                            lightgbm_preds_files_path, self.training_tf_name, cell_line, chrom),
                    "r") as infile:
                model_files = infile['model_files'][...]
                preds = infile['preds'][...]
                labels = np.array(df_test_regions_label.loc[df_test_regions_label['chr'] == chrom, :][cell_line])
                list_preds_matrix.append(preds)
                preds_binary = np.mean(1. / (1. + np.exp(-preds)), axis=1)
                # preds_binary_2 = 1. / (1. + np.exp(-np.mean(preds, axis=1)))
                list_preds_binary.append(preds_binary)
                # list_preds_binary_2.append(preds_binary_2)
                list_labels.append(labels)
        labels = np.hstack(list_labels)
        preds = np.hstack(list_preds_binary)
        # preds_2 = np.hstack(list_preds_binary_2)
        preds_matrix = np.vstack(list_preds_matrix)
        ignore_index = np.where(labels == "A")[0]
        preds_matrix = np.delete(preds_matrix, ignore_index, axis=0)
        preds = np.delete(preds, ignore_index, axis=0)
        label_b_u = np.delete(labels, ignore_index, axis=0)
        label_b_u = np.array(list(map(lambda x: 1 if x == "B" else 0, label_b_u)))
        with open("%s/%s.%s_performance.autoencoder.txt" % (dir_out, self.training_tf_name, cell_line), "w") as outfile:
            fpr, tpr, thresholds = metrics.roc_curve(label_b_u, preds, pos_label=1)
            auc = metrics.auc(fpr, tpr)
            auprc = average_precision_score(label_b_u, preds)
            outfile.write("average model: auc:%.6f auprc:%.6f\n" % (auc, auprc))
            for i in range(preds_matrix.shape[1]):
                fpr, tpr, thresholds = metrics.roc_curve(label_b_u, preds_matrix[:, i], pos_label=1)
                auc = metrics.auc(fpr, tpr)
                auprc = average_precision_score(label_b_u, preds_matrix[:, i])
                outfile.write("%s model: auc:%.6f auprc:%.6f\n" % (
                    model_files[i].decode().split("/")[-1].replace('_model.pkl', ''), auc, auprc))


def prepare_lightgbm_binary_data_dnase_median(self, cell_line, chrom_set_name, dir_out, reference=None,
                                              selected_bin_index_file=None):
    subset_index = -1
    chrom = "chr19"
    # TODO change to 50bp or 100bp
    with h5py.File("%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_all_cell_types.h5" % (
            self.dnase_feature_path, chrom), "r") as infile:
        needed_feature_names = list(infile['feature_names'][...])
        needed_feature_names = list(map(lambda x: "median_" + x.decode('UTF-8'), needed_feature_names))
    list_scores = []
    labels = []
    chrom_set = self.chrom_sets[chrom_set_name]
    for chrom in chrom_set:
        df_temp = self.df_all_regions_label.loc[self.df_all_regions_label['chr'] == chrom, :]
        with h5py.File("%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_all_cell_types.h5" % (
                self.dnase_feature_path, chrom),
                       "r") as infile:
            if selected_bin_index_file is None:
                scores = infile[chrom][:, :, :]
            else:
                selected_bin_index = np.load(selected_bin_index_file)
                scores = infile[chrom][:, selected_bin_index, :]
            samples = list(infile['samples'][...])
            samples = list(map(lambda x: x.decode('UTF-8'), samples))
            for cell_line_name in samples:
                cell_line_index = np.where(np.array(samples) == cell_line_name)[0][0]
                with open("%s/%s_variable_bp_quantile_map.pkl" % (self.quantile_transformer_path, cell_line_name),
                          'rb') as fin:
                    qt = pickle.load(fin, encoding='latin1')
                _ = qt.transform(scores[cell_line_index, :, :])
            scores_median = np.median(scores, axis=0)
        ignore_index = np.where(df_temp[cell_line] == "A")[0]
        scores_median = np.delete(scores_median, ignore_index, axis=0)
        label_b_u = np.delete(np.array(df_temp[cell_line]), ignore_index, axis=0)
        temp_label = list(map(lambda x: 1 if x == "B" else 0, label_b_u))
        labels += temp_label
        list_scores.append(scores_median)
        # print(cell_line, chrom, subset_index)
    all_score = np.vstack(list_scores)

    if reference is not None:
        reference = lgb.Dataset(glob.glob("%s/median/lightGBM.*.*.%d.bin" % (reference, subset_index))[0])
    train_data = lgb.Dataset(all_score, feature_name=list(needed_feature_names), label=labels, reference=reference)
    train_data.save_binary("%s/median/lightGBM.%s.%s.%d.bin" % (dir_out, cell_line, chrom_set_name, subset_index))


def merge_lightgbm_binary_data_median(self, cell_line, chrom_set_name, dir_out):
    # TODO change to 50bp or 100bp
    # with h5py.File("%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_all_cell_types.h5" % (
    #         self.dnase_feature_path, chrom), "r") as infile:
    #     all_feature_names += list(infile['feature_names'][...])
    # chrom = "chr22"
    for cell_line in [cell_line]:
        subset_index = -1
        train_data_all = lgb.Dataset("%s/lightGBM_all.%s.%s.bin" % (dir_out, cell_line, chrom_set_name)).construct()
        train_data = lgb.Dataset("%s/median/lightGBM.%s.%s.%d.bin" %
                                 (dir_out, cell_line, chrom_set_name, subset_index)).construct()

        train_data_all.add_features_from(train_data)
        # print(subset_index)
        train_data_all.save_binary("%s/median/lightGBM_all.%s.%s.bin" % (dir_out, cell_line, chrom_set_name))


if __name__ == '__main__':
    fire.Fire(LightGBMModel)
