#!/usr/bin/env python3
# DNase features on 5' cuts
import os
import fire
import h5py
import pickle
import pyBigWig
import pandas as pd
import numpy as np
from sklearn.preprocessing import QuantileTransformer

import yaml

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


class DNaseFeatures(object):

    def __init__(self, config_file):
        with open(config_file, "r") as infile:
            config = yaml.load(infile, Loader=Loader)
        self.config = config
        self.cell_types = config['cell_types']
        self.chrom_all = config['chrom_all']
        self.chrom_sets = config['chrom_sets']
        self.batch = config['batch']
        self.dic_chrom_length = {}
        with open(config['chrom_size_file'], "r") as infile:
            for line in infile:
                line = line.strip().split("\t")
                if line[0] in self.chrom_all:
                    self.dic_chrom_length[line[0]] = int(line[1])
        # no header
        self.df_all_regions = pd.read_csv(config['regions_all_file'], sep="\t", header=None)
        self.df_all_regions.columns = ['chr', 'start', 'stop']
        self.ATAC_long_short = config['ATAC_long_short']
        self.CA_cluster_mode = config['CA_cluster_mode']
        if self.CA_cluster_mode:
            self.cluster_info_file = config['cluster_info_file']
            self.df_cluster_info_file = pd.read_csv(config['cluster_info_file'], sep=",", header=0)
            self.selected_cluster = config['selected_cluster']


    # per cell_line per chrom
    def generate_dnase_features_from_bw_to_h5(self, chrom, cell_line, dnase_bw_file_path, outfile_path):

        df_temp = self.df_all_regions[self.df_all_regions['chr'] == chrom].copy()
        # 50bp resolution [-700,700] + -
        temp_50 = np.zeros((np.floor(self.dic_chrom_length[chrom] / 50).astype(int) - 3, 14 * 2 * 2))
        # 100bp resolution [-700,700] + -
        temp_100 = np.zeros((np.floor(self.dic_chrom_length[chrom] / 50).astype(int) - 3, 7 * 2 * 2))
        # 25 50 100 *4
        temp_var = np.zeros((np.floor(self.dic_chrom_length[chrom] / 50).astype(int) - 3, 12 * 2 * 2))

        all_bins_50 = list(range(0, 1400 + 50, 50))

        all_bins_100 = list(range(0, 1400 + 100, 100))

        bins = list(range(0, 100, 25)) + list(range(100, 300, 50)) + list(range(300, 700 + 100, 100))
        all_bins_var = list(map(lambda x: 700 - x, bins[::-1][:-1])) + list(map(lambda x: 700 + x, bins))
        # for ind_cell_line,cell_line in enumerate(cell_lines):
        # print(cell_line, chrom)
        for strand in ['forward', 'reverse']:
            # print(strand)
            if strand == 'forward':
                feature_index_50 = 0
                feature_index_100 = 0
                feature_index_var = 0
            else:
                feature_index_50 = int(temp_50.shape[1] / 2)
                feature_index_100 = int(temp_100.shape[1] / 2)
                feature_index_var = int(temp_var.shape[1] / 2)

            bw = pyBigWig.open("%s/DNASE.%s.merge.binSize.1.%s.bw" % (dnase_bw_file_path, cell_line, strand), 'r')

            n = 0
            for j in range(0, self.dic_chrom_length[chrom], self.batch):
                values = bw.values(chrom, max(0, j - 600), min(j + self.batch + 800, self.dic_chrom_length[chrom]),
                                   numpy=True)
                #                     values=bw.values(chrom, 0, min(j+batch+800,dic_chrom_length[chrom]),numpy=True)
                values = np.nan_to_num(values)
                values = np.pad(values, (max(0, -j + 600), max(0, j + self.batch + 800 - self.dic_chrom_length[chrom])),
                                'constant',
                                constant_values=0)
                # print("1 batch")
                for i in range(j, min(j + self.batch, self.dic_chrom_length[chrom]), 50):

                    if i + 200 >= self.dic_chrom_length[chrom]:
                        continue
                    start = i - j + 600 - 600
                    end = i - j + 600 + 800
                    temp_1400 = values[start:end]
                    if np.max(temp_1400) == 0:
                        temp_50[n, feature_index_50:(feature_index_50 + int(temp_50.shape[1] / 2))] = 0
                        temp_100[n, feature_index_100:(feature_index_100 + int(temp_100.shape[1] / 2))] = 0
                        temp_var[n, feature_index_var:(feature_index_var + int(temp_var.shape[1] / 2))] = 0
                    else:
                        temp_1400_cut_index = np.arange(1400)[temp_1400.astype(bool)]
                        #
                        temp_50[n, feature_index_50:(feature_index_50 + int(temp_50.shape[1] / 2))] \
                            = np.histogram(temp_1400_cut_index, bins=all_bins_50, density=True)[0] * len(
                            temp_1400_cut_index)

                        temp_100[n, feature_index_100:(feature_index_100 + int(temp_100.shape[1] / 2))] \
                            = np.histogram(temp_1400_cut_index, bins=all_bins_100, density=True)[0] * len(
                            temp_1400_cut_index)

                        temp_var[n, feature_index_var:(feature_index_var + int(temp_var.shape[1] / 2))] \
                            = np.histogram(temp_1400_cut_index, bins=all_bins_var, density=True)[0] * len(
                            temp_1400_cut_index)
                    n += 1

        with h5py.File("%s/DNASE_bam_5_mer_50bp_all_samples_lightGBM_%s.%s.h5" % (outfile_path, chrom, cell_line),
                       "w") as outfile:
            feature_names = []
            for strand in ['forward', 'reverse']:
                for ind_feature_name, feature_name in enumerate(all_bins_50[:-1]):
                    feature_names.append(
                        "DNase_5_end_%d_%d_%s" % (feature_name - 700, all_bins_50[ind_feature_name + 1] - 700, strand))

            outfile.create_dataset("%s_starts" % chrom, data=df_temp['start'].tolist(), shape=(df_temp.shape[0],),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            outfile.create_dataset(chrom, data=temp_50[(df_temp['start'] / 50).values.astype(int), :], dtype=np.float32,
                                   shape=(df_temp.shape[0], 14 * 2 * 2),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=4)
            outfile.create_dataset("feature_names", data=np.array(feature_names, dtype='S'),
                                   shape=(len(feature_names),),
                                   dtype='S200', compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)

        with h5py.File("%s/DNASE_bam_5_mer_100bp_all_samples_lightGBM_%s.%s.h5" % (outfile_path, chrom, cell_line),
                       "w") as outfile:
            feature_names = []
            for strand in ['forward', 'reverse']:
                for ind_feature_name, feature_name in enumerate(all_bins_100[:-1]):
                    feature_names.append(
                        "DNase_5_end_%d_%d_%s" % (feature_name - 700, all_bins_100[ind_feature_name + 1] - 700, strand))

            outfile.create_dataset("%s_starts" % chrom, data=df_temp['start'].tolist(), shape=(df_temp.shape[0],),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            outfile.create_dataset(chrom, data=temp_100[(df_temp['start'] / 50).values.astype(int), :],
                                   dtype=np.float32,
                                   shape=(df_temp.shape[0], 7 * 2 * 2),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=4)
            outfile.create_dataset("feature_names", data=np.array(feature_names, dtype='S'),
                                   shape=(len(feature_names),),
                                   dtype='S200', compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)

        with h5py.File(
                "%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s.%s.h5" % (outfile_path, chrom, cell_line),
                "w") as outfile:
            feature_names = []
            for strand in ['forward', 'reverse']:
                for ind_feature_name, feature_name in enumerate(all_bins_var[:-1]):
                    feature_names.append(
                        "DNase_5_end_%d_%d_%s" % (feature_name - 700, all_bins_var[ind_feature_name + 1] - 700, strand))
            outfile.create_dataset("%s_starts" % chrom, data=df_temp['start'].tolist(), shape=(df_temp.shape[0],),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            outfile.create_dataset(chrom, data=temp_var[(df_temp['start'] / 50).values.astype(int), :],
                                   dtype=np.float32,
                                   shape=(df_temp.shape[0], 12 * 2 * 2),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=4)
            outfile.create_dataset("feature_names", data=np.array(feature_names, dtype='S'),
                                   shape=(len(feature_names),),
                                   dtype='S200', compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)

    # @staticmethod
    def combine_dnase_features_h5_for_all_cell_types(self, chrom, outfile_path):
        if not self.CA_cluster_mode:
            if self.ATAC_long_short:
                cell_types = ["%s_%s" % (cell_type, frag_len) for cell_type in self.cell_types
                              for frag_len in ["short", 'long']]
            else:
                cell_types = self.cell_types
        else:
            if not self.selected_cluster:
                cell_types = list(self.df_cluster_info_file['cluster'].unique())
            else:
                cell_types = self.selected_cluster
        # with h5py.File("%s/DNASE_bam_5_mer_50bp_all_samples_lightGBM_%s_all_cell_types.h5" % (outfile_path, chrom),
        #                "w") as outfile:
        #     temp = []
        #
        #     for cell_line in cell_types:
        #         # print(cell_line)
        #         if not self.CA_cluster_mode:
        #             with h5py.File(
        #                     "%s/DNASE_bam_5_mer_50bp_all_samples_lightGBM_%s.%s.h5" % (outfile_path, chrom, cell_line),
        #                     "r") as infile:
        #                 starts = infile["%s_starts" % chrom][...]
        #                 feature_names = infile["feature_names"][...]
        #                 temp.append(infile[chrom][...])
        #         else:
        #             cluster_sum = None
        #             for cell_line_id in self.df_cluster_info_file[self.df_cluster_info_file['cluster']==int(cell_line)]['sample'].values:
        #                 h5_file_in = "%s/DNASE_bam_5_mer_50bp_all_samples_lightGBM_%s.%s.h5" % (outfile_path, chrom,
        #                                                                                      str(cell_line_id))
        #                 if os.path.exists(h5_file_in):
        #                     with h5py.File(h5_file_in,"r") as infile:
        #                         starts = infile["%s_starts" % chrom][...]
        #                         feature_names = infile["feature_names"][...]
        #                         if cluster_sum is None:
        #                             cluster_sum = infile[chrom][...]
        #                         else:
        #                             cluster_sum += infile[chrom][...]
        #             temp.append(cluster_sum)
        #     outfile.create_dataset("%s_starts" % chrom, data=starts, shape=(len(starts),),
        #                            compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
        #     outfile.create_dataset("feature_names", data=feature_names,
        #                            shape=(len(feature_names),),
        #                            dtype='S200', compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
        #     outfile.create_dataset("samples", data=np.array(cell_types, dtype='S'),
        #                            shape=(len(cell_types),),
        #                            dtype='S200', compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
        #     data_temp = np.zeros((len(cell_types), len(starts), 14 * 2 * 2))
        #     for ind_cell_line, cell_line in enumerate(cell_types):
        #         data_temp[ind_cell_line, :, :] = temp[ind_cell_line]
        #     outfile.create_dataset(chrom, data=data_temp, dtype=np.float32,
        #                            shape=data_temp.shape,
        #                            compression='gzip', shuffle=True, fletcher32=True, compression_opts=4)
        #
        # with h5py.File("%s/DNASE_bam_5_mer_100bp_all_samples_lightGBM_%s_all_cell_types.h5" % (outfile_path, chrom),
        #                "w") as outfile:
        #     temp = []
        #     for cell_line in cell_types:
        #         # print(cell_line)
        #         if not self.CA_cluster_mode:
        #             with h5py.File(
        #                     "%s/DNASE_bam_5_mer_100bp_all_samples_lightGBM_%s.%s.h5" % (outfile_path, chrom, cell_line),
        #                     "r") as infile:
        #                 starts = infile["%s_starts" % chrom][...]
        #                 feature_names = infile["feature_names"][...]
        #                 temp.append(infile[chrom][...])
        #         else:
        #             cluster_sum = None
        #             for cell_line_id in self.df_cluster_info_file[self.df_cluster_info_file['cluster']==int(cell_line)]['sample'].values:
        #                 h5_file_in = "%s/DNASE_bam_5_mer_100bp_all_samples_lightGBM_%s.%s.h5" % (outfile_path, chrom,
        #                                                                                      str(cell_line_id))
        #                 if os.path.exists(h5_file_in):
        #                     with h5py.File(h5_file_in,"r") as infile:
        #                         starts = infile["%s_starts" % chrom][...]
        #                         feature_names = infile["feature_names"][...]
        #                         if cluster_sum is None:
        #                             cluster_sum = infile[chrom][...]
        #                         else:
        #                             cluster_sum += infile[chrom][...]
        #             temp.append(cluster_sum)
        #
        #     outfile.create_dataset("%s_starts" % chrom, data=starts, shape=(len(starts),),
        #                            compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
        #     outfile.create_dataset("feature_names", data=feature_names,
        #                            shape=(len(feature_names),),
        #                            dtype='S200', compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
        #     outfile.create_dataset("samples", data=np.array(cell_types, dtype='S'),
        #                            shape=(len(cell_types),),
        #                            dtype='S200', compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
        #     data_temp = np.zeros((len(cell_types), len(starts), 7 * 2 * 2))
        #     for ind_cell_line, cell_line in enumerate(cell_types):
        #         data_temp[ind_cell_line, :, :] = temp[ind_cell_line]
        #     outfile.create_dataset(chrom, data=data_temp, dtype=np.float32,
        #                            shape=data_temp.shape,
        #                            compression='gzip', shuffle=True, fletcher32=True, compression_opts=4)

        with h5py.File(
                "%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_all_cell_types.h5" % (outfile_path, chrom),
                "w") as outfile:
            temp = []
            for cell_line in cell_types:
                # print(cell_line)
                if not self.CA_cluster_mode:
                    with h5py.File(
                            "%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s.%s.h5" % (
                                    outfile_path, chrom, cell_line),
                            "r") as infile:
                        starts = infile["%s_starts" % chrom][...]
                        feature_names = infile["feature_names"][...]
                        temp.append(infile[chrom][...])
                else:
                    cluster_sum = None
                    for cell_line_id in self.df_cluster_info_file[self.df_cluster_info_file['cluster']==int(cell_line)]['sample'].values:
                        h5_file_in = "%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s.%s.h5" % (outfile_path, chrom,
                                                                                             str(cell_line_id))
                        if os.path.exists(h5_file_in):
                            with h5py.File(h5_file_in,"r") as infile:
                                starts = infile["%s_starts" % chrom][...]
                                feature_names = infile["feature_names"][...]
                                if cluster_sum is None:
                                    cluster_sum = infile[chrom][...]
                                else:
                                    cluster_sum += infile[chrom][...]
                    temp.append(cluster_sum)
            outfile.create_dataset("%s_starts" % chrom, data=starts, shape=(len(starts),),
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            outfile.create_dataset("feature_names", data=feature_names,
                                   shape=(len(feature_names),),
                                   dtype='S200', compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            outfile.create_dataset("samples", data=np.array(cell_types, dtype='S'),
                                   shape=(len(cell_types),),
                                   dtype='S200', compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
            data_temp = np.zeros((len(cell_types), len(starts), 12 * 2 * 2))
            for ind_cell_line, cell_line in enumerate(cell_types):
                data_temp[ind_cell_line, :, :] = temp[ind_cell_line]
            outfile.create_dataset(chrom, data=data_temp, dtype=np.float32,
                                   shape=data_temp.shape,
                                   compression='gzip', shuffle=True, fletcher32=True, compression_opts=4)

    def generate_quantile_transformer(self, cell_line, dir_dnase_feature, dir_out):
        list_scores = []
        # df_final_all_chr.shape
        # (60519747, 4)
        df_all_regions_sampled = self.df_all_regions.sample(n=1000000, replace=False, random_state=6)
        for chrom in self.chrom_all:
            df_temp = df_all_regions_sampled.loc[df_all_regions_sampled['chr'] == chrom, :]
            filename = '%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_all_cell_types.h5' % (
                dir_dnase_feature, chrom)
            with h5py.File(filename, "r") as infile:
                starts = infile['%s_starts' % chrom][...]
                selected_start = np.array(df_temp['start'].tolist())
                selected_index = np.where(np.isin(starts, selected_start))[0]
                samples = list(infile['samples'][...])
                samples = list(map(lambda x: x.decode('UTF-8'), samples))
                cell_line_index = np.where(np.array(samples) == str(cell_line))[0][0]
                scores_dnase = infile[chrom][cell_line_index, selected_index, :]
            scores = scores_dnase
            # ignore_index = np.where(df_temp[cell_line] == "A")[0]
            # scores = np.delete(scores, ignore_index, axis=0)
            # label_B_U = np.delete(np.array(df_temp[cell_line]), ignore_index, axis=0)
            # temp_label = map(lambda x: 1 if x == "B" else 0, label_B_U)
            # labels += temp_label
            list_scores.append(scores)
            # print(cell_line, chrom)
        all_score = np.vstack(list_scores)
        # print('all_score', all_score.shape)
        qt = QuantileTransformer(n_quantiles=1000, random_state=6, output_distribution='uniform',
                                 ignore_implicit_zeros=False, subsample=100000, copy=False)
        qt.fit(all_score)
        with open("%s/%s_variable_bp_quantile_map.pkl" % (dir_out, cell_line), 'wb') as fout:
            pickle.dump(qt, fout)

        # # 50bp
        # list_scores = []
        # for chrom in self.chrom_all:
        #     df_temp = df_all_regions_sampled.loc[df_all_regions_sampled['chr'] == chrom, :]
        #     filename = '%s/DNASE_bam_5_mer_50bp_all_samples_lightGBM_%s_all_cell_types.h5' % (dir_dnase_feature, chrom)
        #     with h5py.File(filename, "r") as infile:
        #         starts = infile['%s_starts' % chrom][...]
        #         selected_start = np.array(df_temp['start'].tolist())
        #         selected_index = np.where(np.isin(starts, selected_start))[0]
        #         samples = list(infile['samples'][...])
        #         samples = list(map(lambda x: x.decode('UTF-8'), samples))
        #         cell_line_index = np.where(np.array(samples) == cell_line)[0][0]
        #         scores_dnase = infile[chrom][cell_line_index, selected_index, :]
        #     scores = scores_dnase
        #     # ignore_index = np.where(df_temp[cell_line] == "A")[0]
        #     # scores = np.delete(scores, ignore_index, axis=0)
        #     # label_B_U = np.delete(np.array(df_temp[cell_line]), ignore_index, axis=0)
        #     # temp_label = map(lambda x: 1 if x == "B" else 0, label_B_U)
        #     # labels += temp_label
        #     list_scores.append(scores)
        #     # print(cell_line, chrom)
        # all_score = np.vstack(list_scores)
        # # print('all_score', all_score.shape)
        # qt = QuantileTransformer(n_quantiles=1000, random_state=6, output_distribution='uniform',
        #                          ignore_implicit_zeros=False, subsample=1000000, copy=False)
        # qt.fit(all_score)
        # with open("%s/%s_50bp_quantile_map.pkl" % (dir_out, cell_line), 'wb') as fout:
        #     pickle.dump(qt, fout)
        #
        # # 100bp
        # list_scores = []
        # for chrom in self.chrom_all:
        #     df_temp = df_all_regions_sampled.loc[df_all_regions_sampled['chr'] == chrom, :]
        #     filename = '%s/DNASE_bam_5_mer_100bp_all_samples_lightGBM_%s_all_cell_types.h5' % (dir_dnase_feature, chrom)
        #     with h5py.File(filename, "r") as infile:
        #         starts = infile['%s_starts' % chrom][...]
        #         selected_start = np.array(df_temp['start'].tolist())
        #         selected_index = np.where(np.isin(starts, selected_start))[0]
        #         samples = list(infile['samples'][...])
        #         samples = list(map(lambda x: x.decode('UTF-8'), samples))
        #         cell_line_index = np.where(np.array(samples) == cell_line)[0][0]
        #         scores_dnase = infile[chrom][cell_line_index, selected_index, :]
        #     scores = scores_dnase
        #     # ignore_index = np.where(df_temp[cell_line] == "A")[0]
        #     # scores = np.delete(scores, ignore_index, axis=0)
        #     # label_B_U = np.delete(np.array(df_temp[cell_line]), ignore_index, axis=0)
        #     # temp_label = map(lambda x: 1 if x == "B" else 0, label_B_U)
        #     # labels += temp_label
        #     list_scores.append(scores)
        #     # print(cell_line, chrom)
        # all_score = np.vstack(list_scores)
        # # print('all_score', all_score.shape)
        # qt = QuantileTransformer(n_quantiles=1000, random_state=6, output_distribution='uniform',
        #                          ignore_implicit_zeros=False, subsample=1000000, copy=False)
        # qt.fit(all_score)
        # with open("%s/%s_100bp_quantile_map.pkl" % (dir_out, cell_line), 'wb') as fout:
        #     pickle.dump(qt, fout)

    @staticmethod
    def generate_dnase_feature_median(chrom, dir_dnase_feature, dir_quantile_transformer, dir_out,
                                      selected_bin_index_file=None):
        # TODO change to 50bp or 100bp
        with h5py.File("%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_all_cell_types.h5" % (
                dir_dnase_feature, chrom), "r") as infile:
            feature_names = list(infile['feature_names'][...])
            feature_names = list(map(lambda x: "median_" + x.decode('UTF-8'), feature_names))
            starts = list(infile["%s_starts" % chrom][...])
            if selected_bin_index_file is None:
                scores = infile[chrom][:, :, :]
            else:
                selected_bin_index = np.load(selected_bin_index_file)
                scores = infile[chrom][:, selected_bin_index, :]
            samples = list(infile['samples'][...])
            samples = list(map(lambda x: x.decode('UTF-8'), samples))
            for cell_line_name in samples:
                cell_line_index = np.where(np.array(samples) == cell_line_name)[0][0]
                with open("%s/%s_variable_bp_quantile_map.pkl" % (dir_quantile_transformer, cell_line_name),
                          'rb') as fin:
                    qt = pickle.load(fin, encoding='latin1')
                _ = qt.transform(scores[cell_line_index, :, :])
            scores_median = np.median(scores, axis=0)
            with h5py.File("%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_median.h5" % (
                    dir_out, chrom), "w") as outfile:
                outfile.create_dataset("%s_starts" % chrom, data=starts, shape=(len(starts),),
                                       compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
                outfile.create_dataset("feature_names", data=np.array(feature_names, dtype='S'),
                                       shape=(len(feature_names),),
                                       dtype='S200', compression='gzip', shuffle=True, fletcher32=True,
                                       compression_opts=9)
                outfile.create_dataset(chrom, data=scores_median, dtype=np.float32,
                                       shape=scores_median.shape,
                                       compression='gzip', shuffle=True, fletcher32=True, compression_opts=4)


if __name__ == '__main__':
    fire.Fire(DNaseFeatures)
