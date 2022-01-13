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
import pandas as pd
import numpy as np

from sklearn.metrics import average_precision_score
from sklearn.preprocessing import QuantileTransformer


#python3 scripts/merge_train_test.py merge_ave \
#--input_h5file='hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_chr19_all_cell_types.SQTKNNN.h5' \
#--output_h5file='hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_chr19_all_cell_types.SQTKNNN.ave.h5' \
#--chrom='chr19' \
#--train_id_file='HepG2_Epithelium_Liver.train.id.list.txt'
#--test_id_file='MCF-7_Epithelium_Breast.test.id.list.txt'


def merge_ave(input_h5file, output_h5file, chrom, train_id_file, test_id_file):

	print('read h5 file')
	with h5py.File(input_h5file, "r") as infile:
		print('read h5 file score')
		score100 = infile[chrom][...]
		print('read h5 file samples')
		#samples = np.array(infile['samples'][...]).astype(np.int64)
		samples = np.array(infile['samples'][...]).astype(str)
		print('read h5 file featurenames')
		features = infile['feature_names'][...]
		print('read h5 file starts')
		starts = infile['%s_starts' % chrom][...]

	infile.close()


	### read samples & train & test ids 
	train_id = np.array(pd.read_csv(train_id_file, sep="\t", header=None))
	test_id = np.array(pd.read_csv(test_id_file, sep="\t", header=None))

	print(train_id)
	print(test_id)
	### get where are the train & test ids 
	train_id_pos = np.intersect1d(samples, train_id, return_indices=True)[1]
	test_id_pos = np.intersect1d(samples, test_id, return_indices=True)[1]

	print(train_id_pos)
	print(test_id_pos)

	### get average of train & test
	print(score100[train_id_pos,:,:])
	print(score100)
	print(score100.shape)
	print((score100[train_id_pos,:,:]).shape)
	score100_train = np.mean(score100[train_id_pos,:,:], axis=0).reshape(1,score100.shape[1],score100.shape[2])
	print((score100[test_id_pos,:,:]).shape)
	score100_test = np.mean(score100[test_id_pos,:,:], axis=0).reshape(1,score100.shape[1],score100.shape[2])
	score100_train_test_ave = np.concatenate((score100_train, score100_test), axis=0)
	print(score100_train_test_ave.shape)
	samples_new = np.array(['train', 'test']).astype('bytes_')

	score100_train_test_ave[np.isnan(score100_train_test_ave)] = 0	
	print(output_h5file)
	print(score100_train_test_ave)
	print(np.mean(score100_train_test_ave))
	print(np.max(score100_train_test_ave))

	### write output
	with h5py.File(output_h5file, "w") as outfile:
		print('score')
		outfile.create_dataset(chrom, data=np.array(score100_train_test_ave), compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
		print('starts')
		outfile.create_dataset('%s_starts' % chrom, data=starts, compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
		print('feature_names')
		outfile.create_dataset('feature_names', data=features, compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
		print('samples')
		outfile.create_dataset('samples', data=samples_new, compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)

	outfile.close()



if __name__ == '__main__':
	fire.Fire(merge_ave)


