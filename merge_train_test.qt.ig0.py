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

	with h5py.File(input_h5file, "r") as infile:
		score100 = infile[chrom][...]
		samples = np.array(infile['samples'][...]).astype(np.int64)
		features = infile['feature_names'][...]
		starts = infile['%s_starts' % chrom][...]

	infile.close()

	#qt = QuantileTransformer(n_quantiles=1000, random_state=6, output_distribution='uniform', ignore_implicit_zeros=False, subsample=100000, copy=False)
	#all_score = score100[0,:,:]
	#for i in range(1,score100.shape[0]):
	#	all_score = np.concatenate((all_score, score100[i,:,:]), axis=1)
	#qt.fit(all_score)
	#new_data = qt.fit_transform(all_score)

		

	#print(all_score.shape)
	#score100_new = score100
	#qt = QuantileTransformer(n_quantiles=1000, random_state=6, output_distribution='uniform', ignore_implicit_zeros=False, subsample=100000, copy=False)
	#qt.fit(score100[0,:,:])
	#for i in range(0,score100.shape[0]):
	#	print(i)
	#	score100_i = score100[i,:,:]
	#	print(score100_i.shape)
	#	qt = QuantileTransformer(n_quantiles=1000, random_state=6, output_distribution='uniform', ignore_implicit_zeros=False, subsample=100000, copy=False)
	#	qt.fit(score100_i)
	#	score100_i_new = qt.fit_transform(score100_i)
	#	score100_new[i,:,:] = score100_i_new
	#	#score100_new[i,:,:] = all_score[:,0:48]
	#	output_name_i = str(i)+'check.txt'
	#	np.savetxt(output_name_i, score100_i_new[:,11], fmt='%10.2f', delimiter='\t')

	#score100 = score100_new

	### read samples & train & test ids 
	train_id = np.array(pd.read_csv(train_id_file, sep="\t", header=None))
	test_id = np.array(pd.read_csv(test_id_file, sep="\t", header=None))

	### get where are the train & test ids 
	train_id_pos = np.intersect1d(samples, train_id, return_indices=True)[1]
	test_id_pos = np.intersect1d(samples, test_id, return_indices=True)[1]

	### get average of train & test
	score100_train = np.mean(score100[train_id_pos,:,:], axis=0)#.reshape(1,score100.shape[1],score100.shape[2])
	score100_test = np.mean(score100[test_id_pos,:,:], axis=0)#.reshape(1,score100.shape[1],score100.shape[2])
	print('score100_train.shape')
	print(score100_train.shape)

	qt = QuantileTransformer(n_quantiles=1000, random_state=6, output_distribution='uniform', ignore_implicit_zeros=False, subsample=100000, copy=False)
	qt.fit(score100_train)
	score100_train_QT = np.copy(score100_train)
	for i in range(0,score100_train.shape[1]):
		print(i)
		score100_train_raw_i = score100_train[:,i]
		used_row_i = score100_train_raw_i!=0
		score100_train_QT_i = qt.fit_transform(score100_train_raw_i[used_row_i].reshape(np.sum(used_row_i),1))
		score100_train_QT[used_row_i,i] = score100_train_QT_i[:,0]
	score100_train_QT = score100_train_QT.reshape(1,score100.shape[1],score100.shape[2])
	print('score100_train_QT.shape')
	print(score100_train_QT.shape)

	#print('check 0 number')
	#print(score100_train.shape)
	#print(np.sum(score100_train[:,11]==0))
	#print(np.sum(score100_test[:,11]==0))

	qt = QuantileTransformer(n_quantiles=1000, random_state=6, output_distribution='uniform', ignore_implicit_zeros=False, subsample=100000, copy=False)
	qt.fit(score100_test)
	score100_test_QT = np.copy(score100_test)
	for i in range(0,score100_test.shape[1]):
		print(i)
		score100_test_raw_i = score100_test[:,i]
		used_row_i = score100_test_raw_i!=0
		score100_test_QT_i = qt.fit_transform(score100_test_raw_i[used_row_i].reshape(np.sum(used_row_i),1))
		score100_test_QT[used_row_i,i] = score100_test_QT_i[:,0]
	score100_test_QT = score100_test_QT.reshape(1,score100.shape[1],score100.shape[2])


	score100_train_test_ave = np.concatenate((score100_train_QT, score100_test_QT), axis=0)


	#print('check 0 number after QT')
	#print(score100_train_QT.shape)
	#print(np.sum(score100_train_QT[0,:,11]==0))
	#print(np.sum(score100_test_QT[0,:,11]==0))


	samples_new = np.array(['train', 'test']).astype('bytes_')

	### write output
	with h5py.File(output_h5file, "w") as outfile:
		outfile.create_dataset(chrom, data=np.array(score100_train_test_ave), compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
		outfile.create_dataset('%s_starts' % chrom, data=starts, compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
		outfile.create_dataset('feature_names', data=features, compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
		outfile.create_dataset('samples', data=samples_new, compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)

	outfile.close()



if __name__ == '__main__':
	fire.Fire(merge_ave)


