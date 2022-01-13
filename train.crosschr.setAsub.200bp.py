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


def train_lgb_model(self, DNase_h5, DNase_h5_test, motif_h5, motif_h5_test, sigtype, cell_line, chrom, chrom_test, TFname, label_file_train_topN, label_file_train, label_file_test, dir_out, train_ct, test_ct, bins):
	print(DNase_h5)
	print(motif_h5)
	print(cell_line)
	print(chrom)
	print(TFname)
	print(label_file_train)
	print(dir_out)
	### functions & params
	def lgb_auprc_score(y_hat, data):
		y_true = data.get_label()
		# TODO try not to round yhat
		# y_hat = np.round(y_hat)  # scikits f1 doesn't like probabilities
		return 'auprc', average_precision_score(y_true, y_hat), True

	def focal_isoform_binary_object(pred, dtrain, alpha=0.5, beta=0.0, gamma=2.0):
		#	 alpha controls weight of positives
		#	 (0,1) less
		#	 >1 more or(0-0.5 less, 0.5-1 more)
		#	 beta controls the shift of loss function
		#	 >0 to left(less weight to well-trained samples)
		#	 gamma controls the steepness of loss function
		#	 >0
		label = dtrain.get_label()
		x = beta + (2.0 * label - 1) * gamma * pred
		p = 1. / (1. + np.exp(-x))
		# grad = (1 + (alpha - 1) * label) * (2 * label - 1) * (p - 1)
		grad = (1 - label + (label * 2 - 1) * alpha) * (2 * label - 1) * (p - 1)
		# hess = (1 + (alpha - 1) * label) * gamma * (1 - p) * p
		hess = (1 - label + (label * 2 - 1) * alpha) * gamma * (1 - p) * p
		return grad, hess

	params = {
		'boosting_type': 'gbdt',
		'objective': 'binary',
		'metric': ["auc"],
		'metric_freq': 10,
		'num_leaves': 63,
		'num_threads': 2,
		'learning_rate': 0.05,
		'feature_fraction': 1,
		'boost_from_average': False,
		'verbose': 1
	}

	print('read DNase')
	DNase_h5a = DNase_h5+chrom+'_all_cell_types.'+sigtype+'.ave'+'.h5'
	hf = h5py.File(DNase_h5a, 'r')
	dnase_a = np.array(hf.get(chrom))[:,:,np.concatenate((np.arange(8,16), np.arange(32,40)), axis=0)]
	print(dnase_a.shape)
	dnase_a = np.sum(dnase_a, axis=2)
	print(dnase_a.shape)
	dnase_a = dnase_a.reshape((dnase_a.shape[0], dnase_a.shape[1], 1))
	print(dnase_a.shape)
	print(np.max(dnase_a))
	dnase_samples = np.array([x.decode() for x in hf.get('samples')])
	hf.close()


	print('read motif score')
	motif_h5a = motif_h5+chrom+'_motif_features_lightGBM.h5'
	hf = h5py.File(motif_h5a, 'r')
	a = np.array(hf.get('scores'))
	print(np.max(a))
	scale_motif = np.max(a)
	a = a/scale_motif
	hf.close()


	print('read label tsv files')
	d = np.array(pd.read_csv(label_file_train, sep="\t", header=0))
	df_temp = d[d[:,0] == chrom, :]
	d1 = np.array(pd.read_csv(label_file_train_topN, sep="\t", header=0))
	df_temp1 = d1[d1[:,0] == chrom, :]
	
	print('for each traininig chr')

	train_chr = ['chr1', 'chr3']#, 'chr15', 'chr16', 'chr17', 'chr20']
	#['chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr9', 'chr10', 'chr11']
	print(dnase_a.shape)
	print(a.shape)
	print(df_temp.shape)
	print(df_temp1.shape)
	for chrom in train_chr:
		#, 'chr7', 'chr9', 'chr20', 'chr22']:
		DNase_h5a = DNase_h5+chrom+'_all_cell_types.'+sigtype+'.ave'+'.h5'
		hf = h5py.File(DNase_h5a, 'r')
		dnase_a_i = np.array(hf.get(chrom))
		dnase_a_i = np.sum(dnase_a_i, axis=2)
		dnase_a_i = dnase_a_i.reshape((dnase_a_i.shape[0], dnase_a_i.shape[1], 1))
		dnase_a = np.concatenate((dnase_a, dnase_a_i), axis=1)
		print(dnase_a.shape)
		print(np.max(dnase_a))
		hf.close()
		
	for chrom in train_chr:
		print('read motif score')
		motif_h5a = motif_h5+chrom+'_motif_features_lightGBM.h5'
		hf = h5py.File(motif_h5a, 'r')
		a_i = np.array(hf.get('scores'))
		a_i = a_i/scale_motif
		a = np.concatenate((a, a_i), axis=0)
		print(a.shape)
		hf.close()

	for chrom in train_chr:
		df_temp_i = d[d[:,0] == chrom, :]
		df_temp1_i = d1[d1[:,0] == chrom, :]
		df_temp = np.concatenate((df_temp, df_temp_i), axis=0)
		df_temp1 = np.concatenate((df_temp1, df_temp1_i), axis=0)
		print(df_temp.shape)
		print(df_temp1.shape)


	del d
	del d1

	print('get train data')
	used_rows = df_temp[:,3]==df_temp1[:,3]
	df_temp1_new = df_temp1[used_rows,:]

	del df_temp
	del df_temp1
	#print(dnase_samples)
	#print(str(cell_line))
	train_data = np.concatenate((dnase_a[np.where(dnase_samples=='train')[0][0],:,:], a), axis=1)[used_rows,:]
	#train_data = a[used_rows,:]#np.concatenate((a, a), axis=1)[used_rows,:]
	del dnase_a
	del a

	print('get weight & label for training')
	weight = (df_temp1_new[:,3]!='A')*1
	labels = (df_temp1_new[:,3]=='B')*1
	del df_temp1_new
	print('training peak number:')
	print(np.sum(labels))

	print('random split')
	np.random.seed(2019)
	random_indices = np.random.choice(train_data.shape[0], size=int(train_data.shape[0]*0.2), replace=False)
	train_data_1 = np.delete(train_data, random_indices, 0)
	train_data_2 = train_data[random_indices,:]
	del train_data
	labels_1 = np.delete(labels, random_indices)
	labels_2 = labels[random_indices]
	del labels
	weight_1 = np.delete(weight, random_indices)
	weight_2 = weight[random_indices]
	del weight

	print('initialize lgb dataset')
	train_data_lgb = lgb.Dataset(train_data_1, label=labels_1, weight=weight_1)
	test_data_lgb = lgb.Dataset(train_data_2, label=labels_2, weight=weight_2)

	del labels_1
	del labels_2
	del weight_1
	del weight_2

	print('lgb beta')
	beta = - np.log10(2 * train_data_1.shape[0]/np.where(train_data_lgb.get_label() > 0)[0].shape[0] - 1)

	del train_data_1
	del train_data_2

	print('train model')
	evals_result = {}
	gbm = lgb.train(params=params,
		train_set=train_data_lgb,
		fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=0.5, beta=beta, gamma=1),
		feval=lgb_auprc_score,
		valid_sets=[test_data_lgb], 
		early_stopping_rounds=50,
		evals_result=evals_result,
		num_boost_round=200,
		keep_training_booster=False)

	print('save model')
	#with open(dir_out+TFname+'.'+str(cell_line)+'.'+chrom+'.model.pkl', 'wb') as fout:
	#	pickle.dump(gbm, fout)

	del train_data_lgb
	del test_data_lgb

	#label_file_test = 'peaks/5fc/FOXA1.test.5fc.label.tsv'
	#label_file_test = 'peaks/5fc/CTCF.test.5fc.label.tsv'
	print('get test data mat')

	# read bin file
	bin_mat = pd.read_csv(bins, header=None, sep='\t').values

	### get all chr
	chr_all = np.unique(bin_mat[:,0])

	print('test motif')
	#for chrom_test in ['chr1', 'chr12', 'chr4', 'chr6']:
	#for i in range(1,23):
	for chrom_test in chr_all:
		#chrom_test = 'chr'+str(i)
		motif_h5_test1 = motif_h5_test+chrom_test+'_motif_features_lightGBM.h5'
		hf = h5py.File(motif_h5_test1, 'r')
		at = np.array(hf.get('scores'))
		at = at/scale_motif
		hf.close()

		DNase_h5_test1 = DNase_h5_test+chrom_test+'_all_cell_types.'+sigtype+'.ave'+'.h5'
		hf = h5py.File(DNase_h5_test1, 'r')
		dnase_at = np.array(hf.get(chrom_test))
		dnase_at = np.sum(dnase_at, axis=2)
		dnase_at = dnase_at.reshape((dnase_at.shape[0], dnase_at.shape[1], 1))
		hf.close()

		test_data = np.concatenate((dnase_at[np.where(dnase_samples=='test')[0][0],:,:], at), axis=1)
		#test_data = at#np.concatenate((dnase_at[np.where(dnase_samples=='test')[0][0],:,:], at), axis=1)
		ypred = gbm.predict(test_data)
		outputname = dir_out+TFname+'.'+train_ct+'.'+test_ct+'.'+chrom_test+'.lgb.pred.txt'
		ypred = ypred.reshape(ypred.shape[0],1)
		np.savetxt(outputname, ypred, fmt='%10.2f', delimiter='\t')
		
		#d = np.array(pd.read_csv(label_file_test, sep="\t", header=0))
		#df_temp = d[d[:,0] == chrom_test, :]
		#labels = (df_temp[:,3]=='B')*1
		#del df_temp
		#del d
		del at
		del dnase_at
		del test_data
		print('final')
		#print(average_precision_score(labels, ypred))

		print('write output')


if __name__ == '__main__':
    fire.Fire(train_lgb_model)


