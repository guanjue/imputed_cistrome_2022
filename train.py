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


def train_lgb_model(self, DNase_h5, motif_h5, cell_line, chrom, TFname, label_file_train, label_file_test, dir_out):
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
	hf = h5py.File(DNase_h5, 'r')
	dnase_a = np.array(hf.get(chrom))
	print(np.max(dnase_a))
	dnase_samples = np.array([x.decode() for x in hf.get('samples')])
	hf.close()
	print('read motif score')
	hf = h5py.File(motif_h5, 'r')
	a = np.array(hf.get('scores'))
	print(np.max(a))
	a = a/np.max(a)
	hf.close()

	print('read label tsv files')
	d = np.array(pd.read_csv(label_file_train, sep="\t", header=0))
	df_temp = d[d[:,0] == chrom, :]
	del d

	print('get train data')
	train_data = np.concatenate((dnase_a[np.where(dnase_samples==str(cell_line))[0][0],:,:], a), axis=1)

	print('get weight & label for training')
	weight = (df_temp[:,3]!='A')*1
	labels = (df_temp[:,3]=='B')*1

	print('random split')
	np.random.seed(2019)
	random_indices = np.random.choice(train_data.shape[0], size=int(train_data.shape[0]*0.2), replace=False)
	train_data_1 = np.delete(train_data, random_indices, 0)
	train_data_2 = train_data[random_indices,:]
	labels_1 = np.delete(labels, random_indices)
	labels_2 = labels[random_indices]
	weight_1 = np.delete(weight, random_indices)
	weight_2 = weight[random_indices]

	print('initialize lgb dataset')
	train_data_lgb = lgb.Dataset(train_data_1, label=labels_1, weight=weight_1)
	test_data_lgb = lgb.Dataset(train_data_2, label=labels_2, weight=weight_2)

	print('lgb beta')
	beta = - np.log10(2 * train_data_1.shape[0]/np.where(train_data_lgb.get_label() > 0)[0].shape[0] - 1)

	print('train model')
	evals_result = {}
	gbm = lgb.train(params=params,
		train_set=train_data_lgb,
		fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=0.5, beta=beta, gamma=1),
		feval=lgb_auprc_score,
		valid_sets=[test_data_lgb], 
		early_stopping_rounds=50,
		evals_result=evals_result,
		num_boost_round=20,
		keep_training_booster=False)

	print('save model')
	#with open(dir_out+TFname+'.'+str(cell_line)+'.'+chrom+'.model.pkl', 'wb') as fout:
	#	pickle.dump(gbm, fout)


	
	#label_file_test = 'peaks/5fc/FOXA1.test.5fc.label.tsv'
	#label_file_test = 'peaks/5fc/CTCF.test.5fc.label.tsv'
	test_data = np.concatenate((dnase_a[np.where(dnase_samples=='test')[0][0],:,:], a), axis=1)
	d = np.array(pd.read_csv(label_file_test, sep="\t", header=0))
	print(d.shape)
	df_temp = d[d[:,0] == chrom, :]
	print(np.sum(d[:,0] == chrom))
	print(df_temp.shape)
	ypred = gbm.predict(test_data)
	labels = (df_temp[:,3]=='B')*1
	print(ypred.shape)
	print(labels.shape)
	print('final')
	print(average_precision_score(labels, ypred))

	print('write output')
	outputname = dir_out+TFname+'.'+str(cell_line)+'.'+chrom+'.lgb.pred.txt'
	ypred = ypred.reshape(ypred.shape[0],1)
	np.savetxt(outputname, ypred, fmt='%10.2f', delimiter='\t')


if __name__ == '__main__':
    fire.Fire(train_lgb_model)


