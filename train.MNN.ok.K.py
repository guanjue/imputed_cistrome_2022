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
from sklearn.cluster import KMeans

from sklearn.metrics import average_precision_score
from sklearn.preprocessing import StandardScaler

def train_lgb_model(self, DNase_h5, DNase_h5_test, motif_h5, motif_h5_test, cell_line, chrom, chrom_test, TFname, label_file_train, label_file_test, km_n, dir_out):
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

	print('read DNase')
	hf_test = h5py.File(DNase_h5_test, 'r')
	dnase_a_test = np.array(hf_test.get(chrom_test))
	print(np.max(dnase_a_test))
	dnase_samples_test = np.array([x.decode() for x in hf_test.get('samples')])
	hf_test.close()

	print('read motif score')
	hf = h5py.File(motif_h5, 'r')
	a = np.array(hf.get('scores'))
	print(np.max(a))
	a = a/np.max(a)
	hf.close()

	print('read motif score')
	hf = h5py.File(motif_h5_test, 'r')
	a_test = np.array(hf.get('scores'))
	a_test = a_test/np.max(a_test)
	hf.close()

	print('read label tsv files')
	d = np.array(pd.read_csv(label_file_train, sep="\t", header=0))
	df_temp = d[d[:,0] == chrom, :]
	#labels00 = (df_temp[:,3]=='B')*1
	del d

	print('get train data')
	train_data = np.concatenate((dnase_a[np.where(dnase_samples=='train')[0][0],:,:], a), axis=1)
	test_data = np.concatenate((dnase_a_test[np.where(dnase_samples_test=='test')[0][0],:,:], a_test), axis=1)
	ypred = test_data[:,0]
	np.random.seed(2019)
	#alldat = np.concatenate((train_data[:,np.concatenate((np.arange(0,24), np.arange(24,train_data.shape[1])), axis=0)], test_data[:,np.concatenate((np.arange(0,24), np.arange(24,test_data.shape[1])), axis=0)]), axis=0)
	alldat = np.concatenate((train_data, test_data), axis=0)
	scaler = StandardScaler()
	alldat = scaler.fit_transform(alldat)
	#alldat = np.concatenate((train_data[:,np.concatenate((np.arange(0,24), np.arange(24,48)), axis=0)], test_data[:,np.concatenate((np.arange(0,24), np.arange(24,48)), axis=0)]), axis=0)
	kmeans = KMeans(n_clusters=km_n, random_state=0).fit(alldat)
	#kmeans = KMeans(n_clusters=km_n, random_state=0).fit(train_data[:,np.concatenate((np.arange(0,48), np.arange(160,199)), axis=0)])
	#kmeans_c = kmeans.cluster_centers_
	print('check test')

	def get_dist(c,a):
		da = []
		for i in range(0,a.shape[0]):
			if i%10000:
				print(i)
			dj = []
			for j in range(0,c.shape[0]):
				dj.append(np.linalg.norm(a[i]-c[j]))
			da.append(dj)
		da = np.array(da)
		return(da)

	
	#dist_c = get_dist(kmeans_c, test_data[:,np.concatenate((np.arange(0,48), np.arange(160,199)), axis=0)])
	#km_test = np.argmin(dist_c, axis=1)

	d = np.array(pd.read_csv(label_file_test, sep="\t", header=0))
	df_temp = d[d[:,0] == chrom_test, :]
	labels0 = (df_temp[:,3]=='B')*1
	print('kmeanspred.shape')
	print(test_data.shape)
	#km_train = kmeans.labels_
	km_train = kmeans.labels_[0:train_data.shape[0]]
	km_test = kmeans.labels_[train_data.shape[0]:(train_data.shape[0]+test_data.shape[0])]

	for kkk in range(0,km_n):
		used_id0 = km_train == kkk
		print('get weight & label for training')
		train_data1 = train_data[used_id0,:]
		weight = (df_temp[used_id0,3]!='A')*1
		labels = (df_temp[used_id0,3]=='B')*1

		print(train_data1.shape)
		print(weight.shape)
		print(labels.shape)

		print('random split')
		np.random.seed(2019)
		random_indices = np.random.choice(train_data1.shape[0], size=int(train_data1.shape[0]*0.2), replace=False)
		train_data_1 = np.delete(train_data1, random_indices, 0)
		train_data_2 = train_data1[random_indices,:]
		labels_1 = np.delete(labels, random_indices)
		labels_2 = labels[random_indices]
		weight_1 = np.delete(weight, random_indices)
		weight_2 = weight[random_indices]

		print('initialize lgb dataset')
		train_data_lgb = lgb.Dataset(train_data_1, label=labels_1, weight=weight_1)
		test_data_lgb = lgb.Dataset(train_data_2, label=labels_2, weight=weight_2)

		print('lgb beta')
		if np.where(train_data_lgb.get_label() > 0)[0].shape[0]!=0:
			beta = - np.log10(2 * train_data1.shape[0]/np.where(train_data_lgb.get_label() > 0)[0].shape[0] - 1)

			print('train model')
			evals_result = {}
			gbm = lgb.train(params=params,
				train_set=train_data_lgb,
				fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=0.5, beta=beta, gamma=1),
				feval=lgb_auprc_score,
				valid_sets=[test_data_lgb], 
				early_stopping_rounds=100,
				evals_result=evals_result,
				num_boost_round=200,
				keep_training_booster=False)

			#print('save model')
			#with open(dir_out+TFname+'.'+str(cell_line)+'.'+chrom+'.model.pkl', 'wb') as fout:
			#	pickle.dump(gbm, fout)


	
			#label_file_test = 'peaks/5fc/FOXA1.test.5fc.label.tsv'
			#label_file_test = 'peaks/5fc/CTCF.test.5fc.label.tsv'
			#test_data = np.concatenate((dnase_a[np.where(dnase_samples=='test')[0][0],:,:], a), axis=1)
			#d = np.array(pd.read_csv(label_file_test, sep="\t", header=0))
			#df_temp = d[d[:,0] == chrom, :]
			print(ypred.shape)
			print(gbm.predict(test_data).shape)
			used_id0_test = km_test == kkk
			if np.sum(used_id0_test)>0:
				ypred[used_id0_test] = gbm.predict(test_data)[used_id0_test]
		else:
			used_id0_test = km_test == kkk
			ypred[used_id0_test]=(-10)

		print('final')
		print(kkk)
		print(average_precision_score(labels0, ypred))
	#print(average_precision_score(labels00, ypred))
	print(average_precision_score(labels0, ypred))

	#print('write output')
	
	outputname = dir_out+TFname+'.'+str(cell_line)+'.'+chrom+'.lgb.pred.txt'
	ypred = ypred.reshape(ypred.shape[0],1)
	np.savetxt(outputname, ypred, fmt='%10.2f', delimiter='\t')


if __name__ == '__main__':
    fire.Fire(train_lgb_model)


