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



def pred_lgb_model(self, DNase_h5, motif_h5, cell_line, chrom, TFname, label_file_test, dir_out, trained_lgb_model):
	print('read DNase')
	hf = h5py.File(DNase_h5, 'r')
	dnase_a = np.array(hf.get(chrom))
	print(dnase_a.shape)
	dnase_samples = np.array([x.decode() for x in hf.get('samples')])
	hf.close()

	qt = QuantileTransformer(n_quantiles=1000, random_state=6, output_distribution='uniform', ignore_implicit_zeros=False, subsample=100000, copy=False)
	all_score = np.concatenate((dnase_a[0,:,:], dnase_a[1,:,:]), axis=1)
	qt.fit(all_score)
	new_data = qt.fit_transform(all_score)

	dnase_a_qt = np.copy(dnase_a)
	dnase_a_qt[0,:,:] = new_data[:,0:48]
	dnase_a_qt[1,:,:] = new_data[:,48:96]
	dnase_a = dnase_a_qt


	print('read motif score')
	hf = h5py.File(motif_h5, 'r')
	a = np.array(hf.get('scores'))
	hf.close()

	print('read label tsv files')
	d = np.array(pd.read_csv(label_file_test, sep="\t", header=0))
	df_temp = d[d[:,0] == chrom, :]
	del d

	print('get train data')
	test_data = np.concatenate((dnase_a[np.where(dnase_samples==str(cell_line))[0][0],:,:], a), axis=1)

	print('get weight & label for training')
	weight = (df_temp[:,3]!='A')*1
	labels = (df_temp[:,3]=='B')*1

	print('predict')
	with open(trained_lgb_model, 'rb') as pickle_file:
	    gbm = pickle.load(pickle_file)

	pickle_file.close()

	ypred = gbm.predict(test_data)

	print(average_precision_score(labels, ypred))

	print('write output')
	outputname = dir_out+TFname+'.'+str(cell_line)+'.'+chrom+'.lgb.pred.txt'
	ypred = ypred.reshape(ypred.shape[0],1)
	np.savetxt(outputname, ypred, fmt='%10.2f', delimiter='\t')


if __name__ == '__main__':
    fire.Fire(pred_lgb_model)



