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

import scipy
import statsmodels.stats.multitest as multi
from patsy import dmatrix
import statsmodels.api as sm

import random
import numpy as np
from scipy.spatial.distance import cdist  # $scipy/spatial/distance.py
	# http://docs.scipy.org/doc/scipy/reference/spatial.html
from scipy.sparse import issparse  # $scipy/sparse/csr.py


#chrom_j='chr19'
#time python scripts/SQTKNNN.py SQTKNNN --chrom=$chrom_j --CP_id_file=$chrom_j'.50bp.CPid.bed' \
#--input_h5file='DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_'$chrom_j'_all_cell_types.h5' \
#--output_h5file='DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_'$chrom_j'_all_cell_types.SQTKNNN.h5' \
#--DNase_ave_f_DR_file='/homes1/gxiang/projects/impute_cistrome/get_KNNN_sig/DNase.f.'$chrom_j'.ave.h5' \
#--DNase_ave_r_DR_file='/homes1/gxiang/projects/impute_cistrome/get_KNNN_sig/DNase.r.'$chrom_j'.ave.h5' \
#--DNase_ave_f_file='/homes1/gxiang/projects/impute_cistrome/get_KNNN_sig/DNase.f.'$chrom_j'.ave.spline.h5' \
#--DNase_ave_r_file='/homes1/gxiang/projects/impute_cistrome/get_KNNN_sig/DNase.r.'$chrom_j'.ave.spline.h5'


#...............................................................................
def kmeans( X, centres, delta=.001, maxiter=100, metric="euclidean", p=2, verbose=1 ):
	""" centres, Xtocentre, distances = kmeans( X, initial centres ... )
	in:
		X N x dim  may be sparse
		centres k x dim: initial centres, e.g. random.sample( X, k )
		delta: relative error, iterate until the average distance to centres
			is within delta of the previous average distance
		maxiter
		metric: any of the 20-odd in scipy.spatial.distance
			"chebyshev" = max, "cityblock" = L1, "minkowski" with p=
			or a function( Xvec, centrevec ), e.g. Lqmetric below
		p: for minkowski metric -- local mod cdist for 0 < p < 1 too
		verbose: 0 silent, 2 prints running distances
	out:
		centres, k x dim
		Xtocentre: each X -> its nearest centre, ints N -> k
		distances, N
	see also: kmeanssample below, class Kmeans below.
	"""
	if not issparse(X):
		X = np.asanyarray(X)  # ?
	centres = centres.todense() if issparse(centres) \
		else centres.copy()
	N, dim = X.shape
	k, cdim = centres.shape
	if dim != cdim:
		raise ValueError( "kmeans: X %s and centres %s must have the same number of columns" % (X.shape, centres.shape ))
	#if verbose:
	#	print("kmeans: X %s  centres %s  delta=%.2g  maxiter=%d  metric=%s" % (X.shape, centres.shape, delta, maxiter, metric))
	allx = np.arange(N)
	prevdist = 0
	for jiter in range( 1, maxiter+1 ):
		D = cdist_sparse( X, centres, metric=metric, p=p )  # |X| x |centres|
		xtoc = D.argmin(axis=1)  # X -> nearest centre
		distances = D[allx,xtoc]
		avdist = distances.mean()  # median ?
		#if verbose >= 2:
			#print("kmeans: av |X - nearest centre| = %.4g" % avdist)
		if (1 - delta) * prevdist <= avdist <= prevdist or jiter == maxiter:
			break
		prevdist = avdist
		for jc in range(k):  # (1 pass in C)
			c = np.where( xtoc == jc )[0]
			if len(c) > 0:
				centres[jc] = X[c].mean( axis=0 )
	#if verbose:
		#print("kmeans: %d iterations  cluster sizes:" % jiter, np.bincount(xtoc))
	if verbose >= 2:
		r50 = np.zeros(k)
		r90 = np.zeros(k)
		for j in range(k):
			dist = distances[ xtoc == j ]
			if len(dist) > 0:
				r50[j], r90[j] = np.percentile( dist, (50, 90) )
		#print("kmeans: cluster 50 % radius", r50.astype(int))
		#print("kmeans: cluster 90 % radius", r90.astype(int))
			# scale L1 / dim, L2 / sqrt(dim) ?
	return centres, xtoc, distances

#...............................................................................
def kmeanssample( X, k, nsample=0, **kwargs ):
	""" 2-pass kmeans, fast for large N:
		1) kmeans a random sample of nsample ~ sqrt(N) from X
		2) full kmeans, starting from those centres
	"""
		# merge w kmeans ? mttiw
		# v large N: sample N^1/2, N^1/2 of that
		# seed like sklearn ?
	N, dim = X.shape
	if nsample == 0:
		nsample = max( 2*np.sqrt(N), 10*k )
	Xsample = randomsample( X, int(nsample) )
	pass1centres = randomsample( X, int(k) )
	samplecentres = kmeans( Xsample, pass1centres, **kwargs )[0]
	return kmeans( X, samplecentres, **kwargs )

def cdist_sparse( X, Y, **kwargs ):
	""" -> |X| x |Y| cdist array, any cdist metric
		X or Y may be sparse -- best csr
	"""
		# todense row at a time, v slow if both v sparse
	sxy = 2*issparse(X) + issparse(Y)
	if sxy == 0:
		return cdist( X, Y, **kwargs )
	d = np.empty( (X.shape[0], Y.shape[0]), np.float64 )
	if sxy == 2:
		for j, x in enumerate(X):
			d[j] = cdist( x.todense(), Y, **kwargs ) [0]
	elif sxy == 1:
		for k, y in enumerate(Y):
			d[:,k] = cdist( X, y.todense(), **kwargs ) [0]
	else:
		for j, x in enumerate(X):
			for k, y in enumerate(Y):
				d[j,k] = cdist( x.todense(), y.todense(), **kwargs ) [0]
	return d

def randomsample( X, n ):
	""" random.sample of the rows of X
		X may be sparse -- best csr
	"""
	sampleix = random.sample( range( X.shape[0] ), int(n) )
	return X[sampleix]

def nearestcentres( X, centres, metric="euclidean", p=2 ):
	""" each X -> nearest centre, any metric
			euclidean2 (~ withinss) is more sensitive to outliers,
			cityblock (manhattan, L1) less sensitive
	"""
	D = cdist( X, centres, metric=metric, p=p )  # |X| x |centres|
	return D.argmin(axis=1)

def Lqmetric( x, y=None, q=.5 ):
	# yes a metric, may increase weight of near matches; see ...
	return (np.abs(x - y) ** q) .mean() if y is not None \
		else (np.abs(x) ** q) .mean()


def cluster_DR(CP_id_new, score_mat25_ave_spline, q_thresh, r_thresh, Kmeans_n):
	score_mat25_ave_spline_DR = np.copy(score_mat25_ave_spline[:,0:Kmeans_n])
	for i in np.unique(CP_id_new):
		print(i)
		### get score CPi
		used_row_i = CP_id_new==i
		score_mat25_ave_spline_i = score_mat25_ave_spline[used_row_i,:]
		### kmeans cluster based on r
		randomcentres = randomsample( np.transpose(score_mat25_ave_spline_i), Kmeans_n )
		centres, idx, dist = kmeans( np.transpose(score_mat25_ave_spline_i), randomcentres, delta=0.001, maxiter=10, metric='correlation', verbose=2 )
		### get unique k
		idx_u = np.unique(idx)
		### get DR mat
		for j in idx_u:
			used_col_j = idx==j
			if np.sum(used_col_j)>1:
				score_mat25_ave_spline_DR[used_row_i,j] = np.mean(score_mat25_ave_spline_i[:,idx==j], axis=1)
			else:
				score_mat25_ave_spline_DR[used_row_i,j] = score_mat25_ave_spline_i[:,idx==j][:,0]
	return(score_mat25_ave_spline_DR)


def corr2d_coeff(A, B):
	# Rowwise mean of input arrays & subtract from input arrays themeselves
	A_mA = A - A.mean(0)[None, :]
	B_mB = B - B.mean(0)[None, :]
	ssA = (A_mA**2).sum(0)
	ssB = (B_mB**2).sum(0)
	return(np.round(np.dot(A_mA.T, B_mB) / np.sqrt(np.dot(ssA[:, None],ssB[None])), 3))

def fdr_adj_binary(x, fdr_thresh):
	return(multi.fdrcorrection( scipy.stats.norm.sf( (x-np.mean(x))/np.std(x) ), alpha=fdr_thresh)[0])

def spline_cpkcbg(ref, tar, cpkbg, sn):
	ref_log2 = np.log2(ref+sn)
	tar_log2 = np.log2(tar+sn)
	x = ref_log2[cpkbg]/2 + tar_log2[cpkbg]/2
	y = ref_log2[cpkbg] - tar_log2[cpkbg]
	### get all data
	xall = ref_log2/2 + tar_log2/2
	### fit model
	x_cubic = dmatrix('bs(x, knots=np.arange(np.min(x), np.max(x)+0.001, step = np.max(x)/5))', {'x': x})
	fit_cubic = sm.GLM(y, x_cubic).fit()
	### predict in all data
	line_cubic = fit_cubic.predict(dmatrix('bs(xall, knots=np.arange(np.min(x), np.max(x)+0.001, step = np.max(x)/5))', {'xp': xall}))
	line_cubic[np.isnan(line_cubic)] = 0
	tar_spline_norm_log2 = tar_log2 + line_cubic
	tar_spline_norm_log2[tar_log2==np.log2(sn)] = np.log2(sn)
	tar_spline_norm = 2**tar_spline_norm_log2-sn
	tar_spline_norm[tar_spline_norm<0] = 0
	return(tar_spline_norm)

def get_cpk_cbg_array(score_mat100_ave, fdr_thresh, q_thresh):
	### get fdr binary mat
	score_mat100_ave_i_zpfdr_binary_mat = score_mat100_ave[:,0].reshape(score_mat100_ave.shape[0], 1)
	for i in range(0,score_mat100_ave.shape[1]):
		score_mat100_ave_i = score_mat100_ave[:,i]
		score_mat100_ave_i_zpfdr_binary = fdr_adj_binary(score_mat100_ave_i, fdr_thresh)
		score_mat100_ave_i_zpfdr_binary_mat = np.concatenate((score_mat100_ave_i_zpfdr_binary_mat, score_mat100_ave_i_zpfdr_binary.reshape(score_mat100_ave_i_zpfdr_binary.shape[0], 1)), axis=1)
	### get cpk cbg vector
	score_mat100_ave_i_zpfdr_binary_mat = score_mat100_ave_i_zpfdr_binary_mat[:,1:]
	cpk = np.sum(score_mat100_ave_i_zpfdr_binary_mat, axis=1) >= (score_mat100_ave_i_zpfdr_binary_mat.shape[1]*q_thresh)
	cbg = np.sum(score_mat100_ave_i_zpfdr_binary_mat, axis=1) == 0
	return(cpk, cbg)

def QTnorm(tar, ref):
	return(np.sort(ref)[(tar.argsort()).argsort()])

### KNN norm core
def KNN_norm(CP_id_new, score_mat25, score_mat25_SQT, score_mat25_ave_spline, q_thresh, r_thresh):
	score_mat25_SQT_KNN = np.copy(score_mat25_SQT)
	KNN_id = score_mat25_SQT[:,0].reshape(score_mat25_SQT.shape[0],1)
	for i in np.unique(CP_id_new):
		#print(i)
		### get score CPi
		used_row_i = CP_id_new==i
		score_mat25_SQT_i = score_mat25_SQT[used_row_i,:]
		score_mat25_ave_spline_i = score_mat25_ave_spline[used_row_i,:]
		### get r
		sample2ave_r_i = corr2d_coeff(score_mat25_SQT_i, score_mat25_ave_spline_i)
		KNN_id_i = np.argmax(sample2ave_r_i, axis=1)
		KNN_id[used_row_i,0] = KNN_id_i
		### remove na and negative rs
		sample2ave_r_i[np.isnan(sample2ave_r_i)] = 0
		sample2ave_r_i[sample2ave_r_i<0] = 0
		###
		### check if sample2ave_r_i is NA (due to all 0 bins)
		#if not np.isnan(np.max(sample2ave_r_i)):
		#if not 100==1:
		### get top match
		sample_knn_ave = []
		for j in range(0,sample2ave_r_i.shape[0]):
			sample2ave_r_ij = sample2ave_r_i[j,:]
			### get knn ave
			used_ave = (sample2ave_r_ij>np.quantile(sample2ave_r_ij, q_thresh)) & (sample2ave_r_ij>r_thresh)
			if np.sum(used_ave)==0:
				used_ave = (sample2ave_r_ij>np.quantile(sample2ave_r_ij, q_thresh))
			if np.sum(used_ave)>0:
				sample2ave_r_ij_knn = sample2ave_r_ij[used_ave]
				score_mat25_ave_spline_i_knn = score_mat25_ave_spline_i[:,used_ave]
				score_mat25_SQT_i[:,j] = score_mat25_SQT_i[:,j] * np.mean(1-sample2ave_r_ij_knn) + np.dot(score_mat25_ave_spline_i_knn, sample2ave_r_ij_knn) / np.sum(used_ave)
		### replace OD mat
		score_mat25_SQT_KNN[used_row_i,:] = score_mat25_SQT_i
	### keep the 0s
	score_mat25_SQT_KNN[score_mat25==0] = 0
	score_mat25_SQT_KNN[np.isnan(score_mat25_SQT_KNN)] = score_mat25[np.isnan(score_mat25_SQT_KNN)]
	score_mat25_SQT_KNN = np.round(score_mat25_SQT_KNN, 2)
	return(score_mat25_SQT_KNN, KNN_id)

def get_h5mat(c12,c13, bin_size):
	row_n = c12.shape[0]
	### get center 25bp sig
	c12 = np.round(c12/bin_size, 3).reshape(row_n,1)
	c13 = np.round(c13/bin_size, 3).reshape(row_n,1)
	c14 = np.roll(c12, -1).reshape(row_n,1)
	c15 = np.roll(c13, -1).reshape(row_n,1)
	### get center 50bp & 100bp
	cx50 = np.round((c12+c13)/2, 3)
	cx100 = np.round((c12+c13+c14+c15)/4, 3)
	### 100u
	c0 = np.roll(cx100, 14).reshape(row_n,1)
	c1 = np.roll(cx100, 12).reshape(row_n,1)
	c2 = np.roll(cx100, 10).reshape(row_n,1)
	c3 = np.roll(cx100, 8).reshape(row_n,1)
	### 50u
	c4 = np.roll(cx50, 6).reshape(row_n,1)
	c5 = np.roll(cx50, 5).reshape(row_n,1)
	c6 = np.roll(cx50, 4).reshape(row_n,1)
	c7 = np.roll(cx50, 3).reshape(row_n,1)
	### 25u
	c8 = np.roll(c12, 2).reshape(row_n,1)
	c9 = np.roll(c13, 2).reshape(row_n,1)
	c10 = np.roll(c12, 1).reshape(row_n,1)
	c11 = np.roll(c13, 1).reshape(row_n,1)
	### 50d
	c16 = np.roll(cx50, -2).reshape(row_n,1)
	c17 = np.roll(cx50, -3).reshape(row_n,1)
	c18 = np.roll(cx50, -4).reshape(row_n,1)
	c19 = np.roll(cx50, -5).reshape(row_n,1)
	### 100d
	c20 = np.roll(cx100, -6).reshape(row_n,1)
	c21 = np.roll(cx100, -8).reshape(row_n,1)
	c22 = np.roll(cx100, -10).reshape(row_n,1)
	c23 = np.roll(cx100, -12).reshape(row_n,1)
	### h5mat
	h5mat = np.concatenate((c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23), axis=1)
	del(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20,c21,c22,c23)
	return(h5mat)



class SQTKNNN_norm(object):

	def SQTKNNN(self, input_h5file, output_h5file, DNase_ave_f_DR_file, DNase_ave_f_file, DNase_ave_r_DR_file, DNase_ave_r_file, CP_id_file, chrom, fdr_thresh = 0.1, q_thresh = 0.95, r_thresh = 0.6, bin_size = 25, knn_n = 1, sn = 1):
		### get info mat
		#CP_id_file = chrom+'.50bp.CPid.bed'
		#input_h5file = "../cistrome_impute_results_human8/hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_all_cell_types.h5" % (chrom)
		#output_h5file = "../cistrome_impute_results_human8/hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s_all_cell_types.SQTKNN.h5" % (chrom)
		#DNase_ave_f_DR_file = "DNase.f.%s.ave.h5" % (chrom)
		#DNase_ave_f_file = "DNase.f.%s.ave.spline.h5" % (chrom)
		#DNase_ave_r_DR_file = "DNase.r.%s.ave.h5" % (chrom)
		#DNase_ave_r_file = "DNase.r.%s.ave.spline.h5" % (chrom)
		center_25bp_all = [12,13,36,37]

		### read CP id
		CP_id = pd.read_csv(CP_id_file, sep="\t", header=None)[3]
		CP_id[CP_id=='.'] = '0'
		id_old = '0'
		k = 0
		CP_id_new = []
		for CP_id_i in CP_id:
			CP_id_i = [x.strip() for x in CP_id_i.split(',')][0]
			id_new = CP_id_i
			if id_new!=id_old:
				id_old = id_new
				k = k+1
			CP_id_new.append(k)

		###
		CP_id_new = np.array(CP_id_new)


		### read ave Dimension reduced mat
		with h5py.File(DNase_ave_f_DR_file, "r") as infile:
			score_mat25_ave_spline_DR = infile["ave_score"][...]

		infile.close()

		### read ave mat
		with h5py.File(DNase_ave_f_file, "r") as infile:
			score_mat25_ave_spline = infile["ave_score"][...]

		infile.close()

		### read new signal mat
		with h5py.File(input_h5file, "r") as infile:
			score_new = infile[chrom][...]
			features = infile['feature_names'][...]
			starts = infile['%s_starts' % chrom][...]
			samples = infile['samples'][...]

		infile.close()

		### find the highest correlated ave for each sample
		score_new_SQT_KNN = np.copy(score_new)
		KNN_id12 = np.copy(score_new[0,:,0:8])
		### get quantile q based on KNN_n
		q_thresh = 1-knn_n/score_mat25_ave_spline_DR.shape[1]


		for i in range(0, score_new.shape[0]):
			print(i)
			score_mat25a_KNN = np.copy(score_new[0,:,0:2])
			### reverse use column 0 & 1
			for j in range(0,2):
				print('bin set a')
				print('find the closest average sample')
				center_25bp_cola = center_25bp_all[j]
				score_mat25a_test = score_new[i,:,center_25bp_cola].reshape(score_new.shape[1], 1)*25
				sample2ave_r = corr2d_coeff(score_mat25a_test+sn, score_mat25_ave_spline+sn)
				sample2ave_r_max_id = np.argmax(sample2ave_r, axis=1)
				print('initialize QTnorm sample to Spline ave')
				score_mat25_SQTa = np.copy(score_mat25a_test)
				print('get QTnorm sample to Spline ave')
				score_mat25_SQTa = QTnorm(score_mat25_SQTa[:,0], score_mat25_ave_spline[:,sample2ave_r_max_id[0]])
				###### KNN norm
				print('KNN norm for each CP')
				score_mat25_SQT_KNNa = np.copy(score_mat25_SQTa)
				print('get KNN norm matrix')
				score_mat25_SQT_KNNa, KNN_ida = KNN_norm(CP_id_new, score_mat25a_test, score_mat25_SQTa.reshape(score_mat25a_test.shape[0], 1), score_mat25_ave_spline_DR, q_thresh, r_thresh)
				score_mat25a_KNN[:,j] = score_mat25_SQT_KNNa[:,0]
				print('get KNN max id')
				KNN_id12[:,i+j] = KNN_ida[:,0]
			print('get 48 dimension mat')
			score_mat25_SQT_KNNa_mat_f = get_h5mat(score_mat25a_KNN[:,0],score_mat25a_KNN[:,1], bin_size)
			score_new_SQT_KNN[i,:,0:24] = score_mat25_SQT_KNNa_mat_f

		del(score_mat25_ave_spline)
		del(score_mat25_ave_spline_DR)


		### read ave Dimension reduced mat
		with h5py.File(DNase_ave_r_DR_file, "r") as infile:
			score_mat25_ave_spline_DR = infile["ave_score"][...]

		infile.close()

		### read ave mat
		with h5py.File(DNase_ave_r_file, "r") as infile:
			score_mat25_ave_spline = infile["ave_score"][...]

		infile.close()


		for i in range(0, score_new.shape[0]):
			print(i)
			score_mat25a_KNN = np.copy(score_new[i,:,0:2])
			### reverse use column 2 & 3
			for j in range(2,4):
				print('bin set a')
				print('find the closest average sample')
				center_25bp_cola = center_25bp_all[j]
				score_mat25a_test = score_new[i,:,center_25bp_cola].reshape(score_new.shape[1], 1)*25
				sample2ave_r = corr2d_coeff(score_mat25a_test+sn, score_mat25_ave_spline+sn)
				sample2ave_r_max_id = np.argmax(sample2ave_r, axis=1)
				print('initialize QTnorm sample to Spline ave')
				score_mat25_SQTa = np.copy(score_mat25a_test)
				print('get QTnorm sample to Spline ave')
				score_mat25_SQTa = QTnorm(score_mat25_SQTa[:,0], score_mat25_ave_spline[:,sample2ave_r_max_id[0]])
				###### KNN norm
				print('KNN norm for each CP')
				score_mat25_SQT_KNNa = np.copy(score_mat25_SQTa)
				print('get KNN norm matrix')
				score_mat25_SQT_KNNa, KNN_ida = KNN_norm(CP_id_new, score_mat25a_test[:,0].reshape(score_mat25a_test.shape[0], 1), score_mat25_SQTa.reshape(score_mat25a_test.shape[0], 1), score_mat25_ave_spline_DR, q_thresh, r_thresh)
				score_mat25a_KNN[:,j-2] = score_mat25_SQT_KNNa[:,0]
				KNN_id12[:,i+j+4-2] = KNN_ida[:,0]
			print('get 48 dimension mat')
			score_mat25_SQT_KNNa_mat_r = get_h5mat(score_mat25a_KNN[:,0],score_mat25a_KNN[:,1], bin_size)
			score_new_SQT_KNN[i,:,24:48] = score_mat25_SQT_KNNa_mat_r


		del(score_mat25_ave_spline)
		del(score_mat25_ave_spline_DR)

		score_new_SQT_KNN = np.round(score_new_SQT_KNN, 3)


		if os.path.isfile(output_h5file+'KNNid.h5'):
			os.remove(output_h5file+'KNNid.h5')
		with h5py.File(output_h5file+'KNNid.h5', "w") as outfile:
			outfile.create_dataset(chrom, data=np.array(KNN_id12), compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
		outfile.close()


		print('write SQTKNN normalized h5 file')
		if os.path.isfile(output_h5file):
			os.remove(output_h5file) 
		with h5py.File(output_h5file, "w") as outfile:
			outfile.create_dataset(chrom, data=np.array(score_new_SQT_KNN), compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
			outfile.create_dataset('%s_starts' % (chrom), data=starts, compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
			outfile.create_dataset('feature_names', data=features, compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
			outfile.create_dataset('samples', data=samples, compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)

		outfile.close()




if __name__ == '__main__':
	fire.Fire(SQTKNNN_norm)



