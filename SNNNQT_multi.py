
import gc
import os
import pickle

import fire
import h5py
import numpy as np
import pandas as pd


def SQTKNNN_multi(chrom, TF, h5file, CP_id_file, l_list, output_h5file, output_label, test_samples, label_folder):
	#chrom = 'chr13'
	#TF = 'GABPA'
	#h5file = 'DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_'+chrom+'_all_cell_types.SQTKNNN5.h5'
	#CP_id_file = chrom+'.50bp.CPid.bed'
	#l_list = TF+'.train.labels.list.txt'
	#output_h5file = 'DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_'+chrom+'_all_cell_types.SQTKNNN5.ave.train_multi.test.h5'
	#output_label = TF+'.train_multi.train.labels.tsv'
	#test_samples = 'liver'

	center_25bp_all = [12,13,36,37]
	#center_25bp_all = [0,1]


	def corr2d_coeff(A, B):
		# Rowwise mean of input arrays & subtract from input arrays themeselves
		A_mA = A - A.mean(0)[None, :]
		B_mB = B - B.mean(0)[None, :]
		ssA = (A_mA**2).sum(0)
		ssB = (B_mB**2).sum(0)
		return(np.round(np.dot(A_mA.T, B_mB) / np.sqrt(np.dot(ssA[:, None],ssB[None])), 3))



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

	### read labels
	label_list = np.array(pd.read_csv(l_list, sep="\t", header=None))

	l1a = np.array(pd.read_csv(label_folder+TF+'.'+label_list[0,0]+'.train.labels.tsv', sep="\t", header=0))
	l_i_multi_new = l1a[l1a[:,0]==chrom,:]
	l_i_multi = l_i_multi_new[:,3]
	l_i_multi = l_i_multi.reshape(l_i_multi.shape[0],1)

	for i in range(1,label_list.shape[0]):
		print(i)
		lia = np.array(pd.read_csv(label_folder+TF+'.'+label_list[i,0]+'.train.labels.tsv', sep="\t", header=0))
		li_i = lia[lia[:,0]==chrom,:][:,3]
		li_i = li_i.reshape(li_i.shape[0],1)
		l_i_multi = np.concatenate((l_i_multi, li_i), axis=1) 


	### read h5 file
	infile = h5py.File(h5file, "r")
	infile.keys()
	score25 = (infile[chrom][...])
	feature_names = np.array(infile['feature_names'])
	samples = np.array(infile['samples']).astype(str)
	starts = infile['%s_starts' % chrom][...]
	infile.close()


	### get train test samples
	train_id = samples == label_list[0,0]
	for i in range(1,label_list.shape[0]):
		print(i)
		train_id = train_id | (samples == label_list[i,0])

	test_id = samples == test_samples

	### get sigmat
	score25_train = score25[train_id,:,:]
	score25_test = score25[test_id,:,:][0,:,:]
	score25_train25 = np.transpose(np.mean(score25_train[:,:,center_25bp_all], axis=2))
	score25_test25 = np.mean(score25_test[:,center_25bp_all], axis=1)
	score25_train25_multi = np.copy(score25_test)


	cp_used = []
	for i in np.unique(CP_id_new):
		print(i)
		### get score CPi
		used_row_i = CP_id_new==i
		score_mat25_train = score25_train25[used_row_i,:]
		score_mat25_test = score25_test25[used_row_i]
		###
		cp_i_sample_r = corr2d_coeff(score_mat25_test.reshape(score_mat25_test.shape[0],1), score_mat25_train)
		#if np.sum(np.isnan(cp_i_sample_r, where=True))==0:
		cp_i_sample = np.argmax(cp_i_sample_r)
		#else:
		#	cp_i_sample = 0
		###
		### get nulti train h5
		score25_train25_multi[used_row_i, :] = score25_train[cp_i_sample, used_row_i, :]
		cp_used.append(cp_i_sample)
		### get multi label
		l_i_multi_new[used_row_i,3] = l_i_multi[used_row_i, cp_i_sample]

	score25_train25_multi[np.isnan(score25_train25_multi)] = 0

	samples_new = np.array(['train', 'test']).astype('bytes_')


	### get output
	score25_multi_train_test = np.concatenate((score25_train25_multi.reshape(1,score25_train25_multi.shape[0],score25_train25_multi.shape[1]), score25_test.reshape(1,score25_test.shape[0],score25_test.shape[1])), axis=0)

	### write output
	with h5py.File(output_h5file, "w") as outfile:
		print('score')
		outfile.create_dataset(chrom, data=np.array(score25_multi_train_test), compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
		print('starts')
		outfile.create_dataset('%s_starts' % chrom, data=starts, compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
		print('feature_names')
		outfile.create_dataset('feature_names', data=feature_names, compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
		print('samples')
		outfile.create_dataset('samples', data=samples_new, compression='gzip', shuffle=True, fletcher32=True, compression_opts=9)
	outfile.close()

	### write multi train labels
	np.savetxt(output_label, l_i_multi_new, delimiter="\t", fmt="%s", comments='')



if __name__ == '__main__':
	fire.Fire(SQTKNNN_multi)



