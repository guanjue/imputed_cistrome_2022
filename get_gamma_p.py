import fire
import pandas as pd
import numpy as np
from scipy.stats import gamma
from scipy.stats import poisson

def get_gammap(x, mostfreq):
	### remove lower half
	used = (x>np.quantile(x,0.5))
	x = x-mostfreq
	x[x<0] = 0
	### get gamma parameters
	scale = np.var(x[used])/np.mean(x[used])
	shape = np.mean(x[used]) / scale
	### get -log10(p)
	neglog10p = -np.log10(1-gamma.cdf(x, shape, scale=scale))
	neglog10p[neglog10p>16] = 16	
	#neglog10p = gamma.cdf(x, shape, scale=scale)
	#neglog10p = neglog10p.reshape(neglog10p.shape[0], 1)
	return(neglog10p)

def get_sigmoid(x):
	return(1/(1+np.exp(-x)))


def get_poisp(x, mostfreq):
	### remove lower half
	x = x-mostfreq
	x[x<0] = 0
	mean = np.mean(x)
	if mean<1:
		mean = 1
	### get -log10(p)
	pvalue = 1-poisson.cdf(x, mu=mean)
	pvalue[pvalue<1e-16] = 1e-16
	neglog10p = -np.log10(pvalue)
	#neglog10p[neglog10p>16] = 16
	return(neglog10p)



#bins='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.chr1_22.bins.bed'
#data_dir = 'lgbmodel/NFKB1/'
#input_file_start = 'pred.5.chr2.QT_top.setA.cross.NFKB1.L1236_B_Lymphocyte_Blood.Macrophage_None_Blood.train.chr3.'
#input_file_end = '.lgb.pred.txt'
#outputname = 'predict.NFKB1.bed'

def get_predict_bed(data_dir, input_file_start, input_file_end, bins, outputname):
	### get bin file
	bin_mat = pd.read_csv(bins, header=None, sep='\t').values

	### get all chr
	chr_all = np.unique(bin_mat[:,0])

	### get predicted scores
	chr_i = chr_all[0]
	input_file = input_file_start+chr_i+input_file_end
	pred_mat = pd.read_csv(data_dir + input_file, header=None).values

	for chr_i in chr_all[1:]:
		print(chr_i)
		input_file_i = input_file_start+chr_i+input_file_end
		pred_mat_i = pd.read_csv(data_dir + input_file_i, header=None).values
		pred_mat = np.concatenate((pred_mat, pred_mat_i), axis=0)

	### get most frequent element
	count = np.unique(pred_mat[:,0], return_counts=True)
	mostfreq = count[0][np.argmax(count[1])]

	### get -log10(gamma_p)
	#pred_mat_p_gamma = get_gammap(pred_mat, mostfreq)
	pred_mat_p_gamma = get_poisp(pred_mat, mostfreq)
	#pred_mat_p_gamma = get_sigmoid(pred_mat)

	### write output
	output_mat = np.concatenate((bin_mat, pred_mat_p_gamma), axis=1)
	np.savetxt(outputname, output_mat, fmt='%s', delimiter='\t')


if __name__ == '__main__':
	fire.Fire(get_predict_bed)



