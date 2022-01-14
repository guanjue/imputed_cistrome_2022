

library(rjson)
library(rhdf5)

setwd('/homes1/gxiang/projects/impute_cistrome/get_motif_difscore/')

h5ls('/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human/hdf5s/motif/chr22_motifs_top4_scores.h5')
s = h5read('/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human/hdf5s/motif/chr22_motifs_top4_scores.h5','scores')
names = h5read('/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human/hdf5s/motif/chr22_motifs_top4_scores.h5','motif_names')

bins = read.table('/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.chr22.bins.VIS_cCREcount.bed', header=F)
bins = bins[-dim(bins)[1],]
bins_binary = (bins[,4]!=0)*1

TFname = apply(names, 1, function(x) unlist(strsplit(x, '_'))[1])
TF_train = read.table('/homes1/gxiang/projects/impute_cistrome/get_motif_difscore/TF_list.txt', header=F)

TF_done_max = c()
motif_difp_max = c()


TF_done = c()
motif_difp = c()

for (i in 1:dim(TF_train)[1]){
#for (i in 1:5){
	TF_i = TF_train[i,1]
	print(TF_i)
	if (file.exists(paste('/homes1/gxiang/projects/impute_cistrome/get_motif_difscore/TFs/', TF_i, '/', TF_i,'_cofactor_motif_list.json', sep=''))){
		motif_list = fromJSON(file = paste('/homes1/gxiang/projects/impute_cistrome/get_motif_difscore/TFs/', TF_i, '/', TF_i,'_cofactor_motif_list.json', sep=''))
		motif_difp_i = c()
		for (j in 1:length(motif_list)){
			motif_list_j = motif_list[j]
			used_col = which(names==motif_list_j)
			ds = cbind(s[1,,used_col])
			colnames(ds) = names[used_col]
			if (file.exists(paste('/homes1/gxiang/projects/impute_cistrome/get_motif_difscore/peaks/top5k/', TF_i, '.top5k.bincount.chr22.txt', sep=''))){
				bins_pk = scan(paste('/homes1/gxiang/projects/impute_cistrome/get_motif_difscore/peaks/top5k/', TF_i, '.top5k.bincount.chr22.txt', sep=''))[-1]
				bins_pk_binary = (bins_pk!=0)*1
				if (sum(bins_pk_binary==1)>=2){
				if (sum(bins_pk_binary==1)>=500){
				d_pk1 = ds[bins_pk_binary==1,1][sample(sum(bins_pk_binary==1), 500)]
				} else{
				d_pk1 = ds[bins_pk_binary==1,1]
				}
				if (sum(bins_binary==1)>=500){
				d_all1 = ds[bins_binary==1,1][sample(sum(bins_binary==1), 500)]
				} else{
				d_all1 = ds[bins_binary==1,1]
				}
				ttest_result = wilcox.test(d_pk1, d_all1, alternative = "greater")
				motif_difp_i = c(motif_difp_i, -log10(ttest_result$p.value))
				} else{
				motif_difp_i = c(motif_difp_i, 0)
				}
			}
		}
		motif_difp = rbind(motif_difp, motif_difp_i)
		TF_done = c(TF_done, TF_i)
	}
}



rownames(motif_difp) = TF_done
write.table(motif_difp,'motif_difp.txt',quote=F, sep='\t', row.names=T, col.names=F)


