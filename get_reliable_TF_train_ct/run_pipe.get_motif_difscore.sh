#cut -f3 /liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls | sort -u > TF_list.txt
working_dir='/homes1/gxiang/projects/impute_cistrome/get_motif_difscore/'

for tf_i in $(cat TF_list.txt)
do
	echo $tf_i
	###sbatch shell_pipe.get_TOP_pkbin.sh $working_dir $tf_i
	sbatch shell_pipe.get_true_label.topN.sh $working_dir $tf_i
	sbatch shell_pipe.motif_scan.sh $working_dir $tf_i
	###sbatch get_TF_bincount.sh $working_dir $tf_i
done
