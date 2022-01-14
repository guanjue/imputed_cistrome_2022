#cut -f3,9 /liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls | sort -u | awk -F '\t' -v OFS='\t' '{print $1"_"$2}' > TF_CL_list.txt
#cut -f7,8,9 /homes1/gxiang/projects/DNase_cluster/info/cistrome.metadata.Homo_sapiens.DNase.txt | awk -F '\t' -v OFS='\t' '{print $1"_"$2"_"$3}' | sort -u > DNase.CL_CT_TIS.txt

tail -n+2 /liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls | cut -f3,9,10,11 | sort -u | awk -F '\t' -v OFS='\t' '{print $1"_"$2"_"$3"_"$4}' > TF_CL_CT_TIS_list.txt
cut -f7,8,9,27 /homes1/gxiang/projects/DNase_cluster/info/cistrome.metadata.Homo_sapiens.DNase.txt | awk -F '\t' -v OFS='\t' '{if ($4>0.1) print $1"_"$2"_"$3}' | sort -u> DNase.CL_CT_TIS.txt


working_dir='/homes1/gxiang/projects/impute_cistrome/get_motif_difscore/'

# run TF cl merged first run_pipe.get_motif_difscore.sh

#for tf_i in $(cat TF_CL_list.noMYC.3.txt)
for tf_i in $(cat TF_CL_CT_TIS_list.3.txt)
do
	echo $tf_i
	sbatch shell_pipe.get_true_label.topN.CL.sh $working_dir $tf_i
done
