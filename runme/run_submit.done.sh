###### Step0: get input parameters. Use these scripts to get input parameters
### selecting imputed TFs
#target_TF_list=TF_select.tmp.txt
### selecting impute cell types
#test_ct=NONE_B_LYMPHOCYTE_BLOOD
### filter to select reliable TF list
#reliable_TFs_list=/liulab/gjxiang/projects/impute_cistrome/get_motif_difscore/reliable_TFs.txt
#time bash get_thisRun_list.sh $target_TF_list $test_ct $reliable_TFs_list
### get impute_cistrome pipeline input parameters: impute TF; training cell type; testing cell type
#cat TF_train_test.thisRun.list.txt


###### run impute_cistrome pipeline ######
### set impute TF; training cell type; testing cell type; User defined sample id for DNase data
tf=GATA1
train_ct=K562
test_ct=NONE_B_LYMPHOCYTE_BLOOD
user_train_DNase_list=F
user_test_DNase_list=F

### set working directory
working_dir0='/liulab/gjxiang/projects/impute_cistrome/B_LYMPHOCYTE'
### set script directory
script_dir='/liulab/gjxiang/projects/impute_cistrome/scripts_done'
### set input file directory
DNase_h5_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human/hdf5s/DNase'
TF_json_folder='/liulab/gjxiang/projects/impute_cistrome/get_motif_difscore/TFs'
pksumit_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_summits/'
pk_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_peaks/'
motif_h5_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human/hdf5s/motif'

### DNase-seq metadata
DNase_info=$script_dir'/config_file/cistrome.metadata.Homo_sapiens.DNase.withheader.txt'
### TF ChIP-seq metadata
TF_info=$script_dir'/config_file/cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls'
### hg38 200bp bins with 50bp sliding window
bins=$script_dir'/config_file/hg38.200_50slide.bins.bed'


### submit steps as sbatch jobs ######
bash submit.done.sh $tf $train_ct $test_ct $user_train_DNase_list $user_test_DNase_list $working_dir0 $script_dir $DNase_h5_folder $TF_json_folder $pksumit_folder $pk_folder $DNase_info $TF_info $bins $motif_h5_folder


###### output files ######
#cd $working_dir0
### For each run, the output folder will be named as $tf'_'$train_ct'_'$test_ct , e.g GATA1_K562_NONE_B_LYMPHOCYTE_BLOOD
### DNase 5' end count h5 files
#ls -ltrh GATA1_K562_NONE_B_LYMPHOCYTE_BLOOD/hdf5s/DNase/K562_NONE_B_LYMPHOCYTE_BLOOD/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_chr*_all_cell_types.h5
### DNase 5' end count h5 files Quantile transformed
#ls -ltrh GATA1_K562_NONE_B_LYMPHOCYTE_BLOOD/hdf5s/DNase/K562_NONE_B_LYMPHOCYTE_BLOOD/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_chr*_all_cell_types.QT.ave.h5
### motif score h5 files
#ls -ltrh GATA1_K562_NONE_B_LYMPHOCYTE_BLOOD/train/GATA1/selected_motif_hdf5/chr*_motif_features_lightGBM.h5
### lgb predicted scores
#ls -ltrh GATA1_K562_NONE_B_LYMPHOCYTE_BLOOD/lgbmodel/GATA1/K562_NONE_B_LYMPHOCYTE_BLOOD/predict.GATA1.bed
### FDR peaks
#ls -ltrh GATA1_K562_NONE_B_LYMPHOCYTE_BLOOD/lgbmodel/GATA1/K562_NONE_B_LYMPHOCYTE_BLOOD/predict.GATA1.fdr.bed
### FDR bedtools merged peaks
#ls -ltrh GATA1_K562_NONE_B_LYMPHOCYTE_BLOOD/lgbmodel/GATA1/K562_NONE_B_LYMPHOCYTE_BLOOD/predict.GATA1.fdr.topPK.done.bed







