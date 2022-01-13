
working_dir='/homes1/gxiang/projects/impute_cistrome/K562_GM12878/'
cd $working_dir

mkdir hdf5s
mkdir hdf5s/DNase

### modified parameters
TF='STAT5A'
train_ct='K562'
test_ct='Macrophage_None_Blood'

TF=$1
train_ct=$2
test_ct=$3
run_DNase=$4
cell_type_col=$5

### fixed parameters
all_chr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
bins='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.chr1_22.bins.bed'

TF_info='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls'


###### get train cell type TF id list
#time bash get_TF_train_idlist.sh $TF $train_ct $cell_type_col $TF_info

###### motif_scan_jobs_ids
#motif_scan_jobs_ids=$(sbatch shell_pipe.motif_scan.sh $working_dir $TF $train_ct | cut -d ' ' -f4)

###### get train labels
#motif_scan_jobs_ids+=":"
#motif_scan_jobs_ids+=$(sbatch get_true_label.sh $working_dir $TF $train_ct | cut -d ' ' -f4)

###### get motif h5 files
i=${all_chr[0]}
chrom='chr'$i
motif_score_jobs_ids=$(sbatch shell_pipe.motif2h5.sh $working_dir $chrom $TF $train_ct | cut -d ' ' -f4)
preprocessing_jobs_ids=$motif_score_jobs_ids
for i in "${all_chr[@]:1}"
do
        chrom='chr'$i
        motif_job_i=$(sbatch shell_pipe.motif2h5.sh $working_dir $chrom $TF $train_ct | cut -d ' ' -f4)
        motif_score_jobs_ids+=":"
        motif_score_jobs_ids+=$motif_job_i
        preprocessing_jobs_ids+=":"
        preprocessing_jobs_ids+=$motif_job_i
done

###### get trained model and predictions
echo $preprocessing_jobs_ids
#pred_job_id=$(sbatch --dependency=afterok:${preprocessing_jobs_ids} Atrain_QT_lgb.setAsub.sh $working_dir $TF $train_ct $test_ct | cut -d ' ' -f4)

#sbatch Atrain_QT_lgb.setAsub.sh $working_dir $TF $train_ct $test_ct
shuffle_job_id=$(sbatch --dependency=afterok:${preprocessing_jobs_ids} Atrain_QT_lgb.setAsub_motifonly.train.sh $working_dir $TF $train_ct $train_ct | cut -d ' ' -f4)
#time bash Atrain_QT_lgb.setAsub.sh $working_dir $TF $train_ct $test_ct

###### convert predict score to 
#pred2score_job_id=$(sbatch --dependency=afterok:${pred_job_id} convert_score_to_output.sh $working_dir $TF $train_ct $test_ct $bins | cut -d ' ' -f4)



