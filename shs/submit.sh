

working_dir0='/homes1/gxiang/projects/impute_cistrome/macrophage/'
### modified parameters
TF=$1
train_ct=$2
test_ct=$3
run_DNase=$4

mkdir -p $working_dir0'/'$TF'_'$train_ct'_'$test_ct
working_dir=$working_dir0'/'$TF'_'$train_ct'_'$test_ct'/'
cd $working_dir

mkdir -p hdf5s
mkdir -p hdf5s/DNase

### fixed parameters
all_chr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
bins='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.chr1_22.bins.bed'

DNase_info='/homes1/gxiang/projects/DNase_cluster/info/cistrome.metadata.Homo_sapiens.DNase.withheader.txt'
DNase_h5_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human/hdf5s/DNase'

TF_info='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls'
TF_json_folder='/homes1/gxiang/projects/impute_cistrome/get_motif_difscore/TFs'

pksumit_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_summits/'
pk_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_peaks/'
train_test_header_file='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/train_test_header_file.txt'

bins='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.bins.bed'

script_dir='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/scripts/shs/'
config_file='/homes1/gxiang/projects/impute_cistrome/macrophage/config.yaml'

cp $config_file ./

###### get train cell type TF id list
time bash $script_dir/get_DNase_train_test_idlist.sh $train_ct $test_ct $DNase_info $DNase_h5_folder
time bash $script_dir/get_TF_train_idlist.sh $TF $train_ct $TF_info $pksumit_folder

###### get motif scan input files
#motif_scan_jobs_ids=$(sbatch shell_pipe.motif_scan.sh $working_dir $TF $train_ct | cut -d ' ' -f4)

###### get train labels
#motif_scan_jobs_ids+=":"
motif_scan_jobs_ids=$(sbatch $script_dir/get_true_label.sh $working_dir $TF $train_ct $pksumit_folder $pk_folder $train_test_header_file $bins $script_dir | cut -d ' ' -f4)
preprocessing_jobs_ids=$motif_scan_jobs_ids

###### get motif h5 files
i=${all_chr[0]}
chrom='chr'$i
motif_score_jobs_ids=$(sbatch $script_dir/shell_pipe.motif2h5.sh $working_dir $chrom $TF $train_ct $TF_json_folder | cut -d ' ' -f4)
preprocessing_jobs_ids+=":"
preprocessing_jobs_ids+=$motif_score_jobs_ids
for i in "${all_chr[@]:1}"
do
	chrom='chr'$i
	motif_job_i=$(sbatch $script_dir/shell_pipe.motif2h5.sh $working_dir $chrom $TF $train_ct $TF_json_folder | cut -d ' ' -f4)
	motif_score_jobs_ids+=":"
	motif_score_jobs_ids+=$motif_job_i
	preprocessing_jobs_ids+=":"
	preprocessing_jobs_ids+=$motif_job_i
done

if [ $run_DNase == "T" ]
then
###### get DNase QT h5 files
i=${all_chr[0]}
chrom='chr'$i
DNase_QTh5_job_ids=$(sbatch $script_dir/shell_pipe.bam2h5.chrN.sh $working_dir $chrom $train_ct $test_ct $DNase_h5_folder | cut -d ' ' -f4)
preprocessing_jobs_ids+=":"
preprocessing_jobs_ids+=$DNase_QTh5_job_ids
for i in "${all_chr[@]:1}"
do
	chrom='chr'$i
	echo $chrom
	DNase_job_i=$(sbatch $script_dir/shell_pipe.bam2h5.chrN.sh $working_dir $chrom $train_ct $test_ct $DNase_h5_folder | cut -d ' ' -f4)
	DNase_QTh5_job_ids+=":"
	DNase_QTh5_job_ids+=$DNase_job_i
	preprocessing_jobs_ids+=":"
	preprocessing_jobs_ids+=$DNase_job_i
done
fi

###### get trained model and predictions
echo $preprocessing_jobs_ids
pred_job_id=$(sbatch --dependency=afterok:${preprocessing_jobs_ids} $script_dir/Atrain_QT_lgb.setAsub.sh $working_dir $TF $train_ct $test_ct $bins | cut -d ' ' -f4)

#sbatch Atrain_QT_lgb.setAsub.sh $working_dir $TF $train_ct $test_ct
#shuffle_job_id=$(sbatch Atrain_QT_lgb.setAsub_motifonly.train.sh $working_dir $TF $train_ct $train_ct | cut -d ' ' -f4)
#time bash Atrain_QT_lgb.setAsub.sh $working_dir $TF $train_ct $test_ct

###### convert predict score to 
pred2score_job_id=$(sbatch --dependency=afterok:${pred_job_id} $script_dir/convert_score_to_output.sh $working_dir $TF $train_ct $test_ct $bins | cut -d ' ' -f4)

echo 'DNase_QTh5_job_ids:' > $train_ct'.'$test_ct'.'$TF'.jobs.txt'
if [ $run_DNase == "T" ]
then
	echo $DNase_QTh5_job_ids >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'
fi

echo 'motif_scan_jobs_ids' >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'
echo $motif_scan_jobs_ids >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'
echo 'motif_score_jobs_ids' >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'
echo $motif_score_jobs_ids >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'
echo 'pred_job_id' >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'
echo $pred_job_id >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'
echo 'pred2score_job_id' >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'
echo $pred2score_job_id >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'
echo 'shuffle_job_id' >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'
echo $shuffle_job_id >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'





