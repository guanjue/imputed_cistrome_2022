
working_dir='/homes1/gxiang/projects/impute_cistrome/macrophage/'
cd $working_dir

mkdir -p hdf5s
mkdir -p hdf5s/DNase

### modified parameters
TF='STAT5A'
train_ct='K562_Erythroblast_Bone_Marrow'
test_ct='None_Macrophage_Blood'

TF=$1
train_ct=$2
test_ct=$3
run_DNase=$4

### fixed parameters
all_chr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
bins='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.chr1_22.bins.bed'

DNase_info='/homes1/gxiang/projects/DNase_cluster/info/cistrome.metadata.Homo_sapiens.DNase.withheader.txt'

TF_info='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls'
TF_json_folder='/homes1/gxiang/projects/impute_cistrome/get_motif_difscore/TFs/'

all_chr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
bins='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.chr1_22.bins.bed'

DNase_info='/homes1/gxiang/projects/DNase_cluster/info/cistrome.metadata.Homo_sapiens.DNase.withheader.txt'
DNase_h5_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human/hdf5s/DNase'

TF_info='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls'
TF_json_folder='/homes1/gxiang/projects/impute_cistrome/get_motif_difscore/TFs'

pksumit_folder='/liulab/jfan/projects/impute_cistrome/cistrome/human_TF_summits/'
pk_folder='/liulab/jfan/projects/impute_cistrome/cistrome/human_TF_peaks/'
train_test_header_file='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/train_test_header_file.txt'

###### get train cell type TF id list
time bash get_DNase_train_test_idlist.sh $train_ct $test_ct $DNase_info $DNase_h5_folder
time bash get_TF_train_idlist.sh $TF $train_ct $TF_info $pksumit_folder

###### get trained model and predictions
pred_job_id=$(sbatch Atrain_QT_lgb.setAsub.sh $working_dir $TF $train_ct $test_ct $bins | cut -d ' ' -f4)

#sbatch Atrain_QT_lgb.setAsub.sh $working_dir $TF $train_ct $test_ct
#shuffle_job_id=$(sbatch Atrain_QT_lgb.setAsub_motifonly.train.sh $working_dir $TF $train_ct $train_ct | cut -d ' ' -f4)
#time bash Atrain_QT_lgb.setAsub.sh $working_dir $TF $train_ct $test_ct

###### convert predict score to 
pred2score_job_id=$(sbatch --dependency=afterok:${pred_job_id} convert_score_to_output.sh $working_dir $TF $train_ct $test_ct $bins | cut -d ' ' -f4)

