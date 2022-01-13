

### modified parameters
TF=$1
train_ct=$2
test_ct=$3
user_train_DNase_list=$4
user_test_DNase_list=$5
working_dir0=$6
script_dir=$7
DNase_h5_folder=$8
TF_json_folder=$9
pksumit_folder=${10}
pk_folder=${11}
DNase_info=${12}
TF_info=${13}
bins=${14}
motif_h5_folder=${15}

### mkdir folders
mkdir -p $working_dir0
cd $working_dir0
mkdir -p $TF'_'$train_ct'_'$test_ct
cd $TF'_'$train_ct'_'$test_ct
working_dir=$working_dir0'/'$TF'_'$train_ct'_'$test_ct'/'
cd $working_dir
mkdir -p hdf5s
mkdir -p hdf5s/DNase

### fixed parameters
all_chr=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X")
config_file=$script_dir'/config_file/config.yaml'
train_test_header_file=$script_dir'/config_file/train_test_header_file.txt'
cp $config_file ./

###### get train cell type TF id list
time bash $script_dir/shs/get_DNase_train_test_idlist.sh $train_ct $test_ct $DNase_info $DNase_h5_folder
time bash $script_dir/shs/get_TF_train_idlist.sh $TF $train_ct $TF_info $pksumit_folder

### check if use user defined DNase list
if [ $user_train_DNase_list == "F" ]
then
	echo no user defined DNase id list
else
	cp $user_train_DNase_list $train_ct'.id.list.txt'
fi

if [ $user_test_DNase_list == "F" ]
then
	echo no user defined DNase id list
else
	cp $user_test_DNase_list $test_ct'.id.list.txt'
fi



###### get train labels
#motif_scan_jobs_ids+=":"
motif_scan_jobs_ids=$(sbatch $script_dir/shs/get_true_label.sh $working_dir $TF $train_ct $pksumit_folder $pk_folder $train_test_header_file $bins $script_dir/shs/ | cut -d ' ' -f4)
preprocessing_jobs_ids=$motif_scan_jobs_ids

###### get motif h5 files
i=${all_chr[0]}
chrom='chr'$i
motif_score_jobs_ids=$(sbatch $script_dir/shs/shell_pipe.motif2h5.sh $working_dir $chrom $TF $train_ct $TF_json_folder $script_dir $motif_h5_folder | cut -d ' ' -f4)
preprocessing_jobs_ids+=":"
preprocessing_jobs_ids+=$motif_score_jobs_ids
for i in "${all_chr[@]:1}"
do
	chrom='chr'$i
	motif_job_i=$(sbatch $script_dir/shs/shell_pipe.motif2h5.sh $working_dir $chrom $TF $train_ct $TF_json_folder $script_dir $motif_h5_folder | cut -d ' ' -f4)
	motif_score_jobs_ids+=":"
	motif_score_jobs_ids+=$motif_job_i
	preprocessing_jobs_ids+=":"
	preprocessing_jobs_ids+=$motif_job_i
done

###### get DNase QT h5 files
i=${all_chr[0]}
chrom='chr'$i
DNase_QTh5_job_ids=$(sbatch $script_dir/shs/shell_pipe.bam2h5.chrN.sh $working_dir $chrom $train_ct $test_ct $DNase_h5_folder $script_dir | cut -d ' ' -f4)
preprocessing_jobs_ids+=":"
preprocessing_jobs_ids+=$DNase_QTh5_job_ids
for i in "${all_chr[@]:1}"
do
	chrom='chr'$i
	echo $chrom
	DNase_job_i=$(sbatch $script_dir/shs/shell_pipe.bam2h5.chrN.sh $working_dir $chrom $train_ct $test_ct $DNase_h5_folder $script_dir | cut -d ' ' -f4)
	DNase_QTh5_job_ids+=":"
	DNase_QTh5_job_ids+=$DNase_job_i
	preprocessing_jobs_ids+=":"
	preprocessing_jobs_ids+=$DNase_job_i
done

###### get trained model and predictions
echo $preprocessing_jobs_ids
pred_job_id=$(sbatch --dependency=afterok:${preprocessing_jobs_ids} $script_dir/shs/Atrain_QT_lgb.setAsub.sh $working_dir $TF $train_ct $test_ct $bins $script_dir | cut -d ' ' -f4)

###### convert predict score to 
pred2score_job_id=$(sbatch --dependency=afterok:${pred_job_id} $script_dir/shs/convert_score_to_output.sh $working_dir $TF $train_ct $test_ct $bins $script_dir | cut -d ' ' -f4)

###### get imputed top peaks
get_impute_top_pks=$(sbatch --dependency=afterok:${pred2score_job_id} $script_dir/shs/get_impute_top_pk.sh $working_dir $script_dir $TF $train_ct'_'$test_ct | cut -d ' ' -f4)

echo 'DNase_QTh5_job_ids:' > $train_ct'.'$test_ct'.'$TF'.jobs.txt'
echo $DNase_QTh5_job_ids >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'
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
echo 'get_impute_top_pks' >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'
echo $get_impute_top_pks >> $train_ct'.'$test_ct'.'$TF'.jobs.txt'





