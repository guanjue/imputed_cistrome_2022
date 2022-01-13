thread_num=4
TF='FOXA1'
declare -a chr_all=("chr19")
declare -a cell_line_id_all=("40606" "8498")
declare -a cell_line_id_train=("40606")
declare -a cell_line_id_pred=("8498")
cell_line_id_ref="40606"
bam_path='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human7/bam_files/'
config_file='config.yaml'
peak_file_path='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human6/cistrome/peak/'

motif_pwm_path='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/motifs/'
chrom_size_file='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human7/cistrome_impute_hg38_config/hg38.chr19.chrom.sizes'
genome_sequence_fa='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human7/cistrome_impute_hg38_config/chr19.fa'
chrom_all=$chr_all
batch=10000000

seqpos_env_path='/liulab/gjxiang/miniconda3/envs/seqpos/'
motif_database_xml='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/HOCOMOCOv11_full_pwm_HUMAN_mono.xml'
genome='hg38'

label_file_train='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human6/label/train/FOXA1.train.labels.tsv'
label_file_test='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human6/label/train/FOXA1.train.labels.tsv'

### DNase

echo 1 bam_filter
mkdir bams
mkdir bams/with_filter
for ct_id_i in "${cell_line_id_all[@]}"
do
	echo $ct_id_i
	input_bam_i=$bam_path/$ct_id_i.bam
	output_bam_i=bams/with_filter/DNASE.$ct_id_i.merge.filter.bam
	samtools view -F 1804 -q 30 --threads $thread_num -b $input_bam_i -o $output_bam_i
done

echo 2 bam sort
for ct_id_i in "${cell_line_id_all[@]}"
do
	echo $ct_id_i
	input_bam_i=bams/with_filter/DNASE.$ct_id_i.merge.filter.bam
	output_bam_i=bams/with_filter/DNASE.$ct_id_i.merge.filter.sorted.bam
	samtools sort $input_bam_i -o $output_bam_i --threads $thread_num
done

echo 3 bam_index
for ct_id_i in "${cell_line_id_all[@]}"
do
	echo $ct_id_i
	input_bam_i=bams/with_filter/DNASE.$ct_id_i.merge.filter.sorted.bam
#	samtools index $input_bam_i
done

echo 4 generate_bw_5_end_1_bp
mkdir bigwigs
for ct_id_i in "${cell_line_id_all[@]}"
do
	echo $ct_id_i
	input_bam_i=bams/with_filter/DNASE.$ct_id_i.merge.filter.sorted.bam
	forward_bw="bigwigs/DNASE."$ct_id_i".merge.binSize.1.forward.bw"
	reverse_bw="bigwigs/DNASE."$ct_id_i".merge.binSize.1.reverse.bw"
	bamCoverage -b $input_bam_i -o $forward_bw --Offset 1 \
	--samFlagExclude 16  --binSize 1 --numberOfProcessors $thread_num --ignoreDuplicates 
	bamCoverage -b $input_bam_i -o $reverse_bw --Offset 1 \
	--samFlagInclude 16  --binSize 1 --numberOfProcessors $thread_num --ignoreDuplicates 
done

echo 5 generate_DNase_features_from_bw_to_h5
mkdir hdf5s
mkdir hdf5s/DNase
for ct_id_i in "${cell_line_id_all[@]}"
do
	echo $ct_id_i
	for chrom_j in "${chr_all[@]}"
	do
		echo $chrom_j
		python scripts/DNase_features.py generate_dnase_features_from_bw_to_h5 --chrom=$chrom_j \
		--cell_line=$ct_id_i --dnase_bw_file_path='./bigwigs' --outfile_path='./hdf5s/DNase/' \
		--config_file=$config_file
	done
done

echo 6 combine_dnase_features_h5_for_all_cell_types
for chrom_j in "${chr_all[@]}"
do
	echo $chrom_j
	python scripts/DNase_features.py combine_dnase_features_h5_for_all_cell_types \
	--chrom=$chrom_j  --outfile_path='./hdf5s/DNase' \
	--config_file=$config_file
done


### motif
echo 7 prepare_motif_top4_feature #(very slow)
mkdir hdf5s/motif
for chrom_j in "${chr_all[@]}"
do
	echo $chrom_j
	python scripts/motif_features.py MotifFeatures prepare_motif_top4_feature --chrom=$chrom_j --num_threads=$thread_num \
	--motif_pwm_path=$motif_pwm_path --motif_feature_path='./hdf5s/motif' \
	--chrom_size_file=$chrom_size_file --chrom_all=$chrom_all --batch=$batch \
	--genome_sequence_fa=$genome_sequence_fa
done

echo 8 top_5k_peaks
for ct_id_i in "${cell_line_id_train[@]}"
do
	echo $ct_id_i
	input_pk=$peak_file_path'/'$ct_id_i'.'$TF'.sort_peaks.narrowPeak.bed'
	output_pk='peaks/top5k/ChIPseq.'$ct_id_i'.'$TF'.conservative.train.top5k.narrowPeak'
	cat $input_pk | head -n 5000 > $output_pk 
done

echo 9 seqpos_scan_motifs_on_peaks
for ct_id_i in "${cell_line_id_train[@]}"
do
	echo $ct_id_i
	input_pk='peaks/top5k/ChIPseq.'$ct_id_i'.'$TF'.conservative.train.top5k.narrowPeak'
	$seqpos_env_path'/bin/python2.7' $seqpos_env_path'/bin/MDSeqPos.py' $input_pk $genome \
	-m $motif_database_xml -d \
	-O './train/'$TF'/motif_scan/cell_line.'$ct_id_i'.tf.'$TF
done

echo 10 decide_cofactor_motif_set
for ct_id_i in "${cell_line_id_train[@]}"
do
	echo $ct_id_i
	input_file='train/'$TF'/motif_scan/cell_line.'$ct_id_i'.tf.'$TF'/mdseqpos_index.html'
	output_file='train/'$TF'/'$TF'_cofactor_motif_list.json'
	topn=10
	python scripts/motif_features.py decide_cofactor_motif_set --outfilename=$output_file --n=$topn $input_file
done

echo 11 prepare_motif_h5_data
for chrom_j in "${chr_all[@]}"
do
	echo $chrom_j
	cofactor_motif_set_file='train/'$TF'/'$TF'_cofactor_motif_list.json'
	rm 'train/'$TF'/selected_motif_hdf5/'$chrom_j'_motif_features_lightGBM.h5'
	python scripts/train_lightgbm_model.py prepare_motif_h5_data --chrom=$chrom_j \
	--config_file=$config_file --cofactor_motif_set_file=$cofactor_motif_set_file \
	--training_tf_name=$TF \
	--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' \
	--motif_feature_path='./hdf5s/motif' \
	--selected_motif_feature_path='train/'$TF'/selected_motif_hdf5/' \
	--step=120 
done

echo 12 train data
mkdir lgbmodel
mkdir lgbmodel/$TF

for ct_id_i in "${cell_line_id_train[@]}"
do
	echo $ct_id_i
	for chrom_j in "${chr_all[@]}"
	do
		echo $chrom_j
		DNase_h5='hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_'$chrom_j'_all_cell_types.h5'
		motif_h5='train/'$TF'/selected_motif_hdf5/'$chrom_j'_motif_features_lightGBM.h5'
		dir_out='lgbmodel/'$TF'/'
		output_model_file=$dir_out'/'$TF'.'$ct_id_i'.'$chrom_j'.model.pkl'
		rm $output_model_file
		time python scripts/train.py train_lgb_model --chrom=$chrom_j --TFname=$TF \
		--cell_line=$ct_id_i --DNase_h5=$DNase_h5 --motif_h5=$motif_h5 --label_file_train=$label_file_train \
		--dir_out=$dir_out
	done
done


echo 13 predict TF binding
for ct_id_i in "${cell_line_id_pred[@]}"
do
	echo $ct_id_i
	for chrom_j in "${chr_all[@]}"
	do
		echo $chrom_j
		DNase_h5='hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_'$chrom_j'_all_cell_types.h5'
		motif_h5='train/'$TF'/selected_motif_hdf5/'$chrom_j'_motif_features_lightGBM.h5'
		dir_out='lgbmodel/'$TF'/'
		time python scripts/pred.py pred_lgb_model --chrom=$chrom_j --TFname=$TF \
		--cell_line=$ct_id_i --DNase_h5=$DNase_h5 --motif_h5=$motif_h5 --label_file_test=$label_file_test \
		--trained_lgb_model='lgbmodel/'$TF'/'$TF'.'$cell_line_id_ref'.'$chrom_j'.model.pkl' \
		--dir_out=$dir_out
	done
done


