#!/bin/bash

#SBATCH --mem=80G
#SBATCH --cpus-per-task=1

source ~/.bashrc

cd $1
#/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human_macrophage_NFKB1

thread_num=4
#declare -a chr_all=("chr5")
chr_all=$2
config_file='config.yaml'

batch=10000000


train_cl_ct_tis=$3
#'L1236_B_Lymphocyte_Blood'
test_cl_ct_tis=$4
DNase_h5_folder=$5
#'Macrophage_None_Blood'

knn_n=5

scripts_dir='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/scripts/'

### DNase

echo 5 generate_DNase_features_from_bw_to_h5
mkdir -p hdf5s
mkdir -p hdf5s/DNase
mkdir -p hdf5s/DNase/$train_cl_ct_tis'_'$test_cl_ct_tis

echo 6 combine_dnase_features_h5_for_all_cell_types
for chrom_j in "${chr_all[@]}"
do
	echo $chrom_j
	train_id_file=$train_cl_ct_tis'.id.list.txt'
	test_id_file=$test_cl_ct_tis'.id.list.txt'
	python $scripts_dir/DNase_features.py combine_dnase_features_h5_for_all_cell_types \
	--train_id_file=$train_id_file \
	--test_id_file=$test_id_file \
	--chrom=$chrom_j  --outfile_path='./hdf5s/DNase/'$train_cl_ct_tis'_'$test_cl_ct_tis --infile_path=$DNase_h5_folder \
	--config_file=$config_file
done



for chrom_j in "${chr_all[@]}"
do
        echo $chrom_j
        input_h5_file='hdf5s/DNase/'$train_cl_ct_tis'_'$test_cl_ct_tis'/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_'$chrom_j'_all_cell_types.h5'
        output_h5file='hdf5s/DNase/'$train_cl_ct_tis'_'$test_cl_ct_tis'/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_'$chrom_j'_all_cell_types.QT.ave.h5'
        train_id_file=$train_cl_ct_tis'.id.list.txt'
        test_id_file=$test_cl_ct_tis'.id.list.txt'
        time python $scripts_dir/merge_train_test.qt.human2mouse.py \
        --input_h5file=$input_h5_file \
        --output_h5file=$output_h5file \
        --chrom=$chrom_j \
        --train_id_file=$train_id_file \
        --test_id_file=$test_id_file
done


