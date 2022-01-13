#!/bin/bash

#SBATCH --mem=100G
#SBATCH --cpus-per-task=1

source ~/.bashrc
conda activate impute_cistrome_2022

cd $1
chr_all=$2
TF=$3
train_cl_ct_tis=$4
TF_json_folder=$5
scrip_dir=$6
motif_h5_folder=$7

thread_num=4

declare -a set_type=("train" "test")

config_file='config.yaml'


chrom_all=$chr_all



#scrip_dir='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/scripts/'


mkdir -p train
mkdir -p train/$TF
mkdir -p train/$TF/selected_motif_hdf5

echo 11 prepare_motif_h5_data
for chrom_j in "${chr_all[@]}"
do
	echo $chrom_j
	cofactor_motif_set_file=$TF_json_folder'/'$TF'/'$TF'_cofactor_motif_list.json'
	rm -f 'train/'$TF'/selected_motif_hdf5/'$chrom_j'_motif_features_lightGBM.h5'
	time python $scrip_dir/train_lightgbm_model.py prepare_motif_h5_data --chrom=$chrom_j \
	--config_file=$config_file --cofactor_motif_set_file=$cofactor_motif_set_file \
	--training_tf_name=$TF \
	--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' \
	--motif_feature_path=$motif_h5_folder \
	--selected_motif_feature_path='train/'$TF'/selected_motif_hdf5/' \
	--step=120 
done

