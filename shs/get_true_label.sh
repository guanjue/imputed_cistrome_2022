#!/bin/bash

#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

source ~/.bashrc
conda activate impute_cistrome_2022

cd $1

TF=$2
train_ct=$3

pksumit_folder=$4
pk_folder=$5
train_test_header_file=$6
bins=$7
script_dir=$8

time bash $script_dir/shell_pipe.get_true_label.qthresh.sh $TF $train_ct $pksumit_folder $pk_folder $train_test_header_file $bins
time bash $script_dir/shell_pipe.get_true_label.topN.sh $TF $train_ct $pksumit_folder $pk_folder $train_test_header_file $bins

