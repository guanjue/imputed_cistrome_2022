#!/bin/bash

#SBATCH --mem=80G
#SBATCH --cpus-per-task=1

source ~/.bashrc

cd $1

scripts_dir='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/scripts/'

TF=$2
train_cl_ct_tis=$3
test_cl_ct_tis=$4
bins=$5

dir_out='lgbmodel/'$TF'/'$train_cl_ct_tis'_'$test_cl_ct_tis'_200bp/'


/usr/bin/time --verbose python $scripts_dir/get_gamma_p.py \
--data_dir=$dir_out \
--input_file_start=$TF'.'$train_cl_ct_tis'.'$test_cl_ct_tis'.' \
--input_file_end='.lgb.pred.txt' \
--bins=$bins \
--outputname=$dir_out'predict.'$TF'.bed'

 
