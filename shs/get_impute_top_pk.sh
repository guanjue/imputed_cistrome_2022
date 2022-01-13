#!/bin/bash

#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

source ~/.bashrc
conda activate impute_cistrome_2022

cd $1
script_dir=$2
tf=$3
train_test_ct=$4

time Rscript $script_dir/get_FDR.R 'lgbmodel/'$tf'/'$train_test_ct'/predict.'$tf'.bed' 'lgbmodel/'$tf'/'$train_test_ct'/predict.'$tf'.fdr.bed'
sort -k4,4nr 'lgbmodel/'$tf'/'$train_test_ct'/predict.'$tf'.fdr.bed' | head -100000 | sort -k1,1 -k2,2n > 'lgbmodel/'$tf'/'$train_test_ct'/predict.'$tf'.fdr.top.bed'
bedtools merge -i 'lgbmodel/'$tf'/'$train_test_ct'/predict.'$tf'.fdr.top.bed' -c 4 -o mean > 'lgbmodel/'$tf'/'$train_test_ct'/predict.'$tf'.fdr.top.M.bed'
sort -k4,4nr 'lgbmodel/'$tf'/'$train_test_ct'/predict.'$tf'.fdr.top.M.bed' | head -10000 | sort -k1,1 -k2,2n > 'lgbmodel/'$tf'/'$train_test_ct'/predict.'$tf'.fdr.topPK.done.bed'
