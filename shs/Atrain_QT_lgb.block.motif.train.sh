#!/bin/bash

#SBATCH --mem=200G
#SBATCH --cpus-per-task=1

source ~/.bashrc

cd $1

scripts_dir='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/scripts/'

TF=$2
train_cl_ct_tis=$3
test_cl_ct_tis=$4

topN=10000
knn_n=5
ct_id_i='train'

echo $train_cl_ct_tis
echo $test_cl_ct_tis

#label_file_train='peaks/5fc/'$TF'.train.qthresh.label.tsv'
label_file_train='peaks/5fc/'$TF'.train.topN.label.tsv'
label_file_test='peaks/5fc/'$TF'.train.topN.label.tsv'
label_file_train_topN='peaks/5fc/'$TF'.train.topN.label.tsv'


chrom_j='chr2'
chrom_j_test='chr1'



echo $chrom_j
DNase_h5='hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_'
DNase_h5_test='hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_'
motif_h5='train/'$TF'/selected_motif_hdf5/'
motif_h5_test='train/'$TF'/selected_motif_hdf5/'
dir_out='lgbmodel/'$TF'/'

mkdir lgbmodel
mkdir $dir_out
output_model_file=$dir_out'/'$TF'.'$ct_id_i'.'$chrom_j'.model'$knn_n'.pkl'

rm $output_model_file

/usr/bin/time --verbose python $scripts_dir/train.crosschr.setAsub.blockmotif.train.py train_lgb_model --chrom=$chrom_j --chrom_test=$chrom_j_test --TFname=$TF \
--cell_line=$ct_id_i --DNase_h5=$DNase_h5 --motif_h5=$motif_h5 --motif_h5_test=$motif_h5_test \
--label_file_train=$label_file_train --label_file_train_topN=$label_file_train_topN \
--label_file_test=$label_file_test --km_n=5 \
--DNase_h5_test=$DNase_h5_test \
--sigtype='QT' \
--train_ct=$train_cl_ct_tis \
--test_ct=$test_cl_ct_tis \
--dir_out=$dir_out'pred.'$knn_n'.'$chrom_j'.QT_top.setA.cross.motifonly.'




