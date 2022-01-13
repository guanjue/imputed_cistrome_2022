#!/bin/bash

#SBATCH --mem=80G
#SBATCH --cpus-per-task=1

source ~/.bashrc

cd $1
#cd /homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human_macrophage_NFKB1


TF=$2
train_cl_ct_tis=$3
thread_num=4
#TF='FOXA1'
#TF='CTCF'

declare -a set_type=("train" "test")

pksumit_folder='/liulab/jfan/projects/impute_cistrome/cistrome/human_TF_summits/'
seqpos_env_path='/liulab/jfan/miniconda3/envs/seqpos/'
motif_database_xml='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/HOCOMOCOv11_full_pwm_HUMAN_mono.xml'
scrip_dir='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/scripts/'
genome='hg38'

mkdir peaks
mkdir peaks/top5k
### motif

echo 8 prepare top5k.narrowPeak.sumit for train and test
tf_id1=$(head -1 $train_cl_ct_tis'.'$TF'.train.id.list.txt')
head -5000 $pksumit_folder$tf_id1'_sort_summits.bed' | cut -f1,2,3 > 'peaks/top5k/ChIPseq.train.'$TF'.conservative.train.top5k.narrowPeak'
for tf_id in $(cat $train_cl_ct_tis'.'$TF'.train.id.list.txt')
do
	head -5000 $pksumit_folder$tf_id'_sort_summits.bed' | cut -f1,2,3 >> 'peaks/top5k/ChIPseq.train.'$TF'.conservative.train.top5k.narrowPeak'
	#bedtools intersect -a 'peaks/top5k/ChIPseq.train.'$TF'.conservative.train.top5k.narrowPeak' -b tmp.bed -wa -u > 'peaks/top5k/ChIPseq.train.'$TF'.conservative.train.top5k.narrowPeak.tmp' \
	#&& mv 'peaks/top5k/ChIPseq.train.'$TF'.conservative.train.top5k.narrowPeak.tmp' 'peaks/top5k/ChIPseq.train.'$TF'.conservative.train.top5k.narrowPeak'
	#rm tmp.bed
done
sort -u 'peaks/top5k/ChIPseq.train.'$TF'.conservative.train.top5k.narrowPeak' > 'peaks/top5k/ChIPseq.train.'$TF'.conservative.train.top5k.narrowPeak.tmp' \
&& mv 'peaks/top5k/ChIPseq.train.'$TF'.conservative.train.top5k.narrowPeak.tmp' 'peaks/top5k/ChIPseq.train.'$TF'.conservative.train.top5k.narrowPeak0'
time python3 $scrip_dir/sample_rows.py --inputfile='peaks/top5k/ChIPseq.train.'$TF'.conservative.train.top5k.narrowPeak0' --outputfile='peaks/top5k/ChIPseq.train.'$TF'.conservative.train.top5k.narrowPeak' --sample_nrows=5000 --seed=2019



echo 9 seqpos_scan_motifs_on_peaks
ct_id_i='train'
echo $ct_id_i
input_pk='peaks/top5k/ChIPseq.'$ct_id_i'.'$TF'.conservative.'$ct_id_i'.top5k.narrowPeak'
$seqpos_env_path'/bin/python2.7' $seqpos_env_path'/bin/MDSeqPos.py' $input_pk $genome \
-m $motif_database_xml -d \
-O './train/'$TF'/motif_scan/cell_line.'$ct_id_i'.tf.'$TF

echo 10 decide_cofactor_motif_set
ct_id_i='train'
echo $ct_id_i
input_file='train/'$TF'/motif_scan/cell_line.'$ct_id_i'.tf.'$TF'/mdseqpos_index.html'
output_file='train/'$TF'/'$TF'_cofactor_motif_list.json'
topn=10
python $scrip_dir/motif_features.py decide_cofactor_motif_set --outfilename=$output_file --n=$topn $input_file


