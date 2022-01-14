#!/bin/bash

#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

source ~/.bashrc

cd $1
TF=$2
#TF='FOXA1'

chip_sampling_num=30
narrowpk_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_peaks/'
pksumit_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_summits/'
seqpos_env_path='/liulab/gjxiang/miniconda3/envs/seqpos/'
motif_database_xml='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/HOCOMOCOv11_full_pwm_HUMAN_mono.xml'
scrip_dir='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/scripts/'
TF_info_file='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls'
genome='hg38'

### get ouput folders
mkdir -p peaks
mkdir -p peaks/top5k
mkdir -p TFs

### get reproducible peaks
#cat 'peaks/top/'$TF'.label.tsv' | awk -F '\t' -v OFS='\t' '{if ($4=="B") print $1,$2,$3}' > 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.rep.bed0'
cat 'peaks/5fc/'$TF'.topN.label.tsv' | awk -F '\t' -v OFS='\t' '{if ($4=="B") print $1,$2,$3}' > 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.rep.bed0'
### sample 5k if more 5k peaks
time python3 $scrip_dir/sample_rows.py --inputfile='peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.rep.bed0' --outputfile='peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.rep.bed' --sample_nrows=5000 --seed=2019

### get top10 motifs
echo 9 seqpos_scan_motifs_on_peaks
ct_id_i='train'
echo $ct_id_i
input_pk='peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.rep.bed'
$seqpos_env_path'/bin/python2.7' $seqpos_env_path'/bin/MDSeqPos.py' $input_pk $genome \
-m $motif_database_xml -d \
-O './TFs/'$TF'/motif_scan/'$TF
###
echo 10 decide_cofactor_motif_set
input_file='TFs/'$TF'/motif_scan/'$TF'/mdseqpos_index.html'
output_file='TFs/'$TF'/'$TF'_cofactor_motif_list.json'
topn=10
python $scrip_dir/motif_features.py decide_cofactor_motif_set --outfilename=$output_file --n=$topn $input_file


