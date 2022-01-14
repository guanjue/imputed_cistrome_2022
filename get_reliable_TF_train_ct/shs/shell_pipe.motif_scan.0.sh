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

### get id list for the TF
cat $TF_info_file | awk -F '\t' -v OFS='\t' -v target_tf=$TF '{if ($3==target_tf) print $1,$36}' > $TF'.chip.id.list.txt'
id_num=$(wc -l $TF'.chip.id.list.txt' | awk -F ' ' '{print $1}')
if (( id_num > $chip_sampling_num ))
then
	sort -k2,2r $TF'.chip.id.list.txt' | head -$chip_sampling_num | cut -f1 > $TF'.chip.id.list.txt1' && mv $TF'.chip.id.list.txt1' $TF'.chip.id.list.txt'
else
	cut -f1 $TF'.chip.id.list.txt' > $TF'.chip.id.list.txt1' && mv $TF'.chip.id.list.txt1' $TF'.chip.id.list.txt'
fi

### get pooled top peaks
echo 8 prepare top5k.narrowPeak
mkdir -p peaks/top5k/$TF
tf_id1=$(head -1 $TF'.chip.id.list.txt')
head -5000 $narrowpk_folder$tf_id1'_sort_peaks.narrowPeak.bed' | cut -f1,2,3 > peaks/top5k/$TF/$tf_id1'_sort_peaks.narrowPeak.top.bed'
cat peaks/top5k/$TF/$tf_id1'_sort_peaks.narrowPeak.top.bed' > 'peaks/top5k/'$TF'.top5k.narrowPeak'
for tf_id in $(tail -n+2 $TF'.chip.id.list.txt')
do
	echo $tf_id
	head -5000 $narrowpk_folder$tf_id'_sort_peaks.narrowPeak.bed' | cut -f1,2,3 > peaks/top5k/$TF/$tf_id'_sort_peaks.narrowPeak.top.bed'
	cat peaks/top5k/$TF/$tf_id'_sort_peaks.narrowPeak.top.bed' >> 'peaks/top5k/'$TF'.top5k.narrowPeak'
done

### get merge peaks
sort -k1,1 -k2,2n 'peaks/top5k/'$TF'.top5k.narrowPeak' > 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.bed'
bedtools merge -i 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.bed' > 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.bed'
cat 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.bed' | awk -F '\t' -v OFS='\t' '{if (($3-$2)<=1000) print $0}' > 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.bed'

### get peak count mat
cp 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.bed' 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.txt'
for tf_id in $(cat $TF'.chip.id.list.txt')
do
	echo $tf_id
	bedtools intersect -a 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.bed' -b peaks/top5k/$TF/$tf_id'_sort_peaks.narrowPeak.top.bed' -c > 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.bed.tmp'
	cut -f4 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.bed.tmp' > 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.bed.tmp1'
	paste 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.txt' 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.bed.tmp1' > 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.txt1'
	mv 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.txt1' 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.txt'
	rm 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.bed.tmp1' 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.bed.tmp'
done

### get reproducible peaks
#cat 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.txt' | awk -F '\t' -v OFS='\t' '{num0=0; for(i=4; i<=NF; i++) if ($i==0) {num0++}; if ((1-num0/(NF-3))>0.5) print $1,$2,$3}' > 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.rep.bed0'
cat 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.txt' | awk -F '\t' -v OFS='\t' '{num0=0; for(i=4; i<=NF; i++) if ($i==0) {num0++}; if ((NF>4) && ((NF-3-num0)>=2)) print $1,$2,$3; else if (NF==4) print $1,$2,$3}' > 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.rep.bed0'
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


