#!/bin/bash

#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

source ~/.bashrc

cd $1
TF=$2
#TF='FOXA1'

chip_sampling_num=30
bins='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.bins.bed'
pk_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_peaks/'
pksumit_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_summits/'
seqpos_env_path='/liulab/gjxiang/miniconda3/envs/seqpos/'
motif_database_xml='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/HOCOMOCOv11_full_pwm_HUMAN_mono.xml'
scrip_dir='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/scripts/'
TF_info_file='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls'
genome='hg38'
top_pk_n=10000
chrom_used='chr22'

echo 'mk ouput folders'
mkdir -p peaks
mkdir -p peaks/top
mkdir -p peaks/top/$TF
mkdir -p TFs

echo 'get id list for the TF'
cat $TF_info_file | awk -F '\t' -v OFS='\t' -v target_tf=$TF '{if ($3==target_tf) print $1,$36}' > $TF'.chip.id.list.txt'
id_num=$(wc -l $TF'.chip.id.list.txt' | awk -F ' ' '{print $1}')
if (( id_num > $chip_sampling_num ))
then
	sort -k2,2r $TF'.chip.id.list.txt' | head -$chip_sampling_num | cut -f1 > $TF'.chip.id.list.txt1' && mv $TF'.chip.id.list.txt1' $TF'.chip.id.list.txt'
else
	cut -f1 $TF'.chip.id.list.txt' > $TF'.chip.id.list.txt1' && mv $TF'.chip.id.list.txt1' $TF'.chip.id.list.txt'
fi

echo 'get pooled top peaks'
echo 'get pooled top peaks'
tf_id1=$(head -1 $TF'.chip.id.list.txt')
head -$top_pk_n $pksumit_folder$tf_id1'_sort_summits.bed' | cut -f1,2,3 > 'peaks/top/'$TF'.summits.bed.tmp'
head -$top_pk_n  $pk_folder$tf_id1'_sort_peaks.narrowPeak.bed' | cut -f1,2,3 > 'peaks/top/'$TF'.narrowPeak.bed.tmp'
### loop all samples
for tf_id in $(tail -n+2 $TF'.chip.id.list.txt')
do
	echo $tf_id
	head -$top_pk_n $pksumit_folder$tf_id'_sort_summits.bed' | cut -f1,2,3 >> 'peaks/top/'$TF'.summits.bed.tmp'
	head -$top_pk_n $pk_folder$tf_id'_sort_peaks.narrowPeak.bed' | cut -f1,2,3 >> 'peaks/top/'$TF'.narrowPeak.bed.tmp'
done

echo 'sort & merge'
sort -k1,1 -k2,2n 'peaks/top/'$TF'.summits.bed.tmp' > 'peaks/top/'$TF'.summits.bed.tmp.sort'
bedtools merge -i 'peaks/top/'$TF'.summits.bed.tmp.sort' > 'peaks/top/'$TF'.summits.bed.tmp.sort.merge'
sort -k1,1 -k2,2n 'peaks/top/'$TF'.narrowPeak.bed.tmp' > 'peaks/top/'$TF'.narrowPeak.bed.tmp.sort'
bedtools merge -i 'peaks/top/'$TF'.narrowPeak.bed.tmp.sort' > 'peaks/top/'$TF'.narrowPeak.bed.tmp.sort.merge'

echo 'get mat'
bedtools intersect -a $bins -b 'peaks/top/'$TF'.summits.bed.tmp.sort.merge' -wa -u > 'peaks/top/'$TF'.summits.mat.txt'
bedtools intersect -a $bins -b 'peaks/top/'$TF'.narrowPeak.bed.tmp.sort.merge' -wa -u > 'peaks/top/'$TF'.narrowPeak.mat.txt'

echo 'get mat'
cp 'peaks/top/'$TF'.summits.mat.txt' 'peaks/top/'$TF'.summits.mat.bed'
cp 'peaks/top/'$TF'.narrowPeak.mat.txt' 'peaks/top/'$TF'.narrowPeak.mat.bed'
echo 'loop all samples'
for tf_id in $(cat $TF'.chip.id.list.txt')
do
	echo $tf_id
	# summit
	head -$top_pk_n $pksumit_folder$tf_id'_sort_summits.bed' | cut -f1,2,3 > 'peaks/top/'$TF'/'$tf_id'.summits.bed'
	bedtools intersect -a 'peaks/top/'$TF'.summits.mat.bed' -b 'peaks/top/'$TF'/'$tf_id'.summits.bed' -c > 'peaks/top/'$TF'.summits.mat.txt.tmp'
	cut -f4 'peaks/top/'$TF'.summits.mat.txt.tmp' | awk -F '\t' -v OFS='\t' '{if ($1!=0) print 1; else print 0}' > 'peaks/top/'$TF'.summits.mat.txt.tmp1'
	paste 'peaks/top/'$TF'.summits.mat.txt' 'peaks/top/'$TF'.summits.mat.txt.tmp1' > 'peaks/top/'$TF'.summits.mat.txt.tmp'
	mv 'peaks/top/'$TF'.summits.mat.txt.tmp' 'peaks/top/'$TF'.summits.mat.txt'
	## pk
	head -$top_pk_n $pk_folder$tf_id'_sort_peaks.narrowPeak.bed' | cut -f1,2,3 > 'peaks/top/'$TF'/'$tf_id'.narrowPeak.bed'
	bedtools intersect -a 'peaks/top/'$TF'.narrowPeak.mat.bed' -b 'peaks/top/'$TF'/'$tf_id'.narrowPeak.bed' -c > 'peaks/top/'$TF'.narrowPeak.mat.txt.tmp'
	cut -f4 'peaks/top/'$TF'.narrowPeak.mat.txt.tmp' | awk -F '\t' -v OFS='\t' '{if ($1!=0) print 1; else print 0}' > 'peaks/top/'$TF'.narrowPeak.mat.txt.tmp1'
	paste 'peaks/top/'$TF'.narrowPeak.mat.txt' 'peaks/top/'$TF'.narrowPeak.mat.txt.tmp1' > 'peaks/top/'$TF'.narrowPeak.mat.txt.tmp'
	mv 'peaks/top/'$TF'.narrowPeak.mat.txt.tmp' 'peaks/top/'$TF'.narrowPeak.mat.txt'
done

###
echo 'get reporducible pks'
cat 'peaks/top/'$TF'.summits.mat.txt' | awk -F '\t' -v OFS='\t' '{sum=0; for(i=4; i<=NF; i++) sum += $i; if ((NF>=4) && (sum>=2)) print $1,$2,$3; else if (NF==4) print $1,$2,$3}' > 'peaks/top/'$TF'.summits.repro.bed'
cat 'peaks/top/'$TF'.narrowPeak.mat.txt' | awk -F '\t' -v OFS='\t' '{sum=0; for(i=4; i<=NF; i++) sum += $i; if ((NF>=4) && (sum>=2)) print $1,$2,$3; else if (NF==4) print $1,$2,$3}' > 'peaks/top/'$TF'.narrowPeak.repro.bed'

###
echo 'get label vector'
bedtools intersect -a $bins -b 'peaks/top/'$TF'.summits.repro.bed' -c > 'peaks/top/'$TF'.summits.repro.bin.bed'
bedtools intersect -a $bins -b 'peaks/top/'$TF'.narrowPeak.repro.bed' -c > 'peaks/top/'$TF'.narrowPeak.repro.bin.bed'
paste 'peaks/top/'$TF'.summits.repro.bin.bed' 'peaks/top/'$TF'.narrowPeak.repro.bin.bed' | awk -F '\t' -v OFS='\t' '{if ($4!=0) print $1,$2,$3,"B"; else if (($4==0) && ($8==0)) print $1,$2,$3,"U"; else print $1,$2,$3,"A"}' > 'peaks/top/'$TF'.label.tsv'

###
cat 'peaks/top/'$TF'.label.tsv' | awk -F '\t' -v OFS='\t' -v chrom=$chrom_used '{if (($1==chrom) && ($4=="B")) print 1; else if ($1==chrom) print 0}' > 'peaks/top/'$TF'.label.binary.'$chrom_used'.tsv'
chrom_used='chr10'
cat 'peaks/top/'$TF'.label.tsv' | awk -F '\t' -v OFS='\t' -v chrom=$chrom_used '{if (($1==chrom) && ($4=="B")) print 1; else if ($1==chrom) print 0}' > 'peaks/top/'$TF'.label.binary.'$chrom_used'.tsv'


