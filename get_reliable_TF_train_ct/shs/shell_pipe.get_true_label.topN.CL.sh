#!/bin/bash

#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

source ~/.bashrc

cd $1
TF_CL_CT_TIS=$2
#

#TF='FOXA1'
#TF='CTCF'

pksumit_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_summits/'
pk_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_peaks/'
header_file='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/train_test_header_file.txt'


bins='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.bins.bed'


chip_sampling_num=10
bins='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.bins.bed'
pk_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_peaks/'
pksumit_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_summits/'
seqpos_env_path='/liulab/gjxiang/miniconda3/envs/seqpos/'
motif_database_xml='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/HOCOMOCOv11_full_pwm_HUMAN_mono.xml'
scrip_dir='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/scripts/'
TF_CL_CT_TIS_info_file='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls'
genome='hg38'
top_pk_n=5000
chrom_used='chr22'

echo 'mk ouput folders'
mkdir -p peaks
mkdir -p peaks/top
mkdir -p peaks/top/$TF_CL_CT_TIS
mkdir -p TFs

echo 'get id list for the TF'
cat $TF_CL_CT_TIS_info_file | awk -F '\t' -v OFS='\t' -v target_tf=$TF_CL_CT_TIS '{if ($3"_"$9"_"$10"_"$11==target_tf) print $1,$36}' > $TF_CL_CT_TIS'.chip.id.list.txt'

id_num=$(wc -l $TF_CL_CT_TIS'.chip.id.list.txt' | awk -F ' ' '{print $1}')

if (( id_num > $chip_sampling_num ))
then
        sort -k2,2r $TF_CL_CT_TIS'.chip.id.list.txt' | head -$chip_sampling_num | cut -f1 > $TF_CL_CT_TIS'.chip.id.list.txt1' && mv $TF_CL_CT_TIS'.chip.id.list.txt1' $TF_CL_CT_TIS'.chip.id.list.txt'
elif (( id_num == 1 ))
then
	cut -f1 $TF_CL_CT_TIS'.chip.id.list.txt' > $TF_CL_CT_TIS'.chip.id.list.txt1'
	cut -f1 $TF_CL_CT_TIS'.chip.id.list.txt' >> $TF_CL_CT_TIS'.chip.id.list.txt1' && mv $TF_CL_CT_TIS'.chip.id.list.txt1' $TF_CL_CT_TIS'.chip.id.list.txt'
else
        cut -f1 $TF_CL_CT_TIS'.chip.id.list.txt' > $TF_CL_CT_TIS'.chip.id.list.txt1' && mv $TF_CL_CT_TIS'.chip.id.list.txt1' $TF_CL_CT_TIS'.chip.id.list.txt'
fi



echo ################
echo training set topN set
tf_id1=$(head -1 $TF_CL_CT_TIS'.chip.id.list.txt')
head -$top_pk_n $pksumit_folder$tf_id1'_sort_summits.bed' | cut -f1,2,3 > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.5fc.bed.tmp'
head -$top_pk_n  $pk_folder$tf_id1'_sort_peaks.narrowPeak.bed' | cut -f1,2,3 > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.5fc.bed.tmp'
# summit mat


bedtools intersect -a $bins -b 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.5fc.bed.tmp' -c > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.txt.tmp'
cut -f4 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.txt.tmp' | awk '{if ($1!=0) print 1; else print 0}' > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.txt.tmp1'
mv 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.txt.tmp1' 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.topN.txt'
rm 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.txt.tmp'
# peak mat
bedtools intersect -a $bins -b 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.5fc.bed.tmp' -c > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.txt.tmp'
cut -f4 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.txt.tmp' | awk '{if ($1!=0) print 1; else print 0}' > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.txt.tmp1'
mv 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.txt.tmp1' 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.topN.txt'
rm 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.txt.tmp'


for tf_id in $(tail -n+2 $TF_CL_CT_TIS'.chip.id.list.txt')
do
	echo $tf_id
	head -$top_pk_n $pksumit_folder$tf_id'_sort_summits.bed' | cut -f1,2,3 > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.5fc.bed.tmp'
	head -$top_pk_n $pk_folder$tf_id'_sort_peaks.narrowPeak.bed' | cut -f1,2,3 > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.5fc.bed.tmp'
	#
	bedtools intersect -a $bins -b 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.5fc.bed.tmp' -c > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.txt.tmp'
	cut -f4 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.txt.tmp' | awk '{if ($1!=0) print 1; else print 0}' > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.txt.tmp1'
	paste 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.topN.txt' 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.txt.tmp1' > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.txt.tmp'
	mv 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.txt.tmp' 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.topN.txt'
	rm 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.txt.tmp1'
	#
	bedtools intersect -a $bins -b 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.5fc.bed.tmp' -c > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.txt.tmp'
	cut -f4 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.txt.tmp' | awk '{if ($1!=0) print 1; else print 0}' > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.txt.tmp1'
	paste 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.topN.txt' 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.txt.tmp1' > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.txt.tmp'
	mv 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.txt.tmp' 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.topN.txt'
	rm 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.txt.tmp1'
done

### get intersect submit & pk bins
cat 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.topN.txt' | awk -F '\t' -v OFS='\t' '{sum=0; for(i=1; i<=NF; i++) sum += $i; print sum}' > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.topN.sum.txt'
cat 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.topN.txt' | awk -F '\t' -v OFS='\t' '{sum=0; for(i=1; i<=NF; i++) sum += $i; print sum}' > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.topN.sum.txt'

paste $bins 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.mat.topN.sum.txt' | awk -F '\t' -v OFS='\t' '{if ($4>=2) print $1,$2,$3,1; else print $1,$2,$3,0}' > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.5fc.pkbin.topN.txt'
paste $bins 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.mat.topN.sum.txt' | awk -F '\t' -v OFS='\t' '{if ($4>=2) print $1,$2,$3,1; else print $1,$2,$3,0}' > 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.5fc.pkbin.topN.txt'

### get label tsv
cp $header_file 'peaks/5fc/'$TF_CL_CT_TIS'.topN.label.tsv'
paste 'peaks/5fc/'$TF_CL_CT_TIS'.sort_summits.5fc.pkbin.topN.txt' 'peaks/5fc/'$TF_CL_CT_TIS'.sort_pk.5fc.pkbin.topN.txt' \
| awk -F '\t' -v OFS='\t' '{if ($4!=0) print $1,$2,$3,"B"; else if ($4==0 && $8!=0) print $1,$2,$3,"A"; else print $1,$2,$3,"U"}' >> 'peaks/5fc/'$TF_CL_CT_TIS'.topN.label.tsv'

chrom_used='chr22'
cat 'peaks/5fc/'$TF_CL_CT_TIS'.topN.label.tsv' | awk -F '\t' -v OFS='\t' -v chrom=$chrom_used '{if (($1==chrom) && ($4=="B")) print 1; else if ($1==chrom) print 0}' > 'peaks/top/'$TF_CL_CT_TIS'.label.binary.'$chrom_used'.tsv'
chrom_used='chr10'
cat 'peaks/5fc/'$TF_CL_CT_TIS'.topN.label.tsv' | awk -F '\t' -v OFS='\t' -v chrom=$chrom_used '{if (($1==chrom) && ($4=="B")) print 1; else if ($1==chrom) print 0}' > 'peaks/top/'$TF_CL_CT_TIS'.label.binary.'$chrom_used'.tsv'
#cat 'peaks/5fc/'$TF_CL_CT_TIS'.topN.label.tsv' | awk -F '\t' -v OFS='\t' -v chrom=$chrom_used '{if (($1==chrom) && ($4=="B")) print 2; if (($1==chrom) && ($4=="A")) print 1; else if ($1==chrom) print 0}' > 'peaks/top/'$TF_CL_CT_TIS'.label.trinary.'$chrom_used'.tsv'






