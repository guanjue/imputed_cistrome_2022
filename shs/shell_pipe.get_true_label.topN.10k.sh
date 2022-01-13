TF=$1
thread_num=4
#TF='FOXA1'
#TF='CTCF'
declare -a chr_all=("chr1" "chr2")

pksumit_folder='/liulab/jfan/projects/impute_cistrome/cistrome/human_TF_summits/'
pk_folder='/liulab/jfan/projects/impute_cistrome/cistrome/human_TF_peaks/'
train_test_header_file='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/train_test_header_file.txt'

train_cl_ct_tis=$2

bins='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.bins.bed'


pksumit_folder=$3
pk_folder=$4
train_test_header_file=$5
bins=$6

mkdir -p peaks/5fc

echo 8 get bins intersect with submits

mkdir -p peaks
mkdir -p peaks/5fc

echo ################
echo training set topN set
top_pk_n=10000
tf_id1=$(head -1 $train_cl_ct_tis'.'$TF'.train.id.list.txt')
head -$top_pk_n $pksumit_folder$tf_id1'_sort_summits.bed' | cut -f1,2,3 > 'peaks/5fc/'$TF'.train_sort_summits.5fc.bed.tmp'
head -$top_pk_n  $pk_folder$tf_id1'_sort_peaks.narrowPeak.bed' | cut -f1,2,3 > 'peaks/5fc/'$TF'.train_sort_pk.5fc.bed.tmp'
# summit mat
bedtools intersect -a $bins -b 'peaks/5fc/'$TF'.train_sort_summits.5fc.bed.tmp' -c > 'peaks/5fc/'$TF'.train_sort_summits.mat.txt.tmp'
cut -f4 'peaks/5fc/'$TF'.train_sort_summits.mat.txt.tmp' | awk '{if ($1!=0) print 1; else print 0}' > 'peaks/5fc/'$TF'.train_sort_summits.mat.txt.tmp1'
mv 'peaks/5fc/'$TF'.train_sort_summits.mat.txt.tmp1' 'peaks/5fc/'$TF'.train_sort_summits.mat.topN.txt'
rm 'peaks/5fc/'$TF'.train_sort_summits.mat.txt.tmp'

# peak mat
bedtools intersect -a $bins -b 'peaks/5fc/'$TF'.train_sort_pk.5fc.bed.tmp' -c > 'peaks/5fc/'$TF'.train_sort_pk.mat.txt.tmp'
cut -f4 'peaks/5fc/'$TF'.train_sort_pk.mat.txt.tmp' | awk '{if ($1!=0) print 1; else print 0}' > 'peaks/5fc/'$TF'.train_sort_pk.mat.txt.tmp1'
mv 'peaks/5fc/'$TF'.train_sort_pk.mat.txt.tmp1' 'peaks/5fc/'$TF'.train_sort_pk.mat.topN.txt'
rm 'peaks/5fc/'$TF'.train_sort_pk.mat.txt.tmp'


for tf_id in $(tail -n+2 $train_cl_ct_tis'.'$TF'.train.id.list.txt')
do
	echo $tf_id
	head -$top_pk_n $pksumit_folder$tf_id'_sort_summits.bed' | cut -f1,2,3 > 'peaks/5fc/'$TF'.train_sort_summits.5fc.bed.tmp'
	head -$top_pk_n $pk_folder$tf_id'_sort_peaks.narrowPeak.bed' | cut -f1,2,3 > 'peaks/5fc/'$TF'.train_sort_pk.5fc.bed.tmp'
	#
	bedtools intersect -a $bins -b 'peaks/5fc/'$TF'.train_sort_summits.5fc.bed.tmp' -c > 'peaks/5fc/'$TF'.train_sort_summits.mat.txt.tmp'
	cut -f4 'peaks/5fc/'$TF'.train_sort_summits.mat.txt.tmp' | awk '{if ($1!=0) print 1; else print 0}' > 'peaks/5fc/'$TF'.train_sort_summits.mat.txt.tmp1'
	paste 'peaks/5fc/'$TF'.train_sort_summits.mat.topN.txt' 'peaks/5fc/'$TF'.train_sort_summits.mat.txt.tmp1' > 'peaks/5fc/'$TF'.train_sort_summits.mat.txt.tmp'
	mv 'peaks/5fc/'$TF'.train_sort_summits.mat.txt.tmp' 'peaks/5fc/'$TF'.train_sort_summits.mat.topN.txt'
	rm 'peaks/5fc/'$TF'.train_sort_summits.mat.txt.tmp1'
	#
	bedtools intersect -a $bins -b 'peaks/5fc/'$TF'.train_sort_pk.5fc.bed.tmp' -c > 'peaks/5fc/'$TF'.train_sort_pk.mat.txt.tmp'
	cut -f4 'peaks/5fc/'$TF'.train_sort_pk.mat.txt.tmp' | awk '{if ($1!=0) print 1; else print 0}' > 'peaks/5fc/'$TF'.train_sort_pk.mat.txt.tmp1'
	paste 'peaks/5fc/'$TF'.train_sort_pk.mat.topN.txt' 'peaks/5fc/'$TF'.train_sort_pk.mat.txt.tmp1' > 'peaks/5fc/'$TF'.train_sort_pk.mat.txt.tmp'
	mv 'peaks/5fc/'$TF'.train_sort_pk.mat.txt.tmp' 'peaks/5fc/'$TF'.train_sort_pk.mat.topN.txt'
	rm 'peaks/5fc/'$TF'.train_sort_pk.mat.txt.tmp1'
done

### get intersect submit & pk bins
cat 'peaks/5fc/'$TF'.train_sort_summits.mat.topN.txt' | awk -F '\t' -v OFS='\t' '{sum=0; for(i=1; i<=NF; i++) sum += $i; print sum}' > 'peaks/5fc/'$TF'.train_sort_summits.mat.topN.sum.txt'
cat 'peaks/5fc/'$TF'.train_sort_pk.mat.topN.txt' | awk -F '\t' -v OFS='\t' '{sum=0; for(i=1; i<=NF; i++) sum += $i; print sum}' > 'peaks/5fc/'$TF'.train_sort_pk.mat.topN.sum.txt'

paste $bins 'peaks/5fc/'$TF'.train_sort_summits.mat.topN.sum.txt' | awk -F '\t' -v OFS='\t' '{if ($4>=2) print $1,$2,$3,1; else print $1,$2,$3,0}' > 'peaks/5fc/'$TF'.train_sort_summits.5fc.pkbin.topN.txt'
paste $bins 'peaks/5fc/'$TF'.train_sort_pk.mat.topN.sum.txt' | awk -F '\t' -v OFS='\t' '{if ($4>=2) print $1,$2,$3,1; else print $1,$2,$3,0}' > 'peaks/5fc/'$TF'.train_sort_pk.5fc.pkbin.topN.txt'

### get label tsv
cp $train_test_header_file 'peaks/5fc/'$TF'.train.topN.label.tsv'
paste 'peaks/5fc/'$TF'.train_sort_summits.5fc.pkbin.topN.txt' 'peaks/5fc/'$TF'.train_sort_pk.5fc.pkbin.topN.txt' \
| awk -F '\t' -v OFS='\t' '{if ($4!=0) print $1,$2,$3,"B"; else if ($4==0 && $8!=0) print $1,$2,$3,"A"; else print $1,$2,$3,"U"}' >> 'peaks/5fc/'$TF'.train.topN.label.tsv'






