#!/bin/bash

#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

source ~/.bashrc

cd $1
TF=$2
ct1=$3
ct2=$4

echo $TF
TF_train_num=$(wc -l $TF'_'$ct1'_'$ct2'/peaks/5fc/'$TF'.train_sort_summits.5fc.bed.tmp' |awk '{print $1}')
if (( $TF_train_num != 0 ))
then
	#cat $TF'_'$ct1'_'$ct2'/peaks/5fc/'$TF'.train.topN.label.tsv' | awk -F '\t' -v OFS='\t' '{if ($4=="B") print $1,$2,$3-150}' > $TF'_'$ct1'_'$ct2'/peaks/5fc/'$TF'.train.topN.pk.forIDEAS.bed'
	#bedtools merge -i $TF'_'$ct1'_'$ct2'/peaks/5fc/'$TF'.train.topN.pk.forIDEAS.bed' > $TF'_'$ct1'_'$ct2'/peaks/5fc/'$TF'.train.topN.pk.forIDEAS.M.bed'
	#rm $TF'_'$ct1'_'$ct2'/peaks/5fc/'$TF'.train.topN.pk.forIDEAS.bed'
	#bedtools intersect -a /homes1/gxiang/projects/impute_cistrome/macrophage/VS_cCRE/encode_cCRE.bed -b $TF'_'$ct1'_'$ct2'/peaks/5fc/'$TF'.train.topN.pk.forIDEAS.M.bed' -c > $TF'_'$ct1'_'$ct2'/peaks/5fc/'$TF'.train.topN.pk.forIDEAS.cCRE.bed'
	cut -f4 $TF'_'$ct1'_'$ct2'/peaks/5fc/'$TF'.train.topN.pk.forIDEAS.cCRE.bed' | awk -F '\t' -v OFS='\t' '{if ($1!=0) print 10; else print 0}' > $TF'_'$ct1'_'$ct2'/peaks/5fc/'$TF'.train.topN.pk.forIDEAS.cCRE.txt'
fi

