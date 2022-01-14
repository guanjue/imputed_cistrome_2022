#!/bin/bash

#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

source ~/.bashrc

cd $1
TF=$2

echo start
bins='/homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.chr1_22.bins.bed'
echo get bin count
bedtools intersect -a $bins -b 'peaks/top5k/'$TF'.top5k.narrowPeak.sort.merge.widthfiltered.mat.rep.bed0' -c > 'peaks/top5k/'$TF'.top5k.bincount.bed'

echo cut 4th column
cut -f4 'peaks/top5k/'$TF'.top5k.bincount.bed' > 'peaks/top5k/'$TF'.top5k.bincount.txt'
rm 'peaks/top5k/'$TF'.top5k.bincount.bed'

echo get chr1
#paste $bins 'peaks/top5k/'$TF'.top5k.bincount.txt' | awk -F '\t' -v OFS='\t' '{if ($1=="chr1") print $4}' > 'peaks/top5k/'$TF'.top5k.bincount.chr1.txt'
paste $bins 'peaks/top5k/'$TF'.top5k.bincount.txt' | awk -F '\t' -v OFS='\t' '{if ($1=="chr22") print $4}' > 'peaks/top5k/'$TF'.top5k.bincount.chr22.txt'
paste $bins 'peaks/top5k/'$TF'.top5k.bincount.txt' | awk -F '\t' -v OFS='\t' '{if ($1=="chr10") print $4}' > 'peaks/top5k/'$TF'.top5k.bincount.chr10.txt'


### get ENCODE cCRE
#wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/encode3/ccre/encodeCcreCombined.bb
### bigBed2Bed
#~/softwares/ucsc/bigBedToBed encodeCcreCombined.bb encodeCcreCombined.bed
### get chr1 bins
#bedtools intersect -a /homes1/gxiang/projects/impute_cistrome/cistrome_impute_results_human8/hg38.200_50slide.chr1_22.bins.bed -b encodeCcreCombined.bed -c > bin.cre.count.bed
#cat bin.cre.count.bed | awk -F '\t' -v OFS='\t' '{if ($1=="chr1") print $0}' > bin.cre.count.chr1.bed
#cat bin.cre.count.bed | awk -F '\t' -v OFS='\t' '{if ($1=="chr22") print $0}' > bin.cre.count.chr22.bed
#cat bin.cre.count.bed | awk -F '\t' -v OFS='\t' '{if ($1=="chr10") print $0}' > bin.cre.count.chr10.bed

#time ~/softwares/R/R-4.0.0/bin/Rscript get_motif_score.R

