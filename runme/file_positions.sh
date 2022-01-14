############ Input files
### genome info fasta chrom size
genome_info_hg38_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_hg38_config'

### whole genome 200bp slide bin (50 sliding window)
##### human hg38
bin_hg38='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/hg38_regions_all_sorted.bed'
##### mouse mm10
bin_mm10='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/mm10_regions_all_sorted.bed'

### motif h5 files
##### human hg38
motif_hg38_h5_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human/hdf5s/motif'
##### mouse mm10
motif_mm10_h5_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_mouse/hdf5s/motif'

### DNase bam files
##### human hg38
DNase_bam_hg38_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human/bams/with_filter'
##### mouse mm10
DNase_bam_mm10_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_mouse/bams/with_filter'

### DNase 5' counts bigwig files
##### human hg38
DNase_bigwig_hg38_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human/bigwigs'
##### mouse mm10
DNase_bigwig_mm10_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_mouse/bigwigs'

### DNase signal h5 files
##### human hg38
DNase_hg38_h5_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human/hdf5s/DNase'
##### mouse mm10
DNase_mm10_h5_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_mouse/hdf5s/DNase'


### TF json file
TF_json_folder='/liulab/gjxiang/projects/impute_cistrome/get_motif_difscore/TFs'

### TF ChIP-seq peak sumit files
pksumit_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_summits'

### TF ChIP-seq peak files
pk_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_peaks'

### reliable TF cell type list
reliable_TF_ct_list='/liulab/gjxiang/projects/impute_cistrome/get_motif_difscore/reliable_TFs.txt'



############ packages
### Seqpos package
seqpos_package_folder='/liulab/gjxiang/miniconda3/envs/seqpos/'

### giggle package
giggle_package_folder='/liulab/jfan/test/giggle/bin/giggle'



############ Previous files, useful but not necessary for running the pipeline
### motif pwm files
motif_pwm_path='/liulab/gjxiang/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/motifs'

### Cistrome RPs
Cistrome_RPs='/liulab/gjxiang/projects/impute_cistrome/cistrome_TF_peak_RPs'

############
### Jingyu's previous impute_cistrome files
impute_cistrome_jfan_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human'

############
### ENCODE DREAM challenge files
ENCODE_DREAM_folder='/liulab/gjxiang/projects/impute_cistrome/ENCODE_DREAM'



