#!/usr/bin/env python3
import os
import yaml

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

config = dict(DNase_bam_file_path='/liulab/jfan/projects/impute_cistrome/cistrome/DNase_bam_mouse',
              cell_types=['45133', '78618', '78205','65405_ATAC', '48713'],
              reference_cell_type='65405_ATAC',
              tf_names=[],
              chrom_all=list(map(lambda x: "chr" + str(x), list(range(1, 20)) + ["X"])),
              genome="mm10",
              chrom_size_file='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/mm10.chrom'
                              '.sizes',
              regions_all_file='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config'
                               '/mm10_regions_all_sorted.bed',
              batch=10000000,
              genome_sequence_fa='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/mm10'
                                 '.fa',
              motif_pwm_path='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/motifs/',
              peak_file_path='/liulab/jfan/projects/impute_cistrome/cistrome/peak_mouse',
              motif_database_xml='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config'
                                 '/HOCOMOCOv11_full_pwm_HUMAN_mono.xml',
              chrom_sets=dict(
                  chrom_setA=['chr2', 'chr4', 'chr6', 'chr7', 'chr12', 'chr15', 'chr16', 'chr17',
                              'chrX'],
                  chrom_setB=['chr1', 'chr3', 'chr5', 'chr9', 'chr10', 'chr11', 'chr13', 'chr14', 'chr18'],
                  chrom_set_test=["chr8", "chr19"]),
              training_cell_types_regions_label_path="/liulab/jfan/projects/impute_cistrome/cistrome_impute_results_mouse/label"
                                                     "/train",
              training_cell_types_regions_label_name="train.labels.tsv",
              test_cell_types_regions_label_path="/liulab/jfan/projects/impute_cistrome/cistrome_impute_results_mouse/label"
                                                 "/test",
              test_cell_types_regions_label_name="train.labels.tsv",
              seqpos_env_path="/liulab/jfan/miniconda3/envs/seqpos",
              region_topic_model_h5='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/'
                                                'hg38_regions_topic_all_chroms.h5'
              )

print(os.getcwd())
with open('config.yaml', 'w') as outfile:
    yaml.dump(config, outfile, Dumper=Dumper)

with open('config.yaml', 'r') as infile:
    print(yaml.load(infile, Loader=Loader))

cluster_config = dict(
    __default__=dict(
        job_name="default",
        time="1:00:00",
        n=1,
        mem_mb=1000,
        partition="defq"),
    bam_sort=dict(
                job_name="bam_sort",
                time="2:00:00",
                n=10,
                mem_mb=10000,
                partition="defq"),
    combine_bam_file_for_same_cell_line=dict(
        job_name="combine_bam_file_for_same_cell_line",
        time="6:00:00",
        n=4,
        mem_mb=1000,
        partition="defq"),
    generate_bw_5_end_1_bp=dict(
        job_name="generate_bw_5_end_1_bp",
        time="8:00:00",
        n=8,
        mem_mb=6000,
        partition="defq"),
    generate_DNase_features_from_bw_to_h5=dict(
        job_name="generate_DNase_features_from_bw_to_h5",
        time="8:00:00",
        n=8,
        mem_mb=20000,
        partition="defq"),
    combine_dnase_features_h5_for_all_cell_types=dict(
        job_name="combine_dnase_features_h5_for_all_cell_types",
        time="4:00:00",
        n=8,
        mem_mb=80000,
        partition="defq"),
    generate_quantile_transformer=dict(
        job_name="generate_quantile_transformer",
        time="2:00:00",
        n=2,
        mem_mb=10000,
        partition="defq"),
    generate_dnase_feature_median=dict(
        job_name="generate_dnase_feature_median",
        time="1:00:00",
        n=2,
        mem_mb=80000,
        partition="defq"),
    prepare_motif_top4_feature=dict(
        job_name="prepare_motif_top4_feature",
        time="3:00:00",
        n=32,
        mem_mb=120000,
        partition="defq"),
    seqpos_scan_motifs_on_peaks=dict(
        job_name="seqpos_scan_motifs_on_peaks",
        time="10:00:00",
        n=1,
        mem_mb=10000,
        partition="defq"),
    prepare_motif_h5_data=dict(
        job_name="prepare_motif_h5_data",
        time="00:30:00",
        n=1,
        mem_mb=120000,
        partition="defq"),
    prepare_lightgbm_binary_data_motif_feature_reference=dict(
        job_name="prepare_lightgbm_binary_data_motif_feature_reference",
        time="1:00:00",
        n=1,
        mem_mb=240000,
        partition="defq"),
    prepare_lightgbm_binary_data_motif_feature_subset=dict(
        job_name="prepare_lightgbm_binary_data_motif_feature_subset",
        time="4:00:00",
        n=1,
        mem_mb=240000,
        partition="defq"),
    prepare_lightgbm_binary_data_dnase_feature_reference=dict(
        job_name="prepare_lightgbm_binary_data_dnase_feature_reference",
        time="4:00:00",
        n=1,
        mem_mb=20000,
        partition="defq"),
    prepare_lightgbm_binary_data_dnase_feature_all=dict(
        job_name="prepare_lightgbm_binary_data_dnase_feature_all",
        time="4:00:00",
        n=1,
        mem_mb=40000,
        partition="defq"),
    merge_lightgbm_binary_data=dict(
        job_name="merge_lightgbm_binary_data",
        time="4:00:00",
        n=1,
        mem_mb=120000,
        partition="defq"),
    train_models=dict(
        job_name="train_models",
        time="4:00:00",
        n=16,
        mem_mb=240000,
        partition="defq"),
    train_models_hyperopt=dict(
            job_name="train_models",
            time="16:00:00",
            n=32,
            mem_mb=240000,
            partition="defq"),
    make_prediction_chrom=dict(
        job_name="make_prediction_chrom",
        time="2:00:00",
        n=1,
        mem_mb=60000,
        partition="defq"),
    make_prediction_hyperopt=dict(
            job_name="make_prediction_chrom",
            time="2:00:00",
            n=1,
            mem_mb=60000,
            partition="defq"),
    evaluation=dict(
        job_name="evaluation",
        time="2:00:00",
        n=1,
        mem_mb=60000,
        partition="defq"),
    evaluation_hyperopt=dict(
            job_name="evaluation",
            time="2:00:00",
            n=1,
            mem_mb=60000,
            partition="defq"),
    evaluation_leafs=dict(
            job_name="evaluation",
            time="6:00:00",
            n=1,
            mem_mb=120000,
            partition="defq"),
    make_prediction_chrom_leaf=dict(
        job_name="make_prediction_chrom",
        time="12:00:00",
        n=1,
        mem_mb=120000,
        partition="defq"),
    make_prediction_autoencoder=dict(
        job_name="evaluation",
        time="2:00:00",
        n=1,
        mem_mb=60000,
        partition="defq"),
)

with open('cluster.yaml', 'w') as outfile:
    yaml.dump(cluster_config, outfile, Dumper=Dumper)
