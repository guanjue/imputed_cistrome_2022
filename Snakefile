import glob
import os

import numpy as np
import pandas as pd
import h5py

configfile: "config.yaml"



localrules: all, top_5k_peaks, decide_cofactor_motif_set, motif_h5_data_ready, copy_binary_files_reference,
          copy_binary_files_motif_reference, TF_ChIP_sample_CA_cluster_matching, combine_TF_ChIP_sample_CA_cluster_results,
            train_models_all

ruleorder: prepare_lightgbm_binary_data_dnase_feature_all > copy_binary_files_reference
ruleorder: prepare_lightgbm_binary_data_dnase_feature_all_autoencoder > copy_binary_files_reference_autoencoder

def cell_line_for_evaluation(tf_name):
    df_all_regions_label = pd.read_csv(
        "%s/%s.%s" % (config['test_cell_types_regions_label_path'], tf_name,
                      config['test_cell_types_regions_label_name']),
        sep="\t", header=0, nrows=10)
    return list(df_all_regions_label.columns[3:])


def cell_line_for_train(tf_name):
    if config['CA_cluster_mode']:
        label_file = checkpoints.prepare_label_file.get(tf_name=tf_name).output[0]
    else:
        label_file = '%s/%s.%s' % (config['training_cell_types_regions_label_path'], tf_name,
                                   config['training_cell_types_regions_label_name'])
    df_all_regions_label = pd.read_csv(
        label_file,
        sep="\t", header=0, nrows=10)
    return list(df_all_regions_label.columns[3:])


CELL_LINE = []
TF_NAME = []
CHROM_SET_NAME = []

def tf_name_to_train(list_tf_names):
    CELL_LINE = []
    TF_NAME = []
    CHROM_SET_NAME = []
    for tf_name in list_tf_names:
        cell_lines = cell_line_for_train(tf_name)
        CELL_LINE += cell_lines * 2
        TF_NAME += [tf_name] * (2 * len(cell_lines))
        CHROM_SET_NAME += ['chrom_setA'] * len(cell_lines) + ['chrom_setB'] * len(cell_lines)
    return CELL_LINE, TF_NAME, CHROM_SET_NAME

# train models for a list of TFs
# example: list_tf_names = ['CTCF', 'MYB']
list_tf_names = []

CELL_LINE, TF_NAME, CHROM_SET_NAME = tf_name_to_train(list_tf_names)



rule all:
    input:
        # "train/data/quantile_transformer/62395_variable_bp_quantile_map.pkl"
        #  expand(
        #      "train/NR3C1/predictions/NR3C1.48.{chrom}_preds.h5",
        #      chrom=config['chrom_all']
        #  )
         # 'train/TCF7/predictions/TCF7.44921.chr21_preds.h5'
         expand(
             "train/{tf_name}/models/{tf_name}.{training_cell_line}.{training_chrom_set_name}_model.pkl",
             zip,
             tf_name=TF_NAME,
             training_cell_line=CELL_LINE,
             training_chrom_set_name=CHROM_SET_NAME
         )
         # expand("train/{tf_name}/evaluations/{tf_name}.{cell_line}_performance.txt",
         #        zip,
         #        tf_name=TF_NAME,
         #        cell_line=CELL_LINE
         #        )
         # expand(
         #     "train/BCL6/models/BCL6.{training_cell_line}.{training_chrom_set_name}_model.pkl",
         #     training_cell_line=[78454,3139,3203,3142,44907,40606],
         #     training_chrom_set_name=['chrom_setA']*6+['chrom_setB']*6)
         # expand("hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
         #        chrom=config['chrom_all']),
         # expand("train/{tf_name}/selected_motif_hdf5/{chrom}_motif_features_lightGBM.h5", chrom = config['chrom_all'],
         #        tf_name=["NFATC3"])
         # expand("hdf5s/motif/{chrom}_motifs_top4_scores.h5", chrom = config['chrom_all'] )
         # expand("train/{tf_name}/evaluations/{tf_name}.{cell_line}_performance.hyperopt.txt",
         #        zip,
         #        tf_name=TF_NAME,
         #        cell_line=CELL_LINE
         #        )
        # expand("train/{tf_name}/evaluations/{tf_name}.{cell_line}_confusion_matrix.txt",
        #         zip,
        #         tf_name=TF_NAME,
        #         cell_line=CELL_LINE
        #         )
        # expand("train/{tf_name}/evaluations/{tf_name}.{cell_line}_models_fro_norm.txt",
        #         zip,
        #         tf_name=TF_NAME,
        #         cell_line=CELL_LINE
        #         )
         # expand("train/{tf_name}/evaluations/{tf_name}.{cell_line}_performance.autoencoder.txt",
         #        zip,
         #        tf_name=TF_NAME,
         #        cell_line=CELL_LINE
         #        )
         # expand("train/{tf_name}/evaluations/{tf_name}.{cell_line}_performance.txt",
         #        zip,
         #        tf_name=TF_NAME,
         #        cell_line=CELL_LINE
         #        )
         # "peaks/top5k/ChIPseq.HepG2.FOXA1.conservative.train.top5k.narrowPeak"
         # expand("train/quantile_transformer/{cell_line}_variable_bp_quantile_map.pkl",
         #        cell_line=config['cell_types'])
         # expand('hdf5s/DNase/median/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_median.h5',
         #        chrom=config['chrom_all'])
         # 'train/EGR1_cofactor_motif_list.json'
         # "peaks/cell_line.HepG2.tf.FOXA2/mdseqpos_index.html"
         # expand("hdf5s/DNASE_bam_5_mer_50bp_all_samples_lightGBM_{chrom}_all_cell_types.h5", chrom=config['chrom_all'])
         # expand("peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/{sample_id}.giggle.txt",
         #        zip,
         #        tf_name=["RUNX1"]*len(all_TF_ChIP_sample_id('RUNX1')),
         #        sample_id=all_TF_ChIP_sample_id('RUNX1')
         #        )
         # "peaks/TF_ChIP_sample_CA_cluster_matching/RUNX1/RUNX1_ChIP_sample_CA_cluster_assigned.xls"
         # "peaks/TF_ChIP_sample_CA_cluster_matching/GATA3/GATA3_ChIP_sample_CA_cluster_assigned.xls"
        # "{training_cell_types_regions_label_path}/{tf_name}.{training_cell_types_regions_label_name}".format(
        #           training_cell_types_regions_label_path = config['training_cell_types_regions_label_path'],
        #           training_cell_types_regions_label_name = config['training_cell_types_regions_label_name'],
        #           tf_name = 'TCF12'
        # )
        # expand("bams/with_filter/DNASE.{cell_line}.merge.filter.sorted.bam.bai",
        #        cell_line=[40595, 6712, 6715, 6716, 6717, 6718, 40574, 40575, 40562, 40594, 40585, 40566, 35745, 40582, 40586, 40560, 40593, 40584, 40564, 40568, 40577, 40561, 40565, 40569, 40563, 40559, 40428, 40426, 40416, 40430, 40420, 40421, 40417, 40415, 40410, 40423, 40411, 40422, 40427, 40418, 40429, 40401, 40405, 40406, 40419, 40412, 6721, 6722, 6724, 6723, 6719, 6720, 40408, 45067, 40601, 40604, 39138, 35746, 5217, 5216, 5215, 35742, 35747, 35744, 40882, 40875, 40879, 40880, 40877, 40878, 40876, 40881, 33603, 33612, 33596, 44961]
        #
        #        )
        # expand("train/{tf_name}/models/{tf_name}.models.done",
        #        # tf_name=['MYB','TAL1','GATA3','RUNX1','CEBPA']
        #        tf_name=['ESR1','NR3C1','NR4A1','RXRA','NR1H2']
        #        )
        # expand("peaks/all_TF_samples_CA_cluster_matching/{sample_id}.giggle.txt",
        #        sample_id=[i for i in list(map(lambda x: x.split('_')[0], os.listdir("/liulab/jfan/projects/impute_cistrome/cistrome/human_TF_peaks/")))
        #                   if i != "cistromeDB"],
        #        )
        # "train/TCF12/evaluations/TCF12.72_performance.txt",
        # "train/MYB/evaluations/MYB.47_performance.txt",
        # "train/RUNX1/evaluations/RUNX1.47_performance.txt",
        # "train/GATA3/evaluations/GATA3.47_performance.txt",
        # "train/TAL1/evaluations/TAL1.47_performance.txt",
         # "train/CEBPA/evaluations/CEBPA.15_performance.txt",
        # "train/TCF7/evaluations/TCF7.61485_performance.txt",
        # expand("/liulab/jfan/projects/impute_cistrome/cistrome/DNase_mapping/human_DNase_bam/{sample_id}.bam",
        #         sample_id=[40595, 6712, 6715, 6716, 6717, 6718, 40574, 40575, 40562, 40594, 40585, 40566, 35745, 40582, 40586, 40560, 40593, 40584, 40564, 40568, 40577, 40561, 40565, 40569, 40563, 40559, 40428, 40426, 40416, 40430, 40420, 40421, 40417, 40415, 40410, 40423, 40411, 40422, 40427, 40418, 40429, 40401, 40405, 40406, 40419, 40412, 6721, 6722, 6724, 6723, 6719, 6720, 40408, 45067, 40601, 40604, 39138, 35746, 5217, 5216, 5215, 35742, 35747, 35744, 40882, 40875, 40879, 40880, 40877, 40878, 40876, 40881, 33603, 33612, 33596, 44961]
        #        )



rule TF_ChIP_sample_CA_cluster_matching_all:
    input:
     # lambda wildcards: "{peak_path}/ChIPseq.{cell_line}.{tf_name}.conservative.train.narrowPeak.gz".format(
         lambda wildcards: "{peak_path}/{sample_id}_sort_peaks.narrowPeak.bed".format(
             peak_path=config['cistromeDB_narrow_peak_path'],
             sample_id=wildcards.sample_id)
    output:
          "peaks/all_TF_samples_CA_cluster_matching/{sample_id}.giggle.txt"
    # group: "peak_operation"
    resources:
             mem_mb=1000
    params:
          runtime="1h",
          giggle_path=config['giggle_path'],
          giggle_index_path=config['giggle_index_path']
    shell:
        '''
        cut -f 1,2,3 {input} | head -n 10000 > peaks/all_TF_samples_CA_cluster_matching/{wildcards.sample_id}.bed || true
        bedtools sort -i peaks/all_TF_samples_CA_cluster_matching/{wildcards.sample_id}.bed > \
            peaks/all_TF_samples_CA_cluster_matching/{wildcards.sample_id}.sorted.bed
        bgzip peaks/all_TF_samples_CA_cluster_matching/{wildcards.sample_id}.sorted.bed -c > \
            peaks/all_TF_samples_CA_cluster_matching/{wildcards.sample_id}.sorted.bed.gz
        {params.giggle_path} search -i {params.giggle_index_path} \
            -q peaks/all_TF_samples_CA_cluster_matching/{wildcards.sample_id}.sorted.bed.gz \
            -s > {output}
        '''





rule bwa_mapping:
    input:
         "/liulab/jfan/projects/impute_cistrome/cistrome/DNase_mapping/human_DNase_fastq/{sample_id}.fastq.gz"
    output:
          "/liulab/jfan/projects/impute_cistrome/cistrome/DNase_mapping/human_DNase_bam/{sample_id}.bam"
    log:
       "logs/bwa_mapping/{sample_id}.log"
    benchmark:
             "benchmarks/bwa_mapping/{sample_id}.tsv"
    params:
          runtime="6h"
    resources:
             mem_mb=12000
    # group: "bam_file_operation"
    threads: 16
    shell:
         # "(bwa mem -t 16 -p /liulab/jfan/projects/impute_cistrome/cistrome/bwa_index/hg38/bwa_indices/hg38/hg38.fa "
         # "{input} | samtools view  -Sb -o {output} - ) > {log} 2>&1 "
         "(bwa aln -q 5 -l 32 -k 2 -t 16 /liulab/jfan/projects/impute_cistrome/cistrome/bwa_index/hg38/bwa_indices/hg38/hg38.fa "
         "{input} | bwa samse /liulab/jfan/projects/impute_cistrome/cistrome/bwa_index/hg38/bwa_indices/hg38/hg38.fa - "
         "{input} | samtools view -Sb -o {output} - ) > {log} 2>&1 "




rule combine_bam_file_for_same_cell_line:
    input:
         lambda wildcards: glob.glob(
             "{bam_path}/DNASE.{cell_line}.biorep*.techrep*.bam".format(bam_path=config['bam_file_path'],
                                                                        cell_line=wildcards.cell_line))
    output:
          "bams/no_filter/DNASE.{cell_line}.merge.bam"
    log:
       "logs/combine_bam_file_for_same_cell_line/{cell_line}.log"
    benchmark:
             "benchmarks/combine_bam_file_for_same_cell_line/{cell_line}.tsv"
    params:
          runtime="6h"
    resources:
             mem_mb=1000
    # group: "bam_file_operation"
    threads: 4
    shell:
         "(samtools merge {output} {input} --threads $(({threads}-1)) -f) > {log} 2>&1 "

rule bam_filter:
    input:
         # "bams/no_filter/DNASE.{cell_line}.merge.bam"
         # "/liulab/jfan/projects/impute_cistrome/cistrome/DNase_bam/all_bam/{cell_line}.bam"
          lambda wildcards: "{bam_path}/{cell_line}.bam".format(bam_path=config['bam_file_path'],cell_line=wildcards.cell_line)
    output:
          temp("bams/with_filter/DNASE.{cell_line}.merge.filter.bam")
    log:
       "logs/bam_filter/{cell_line}.log"
    benchmark:
             "benchmarks/bam_filter/{cell_line}.tsv"
    params:
          runtime="6h"
    resources:
             mem_mb=1000
    # group: "bam_file_operation"
    threads: 4
    shell:
         "samtools view -F 1804 -q 30 --threads {threads} -b {input} -o {output} > {log} 2>&1 "


rule bam_sort:
    input:
         "bams/with_filter/DNASE.{cell_line}.merge.filter.bam"
    output:
          "bams/with_filter/DNASE.{cell_line}.merge.filter.sorted.bam"
    log:
       "logs/bam_sort/{cell_line}.log"
    benchmark:
             "benchmarks/bam_sort/{cell_line}.tsv"
    params:
          runtime="6h"
    resources:
             mem_mb=1000
    # group: "bam_file_operation"
    threads: 4
    shell:
         "samtools sort {input} -o {output} --threads {threads} > {log} 2>&1 "

rule bam_index:
    input:
         "bams/with_filter/DNASE.{cell_line}.merge.filter.sorted.bam"
    output:
          "bams/with_filter/DNASE.{cell_line}.merge.filter.sorted.bam.bai"
    params:
          runtime="1h"
    resources:
             mem_mb=1000
    # group: "bam_file_operation"
    shell:
         "samtools index {input}"

# check SE or PE : { samtools view -H DNASE.HepG2.merge.filter.bam ;
# samtools view DNASE.HepG2.merge.filter.bam | head -n 1000; } | samtools view -c -f 1
rule generate_bw_5_end_1_bp:
    input:
         bam="bams/with_filter/DNASE.{cell_line}.merge.filter.sorted.bam",
         bai="bams/with_filter/DNASE.{cell_line}.merge.filter.sorted.bam.bai"
    output:
          forward_bw="bigwigs/DNASE.{cell_line}.merge.binSize.1.forward.bw",
          reverse_bw="bigwigs/DNASE.{cell_line}.merge.binSize.1.reverse.bw"
    threads: 8
    params:
          runtime="8h"
    resources:
             mem_mb=6000
    benchmark:
             "benchmarks/generate_bw_5_end_1_bp/{cell_line}.tsv"
    log:
       "logs/generate_bw_5_end_1_bp/{cell_line}.log"
    shell:
         """
         bamCoverage -b {input.bam} -o {output.forward_bw} --Offset 1 \
          --samFlagExclude 16  --binSize 1 --numberOfProcessors {threads} --ignoreDuplicates > {log} 2>&1 
         bamCoverage -b {input.bam} -o {output.reverse_bw} --Offset 1 \
          --samFlagInclude 16  --binSize 1 --numberOfProcessors {threads} --ignoreDuplicates > {log} 2>&1 
         """

# TODO remove unused h5 files after comparision
rule generate_DNase_features_from_bw_to_h5:
    input:
         forward_bw="bigwigs/DNASE.{cell_line}.merge.binSize.1.forward.bw",
         reverse_bw="bigwigs/DNASE.{cell_line}.merge.binSize.1.reverse.bw"
    output:
          # temp("hdf5s/DNase/DNASE_bam_5_mer_50bp_all_samples_lightGBM_{chrom}.{cell_line}.h5"),
          # temp("hdf5s/DNase/DNASE_bam_5_mer_100bp_all_samples_lightGBM_{chrom}.{cell_line}.h5"),
          "hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}.{cell_line}.h5"
    threads: 1
    resources:
             mem_mb=10000
    benchmark:
             "benchmarks/generate_DNase_features_from_bw_to_h5/{chrom}.{cell_line}.tsv"
    log:
       "logs/generate_DNase_features_from_bw_to_h5/{chrom}.{cell_line}.log"
    params:
          runtime="8h",
          config_file="./config.yaml"
    shell:
         "python scripts/DNase_features.py generate_dnase_features_from_bw_to_h5 --chrom={wildcards.chrom} "
         "--cell_line={wildcards.cell_line} --dnase_bw_file_path='./bigwigs' --outfile_path='./hdf5s/DNase/' "
         "--config_file='{params.config_file}' > {log} 2>&1"

def all_CA_sample_id_for_CA_clusters():
        # file_cluster_assigned = "peaks/TF_ChIP_sample_CA_cluster_matching/%s/" \
        #                         "%s_ChIP_sample_CA_cluster_assigned.xls" % (wildcards.tf_name, wildcards.tf_name)
        df_cluster_info_file = pd.read_csv(config['cluster_info_file'], sep=",", header=0)
        # df_cluster_info_file = df_cluster_info_file[df_cluster_info_file['included']=='Y']
        if config['selected_cluster']:
            int_clusters = list(map(int, config['selected_cluster']))
            return df_cluster_info_file[df_cluster_info_file['cluster'].isin(int_clusters)]['sample'].tolist()
        else:
            return df_cluster_info_file['sample'].tolist()

# checkpoint combine_dnase_features_h5_for_all_cell_types:
rule combine_dnase_features_h5_for_all_cell_types:
    input:
         # expand("hdf5s/DNase/DNASE_bam_5_mer_50bp_all_samples_lightGBM_{{chrom}}.{cell_line}.h5",
         #        cell_line=config['cell_types']),
         # expand("hdf5s/DNase/DNASE_bam_5_mer_100bp_all_samples_lightGBM_{{chrom}}.{cell_line}.h5",
         #        cell_line=config['cell_types']),
         expand("hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{{chrom}}.{cell_line}.h5",
                cell_line=config['cell_types']  if not config['ATAC_long_short']
                else ["%s_%s" % (cell_type,frag_len) for cell_type in config['cell_types'] for frag_len in ["short",'long']])
                if not config['CA_cluster_mode'] else
                expand("hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{{chrom}}.{cell_line}.h5",
                    cell_line=all_CA_sample_id_for_CA_clusters()
                   )
    output:
          # "hdf5s/DNase/DNASE_bam_5_mer_50bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
          # "hdf5s/DNase/DNASE_bam_5_mer_100bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
          "hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_all_cell_types.h5"
    threads: 1
    benchmark:
             "benchmarks/combine_dnase_features_h5_for_all_cell_types/{chrom}.tsv"
    log:
       "logs/combine_dnase_features_h5_for_all_cell_types/{chrom}.log"
    resources:
             mem_mb=10000
    params:
          runtime="4h",
          config_file="./config.yaml"
    shell:
         "python scripts/DNase_features.py combine_dnase_features_h5_for_all_cell_types "
         "--chrom={wildcards.chrom}  --outfile_path='./hdf5s/DNase' "
         "--config_file='{params.config_file}' > {log} 2>&1"

# checkpoint dnase_features_h5_for_all_cell_types_ready:
#     input:
#         # expand("hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_all_cell_types.h5", chrom=config['chrom_all']),
#         "hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_chr22_all_cell_types.h5"
#     output:
#         temp("train/{tf_name}/selected_motif_hdf5/dnase_features_h5_for_all_cell_types_ready.done")
#         # temp(touch("hdf5s/DNase/selected_motif_hdf5/dnase_features_h5_for_all_cell_types_ready.done"))
#     resources:
#         mem_mb=10
#     threads: 1
#     # shell:
#     #     "touch {output}"
#     run:
#         with open(output[0], 'w') as outfile:
#             with h5py.File(input[0], "r") as infile:
#                 feature_names_num = infile['feature_names'].shape[0]
#             outfile.write("%s" % feature_names_num)


rule generate_quantile_transformer:
    input:
         # expand("hdf5s/DNase/DNASE_bam_5_mer_50bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
         #        chrom=config['chrom_all']),
         # expand("hdf5s/DNase/DNASE_bam_5_mer_100bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
         #        chrom=config['chrom_all']),
         expand("hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
                chrom=config['chrom_all']),
    output:
          # "train/quantile_transformer/{cell_line}_50bp_quantile_map.pkl",
          # "train/quantile_transformer/{cell_line}_100bp_quantile_map.pkl",
          "train/data/quantile_transformer/{cell_line}_variable_bp_quantile_map.pkl"
    threads: 1
    benchmark:
             "benchmarks/generate_quantile_transformer/{cell_line}.tsv"
    log:
       "logs/generate_quantile_transformer/{cell_line}.log"
    resources:
             mem_mb=10000
    params:
          runtime="2h",
          config_file="./config.yaml"
    shell:
         "python scripts/DNase_features.py generate_quantile_transformer --cell_line={wildcards.cell_line} "
         "--dir_dnase_feature='./hdf5s/DNase' --dir_out='./train/data/quantile_transformer' "
         "--config_file='{params.config_file}' > {log} 2>&1"

rule generate_dnase_feature_median:
    input:
         # expand("hdf5s/DNase/DNASE_bam_5_mer_50bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
         #        chrom=config['chrom_all']),
         # expand("hdf5s/DNase/DNASE_bam_5_mer_100bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
         #        chrom=config['chrom_all']),
         "hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
         expand("train/data/quantile_transformer/{cell_line}_variable_bp_quantile_map.pkl",
                cell_line=config['cell_types'])
    output:
          "hdf5s/DNase/median/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_median.h5"
    threads: 1
    benchmark:
             "benchmarks/generate_dnase_feature_median/{chrom}.tsv"
    log:
       "logs/generate_dnase_feature_median/{chrom}.log"
    resources:
             mem_mb=80000
    params:
          runtime="1h",
          config_file="./config.yaml"
    shell:
         "python scripts/DNase_features.py generate_dnase_feature_median "
         "--chrom={wildcards.chrom}  --dir_dnase_feature='./hdf5s/DNase' --selected_bin_index_file=None "
         "--dir_quantile_transformer='train/data/quantile_transformer' --dir_out='hdf5s/DNase/median' "
         "--config_file='{params.config_file}' > {log} 2>&1"

rule prepare_motif_top4_feature:
    output:
          # expand("hdf5s/motif/{chrom}_motifs_top4_scores.h5", chrom=config['chrom_all'])
          "hdf5s/motif/{chrom}_motifs_top4_scores.h5"
    threads: 16
    benchmark:
             "benchmarks/prepare_motif_top4_feature/{chrom}.tsv"
    log:
       "logs/prepare_motif_top4_feature/{chrom}.log"
    resources:
             mem_mb=160000
    params:
          runtime="12h",
          chrom_size_file=config['chrom_size_file'],
          chrom_all=config['chrom_all'],
          genome_sequence_fa=config['genome_sequence_fa'],
          batch=config['batch'],
          motif_pwm_path=config['motif_pwm_path'],
    shell:
         "python scripts/motif_features.py MotifFeatures prepare_motif_top4_feature --chrom={wildcards.chrom} --num_threads={threads} "
         "--motif_pwm_path='{params.motif_pwm_path}' --motif_feature_path='./hdf5s/motif' "
         "--chrom_size_file='{params.chrom_size_file}' --chrom_all='{params.chrom_all}' --batch={params.batch} "
         "--genome_sequence_fa='{params.genome_sequence_fa}' > {log} 2>&1"

# def all_TF_ChIP_sample_id(wildcards):
#     df_TF_ChIP_sample_meta = pd.read_csv(config['TF_sample_info_file'], sep="\t", header=0)
#     return df_TF_ChIP_sample_meta[df_TF_ChIP_sample_meta['factor']==wildcards.tf_name]['DCid'].tolist()

def all_TF_ChIP_sample_id(tf_name):
    df_TF_ChIP_sample_meta = pd.read_csv(config['TF_sample_info_file'], sep="\t", header=0)
    return df_TF_ChIP_sample_meta[df_TF_ChIP_sample_meta['factor']==tf_name]['DCid'].tolist()

# sample_id=all_TF_ChIP_sample_id(wildcards)
if config['CA_cluster_mode']:
    rule TF_ChIP_sample_CA_cluster_matching:
        input:
         # lambda wildcards: "{peak_path}/ChIPseq.{cell_line}.{tf_name}.conservative.train.narrowPeak.gz".format(
             lambda wildcards: "{peak_path}/{sample_id}_sort_peaks.narrowPeak.bed".format(
                 peak_path=config['cistromeDB_narrow_peak_path'],
                 sample_id=wildcards.sample_id)
        output:
              "peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/{sample_id}.giggle.txt"
        # group: "peak_operation"
        resources:
                 mem_mb=1000
        params:
              runtime="1h",
              giggle_path=config['giggle_path'],
              giggle_index_path=config['giggle_index_path']
        shell:
            '''
            cut -f 1,2,3 {input} | head -n 10000 > peaks/TF_ChIP_sample_CA_cluster_matching/{wildcards.tf_name}/{wildcards.sample_id}.bed || true
            bedtools sort -i peaks/TF_ChIP_sample_CA_cluster_matching/{wildcards.tf_name}/{wildcards.sample_id}.bed > \
                peaks/TF_ChIP_sample_CA_cluster_matching/{wildcards.tf_name}/{wildcards.sample_id}.sorted.bed
            bgzip peaks/TF_ChIP_sample_CA_cluster_matching/{wildcards.tf_name}/{wildcards.sample_id}.sorted.bed -c > \
                peaks/TF_ChIP_sample_CA_cluster_matching/{wildcards.tf_name}/{wildcards.sample_id}.sorted.bed.gz
            {params.giggle_path} search -i {params.giggle_index_path} \
                -q peaks/TF_ChIP_sample_CA_cluster_matching/{wildcards.tf_name}/{wildcards.sample_id}.sorted.bed.gz \
                -s > {output}
            '''

    checkpoint combine_TF_ChIP_sample_CA_cluster_results:
        input:
             lambda wildcards: expand("peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/{sample_id}.giggle.txt",
                                   tf_name=wildcards.tf_name, sample_id=all_TF_ChIP_sample_id(wildcards.tf_name))

        output:
              "peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/{tf_name}_ChIP_sample_CA_cluster_assigned.xls"
        # group: "peak_operation"
        resources:
                 mem_mb=1000
        params:
              runtime="1h",
              TF_sample_info_file=config['TF_sample_info_file'],
              cluster_info_file=config['cluster_info_file']
        run:
            df_TF_ChIP_sample_meta = pd.read_csv(params.TF_sample_info_file, sep="\t", header=0)
            df_temp = df_TF_ChIP_sample_meta[df_TF_ChIP_sample_meta['factor']==wildcards.tf_name].copy()
            dic_sampel_id_cluster_id = {}
            for filename in input:
                df_giggle = pd.read_csv(filename,sep='\t',header=0)
                column_names=df_giggle.columns.tolist()
                column_names[0]='filename'
                df_giggle=df_giggle.reset_index()
                del df_giggle['combo_score']
                df_giggle.columns=column_names
                df_giggle.sort_values('combo_score',inplace=True,ascending=False)
                dic_sampel_id_cluster_id[int(filename.split('/')[-1].split('.')[0])]=int(df_giggle.iloc[0,:]['filename'].split("_")[2])
            df_temp.loc[:,'cluster'] = df_temp['DCid'].apply(lambda x:dic_sampel_id_cluster_id[x])
            # print(df_temp.shape)
            df_cluster_info_file = pd.read_csv(params.cluster_info_file, sep=",", header=0)
            df_cluster_info_file = df_cluster_info_file[df_cluster_info_file['included']=='Y'][['cluster','annotation']]
            df_cluster_info_file.drop_duplicates(inplace=True)
            # print(output[0])
            # print(df_temp.shape)
            df_temp.merge(df_cluster_info_file[['cluster','annotation']], on='cluster', how='left').to_csv(output[0],
                                                                                      sep='\t', header=True, index=False)
    def all_TF_ChIP_sample_id_for_same_CA_cluster(wildcards):
        # file_cluster_assigned = "peaks/TF_ChIP_sample_CA_cluster_matching/%s/" \
        #                         "%s_ChIP_sample_CA_cluster_assigned.xls" % (wildcards.tf_name, wildcards.tf_name)
        file_cluster_assigned = checkpoints.combine_TF_ChIP_sample_CA_cluster_results.get(tf_name=wildcards.tf_name).output[0]
        df_cluster_assigned = pd.read_csv(file_cluster_assigned, sep="\t", header=0)
        return df_cluster_assigned[df_cluster_assigned['cluster']==int(wildcards.cluster_id)]['DCid'].tolist()

    rule combine_TF_ChIP_sample_peaks_for_same_CA_cluster:
        input:
             # "peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/{tf_name}_ChIP_sample_CA_cluster_assigned.xls",
             lambda wildcards: expand("peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/{sample_id}.bed",
                 tf_name=wildcards.tf_name,
                 sample_id=all_TF_ChIP_sample_id_for_same_CA_cluster(wildcards))
        output:
              "peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/tf_name.{tf_name}.cluster_id.{cluster_id}.peak.bed"
        benchmark:
             "benchmarks/combine_TF_ChIP_sample_peaks_for_same_CA_cluster/{tf_name}.{cluster_id}.tsv"
        log:
           "logs/combine_TF_ChIP_sample_peaks_for_same_CA_cluster/{tf_name}.{cluster_id}.log"
        resources:
                 mem_mb=10000
        params:
              runtime="1h",
        run:
            # df_TF_ChIP_sample_CA_cluster_assigned = pd.read_csv(input[0], sep="\t", header=0)
            # df_temp = df_TF_ChIP_sample_CA_cluster_assigned[df_TF_ChIP_sample_CA_cluster_assigned['cluster']==int(wildcards.cluster_id)]
            path = 'peaks/TF_ChIP_sample_CA_cluster_matching/%s' % wildcards.tf_name
            # str_bedfiles = " ".join(list(map(lambda x: "%s/%d.bed" % (path, x), df_temp['DCid'].unique().tolist())))
            str_bedfiles = " ".join(list(map(lambda x: "%s/%d.bed" % (path, x), all_TF_ChIP_sample_id_for_same_CA_cluster(wildcards))))
            # print(str_bedfiles)
            _ = os.system("cat %s > %s/%d.cat.bed" % (str_bedfiles, path, int(wildcards.cluster_id)))
            _ = os.system("bedtools sort -i %s/%d.cat.bed > %s/%d.sorted.bed" % (path, int(wildcards.cluster_id), path, int(wildcards.cluster_id)))
            _ = os.system("bedtools merge -i %s/%d.sorted.bed > %s" % (path, int(wildcards.cluster_id), output[0]))


    def all_CA_cluster(wildcards):
        # file_cluster_assigned = "peaks/TF_ChIP_sample_CA_cluster_matching/%s/" \
                                # "%s_ChIP_sample_CA_cluster_assigned.xls" % (wildcards.tf_name, wildcards.tf_name)
        file_cluster_assigned = checkpoints.combine_TF_ChIP_sample_CA_cluster_results.get(tf_name=wildcards.tf_name).output[0]
        df_cluster_assigned = pd.read_csv(file_cluster_assigned, sep="\t", header=0)
        return df_cluster_assigned['cluster'].unique().tolist()

    checkpoint prepare_label_file:
        input:
             # "peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/{tf_name}_ChIP_sample_CA_cluster_assigned.xls",
             lambda wildcards: expand("peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/tf_name.{tf_name}.cluster_id.{cluster_id}.peak.bed",
                 tf_name=wildcards.tf_name,
                 cluster_id=all_CA_cluster(wildcards))
        output:
              "{training_cell_types_regions_label_path}/{{tf_name}}.{training_cell_types_regions_label_name}".format(
                  training_cell_types_regions_label_path = config['training_cell_types_regions_label_path'],
                  training_cell_types_regions_label_name = config['training_cell_types_regions_label_name'])
        benchmark:
             "benchmarks/prepare_label_file/{tf_name}.tsv"
        log:
           "logs/prepare_label_file/{tf_name}.log"
        resources:
                 mem_mb=10000
        params:
              runtime="1h",
              regions_all_file=config['regions_all_file'],
              training_cell_types_regions_label_path=config['training_cell_types_regions_label_path']
        run:
            if not os.path.exists('{training_cell_types_regions_label_path}/{tf_name}'.format(
                                     tf_name=wildcards.tf_name,
                                     training_cell_types_regions_label_path=params.training_cell_types_regions_label_path
                    )):
                os.mkdir('{training_cell_types_regions_label_path}/{tf_name}'.format(
                                         tf_name=wildcards.tf_name,
                                         training_cell_types_regions_label_path=params.training_cell_types_regions_label_path
                        ))
            list_DNase_sample_id=[]
            for ind,cluster_id in enumerate(all_CA_cluster(wildcards)):
                list_DNase_sample_id.append(cluster_id)
                if ind==0:
                    _ = os.system('bedtools intersect -c -a {regions_all_file} \
                        -b  peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/tf_name.{tf_name}.cluster_id.{cluster_id}.peak.bed\
                         > {training_cell_types_regions_label_path}/{tf_name}/{tf_name}_train_temp_{ind}.bed'.format(
                                     regions_all_file=params.regions_all_file,
                                     tf_name=wildcards.tf_name,
                                     cluster_id=cluster_id,
                                     training_cell_types_regions_label_path=params.training_cell_types_regions_label_path,
                                     ind=ind
                    ))
                else:
                    _ = os.system('bedtools intersect -c -a {training_cell_types_regions_label_path}/{tf_name}/{tf_name}_train_temp_{prev_ind}.bed \
                        -b  peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/tf_name.{tf_name}.cluster_id.{cluster_id}.peak.bed\
                         > {training_cell_types_regions_label_path}/{tf_name}/{tf_name}_train_temp_{ind}.bed'.format(
                                     regions_all_file=params.regions_all_file,
                                     tf_name=wildcards.tf_name,
                                     cluster_id=cluster_id,
                                     training_cell_types_regions_label_path=params.training_cell_types_regions_label_path,
                                     ind=ind,
                                     prev_ind=ind-1
                    ))
            df_temp=pd.read_csv('{training_cell_types_regions_label_path}/{tf_name}/{tf_name}_train_temp_{ind}.bed'.format(
                                     tf_name=wildcards.tf_name,
                                     training_cell_types_regions_label_path=params.training_cell_types_regions_label_path,
                                     ind=len(list_DNase_sample_id)-1
                                ),sep="\t",header=None)
            df_temp.columns=['chr', 'start', 'stop']+list(map(lambda x: str(x), list_DNase_sample_id))
            for sample_id in list_DNase_sample_id:
                df_temp[str(sample_id)]=df_temp[str(sample_id)].apply(lambda x: 'U' if x==0 else "B")
            df_temp.to_csv(output[0], sep="\t", header=True,index=False)
            for ind in range(len(list_DNase_sample_id)):
                os.remove('{training_cell_types_regions_label_path}/{tf_name}/{tf_name}_train_temp_{ind}.bed'.format(
                                     tf_name=wildcards.tf_name,
                                     training_cell_types_regions_label_path=params.training_cell_types_regions_label_path,
                                     ind=ind
                    ))
            os.rmdir('{training_cell_types_regions_label_path}/{tf_name}'.format(
                                     tf_name=wildcards.tf_name,
                                     training_cell_types_regions_label_path=params.training_cell_types_regions_label_path
                    ))






rule top_5k_peaks:
    input:
         # lambda wildcards: "{peak_path}/ChIPseq.{cell_line}.{tf_name}.conservative.train.narrowPeak.gz".format(
         lambda wildcards: "peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/tf_name.{tf_name}.cluster_id.{cell_line}.peak.bed".format(
             cell_line=wildcards.cell_line, tf_name=wildcards.tf_name)
            if config['CA_cluster_mode']
            else "{peak_path}/{cell_line}.{tf_name}.sort_peaks.narrowPeak.bed".format(
             peak_path=config['peak_file_path'],
             cell_line=wildcards.cell_line, tf_name=wildcards.tf_name)
    output:
          "peaks/top5k/ChIPseq.{cell_line}.{tf_name}.conservative.train.top5k.narrowPeak"
    # group: "peak_operation"
    resources:
             mem_mb=1000
    params:
          runtime="1h",
    shell:
         # "zcat {input} | head -n 5000 > {output} || true "
         # "cat {input} | head -n 5000 > {output} || true "
         "cat {input} > {output} || true "

rule seqpos_scan_motifs_on_peaks:
    input:
         "peaks/top5k/ChIPseq.{cell_line}.{tf_name}.conservative.train.top5k.narrowPeak"
    output:
          "train/{tf_name}/motif_scan/cell_line.{cell_line}.tf.{tf_name}/mdseqpos_index.html"
    benchmark:
             "benchmarks/seqpos_scan_motifs_on_peaks/{tf_name}.{cell_line}.tsv"
    log:
       "logs/seqpos_scan_motifs_on_peaks/{tf_name}.{cell_line}.log"
    # group: "peak_operation"
    threads: 1
    resources:
             mem_mb=10000
    params:
          runtime="10h",
          motif_database_xml=config['motif_database_xml'],
          genome=config['genome'],
          seqpos_env_path=config['seqpos_env_path']
    shell:
         "{params.seqpos_env_path}/bin/python2.7 {params.seqpos_env_path}/bin/MDSeqPos.py {input} {params.genome} "
         "-m {params.motif_database_xml} -d "
         "-O ./train/{wildcards.tf_name}/motif_scan/cell_line.{wildcards.cell_line}.tf.{wildcards.tf_name} > {log} 2>&1"


def training_cell_line(wildcards):
    if config['CA_cluster_mode']:
        return list(glob_wildcards("peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/"
                                   "tf_name.{tf_name}.cluster_id.{{cell_line}}.peak.bed".format(
                    tf_name=wildcards.tf_name
                ))[0])
        "peaks/TF_ChIP_sample_CA_cluster_matching/{tf_name}/tf_name.{tf_name}.cluster_id.{cluster_id}.peak.bed"
    else:
        # return list(glob_wildcards("{peak_path}/ChIPseq.{{cell_line}}.{tf_name}.conservative.train.narrowPeak.gz".format(
        return list(glob_wildcards("{peak_path}/{{cell_line}}.{tf_name}.sort_peaks.narrowPeak.bed".format(
            peak_path=config['peak_file_path'],
            tf_name=wildcards.tf_name
        ))[0])


rule decide_cofactor_motif_set:
    input:
         lambda wildcards: expand("train/{tf_name}/motif_scan/cell_line.{cell_line}.tf.{tf_name}/mdseqpos_index.html",
                                  cell_line=cell_line_for_train(wildcards.tf_name),
                                  # cell_line=training_cell_line(wildcards),
                                  tf_name=wildcards.tf_name),
    output:
          "train/{tf_name}/{tf_name}_cofactor_motif_list.json"
    threads: 1
    # group: "cofactor_motif_set_operation"
    resources:
             mem_mb=1000
    params:
          runtime="1h",
          topn=10
    shell:
         "python scripts/motif_features.py decide_cofactor_motif_set --outfilename={output} --n={params.topn} {input}"

rule prepare_motif_h5_data:
    input:
         cofactor_motif_set_file="train/{tf_name}/{tf_name}_cofactor_motif_list.json",
         motif_top4_feature="hdf5s/motif/{chrom}_motifs_top4_scores.h5"
    output:
          "train/{tf_name}/selected_motif_hdf5/{chrom}_motif_features_lightGBM.h5"
    threads: 1
    # group: "cofactor_motif_set_operation"
    resources:
             # mem_mb=60000
             mem_mb=120000
    benchmark:
             "benchmarks/prepare_motif_h5_data/{tf_name}.{chrom}.tsv"
    params:
          runtime="6h",
          config_file="./config.yaml"
    log:
       "logs/prepare_motif_h5_data/{tf_name}.{chrom}.log"
    shell:
         "python scripts/train_lightgbm_model.py prepare_motif_h5_data --chrom={wildcards.chrom} "
         "--config_file='{params.config_file}' --cofactor_motif_set_file='{input.cofactor_motif_set_file}' "
         "--training_tf_name='{wildcards.tf_name}' "
         "--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' "
         "--motif_feature_path='./hdf5s/motif' "
         "--selected_motif_feature_path='train/{wildcards.tf_name}/selected_motif_hdf5/' "
         "--step=120 > {log} 2>&1 "

rule prepare_dnase_autoencoder_h5_data:
    input:
         "/n/scratchlfs/xiaoleliu_lab/cchen/Cistrome_imputation/encode/data/DNase_scanning/scan_result/"
         "DNASE.{cell_line}.merge.binSize.1.corrected_sorted_hg19_25bpbin_bwaverage_transformed_{chrom}_scanned_with_autoencoder.hdf5"
    output:
          "hdf5s/DNase/DNASE_autoencoder_lightGBM.{chrom}.{cell_line}.h5"
    threads: 1
    resources:
             mem_mb=20000
    benchmark:
             "benchmarks/prepare_dnase_autoencoder_h5_data/{chrom}.{cell_line}.tsv"
    params:
          runtime="1h",
          config_file="./config.yaml"
    log:
       "logs/prepare_dnase_autoencoder_h5_data/{chrom}.{cell_line}.log"
    shell:
         "python scripts/train_lightgbm_model.py prepare_dnase_autoencoder_h5_data --chrom={wildcards.chrom} "
         "--cell_line={wildcards.cell_line} --outfile_path='./hdf5s/DNase' "
         "--config_file='{params.config_file}' "
         "--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' "
         "--motif_feature_path='./hdf5s/motif' --training_tf_name='HNF4A' "
         "--step=120 > {log} 2>&1 "

checkpoint motif_h5_data_ready:
    input:
         # expand("train/{{tf_name}}/selected_motif_hdf5/{chrom}_motif_features_lightGBM.h5", chrom=config['chrom_all']),
         # "train/{tf_name}/selected_motif_hdf5/chr19_motif_features_lightGBM.h5",
         expand("train/{{tf_name}}/selected_motif_hdf5/{chrom}_motif_features_lightGBM.h5",
                chrom=config['chrom_all'])
    output:
          "train/{tf_name}/selected_motif_hdf5/motif_h5_data_ready.done"
          # temp(touch("train/{tf_name}/selected_motif_hdf5/motif_h5_data_ready.done"))
    # shell:
    #     "touch {output}"
    run:
        with open(output[0], 'w') as outfile:
            with h5py.File(input[0], "r") as infile:
                feature_names_num = infile['feature_names'].shape[0]
            outfile.write("%s" % feature_names_num)

#
# def selected_cell_line_for_reference(wildcards, n=1):
#     df_all_regions_label = pd.read_csv("%s/%s.%s" % (
#         config['training_cell_types_regions_label_path'], wildcards.tf_name,
#         config['training_cell_types_regions_label_name']), sep="\t", header=0, nrows=10)
#     if n == 1:
#         return [df_all_regions_label.columns[3]]
#     else:
#         return list(df_all_regions_label.columns[3:])


rule prepare_lightgbm_binary_data_motif_feature_reference:
    input:
         expand("train/{{tf_name}}/selected_motif_hdf5/{chrom}_motif_features_lightGBM.h5",
                chrom=config['chrom_sets']['chrom_setA'])
    output:
          "train/{tf_name}/binary_files/reference/lightGBM.motif.chrom_setA.{subset_index}.bin"
    threads: 1
    benchmark:
             "benchmarks/prepare_lightgbm_binary_data_motif_feature_reference/{tf_name}.chrom_setA.{subset_index}.tsv"
    log:
       "logs/prepare_lightgbm_binary_data_motif_feature_reference/{tf_name}.chrom_setA.{subset_index}.log"
    resources:
             mem_mb=40000
    params:
          runtime="4h",
          config_file="./config.yaml",
    shell:
         "python scripts/train_lightgbm_model.py prepare_lightgbm_binary_data_motif_feature_subset "
         "--chrom_set_name='chrom_setA' --subset_index={wildcards.subset_index} "
         "--dir_out='./train/{wildcards.tf_name}/binary_files/reference' "
         "--config_file='{params.config_file}' "
         "--training_tf_name='{wildcards.tf_name}' "
         "--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' "
         "--motif_feature_path='./hdf5s/motif' "
         "--selected_motif_feature_path='train/{wildcards.tf_name}/selected_motif_hdf5/' "
         "--step=120 > {log} 2>&1"

rule copy_binary_files_motif_reference:
    input:
         "train/{tf_name}/binary_files/reference/lightGBM.motif.chrom_setA.{subset_index}.bin"
    output:
          "train/{tf_name}/binary_files/lightGBM.motif.chrom_setA.{subset_index}.bin",
          # touch(
          #     "train/data/dnase_feature_binary_files/reference/"
          #     "copy_binary_files_reference.{cell_line}.chrom_set_test.done")
    threads: 1
    resources:
             mem_mb=1000
    params:
          runtime="1h",
    shell:
         "cp {input} {output[0]}"

rule prepare_lightgbm_binary_data_motif_feature_subset:
    input:
         lambda wildcards: expand("train/{{tf_name}}/selected_motif_hdf5/{chrom}_motif_features_lightGBM.h5",
                                  chrom=config['chrom_sets'][wildcards.chrom_set_name]),
         lambda wildcards: expand("train/{{tf_name}}/binary_files/lightGBM.motif.chrom_setA.{all_subset_index}.bin",
                                  all_subset_index=subset_index_list(wildcards))
    output:
          "train/{tf_name}/binary_files/lightGBM.motif.{chrom_set_name}.{subset_index}.bin"
    threads: 1
    benchmark:
             "benchmarks/prepare_lightgbm_binary_data_motif_feature_subset/{tf_name}.{chrom_set_name}.{subset_index}.tsv"
    log:
       "logs/prepare_lightgbm_binary_data_motif_feature_subset/{tf_name}.{chrom_set_name}.{subset_index}.log"
    resources:
             mem_mb=40000
    params:
          runtime="4h",
          config_file="./config.yaml",
    shell:
         "python scripts/train_lightgbm_model.py prepare_lightgbm_binary_data_motif_feature_subset "
         "--chrom_set_name={wildcards.chrom_set_name} --subset_index={wildcards.subset_index} "
         "--dir_out='./train/{wildcards.tf_name}/binary_files' "
         "--config_file='{params.config_file}' --reference='./train/{wildcards.tf_name}/binary_files/reference' "
         "--training_tf_name='{wildcards.tf_name}' "
         "--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' "
         "--motif_feature_path='./hdf5s/motif' "
         "--selected_motif_feature_path='train/{wildcards.tf_name}/selected_motif_hdf5/' "
         "--step=120 > {log} 2>&1"

if config['ATAC_long_short']:
    rule prepare_lightgbm_binary_data_dnase_feature_reference:
        input:
             # expand("hdf5s/DNase/DNASE_bam_5_mer_50bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
             #        chrom=config['chrom_all']),
             # expand("hdf5s/DNase/DNASE_bam_5_mer_100bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
             #        chrom=config['chrom_all']),
             expand("hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
                    chrom=config['chrom_sets']['chrom_set_test']),
             # expand("hdf5s/DNase/median/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_median.h5",
             #        chrom=config['chrom_sets']['chrom_set_test']),
             "train/data/quantile_transformer/{cell_line}_short_variable_bp_quantile_map.pkl",
             "train/data/quantile_transformer/{cell_line}_long_variable_bp_quantile_map.pkl"
        output:
              "train/data/dnase_feature_binary_files/reference/lightGBM.dnase.{cell_line}.chrom_set_test.bin"
        threads: 1
        benchmark:
                 "benchmarks/prepare_lightgbm_binary_data_dnase_feature_reference/{cell_line}.chrom_set_test.tsv"
        log:
           "logs/prepare_lightgbm_binary_data_dnase_feature_reference/{cell_line}.chrom_set_test.log"
        resources:
                 mem_mb=20000
        params:
              runtime="4h",
              config_file="./config.yaml",
              ATAC_long_short=config['ATAC_long_short']
        shell:
             "python scripts/train_lightgbm_model.py prepare_lightgbm_binary_dnase_feature --cell_line={wildcards.cell_line} "
             "--chrom_set_name='chrom_set_test' --dir_dnase_feature_median='hdf5s/DNase/median' "
             "--dir_out='./train/data/dnase_feature_binary_files/reference' --ATAC_long_short='{params.ATAC_long_short}' "
             "--config_file='{params.config_file}' "
             "--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' "
             "--motif_feature_path='./hdf5s/motif' "
             "--step=120 > {log} 2>&1"


    rule prepare_lightgbm_binary_data_dnase_feature_all:
        input:
             # expand("hdf5s/DNase/DNASE_bam_5_mer_50bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
             #        chrom=config['chrom_all']),
             # expand("hdf5s/DNase/DNASE_bam_5_mer_100bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
             #        chrom=config['chrom_all']),
             lambda wildcards: expand(
                 "hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
                 chrom=config['chrom_sets'][wildcards.chrom_set_name]),
             # lambda wildcards: expand(
                 # "hdf5s/DNase/median/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_median.h5",
                 # chrom=config['chrom_sets'][wildcards.chrom_set_name]),
             # "train/data/dnase_feature_binary_files/reference/copy_binary_files_reference.{}.chrom_set_test.done".format(
             #     config['reference_cell_type']),
             "train/data/dnase_feature_binary_files/lightGBM.dnase.{}.chrom_set_test.bin".format(
                 config['reference_cell_type']),
             "train/data/quantile_transformer/{cell_line}_short_variable_bp_quantile_map.pkl",
             "train/data/quantile_transformer/{cell_line}_long_variable_bp_quantile_map.pkl"
        output:
              "train/data/dnase_feature_binary_files/lightGBM.dnase.{cell_line}.{chrom_set_name}.bin"
        threads: 1
        benchmark:
                 "benchmarks/prepare_lightgbm_binary_data_dnase_feature_all/{cell_line}.{chrom_set_name}.tsv"
        log:
           "logs/prepare_lightgbm_binary_data_dnase_feature_all/{cell_line}.{chrom_set_name}.log"
        resources:
                 mem_mb=40000
        params:
              runtime="4h",
              config_file="./config.yaml",
              reference_cell_type=config['reference_cell_type'],
              ATAC_long_short=config['ATAC_long_short']
        shell:
             "python scripts/train_lightgbm_model.py prepare_lightgbm_binary_dnase_feature --cell_line={wildcards.cell_line} "
             "--chrom_set_name={wildcards.chrom_set_name} --dir_dnase_feature_median='hdf5s/DNase/median' "
             "--dir_out='./train/data/dnase_feature_binary_files' --ATAC_long_short='{params.ATAC_long_short}' "
             "--config_file='{params.config_file}' --reference='./train/data/dnase_feature_binary_files/reference/"
             "lightGBM.dnase.{params.reference_cell_type}.chrom_set_test.bin' "
             "--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' "
             "--motif_feature_path='./hdf5s/motif' "
             "--step=120 > {log} 2>&1"
else:
    rule prepare_lightgbm_binary_data_dnase_feature_reference:
        input:
             # expand("hdf5s/DNase/DNASE_bam_5_mer_50bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
             #        chrom=config['chrom_all']),
             # expand("hdf5s/DNase/DNASE_bam_5_mer_100bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
             #        chrom=config['chrom_all']),
             expand("hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
                    chrom=config['chrom_sets']['chrom_set_test']),
             # expand("hdf5s/DNase/median/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_median.h5",
             #        chrom=config['chrom_sets']['chrom_set_test']),
             "train/data/quantile_transformer/{cell_line}_variable_bp_quantile_map.pkl"
        output:
              "train/data/dnase_feature_binary_files/reference/lightGBM.dnase.{cell_line}.chrom_set_test.bin"
        threads: 1
        benchmark:
                 "benchmarks/prepare_lightgbm_binary_data_dnase_feature_reference/{cell_line}.chrom_set_test.tsv"
        log:
           "logs/prepare_lightgbm_binary_data_dnase_feature_reference/{cell_line}.chrom_set_test.log"
        resources:
                 mem_mb=20000
        params:
              runtime="4h",
              config_file="./config.yaml",
        shell:
             "python scripts/train_lightgbm_model.py prepare_lightgbm_binary_dnase_feature --cell_line={wildcards.cell_line} "
             "--chrom_set_name='chrom_set_test' --dir_dnase_feature_median='hdf5s/DNase/median' "
             "--dir_out='./train/data/dnase_feature_binary_files/reference' "
             "--config_file='{params.config_file}' "
             "--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' "
             "--motif_feature_path='./hdf5s/motif' "
             "--step=120 > {log} 2>&1"


    rule prepare_lightgbm_binary_data_dnase_feature_all:
        input:
             # expand("hdf5s/DNase/DNASE_bam_5_mer_50bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
             #        chrom=config['chrom_all']),
             # expand("hdf5s/DNase/DNASE_bam_5_mer_100bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
             #        chrom=config['chrom_all']),
             lambda wildcards: expand(
                 "hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
                 chrom=config['chrom_sets'][wildcards.chrom_set_name]),
             # lambda wildcards: expand(
                 # "hdf5s/DNase/median/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_median.h5",
                 # chrom=config['chrom_sets'][wildcards.chrom_set_name]),
             # "train/data/dnase_feature_binary_files/reference/copy_binary_files_reference.{}.chrom_set_test.done".format(
             #     config['reference_cell_type']),
             "train/data/dnase_feature_binary_files/lightGBM.dnase.{}.chrom_set_test.bin".format(
                 config['reference_cell_type']),
             "train/data/quantile_transformer/{cell_line}_variable_bp_quantile_map.pkl"
        output:
              "train/data/dnase_feature_binary_files/lightGBM.dnase.{cell_line}.{chrom_set_name}.bin"
        threads: 1
        benchmark:
                 "benchmarks/prepare_lightgbm_binary_data_dnase_feature_all/{cell_line}.{chrom_set_name}.tsv"
        log:
           "logs/prepare_lightgbm_binary_data_dnase_feature_all/{cell_line}.{chrom_set_name}.log"
        resources:
                 mem_mb=40000
        params:
              runtime="4h",
              config_file="./config.yaml",
              reference_cell_type=config['reference_cell_type'],
        shell:
             "python scripts/train_lightgbm_model.py prepare_lightgbm_binary_dnase_feature --cell_line={wildcards.cell_line} "
             "--chrom_set_name={wildcards.chrom_set_name} --dir_dnase_feature_median='hdf5s/DNase/median' "
             "--dir_out='./train/data/dnase_feature_binary_files' "
             "--config_file='{params.config_file}' --reference='./train/data/dnase_feature_binary_files/reference/"
             "lightGBM.dnase.{params.reference_cell_type}.chrom_set_test.bin' "
             "--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' "
             "--motif_feature_path='./hdf5s/motif' "
             "--step=120 > {log} 2>&1"


def subset_index_list(wildcards):
    filename = checkpoints.motif_h5_data_ready.get(tf_name=wildcards.tf_name).output[0]
    with open(filename, "r") as infile:
        all_feature_names_len = int(infile.readline().strip())
    return list(range(1, int(np.ceil(all_feature_names_len / 120) + 1)))


rule copy_binary_files_reference:
    input:
         "train/data/dnase_feature_binary_files/reference/lightGBM.dnase.{cell_line}.chrom_set_test.bin"
    output:
          "train/data/dnase_feature_binary_files/lightGBM.dnase.{cell_line}.chrom_set_test.bin",
          # touch(
          #     "train/data/dnase_feature_binary_files/reference/"
          #     "copy_binary_files_reference.{cell_line}.chrom_set_test.done")
    threads: 1
    resources:
             mem_mb=1000
    params:
          runtime="1h",
    shell:
         "cp {input} {output[0]}"

rule merge_lightgbm_binary_data:
    input:
         "train/{tf_name}/selected_motif_hdf5/chr19_motif_features_lightGBM.h5",
         "train/data/dnase_feature_binary_files/lightGBM.dnase.{cell_line}.{chrom_set_name}.bin",
         lambda wildcards: expand(
             "train/{{tf_name}}/binary_files/lightGBM.motif.{{chrom_set_name}}.{subset_index}.bin",
             subset_index=subset_index_list(wildcards)
         ),
         "{training_cell_types_regions_label_path}/{{tf_name}}.{training_cell_types_regions_label_name}".format(
                  training_cell_types_regions_label_path = config['training_cell_types_regions_label_path'],
                  training_cell_types_regions_label_name = config['training_cell_types_regions_label_name'])
    output:
          "train/{tf_name}/binary_files/lightGBM.all.{cell_line}.{chrom_set_name}.bin"
    threads: 1
    wildcard_constraints:
                        subset_index="\d+"
    benchmark:
             "benchmarks/merge_lightgbm_binary_data/{tf_name}.{cell_line}.{chrom_set_name}.tsv"
    log:
       "logs/merge_lightgbm_binary_data/{tf_name}.{cell_line}.{chrom_set_name}.log"
    resources:
             mem_mb=140000
    params:
          runtime="4h",
          config_file="./config.yaml"
    shell:
         "python scripts/train_lightgbm_model.py merge_lightgbm_binary_data --cell_line={wildcards.cell_line} "
         "--dir_out='./train/{wildcards.tf_name}/binary_files' --chrom_set_name='{wildcards.chrom_set_name}' "
         "--lightgbm_dnase_binary_files_path='./train/data/dnase_feature_binary_files' "
         "--lightgbm_motif_binary_files_path='./train/{wildcards.tf_name}/binary_files' "
         "--config_file='{params.config_file}' "
         "--training_tf_name='{wildcards.tf_name}' "
         "--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' "
         "--motif_feature_path='./hdf5s/motif' "
         "--selected_motif_feature_path='train/{wildcards.tf_name}/selected_motif_hdf5/' "
         "--step=120 > {log} 2>&1"


rule train_models:
    input:
         # lambda wildcards: expand(
         #     "train/data/dnase_feature_binary_files/lightGBM.dnase.{cell_line}.{chrom_set_name}.bin",
         #     cell_line=selected_cell_line_for_reference(wildcards, 2),
         #     chrom_set_name=list(set(config['chrom_sets'].keys()))),
         # lambda wildcards: expand("train/{{tf_name}}/binary_files/lightGBM.motif.{chrom_set_name}.{subset_index}.bin",
         #                          subset_index=subset_index_list(wildcards),
         #                          chrom_set_name=list(set(config['chrom_sets'].keys()))
         #                          )
         lambda wildcards: expand(
             "train/{{tf_name}}/binary_files/lightGBM.all.{cell_line}.{chrom_set_name}.bin",
             # cell_line=selected_cell_line_for_reference(wildcards, 2),
             cell_line=cell_line_for_train(wildcards.tf_name),
             chrom_set_name=list(set(config['chrom_sets'].keys()))),
    output:
          "train/{tf_name}/models/{tf_name}.{training_cell_line}.{training_chrom_set_name}_model.pkl",
          "train/{tf_name}/models/{tf_name}.{training_cell_line}.{training_chrom_set_name}_evals_result.pkl",
    threads: 16
    benchmark:
             "benchmarks/train_models/{tf_name}.{training_cell_line}.{training_chrom_set_name}.tsv"
    log:
       "logs/train_models/{tf_name}.{training_cell_line}.{training_chrom_set_name}.log"
    resources:
             mem_mb=120000
    params:
          runtime="4h",
          config_file="./config.yaml"
    shell:
         "python scripts/train_lightgbm_model.py train_models --cell_line={wildcards.training_cell_line} "
         "--dir_out='./train/{wildcards.tf_name}/models' --chrom_set_name={wildcards.training_chrom_set_name} "
         "--lightgbm_dnase_binary_files_path='./train/data/dnase_feature_binary_files' "
         "--lightgbm_motif_binary_files_path='./train/{wildcards.tf_name}/binary_files' --num_threads={threads} "
         "--config_file='{params.config_file}' "
         "--training_tf_name='{wildcards.tf_name}' "
         "--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' "
         "--motif_feature_path='./hdf5s/motif' "
         "--selected_motif_feature_path='train/{wildcards.tf_name}/selected_motif_hdf5/' "
         "--step=120  > {log} 2>&1"



rule make_prediction_chrom:
    input:
         # expand("hdf5s/DNase/DNASE_bam_5_mer_50bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
         #         chrom=config['chrom_all']),
         # expand("hdf5s/DNase/DNASE_bam_5_mer_100bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
         #     chrom=config['chrom_all']),
         "hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_all_cell_types.h5",
         "train/{tf_name}/selected_motif_hdf5/{chrom}_motif_features_lightGBM.h5",
         # "hdf5s/DNase/median/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_{chrom}_median.h5",
         "train/data/quantile_transformer/{cell_line}_variable_bp_quantile_map.pkl",
         lambda wildcards: expand(
             "train/{{tf_name}}/models/{{tf_name}}.{training_cell_line}.{training_chrom_set_name}_model.pkl",
             # training_cell_line=selected_cell_line_for_reference(wildcards, 2),
             training_cell_line=cell_line_for_train(wildcards.tf_name),
             training_chrom_set_name=list(set(config['chrom_sets'].keys()) - {'chrom_set_test'}))
    output:
          'train/{tf_name}/predictions/{tf_name}.{cell_line}.{chrom}_preds.h5'
    threads: 1
    benchmark:
             "benchmarks/make_prediction_chrom/{tf_name}.{cell_line}.{chrom}.tsv"
    log:
       "logs/make_prediction_chrom/{tf_name}.{cell_line}.{chrom}.log"
    resources:
             mem_mb=30000
    params:
          runtime="2h",
          config_file="./config.yaml"
    shell:
         "python scripts/train_lightgbm_model.py make_prediction --cell_line={wildcards.cell_line} "
         "--dir_out='./train/{wildcards.tf_name}/predictions' --chrom={wildcards.chrom} "
         "--lightgbm_model_files_path='./train/{wildcards.tf_name}/models' "
         "--dir_dnase_feature_median='hdf5s/DNase/median' "
         "--config_file='{params.config_file}' "
         "--training_tf_name='{wildcards.tf_name}' "
         "--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' "
         "--motif_feature_path='./hdf5s/motif' "
         "--selected_motif_feature_path='train/{wildcards.tf_name}/selected_motif_hdf5/' "
         "--step=120 > {log} 2>&1"

rule evaluation:
    input:
         lambda wildcards: expand(
             "train/{{tf_name}}/predictions/{{tf_name}}.{{cell_line}}.{chrom}_preds.h5",
             chrom=config['chrom_all'])
    output:
          "train/{tf_name}/evaluations/{tf_name}.{cell_line}_performance.txt",
          "train/{tf_name}/evaluations/{tf_name}.{cell_line}_confusion_matrix.txt"
    threads: 1
    benchmark:
             "benchmarks/evaluation/{tf_name}.{cell_line}.tsv"
    log:
       "logs/evaluation/{tf_name}.{cell_line}.log"
    resources:
             mem_mb=20000
    params:
          runtime="2h",
          config_file="./config.yaml"
    priority: 50
    shell:
         "python scripts/train_lightgbm_model.py evaluation --cell_line={wildcards.cell_line} "
         "--dir_out='./train/{wildcards.tf_name}/evaluations' "
         "--lightgbm_preds_files_path='./train/{wildcards.tf_name}/predictions' "
         "--config_file='{params.config_file}' "
         "--training_tf_name='{wildcards.tf_name}' "
         "--quantile_transformer_path='./train/data/quantile_transformer' --dnase_feature_path='./hdf5s/DNase/' "
         "--motif_feature_path='./hdf5s/motif' "
         "--selected_motif_feature_path='train/{wildcards.tf_name}/selected_motif_hdf5/' "
         "--step=120 > {log} 2>&1"

