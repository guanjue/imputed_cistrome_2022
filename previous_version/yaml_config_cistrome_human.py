#!/usr/bin/env python3
import os
import yaml
import glob

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

config = dict(
              # the path to store all DNase or ATAC-seq bam files
              bam_file_path='/liulab/jfan/projects/impute_cistrome/cistrome/ATAC_bam',
              # the element in cell_types should be the sample name and there should be a
              # {cell_type}.bam in bam_file_path.
              # cell_types=['45056', '45057', '45058', '45059', '45060', '45061', '45062', '45063', '40968', '45065', '45066', '45067', '45068', '45069', '45070', '45071', '45072', '45073', '45074', '45075', '45076', '45077', '45078', '45079', '45080', '45081', '45082', '45083', '100380', '45085', '45086', '3138', '40966', '44942', '102440', '102441', '40967', '45064', '40919', '63547', '40970', '41022', '41023', '41024', '41025', '100418', '100419', '100420', '41029', '100422', '41031', '100424', '100425', '100426', '63563', '41036', '41037', '41038', '41039', '41040', '41041', '41042', '41043', '41044', '41045', '41046', '94295', '41048', '41049', '41050', '94299', '94300', '41053', '41054', '41055', '41056', '41057', '41058', '41059', '41060', '41061', '41062', '41063', '41064', '41065', '41066', '41067', '41068', '41069', '41070', '3091', '63604', '55174', '40853', '62484', '45013', '1801', '43138', '43139', '43141', '43142', '43143', '43144', '43145', '43146', '43147', '43148', '43149', '43150', '43151', '41104', '43153', '41106', '43155', '41108', '41109', '41110', '43159', '40884', '64197', '45084', '100381', '63668', '63669', '63670', '63671', '3103', '44947', '78197', '3104', '62497', '40849', '62542', '62498', '63698', '62499', '62468', '62500', '80090', '62501', '39138', '63716', '62502', '40925', '62503', '63724', '40516', '45017', '63735', '63736', '63752', '64557', '63768', '78113', '63781', '8488', '8489', '40518', '45019', '8495', '8496', '8498', '44951', '62467', '78138', '78139', '63807', '102723', '102724', '80197', '80198', '40810', '63818', '40792', '44952', '49079', '40929', '78182', '40793', '78195', '78196', '64574', '44943', '64577', '41026', '41027', '41028', '63906', '63907', '63908', '41030', '45023', '64527', '100334', '100423', '44955', '41032', '64585', '40692', '78266', '78268', '63562', '51646', '51647', '41035', '40694', '63944', '63945', '44956', '40695', '63949', '40696', '41429', '41430', '40697', '63964', '63965', '55186', '100127', '63969', '40699', '40511', '40797', '44957', '40700', '89511', '63985', '63986', '47604', '47606', '62552', '40934', '63575', '40866', '64598', '44958', '41047', '64014', '94296', '44944', '62553', '41051', '44959', '1798', '41052', '94763', '94764', '64049', '64050', '64055', '6712', '64057', '6715', '6716', '6717', '6718', '6719', '6720', '6721', '6722', '6723', '6724', '3169', '92746', '64075', '92748', '92749', '92750', '92752', '92753', '55190', '45029', '40733', '45040', '80486', '80487', '80488', '80489', '80490', '64110', '78447', '78448', '78449', '78450', '100979', '78452', '78453', '78454', '51645', '3090', '44972', '44894', '62625', '64141', '64142', '3181', '90770', '62489', '45031', '40803', '78492', '78495', '44904', '78500', '64167', '64168', '49064', '40734', '58042', '58043', '58044', '58045', '44896', '58050', '58051', '58052', '58053', '64198', '64199', '44975', '40941', '63266', '40949', '40873', '78627', '44965', '44897', '78995', '40601', '49076', '40847', '45000', '45034', '44966', '64252', '64253', '64254', '40406', '40943', '78636', '64983', '64276', '44967', '3204', '78981', '78625', '78626', '80675', '64586', '64608', '92972', '92973', '1800', '92975', '92978', '44968', '100148', '7647', '33596', '100149', '62603', '93000', '33612', '55189', '64337', '3215', '93020', '78686', '40421', '41105', '43154', '41107', '43156', '44960', '43157', '64384', '40427', '43158', '64392', '64393', '3167', '62354', '62355', '62356', '62357', '62358', '62359', '64408', '62361', '62362', '62363', '62364', '62365', '62366', '35743', '35744', '35745', '35746', '35747', '62372', '62373', '62374', '62375', '62376', '62377', '62379', '62380', '62381', '62382', '62383', '62384', '62385', '62386', '62387', '62388', '62389', '62390', '62391', '62392', '62393', '62394', '62395', '62396', '62397', '62398', '64447', '62400', '62401', '62402', '62403', '62404', '62405', '62406', '62407', '62408', '62409', '5215', '91639', '44905', '91640', '45042', '44974', '5216', '62444', '40786', '62446', '62447', '62448', '62450', '62451', '62452', '62453', '62454', '62455', '62456', '44884', '62458', '62459', '62460', '62461', '62462', '62463', '5120', '62465', '62466', '64515', '64516', '62469', '62470', '62471', '62472', '62473', '62475', '62476', '49067', '62478', '62479', '62480', '62481', '62482', '62483', '3092', '3093', '3094', '3095', '3096', '3097', '3098', '3099', '3100', '3101', '3102', '62495', '62496', '3105', '3106', '3107', '3108', '3109', '3110', '3111', '3112', '3113', '3114', '3115', '3116', '3117', '3118', '3119', '3120', '3121', '3122', '3123', '3124', '3125', '3126', '3127', '3128', '3129', '3131', '3132', '3133', '3134', '3135', '3136', '3137', '64578', '3139', '3140', '3141', '3142', '3143', '3144', '3145', '3146', '64587', '3148', '3149', '3150', '3151', '3152', '3153', '3154', '3155', '3156', '3157', '3158', '3159', '3160', '3161', '3162', '3163', '3164', '3165', '3166', '64607', '3168', '5217', '3170', '3171', '3172', '3173', '3174', '3175', '3176', '3177', '3178', '3179', '3180', '101485', '3182', '3183', '3184', '3185', '3186', '3187', '3188', '3189', '3190', '3191', '3192', '3193', '3194', '3195', '3196', '3197', '3198', '3199', '3200', '3201', '3202', '3203', '78980', '3205', '3206', '3207', '3208', '3209', '3210', '3211', '3212', '3213', '3214', '64655', '3216', '3217', '3218', '3219', '42138', '42139', '42140', '42141', '42142', '42143', '42144', '42145', '42146', '49068', '40818', '64687', '62611', '40819', '42165', '42166', '42167', '42168', '62370', '40886', '62664', '40816', '44935', '44921', '64729', '44976', '63054', '40826', '50400', '50401', '44923', '62695', '44924', '44925', '92099', '40830', '44908', '44927', '49075', '44928', '44929', '64965', '40834', '44931', '40815', '64793', '64795', '64796', '40837', '44964', '93474', '40838', '40839', '40840', '62773', '40841', '7480', '40842', '40843', '40885', '40894', '55180', '40515', '44941', '64056', '79429', '40846', '49078', '55183', '45033', '40963', '55184', '44936', '55185', '44946', '62831', '78403', '55187', '64062', '55188', '64063', '44973', '44949', '49080', '40854', '44883', '64905', '64906', '64407', '78402', '64910', '62360', '79427', '41033', '44953', '79428', '40858', '44613', '40859', '44614', '62886', '62887', '5020', '5021', '63048', '40693', '41034', '35742', '62903', '65039', '91578', '62367', '1890', '62911', '89536', '62368', '92747', '64964', '89850', '62369', '79497', '64969', '40395', '40396', '40397', '40398', '40399', '40401', '44963', '40405', '7638', '7639', '40408', '63267', '40410', '40411', '40412', '62941', '40414', '40415', '40416', '40417', '40418', '40419', '40420', '49062', '40422', '40423', '40426', '49063', '40428', '40429', '40430', '44907', '40434', '91635', '91636', '91637', '7670', '49065', '40440', '91641', '91642', '44937', '49066', '3147', '65024', '65025', '40875', '62982', '62983', '91638', '40458', '40948', '49069', '33603', '40468', '49070', '62999', '40474', '49071', '64565', '7713', '44993', '63012', '63013', '40881', '44948', '78994', '40882', '40494', '7727', '40498', '40499', '40500', '40501', '40502', '40503', '40504', '40505', '40506', '40507', '40508', '40509', '40510', '44981', '40512', '40513', '40514', '44611', '44612', '40517', '79430', '63047', '40520', '40521', '40522', '40523', '40524', '40525', '40526', '40527', '40528', '40529', '40530', '40531', '49081', '40890', '63070', '40559', '40560', '40561', '40562', '40563', '40564', '40565', '102006', '102007', '40568', '40569', '40836', '62399', '62464', '40574', '40575', '40577', '40582', '89735', '40584', '40585', '40586', '44994', '40593', '40594', '40595', '40596', '78446', '40598', '40599', '40600', '44996', '40602', '40603', '40604', '40605', '40606', '92101', '101659', '44998', '40614', '40615', '1898', '89734', '92104', '45032', '63164', '40566', '40646', '45004', '1899', '40865', '62457', '40910', '45008', '1900', '91883', '45010', '7922', '89843', '89844', '89845', '89846', '89847', '89848', '89849', '40698', '89851', '89852', '40701', '40702', '40703', '40704', '40705', '40706', '40707', '40708', '40709', '40710', '1799', '40712', '40713', '40714', '40715', '40716', '40717', '40718', '40719', '40720', '40721', '40722', '40723', '40724', '40725', '40726', '40727', '40728', '40729', '40730', '40731', '40732', '44898', '100126', '40735', '40736', '40737', '40738', '40739', '40740', '40741', '40742', '40743', '40744', '40745', '40746', '40747', '40748', '40749', '40750', '40751', '40752', '40753', '40754', '40755', '40756', '40757', '40758', '40759', '40760', '40813', '40762', '40763', '40764', '40765', '40766', '40767', '40768', '40928', '40772', '40773', '40774', '40775', '40776', '40777', '40783', '40784', '40785', '44882', '40787', '40788', '40789', '40790', '40791', '44888', '44889', '40794', '40795', '40796', '44893', '40798', '40799', '40800', '40801', '40802', '44899', '40804', '40805', '40806', '40807', '40808', '40809', '44906', '40811', '40812', '44909', '40814', '44911', '44912', '40817', '44914', '44915', '40820', '40821', '40822', '40823', '40824', '40825', '44922', '40827', '40828', '40829', '44926', '40831', '40832', '40833', '44930', '40835', '44932', '44933', '44934', '55175', '55176', '55177', '55178', '55179', '40844', '40845', '55182', '62445', '40848', '44945', '40850', '40851', '40852', '45038', '44950', '40855', '40856', '40857', '44954', '45039', '40860', '40861', '40862', '40863', '40864', '44961', '44962', '40867', '40868', '40869', '40870', '40871', '40872', '44969', '40874', '44971', '40876', '40877', '40878', '40879', '40880', '44977', '44978', '44979', '44980', '63413', '44982', '44983', '44984', '44985', '44986', '44987', '44988', '44989', '44990', '44991', '44992', '40897', '64579', '44995', '92100', '44997', '92102', '44999', '40904', '45001', '45002', '45003', '40908', '45005', '45006', '45007', '40711', '45009', '40914', '45011', '45012', '40917', '45014', '45015', '45016', '40921', '45018', '40923', '45020', '45021', '45022', '40927', '45024', '45025', '45026', '45027', '45028', '40933', '45030', '64849', '40936', '40937', '44940', '45035', '45036', '45037', '40942', '100335', '40944', '40946', '45043', '45044', '45045', '45046', '45047', '45048', '45049', '45050', '45051', '45052', '45053', '45054', '45055'],
              cell_types=['65259','80029','80143','66044','78923','79678','65871','65805'],
              # the reference_cell_type used as a template to construct all other feature dataset, choose one good
              # quality data.
              reference_cell_type='80029',
              # TFs to be trained, should be useless now...
              tf_names=['TCF7'],
              # all chrom names used in the project, excluding chrY for avoiding strange phenomenon
              chrom_all=list(map(lambda x: "chr" + str(x), list(range(1, 23)) + ["X"])),
              # the genome version
              genome="hg38",
              # chrom size file of the genome version, head 2 lines:
              # chr1	248956422
              # chr2	242193529
              chrom_size_file='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/hg38.chrom'
                              '.sizes',
              # the genome bins using for training and testing, here 200bp window with 50bp step is used. Head 2 lines:
              # chr1	0	200
              # chr1	50	250
              regions_all_file='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config'
                               '/hg38_regions_all_sorted.bed',
              # batch size used when precomputed the motif score across the genome.
              batch=10000000,
              # genome sequence for motif scan
              genome_sequence_fa='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/hg38'
                                 '.fa',
              # file path for motif files, example:
              # head /liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/motifs/PAX6_HUMAN.H11MO.0.C.pwm
              # >PAX6_HUMAN.H11MO.0.C
              # -0.6887818819947552	-1.975404324094267	-1.1886155146558977	1.1166043285356844
              motif_pwm_path='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/motifs/',
              # MDseqpos motif library
              motif_database_xml='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config'
                                 '/HOCOMOCOv11_full_pwm_HUMAN_mono.xml',
              # divide chrom into 3 sets: A, B for training and chrom_set_test is used for cross validation
              chrom_sets=dict(
                  chrom_setA=['chr2', 'chr4', 'chr6', 'chr7', 'chr12', 'chr13', 'chr15', 'chr16', 'chr17', 'chr20',
                              'chrX'],
                  chrom_setB=['chr3', 'chr5', 'chr9', 'chr10', 'chr11', 'chr14', 'chr18', 'chr19', 'chr22'],
                  chrom_set_test=["chr1", "chr8", "chr21"]),
              # file path for training label
              # example:
              # chr	start	stop	80029	65259	66044	79678	78923
              # chr1	0	200	U	U	U	U	U
              # chr1	50	250	U	U	U	U	U
              training_cell_types_regions_label_path="/liulab/jfan/projects/impute_cistrome/cistrome_impute_results_human_ATAC_short_long/label"
                                                     "/train",
              training_cell_types_regions_label_name="train.labels.tsv",
              test_cell_types_regions_label_path="/liulab/jfan/projects/impute_cistrome/cistrome_impute_results_human_ATAC_short_long/label"
                                                 "/test",
              test_cell_types_regions_label_name="test.labels.tsv",
              seqpos_env_path="/liulab/jfan/miniconda3/envs/seqpos",
              peak_file_path='/liulab/jfan/projects/impute_cistrome/cistrome/peak',
              # to be removed
              region_topic_model_h5='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/'
                                    'hg38_regions_topic_all_chroms.h5',
              # whether to merge the DNase-seq or ATAC-seq samples that belong to the same cluster for training and predicting
              CA_cluster_mode=False,
              # sample cluster result:
              # sample,cluster,Tissue,celltypes,celllines,included,annotation
              # 62388,0,Kidney,Clear Cell Carcinoma,Caki-2,Y,Renal Cell
              # 62387,0,Kidney,Clear Cell Carcinoma,Caki-2,Y,Renal Cell
              cluster_info_file="/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/"
                                "human_peaks_RP_DNaseseq_with_log_QN_hierarchical_sample_80_clusters_reannotated_missing_bad_samples_removed.csv",
              # ignore all other clusters but the selected cluster, leave blank to include all clusters
              selected_cluster=[],
              # TF ChIP-seq sample information
              TF_sample_info_file='/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/'
                                  'cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls',
              # 82348_sort_peaks.narrowPeak.bed
              cistromeDB_narrow_peak_path='/liulab/jfan/projects/impute_cistrome/cistrome/human_TF_peaks',
              giggle_path='/liulab/jfan/test/giggle/bin/giggle',
              giggle_index_path='/liulab/jfan/projects/impute_cistrome/cistrome/giggle_index/hg38_DNase_cluster_dup_removed_merged',
              # whether divide the reads based on the fragment size, this mode is designed for ATAC-seq data
              ATAC_long_short=True
              )


with open('config.yaml', 'w') as outfile:
    yaml.dump(config, outfile, Dumper=Dumper)

with open('config.yaml', 'r') as infile:
    print(yaml.load(infile, Loader=Loader))


# resource allocated to each type of job
cluster_config = dict(
    __default__=dict(
        job_name="default",
        time="1:00:00",
        n=1,
        mem_mb=1000,
        partition="defq"),
    TF_ChIP_sample_CA_cluster_matching_all=dict(
        job_name="TF_ChIP_sample_CA_cluster_matching_all",
        time="0:10:00",
        n=1,
        mem_mb=8000,
        partition="defq"),
    bwa_mapping=dict(
                    job_name="bwa_mapping",
                    time="12:00:00",
                    n=16,
                    mem_mb=60000,
                    partition="defq"),
    bam_filter=dict(
                    job_name="bam_filter",
                    time="4:00:00",
                    n=1,
                    mem_mb=16000,
                    partition="defq"),
    bam_sort=dict(
                job_name="bam_sort",
                time="2:00:00",
                n=10,
                mem_mb=10000,
                partition="defq"),
    bam_index=dict(
                job_name="bam_index",
                time="1:00:00",
                n=1,
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
        time="6:00:00",
        n=8,
        mem_mb=60000,
        partition="defq"),
    generate_quantile_transformer=dict(
        job_name="generate_quantile_transformer",
        time="6:00:00",
        n=2,
        mem_mb=20000,
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
    combine_TF_ChIP_sample_peaks_for_same_CA_cluster=dict(
        job_name="combine_TF_ChIP_sample_peaks_for_same_CA_cluster",
        time="1:00:00",
        n=1,
        mem_mb=16000,
        partition="defq"),
    prepare_label_file=dict(
        job_name="prepare_label_file",
        time="4:00:00",
        n=1,
        mem_mb=32000,
        partition="defq"),
    seqpos_scan_motifs_on_peaks=dict(
        job_name="seqpos_scan_motifs_on_peaks",
        time="10:00:00",
        n=1,
        mem_mb=10000,
        partition="defq"),
    prepare_motif_h5_data=dict(
        job_name="prepare_motif_h5_data",
        time="01:00:00",
        n=1,
        mem_mb=120000,
        partition="defq"),
    prepare_lightgbm_binary_data_motif_feature_reference=dict(
        job_name="prepare_lightgbm_binary_data_motif_feature_reference",
        time="1:00:00",
        n=1,
        mem_mb=60000,
        partition="defq"),
    prepare_lightgbm_binary_data_motif_feature_subset=dict(
        job_name="prepare_lightgbm_binary_data_motif_feature_subset",
        time="4:00:00",
        n=1,
        mem_mb=60000,
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
        mem_mb=60000,
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
        time="4:00:00",
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
