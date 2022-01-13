from sklearn.metrics import f1_score
from sklearn.metrics import average_precision_score
import lightgbm as lgb
import numpy as np
import pickle, os, sys


def logistic_obj(y_hat, dtrain, imbalance_alpha=None):
    if imbalance_alpha is None:
        imbalance_alpha = 5.0
    y = dtrain.get_label()
    #     p = y_hat
    p = 1. / (1. + np.exp(-y_hat))
    grad = (imbalance_alpha - 1) * p * y + p - imbalance_alpha * y
    hess = ((imbalance_alpha - 1) * y + 1) * (p * (1.0 - p))
    return grad, hess


def err_rate(y_hat, dtrain):
    y = dtrain.get_label()
    # y_hat = 1.0 / (1.0 + np.exp(-y_hat))
    y_hat = np.clip(y_hat, 10e-7, 1 - 10e-7)
    loss_fn = y * np.log(y_hat)
    loss_fp = (1.0 - y) * np.log(1.0 - y_hat)
    return 'error', np.sum(-(5 * loss_fn + loss_fp)) / len(y), False


def focal_isoform_binary_object(pred, dtrain, alpha=1, beta=0, gamma=2):
    #     alpha controls weight of positives
    #     (0,1) less
    #     >1 more
    #     beta controls the shift of loss function
    #     >0 to left(less weight to well-trained samples)
    #     gamma controls the steepness of loss function
    #     >0
    label = dtrain.get_label()
    x = beta + (2.0 * label - 1) * gamma * pred
    p = 1. / (1. + np.exp(-x))
    grad = (1 + (alpha - 1) * label) * (2 * label - 1) * (p - 1)
    hess = (1 + (alpha - 1) * label) * gamma * (1 - p) * p
    return grad, hess


# dir_out="/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/test"
dir_out = '/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/test/whole_genome'

dir_motif_feature = "/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/test"
dir_DNase_feature = '/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/DNASE/fold_coverage_wiggles/3M_features'


def lgb_f1_score(y_hat, data):
    y_true = data.get_label()
    y_hat = np.round(y_hat)  # scikits f1 doesn't like probabilities
    return 'f1', f1_score(y_true, y_hat), True


def lgb_auprc_score(y_hat, data):
    y_true = data.get_label()
    y_hat = np.round(y_hat)  # scikits f1 doesn't like probabilities
    return 'auprc', average_precision_score(y_true, y_hat), True


tf = "EGR1"

params = {
    'boosting_type': 'gbdt',

    #             'boosting_type': 'dart',
    #             'drop_rate':0.3,
    #             'max_drop':50,
    #             'skip_drop':0.5,
    #             'drop_seed':6,
    # 'pos_bagging_fraction': 1,
    # 'neg_bagging_fraction': 0.01,
    # 'bagging_freq': 10000,
    # 'bagging_seed': 6,
    'objective': 'binary',
    #             'objective': focal_binary_object,
    #             'metric': ['binary_error','binary_logloss',"auc"],
    'metric': ["auc"],
    #             'is_unbalance': True,
    #             "scale_pos_weight": 100,
    'metric_freq': 10,
    'num_leaves': 63,
    #             'max_bin':255,
    'num_threads': 16,
    'learning_rate': 0.1,
    'feature_fraction': 1,
    'boost_from_average': False,
    #             'early_stopping_round':20,

    'verbose': 1
}

import gc

cell_line, set_type = 'HCT116', 'setA'

other_set_type = list(set(["setA", "setB"]) - set([set_type]))[0]
other_cell_lines = list(set(['GM12878', 'HCT116', 'H1-hESC', 'MCF-7']) - set([cell_line]))
#         train_data=lgb.Dataset("%s/%s_%s_%s_lightGBM.bin" % (dir_out,tf,cell_line,set_type))
train_data = None
gc.collect()
train_data = lgb.Dataset("%s/%s_%s_%s_all_quantile_variable_bp_lightGBM.bin" %
                         (dir_out, tf, cell_line, set_type)).construct()

print
cell_line, set_type
#         train_data=train_data.construct()
list_validation_data = []
gc.collect()
for other_cell_line in other_cell_lines:
    #             validation_data = lgb.Dataset("%s/%s_%s_%s_lightGBM.bin" % \
    #                                           (dir_out,tf,other_cell_line,other_set_type), reference=train_data)
    #             validation_data=lgb.Dataset("%s/%s_%s_%s_all_lightGBM.bin" %
    #                                    (dir_out, tf, other_cell_line, other_set_type), reference=train_data)
    validation_data = None
    gc.collect()
    validation_data = lgb.Dataset("%s/%s_%s_%s_all_quantile_variable_bp_lightGBM.bin" %
                                  (dir_out, tf, other_cell_line, "set_test"), reference=train_data).construct()
    #             validation_data = train_data.create_valid("%s/%s_%s_%s_lightGBM.bin" % \
    #                                           (dir_out,tf,other_cell_line,other_set_type))
    list_validation_data.append(validation_data)
    print
    other_cell_line, other_set_type

evals_result = {}
print
"training starts"
gbm = lgb.train(params=params,
                train_set=train_data,
                fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=1, beta=-1.5, gamma=1.01),
                #                 fobj=lambda x,y:logistic_obj(x,y,imbalance_alpha=1.0),
                valid_sets=[train_data] + list_validation_data,
                valid_names=['train'] + ["%s_%s" % (other_cell_line, other_set_type) \
                                         for other_cell_line in other_cell_lines],
                feval=lgb_auprc_score,
                early_stopping_rounds=20,
                evals_result=evals_result,
                keep_training_booster=True)

with open("%s/%s_HCT116_test_model_genome_variable_bp_with_forward_strand.pkl" % (dir_out, tf), 'wb') as fout:
    pickle.dump(gbm, fout)
#         # load model with pickle to predict
#         with open('model.pkl', 'rb') as fin:
#             pkl_bst = pickle.load(fin)
#         with open("%s/%s_%s_%s_evals_result.pkl" % (dir_out,tf,cell_line,set_type), 'wb') as outfile_evals_result:
with open("%s/%s_%s_%s_evals_result_genome_variable_bp.pkl" % (dir_out, tf, cell_line, set_type),
          'wb') as outfile_evals_result:
    pickle.dump(evals_result, outfile_evals_result, pickle.HIGHEST_PROTOCOL)

import pandas as pd
import h5py

data_dir = "/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM"

tf = "EGR1"
# 4 cell types, 20 chroms
train_path_file = "ChIPseq/labels/%s.train.labels.tsv.gz" % tf

df_data = pd.read_table("%s/%s" % (data_dir, train_path_file), sep="\t", header=0)

cell_line = "HCT116"
chrom_setA = ['chr2', 'chr4', 'chr6', 'chr7', 'chr12', 'chr13', 'chr15', 'chr16', 'chr17', 'chr20', 'chrX']
chrom_setB = ['chr3', 'chr5', 'chr9', 'chr10', 'chr11', 'chr14', 'chr18', 'chr19', 'chr22']
chrom_set_test = ["chr1", "chr8", "chr21"]

dir_out = '/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/test/whole_genome'

dir_motif_feature = "/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/test"
dir_DNase_feature = '/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/DNASE/fold_coverage_wiggles/3M_features'

model_file = "%s/%s_HCT116_test_model_genome_variable_bp_with_forward_strand.pkl" % (dir_out, tf)
with open(model_file, 'rb') as fin:
    gbm = pickle.load(fin)

new_train_data = []
for chrom in chrom_setA + chrom_setB:
    print
    chrom,
    df_temp = df_data.loc[df_data['chr'] == chrom, :]
    feature_names = []
    dnase_file_name = "%s/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_%s.h5" % (dir_DNase_feature, chrom)
    with h5py.File(dnase_file_name, "r") as infile:
        celline_index = np.where(infile['samples'][...] == cell_line)[0][0]
        scores_DNase = infile[chrom][celline_index, :, :]
        feature_names += list(infile['feature_names'][...])
        with open("%s/../%s_variable_bp_quantile_map.pkl" %
                  (dir_out, cell_line), 'rb') as fin:
            qt = pickle.load(fin)
        _ = qt.transform(scores_DNase)
    with h5py.File("%s/%s_motif_features_lightGBM.h5" % (dir_motif_feature, chrom), "r") as infile:
        starts = infile["starts"][...]
        scores_motif = infile['scores'][...]
        feature_names += list(infile['feature_names'][...])
    #         df_lightGBM_2=pd.DataFrame(scores,columns=infile['feature_names'][...])
    #     df_lightGBM_all=pd.concat([df_lightGBM,df_lightGBM_2],axis=1)
    scores = np.hstack((scores_DNase, scores_motif))
    ignore_index = np.where(df_temp[cell_line] == "A")[0]
    scores = np.delete(scores, ignore_index, axis=0)
    label_B_U = np.delete(np.array(df_temp[cell_line]), ignore_index, axis=0)
    starts = np.delete(np.array(starts), ignore_index, axis=0)
    temp_label = map(lambda x: 1 if x == "B" else 0, label_B_U)
    temp_label = np.array(temp_label)
    preds = None
    preds = np.zeros((len(temp_label), 1))
    params = {
        'num_threads': 8,
    }
    gbm = gbm.reset_parameter(params)
    ypred = gbm.predict(scores)
    preds[:, 0] = ypred
    preds = 1. / (1. + np.exp(-preds))
    preds = preds[:, 0]
    fp_index = np.where(np.logical_and(temp_label == 0, preds > 0.5))[0]
    p_index = np.where(temp_label == 1)[0]
    selected_index = np.hstack((fp_index, p_index))
    new_train_data.append(scores[selected_index, :])

data = np.hstack(new_train_data)

cell_line = "HCT116"
chrom_setA = ['chr2', 'chr4', 'chr6', 'chr7', 'chr12', 'chr13', 'chr15', 'chr16', 'chr17', 'chr20', 'chrX']
chrom_setB = ['chr3', 'chr5', 'chr9', 'chr10', 'chr11', 'chr14', 'chr18', 'chr19', 'chr22']
chrom_set_test = ["chr1", "chr8", "chr21"]

dir_out = '/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/test/whole_genome'

dir_motif_feature = "/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/test"
dir_DNase_feature = '/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/DNASE/fold_coverage_wiggles/3M_features'

temp = []
temp_2 = []
for chrom in chrom_setA + chrom_setB:
    temp.append(np.load("%s/retrain_data_%s.npy" % (dir_out, chrom)))
    temp_2.append(np.load("%s/retrain_label_%s.npy" % (dir_out, chrom)))

data = np.vstack(temp)
label = np.hstack(temp_2)
model_file = "%s/%s_HCT116_test_model_genome_variable_bp_with_forward_strand.pkl" % (dir_out, tf)
with open(model_file, 'rb') as fin:
    gbm = pickle.load(fin)

tf = "EGR1"

params = {
    'boosting_type': 'gbdt',

    #             'boosting_type': 'dart',
    #             'drop_rate':0.3,
    #             'max_drop':50,
    #             'skip_drop':0.5,
    #             'drop_seed':6,
    # 'pos_bagging_fraction': 1,
    # 'neg_bagging_fraction': 0.01,
    # 'bagging_freq': 10000,
    # 'bagging_seed': 6,
    'objective': 'binary',
    #             'objective': focal_binary_object,
    #             'metric': ['binary_error','binary_logloss',"auc"],
    'metric': ["auc"],
    #             'is_unbalance': True,
    #             "scale_pos_weight": 100,
    'metric_freq': 10,
    'num_leaves': 63,
    #             'max_bin':255,
    'num_threads': 16,
    'learning_rate': 0.1,
    'feature_fraction': 1,
    'boost_from_average': False,
    #             'early_stopping_round':20,

    'verbose': 1
}

import gc

cell_line, set_type = 'HCT116', 'setA'

other_set_type = list(set(["setA", "setB"]) - set([set_type]))[0]
other_cell_lines = list(set(['GM12878', 'HCT116', 'H1-hESC', 'MCF-7']) - set([cell_line]))
#         train_data=lgb.Dataset("%s/%s_%s_%s_lightGBM.bin" % (dir_out,tf,cell_line,set_type))
train_data = None
gc.collect()

#         train_data=train_data.construct()
list_validation_data = []
gc.collect()

init_score = gbm.predict(data)
train_data = lgb.Dataset(data, label, init_score=init_score, free_raw_data=False, reference=validation_data).construct()

for other_cell_line in other_cell_lines:
    #             validation_data = lgb.Dataset("%s/%s_%s_%s_lightGBM.bin" % \
    #                                           (dir_out,tf,other_cell_line,other_set_type), reference=train_data)
    #             validation_data=lgb.Dataset("%s/%s_%s_%s_all_lightGBM.bin" %
    #                                    (dir_out, tf, other_cell_line, other_set_type), reference=train_data)
    validation_data = None
    gc.collect()
    validation_data = lgb.Dataset("%s/%s_%s_%s_all_quantile_variable_bp_lightGBM.bin" %
                                  (dir_out, tf, other_cell_line, "set_test"), reference=train_data).construct()
    #             validation_data = train_data.create_valid("%s/%s_%s_%s_lightGBM.bin" % \
    #                                           (dir_out,tf,other_cell_line,other_set_type))
    list_validation_data.append(validation_data)
    print
    other_cell_line, other_set_type

evals_result = {}
print
"training starts"

gbm_new = lgb.train(params=params,
                    train_set=train_data,
                    fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=1, beta=0, gamma=1),
                    #                 fobj=lambda x,y:logistic_obj(x,y,imbalance_alpha=1.0),
                    valid_sets=[train_data] + list_validation_data,
                    valid_names=['train'] + ["%s_%s" % (other_cell_line, other_set_type) \
                                             for other_cell_line in other_cell_lines],
                    feval=lgb_auprc_score,
                    early_stopping_rounds=20,
                    evals_result=evals_result,
                    keep_training_booster=True,
                    # init_model=gbm
                    )

gbm_new = lgb.train(params=params,
                    train_set=train_data,
                    fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=1, beta=0, gamma=1),
                    #                 fobj=lambda x,y:logistic_obj(x,y,imbalance_alpha=1.0),
                    # valid_sets=[train_data]  ,
                    # valid_names=['train'] ,
                    feval=lgb_auprc_score,
                    # early_stopping_rounds=20,
                    evals_result=evals_result,
                    keep_training_booster=True,
                    init_model=gbm)

with open("%s/%s_HCT116_retrained_model_genome_variable_bp_with_forward_strand.pkl" % (dir_out, tf), 'wb') as fout:
    pickle.dump(gbm_new, fout)

import numpy as np

cell_line = "liver"
chrom = "chr22"
labels = np.array(df_test_regions_label.loc[df_test_regions_label['chr'] == chrom, :][cell_line])

import pandas as pd

training_tf_name = "EGR1"
if dir_out is None:
    dir_out = "./train/%s/evaluations/" % training_tf_name
df_test_regions_label = pd.read_csv(
    "%s/%s.%s" % (
        '/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/test/pipeline_py3/label/test',
        training_tf_name,
        'train.labels.tsv.gz'), sep="\t", header=0)
list_preds_binary = []
# list_preds_binary_2 = []
list_labels = []
list_preds_matrix = []
for chrom in list(map(lambda x: "chr" + str(x), list(range(1, 23)) + ["X"])):
    with h5py.File(
            '%s/%s.%s.%s_preds.h5' % ('./train/EGR1/predictions', training_tf_name, cell_line, chrom),
            "r") as infile:
        model_files = infile['model_files'][...]
        preds = infile['preds'][...]
        labels = np.array(df_test_regions_label.loc[df_test_regions_label['chr'] == chrom, :][cell_line])
        list_preds_matrix.append(preds)
        preds_binary = np.mean(1. / (1. + np.exp(-preds)), axis=1)
        # preds_binary_2 = 1. / (1. + np.exp(-np.mean(preds, axis=1)))
        list_preds_binary.append(preds_binary)
        # list_preds_binary_2.append(preds_binary_2)
        list_labels.append(labels)
labels = np.hstack(list_labels)
preds = np.hstack(list_preds_binary)
# preds_2 = np.hstack(list_preds_binary_2)
preds_matrix = np.vstack(list_preds_matrix)
ignore_index = np.where(labels == "A")[0]
preds_matrix = np.delete(preds_matrix, ignore_index, axis=0)
preds = np.delete(preds, ignore_index, axis=0)
label_b_u = np.delete(labels, ignore_index, axis=0)
label_b_u = np.array(map(lambda x: 1 if x == "B" else 0, label_b_u))
print(label_b_u.shape)
print(label_b_u[:10])
print(preds.shape)
print(preds[:10])
with open("%s/%s.%s_performance.txt" % (dir_out, self.training_tf_name, cell_line), "w") as outfile:
    fpr, tpr, thresholds = metrics.roc_curve(label_b_u, preds, pos_label=1)
    auc = metrics.auc(fpr, tpr)
    auprc = average_precision_score(label_b_u, preds)
    outfile.write("average model: auc:%.6f auprc:%.6f\n" % (auc, auprc))
    for i in range(preds_matrix.shape[1]):
        fpr, tpr, thresholds = metrics.roc_curve(label_b_u, preds_matrix[:, i], pos_label=1)
        auc = metrics.auc(fpr, tpr)
        auprc = average_precision_score(label_b_u, preds_matrix[:, i])
        outfile.write("%s model: auc:%.6f auprc:%.6f\n" % (model_files, auc, auprc))

import numpy as np
import lightgbm as lgb
import pickle

a = lgb.Dataset('lightGBM.all.MCF-7.chrom_set_test.bin')
a = lgb.Dataset('train/data/dnase_feature_binary_files/lightGBM.dnase.GM12878.chrom_setA.bin')

with open("train/EGR1/models/EGR1.HCT116.chrom_setA_model.pkl", 'rb') as fin:
    gbm = pickle.load(fin, encoding='latin1')

data = np.random.random((500, 10))
y = [1] * 250 + [0] * 250
lgb_train_reference = lgb.Dataset(data, y, free_raw_data=True)
lgb_train_reference.save_binary("lgb_train_data_reference.bin")

data = np.random.random((500, 10))
y = [1] * 250 + [0] * 250
lgb_train_1 = lgb.Dataset(data, y, free_raw_data=True, feature_name=list(map(lambda x: "c" + str(x), range(10))),
                          reference=lgb_train_reference)
lgb_train_1.save_binary("lgb_train_data_1.bin")

data = np.random.random((500, 10))
y = [1] * 250 + [0] * 250
lgb_train_2 = lgb.Dataset(data, y, free_raw_data=True, reference=lgb_train_reference)
lgb_train_2.save_binary("lgb_train_data_2.bin")

lgb_train_from_file_1 = lgb.Dataset('lgb_train_data_1.bin', params={'ignore_column': 'name:c0,c1'}
                                    )
lgb_train_from_file_1 = lgb.Dataset('lgb_train_data_1.txt', params={'ignore_column': '0,1'}
                                    )
lgb_train_from_file_2 = lgb.Dataset('lgb_train_data_2.bin', reference=lgb_train_from_file_1)

subset_index_1 = np.random.choice(np.arange(500), 300, replace=False)
subset_data_1 = lgb_train_from_file_1.subset(subset_index_1)
subset_index_2 = np.random.choice(np.arange(500), 200, replace=False)
subset_data_2 = lgb_train_from_file_2.subset(subset_index_2)
params = {
    # 'ignore_column': '0,1',
    'objective': 'binary',
    'verbose': 1,
    'learning_rates':10,
    'metric':"auc"
}

gbm = lgb.train(params=params,
                train_set=lgb_train_reference,
                valid_sets=[lgb_train_reference],
                num_boost_round=100,
                early_stopping_rounds=1
                )

a = gbm.predict(data[:4, :], raw_score=True, pred_leaf=True, pred_contrib=False)
b = gbm.predict(data[:4, :], raw_score=True, pred_leaf=False, pred_contrib=False)

temp = np.zeros(a.shape)
for i in range(a.shape[0]):
    for j in range(a.shape[1]):
        temp[i, j] = gbm.get_leaf_output(j, a[i, j])

gbm = lgb.train(params=params,
                train_set=subset_data_1,
                valid_sets=[subset_data_1, subset_data_2],
                num_boost_round=2,
                )

gbm = lgb.train(params=params,
                train_set=lgb_train_from_file_1,
                # valid_sets=[lgb_train_from_file_1, lgb_train_from_file_2],
                num_boost_round=2,
                )

# def merge_lightgbm_binary_data(training_tf_name, cell_line, chrom_set_name, step=120, lightgbm_dnase_binary_files_path=None,
# #                                lightgbm_motif_binary_files_path=None):

import numpy as np
import lightgbm as lgb
import pandas as pd
import h5py

training_tf_name = 'EGR1'
cell_line = 'GM12878'
chrom_set_name = 'chrom_set_test'
step = 120
lightgbm_dnase_binary_files_path = None
lightgbm_motif_binary_files_path = None
if lightgbm_motif_binary_files_path is None:
    lightgbm_motif_binary_files_path = "./train/%s/binary_files" % training_tf_name

if lightgbm_dnase_binary_files_path is None:
    lightgbm_dnase_binary_files_path = "./train/data/dnase_feature_binary_files"

all_feature_names = []
chrom = "chr22"
selected_motif_feature_path = 'train/%s/selected_motif_hdf5/' % training_tf_name
# TODO change to 50bp or 100bp
with h5py.File("%s/%s_motif_features_lightGBM.h5" % (selected_motif_feature_path, chrom),
               "r") as infile:
    all_feature_names += list(infile['feature_names'][...])

all_feature_names = list(map(lambda x: x.decode('UTF-8'), all_feature_names))
train_data_all = lgb.Dataset("%s/lightGBM.dnase.%s.%s.bin" %
                             (lightgbm_dnase_binary_files_path, cell_line, chrom_set_name)).construct()
for subset_index in range(int(np.ceil(len(all_feature_names) / step))):
    train_data = lgb.Dataset("%s/lightGBM.motif.%s.%d.bin" %
                             (lightgbm_motif_binary_files_path, chrom_set_name, subset_index + 1)).construct()
    train_data_all.add_features_from(train_data)

temp = []
df_all_regions_label = pd.read_csv(
    "%s/%s.%s" % (
        "/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM/test/pipeline_py3/label/train",
        training_tf_name,
        'train.labels.tsv'),
    sep="\t", header=0)
chrom_sets = dict(
    chrom_setA=['chr2', 'chr4', 'chr6', 'chr7', 'chr12', 'chr13', 'chr15', 'chr16', 'chr17', 'chr20',
                'chrX'],
    chrom_setB=['chr3', 'chr5', 'chr9', 'chr10', 'chr11', 'chr14', 'chr18', 'chr19', 'chr22'],
    chrom_set_test=["chr1", "chr8", "chr21"])
chrom_set = chrom_sets[chrom_set_name]
temp = []
for chrom in chrom_set:
    df_temp = df_all_regions_label.loc[df_all_regions_label['chr'] == chrom, :]
    temp.append(df_temp)

df_all_temp = pd.concat(temp, ignore_index=True)
# selected_index = np.where(df_all_temp[cell_line] != "A")[0]
# ignore_index = np.where(df_all_temp[cell_line] == "A")[0]
# label_b_u = np.delete(np.array(df_all_temp[cell_line]), ignore_index, axis=0)
# labels = list(map(lambda x: 1 if x == "B" else 0, label_b_u))
# train_data_all_subset = train_data_all.subset(selected_index)
# train_data_all_subset.set_label(labels)
# return train_data_all_subset
weight = (np.array(df_all_temp[cell_line]) != "A").astype(int)
train_data_all.set_weight(weight)
labels = (np.array(df_all_temp[cell_line]) == "B").astype(int)
train_data_all.set_label(labels)
train_data_all.save_binary('2.bin')

lgb_train_from_file_1 = lgb.Dataset('1.bin')
lgb_train_from_file_2 = lgb.Dataset('2.bin')

gbm = lgb.train(params=params,
                train_set=lgb_train_from_file_1,
                valid_sets=[lgb_train_from_file_1, lgb_train_from_file_2],
                num_boost_round=2,
                )

import os
import subprocess


def addPy2Paths_Config(config):
    """ADDS the python2 paths to config"""
    conda_root = subprocess.check_output('conda info --root', shell=True).decode('utf-8').strip()
    conda_path = os.path.join(conda_root, 'pkgs')
    config["python2_pythonpath"] = os.path.join(conda_root, 'envs', 'chips_py2', 'lib', 'python2.7', 'site-packages')
    if not "python2" in config or not config["python2"]:
        config["python2"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'python2.7')
    if not "mdseqpos_path" in config or not config["mdseqpos_path"]:
        config["mdseqpos_path"] = os.path.join(conda_root, 'envs', 'chips_py2', 'bin', 'MDSeqPos.py')


config = {}
addPy2Paths_Config(config)

/ n / xiaoleliu_lab / chips / miniconda3 / envs / chips_py2
/ n / xiaoleliu_lab / chips / miniconda3 / envs / chips_py2 / bin / MDSeqPos.py
/ n / xiaoleliu_lab / chips / miniconda3 / envs / chips_py2 / bin / python2
.7

a = h5py.File(
    '/n/scratchlfs/xiaoleliu_lab/cchen/Cistrome_imputation/encode/data/DNase_scanning/scan_result/DNASE.PC-3.merge.binSize.1.corrected_sorted_hg19_25bpbin_bwaverage_transformed_chr5_scanned_with_autoencoder.hdf5',
    "r")
b = h5py.File("hdf5s/DNase/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_chr5_all_cell_types.h5", "r")
import h5py
import glob

for filename in glob.glob("train/REST/selected_motif_hdf5/*.h5"):
    print(filename)
    a = h5py.File(filename, "r")
    print(a.keys())

import os
import glob

for filename in glob.glob("bams/with_filter/DNASE.*.merge.filter.bam"):
    print(filename)
    a = os.popen("{ samtools view -H %s ; samtools view %s | head -n 1000; } | samtools view -c -f 1" % (
    filename, filename)).read().strip()
    print(a)

import lightgbm as lgb
import glob
import os

for filename in glob.glob("train/HNF4A/binary_files/lightGBM.all*"):
    print(filename)
    try:
        a = lgb.Dataset(filename).construct()
    except:
        print("error")
        os.remove(filename)


a = lgb.Dataset("train/data/dnase_feature_binary_files/lightGBM.autoencoder.dnase.HepG2.chrom_setA.bin").construct()
b = lgb.Dataset("train/data/dnase_feature_binary_files/lightGBM.dnase.A549.chrom_setA.bin").construct()
a.num_data()
b.num_data()

/ n / xiaoleliu_lab / chips / miniconda3 / envs / chips_py2 / bin / python2
.7 / n / xiaoleliu_lab / chips / miniconda3 / envs / chips_py2 / bin / MDSeqPos.py \
    peaks / top5k / ChIPseq.HepG2.FOXA1.conservative.train.top5k.narrowPeak
hg19 \
- m / n / scratchlfs / xiaoleliu_lab / Jingyu / impute_cistrome / ENCODE_DREAM / annotations / HOCOMOCOv11_full_pwm_HUMAN_mono.xml \
- d - O. / train / FOXA1 / motif_scan / cell_line.HepG2.tf.FOXA1

shell(
    "{params.pypath} {config[mdseqpos_path]} {input} {params.genome} -m cistrome.xml -d -O analysis/motif/{params.runName}/results 1>>{log}")

import gc
import os
import pickle

import fire
import h5py

from sklearn.preprocessing import StandardScaler
from numpy import linalg as LA
import sklearn.metrics as metrics
import json
import lightgbm as lgb

import numpy as np

import pandas as pd
import glob
from sklearn.preprocessing import QuantileTransformer
import yaml

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from sklearn.metrics import average_precision_score


def focal_isoform_binary_object(pred, dtrain, alpha=1.0, beta=0.0, gamma=2.0):
    #     alpha controls weight of positives
    #     (0,1) less
    #     >1 more
    #     beta controls the shift of loss function
    #     >0 to left(less weight to well-trained samples)
    #     gamma controls the steepness of loss function
    #     >0
    label = dtrain.get_label()
    x = beta + (2.0 * label - 1) * gamma * pred
    p = 1. / (1. + np.exp(-x))
    grad = (1 + (alpha - 1) * label) * (2 * label - 1) * (p - 1)
    hess = (1 + (alpha - 1) * label) * gamma * (1 - p) * p
    return grad, hess


def lgb_auprc_score(y_hat, data):
    y_true = data.get_label()
    # TODO try not to round yhat
    y_hat = np.round(y_hat)  # scikits f1 doesn't like probabilities
    return 'auprc', average_precision_score(y_true, y_hat), True


training_tf_name = "EGR1"

lightgbm_motif_binary_files_path = "./train/%s/binary_files" % training_tf_name

lightgbm_dnase_binary_files_path = "./train/data/dnase_feature_binary_files"

dir_out = "./train/%s/models/" % training_tf_name

params = {
    'boosting_type': 'gbdt',
    # 'ignore_column': list(range(500)),
    # 'boosting_type': 'dart',
    # 'drop_rate': 0.3,
    # 'max_drop': 50,
    # 'skip_drop': 0.5,
    # 'drop_seed': 6,
    # 'pos_bagging_fraction': 1,
    # 'neg_bagging_fraction': 0.01,
    # 'bagging_freq': 10000,
    # 'bagging_seed': 6,
    'objective': 'binary',
    # 'objective': focal_binary_object,
    # 'metric': ['binary_error', 'binary_logloss', "auc"],
    'metric': ["auc"],
    # 'is_unbalance': True,
    # "scale_pos_weight": 100,
    'metric_freq': 10,
    'num_leaves': 31,
    'min_data_in_leaf': 20,
    # 'max_bin': 255,
    'num_threads': 32,
    'learning_rate': 0.01,
    'feature_fraction': 1,
    'boost_from_average': False,
    'verbose': 1
}

chrom_set_name = 'chrom_setA'
other_cell_lines = ['MCF-7', 'H1-hESC', 'GM12878']
cell_line = 'HCT116'
train_data = lgb.Dataset(
    "%s/lightGBM.all.%s.%s.bin" % (lightgbm_motif_binary_files_path, cell_line, chrom_set_name))
list_validation_data = []
for other_cell_line in other_cell_lines:
    # validation_data = self.merge_lightgbm_binary_data(other_cell_line, "chrom_set_test",
    #                                                   lightgbm_dnase_binary_files_path,
    #                                                   lightgbm_motif_binary_files_path)
    validation_data = lgb.Dataset("%s/lightGBM.all.%s.%s.bin" % (
        lightgbm_motif_binary_files_path, other_cell_line, "chrom_set_test"), reference=train_data)
    list_validation_data.append(validation_data)
os.remove('train/REST/binary_files/lightGBM.all.HeLa-S3.chrom_set_test.bin')

with open("train/EGR1/models/EGR1.HCT116.chrom_setA_model.pkl", 'rb') as fin:
    init_gbm = pickle.load(fin, encoding='latin1')

evals_result = {}
gbm = lgb.train(params=params,
                  train_set=train_data,
                  num_boost_round=1000,
                  fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=1, beta=-1.5, gamma=1.01),
                  #                 fobj=lambda x,y:logistic_obj(x,y,imbalance_alpha=1.0),
                  valid_sets=[train_data] + list_validation_data,
                  valid_names=['train'] + ["%s_%s" % (other_cell_line, "set_test") \
                                           for other_cell_line in other_cell_lines],
                  feval=lgb_auprc_score,
                  early_stopping_rounds=20,
                  evals_result=evals_result,
                  keep_training_booster=False,
                  learning_rates=lambda x: 0.2 * (0.98 ** x),
                  # init_model=gbm
                  # init_model=init_gbm
                  )

with open("%s/%s.%s.%s_model.pkl" % (
        dir_out, self.training_tf_name, cell_line, chrom_set_name), 'wb') as fout:
    pickle.dump(gbm, fout)
with open("%s/%s.%s.%s_evals_result.pkl" % (
        dir_out, self.training_tf_name, cell_line, chrom_set_name),
          'wb') as outfile_evals_result:
    pickle.dump(evals_result, outfile_evals_result, pickle.HIGHEST_PROTOCOL)




import os, torch, pybedtools


def chrom_loc2motif_scores(fasta_instance, chrom, start, end, strand, motif_pwms, num_threads=1):
    torch.set_num_threads(num_threads)
    sequence = pybedtools.BedTool.seq((chrom,
                                       max(start - int(motif_pwms.shape[1] / 2), 0),
                                       end + motif_pwms.shape[1] + 200 - 50 - 1 - int(motif_pwms.shape[1] / 2)
                                       ), fasta_instance)
    sequence = sequence.replace("a", "A")
    sequence = sequence.replace("c", "C")
    sequence = sequence.replace("g", "G")
    sequence = sequence.replace("t", "T")
    sequence = sequence.replace("N", "A")
    sequence = sequence.replace("n", "A")
    if strand == "+":
        sequence = sequence.replace("A", "0")
        sequence = sequence.replace("C", "1")
        sequence = sequence.replace("G", "2")
        sequence = sequence.replace("T", "3")
    else:
        motif_pwms = np.flip(motif_pwms, axis=1).copy()
        sequence = sequence.replace("A", "3")
        sequence = sequence.replace("C", "2")
        sequence = sequence.replace("G", "1")
        sequence = sequence.replace("T", "0")
    sequence = np.asarray([0] * max(int(motif_pwms.shape[1] / 2) - start, 0) + list(sequence), dtype='int')
    sequence_one_hot = torch.nn.functional.one_hot(torch.from_numpy(sequence), num_classes=4)
    sequence_one_hot = sequence_one_hot.reshape((1, 1, -1, 4)).float()
    motif_filters = torch.from_numpy(motif_pwms.reshape((motif_pwms.shape[0], 1, motif_pwms.shape[1], 4))).float()
    scores = torch.nn.functional.conv2d(sequence_one_hot, motif_filters, padding=0, stride=1)
    #     scores=scores.reshape((motif_pwms.shape[0], -1)).numpy()
    scores = scores.reshape((motif_pwms.shape[0], -1))
    return scores


import numpy as np
import os




tf_length_max = 0
motif_names = []
motif_pwms = []
motif_lengths = []

for motif_file in os.listdir('./motifs'):
    tf_matrix = np.loadtxt("./motifs/"+motif_file, comments=["#", ">"], delimiter='\t',
                           dtype='float32')
    tf_length = tf_matrix.shape[0]
    motif_lengths.append(tf_length)
    if tf_length_max < tf_length:
        tf_length_max = tf_length
    if motif_file.endswith('.pwm'):
        tf_name = motif_file.replace('.pwm', '')
    else:
        tf_name = motif_file
    motif_names.append(tf_name)
    motif_pwms.append(tf_matrix)


for ind, pwm in enumerate(motif_pwms):
    zero_padding_length = tf_length_max - pwm.shape[0]
    motif_pwms[ind] = np.vstack((
        np.zeros((int(zero_padding_length / 2), 4)),
        pwm,
        np.zeros((zero_padding_length - int(zero_padding_length / 2), 4))
    ))

motif_pwms = np.array(motif_pwms)




fasta_instance = pybedtools.example_filename('/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/hg38.fa')


motif_scores = chrom_loc2motif_scores(fasta_instance, 'chr1', 0,
                                                           20000, '+',
                                                           motif_pwms, 1)

motif_scores_unfold = motif_scores.unfold(dimension=1, size=200, step=50)

def torch_topk(tensor, k, dim, largest=None, sorted=None):
    if largest is None:
        largest = True
    if sorted is None:
        sorted = True
    return torch.topk(tensor, k, dim, largest, sorted)[0]

from functools import partial

torch_topk_new = partial(torch_topk, k=3, dim=2, largest=True, sorted=True)

a=torch_topk_new(motif_scores_unfold)



with open('/liulab/jfan/projects/impute_cistrome/cistrome/cistrome_impute_hg38_config/hg38_regions_all.bed','w') as outfile:
    for chrom in chrom_all:
        for i in range(0, dic_chrom_length[chrom]-200,50):
            outfile.write("%s\t%d\t%d\n" % (chrom, i, i+200))
