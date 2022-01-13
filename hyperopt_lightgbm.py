from hyperopt import STATUS_OK
from hyperopt import hp
from timeit import default_timer as timer
import numpy as np
import lightgbm as lgb
from hyperopt import tpe
from hyperopt import Trials

from hyperopt import fmin

from sklearn.metrics import average_precision_score
from hyperopt.pyll.stochastic import sample
from sklearn.metrics import roc_auc_score, f1_score, precision_recall_curve, auc


def focal_isoform_binary_object(pred, dtrain, alpha=0.5, beta=0.0, gamma=2.0):
    #     alpha controls weight of positives
    #     (0,1) less
    #     >1 more or(0-0.5 less, 0.5-1 more)
    #     beta controls the shift of loss function
    #     >0 to left(less weight to well-trained samples)
    #     gamma controls the steepness of loss function
    #     >0
    label = dtrain.get_label()
    x = beta + (2.0 * label - 1) * gamma * pred
    p = 1. / (1. + np.exp(-x))
    # grad = (1 + (alpha - 1) * label) * (2 * label - 1) * (p - 1)
    grad = (1 - label + (label * 2 - 1) * alpha) * (2 * label - 1) * (p - 1)
    # hess = (1 + (alpha - 1) * label) * gamma * (1 - p) * p
    hess = (1 - label + (label * 2 - 1) * alpha) * gamma * (1 - p) * p
    return grad, hess


def lgb_auprc_score_sklearn(y_hat, data):
    y_true = data.get_label()
    precision, recall, _ = precision_recall_curve(y_true, y_hat)
    return 'auprc', auc(recall, precision), True


def lgb_auprc_score(y_hat, data):
    y_true = data.get_label()
    # TODO try not to round yhat
    # y_hat = np.round(y_hat)  # scikits f1 doesn't like probabilities
    return 'auprc', average_precision_score(y_true, y_hat), True


params = {
    'boosting_type': 'gbdt',
    # 'ignore_column': list(range(500)),
    # 'boosting_type': 'dart',
    # 'drop_rate': 0.3,
    # 'max_drop': 50,
    # 'skip_drop': 0.5,
    # 'drop_seed': 6,
    # 'pos_bagging_fraction': 0.001,
    # 'neg_bagging_fraction': 0.001,
    # 'bagging_freq': 10000,
    # 'bagging_seed': 6,
    'objective': 'binary',
    # 'objective': focal_binary_object,
    # 'metric': ['binary_error', 'binary_logloss', "auc"],
    'metric': ["auc"],
    # 'first_metric_only': True,
    # 'is_unbalance': True,
    # "scale_pos_weight": 100,
    'metric_freq': 10,
    'num_leaves': 31,
    'min_data_in_leaf': 20,
    # 'max_bin': 255,
    'num_threads': 32,
    'learning_rate': 0.1,
    'feature_fraction': 1,
    'boost_from_average': False,
    'verbose': 1
}

evals_result = {}
gbm = lgb.train(params=params,
                train_set=dtrain_subset,
                num_boost_round=1000,
                fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=0.5, beta=-1.5, gamma=1.01),
                #                 fobj=lambda x,y:logistic_obj(x,y,imbalance_alpha=1.0),
                valid_sets=[dtrain_subset,dtrain_subset_2],
                valid_names=['test'],
                feval=lgb_auprc_score,
                # early_stopping_rounds=1,
                evals_result=evals_result,
                keep_training_booster=False,
                learning_rates=lambda x: 0.2 * (0.98 ** x),
                callbacks=[early_stopping(1, first_metric_only=False, verbose=True)]
                # init_model=gbm
                # init_model=init_gbm
                )

N_FOLDS = 6

# dtrain = lgb.Dataset('train/GABPA/binary_files/lightGBM.all.MCF-7.chrom_set_test.bin').construct()
dtrain = lgb.Dataset('train/JUND/binary_files/lightGBM.all.HepG2.chrom_set_test.bin').construct()
subset_index = np.random.choice(np.arange(dtrain.num_data()), int(dtrain.num_data() / 10), replace=False)
dtrain_subset = dtrain.subset(subset_index).construct()

# subset_index = np.random.choice(np.arange(dtrain.num_data()), int(dtrain.num_data() / 10), replace=False)
# dtrain_subset_2 = dtrain.subset(subset_index).construct()

subset_index = np.random.choice(np.arange(dtrain_subset.num_data()), int(dtrain_subset.num_data() / 10), replace=False)
dtrain_subset_2 = dtrain_subset.subset(subset_index).construct()

from scripts.early_stopping_avg import early_stopping



cv_results = lgb.cv(params, dtrain_subset, num_boost_round=10000,
                    nfold=6,
                    fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=0.5, beta=1, gamma=1.01),
                    feval=lgb_auprc_score,
                    # early_stopping_rounds=20,
                    seed=6, verbose_eval=1,
                    callbacks=[early_stopping(1, first_metric_only=False, verbose=True)])

start = timer()

# cv_results = lgb.cv(params, dtrain_subset, num_boost_round=100,
#                     nfold=2,
#                     fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=0.5, beta=1, gamma=1.01),
#                     feval=lgb_auprc_score,
#                     # early_stopping_rounds=20,
#                     seed=6, verbose_eval=1,
#                     # callbacks=[early_stopping(1, first_metric_only=False, verbose=True)]
#                     )

evals_result = {}
gbm = lgb.train(params=params,
                # train_set=dtrain_subset,
                train_set=dtrain,
                num_boost_round=100,
                fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=0.5, beta=-1.5, gamma=1.01),
                #                 fobj=lambda x,y:logistic_obj(x,y,imbalance_alpha=1.0),
                # valid_sets=[dtrain_subset, dtrain_subset_2],
                # valid_names=['train', 'test'],
                # valid_sets=[dtrain_subset],
                valid_sets=[dtrain],
                valid_names=['train'],
                feval=lgb_auprc_score,
                early_stopping_rounds=1,
                evals_result=evals_result,
                keep_training_booster=False,
                learning_rates=lambda x: 0.2 * (0.98 ** x),
                # callbacks=[lgb.early_stopping(1, first_metric_only=False, verbose=True)]
                # init_model=gbm
                # init_model=init_gbm
                )

run_time = timer() - start
print(run_time)




112.17675656080246
57.596633109264076
37.916644868440926
16.787454077973962


# def hyperopt_objective(argsDict, dtrain=dtrain, n_folds=6):
def hyperopt_objective(argsDict, dtrain=dtrain_subset, n_folds=3):
    """Objective function for Gradient Boosting Machine Hyperparameter Optimization"""
    # Keep track of evals
    global ITERATION
    ITERATION += 1
    # Make sure parameters that need to be integers are integers
    for parameter_name in ['num_leaves',
                           'min_data_in_leaf',
                           # 'max_depth'
                           ]:
        argsDict[parameter_name] = int(argsDict[parameter_name])
    start = timer()
    params = {
        'boosting_type': 'gbdt',
        # 'ignore_column': list(range(500)),
        # 'boosting_type': 'dart',
        # 'drop_rate': 0.3,
        # 'max_drop': 50,
        # 'skip_drop': 0.5,
        # 'drop_seed': 6,
        # 'pos_bagging_fraction': 0.001,
        # 'neg_bagging_fraction': 0.001,
        # 'bagging_freq': 10000,
        # 'bagging_seed': 6,
        'objective': 'binary',
        # 'objective': focal_binary_object,
        # 'metric': ['binary_error', 'binary_logloss', "auc"],
        # 'metric': ["auc"],
        # 'first_metric_only': True,
        # 'is_unbalance': True,
        # "scale_pos_weight": 100,
        # 'feature_fraction_bynode': True,
        'metric_freq': 10,
        'num_leaves': argsDict['num_leaves'],
        'min_data_in_leaf': argsDict['min_data_in_leaf'],
        # 'min_data_in_leaf': 20,
        # 'max_depth': argsDict['max_depth'],
        # 'min_sum_hessian_in_leaf': argsDict['min_sum_hessian_in_leaf'],
        # 'bagging_fraction': argsDict['bagging_fraction'],
        # 'feature_fraction': argsDict['feature_fraction'],
        # 'lambda_l1': argsDict['lambda_l1'],
        # 'lambda_l2': argsDict['lambda_l2'],
        # 'max_bin': 255,
        'num_threads': 32,
        # 'learning_rate': argsDict['learning_rate'],
        'learning_rate': 0.1,
        'bagging_freq': 1,
        'boost_from_average': False,
        'verbose': -1
    }
    # Perform n_folds cross validation
    # cv_results = lgb.cv(params, dtrain, num_boost_round=10000, nfold=n_folds,
    #                     early_stopping_rounds=20, metrics='auc', seed=50)
    cv_results = lgb.cv(params, dtrain, num_boost_round=300, nfold=n_folds,
                        fobj=lambda x, y: focal_isoform_binary_object(x, y,
                                                                      # alpha=float(
                                                                      #     np.clip(argsDict['alpha'], 0.001, 0.999)),
                                                                      alpha=1. / (1. + np.exp(
                                                                          -argsDict['alpha_isoform'])),
                                                                      beta=argsDict['beta'],
                                                                      gamma=argsDict['gamma']),
                        feval=lgb_auprc_score,
                        early_stopping_rounds=20, seed=6,
                        # verbose_eval=10
                        )
    run_time = timer() - start
    # Extract the best score
    best_score = np.max(cv_results['auprc-mean'])
    # Loss must be minimized
    loss = 1 - best_score
    # Boosting rounds that returned the highest cv score
    n_estimators = int(np.argmax(cv_results['auprc-mean']) + 1)
    print('auprc:{} ITERATION:{} n_estimators:{} run_time:{}'.format(best_score, ITERATION, n_estimators, run_time),
          end="\n")
    # Dictionary with information for evaluation
    return {'loss': loss,
            'params': argsDict,
            'iteration': ITERATION,
            'estimators': n_estimators,
            'train_time': run_time,
            'status': STATUS_OK}
    # return loss


# Define the search space
space = {
    # 'class_weight': hp.choice('class_weight', [None, 'balanced']),
    'num_leaves': hp.qloguniform('num_leaves', np.log(15), np.log(1023), 5),
    # 'max_depth': hp.quniform('max_depth', 3, 63, 1),
    # 'min_sum_hessian_in_leaf': hp.loguniform('min_sum_hessian_in_leaf', np.log(0.001), np.log(1)),
    # 'learning_rate': hp.loguniform('learning_rate', np.log(0.001), np.log(0.2)),
    'min_data_in_leaf': hp.quniform('min_data_in_leaf', 20, 1000, 5),
    # 'lambda_l1': hp.uniform('lambda_l1', 0.0, 1.0),
    # 'lambda_l2': hp.uniform('lambda_l2', 0.0, 1.0),
    # 'bagging_fraction': hp.uniform('bagging_fraction', 0.4, 1.0),
    # 'feature_fraction': hp.uniform('feature_fraction', 0.4, 1.0),
    # 'alpha': hp.loguniform('alpha', np.log(1), np.log(100)),
    # 'alpha': hp.normal('alpha', 0.5, 0.15),
    'alpha_isoform': hp.normal('alpha_isoform', 0, 3),
    'beta': hp.uniform('beta', -10, 10),
    'gamma': hp.loguniform('gamma', np.log(1), np.log(20)),
}

# sample(space)

# Keep track of results
from hyperopt.fmin import generate_trials_to_calculate
bayes_trials = generate_trials_to_calculate(
            [{'num_leaves': 63, 'min_data_in_leaf': 20, 'alpha_isoform': 0, 'beta': -1.5, 'gamma': 1.01}])
# bayes_trials = Trials()

# Global variable
global ITERATION

ITERATION = 0

# Run optimization
best = fmin(fn=hyperopt_objective, space=space, algo=tpe.suggest,
            max_evals=100, trials=bayes_trials, rstate=np.random.RandomState(6))

# Sort the trials with lowest loss (highest AUC) first

bayes_trials_results = sorted(bayes_trials.results, key=lambda x: x['loss'])
bayes_trials_results[:10]
bayes_trials_small=bayes_trials

subset_index = np.random.choice(np.arange(dtrain.num_data()), int(dtrain.num_data() / 100), replace=False)
dtrain_subset = dtrain.subset(subset_index).construct()

params = {
    'boosting_type': 'gbdt',
    # 'ignore_column': list(range(500)),
    # 'boosting_type': 'dart',
    # 'drop_rate': 0.3,
    # 'max_drop': 50,
    # 'skip_drop': 0.5,
    # 'drop_seed': 6,
    # 'pos_bagging_fraction': 0.001,
    # 'neg_bagging_fraction': 0.001,
    # 'bagging_freq': 10000,
    # 'bagging_seed': 6,
    'objective': 'binary',
    # 'objective': focal_binary_object,
    # 'metric': ['binary_error', 'binary_logloss', "auc"],
    # 'metric': ["auc"],
    # 'first_metric_only': True,
    # 'is_unbalance': True,
    # "scale_pos_weight": 100,
    'metric_freq': 10,
    'num_leaves': 31,
    'min_data_in_leaf': 20,
    # 'max_bin': 255,
    'num_threads': 32,
    'learning_rate': 0.1,
    'feature_fraction': 1,
    'boost_from_average': False,
    'verbose': 1
}

cv_results = lgb.cv(params, dtrain_subset, num_boost_round=10000, nfold=3,
                    fobj=lambda x, y: focal_isoform_binary_object(x, y,
                                                                  alpha=1. / (1. + np.exp(
                                                                      -argsDict['alpha_isoform'])),
                                                                  beta=0, gamma=1),
                    feval=lgb_auprc_score,
                    early_stopping_rounds=20, seed=6, verbose_eval=10)
max(cv_results["auprc-mean"])

argsDict = best

params = {
    'boosting_type': 'gbdt',
    # 'ignore_column': list(range(500)),
    # 'boosting_type': 'dart',
    # 'drop_rate': 0.3,
    # 'max_drop': 50,
    # 'skip_drop': 0.5,
    # 'drop_seed': 6,
    # 'pos_bagging_fraction': 0.001,
    # 'neg_bagging_fraction': 0.001,
    # 'bagging_freq': 10000,
    # 'bagging_seed': 6,
    'objective': 'binary',
    # 'objective': focal_binary_object,
    # 'metric': ['binary_error', 'binary_logloss', "auc"],
    # 'metric': ["auc"],
    # 'first_metric_only': True,
    # 'is_unbalance': True,
    # "scale_pos_weight": 100,
    # 'feature_fraction_bynode': False,
    'metric_freq': 10,
    'num_leaves': int(argsDict['num_leaves']),
    'min_data_in_leaf': int(argsDict['min_data_in_leaf']),
    # 'max_depth': argsDict['max_depth'],
    # 'min_sum_hessian_in_leaf': 0.001,
    # 'min_sum_hessian_in_leaf': argsDict['min_sum_hessian_in_leaf'],
    # 'bagging_fraction': argsDict['bagging_fraction'],
    # 'feature_fraction': argsDict['feature_fraction'],
    # 'lambda_l1': argsDict['lambda_l1'],
    # 'lambda_l2': argsDict['lambda_l2'],
    # 'max_bin': 255,
    'num_threads': 32,
    # 'learning_rate': argsDict['learning_rate'],
    'learning_rate': 0.1,
    'bagging_freq': 1,
    'boost_from_average': False,
    'verbose': 1
}

cv_results = lgb.cv(params, dtrain_subset, num_boost_round=10000, nfold=6,
                    fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=argsDict['alpha'], beta=argsDict['beta'],
                                                                  gamma=argsDict['gamma']),
                    feval=lgb_auprc_score,
                    early_stopping_rounds=20, seed=6, verbose_eval=10)

max(cv_results["auprc-mean"])

evals_result = {}
gbm = lgb.train(params=params,
                train_set=dtrain,
                num_boost_round=1000,
                fobj=lambda x, y: focal_isoform_binary_object(x, y, alpha=argsDict['alpha'], beta=argsDict['beta'],
                                                              gamma=argsDict['gamma']),
                #                 fobj=lambda x,y:logistic_obj(x,y,imbalance_alpha=1.0),
                feval=lgb_auprc_score,
                early_stopping_rounds=20,
                evals_result=evals_result,
                keep_training_booster=False,
                learning_rates=lambda x: max(0.2 * (0.98 ** x), 0.01),
                # init_model=gbm
                # init_model=init_gbm
                )

a = [
    {'loss': 0.538821527627692,
     'params': {'alpha': 84.25023221444668, 'beta': 3.82429677823259, 'gamma': 4.766755613108209,
                'min_sum_hessian_in_leaf': 0.01577482235883418, 'num_leaves': 23}, 'iteration': 83, 'estimators': 163,
     'train_time': 19.835279636085033, 'status': 'ok'},
    {'loss': 0.5432677259278311,
     'params': {'alpha': 52.6345975934508, 'beta': 2.3866439100496337, 'gamma': 12.773093162589662,
                'min_sum_hessian_in_leaf': 0.002075441584863259, 'num_leaves': 44}, 'iteration': 86, 'estimators': 144,
     'train_time': 22.59964281693101, 'status': 'ok'},
    {'loss': 0.5435050408889868,
     'params': {'alpha': 69.21035509824684, 'beta': 3.7776317895399765, 'gamma': 2.9957990208563214,
                'min_sum_hessian_in_leaf': 0.0710012500578193, 'num_leaves': 35}, 'iteration': 33, 'estimators': 115,
     'train_time': 15.241986678913236, 'status': 'ok'},
    {'loss': 0.5437308590626344,
     'params': {'alpha': 53.552913631842344, 'beta': 4.293668047714721, 'gamma': 1.0731705381626415,
                'min_sum_hessian_in_leaf': 0.03681059881682353, 'num_leaves': 32}, 'iteration': 70, 'estimators': 96,
     'train_time': 13.097969220019877, 'status': 'ok'},
    {'loss': 0.5439951072227043,
     'params': {'alpha': 24.70686311576184, 'beta': -0.17229772336135296, 'gamma': 2.0533809716399585,
                'min_sum_hessian_in_leaf': 0.11494351628474447, 'num_leaves': 87}, 'iteration': 51, 'estimators': 151,
     'train_time': 27.601539500057697, 'status': 'ok'},
    {'loss': 0.5469903253986683,
     'params': {'alpha': 96.97876566546475, 'beta': 2.1078888940414195, 'gamma': 9.154680014826484,
                'min_sum_hessian_in_leaf': 0.006052604450051556, 'num_leaves': 49}, 'iteration': 91, 'estimators': 215,
     'train_time': 28.905053738504648, 'status': 'ok'},
    {'loss': 0.5472778295324561,
     'params': {'alpha': 32.319763973183505, 'beta': 2.0931674022191302, 'gamma': 1.2093211495250065,
                'min_sum_hessian_in_leaf': 0.3549295928774973, 'num_leaves': 131}, 'iteration': 18, 'estimators': 93,
     'train_time': 23.846948521211743, 'status': 'ok'},
    {'loss': 0.5473851773961442,
     'params': {'alpha': 82.84666453971768, 'beta': 4.774933368779824, 'gamma': 15.444269469197687,
                'min_sum_hessian_in_leaf': 0.008687681141799148, 'num_leaves': 64}, 'iteration': 100, 'estimators': 123,
     'train_time': 21.544361477717757, 'status': 'ok'},
    {'loss': 0.5484552131429545,
     'params': {'alpha': 43.75631201693451, 'beta': 3.390351610792944, 'gamma': 8.055598985363558,
                'min_sum_hessian_in_leaf': 1.307349821660981, 'num_leaves': 343}, 'iteration': 47, 'estimators': 96,
     'train_time': 27.99754294939339, 'status': 'ok'},
    {'loss': 0.54889556409551,
     'params': {'alpha': 20.952739599883802, 'beta': -0.6283707567139105, 'gamma': 15.368059424933069,
                'min_sum_hessian_in_leaf': 0.004890346719512667, 'num_leaves': 352}, 'iteration': 8, 'estimators': 186,
     'train_time': 101.53110141307116, 'status': 'ok'}]





