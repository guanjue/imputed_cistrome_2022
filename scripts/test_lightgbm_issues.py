
import warnings
from operator import gt, lt
import numpy as np
import lightgbm as lgb
from lightgbm.callback import EarlyStopException
from lightgbm.callback import _format_eval_result
from lightgbm.compat import range_



from sklearn.metrics import average_precision_score


def lgb_auprc_score(y_hat, data):
    y_true = data.get_label()
    # TODO try not to round yhat
    # y_hat = np.round(y_hat)  # scikits f1 doesn't like probabilities
    return 'auprc', average_precision_score(y_true, y_hat), True



def early_stopping(stopping_rounds, first_metric_only=False, verbose=True):
    best_score = []
    best_score_avg = []
    best_iter = []
    best_iter_avg = []
    best_score_list = []
    best_score_avg_list = []
    cmp_op = []
    enabled = [True]

    def _init(env):
        enabled[0] = not any((boost_alias in env.params
                              and env.params[boost_alias] == 'dart') for boost_alias in ('boosting',
                                                                                         'boosting_type',
                                                                                         'boost'))
        if not enabled[0]:
            warnings.warn('Early stopping is not available in dart mode')
            return
        if not env.evaluation_result_list:
            raise ValueError('For early stopping, '
                             'at least one dataset and eval metric is required for evaluation')
        if verbose:
            msg = "Training until validation scores don't improve for {} rounds."
            print(msg.format(stopping_rounds))
        for eval_ret in env.evaluation_result_list:
            best_iter.append(0)
            best_score_list.append(None)
            if eval_ret[3]:
                best_score.append(float('-inf'))
                # best_score_avg = float('-inf')
                cmp_op.append(gt)
            else:
                best_score.append(float('inf'))
                # best_score_avg = float('inf')
                cmp_op.append(lt)
        best_score_avg.append(None)
        best_iter_avg.append(None)
        best_score_avg_list.append(None)

    def _callback(env):
        if not cmp_op:
            _init(env)
        if not enabled[0]:
            return
        print(env.evaluation_result_list)
        for i in range_(len(env.evaluation_result_list)):
            score = env.evaluation_result_list[i][2]
            if best_score_list[i] is None or cmp_op[i](score, best_score[i]):
                best_score[i] = score
                best_iter[i] = env.iteration
                best_score_list[i] = env.evaluation_result_list
            elif env.iteration - best_iter[i] >= stopping_rounds:
                if verbose:
                    print('Early stopping, best iteration is:\n[%d]\t%s' % (
                        best_iter[i] + 1, '\t'.join([_format_eval_result(x) for x in best_score_list[i]])))
                raise EarlyStopException(best_iter[i], best_score_list[i])
            if env.iteration == env.end_iteration - 1:
                if verbose:
                    print('Did not meet early stopping. Best iteration is:\n[%d]\t%s' % (
                        best_iter[i] + 1, '\t'.join([_format_eval_result(x) for x in best_score_list[i]])))
                raise EarlyStopException(best_iter[i], best_score_list[i])
            if first_metric_only:  # the only first metric is used for early stopping
                break

    _callback.order = 30
    return _callback



data = np.random.random((500, 10))
y = [1] * 250 + [0] * 250
lgb_train = lgb.Dataset(data, y, free_raw_data=True)


data = np.random.random((500, 10))
y = [1] * 250 + [0] * 250
lgb_test = lgb.Dataset(data, y, free_raw_data=True)


params = {
    'objective': 'binary',
    'verbose': 1,
    # 'metric': ['binary_logloss', 'auc']
    'metric': ['auc', 'binary_logloss']
}
evals_result = {}
gbm = lgb.train(params=params,
                train_set=lgb_train,
                # valid_sets=[lgb_train, lgb_test],
                valid_sets=[lgb_train,lgb_test],
                feval=lgb_auprc_score,
                num_boost_round=1000,
                evals_result=evals_result,
                callbacks=[early_stopping(1, first_metric_only=False, verbose=True)]
                )






def early_stopping(stopping_rounds, first_metric_only=False, verbose=True):
    """Create a callback that activates early stopping.
    Note
    ----
    Activates early stopping.
    The model will train until the validation score stops improving.
    Validation score needs to improve at least every ``early_stopping_rounds`` round(s)
    to continue training.
    Requires at least one validation data and one metric.
    If there's more than one, will check all of them. But the training data is ignored anyway.
    To check only the first metric set ``first_metric_only`` to True.
    Parameters
    ----------
    stopping_rounds : int
       The possible number of rounds without the trend occurrence.
    first_metric_only : bool, optional (default=False)
       Whether to use only the first metric for early stopping.
    verbose : bool, optional (default=True)
        Whether to print message with early stopping information.
    Returns
    -------
    callback : function
        The callback that activates early stopping.
    """
    best_score = []
    best_score_avg = [None]
    best_iter = []
    best_iter_avg = [None]
    best_score_list = []
    best_score_avg_list = [None]
    cmp_op = []
    enabled = [True]
    first_metric = ['']
    def _init(env):
        enabled[0] = not any((boost_alias in env.params
                              and env.params[boost_alias] == 'dart') for boost_alias in ('boosting',
                                                                                         'boosting_type',
                                                                                         'boost'))
        if not enabled[0]:
            warnings.warn('Early stopping is not available in dart mode')
            return
        if not env.evaluation_result_list:
            raise ValueError('For early stopping, '
                             'at least one dataset and eval metric is required for evaluation')
        if verbose:
            msg = "Training until validation scores don't improve for {} rounds"
            print(msg.format(stopping_rounds))
        # split is needed for "<dataset type> <metric>" case (e.g. "train l1")
        first_metric[0] = env.evaluation_result_list[0][1].split(" ")[-1]
        for eval_ret in env.evaluation_result_list:
            best_iter.append(0)
            best_score_list.append(None)
            if eval_ret[3]:
                best_score.append(float('-inf'))
                cmp_op.append(gt)
            else:
                best_score.append(float('inf'))
                cmp_op.append(lt)
    def _final_iteration_check(env, eval_name_splitted, i):
        if env.iteration == env.end_iteration - 1:
            if verbose:
                print('Did not meet early stopping. Best iteration is:\n[%d]\t%s' % (
                    best_iter[i] + 1, '\t'.join([_format_eval_result(x) for x in best_score_list[i]])))
                if first_metric_only:
                    print("Evaluated only: {}".format(eval_name_splitted[-1]))
            raise EarlyStopException(best_iter[i], best_score_list[i])
    def _callback(env):
        if not cmp_op:
            _init(env)
        if not enabled[0]:
            return
        # score = np.mean([env.evaluation_result_list[i][2] for i in range(len(env.evaluation_result_list))])
        # if best_score_avg[0] is None or (score > best_score_avg[0]):
        #     best_score_avg[0] = score
        #     best_iter_avg[0] = env.iteration
        #     best_score_avg_list[0] = env.evaluation_result_list
        #     print(env.iteration, env.evaluation_result_list, score)
        # if env.iteration == env.end_iteration - 1:
        #     if verbose:
        #         print('Did not meet early stopping. Best iteration is:\n[%d]\t%s' % (
        #             best_iter_avg[0] + 1, '\t'.join([_format_eval_result(x) for x in best_score_avg_list[0]])))
        #     raise EarlyStopException(best_iter_avg[0], best_score_avg_list[0])
        # if env.iteration - best_iter_avg[0] >= stopping_rounds:
        #     if verbose:
        #         print('Early stopping, best iteration is:\n[%d]\t%s' % (
        #             best_iter_avg[0] + 1, '\t'.join([_format_eval_result(x) for x in best_score_avg_list[0]])))
        #     raise EarlyStopException(best_iter_avg[0], best_score_avg_list[0])
        print(env.evaluation_result_list)
        for i in range_(len(env.evaluation_result_list)):
            score = env.evaluation_result_list[i][2]
            if best_score_list[i] is None or cmp_op[i](score, best_score[i]):
                best_score[i] = score
                best_iter[i] = env.iteration
                best_score_list[i] = env.evaluation_result_list
            # split is needed for "<dataset type> <metric>" case (e.g. "train l1")
            eval_name_splitted = env.evaluation_result_list[i][1].split(" ")
            if first_metric_only and first_metric[0] != eval_name_splitted[-1]:
                print(env.evaluation_result_list[i],"first_metric_only")
                continue  # use only the first metric for early stopping
            if ((env.evaluation_result_list[i][0] == "cv_agg" and eval_name_splitted[0] == "train"
                 or env.evaluation_result_list[i][0] == env.model._train_data_name)):
                print(env.evaluation_result_list[i], "no training")
                _final_iteration_check(env, eval_name_splitted, i)
                continue  # train data for lgb.cv or sklearn wrapper (underlying lgb.train)
            elif env.iteration - best_iter[i] >= stopping_rounds:
                if verbose:
                    print('Early stopping, best iteration is:\n[%d]\t%s' % (
                        best_iter[i] + 1, '\t'.join([_format_eval_result(x) for x in best_score_list[i]])))
                    if first_metric_only:
                        print("Evaluated only: {}".format(eval_name_splitted[-1]))
                raise EarlyStopException(best_iter[i], best_score_list[i])
            _final_iteration_check(env, eval_name_splitted, i)
    _callback.order = 30
    return _callback
