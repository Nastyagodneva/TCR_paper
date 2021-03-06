###################################################################################################
# File: RandomizedSearchCV-SHAP_LGBM.py
# Version: 0.0
# Date: 12.8.2018
# Noam Bar, noam.bar@weizmann.ac.il
#
# 
# Python version: 2.7
###################################################################################################

from __future__ import print_function
from sklearn.model_selection import RandomizedSearchCV
from sklearn import metrics
# from Analyses.AnalysisHelperFunctions import make_dir_if_not_exists
import lightgbm as lgb
import shap
import os
import pickle
from scipy.stats.stats import spearmanr, pearsonr
from sklearn.metrics import r2_score, precision_recall_curve
from sklearn.model_selection import GroupKFold
import argparse
import numpy as np
import pandas as pd
from addloglevels import sethandlers
from queue.qp import qp, fakeqp
from datetime import datetime
import Utils
import sys
from os import listdir, mkdir, makedirs
from os.path import isfile, join, isdir, exists


# lightgbm params
# LEARNING_RATE = [0.1, 0.05, 0.02, 0.015, 0.01, 0.0075, 0.005, 0.002, 0.001, 0.0005] #Noam's
LEARNING_RATE = [0.1, 0.05, 0.01, 0.005, 0.001, 0.0005] #Mine
# NUM_LEAVES = range(2, 15) #Noam's
NUM_LEAVES = range(2, 15,2) #mine
# MAX_DEPTH = [-1, 2, 3, 4, 5, 10, 20, 40, 50] #Noam's
MAX_DEPTH = [-1, 2, 3, 4, 5, 10]
# MIN_DATA_IN_LEAF = range(1, 44, 2) #Noam's
MIN_DATA_IN_LEAF = range(1, 41, 5)
FEATURE_FRACTION = [i/10. for i in range(4,11)]
METRIC = ['l2']
EARLY_STOPPING_ROUNDS = [None]
N_THREADS = [1]
VERBOSE = [-1]
SILENT = [True]

# N_ESTIMATORS = range(50, 750, 50) #Noam's
N_ESTIMATORS = range(50, 850, 100)
BAGGING_FRACTION = [i/10. for i in range(6,11)]
BAGGING_FREQ = [0, 1, 2, 3, 5, 8, 10]

PREDICTOR_PARAMS = {'learning_rate':LEARNING_RATE, 'max_depth':MAX_DEPTH, 
                    'feature_fraction':FEATURE_FRACTION, 'num_leaves':NUM_LEAVES,
                      'min_data_in_leaf':MIN_DATA_IN_LEAF, 'metric':METRIC,
                      'early_stopping_rounds':EARLY_STOPPING_ROUNDS, 'n_estimators':N_ESTIMATORS,
                      'bagging_fraction':BAGGING_FRACTION, 'bagging_freq':BAGGING_FREQ, 
                      'num_threads':N_THREADS, 'verbose':VERBOSE, 'silent':SILENT}


def randomized_search_cv_and_shap(command_args, Y, idx):
    print ('randomized_search_cv_and_shap', str(idx))
    X = Utils.Load(command_args.path_to_X)
    
    shap_values_dic = {}
    results_df = pd.DataFrame(index=Y.columns, columns=['Size', 'Coefficient_of_determination', 
                                                         'pearson_r', 'pearson_p', 'spearman_r', 
                                                         'spearman_p'])
    predictions_df = pd.DataFrame(index=Y.index, columns=Y.columns)
    for y_name in Y.columns:
        
        y = Y[y_name].dropna().astype(float)
        if y.unique().shape[0] == 2:
            classification_problem = True
        else:
            classification_problem = False
            
        print (y_name, y.shape)
        print (shap.__path__)
        X_temp = X.loc[y.index].copy()
        groups = np.array(range(X_temp.shape[0]))
        group_kfold = GroupKFold(n_splits=command_args.k_folds)
        shap_values = pd.DataFrame(np.nan, index=X_temp.index, columns=X_temp.columns)
        final_pred = pd.DataFrame(index=X_temp.index, columns=[y_name])
        
        try:
            for train_index, test_index in group_kfold.split(X_temp, y, groups):
                X_train, X_test = X_temp.iloc[train_index,:], X_temp.iloc[test_index,:]
                y_train, y_test = y.iloc[train_index], y.iloc[test_index]
                
                if classification_problem:
                    gbm = lgb.LGBMClassifier()
                else:
                    gbm = lgb.LGBMRegressor()
                rscv = RandomizedSearchCV(gbm, PREDICTOR_PARAMS, n_iter=command_args.n_random)
                rscv.fit(X_train, y_train)
                # use best predictor according to random hyper-parameter search
                if classification_problem:
                    y_pred = rscv.best_estimator_.predict_proba(X_test)
                    final_pred.loc[X_test.index, :] = np.expand_dims(y_pred[:,1], 1)
                else:
                    y_pred = rscv.best_estimator_.predict(X_test)
                    final_pred.loc[X_test.index, :] = np.expand_dims(y_pred, 1)
                # run SHAP
                explainer = shap.TreeExplainer(rscv.best_estimator_)
                try:
                    shap_values.loc[X_test.index, :] = explainer.shap_values(X_test) # changed on 3.10.2018, last column is the bias column
                except:
                    shap_values.loc[X_test.index, :] = explainer.shap_values(X_test)[:, :-1]
            shap_values_dic[y_name] = shap_values
            results_df = _evaluate_performance(y_name, final_pred.values.ravel(), y, results_df, classification_problem)
            predictions_df.loc[final_pred.index, y_name] = final_pred.values.ravel()
        except:
            print ('RandomizedSearchCV failed with metabolite %s'%y_name)
            continue
    _save_temporary_files(command_args, idx, shap_values_dic, results_df, predictions_df)
    return

def _save_temporary_files(command_args, idx, shap_values_dic, results_df, predictions_df):
    with open(command_args.output_dir + '/temp_resdf_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(results_df, fout)
    with open(command_args.output_dir + '/temp_shap_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(shap_values_dic, fout)
    with open(command_args.output_dir + '/temp_pred_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(predictions_df, fout)
    return

def _evaluate_performance(y_name, y_pred, y_test, results_df, classification_problem):
    results_df.loc[y_name, 'Size'] = y_pred.shape[0]
    if classification_problem:
        # Prevalence
        results_df.loc[y_name, 'prevalence'] = float(y_test.sum()) / y_test.shape[0]
        # AUC
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred, pos_label=1)
        results_df.loc[y_name, 'AUC'] = metrics.auc(fpr, tpr)
        # PR
        precision, recall, _ = precision_recall_curve(y_test, y_pred)
        results_df.loc[y_name, 'Precision_Recall'] = metrics.auc(recall, precision)
    else:
        results_df.loc[y_name, 'Coefficient_of_determination'] = r2_score(y_true = y_test, y_pred = y_pred)
        results_df.loc[y_name, 'pearson_r'], results_df.loc[y_name, 'pearson_p'] = pearsonr(y_pred, y_test)
        results_df.loc[y_name, 'spearman_r'], results_df.loc[y_name, 'spearman_p'] = spearmanr(y_pred, y_test)
    return results_df


def concat_outputs(command_args):
    print ('concat_outputs')
    all_temp_files = os.listdir(command_args.output_dir)
    resdf_files = [command_args.output_dir + f for f in all_temp_files if f.startswith('temp_resdf_')]
    shap_files = [command_args.output_dir + f for f in all_temp_files if f.startswith('temp_shap_')]
    pred_files = [command_args.output_dir + f for f in all_temp_files if f.startswith('temp_pred_')]
    _concat_files(resdf_files, command_args.output_dir + '/results_df.pkl', how='dataframe')
    _concat_files(shap_files, command_args.output_dir + '/shap_values.pkl', how='dic')
    _concat_files(pred_files, command_args.output_dir + '/predictions_df.pkl', how='dataframe', axis=1)
    return
    
    
def _concat_files(files, final_path, how='dataframe', axis=0):
    if how == 'dataframe':
        final_file = pd.DataFrame()
        for f in files:
            final_file = pd.concat((final_file, pd.read_pickle(f)), axis=axis)
            os.remove(f)
        with open(final_path, 'wb') as fout:
            pickle.dump(final_file, fout)
    elif how == 'dic':
        final_file = {}
        for f in files:
            final_file.update(pd.read_pickle(f))
            os.remove(f)
        with open(final_path, 'wb') as fout:
            pickle.dump(final_file, fout)
    return
    # TODO: add also csv option...
     



def upload_these_jobs(q, command_args):
    print ('upload_these_jobs')
    waiton = []
    Y = Utils.Load(command_args.path_to_Y)
    
    for idx in range(0, Y.shape[1], command_args.n_cols_per_job):
        waiton.append(q.method(randomized_search_cv_and_shap, (command_args, 
                                                               Y.iloc[:, idx:idx+command_args.n_cols_per_job], 
                                                               idx)))
    res = q.waitforresults(waiton)
    # merge the temp results files
    concat_outputs(command_args)
    return res




def main():
    print ('main')
    parser = argparse.ArgumentParser()
    parser.add_argument('output_dir', help='Path to output directory', type=str, default=None)
    parser.add_argument('-n_cols_per_job', help='Number of columns per job', type=int, default=10)
    parser.add_argument('-n_random', help='Number of random samples', type=int, default=20)
    parser.add_argument('-path_to_X', '--path_to_X', help='Path to features data - X', type=str, default='/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/SHAP/dataframes/mar17_features+diet+MPA_species.dat')
    parser.add_argument('-path_to_Y', help='Path to labels - Y', type=str, default='/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/Metabolon/SHAP/dataframes/mar17_metabolomics_unnormed.dat')
    parser.add_argument('-k_folds', help='Number of folds', type=int, default=10)
    parser.add_argument('-only_concat', help='Whether to only concatenate the output files', type=bool, default=False)
    parser.add_argument('-mem_def', help='Number of folds', type=int, default=2)
    command_args = parser.parse_args()
    
    if (not os.path.exists(command_args.path_to_X)) or (not os.path.exists(command_args.path_to_Y)):
        print ("X or Y doesn't exist!"); return
    if command_args.n_cols_per_job < 1 or command_args.n_cols_per_job > 1000:
        print ("n_cols_per_job must be between 1 and 1000"); return

    if command_args.only_concat:
        concat_outputs(command_args); return
        
#     make_dir_if_not_exists(command_args.output_dir)
    if not isdir(command_args.output_dir):
        makedirs(command_args.output_dir)
    
    with open(command_args.output_dir + '/args' + str(datetime.now()), 'w') as handle:
        for arg in vars(command_args):
            handle.write(str(arg) + '\t' + str(getattr(command_args, arg)) + '\n')

    with qp(jobname = 'SHAP-RSCV', q=['himem7.q'], mem_def = str(command_args.mem_def) + 'G', 
            trds_def = 2, tryrerun=True, max_u = 350, delay_batch=15) as q:
#         os.chdir("/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/temp_q_dir/")
        tempDir="/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/temp_q_dir"
        if not isdir(tempDir):
            makedirs(tempDir)
        os.chdir(tempDir)
        q.startpermanentrun()
        upload_these_jobs(q, command_args)
        
        
if __name__ == "__main__":
    sethandlers()
    main()