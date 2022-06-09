######################################################################################################
# File: PredictionModuleShani.py
# Version: 0.0
# Date: 08.11.18
# Updated: 28.11.18
# Shani Ben-Ari Fuchs, shani.ben-ari@weizmann.ac.il
# adopted from Noam Bar
#
# 
# Python version: 2.7
###################################################################################################
#
# IMPORTANT NOTES!
# 1. USE -pathToModelParams parameter to specific the path to excel file containing model parameters! see an example in:
# '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/predictions2/Model_params_files/XGB_model_params_1.xlsx'
#
#
# add:
# 1. scoring for selecting hyper-parameter - enable different metric for regression and classification 
# (currently the user-defined parameter apply only to binary classfication
# 2. multiclass - improve performance evaluation
# 3. bootstrapping ?)
# 4. plots (?)
# 
#
#
from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from sklearn.model_selection import RandomizedSearchCV, GridSearchCV
from sklearn import metrics
# from Analyses.AnalysisHelperFunctions import make_dir_if_not_exists
import lightgbm as lgb
from xgboost import XGBClassifier, XGBRegressor
import shap
import os
import pickle
from scipy.stats.stats import spearmanr, pearsonr
from sklearn.metrics import r2_score, precision_recall_curve, f1_score
from sklearn.model_selection import GroupKFold,StratifiedKFold
import argparse
import numpy as np
import pandas as pd
from addloglevels import sethandlers
from SegalQueue.qp import qp, fakeqp
from datetime import datetime
import Utils
import sys
from os import listdir, mkdir, makedirs
from os.path import isfile, join, isdir, exists
import json
from sklearn.linear_model import LogisticRegression, LinearRegression, Lasso
from sklearn.preprocessing import StandardScaler, scale
from sklearn.feature_selection import f_regression, f_classif, SelectKBest, mutual_info_regression, mutual_info_classif,chi2


def model_and_shap(command_args, Y, idx):
    print ('model_and_shap', str(idx))
    try:
        X = Utils.Load(command_args.path_to_X)
    except:
        X = pd.read_excel(command_args.path_to_X).set_index('BD')
    
    shap_values_dic = {}
    coef_values_dic = {}
    model_params_dic = {}
    results_df = pd.DataFrame(index=Y.columns, columns=['Size', 'Coefficient_of_determination',
                                                         'pearson_r', 'pearson_p', 'spearman_r',
                                                         'spearman_p'])
    predictions_df = pd.DataFrame(index=Y.index, columns=Y.columns)
    
    
    #get predictor_params:
    
    if command_args.pathToModelParams is not None:
        params_df=pd.read_excel(command_args.pathToModelParams)
        predictor_params={}
        for n in range(params_df.shape[0]):
            param=str(params_df.iloc[n,0])
            values=params_df.iloc[n,1:].dropna().tolist()
#             print ('values: ', values)
            values=[int(x) if ('float' in str(type(x)) and ((round(x)==x) or (round(x)==x+1))) else x for x in values]
            predictor_params[param]=values
    else:
        predictor_params=command_args.predictor_params   
      
    #organize Y:
    for y_name in Y.columns:
        try:
            Ys = Y[y_name].loc[X.index].dropna().astype(float)
        except:
            print ('didnt convert y to float')
            Ys = Y[y_name].loc[X.index].dropna()
       
        # define probem type (classification/multiclass classification or regression)
        # and scoring metric for hyperparameter tuning. defaults are: 'roc_auc' for binary, 'f1_macro' for multiclass classification,
        # explained variance for regression
        #other enabled metrics: 'kappa', '
        nCat = Ys.unique().shape[0]
        if nCat == 2:
            classification_problem = True
            if command_args.score_binclass is None:
                scoring = 'roc_auc'
            elif command_args.score_binclass == 'kappa':
                kappa=metrics.make_scorer(metrics.cohen_kappa_score)
                scoring=kappa
            else:
                scoring=command_args.score_binclass
            # make sure the values are 0 and 1:
            y = pd.Series(index=Ys.index, data=np.where(Ys == Ys.unique().min(), 0, 1))
        else:
            if command_args.multiclass and nCat < 7:
                classification_problem = True
                if command_args.score_multiclass is None:
                    f1macro = metrics.make_scorer(f1_score, average='macro')
                    scoring = f1macro
                elif command_args.score_multiclass == 'kappa':
                    kappa=metrics.make_scorer(metrics.cohen_kappa_score)
                    scoring=kappa                
                else:
                    scoring=command_args.score_multiclass
                y = Ys
            else:
                classification_problem = False
                if command_args.score_reg is None:
                    expVar = metrics.make_scorer(metrics.explained_variance_score)
                    scoring = expVar
                else:
                    scoring=command_args.score_reg
                y = Ys
            
        X_temp = X.loc[y.index].copy()
        print (y_name, y.shape)
        print ('X shape is: ', X_temp.shape)
        print (shap.__path__)
        #define 'leave-1-out if required:
        if command_args.k_folds<1:
            command_args.k_folds=X_temp.shape[0]
            print ('using leave-1 out design')
        shap_values = pd.DataFrame(np.nan, index=X_temp.index, columns=X_temp.columns)
        coef_values = pd.DataFrame(np.nan, index=X_temp.index, columns=X_temp.columns)
        final_pred = pd.DataFrame(index=X_temp.index, columns=[y_name])
        selected_model_params_df = pd.DataFrame()
        
        if classification_problem:
            if command_args.model == 'XGB' : model = XGBClassifier
            elif command_args.model == 'LGBM' : model = lgb.LGBMClassifier
            else: model = LogisticRegression
            
            if command_args.k_folds>0:
                group_kfold = StratifiedKFold(n_splits=command_args.k_folds)
                groups=None
            
        else:
            groups = np.array(range(X_temp.shape[0]))
            group_kfold = GroupKFold(n_splits=command_args.k_folds)
            if command_args.model == 'XGB' : model = XGBRegressor
            if command_args.model == 'LGBM' : model = lgb.LGBMRegressor
            if command_args.model == 'lin' : model = Lasso
        print ('model to be used: ' ,model) 
        split = 0
#         try:
        for train_index, test_index in group_kfold.split(X_temp, y, groups):
            split = split + 1
            print('split:', split)
            X_train, X_test = X_temp.iloc[train_index, :], X_temp.iloc[test_index, :]
            y_train, y_test = y.iloc[train_index], y.iloc[test_index]
            #compare train-test data:
            if classification_problem:
                print ('y_train nornalized value counts; ')
                print (y_train.value_counts(normalize=True))
                print ('y_test normalized value counts; ')
                print (y_test.value_counts(normalize=True))
            else:
                print ('y_train: mean=%s, std=%s' %(y_train.mean(),y_train.std())) 
                print ('y_test: mean=%s, std=%s' %(y_test.mean(),y_test.std()))          
            
            #apply univariate feature selection if needed:
            if (command_args.nPreSelectedFeatures is not None) and (command_args.nPreSelectedFeatures <X.shape[1]) :
                k=command_args.nPreSelectedFeatures
                print ('choosing %s features' %k)
                
                if command_args.scaleForFeatureSelection:
                    for col in X.columns:
                        X[col] = X[col].fillna(X[col].median())
                        scaler = StandardScaler()
                        scaler.fit(X_train)
                        X_train = pd.DataFrame(index=X_train.index,columns=X_train.columns, data=scaler.transform(X_train))
                        X_test = pd.DataFrame(index=X_test.index,columns=X_test.columns, data=scaler.transform(X_test))
                    print (' data was scaled using standardScaler')
                
                if classification_problem:
                    if command_args.featureSelectionMethod_clas is None:
                        featureSelectionMethod=f_classif
                    else:
                        print ('feature selection method= ',command_args.featureSelectionMethod_clas)
                        if command_args.featureSelectionMethod_clas == 'f_classif':
                            featureSelectionMethod=f_classif
                        elif command_args.featureSelectionMethod_clas == 'mutual_info_classif':
                            featureSelectionMethod=mutual_info_classif
                        elif command_args.featureSelectionMethod_clas == 'chi2':
                            featureSelectionMethod= chi2
                        else:
                            print ('feature selection method unknown')
                else:
                    if command_args.featureSelectionMethod_reg is None:
                        featureSelectionMethod=f_regression
                    else:
                        print ('feature selection method= ',command_args.featureSelectionMethod_reg)
                        if command_args.featureSelectionMethod_reg == 'f_regression':
                            featureSelectionMethod=f_regression
                        elif command_args.featureSelectionMethod_reg == 'mutual_info_regression':
                            featureSelectionMethod=mutual_info_regression
                        else:
                            print ('feature selection method unknown')
   
                selector = SelectKBest(featureSelectionMethod, k)
                selector.fit(X_train, y_train)
                idxs_selected = selector.get_support(indices=True)
                X_train_new = X_train.iloc[:, idxs_selected]
                X_test_new = X_test.iloc[:, idxs_selected]
                X_train=X_train_new
                X_test=X_test_new
            selected_features=X_test.columns.tolist()
          
            
            # for linear models, scale data:
            if command_args.model == 'lin':
                for col in X.columns:
                    X[col] = X[col].fillna(X[col].median())
                scaler = StandardScaler()
                scaler.fit(X_train)
                X_train = pd.DataFrame(index=X_train.index,columns=X_train.columns, data=scaler.transform(X_train))
                X_test = pd.DataFrame(index=X_test.index,columns=X_test.columns, data=scaler.transform(X_test))
#                 print (' data was scaled using standardScaler')
#                 print (' nans were filled with column median value')
                
                
            # build the parameter search object:
#             parser.add_argument('-parameter_search_type', help='randomSearch or gridSearch or None', type=str, default=None)   
#                 fit_dict={"eval_set":[(X_test, y_test)],
#           "early_stopping_rounds":5,
#           "eval_metric":"logloss","verbose":5}
            print ('scoring metric used: ', scoring)
            if command_args.parameter_search_type == 'randomSearch':
                modelObject = RandomizedSearchCV(model(), param_distributions=predictor_params, n_iter=command_args.n_random,
                n_jobs=8, cv=5, scoring=scoring)
            elif command_args.parameter_search_type == 'gridSearch':
                modelObject = GridSearchCV(model(), param_grid=predictor_params,
                                            n_jobs=8, cv=5, scoring=scoring)
            else:
                model_params = predictor_params
                modelObject = model(**model_params)
            
            modelObject.fit(X_train, y_train)
            
            #save model if required:
            if command_args.saveModel:
                #savwe split model:
                modelFileName = '%s_split%s_model.pkl' % (y_name, split)
                f1 = command_args.output_dir + modelFileName
                pickle.dump(modelObject, open(f1, "wb"))  
                # save split scaler:
                scalerFileName = '%s_split%s_scaler.pkl' % (y_name, split)
                f2 = command_args.output_dir + scalerFileName
                try:
                    pickle.dump(scaler, open(f2, "wb"))
                except:
                    print ('no scaler object to save')  
                # save split selected_features
                selectedFeaturesFileName = '%s_split%s_selectedFeatures.pkl' % (y_name, split)
                f3 = command_args.output_dir + selectedFeaturesFileName
                pickle.dump(selected_features, open(f3, "wb"))  # save split model
                       
            # use best predictor according to random hyper-parameter search
            if (command_args.parameter_search_type != 'randomSearch') and (command_args.parameter_search_type != 'gridSearch'):
                if (classification_problem) and not command_args.multiclass:
                    y_pred = modelObject.predict_proba(X_test)
                    final_pred.loc[X_test.index, :] = np.expand_dims(y_pred[:, 1], 1)
                else:
                    y_pred = modelObject.predict(X_test)
                    final_pred.loc[X_test.index, :] = np.expand_dims(y_pred, 1)
                # document best parameters chosen:
                selected_model_params_df = pd.DataFrame()
                # run SHAP or get coeffs
                if not command_args.useShortVersion:
                    if command_args.model == 'lin':
                        # get p-values if linear model is used:
                        try:
                            coef_values.loc[X_test.index, selected_features] = modelObject.coef_ # changed on 3.10.2018, last column is the bias column
                        except:
                            coef_values.loc[X_test.index, selected_features] = modelObject.coef_.mean()
                    else:
                        explainer = shap.TreeExplainer(modelObject)
                        try:
                            shap_values.loc[X_test.index, selected_features] = explainer.shap_values(X_test)  # changed on 3.10.2018, last column is the bias column
                        except:
                            shap_values.loc[X_test.index, selected_features] = explainer.shap_values(X_test)[:, :-1]
            else:
                if (classification_problem) and not command_args.multiclass:
                    y_pred = modelObject.best_estimator_.predict_proba(X_test)
                    final_pred.loc[X_test.index, :] = np.expand_dims(y_pred[:, 1], 1)
                else:
                    y_pred = modelObject.best_estimator_.predict(X_test)
                    final_pred.loc[X_test.index, :] = np.expand_dims(y_pred, 1)
                
                
                if not command_args.useShortVersion:
                    # document best parameters chosen:
                    df = pd.DataFrame(index=[split], data=modelObject.best_params_)
                    selected_model_params_df = pd.concat([selected_model_params_df, df])
                    # run SHAP or get coeffs
                    if command_args.model == 'lin':
                        # get p-values if linear model is used:
                        try:
                            coef_values.loc[X_test.index, selected_features] = modelObject.best_estimator_.coef_  # changed on 3.10.2018, last column is the bias column
                        except:
                            coef_values.loc[X_test.index, selected_features] = modelObject.best_estimator_.coef_.mean() 
                    else:
                        explainer = shap.TreeExplainer(modelObject.best_estimator_)
                        try:
                            shap_values.loc[X_test.index, selected_features] = explainer.shap_values(X_test)  # changed on 3.10.2018, last column is the bias column
                        except:
                            try:
                                shap_values.loc[X_test.index, selected_features] = explainer.shap_values(X_test)[:, :-1]
                            except:
                                print ('couldnt calculate shap values')
                            
                         
                                              
        shap_values_dic[y_name] = shap_values.loc[:,shap_values.sum()>0]
        coef_values_dic[y_name] = coef_values.loc[:,coef_values.sum()>0]
        model_params_dic[y_name] = selected_model_params_df
        results_df = _evaluate_performance(y_name, final_pred.values.ravel(), y, results_df, classification_problem, command_args.binary_threshold, nCat)
        predictions_df.loc[final_pred.index, y_name] = final_pred.values.ravel()
#         except:
#             print ('prediction failed with phenotype %s'%y_name)
#             continue
    _save_temporary_files(command_args, idx, shap_values_dic,coef_values_dic, model_params_dic, results_df, predictions_df)
    return


# def prediction_bs(q2,command_args):
    

def _save_temporary_files(command_args, idx, shap_values_dic, coef_values_dic,model_params_dic, results_df, predictions_df):
    with open(command_args.output_dir + '/temp_resdf_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(results_df, fout)
    with open(command_args.output_dir + '/temp_shap_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(shap_values_dic, fout)
    with open(command_args.output_dir + '/temp_coef_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(coef_values_dic, fout)
    with open(command_args.output_dir + '/temp_model_params_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(model_params_dic, fout)
    with open(command_args.output_dir + '/temp_pred_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(predictions_df, fout)
    return

def _evaluate_performance(y_name, y_pred, y_test, results_df, classification_problem, binary_threshold, nCat):
    results_df.loc[y_name, 'Size'] = y_pred.shape[0]
    if classification_problem:
        print ('nCat: ', nCat)
        if nCat == 2:
            y_pred_binary = [0 if x < binary_threshold else 1 for x in list(y_pred)]
            # Prevalence   
            print ('prevalence: ', float(y_test.sum()) / y_test.shape[0])
            results_df.loc[y_name, 'prevalence'] = float(y_test.sum()) / y_test.shape[0]
            # AUC
            fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred, pos_label=1)
            results_df.loc[y_name, 'AUC'] = metrics.auc(fpr, tpr)
            print ('AUC:', metrics.auc(fpr, tpr))
            # PR
            precision, recall, _ = precision_recall_curve(y_test, y_pred)
            results_df.loc[y_name, 'Precision_Recall'] = metrics.auc(recall, precision)
            print ('PR:', metrics.auc(recall, precision))
            # Cohen's Kappa:
            print ('Binary probability threshold=%s' % binary_threshold)
            results_df.loc[y_name, 'binary_threshold'] = binary_threshold
            kappa = metrics.cohen_kappa_score(y_test, y_pred_binary)
            print ("cohen's kappa: ", kappa)
            results_df.loc[y_name, 'kappa'] = kappa
            # f1 scores:
            f1score = metrics.f1_score(y_test, y_pred_binary)
            f1scoreMacro = metrics.f1_score(y_test, y_pred_binary, average='macro')
            print ('f1score:', f1score)
            print ('f1scoremacro:' , f1scoreMacro)
            results_df.loc[y_name, 'f1score'] = f1score
            results_df.loc[y_name, 'f1scoreMacro'] = f1scoreMacro
        else:
            print ('prevalence:')
            print (y_test.value_counts(normalize=True))
            f1scoreMacro = metrics.f1_score(y_test, y_pred, average='macro')
            print ('f1scoremacro:' , f1scoreMacro)
            results_df.loc[y_name, 'f1scoreMacro'] = f1scoreMacro
            kappa = metrics.cohen_kappa_score(y_test, y_pred)
            print ("cohen's kappa: ", kappa)
            results_df.loc[y_name, 'kappa'] = kappa
            
    else:
        results_df.loc[y_name, 'Coefficient_of_determination'] = r2_score(y_true=y_test, y_pred=y_pred)
        results_df.loc[y_name, 'pearson_r'], results_df.loc[y_name, 'pearson_p'] = pearsonr(y_pred, y_test)
        results_df.loc[y_name, 'spearman_r'], results_df.loc[y_name, 'spearman_p'] = spearmanr(y_pred, y_test)
        results_df.loc[y_name, 'explained_variance_score'] = metrics.explained_variance_score(y_true=y_test, y_pred=y_pred)
        
        print ('pearson r,p:', pearsonr(y_pred, y_test))
        print ('explained_variance_score: ', metrics.explained_variance_score(y_true=y_test, y_pred=y_pred))
    return results_df


def concat_outputs(command_args):
    print ('concat_outputs')
    all_temp_files = os.listdir(command_args.output_dir)
    resdf_files = [command_args.output_dir + f for f in all_temp_files if f.startswith('temp_resdf_')]
    shap_files = [command_args.output_dir + f for f in all_temp_files if f.startswith('temp_shap_')]
    coef_files = [command_args.output_dir + f for f in all_temp_files if f.startswith('temp_coef_')]
    model_params_files = [command_args.output_dir + f for f in all_temp_files if f.startswith('temp_model_params_')]
    pred_files = [command_args.output_dir + f for f in all_temp_files if f.startswith('temp_pred_')]
    _concat_files(resdf_files, command_args.output_dir + '/results_df.pkl', how='dataframe')
    _concat_files(shap_files, command_args.output_dir + '/shap_values.pkl', how='dic')
    _concat_files(coef_files, command_args.output_dir + '/coef_values.pkl', how='dic')
    _concat_files(model_params_files, command_args.output_dir + '/selected_model_params.pkl', how='dic')
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
    try:
        Y = Utils.Load(command_args.path_to_Y)
    except:
        Y = pd.read_excel(command_args.path_to_Y).set_index('BD')  # enable loading dataframe from excel files
    
    # fill missing target values with median target value:
    if command_args.fillMissingY:
        for col in Y.columns:
            try:
                Y[col] = Y[col].fillna(Y[col].median())
            except:
                try:
                    Y[col] = Y[col].fillna(Y[col].mode())
                except:
                    print ('couldnt fill nas')
    
    for idx in range(0, Y.shape[1], command_args.n_cols_per_job):
        print ('start running prediction for idx:', idx)
        waiton.append(q.method(model_and_shap, (command_args, Y.iloc[:, idx:idx + command_args.n_cols_per_job], idx)))
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
    parser.add_argument('-k_folds', help='Number of folds, for leave one out design use 0, default is 10', type=int, default=10)
    parser.add_argument('-only_concat', help='Whether to only concatenate the output files', type=bool, default=False)
    parser.add_argument('-mem_def', help='mem_def', type=int, default=2)
    parser.add_argument('-predictor_params', help='model parameters (insert the dictionary directly. use '' to encapsulate the\
whole dict, and "" for each string in the dicitonary. use only lists as values. NOTE - if the argument pathToModelParams\
is not None, it overrides this parameter', type=json.loads, default={})
    parser.add_argument('-parameter_search_type', help='randomSearch or gridSearch or None. None will use default model parameters', \
    type=str, default=None)
    parser.add_argument('-model', help='which model to use - "LGBM"/"XGB"/"lin"', type=str, default="LGBM")
    parser.add_argument('-score_binclass', help='scoring metric used for binary classification hyperparameter tuning', type=str, default=None)
    parser.add_argument('-score_multiclass', help='scoring metric used for multiclass classification hyperparameter tuning', type=str, default=None)
    parser.add_argument('-score_reg', help='scoring metric used for reg hyperparameter turession', type=str, default=None)
    parser.add_argument('-n_bs', help='number of bootstrapping to conduct', type=int, default=0)
    parser.add_argument('-fillMissingY', help='whether to fill in missing target values with the median value', type=bool, default=True)
    parser.add_argument('-saveModel', help='whether to save the split model', type=bool, default=0)
    parser.add_argument('-binary_threshold', help='float, the probablity value from which onwards the prediction is considered positive for scoring calculations',
    type=float, default=0.5)
    parser.add_argument('-multiclass', help='bolean, execute multiclass classification or not (will default into regression in case of more than two targets',
    type=bool, default=0)
    parser.add_argument('-nPreSelectedFeatures', help='int, number of features to be selected by univariate selection before model\
    application', type=int, default=999999)
    parser.add_argument('-useShortVersion', help='boolean, if True, do not document params and do not calculate shap/coefs', type=bool,default=False)
    parser.add_argument('-pathToModelParams', help='path to excel file containing model params, first column is param\
name and the other columns in the row are the possible values. if this value is not None, then the parameters in this file\
override what is written in the argument predictor_params', type=str,default=None)
    parser.add_argument('-featureSelectionMethod_clas', help='string,use any eligible sklearn\
.feature_selection method', type=str,default=None)
    parser.add_argument('-featureSelectionMethod_reg', help='string,use any eligible sklearn\
.feature_selection method', type=str,default=None)
    parser.add_argument('-scaleForFeatureSelection',help='bolean',type=bool,default=0)
    
    command_args = parser.parse_args()
    model = command_args.model
#     print ("MODEL USED:", model)
    
    # check folder existance: 
    if (not os.path.exists(command_args.path_to_X)) or (not os.path.exists(command_args.path_to_Y)):
        print ("X or Y doesn't exist!"); return
    if command_args.n_cols_per_job < 1 or command_args.n_cols_per_job > 1000:
        print ("n_cols_per_job must be between 1 and 1000"); return
    #     make_dir_if_not_exists(command_args.output_dir)
    if not isdir(command_args.output_dir):
        makedirs(command_args.output_dir)
        print ('made new output dir')
    
    # save parameters to a text file:
    print ('saving parameters to file')
    argFile = command_args.output_dir + 'argFile.txt'
    new_file = open(argFile, mode="w")
    for arg in vars(command_args): 
        new_file.write(str(arg) +": "+ str(getattr(command_args, arg))+"\n") 
    new_file.write('NOTE! if -pathToModelParams is not None, relate to the parameters in this file as model paramters!')
    new_file.close()
    
    if command_args.only_concat:
        concat_outputs(command_args); return
       
    with open(command_args.output_dir + '/args' + str(datetime.now()), 'w') as handle:
        for arg in vars(command_args):
            handle.write(str(arg) + '\t' + str(getattr(command_args, arg)) + '\n')

    with qp(jobname='SHAP-modelObject', q=['himem7.q'], mem_def=str(command_args.mem_def) + 'G',
            trds_def=4, tryrerun=True, max_u=350) as q:
#         os.chdir("/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/temp_q_dir/")
        tempDir = "/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/temp_q_dir"
        if not isdir(tempDir):
            makedirs(tempDir)
        os.chdir(tempDir)
        q.startpermanentrun()
        upload_these_jobs(q, command_args)
    
#     if command_args.n_bs !=0:
#         with qp(jobname='bsPrediction', q=['himem7.q'], mem_def=str(command_args.mem_def) + 'G',
#             trds_def=4, tryrerun=True, max_u=350) as q2:
#             tempDir = "/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/temp_q_dir"
#             if not isdir(tempDir):
#                 makedirs(tempDir)
#             os.chdir(tempDir)
#             q.startpermanentrun()
#             prediction_bs(q2,command_args)
        
        
        
if __name__ == "__main__":
    sethandlers()
    main()
