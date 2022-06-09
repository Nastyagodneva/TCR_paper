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
import yaml

yaml_param_file= '/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/\
predictions2/outcomes/UnplannedPCI_params'
yaml_default_file= '/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/\
predictions2/prediction_module_default_args.yaml'

#get params:
with open(yaml_param_file, 'r') as yaml_file:
    params = yaml.load(yaml_file)
output_dir = params['output_dir']
n_cols_per_job = params['n_cols_per_job']
n_random = params['n_random']
path_to_X = params['path_to_X']
path_to_Y = params['path_to_Y']
k_folds = params['k_folds']
only_concat = params['only_concat']
mem_def = params['mem_def']
predictor_params = params['predictor_params']
parameter_search_type = params['parameter_search_type']
model = params['model']
score_binclass = params['score_binclass']
score_multiclass = params['score_multiclass']
score_reg = params['score_reg']
n_bs = params['n_bs']
fillMissingY = params['fillMissingY']
saveModel = params['saveModel']
binary_threshold = params['binary_threshold']
multiclass = params['multiclass']
nPreSelectedFeatures = params['nPreSelectedFeatures']
useShortVersion = params['useShortVersion']
pathToModelParams = params['pathToModelParams']
featureSelectionMethod_class = params['featureSelectionMethod_class']
featureSelectionMethod_reg = params['featureSelectionMethod_reg']
scaleForFeatureSelection = params['scaleForFeatureSelection']

def model_and_shap(Y, idx):
    print ('model_and_shap', str(idx))
    try:
        X = Utils.Load(path_to_X)
    except:
        X = pd.read_excel(path_to_X).set_index('BD')
    
    shap_values_dic = {}
    coef_values_dic = {}
    model_params_dic = {}
    results_df = pd.DataFrame(index=Y.columns, columns=['Size', 'Coefficient_of_determination',
                                                         'pearson_r', 'pearson_p', 'spearman_r',
                                                         'spearman_p'])
    predictions_df = pd.DataFrame(index=Y.index, columns=Y.columns)
    
    
    #get predictor_params:
    
    if pathToModelParams is not None:
        params_df=pd.read_excel(pathToModelParams)
        predictor_params={}
        for n in range(params_df.shape[0]):
            param=str(params_df.iloc[n,0])
            values=params_df.iloc[n,1:].dropna().tolist()
#             print ('values: ', values)
            values=[int(x) if ('float' in str(type(x)) and ((round(x)==x) or (round(x)==x+1))) else x for x in values]
            predictor_params[param]=values
    # else:
    #     predictor_params=predictor_params
      
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
            if score_binclass is None:
                scoring = 'roc_auc'
            elif score_binclass == 'kappa':
                kappa=metrics.make_scorer(metrics.cohen_kappa_score)
                scoring=kappa
            else:
                scoring=score_binclass
            # make sure the values are 0 and 1:
            y = pd.Series(index=Ys.index, data=np.where(Ys == Ys.unique().min(), 0, 1))
        else:
            if multiclass and nCat < 7:
                classification_problem = True
                if score_multiclass is None:
                    f1macro = metrics.make_scorer(f1_score, average='macro')
                    scoring = f1macro
                elif score_multiclass == 'kappa':
                    kappa=metrics.make_scorer(metrics.cohen_kappa_score)
                    scoring=kappa                
                else:
                    scoring=score_multiclass
                y = Ys
            else:
                classification_problem = False
                if score_reg is None:
                    expVar = metrics.make_scorer(metrics.explained_variance_score)
                    scoring = expVar
                else:
                    scoring=score_reg
                y = Ys
            
        X_temp = X.loc[y.index].copy()
        print (y_name, y.shape)
        print ('X shape is: ', X_temp.shape)
        print (shap.__path__)
        #define 'leave-1-out if required:
        if k_folds<1:
            k_folds=X_temp.shape[0]
            print ('using leave-1 out design')
        shap_values = pd.DataFrame(np.nan, index=X_temp.index, columns=X_temp.columns)
        coef_values = pd.DataFrame(np.nan, index=X_temp.index, columns=X_temp.columns)
        final_pred = pd.DataFrame(index=X_temp.index, columns=[y_name])
        selected_model_params_df = pd.DataFrame()
        
        if classification_problem:
            if model == 'XGB' : model = XGBClassifier
            elif model == 'LGBM' : model = lgb.LGBMClassifier
            else: model = LogisticRegression
            
            if k_folds>0:
                group_kfold = StratifiedKFold(n_splits=k_folds)
                groups=None
            
        else:
            groups = np.array(range(X_temp.shape[0]))
            group_kfold = GroupKFold(n_splits=k_folds)
            if model == 'XGB' : model = XGBRegressor
            if model == 'LGBM' : model = lgb.LGBMRegressor
            if model == 'lin' : model = Lasso
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
            if (nPreSelectedFeatures is not None) and (nPreSelectedFeatures <X.shape[1]) :
                k=nPreSelectedFeatures
                print ('choosing %s features' %k)
                
                if scaleForFeatureSelection:
                    for col in X.columns:
                        X[col] = X[col].fillna(X[col].median())
                        scaler = StandardScaler()
                        scaler.fit(X_train)
                        X_train = pd.DataFrame(index=X_train.index,columns=X_train.columns, data=scaler.transform(X_train))
                        X_test = pd.DataFrame(index=X_test.index,columns=X_test.columns, data=scaler.transform(X_test))
                    print (' data was scaled using standardScaler')
                
                if classification_problem:
                    if featureSelectionMethod_class is None:
                        featureSelectionMethod=f_classif
                    else:
                        print ('feature selection method= ',featureSelectionMethod_class)
                        if featureSelectionMethod_class == 'f_classif':
                            featureSelectionMethod=f_classif
                        elif featureSelectionMethod_class == 'mutual_info_classif':
                            featureSelectionMethod=mutual_info_classif
                        elif featureSelectionMethod_class == 'chi2':
                            featureSelectionMethod= chi2
                        else:
                            print ('feature selection method unknown')
                else:
                    if featureSelectionMethod_reg is None:
                        featureSelectionMethod=f_regression
                    else:
                        print ('feature selection method= ',featureSelectionMethod_reg)
                        if featureSelectionMethod_reg == 'f_regression':
                            featureSelectionMethod=f_regression
                        elif featureSelectionMethod_reg == 'mutual_info_regression':
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
            if model == 'lin':
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
            if parameter_search_type == 'randomSearch':
                modelObject = RandomizedSearchCV(model(), param_distributions=predictor_params, n_iter=n_random,
                n_jobs=8, cv=5, scoring=scoring)
            elif parameter_search_type == 'gridSearch':
                modelObject = GridSearchCV(model(), param_grid=predictor_params,
                                            n_jobs=8, cv=5, scoring=scoring)
            else:
                model_params = predictor_params
                modelObject = model(**model_params)
            
            modelObject.fit(X_train, y_train)
            
            #save model if required:
            if saveModel:
                #savwe split model:
                modelFileName = '%s_split%s_model.pkl' % (y_name, split)
                f1 = output_dir + modelFileName
                pickle.dump(modelObject, open(f1, "wb"))  
                # save split scaler:
                scalerFileName = '%s_split%s_scaler.pkl' % (y_name, split)
                f2 = output_dir + scalerFileName
                try:
                    pickle.dump(scaler, open(f2, "wb"))
                except:
                    print ('no scaler object to save')  
                # save split selected_features
                selectedFeaturesFileName = '%s_split%s_selectedFeatures.pkl' % (y_name, split)
                f3 = output_dir + selectedFeaturesFileName
                pickle.dump(selected_features, open(f3, "wb"))  # save split model
                       
            # use best predictor according to random hyper-parameter search
            if (parameter_search_type != 'randomSearch') and (parameter_search_type != 'gridSearch'):
                if (classification_problem) and not multiclass:
                    y_pred = modelObject.predict_proba(X_test)
                    final_pred.loc[X_test.index, :] = np.expand_dims(y_pred[:, 1], 1)
                else:
                    y_pred = modelObject.predict(X_test)
                    final_pred.loc[X_test.index, :] = np.expand_dims(y_pred, 1)
                # document best parameters chosen:
                selected_model_params_df = pd.DataFrame()
                # run SHAP or get coeffs
                if not useShortVersion:
                    if model == 'lin':
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
                if (classification_problem) and not multiclass:
                    y_pred = modelObject.best_estimator_.predict_proba(X_test)
                    final_pred.loc[X_test.index, :] = np.expand_dims(y_pred[:, 1], 1)
                else:
                    y_pred = modelObject.best_estimator_.predict(X_test)
                    final_pred.loc[X_test.index, :] = np.expand_dims(y_pred, 1)
                
                
                if not useShortVersion:
                    # document best parameters chosen:
                    df = pd.DataFrame(index=[split], data=modelObject.best_params_)
                    selected_model_params_df = pd.concat([selected_model_params_df, df])
                    # run SHAP or get coeffs
                    if model == 'lin':
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
        results_df = _evaluate_performance(y_name, final_pred.values.ravel(), y, results_df, classification_problem, binary_threshold, nCat)
        predictions_df.loc[final_pred.index, y_name] = final_pred.values.ravel()
#         except:
#             print ('prediction failed with phenotype %s'%y_name)
#             continue
    _save_temporary_files(idx, shap_values_dic,coef_values_dic, model_params_dic, results_df, predictions_df)
    return


# def prediction_bs(q2,command_args):
    

def _save_temporary_files(idx, shap_values_dic, coef_values_dic,model_params_dic, results_df, predictions_df):
    with open(output_dir + '/temp_resdf_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(results_df, fout)
    with open(output_dir + '/temp_shap_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(shap_values_dic, fout)
    with open(output_dir + '/temp_coef_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(coef_values_dic, fout)
    with open(output_dir + '/temp_model_params_' + str(idx) + '.pkl', 'wb') as fout:
        pickle.dump(model_params_dic, fout)
    with open(output_dir + '/temp_pred_' + str(idx) + '.pkl', 'wb') as fout:
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


def concat_outputs():
    print ('concat_outputs')
    all_temp_files = os.listdir(output_dir)
    resdf_files = [output_dir + f for f in all_temp_files if f.startswith('temp_resdf_')]
    shap_files = [output_dir + f for f in all_temp_files if f.startswith('temp_shap_')]
    coef_files = [output_dir + f for f in all_temp_files if f.startswith('temp_coef_')]
    model_params_files = [output_dir + f for f in all_temp_files if f.startswith('temp_model_params_')]
    pred_files = [output_dir + f for f in all_temp_files if f.startswith('temp_pred_')]
    _concat_files(resdf_files, output_dir + '/results_df.pkl', how='dataframe')
    _concat_files(shap_files, output_dir + '/shap_values.pkl', how='dic')
    _concat_files(coef_files, output_dir + '/coef_values.pkl', how='dic')
    _concat_files(model_params_files, output_dir + '/selected_model_params.pkl', how='dic')
    _concat_files(pred_files, output_dir + '/predictions_df.pkl', how='dataframe', axis=1)
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
     



def upload_these_jobs(q):
    print ('upload_these_jobs')
    waiton = []
    try:
        Y = Utils.Load(path_to_Y)
    except:
        Y = pd.read_excel(path_to_Y).set_index('BD')  # enable loading dataframe from excel files
    
    # fill missing target values with median target value:
    if fillMissingY:
        for col in Y.columns:
            try:
                Y[col] = Y[col].fillna(Y[col].median())
            except:
                try:
                    Y[col] = Y[col].fillna(Y[col].mode())
                except:
                    print ('couldnt fill nas')
    
    for idx in range(0, Y.shape[1], n_cols_per_job):
        print ('start running prediction for idx:', idx)
        waiton.append(q.method(model_and_shap, (Y.iloc[:, idx:idx + n_cols_per_job], idx)))
    res = q.waitforresults(waiton)
    # merge the temp results files
    concat_outputs()
    return res




def main():
    print ('main')  
    # check folder existance: 
    if (not os.path.exists(path_to_X)) or (not os.path.exists(path_to_Y)):
        print ("X or Y doesn't exist!"); return
    if n_cols_per_job < 1 or n_cols_per_job > 1000:
        print ("n_cols_per_job must be between 1 and 1000"); return
    #     make_dir_if_not_exists(output_dir)
    if not isdir(output_dir):
        makedirs(output_dir)
        print ('made new output dir')
    
    # save parameters to a text file:
    # print ('saving parameters to file')
    # argFile = output_dir + 'argFile.txt'
    # new_file = open(argFile, mode="w")
    # for arg in vars(command_args):
    #     new_file.write(str(arg) +": "+ str(getattr(command_args, arg))+"\n")
    # new_file.write('NOTE! if -pathToModelParams is not None, relate to the parameters in this file as model paramters!')
    # new_file.close()
    
    if only_concat:
        concat_outputs(); return
       
    # with open(output_dir + '/args' + str(datetime.now()), 'w') as handle:
    #     for arg in vars(command_args):
    #         handle.write(str(arg) + '\t' + str(getattr(command_args, arg)) + '\n')

    with qp(jobname='SHAP-modelObject', q=['himem7.q'], mem_def=str(mem_def) + 'G',
            trds_def=2, tryrerun=True, max_u=350, qworker='~/Develop/Python/lib/SegalQueue/qworker.py') as q:
#         os.chdir("/net/mraid08/export/jafar/Microbiome/Analyses/Noamba/temp_q_dir/")
        tempDir = "/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/temp_q_dir"
        if not isdir(tempDir):
            makedirs(tempDir)
        os.chdir(tempDir)
        q.startpermanentrun()
        upload_these_jobs(q)
    
#     if n_bs !=0:
#         with qp(jobname='bsPrediction', q=['himem7.q'], mem_def=str(mem_def) + 'G',
#             trds_def=4, tryrerun=True, max_u=350) as q2:
#             tempDir = "/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/temp_q_dir"
#             if not isdir(tempDir):
#                 makedirs(tempDir)
#             os.chdir(tempDir)
#             q.startpermanentrun()
#             prediction_bs(q2,command_args)
        
        
        
sethandlers()
main()
