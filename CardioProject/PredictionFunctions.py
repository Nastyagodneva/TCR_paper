import sys
sys.path.insert(0, '/net/mraid08/export/genie/workspace/Microbiome/ShaniBA')

from ShaniBA.myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot
from os import listdir, mkdir, makedirs
from os.path import isfile, join, isdir, exists
import pandas as pd
import numpy as np
from scipy import stats
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from ShaniBA.myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot, draw_correlation_scatter
from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
from Bio.SeqUtils import GC
import seaborn as sns
import random
from scipy.stats import pearsonr, fisher_exact, spearmanr
from skbio.diversity.alpha import shannon, simpson, berger_parker_d
# 
# from pop_organize import get_sample_data, get_sample_with_dfs
# from SufficientStatistics import *
from ShaniBA.MyFunctionsShani import *
import math
# from myplots import roundup, rounddown, find_decimal_fold
from skbio.stats.distance import mantel
from scipy.spatial.distance import braycurtis, pdist
from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import *
from ShaniBA.TCR_feature_generation.publicSeqAnalysis import *
from ShaniBA.TCR_microbiome_interactions.TCR_microbiome_interactions_functions import *

import os
from Utils import cacheOnDisk
from SegalQueue.qp import qp, fakeqp
from addloglevels import sethandlers

# ML imports:
from xgboost import XGBClassifier
import lightgbm as lgb
from collections import OrderedDict
from sklearn.model_selection import GroupKFold, StratifiedKFold, KFold
import statsmodels.formula.api as sm
from sklearn.linear_model import LogisticRegression
import shap
from sklearn import metrics, preprocessing
from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel, SelectKBest, chi2, mutual_info_classif, f_classif, f_regression, mutual_info_regression
from sklearn.naive_bayes import GaussianNB



MyPath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF'




#-------------------------------------------------------------------------------------------------
def predictBinaryByDistMat(Y, YName, X, XName, ResultFolder, model, modelName, model_params, n_splits, useCV, stratifiedCV,
                           featureSelectionMethod, k, cdate=None, overwriteResults=False):
    
    from collections import Counter
    print 'k after calling prediction function is %s' % k
    predResultsDF = pd.DataFrame()
    # (1) arrange result folders:
    print 'arranging  result folders:'
    
    predResultsDFFolder = '%s/%s_predictions/%s/predictionDFs' % (ResultFolder, YName, modelName)  # define folder for all result dfs in this model
    if not isdir(predResultsDFFolder):
        makedirs(predResultsDFFolder)
        print 'generating predResultsDFFolder %s' % predResultsDFFolder      
    predResultsfigFolder = '%s/%s_predictions/%s/figs' % (ResultFolder, YName, modelName)  # define folder for figs in this model
    if not isdir(predResultsfigFolder):
        makedirs(predResultsfigFolder)
        print 'predResultsfigFolder %s' % predResultsfigFolder

    # generate result DF name based on model_params:
    d = OrderedDict(sorted(model_params.items(), key=lambda t: t[0]))
    try:
        d2 = d.copy()
        del d2['numthreads']
    except:
        d2 = d.copy()
    predResultsDFName = '_'.join(['%s%s' % (key.replace('_', ''), value) for (key, value) in d2.items()])  # generate a file name based on params
    predResultsDFName = predResultsDFName.replace('.', '-')
    if model_params == {}:
        predResultsDFName = 'defaultParams'
    if useCV:
        predResultsDFName = '%s_CV%s' % (predResultsDFName, n_splits)
    predResultsDFName = '%s_%s' % (XName, predResultsDFName)
    
    predResultsDFfile = '%s/%s.xlsx' % (predResultsDFFolder, predResultsDFName)
    existingDFs = [f for f in listdir(predResultsDFFolder) if isfile(join(predResultsDFFolder, f))]
#     print 'first existing df names:'
#     print  existingDFs[:10]

    if (predResultsDFName not in existingDFs) | overwriteResults:  
        # common processing of X and y:
        # leave only common samples in each df (X and seq1data)
        X = X.loc[[str(x) for x in X.index.tolist() if x in Y.index], :]
        X = X.sort_index()
        print 'X shape is %s_%s' % (X.shape[0], X.shape[1])
        print 'the 100th sample in X is %s' % X.index[100]
        print X.iloc[:3, :3]
        Y = Y.loc[[str(x) for x in Y.index.tolist() if x in X.index]]
        Y = Y.sort_index()
#         oldColName=str(Y.columns.values[0])
#         Y=Y.rename(columns={oldColName:'Class'})
        print 'Y shape is %s' % (Y.shape[0])
        print 'the 100th sample in Y is %s' % Y.index[100]
        print Y.head(3)

        # (2) model fitting and predictions:
        if useCV:
            print 'splitting train_test using cross validation with %s splits...' % n_splits
            
            if stratifiedCV:
                group_kfold = StratifiedKFold(n_splits=n_splits)
                groups = None
            else:
                group_kfold = GroupKFold(n_splits=n_splits)
                groups = np.array(range(X.shape[0]))
                
        
            y_pred_df = pd.DataFrame(index=Y.index, columns=['pred_proba'])
            i = 0
            selected_features_list = []
            for train_index, val_index in group_kfold.split(X, Y, groups):
                print i
                i += 1
                X_train, X_val = X.iloc[train_index], X.iloc[val_index]
                y_train, y_val = Y.loc[X_train.index], Y.loc[X_val.index]
                
                print 'fraction of 1s in train set=%s' % (float(y_train.sum()) / len(y_train))
                print 'fraction of 1s in test set=%s' % (float(y_val.sum()) / len(y_val))
                
                if (featureSelectionMethod is not None) & (featureSelectionMethod != 'usingModel'):
          
                    if k is None:
                        k = 50  
                    if len(X_train.columns) < k:
#                         k = len(X_train.columns)
                        featureSelectionMethod = None
                        print 'k is bigger than the number of columns, therefore no column was filtered by univariate'
                        selected_features = X_train.columns.values.tolist()
                        print selected_features
                    else:
                        print 'filtering by %s with k=%s' % (featureSelectionMethod, k)
                        if featureSelectionMethod == 'random':
                            import random
                            idxs_selected = random.sample(X_train.columns.values.tolist(), k)
                            X_train_new = X_train.loc[:, idxs_selected]
                            X_val_new = X_val.loc[:, idxs_selected]
#                             uniName='random'
                        else:
                            selector = SelectKBest(featureSelectionMethod, k)
                            selector.fit(X_train, y_train)
                            idxs_selected = selector.get_support(indices=True)
                            # Create new dataframe with only desired columns, or overwrite existing
                            X_train_new = X_train.iloc[:, idxs_selected]
                            X_val_new = X_val.iloc[:, idxs_selected]
#                             uniName = '%s' % filterFeaturesByUnivariate
#                             uniName = uniName.split(' ')[1]
#                             for ch in ['<', '>', '.']:
#                                 uniName = uniName.replace(ch, '')
                        X_train = X_train_new
                        X_val = X_val_new
                        print 'k=%s' % k
                        print 'length of idxs_selected=%s' % len(idxs_selected)
                        print 'X_train shape after univaraite feature selection is  %s_%s' % (X_train.shape[0], X_train.shape[1])
                        print 'X_val shape after univaraite feature selection is  %s_%s' % (X_val.shape[0], X_val.shape[1])
                        print 'selected columns are:'
                        selected_features = X_train.columns.values.tolist()
                        print selected_features
#                         XName = '%s_fUni%s%s' % (XName, uniName, k)
                
                
                # creating the model object
                m = model(**model_params)

                # fitting the training
#                 if modelName == 'LGBMClassifier':
#                     m.fit(X_train, y_train,
# #                             eval_set=[(X_val, y_val)],
#                             early_stopping_rounds=None, verbose=-1)
#                 else:
                m.fit(X_train, y_train)
                
                if featureSelectionMethod == 'usingModel':
                    
                    features = X_train.columns.tolist()
                    importances = m.feature_importances_
                    
#                     #plot feature importances:
#                     plt.subplots(figsize=(40,6))
#                     plt.bar(range(len(features)),importances)
#                     plt.xticks ([x+0.4 for x in range(len(features))],features,fontsize=7,rotation=45)
#                     plt.title('Feature importance')
#                     plt.show()
                    
                    # take only important fefatures:
                    indices = (np.argsort(importances)[::-1])[:k]
                    X_train_new = X_train.iloc[:, indices]
                    X_val_new = X_val.iloc[:, indices]
                    
                    selected_features = X_train_new.columns.tolist()
                    print 'selected features are:'
                    print selected_features
                    
                    # fit model again, now using only selected features:
                    m2 = model(**model_params)

                    # fitting the training
#                     if modelName == 'LGBMClassifier':
#                         m2.fit(X_train_new, y_train,
#     #                             eval_set=[(X_val, y_val)],
#                                 early_stopping_rounds=None, verbose=-1)
#                     else:
                    m2.fit(X_train_new, y_train)
                    
                    print 'X_train shape after feature selection using model is  %s_%s' % (X_train_new.shape[0], X_train_new.shape[1])
                    print 'X_val shape after feature selection using model is  %s_%s' % (X_val_new.shape[0], X_val_new.shape[1])
                
                
                    # getting the predictions for the test using m2 (trained on select feature only!), X_train_new
                    # and X_val_new
                    y_pred_proba = m2.predict_proba(X_val_new)
                    y_pred_df.loc[y_val.index, :] = np.expand_dims(y_pred_proba[:, 1], 1)
                    numFeaturesPerSplit = len(X_train_new.columns)
                
                else:
                    # getting the predictions for the test using m1 (trained on all features or features pre-selected
                    # before applying the model:
                    y_pred_proba = m.predict_proba(X_val)
                    y_pred_df.loc[y_val.index, :] = np.expand_dims(y_pred_proba[:, 1], 1)
                    numFeaturesPerSplit = len(X_train.columns)
                
                # appending list of selected features:
                try:
                    selected_features_list = selected_features_list + selected_features
                except:
                    print 'no features were selected'
                    
                
        else:
            test_size = 0.33
            print 'using normal train_test split with test_size=%s' % test_size
            from sklearn.model_selection import train_test_split  
            # seed = 7
            
            X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=test_size, stratify=Y)
            print 'fitting model...'
            model = model
            model.fit(X_train, y_train)
            # make predictions for test data
            print 'predicting...'
            y_pred = model.predict(X_test)
            predictions = [round(value) for value in y_pred]
            # evaluate predictions
            from sklearn.metrics import accuracy_score
            accuracy = accuracy_score(y_test, y_pred)
            print("Accuracy: %.2f%%" % (accuracy * 100.0))
            y_pred_df = pd.DataFrame(index=Y.index, data={'pred':y_pred})
            
            
        # generating selected feature counter and transforming into string:
        selectedFeatureCounter = Counter(selected_features_list)
        try:
            selectedFeatureCounterString = ''
            for k in sorted(selectedFeatureCounter.items(), key=lambda x:x[1], reverse=True):
                s1 = ':'.join([str(k[0]), str(k[1])])
#                 print s1
                selectedFeatureCounterString = ','.join([selectedFeatureCounterString, s1])
#             print selectedFeatureCounterString
        except:
            print 'feature counter is empty'
            selectedFeatureCounterString = ''
            

        # (4) plots:
        print 'generating plots...'

        # this plot shows the probabilities returned by the predictor colored by the class
#         plt.figure(figsize=(3, 2))
#         plt.scatter(range(y_pred_df.shape[0]), y_pred_df.pred_proba, c=Y)
#         plt.show()

        # plot ROC and PR curves

        fpr, tpr, thresholds = metrics.roc_curve(Y + 1, y_pred_df.pred_proba, pos_label=2)
        roc_auc = metrics.auc(fpr, tpr)
        print roc_auc

        plt.figure(figsize=(20, 7))
        plt.subplot(1, 2, 1)
        lw = 2
        plt.plot(fpr, tpr, color='darkorange', lw=lw, label='ROC curve (area = %0.3f)' % roc_auc)
        plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate', fontsize=15)
        plt.ylabel('True Positive Rate', fontsize=15)
        plt.title('ROC curve - ' , fontsize=20)
        plt.legend(loc="lower right", fontsize=15)

        plt.subplot(1, 2, 2)
        precision, recall, _ = metrics.precision_recall_curve(Y, y_pred_df.pred_proba)
        pr_auc = metrics.auc(recall, precision)

        plt.step(recall, precision, color='b', alpha=0.2, where='post')
        plt.fill_between(recall, precision, step='post', alpha=0.5,
                         color='darkorange', label='Precision Recall curve - AUC = {0:0.3f}'.format(pr_auc))
        plt.plot([0, 1], [Y.sum() / Y.shape[0], Y.sum() / Y.shape[0]], color='navy', lw=lw, linestyle='--')
        plt.xlabel('Recall', fontsize=15)
        plt.ylabel('Precision', fontsize=15)
        plt.ylim([0.0, 1.05])
        plt.xlim([0.0, 1.0])
        plt.title('Precision Recall curve', fontsize=20)
        plt.legend(loc="upper right", fontsize=15)

        if roc_auc > 0.55:
            predResultsFigfile = '%s/%s_ROC-PR.png' % (predResultsfigFolder, predResultsDFName)
            plt.savefig(predResultsFigfile, bbox_inches='tight')

        plt.show()

        # (5) generate summarizing df:
        print 'generating summarizing df'
        predResultsDF.loc[0, 'Yname'] = YName
        predResultsDF.loc[0, 'Xname'] = XName
        predResultsDF.loc[0, 'modelName'] = modelName
        predResultsDF.loc[0, 'roc_auc'] = round(roc_auc, 3)
        predResultsDF.loc[0, 'pr_auc'] = round(pr_auc, 3)
        predResultsDF.loc[0, 'perc_pos_target'] = round(float(Y.sum()) / len(Y), 3)
        predResultsDF.loc[0, 'pr_auc_corrected'] = round(pr_auc, 3) - round(float(Y.sum()) / len(Y), 3)
        predResultsDF.loc[0, 'useCV'] = useCV
        predResultsDF.loc[0, 'featureSelectionMethod'] = featureSelectionMethod
#         predResultsDF.loc[0, 'C'] = C
        predResultsDF.loc[0, 'features used'] = selectedFeatureCounterString
        predResultsDF.loc[0, 'total num features used'] = len(selectedFeatureCounter)
        predResultsDF.loc[0, 'Num features used per split'] = numFeaturesPerSplit
        
        if useCV:
            predResultsDF.loc[0, 'n_splits'] = n_splits
        
        
        for (key, value) in OrderedDict(sorted(model_params.items(), key=lambda t: t[0])).items():
            predResultsDF.loc[0, key] = value
        
        nPos = int(Y.sum())
        nNeg = len(Y) - nPos
        
        print nPos, nNeg

        predResultsDF.loc[0, 'nPos'] = nPos
        predResultsDF.loc[0, 'nNeg'] = nNeg
        predResultsDF.loc[0, 'cdate'] = cdate
        
        
        predResultsDF.to_excel(predResultsDFfile)
        print predResultsDFfile
        
    else:
        print 'this prediction already exists'
    
    return predResultsDF
    
#---------------------------------------------------------------------------------------------------------------
'''
the following function is adequate for binary classification using logReg/lightGBM classifier or XGBoost classifier

input:
featureSelectionMethod - None,chi2,f_classif,mutual_info_classif,'random','usingModel'

continue documentation



'''
def predictionPipeline(datasetType, ResultFolder, model, modelName, model_params, n_splits, useCV, stratifiedCV,
                      XName, usePhenotype, useTCRdf, useTCRfeatures, usePCAdf,
                      YName, targetDF,
                      genPhenotypeDF, phenotypeDF, phenotypeDFname,
                      datasetFolder, datasetName, TCRdf, TCRdfName, n_comp, isSparse,
                      getTCRfeatures, TCRfeatureDF, sampleListName, sampleList,
                      filterFeaturesByCorr, featureSelectionMethod, scale=False,
                      k=50, cdate=None, overwriteResults=False):  
    
    
    # (1) get y:
    y = targetDF[YName]
    y = editSampleNames(y)
    print 'y shape is %s' % y.shape
    print 'y head:'
    print y.head()
    print y.tail() 
    
    # (2) get phenotypeDF:
    if usePhenotype:
        if genPhenotypeDF:
            # # generate phenotype df if necessary
            print 'generating phenotypeDF...'
            # ## get the original phenotype file(s) - concat if necessary:
            
            PNPphenotypeFile = '%s/TCR_real_data/NewPhenotypicData/PNP530_phen_new.xlsx' % MyPath
            CardioPhenotypeFile = '%s/TCR_real_data/CardioSamples/phenotypicData/Cardio126phenNew.xlsx' % MyPath

            PNPphenotypeDF = pd.read_excel(PNPphenotypeFile).set_index('BD')
            CardiophenotypeDF = pd.read_excel(CardioPhenotypeFile).set_index('BD')

            if datasetType == 'PNP':
                phenotypeDF = PNPphenotypeDF
            elif datasetType == 'PNP_Cardio':
                phenotypeDF = pd.concat([PNPphenotypeDF, CardiophenotypeDF])
            
            phenotypeDF = phenotypeDF.drop(['RegistrationCode', 'DM', 'Diagnosis', 'eGFR by CKD-EPI', 'Known CAD', 'Previous PCI', 'Previous CABG',
                                              'LVEF', 'Initial Troponin', 'Maximal Troponin'], axis=1)
            if datasetType == 'Cardio':
                phenotypeDF = CardiophenotypeDF

#             phenotypeDF['isCardio']=np.where(phenotypeDF.index.str.replace('BD','').astype(int)>949,1,0)
            print 'phenotypeDF shape is %s_%s' % (phenotypeDF.shape[0], phenotypeDF.shape[1])
            print phenotypeDF.tail()
            
            # ## generate dummy variables:
            # the function gen_dummies was copied to Feature_phenotype_functions.py
            toDummyColList = ['Gender', 'Smoking', 'PCRplate']
            df = phenotypeDF.copy()
            phenotypeDF = gen_dummies(df, toDummyColList)
            print 'pehonotypeDF head after dummy addition:'
            phenotypeDF.head()
            
            # save phenotypeDF:
            if datasetType == 'PNP':
                f2 = '%s/TCR_real_data/NewPhenotypicData/PNP530_phen_new_withDummies.xlsx' % MyPath
            elif datasetType == 'PNP_Cardio':
                f2 = '%s/TCR_real_data/PNP530Cardio126Combined/Phenotypes/PNP530Cardio126_phen_new_withDummies.xls' % MyPath
            elif datasetType == 'PNP_Cardio':
                f2 = '%s/TCR_real_data/CardioSamples/NewPhenotypicData/Cardio126_phen_new_withDummies.xlsx' % MyPath
           
            phenotypeDF.to_excel(f2)
#         else:
#             print 'loading phenotype file...'
#             phenotypeDF = pd.read_excel(phenotypefile).set_index('BD')
#             print phenotypeDF.shape
#             phenotypeDF.tail()
            
        #
        # Xtypes:
        allNum = ['Age', 'BMI', 'Creatinine', 'isCardio', 'nTemplates', 'Gender_Male', 'Smoking_Past', 'Smoking_Yes',
               'PCR_Plate1', 'PCR_Plate10', 'PCR_Plate2', 'PCR_Plate3',
               'PCR_Plate4', 'PCR_Plate5', 'PCR_Plate6', 'PCR_Plate7',
               'PCR_Plate8', 'PCR_Plate9']
        small = [ 'Age', 'BMI', 'Creatinine', 'isCardio', 'Gender_Male', 'Smoking_Past', 'Smoking_Yes']
        smallNoCardio = ['Age', 'BMI', 'Creatinine', 'Gender_Male', 'Smoking_Past', 'Smoking_Yes']
        allNumNoPlate = ['Age', 'BMI', 'Creatinine', 'isCardio', 'nTemplates', 'Gender_Male', 'Smoking_Past', 'Smoking_Yes']
        newSmall = ['Age', 'BMI', 'Gender_Male', 'Smoking_Past', 'Smoking_Yes', 'WBC', 'eGFR_CKD-EPI_new']
        cardioFull = ['Age', 'BMI', 'Gender_Male', 'Smoking_Past', 'Smoking_Yes', 'WBC', 'eGFR_CKD-EPI_new', 'CRP',
                    'Hypertension', 'Known CAD', 'PreviousPCI_binary']
        septFeatures = ['Smoking Status_Past', 'Smoking Status_Yes', 'Gender_Male', 'Age', 'BMI', 'eGFR by CKD-EPI',
                       'HbA1C', 'AST', 'Hemoglobin', 'WBC', 'HDL', 'Total Cholesterol']
        if phenotypeDFname == 'small':
            Xcols = small
        elif phenotypeDFname == 'allNum':
            Xcols = allNum
        elif phenotypeDFname == 'smallNoCardio':
            Xcols = smallNoCardio
        elif phenotypeDFname == 'allNumNoPlate':
            Xcols = allNumNoPlate
        elif phenotypeDFname == 'newSmall':
            Xcols = newSmall
        elif phenotypeDFname == 'cardioFull':
            Xcols = cardioFull
        elif phenotypeDFname == 'septFeatures':
            Xcols = septFeatures
        phenotypeDF = phenotypeDF[Xcols]          
        phenotypeDF = editSampleNames(phenotypeDF)
        
        if YName in phenotypeDF.columns:
            print 'the target %s is included in phenotypeDF columns and therefore is being removed...' % YName
            phenotypeDF = phenotypeDF.drop(YName, axis=1)
        print 'fillnas with column medians in phenotypeDF..'
        for col in phenotypeDF.columns:
            phenotypeDF[col] = phenotypeDF[col].fillna(phenotypeDF[col].median())
        print 'final phenotype shape is %s_%s' % (phenotypeDF.shape[0], phenotypeDF.shape[1])
        print 'phenotypeDF head:'
        print phenotypeDF.head()

        
#     # (3) get TCRdf:
#     if useTCRdf or usePCAdf:
#         # # ** load shared sequence matrix
#         print 'loading shared sequence matrix:'
#         TCRdfShortName = TCRdfName.replace('sharingMatrix_PNP530Cardio126MatchedSamples_', '')
#         TCRdfShortName = TCRdfShortName.replace('__', '_')
# 
#         f1 = '%s/sharingAnalysis/%s' % (datasetFolder, TCRdfName)
#         TCRdf = pd.read_pickle(f1)
#         print 'TCRdf shape is %s_%s' % (TCRdf.shape[0], TCRdf.shape[1])
#         print 'TCRdf head:'
#         print TCRdf.iloc[:4, :4]
    
    # (4) get PCAdf:
    if usePCAdf:
        # # **generate PCs as features:
        if n_comp is None:
            n_comp = 5
        print 'generating PCAdf with %s PCs based on TCRdf:' % n_comp
        fig, ax = plt.subplots()
        TCRdf = editSampleNames(TCRdf)
        TCRdf = TCRdf.loc[[str(x) for x in TCRdf.index.tolist() if x in y.index], :]
        PCAdf, ax, p_ttest_PC1, p_ttest_PC2 = PCAfunc(TCRdf, n_comp, isSparse, ax, groupingByDF=y, groupbyName=YName)
        for c in ['BDindex', YName]:
            try:
                PCAdf = PCAdf.drop(c, axis=1)
            except:
                print 'PCAdf doesnt include %s column' % c
        print 'PCAdf shape is %s_%s' % (PCAdf.shape[0], PCAdf.shape[1])
#         print 'PCAdf HEAD:'
#         print PCAdf.head()
#       
    # (5) get TCRfeature:
    if useTCRfeatures:
        if getTCRfeatures:
            print 'generating new feature DF...this will take a while and will print out many things...'
            # ## calculate features - from PopulationAnalaysis_new_version.ipynb
            data_folder = datasetFolder.replace('%s/' % MyPath, '')
            dfs_folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/SamplesForAnalysis_corrected' % data_folder
            filenames = [f for f in listdir(dfs_folder) if isfile(join(dfs_folder, f))]
            filenames = [f.strip('.tsv') for f in filenames]
            filenames = [f.strip('.xlsx') for f in filenames]
            print 'number of samples to extract features is %s' % len(filenames)
            for n, sample_name in enumerate(filenames): 
                print  n, sample_name
                gen_descriptive_stats(sample_name, data_folder, newColumnList=None)
                gen_geneUsageCount(sample_name, data_folder, newColumnList=None)
            # ## generate feature DFs - based on 'Generate features DF' notebook
            seqTypeList = ['Total', 'Prod', 'nonProd']
            dfs_folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/Total' % data_folder
            FeatureGroups = listdir(dfs_folder)
            for seqType in seqTypeList:
                print seqType
                gen_featureSummaryDF_forSeqType(seqType, data_folder, datasetName, FeatureGroups)
                
            data_folder = datasetFolder.replace('%s/' % MyPath, '')
            gen_merged_feature_df_for_dataset(data_folder, datasetName)
            f2 = '%s/featureSummaryDFs/%s_allFeatures' % (datasetFolder, datasetName)
            TCRfeatureDF = pd.read_pickle(f2)
        else:
            if TCRfeatureDF is None:
                print 'loading featureDF...'
                f2 = '%s/featureSummaryDFs/%s_allFeatures' % (datasetFolder, datasetName)
                TCRfeatureDF = pd.read_pickle(f2)
            else:
                TCRfeatureDF = TCRfeatureDF
        print 'TCRfeatureDF shape is %s_%s' % (TCRfeatureDF.shape[0], TCRfeatureDF.shape[1])
        print 'TCRfeatureDF head:'
        print TCRfeatureDF.iloc[:4, :4]
        TCRfeatureDF = editSampleNames(TCRfeatureDF)
        
        # marked next rows as I use featureDF which was already filtered for too-many nans...    
#         # # ** remove columns with too many nans
#         print 'removing feature columns with less than 5 samples that are not nan...'
#         print 'number of such columns is %s' %len(TCRfeatureDF.loc[:, TCRfeatureDF.notnull().sum() <=5].columns.tolist())
#         TCRfeatureDF=TCRfeatureDF.loc[:, TCRfeatureDF.loc[:, TCRfeatureDF.notnull().sum() > 5].columns.tolist()]
#         print 'TCRfeatureDF shape after removing all columns that dont have at least 5 valid values is %s_%s' %(TCRfeatureDF.shape[0],  
#                                                                                                                 TCRfeatureDF.shape[1])
        # # ** fillna's in featureDF:                                                                                                      TCRfeatureDF.shape[1])
#         # get number of columns with number of nan values?
#         print 'fill nans in featureDF...:'
#         print 'columns with top numbers of nan values:'
#         print TCRfeatureDF.isnull().sum().sort_values(ascending=False).head()
#         print 'value counts of nan numbers per column:'
#         print TCRfeatureDF.isnull().sum().value_counts(dropna=False).sort_values(ascending=False).head()
#         # find how many columns include nans:
#         print 'number of columns with nan values = %s' % len(TCRfeatureDF.loc[:, TCRfeatureDF.isnull().sum() > 0].columns.tolist())
#         print 'total number of columns is %s' % len(TCRfeatureDF.columns)
        # get list of columns that inlcude nans and should be fillna-ed with 0s:
        columnsToFill = TCRfeatureDF.loc[:, TCRfeatureDF.isnull().sum() > 0].columns.tolist()
        columnsToFill_clean = [x for x in columnsToFill if 'norm' not in x]
        print 'number of columns that include nans and their name doesnt include -norm- is %s' % len(columnsToFill_clean)
        print 'these columns will be fillna-ed with 0s'
        # fillna with zeros
        TCRfeatureDF[columnsToFill_clean] = TCRfeatureDF[columnsToFill_clean].fillna(0)
        # check how many columns remained with nans:
        print 'now the number of columns with nan values is %s' % len(TCRfeatureDF.loc[:, TCRfeatureDF.isnull().sum() > 0].columns.tolist())
     
    # (6) edit X- filter with sampleList, filter nan columns, fillna, filter for y rows
    
    # build X from its components:
    print 'building X from its components:'
    XcomponentList = []
    if usePhenotype:
        XcomponentList.append(phenotypeDF)
    if useTCRdf:
        XcomponentList.append(TCRdf)
    if useTCRfeatures:
        XcomponentList.append(TCRfeatureDF)
    if usePCAdf:
        XcomponentList.append(PCAdf)
#         XName = '%s_PCA%s' % (XName, n_comp)
    print 'number of Xcomponents to be used is %s' % len(XcomponentList)
    
    # edit sample names and merge all X components:
    for n, comp in enumerate(XcomponentList):
        # comp = editSampleNames(comp)  ## this row was marked because all components were already edited. unmark if necessary
        if n == 0:
            X = comp
        else:
            X = pd.merge(X, comp, how='inner', right_index=True, left_index=True)
    print 'combined X shape is %s_%s' % (X.shape[0], X.shape[1])
    print len(X.index.tolist())
    
    
    # if model is 'logit', drop columns that contains nans, except for phenotype columns; fillna in phenotype columns
    # and check for rows with nan values:

    if ('log' in modelName) or ('Log' in modelName) or ('NB' in modelName) or\
((featureSelectionMethod is not None) & (featureSelectionMethod != 'usingModel')):
        print 'since logRef model or univaraite feature selection is used, needs to fillnas in phenotype columns...'
        nullSum = X.isnull().sum()
        ignoreColumns = nullSum[nullSum > 0].index.tolist()
        for column in phenotypeDF.columns.values:
            if column in ignoreColumns:
                ignoreColumns.remove(column)
        X = X.drop(ignoreColumns, axis=1)
        print 'X shape after dropping nan columns is %s_%s' % (X.shape[0], X.shape[1])
        for column in phenotypeDF.columns.values:
            X[column] = X[column].fillna(X[column].median())
        X = X.dropna(how='any')
        print 'Xshape after dropping rows with nans=%s_%s' % (X.shape[0], X.shape[1])
    
    # filter X rows by sampleList:
    if sampleListName is not None:
        print 'filtering X rows with sample List = %s' % sampleListName
        X = X.loc[sampleList, :]
        print 'Xshape after filtering rows with sample list = %s_%s:' % (X.shape[0], X.shape[1])
    # common processing of X and y:
    # leave only common samples in each df (X and seq1data)
    X = X.loc[[str(x) for x in X.index.tolist() if x in y.index], :]
    X = X.sort_index()
    print 'X shape is %s_%s' % (X.shape[0], X.shape[1])
    print 'the 100th sample in X is %s' % X.index[100]
    print X.iloc[:3, :3]
    y = y.loc[[str(x) for x in y.index.tolist() if x in X.index]]
    y = y.sort_index()
#         oldColName=str(Y.columns.values[0])
#         Y=Y.rename(columns={oldColName:'Class'})
    print 'Y shape is %s' % (y.shape[0])
    print 'the 100th sample in Y is %s' % y.index[100]
    print y.head(3)
#     # drop samples that are not in y
#     print 'dropping X rows that do not appear in y:'
#     ysamples = y.index.tolist()
#     X = X.loc[ysamples, :]
#     X = X.dropna(how='all')
    print 'Xshape after dropping rows that do not appear in y = %s_%s:' % (X.shape[0], X.shape[1])    
       
    # (7) processing required for feature selection (which will occure within the CV loops): 
    # # naming change occures here:
    if not 'phenotypesOnly' in XName:
        if featureSelectionMethod == 'random':
            uniName = 'random'
        else:
            uniName = '%s' % featureSelectionMethod
            try:
                uniName = uniName.split(' ')[1]
            except:
                uniName = uniName
            for ch in ['<', '>', '.']:
                uniName = uniName.replace(ch, '')
        XName = '%s_FSM%s%s' % (XName, uniName, k)
    
    
    # # need to transfer the following into the cv
    if filterFeaturesByCorr is not None:
        
        print 'filtering out columns with correlation larger than %s' % filterFeaturesByCorr
        colToDrop = []
        for i in range(len(X.columns)):
            for j in range(i + 1, len(X.columns)):
                col1 = X.columns.values[i]
                col2 = X.columns.values[j]
                r, p = MyPearsonr(X[col1], X[col2])
                if r > filterFeaturesByCorr:
                    print col1, col2, r
                    colToDrop.append(col2)
        X = X.drop(colToDrop, axis=1)
        print 'X shape after dropping highly correlated columns=%s_%s' % (X.shape[0], X.shape[1])
        XName = '%s_fCorr%s' % (XName, filterFeaturesByCorr) 
        
    if featureSelectionMethod == chi2:  # using chi2 required scaled data and non-negative values only (thus
                                         # PCA data was omitted (contains negatives) and the data is scaled here
        scale = True
    
    if scale:
        print 'scaling all non-binary features'
        print 'X shape before scaling is  %s_%s' % (X.shape[0], X.shape[1])
        colsToScale = []
        for col in X.columns:
            if len(X[col].value_counts()) > 2:
                X[col] = X[col].fillna(X[col].median())
                colsToScale.append(col)
        X[colsToScale] = preprocessing.scale(X[colsToScale])
        print 'X shape after scaling is  %s_%s (shouldnt change!)' % (X.shape[0], X.shape[1])
        XName = '%s_scaled' % XName
    
    # (8) model prediction: yield df, plots
    print 'k before calling prediction function is %s' % k
    print 'calculating model prediction'
    Y = y
    predResultsDF = predictBinaryByDistMat(Y, YName, X, XName, ResultFolder, model, modelName, model_params, n_splits,
                                         useCV, stratifiedCV, featureSelectionMethod, k, cdate, overwriteResults)
    print predResultsDF
    print 'done!!'


#------------------------------------------------------
'''
the following function takes an excel file which contains the following columns: Model,Features,roc_auc,pr_auc_corrected
and plots the ROC and PR AUC values and change of values for different models 

'''

def plot_modelComparison(resultDF, title, folderToSave=None):
    resultDF['roc_auc'] = resultDF['roc_auc'].astype(float)
    resultDF = resultDF.dropna(how='any')
#     print resultDF
    fig, ax = plt.subplots(figsize=(4, 7))
    colorList = ['r', 'b', 'g', 'y']
    x1 = 0
    x2 = 0.4
    sumDF = pd.DataFrame()
    for n, model in enumerate(resultDF.Model.unique()):
        
        print n, model
        df = resultDF[resultDF.Model == model].drop('Model', axis=1).set_index('Features')
        ROC = list(df.loc[:, 'roc_auc'])
        print ROC
        PR = list(df.loc[:, 'pr_auc_corrected'])
        print PR

        if ROC[1] < ROC[0]:
            marker1 = 'v'
        else:
            marker1 = '^'
        if PR[1] < PR[0]:
            marker2 = 'v'
        else:
            marker2 = '^'

        ax.axvline(x=x1 + 0.05 * (n + 1), ymin=ROC[0], ymax=ROC[1], lw=4, marker=marker1, ms=10, color=colorList[n], label=model)
        ax.axvline(x=x2 + 0.05 * (n + 1), ymin=PR[0], ymax=PR[1], lw=4, marker=marker2, ms=10, color=colorList[n])

        ROCdif = ROC[1] - ROC[0]
        PRdif = PR[1] - PR[0]
        

        sumDF.loc[n, 'model'] = model
        sumDF.loc[n, 'ROCdif'] = ROCdif
        sumDF.loc[n, 'PRdif'] = PRdif

    ax.set_xticks([0.1, 0.5])
    ax.set_xticklabels(['roc_auc', 'pr_auc_corrected'], fontsize=14)
    
#     ax.set_xticks([0.1,0.5],['roc_auc','pr_auc_corrected'],fontsize=20)
    ax.legend()
    ax.grid(False)
    ax.set_xlim(0, 0.7)
    
    # ax.set_xlabel('Metric')
    ax.set_title(title)
    if folderToSave is not None:
        f1 = '%s/%s' % (folderToSave, title)
    
    plt.show()
    print sumDF
    return fig, sumDF


#--------------------------------------------------------------------
'''
the following function takes featureDF (can generate it for matched sample cohorts if necessary) and removes 
redundant columns:
* columns with less than 5 samples which are not nans.
* columns with complete correlation to another column

it outputs the list of columns to drop, and the featureDF after removal of these columns
'''





def removeRedundantFeatures(ss, repeat, filterFeaturesByCorr, sampleList, sampleListName=None, genFeatureDF=True):

    if genFeatureDF:
        # generate featureDF:
        print 'getting TCRfeature DF...'
        if ss is None:
            PNPfeatureFile = '%s/TCR_real_data/featureSummaryDFs/PNP530_allFeatures' % MyPath
        else:
            PNPfeatureFile = '%s/TCR_real_data/PNP530_SubSampled%sdata_rep%s/featureSummaryDFs/\
PNP530_ss%s_rep%s_allFeatures' % (MyPath, ss, repeat, ss, repeat)
        PNPfeatureDF = pd.read_pickle(PNPfeatureFile)
        if ss is None:
            CardiofeatureFile = '%s/TCR_real_data/CardioSamples/featureSummaryDFs/Cardio126_allFeatures' % MyPath
        else:
            CardiofeatureFile = '%s/TCR_real_data/CardioSamples/Cardio126_SubSampled%sdata_rep%s/featureSummaryDFs/\
Cardio126_ss%s_rep%s_allFeatures' % (MyPath, ss, repeat, ss, repeat)
        CardiofeatureDF = pd.read_pickle(CardiofeatureFile)
        
        if ss is None:
            datasetFolder = '%s/TCR_real_data/PNP530Cardio126Combined' % MyPath
            datasetName = 'PNP530Cardio126'
        else:
            datasetFolder = '%s/TCR_real_data/PNP530Cardio126Combined/MatchedSamples/ss%srep%s/' % (MyPath, ss, repeat)
            datasetName = 'MatchedSamples_ss%srep%s' % (ss, repeat)

        TCRfeatureDF = pd.concat([PNPfeatureDF, CardiofeatureDF])
        TCRfeatureDF = editSampleNames(TCRfeatureDF)
        # TCRfeatureDF.head()
        
        if sampleListName is None:
            featureFileXls = '%s/featureSummaryDFs/%s_allFeatures.xlsx' % (datasetFolder, datasetName)
            featureFilePickle = '%s/featureSummaryDFs/%s_allFeatures' % (datasetFolder, datasetName)
        else:
            featureFileXls = '%s/featureSummaryDFs/%s_%s_allFeatures.xlsx' % (datasetFolder, sampleListName, datasetName)
            featureFilePickle = '%s/featureSummaryDFs/%s_%s_allFeatures' % (datasetFolder, sampleListName, datasetName)
            TCRfeatureDF_matched = TCRfeatureDF.loc[sampleList, :]
        print 'TCR feature DF shape before any removal and saving is %s_%s' % (TCRfeatureDF_matched.shape[0], TCRfeatureDF_matched.shape[1])
        TCRfeatureDF_matched.to_excel(featureFileXls)
        TCRfeatureDF_matched.to_pickle(featureFilePickle)
 
    # load featureDF:
    try:
        featureFilePickle = '%s/featureSummaryDFs/%s_%s_allFeatures' % (datasetFolder, sampleListName, datasetName)
        TCRfeatureDF = pd.read_pickle(featureFilePickle)    
    except:
        print 'still trying to load feature file...'
        try:
            featureFilePickle = '%s/featureSummaryDFs/%s_allFeatures' % (datasetFolder, datasetName)
            TCRfeatureDF = pd.read_pickle(featureFilePickle)
        except:
            print 'couldnt load feature file...'
    
    if sampleListName is None:
        sampleListName = ''
    featureFilePickle2 = '%s/featureSummaryDFs/%s%s_allFeatures_noCorrelated' % (datasetFolder, datasetName, '_filteredBy' + sampleListName)
    featureFileExcel2 = '%s/featureSummaryDFs/%s%s_allFeatures_noCorrelated.xlsx' % (datasetFolder, datasetName, '_filteredBy' + sampleListName)

    # # ** remove columns with too many nans
    print 'TCR feature DF shape after loading and before any removal  is %s_%s' % (TCRfeatureDF.shape[0], TCRfeatureDF.shape[1])
    print 'removing feature columns with less than 5 samples that are not nan...'
    nanColumns = TCRfeatureDF.loc[:, TCRfeatureDF.notnull().sum() <= 5].columns.tolist()
    print nanColumns
    print 'number of such columns is %s' % len(TCRfeatureDF.loc[:, TCRfeatureDF.notnull().sum() <= 5].columns)
    TCRfeatureDF = TCRfeatureDF.loc[:, TCRfeatureDF.loc[:, TCRfeatureDF.notnull().sum() > 5].columns.tolist()]
    print 'TCRfeatureDF shape after removing all columns that dont have at least 5 valid values is %s_%s'\
    % (TCRfeatureDF.shape[0], TCRfeatureDF.shape[1])

    # generating list of columns to drop with correlation higher or equal to threshold                                                                                                        
    print 'TCRfeatureDF2 shape before dropping highly correlated columns=%s_%s' % (TCRfeatureDF.shape[0], TCRfeatureDF.shape[1])

    TCRfeatureDF2 = TCRfeatureDF.copy()
    print 'filtering out columns with correlation larger than %s' % filterFeaturesByCorr
    colToDrop = []
    for i in range(len(TCRfeatureDF2.columns)):
        for j in range(i + 1, len(TCRfeatureDF2.columns)):
            col1 = TCRfeatureDF2.columns.values[i]
            col2 = TCRfeatureDF2.columns.values[j]
            r, p = MyPearsonr(TCRfeatureDF2[col1], TCRfeatureDF2[col2])
            if r >= filterFeaturesByCorr:
                print col1, col2, r
#                 print col1, col2, r, len(TCRfeatureDF2[TCRfeatureDF2[col1].notnull()]),len(TCRfeatureDF2[TCRfeatureDF2[col2].notnull()])
                colToDrop.append(col2)
                
    TCRfeatureDF3 = TCRfeatureDF2.drop(colToDrop, axis=1)
    print 'TCRfeatureDF3 shape after dropping highly correlated columns=%s_%s' % (TCRfeatureDF3.shape[0], TCRfeatureDF3.shape[1])
    # TCRfeatureDF2Name = '%s_fCorr%s' % (TCRfeatureDF2Name, filterFeaturesByCorr)
    
    colToDrop2 = list(set(colToDrop))
    colToDrop3 = colToDrop2 + nanColumns
    print 'colToDrop list length with repeats: %s' % len(colToDrop)
    print 'colToDrop list length without repeats: %s' % len(colToDrop2)
    print 'colToDrop list length without repeats and with nan columns: %s' % len(colToDrop3)
    print colToDrop3
    
    print 'saving list of columns to drop and df after dropping...'
    with open('%s/featureSummaryDFs/colToDrop_filteredBy%s_corr%s' % (datasetFolder, sampleListName, filterFeaturesByCorr), 'wb') as fp:
        pickle.dump(colToDrop3, fp)
    
    # save clean featureDFs to pickle and excel
    TCRfeatureDF3.to_pickle(featureFilePickle2)
    TCRfeatureDF3.to_excel(featureFileExcel2)
    
    return TCRfeatureDF3


#-------------------------------------------------------------------------------------



#------------------------------------------------------------------
'''
the following function plots prediction results for a specific dataset (defined by ss and repeat)
it compares results from three different models (logReg, XGB, LGBM) and can compare results with and without TCR sequence
data.

input:
ss = int, number of sequence subsampled in the dataset
repeat- int, repeat number of the dataset
colorList - for plots
datacolumn = string 'roc_auc'/'pr_auc_corrected'
ax- axes info
withTCRdf = true/false - should include TCR seq or not
title - string
ylim, xlim - tuples

usage example:
colorlist=['blue','deepskyblue','darkgreen','lime','red','coral']
colList=[('roc_auc',(0.5,1)),('pr_auc_corrected',(0,0.6))]
withTCRdfList=[True,False/'only' (means use only the combination including TCRdf)
count=1

for dataColumnInfo in colList:
    dataColumn=dataColumnInfo[0]
    ylim=dataColumnInfo[1]
    for withTCRdf in withTCRdfList:
        fig,axes=plt.subplots(nrows=1, ncols=3,figsize=(12,7),sharex=True)
        datasetList=[(9000,2),(12500,1),(12500,2)]
        print count
        count=count+1
        for n,dataset in enumerate(datasetList):
            ax1=axes[n]
            ss=dataset[0]
            repeat=dataset[1]
            title='ss%s_rep%s' %(ss,repeat)
            ax=plot_prediction_results(ss,repeat,colorList,dataColumn,ax1,withTCRdf,xlim=(0,300),ylim=ylim,title=title)
            figName='predPerformance_compareDatasetFeaturenumFeaturetype_%s_withTCRdf%s_%s' %(dataColumn,withTCRdf,cdate)
            folder='%s/TCR_real_data/Predictions/isCardio_predictions' %MyPath
            figFile='%s/%s' %(folder,figName)
            fig.savefig(figFile,dpi=200)
            
##if not working, check the compatibility between XName and resultTypes definition####


'''
def helperPlotFunction(r0, r1, results, dataColumn, phenotypeOnly, colorList, modelName, ax, i, j=None):
    df = results[results['Xname'].str.contains(r1)]
    df['k'] = df['Num features used per split']
    df2 = df[['k', dataColumn]].sort_values(by=['k', dataColumn])
    df2 = df2.drop_duplicates(subset='k', keep='last')
    df3 = df2.copy()
    df3 = pd.concat([phenotypeOnly.T, df2], ignore_index=True)
#         print df3
    color = colorList[i * 2 + j]
    ax.plot(df3['k'], df3[dataColumn], label='-'.join([modelName, r0]), c=color, lw=2)
    ax.scatter(df3['k'], df3[dataColumn], c=color, label='')
    
    return ax


def plot_prediction_results(ss, repeat, colorList, dataColumn, ax, resultFolder=None, withTCRdf=True, title='', ylim=(0, 1), xlim=(0, 300),
                            filedate='14082018', YName='isCardio', datasetName=None, filterColumn=None, filterCrit=None):

    modelNameList = ['LogisticRegression', 'LGBMClassifier', 'XGBClassifier', ]

    for i, modelName in enumerate(modelNameList):
        print modelName
        
        # get result df from file:
        if resultFolder is not None:
            f1 = '%s/%s_predictions/%s/resultSummary_%s.xlsx' % (resultFolder, YName, modelName, filedate)
        else:
            if datasetName == 'Cardio126':
                f1 = '%s/TCR_real_data/CardioSamples/Predictions/%s_predictions/%s/resultSummary_%s.xlsx' % (MyPath, YName, modelName, filedate)
            else:
                f1 = '%s/TCR_real_data/Predictions/%s_predictions/%s/resultSummary_%s.xlsx' % (MyPath, YName, modelName, filedate)
        print 'result file is:'
        print f1
        results = pd.read_excel(f1)

        # get results for phenotypes only:
        if ss is None:
            phenotypeOnly = results[(results['Xname'].str.contains('phenotypesOnly')) & (~results['Xname'].str.contains('ss'))]    
        else:
            phenotypeOnly = results[results['Xname'].str.contains('MatchedSamples_ss%srep%s_phenotypesOnly' % (ss, repeat))]
        phenotypeOnly['k'] = 0
        phenotypeOnly = phenotypeOnly[['k', dataColumn]].sort_values(by=dataColumn)
        phenotypeOnly = pd.DataFrame(phenotypeOnly.iloc[-1])
    #     print phenotypeOnly
    #     print phenotypeOnly.T

        
        # get subset of results: (done only after phenotypeOnly extraction)
        if (filterColumn is not None) & (filterCrit is not None):
            print filterColumn
            print filterCrit
            results[filterColumn] = results[filterColumn].fillna('No')

            results = results[results[filterColumn].str.contains(filterCrit)]
            print results[filterColumn].value_counts(dropna=False)
        
        # get results for phenotype+PCA+feature and phenotype+PCA+tcrDF+features
        # process results dfs: add 'k' column, make sure there is only one row per each 'k', and add 'phenotype only' data as 
        # 'k'=o
        if ss is None:
                resultTypes = [('phenFeatPCA', '%s_phenotypes&PCA&TCRfeatures' % datasetName),
                     ('phenFeatPCAseqs', '%s_phenotypes&TCRdf&PCA&TCRfeatures' % datasetName)]
                
        else:
                resultTypes = [('phenFeatPCA', 'MatchedSamples_ss%srep%s_phenotypes&PCA&TCRfeatures' % (ss, repeat)),
                              ('phenFeatPCAseqs', 'MatchedSamples_ss%srep%s_phenotypes&TCRdf&PCA&TCRfeatures' % (ss, repeat))]
        if withTCRdf == 'both':
            print 'true'
            resultTypes = resultTypes
            
        elif withTCRdf == 'onlyYES':
            print 'only'
            resultTypes = resultTypes[1]
            j = 1  
        elif withTCRdf == 'onlyNO':
            print 'false'
            resultTypes = resultTypes[0]
            j = 0
            print resultTypes
 
        if type(resultTypes) == list:
            for j, Rtype in enumerate(resultTypes):
                r0 = Rtype[0]
                r1 = Rtype[1]
                print j, Rtype[1]
                ax = helperPlotFunction(r0, r1, results, dataColumn, phenotypeOnly, colorList, modelName, ax, i, j)
        else:
            r0 = resultTypes[0]
            r1 = resultTypes[1]
            ax = helperPlotFunction(r0, r1, results, dataColumn, phenotypeOnly, colorList, modelName, ax, i, j)

    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.margins()
    ax.legend() 
    ax.set_xlabel ('Selected feature number')
    ax.set_ylabel(dataColumn)
    ax.set_title(title)

    return ax



#------------------------------------------------------------------

###############################
#   new prediction functions: #
###############################

#-------------------------------------------------------------------


def add_text_at_corner(myplt, text, where='top right', **kwargs):
    legal_pos = ['top right', 'top left', 'bottom right', 'bottom left']
    if where not in legal_pos:
        print "where should be one of: " + ', '.join(legal_pos)
        return
    topbottom = where.split()[0]
    rightleft = where.split()[1]
    if str(type(myplt)) == "<class 'matplotlib.axes._subplots.AxesSubplot'>" or str(type(myplt)) == "<class 'mpl_toolkits.axes_grid1.parasite_axes.AxesHostAxes'>":
        x = myplt.get_xlim()
        y = myplt.get_ylim()
    elif str(type(myplt)) == "<type 'module'>":
        x = myplt.xlim()
        y = myplt.ylim()
    else:
       
        raise
    newaxis = {'left':x[0] + (x[1] - x[0]) * 0.01, 'right':x[1] - (x[1] - x[0]) * 0.01, 'top':y[1] - (y[1] - y[0]) * 0.01, 'bottom':y[0] + (y[1] - y[0]) * 0.01}
    myplt.text(newaxis[rightleft], newaxis[topbottom], text, horizontalalignment=rightleft, verticalalignment=topbottom, **kwargs)




# y=target series (binary), y_pred_df: df summarizing predicted y's (0<y<1)
def plot_ROC_PR_AUC(y, y_pred_df, ax,color1='mediumorchid',color2='navy',ticklabelsize=11,textsize=14,
                   labelsize=12,add_texts=True):

    from mpl_toolkits.axes_grid.inset_locator import inset_axes
    parent_axes = ax
    prevalence=round((y.sum() / y.shape[0]),2)
    
    fpr, tpr, thresholds = metrics.roc_curve(y.astype(float), y_pred_df.loc[y.index, 'pred_proba'])
    roc_auc = metrics.auc(fpr, tpr)
    
    parent_axes.plot(fpr, tpr, color=color1, lw=2, label='ROC curve (area = %0.3f)' % roc_auc)
    parent_axes.plot([0, 1], [0, 1], color=color2, lw=2, linestyle='--')
    parent_axes.set_xlim([-0.01, 1.0])
    parent_axes.set_ylim([0.0, 1.05])
    # parent_axes.set_xticks([0, 0.2, 0.4])
    parent_axes.set_xlabel('False Positive Rate', fontsize=labelsize)
    parent_axes.set_ylabel('True Positive Rate', fontsize=labelsize)
    parent_axes.tick_params(labelsize=ticklabelsize)
    
#     parent_axes.set_title('ACS (TCR features)', fontsize=20)

    inset_axes = inset_axes(parent_axes,
                        width="35%",  # width = 30% of parent_bbox
                        height="35%",  # height : 1 inch
    #                 bbox_to_anchor=(.0, 0., .5, .8), bbox_transform=parent_axes.transAxes)
                        loc=4, borderpad=1.5)
    inset_axes.set_xticks([0, 1])
    inset_axes.set_yticks([0, 1])
    for _, spine in inset_axes.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(1.)
        spine.set_color('black')

    precision, recall, _ = metrics.precision_recall_curve(y.astype(float), y_pred_df.loc[y.index, 'pred_proba'])
    pr_auc = metrics.auc(recall, precision)

    inset_axes.step(recall, precision, color='b', alpha=0.2, where='post')
    inset_axes.fill_between(recall, precision, step='post', alpha=0.5,
                     color=color1, label='Precision Recall curve - AUC = {0:0.3f}'.format(pr_auc))
    inset_axes.plot([0, 1], [y.sum() / y.shape[0], y.sum() / y.shape[0]], color=color2, lw=1, linestyle='--')
    inset_axes.set_xlabel('Recall', fontsize=9, labelpad=-15)
    inset_axes.set_ylabel('Precision', fontsize=9, labelpad=-12)
    inset_axes.set_ylim([0.0, 1.02])
    inset_axes.set_xlim([0.0, 1.0])
    
    if add_texts:
        add_text_at_corner(parent_axes, '\n AUC=%0.3g' % roc_auc, 'top left', fontsize=textsize, weight='bold')
        add_text_at_corner(inset_axes, 'PR=%0.2g' % pr_auc, 'bottom left', fontsize=textsize-2, weight='bold')
        inset_axes.text(0.02,(y.sum() / y.shape[0])+0.02,'prevalence=%s' %prevalence,fontsize=textsize-5)
#     plt.show()
    
    return ax, inset_axes,roc_auc, pr_auc,prevalence


'''
FOR REGRESSION MODELS - CHOOSE TASK='reg', for classification use task='class'
choose the right featureSelectionMethod for regression or classification if using pre-selection method. 


'''
# do not use featureSelectionMethod

def run_old_classification(X, XName, Y, k_list, model, modelName, model_params, task, resultFolder, featureSelectionMethod_list,
                           n_splits=10, AUC_figsize=(4, 3), nTOPfeatures=20, fig2toAnnotate=True):

    print ('klist:', k_list)
    print ('featureSelectionMethod_list:', featureSelectionMethod_list)
    
    # DEFINE NAMES AND FOLDERS:
    import datetime
    YName = pd.DataFrame(Y).columns[0]
    resultFolder = '%s/%s' % (resultFolder, YName)
    figResultFolder = '%s/plots' % resultFolder
    if not isdir(figResultFolder):
        makedirs(figResultFolder)
    dfResultFolder = '%s/resultDFs' % resultFolder
    if not isdir(dfResultFolder):
        makedirs(dfResultFolder) 
               
    modelParamString = ''   
    for key, value in model_params.iteritems():
        item2 = ':'.join([key, str(value)])
        modelParamString = ','.join([modelParamString, item2])      
    
    # GENERATE CV AND MODEL OBJECTS:
    if task == 'class':
        group_kfold = StratifiedKFold(n_splits=n_splits)
        groups = None
    else:
        group_kfold = GroupKFold(n_splits=n_splits)
        groups = np.array(range(X.shape[0]))
    
    print ('klist:', k_list)
    print ('featureSelectionMethod_list:', featureSelectionMethod_list)     
    for featureSelectionMethod in featureSelectionMethod_list:
        print featureSelectionMethod
        
        if (featureSelectionMethod == 'usingModel') & (model == LogisticRegression):
            print 'cant perform logistic regression with importance-based feature selection'
            break
            
        if featureSelectionMethod == None: 
            k_list = [10]
        for k in k_list:
            if task == 'class':
                columns = 'pred_proba'
            else:
                columns = 'pred_value'
            y_pred_df = pd.DataFrame(index=Y.index, columns=[columns])
            i = 0
            selected_features_list = []
            importancesDF_list = []
            
            # LOOP OVER CV FOLDS:
            for train_index, val_index in group_kfold.split(X, Y, groups):
                currentDT = datetime.datetime.now()
                forCurrentDF = currentDT.strftime("%d%m%y_%H%M%S")
                
                print i
                i += 1
                X_train, X_val = X.iloc[train_index], X.iloc[val_index]
                y_train, y_val = Y.loc[X_train.index], Y.loc[X_val.index]
                
                if task == 'class':
                    print 'fraction of 1s in train set=%s' % (float(y_train.sum()) / len(y_train))
                    print 'fraction of 1s in test set=%s' % (float(y_val.sum()) / len(y_val))

                if (featureSelectionMethod is not None) & (featureSelectionMethod != 'usingModel'):
                    if k is None:
                        k = 50  
                    if len(X_train.columns) < k:
                        featureSelectionMethod = None
                        print 'k is bigger than the number of columns, therefore no column was filtered by univariate'
                        selected_features = X_train.columns.values.tolist()
                        print selected_features
                    else:
                        print 'filtering by %s with k=%s' % (featureSelectionMethod, k)
                        if featureSelectionMethod == 'random':
                            import random
                            idxs_selected = random.sample(X_train.columns.values.tolist(), k)
                            X_train_new = X_train.loc[:, idxs_selected]
                            X_val_new = X_val.loc[:, idxs_selected]
                        else:
                            selector = SelectKBest(featureSelectionMethod, k)
                            selector.fit(X_train, y_train)
                            idxs_selected = selector.get_support(indices=True)
                            X_train_new = X_train.iloc[:, idxs_selected]
                            X_val_new = X_val.iloc[:, idxs_selected]
                        X_train = X_train_new
                        X_val = X_val_new
                        print 'k=%s' % k
                        print 'length of idxs_selected=%s' % len(idxs_selected)
                        print 'X_train shape after univaraite/random feature selection is  %s_%s' % (X_train.shape[0], X_train.shape[1])
                        print 'X_val shape after univaraite/random feature selection is  %s_%s' % (X_val.shape[0], X_val.shape[1])
                        print 'selected columns are:'
                        selected_features = X_train.columns.values.tolist()
                        print selected_features


                # creating the model object
                m = model(**model_params)
                m.fit(X_train, y_train)
                
                # get feature importances information from this fold:
                if model != LogisticRegression:
                    features = X_train.columns.tolist()
                    importances = m.feature_importances_
                    importancesDF = pd.DataFrame(index=features, data={i:importances}).T
                    importancesDF_list.append(importancesDF)

                if featureSelectionMethod == 'usingModel':  
                    # take only important fefatures:
                    indices = (np.argsort(importances)[::-1])[:k]
                    X_train_new = X_train.iloc[:, indices]
                    X_val_new = X_val.iloc[:, indices]

                    selected_features = X_train_new.columns.tolist()
                    print 'selected features are:'
                    print selected_features

                    # fit model again, now using only selected features:
                    m2 = model(**model_params)

                    # fitting the training
                    m2.fit(X_train_new, y_train)

                    print 'X_train shape after feature selection using model is  %s_%s' % (X_train_new.shape[0], X_train_new.shape[1])
                    print 'X_val shape after feature selection using model is  %s_%s' % (X_val_new.shape[0], X_val_new.shape[1])


                    # getting the predictions for the test using m2 (trained on select feature only!), X_train_new
                    # and X_val_new
                    if task == 'class':
                        y_pred_proba = m2.predict_proba(X_val_new)
                        y_pred_df.loc[y_val.index, :] = np.expand_dims(y_pred_proba[:, 1], 1)
                        
                    else:
                        y_pred = m2.predict(X_val_new)
                        y_pred_df.loc[y_val.index, :] = np.expand_dims(y_pred, 1)
                    numFeaturesPerSplit = len(X_train_new.columns)

                else:
                    # getting the predictions for the test using m1 (trained on all features or features pre-selected
                    # before applying the model:
                    if task == 'class':
                        y_pred_proba = m.predict_proba(X_val)
                        y_pred_df.loc[y_val.index, :] = np.expand_dims(y_pred_proba[:, 1], 1)
                    else:
                        y_pred = m.predict(X_val)
                        y_pred_df.loc[y_val.index, :] = np.expand_dims(y_pred, 1)
                    numFeaturesPerSplit = len(X_train.columns)

                # appending list of selected features:
                try:
                    selected_features_list = selected_features_list + selected_features
                except:
                    print 'no features were selected'
                    
            print 'done with CV...'
            
            # # PLOTTING RESULTS:
                    
            if type(featureSelectionMethod) == str:
                featureSelectionMethodString = featureSelectionMethod
            else:
                try:
                    featureSelectionMethodString = str(featureSelectionMethod).split(' ')[1].replace('_', '')
                except:
                    featureSelectionMethodString = str(featureSelectionMethod).split(' ')[0].replace('_', '')
            
                       
            figname = '%s_%s_%s_k%s_featSelMeth%s_%s' % (XName, YName, modelName, k,
                                                                      featureSelectionMethodString, modelParamString)
            figname = figname.replace(':', '')
            figname = figname.replace(',', '')
            figname = figname.replace('.', '')
            
            print ('********')
            print ('YName:', YName, 'model:', model, 'k:', k, 'featureSelectionMethod', featureSelectionMethodString)
            print ('model params', modelParamString)
            print ('********')
       
            # generate summarizing feature importance df and plot:
            try:
                sum_importances = pd.concat(importancesDF_list)
                print ('sum_importances_shape:', sum_importances.shape)
                sum_importances = sum_importances.T
                sum_importances['mean_importance'] = sum_importances.mean(axis=1)
                sum_importances['mean_importance'] = sum_importances['mean_importance'].round(4)
                sum_importances = sum_importances.sort_values(by='mean_importance', ascending=False)
    #             print 'sum importances head:'
    #             print sum_importances['mean_importance'].head()
                featureImportanceList = zip(sum_importances.index[:nTOPfeatures],
                                          sum_importances['mean_importance'][:nTOPfeatures])
                top20features = ''
                for x in featureImportanceList:
                    item = ':'.join([x[0], str(x[1])])
                    top20features = ','.join([top20features, item])
            
                # plot feature importances:
                fig1, ax1 = plt.subplots(figsize=(nTOPfeatures * 0.8, 5))
                sum_importances['mean_importance'].head(nTOPfeatures).plot(kind='bar', ax=ax1, color='mediumorchid', edgecolor='black', lw=2)
#                 ax1.set_title('Top %s important features over all CV-folds' %nTOPfeatures)
                ax1.tick_params(axis='x', labelsize='x-large')
                ax1.tick_params(axis='y', labelsize='large')
                fig1.savefig('%s/featImp_%s' % (figResultFolder, figname), bbox_inches='tight')
                plt.show()
            except:
                print 'no feature importance info...'
                top20features = ''
                fig1 = None
                ax1 = None
            fig1tuple = tuple([fig1, ax1])
    
            
            if task == 'class':
                # calculate and plot ROC+PR graphs:
                fig2, ax2 = plt.subplots(figsize=AUC_figsize)
                y = Y
                ax2, roc_auc, pr_auc = plot_ROC_PR_AUC(y, y_pred_df, ax2)
                plt.show()
                fig2.savefig('%s/ROC_PR_plot_%s' % (figResultFolder, figname), bbox_inches='tight')
                text2 = None
                r_pear = None
                p_pear = None
            else:
                fig2, ax2 = plt.subplots(figsize=AUC_figsize)
                data1 = Y.tolist()
                data2 = y_pred_df['pred_value'].tolist()
                data1name = 'Y values'
                data2name = 'Predicted Y values'
                ax = ax2
                title = "%s prediction by TCR features\n(Pearson's correlation test" % YName
                corrType = 'pearson'
                
                ax2, nsamples, r_pear, p_pear, text2 = plot_corr(data1, data2, data1name, data2name, ax, title, corrType,
                                                      toAnnotate=fig2toAnnotate, color='darkblue', s=100, alpha=0.4)
                r_spear, p_spear = MySpearmanr(data1, data2)
                fig2.savefig('%s/corr_PR_plot_%s' % (figResultFolder, figname), bbox_inches='tight')
                plt.show()
            
            fig2tuple = tuple([fig2, ax2, text2, r_pear, p_pear])
                
            
            # GENERATE SUMMARIZING DF:
            Resdf = pd.DataFrame()
            Resdf.loc[0, 'XName'] = XName
            Resdf.loc[0, 'YName'] = YName
            Resdf.loc[0, 'modelName'] = modelName
            Resdf.loc[0, 'k'] = k
            Resdf.loc[0, 'featureSelectionMethod'] = featureSelectionMethod
            Resdf.loc[0, 'top%sfeatures' % nTOPfeatures] = top20features
            Resdf.loc[0, 'forCurrentDF'] = forCurrentDF 
            Resdf.loc[0, 'modelParams'] = modelParamString
            Resdf.loc[0, 'task'] = task
            
            if task == 'class':
                Resdf.loc[0, 'roc_auc'] = round(roc_auc, 3)
                Resdf.loc[0, 'pr_auc'] = round(pr_auc, 3)
                Resdf.loc[0, 'perc_pos_target'] = round(float(Y.sum()) / len(Y), 3)
                Resdf.loc[0, 'pr_auc_corrected'] = round(pr_auc, 3) - round(float(Y.sum()) / len(Y), 3)
            else:
                Resdf.loc[0, 'r_pear'] = round(r_pear, 2)
                Resdf.loc[0, 'p_pear'] = round(p_pear, 4)
                Resdf.loc[0, 'r_spear'] = round(r_spear, 2)
                Resdf.loc[0, 'p_spear'] = round(p_spear, 4) - round(float(Y.sum()) / len(Y), 3)

            Resdf.to_excel('%s/Resdf_%s.xlsx' % (dfResultFolder, figname))
           
            print 'Prediction round is over'
    
    print 'all predictions are done'       
    return Resdf, fig1tuple, fig2tuple



#-----------------------------------------------------------------------
def plot_roc_curve(ax,y_true,y_pred,pos_label,dataset_name='roc curve',color='orange'):
    from sklearn import metrics
    fpr, tpr, thresholds = metrics.roc_curve(y_true,y_pred,pos_label)
    roc_auc = metrics.auc(fpr, tpr)
    print ('ROC_AUC: ',roc_auc)
    
    lw=2
    ax.plot(fpr, tpr, color=color, lw=lw, label='%s (area = %0.3f)' %(dataset_name,roc_auc))
    ax.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False Positive Rate', fontsize='large')
    ax.set_ylabel('True Positive Rate', fontsize='large')
    ax.set_title('ROC curve - ' , fontsize='x-large')
    ax.legend(loc="lower right", fontsize='large')
    
    return ax, roc_auc

def plot_PR_curve(ax,y_true,y_pred_proba,pos_label,dataset_name='Precision-Recall curve',color='b',alpha=0):

    from sklearn import metrics
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    for _, spine in ax.spines.items():
        spine.set_visible(True)
        spine.set_linewidth(1.)
        spine.set_color('black')

    precision, recall, _ = metrics.precision_recall_curve(y_true.astype(float), y_pred_proba.loc[y_true.index, :].iloc[:,0])
    pr_auc = metrics.auc(recall, precision)
    prevalence=round((y_true.sum() / y_true.shape[0]),2)
    ax.step(recall, precision, color=color, where='post',label='%s({0:0.3f})'.format(pr_auc) %dataset_name)
    ax.fill_between(recall, precision, step='post', alpha=alpha,
                     color=color)
    ax.plot([0, 1], [y_true.sum() / y_true.shape[0], y_true.sum() / y_true.shape[0]], color='navy', lw=1, linestyle='--')
    ax.set_xlabel('Recall', fontsize='large', labelpad=-15)
    ax.set_ylabel('Precision', fontsize='large', labelpad=-12)
    ax.set_ylim([0.0, 1.02])
    ax.set_xlim([0.0, 1.0])

    
    return ax,pr_auc
#---------------------------------------------------------

def plot_many_roc(ax,sample_list,color_map,prediction_subdir_list=None,prediction_dir=None,
                 dataset_name_list=None,AUC_diff_min=None,AUC_diff_max=None):

    Y_file='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/TargetDFs/isCardio.dat'
    Y=pd.read_pickle(Y_file).loc[sample_list,:]
    # print Y.head()

    #plot best prediction:
    if prediction_dir is None:
        prediction_dir='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/isCardio/'
    if prediction_subdir_list is None:
        prediction_subdir_list=['XGB_randomSearch_25_byPredictedGender/', 'XGB_randomSearch_25_byPredictedAge/',
               'XGB_randomSearch_25_byPredictedAgeGender/','XGB_randomSearch_25_byRepFeatPCA10percVDJ0999PredictedAgeGender/']
    cmap=plt.get_cmap(color_map, len(prediction_subdir_list)+1) #get seperate colors for graphs
    color_list=[cmap(x+1) for x in range(len(prediction_subdir_list))]
    [cmap(1),cmap(2),cmap(3),cmap(4)]
    # print colorlist
    if dataset_name_list is None:
        dataset_name_list=['Pred. Gender','Pred. Age','Pred. Age+Gender','TCR features + Pred. Age+Gender']
    
    # for calculation of AUC difference between the 'interesting' curve ('max') and the 'control' curve ('min'),
    #define which dataset is the min and the max:
    if AUC_diff_min is None:
        AUC_diff_min='byPredictedAgeGender'
    if AUC_diff_max is None:
        AUC_diff_max='byRepFeat'
        
    ### plot ROC curves for each prediction data set in the prediction_subdir_list:
    for n in range(len(prediction_subdir_list)):
        d=prediction_subdir_list[n]
        color=color_list[n]

        pred=pd.read_pickle(prediction_dir+d+'predictions_df.pkl') #get prediction data
        y_true=Y.isCardio
        y_pred=pred.isCardio.loc[sample_list]
        pos_label=1
        dataset_name=dataset_name_list[n]

        #plot:
        ax, roc_auc=plot_roc_curve(ax,y_true,y_pred,pos_label,dataset_name=dataset_name,color=color)
        
        if AUC_diff_min in d: auc_min=roc_auc
        if AUC_diff_max in d: auc_max=roc_auc
    
    #calculate AUC diff:
    AUC_diff=round(auc_max-auc_min,3)
    
    #edit legend:
    handles, labels = ax.get_legend_handles_labels()
    labels=[l.replace('area = ','') for l in labels]
    ax.legend(handles, labels,loc='lower right')

    return ax,AUC_diff

def plot_many_PR(ax,sample_list,color_map,prediction_subdir_list=None,prediction_dir=None,
                 dataset_name_list=None,AUC_diff_min=None,AUC_diff_max=None):
    
    Y_file='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/TargetDFs/isCardio.dat'
    Y=pd.read_pickle(Y_file).loc[sample_list,:]
    if AUC_diff_min is None:
        AUC_diff_min='byPredictedAgeGender'
    if AUC_diff_max is None:
        AUC_diff_max='byRepFeat'

    ###get default data if necessary:
    if prediction_dir is None:
        prediction_dir='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/isCardio/'
    if prediction_subdir_list is None:
        prediction_subdir_list=['XGB_randomSearch_25_byPredictedGender/', 'XGB_randomSearch_25_byPredictedAge/',
               'XGB_randomSearch_25_byPredictedAgeGender/','XGB_randomSearch_25_byRepFeatPCA10percVDJ0999PredictedAgeGender/'] 
    if dataset_name_list is None:
        dataset_name_list=['Pred. Gender','Pred. Age','Pred. Age+Gender','TCR features + Pred. Age+Gender']
        
    ### define graph colors:   
    cmap=plt.get_cmap(color_map, len(prediction_subdir_list)+1) #get seperate colors for graphs
    color_list=[cmap(x+1) for x in range(len(prediction_subdir_list))]

    ### plot PR graph for each data set:
    for n in range(len(prediction_subdir_list)):
        d=prediction_subdir_list[n]
        color=color_list[n]
        pred=pd.read_pickle(prediction_dir+d+'predictions_df.pkl')
        y_true=Y.isCardio
        y_pred=pred.isCardio.loc[sample_list]
        pos_label=1
        dataset_name=dataset_name_list[n]

        ax,pr_auc=plot_PR_curve(ax,y_true,pred,pos_label,dataset_name=dataset_name,color=color) #don't input alpha to have alpha=0 and not filling.
        
        if AUC_diff_min in d: auc_min=pr_auc
        if AUC_diff_max in d: auc_max=pr_auc
    
    #calculate AUC diff:
    AUC_diff=round(auc_max-auc_min,3)
    #edit legend:
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(loc='best',fontsize='medium')   

    return ax,AUC_diff


def plot_many_corrs(ax,YName,Y_file,sample_list,color_map,prediction_subdir_list=None,prediction_dir=None,
                 dataset_name_list=None):

    Y=pd.read_excel(Y_file).set_index('BD').loc[sample_list,:]
    # print Y.head()

    #plot best prediction:
    if prediction_dir is None:
        prediction_dir='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/Age/'
    if prediction_subdir_list is None:
        prediction_subdir_list=['XGB_randomSearch_25_byOldTCRfeatureMatrix','XGB_randomSearch_25_byOldTCRfeatureMatrix_ageGenderRegCorr','XGB_randomSearch_25_byNewTCRfeatureMatrix_ageRegCorr',
         'XGB_randomSearch_25_byOldVDJ','XGB_randomSearch_25_byOldVDJ_ageGenderRegCorr']
    cmap=plt.get_cmap(color_map, len(prediction_subdir_list)+1) #get seperate colors for graphs
    color_list=[cmap(x+1) for x in range(len(prediction_subdir_list))]
    [cmap(1),cmap(2),cmap(3),cmap(4)]
    # print colorlist
    if dataset_name_list is None:
        dataset_name_list=['Pred. Gender','Pred. Age','Pred. Age+Gender','TCR features + Pred. Age+Gender']

        
    ### plot ROC curves for each prediction data set in the prediction_subdir_list:
    r_list=[]; p_list=[]
    handles_list=[]; labels_list=[]
    for n in range(len(prediction_subdir_list)):
        d=prediction_subdir_list[n]
        color=color_list[n]

        pred=pd.read_pickle(prediction_dir+d+'predictions_df.pkl') #get prediction data
        y_true=pd.DataFrame(Y[YName].rename('y_true'))
        y_pred=pd.DataFrame(pred[YName].loc[sample_list].rename('y_pred'))
        merged=pd.merge(y_true,y_pred,how='inner',left_index=True,right_index=True)
        dataset_name=dataset_name_list[n]
        data1name='True Value (%s)' %YName
        data2name='Predicted Value (%s)' %YName
        

        #plot:
        ax, nsamples,r,p,text,handles, labels=plot_corr(merged['y_true'],merged['y_pred'],data1name,data2name,ax,title=None,corrType='pearson',
              toAnnotate=False,plotTrendLine=True,scatter_kws={'alpha':0.4,'s':100,'color':color}, text_kws={'fontsize':'x-large'},
                                       linecolor=color)
        r_list.append(r)
        p_list.append(p)
        handles_list.append(handles)
        labels_list.append(labels)


    return ax,r_list,p_list,handles_list,labels_list


