###################################################################################################
# File: FeatureCorrection_byLinearRegression.py
# Version: 0.0
# Date: 11.11.18
# Shani Ben-Ari Fuchs, shani.ben-ari@weizmann.ac.il
# 
# Python version: 2.7
###################################################################################################
#
# add:
# 1.


import os
import cPickle as pickle
from scipy.stats.stats import spearmanr, pearsonr
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
import json
from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import * 
from sklearn.linear_model import LogisticRegression, LinearRegression


def check_features_to_correct(X, Xname, phenDF, phen, sampleSubset, command_args):
    
        if sampleSubset is not None:
            X = X.loc[sampleSubset, :]
            phenotypeDF = pd.DataFrame(phenDF.loc[sampleSubset, phen])        
            print phenotypeDF.head()
#               
        print ('%s series type:' % phen, phenDF[phen].dtype)
        
        if (phenotypeDF[phen].nunique() > 2) or (command_args.useCorrForBinary):
            f1 = '%s/%s_%s_%s_%s.xlsx' % (command_args.path_to_association_doc_dir, Xname, phen, command_args.corrTestType,str(command_args.FDR).replace('.',''))       
        elif phenotypeDF[phen].nunique() == 2:
            f1 = '%s/%s_%s_ttest_%s.xlsx' % (command_args.path_to_association_doc_dir, Xname, phen,str(command_args.FDR).replace('.',''))
        f1=f1.replace(' ','')
        if command_args.calc_association:
            print 'checking associations between %s and all TCR features:' % phen
            
            if (phenotypeDF[phen].nunique() > 2) or (command_args.useCorrForBinary):
                print 'number of phenotype unique values is %s and useCorrForBinary=%s\
therefore correlation test is used' % (phenotypeDF.nunique(),command_args.useCorrForBinary)
                feature_phen_assoc = calc_numPhens_TCRfeatures_associations_correlation(X, Xname, phenotypeDF, phen,
                                                    command_args.corrTestType,command_args.FDR)
            elif phenotypeDF[phen].nunique() == 2:
                print 'number of phenotype unique values is %s and therefore ttest is used' % phenotypeDF.nunique()
                feature_phen_assoc = calc_catPhens_TCRfeatures_associations_ttest(X, Xname, phenotypeDF, phen, command_args.FDR)
                                                                        
            feature_phen_assoc.to_excel(f1)
        else:
            feature_phen_assoc = pd.read_excel(f1)
            
        # CHECK WHICH FEATURES TO CORRECT:
        sigDF = feature_phen_assoc[feature_phen_assoc['sig. by FDR=%s' % command_args.FDR] == 1]
        featToCorrect = sigDF['Feature'].unique().tolist()
        print ('number of features to correct is:', len(featToCorrect))
        print ('the most significantly associated features are:', sigDF['Feature'].unique().tolist()[:5])
        
        return featToCorrect
    
    
def correct_feature_by_reg(X, Xname, phen, phenDF, featToCorrect, command_args):
    
    # correct all featues for phen and generate new X:
    # PNP530 and Cardio126 samples:
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP530', 'rb') as fp:
        PNP530 = pickle.load(fp)
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/Cardio126', 'rb') as fp:
        Cardio126 = pickle.load(fp)
    PNP530Cardio126 = PNP530 + Cardio126
    
    X_Corr = pd.DataFrame(index=PNP530Cardio126, columns=X.columns, data=np.nan)
    for n, feature in enumerate(X_Corr.columns):
        if n % 100 == 0:
            print n, feature
        if feature in featToCorrect:
            phen_data = pd.DataFrame(phenDF[phen].loc[PNP530Cardio126].fillna(phenDF[phen].median()))
            if 'norm' in feature:
                feature_data = pd.DataFrame(X[feature].loc[PNP530Cardio126].fillna(X[feature].median()))
            else:
                feature_data = pd.DataFrame(X[feature].loc[PNP530Cardio126].fillna(0)) 
            featureName = feature
            
            res = calc_residuals(feature_data, phen_data, featureName, phen, PNP530)  # use only PNP530 samples to study
            X_Corr = pd.merge(X_Corr, pd.DataFrame(res), how='inner', left_index=True, right_index=True)
            X_Corr = X_Corr.rename(columns={'res':feature + '_' + phen + 'CorrRes'})
            X_Corr = X_Corr.drop(feature, axis=1)
        else:
            X_Corr[feature] = X[feature]
    print ('X_Corr shape after correcting all features is', X_Corr.shape)
    
    Xname = Xname + '_' + phen + 'Corr'
    return X_Corr, Xname    

# regress feature using phenotype data and generate the resulting residuals:
# the function should take feature data and phenotype data for PNP530+Cardio126!!!!

def calc_residuals(feature_data, phen_data, featureName, phen, sampleSubset):
    lm = LinearRegression()  # intercept=True, normalize=False
    lm.fit(phen_data.loc[sampleSubset, :], feature_data.loc[sampleSubset, :])
#     print ('coef: ', coef_)
#     print ('intercept: ', intercept_)
#     print ('score: ', lm.score(phen_data.loc[sampleSubset, :], feature_data.loc[sampleSubset, :]))
    y_pred = pd.DataFrame(index=phen_data.index.tolist(), data=lm.predict(phen_data))
    
    merged = pd.merge(feature_data, y_pred, how='inner', left_index=True, right_index=True)
#     merged = pd.merge(merged, phen_data, how='inner', left_index=True, right_index=True)
    merged['res'] = merged[featureName] - merged[0]
    
#     f1 = '%s/TCR_real_data/TCR_phenotype_relations/resDFs/merged_%s_%s.xlsx' % (MyPath, phen, featureName)
#     merged.to_excel(f1)
#     print ('sum of residuals for the fitted data only:', merged.loc[sampleSubset,'res'].sum())
    
#     print merged.head(10)
#     print merged[featureName].median()
    return merged['res']




def main():
    print ('main')
    parser = argparse.ArgumentParser()
    parser.add_argument('path_to_XtoRegress', help='X matrix that needs correction', type=str, default=None)
    parser.add_argument('-path_to_phenotype_DF', help='path to phenotype DF', type=str, default=None)
    parser.add_argument('-phen_list', help='list of phenotypes to correct for, order matters, encapsulate with "',
    type=json.loads, default=None)
    parser.add_argument('-path_to_correctedX_dir', help='path_to_correctedX_dir', type=str,
    default='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Predictions/featureDFs/')
    parser.add_argument('-path_to_association_doc_dir', help='path_to_association_doc_dir', type=str,
    default='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/TCR_phenotype_relations/')
    parser.add_argument('-checkAssociationsOnPNP530Only', help='heckAssociationsOnPNP530Only', type=bool, default=True)
    parser.add_argument('-calc_association', help='calc_association - True/False', type=bool, default=True)
    parser.add_argument('-FDR', help='association test FDR', type=float, default=0.25)
    parser.add_argument('-corrTestType', help='corrTestType - pearson or spearman', type=str, default='pearson')
    parser.add_argument('-useCorrForBinary',help='True or False', type=bool, default=True)
    command_args = parser.parse_args()
    
    # check folder existance: 
    if (not os.path.exists(command_args.path_to_XtoRegress)) or (not os.path.exists(command_args.path_to_phenotype_DF)):
        print ("XtoRegress or phenotype_DF doesn't exist!"); return
    if not isdir(command_args.path_to_correctedX_dir):
        makedirs(command_args.path_to_correctedX_dir)
        print ('made new correctedX_dir')
    if not isdir(command_args.path_to_association_doc_dir):
        makedirs(command_args.path_to_association_doc_dir)
        print ('made new association_doc_dir')
        
    
    # definitions:
    Xfile = command_args.path_to_XtoRegress
    Xname = Xfile.split('/')[-1]
    Xdir = Xfile.replace(Xname, '')
    Xname = Xname.replace('.dat', '')
    print ('Xfile: ', Xfile)
    print ('Xdir: ', Xdir)
    print ('Xname: ', Xname)
    
    try:
        X = pd.read_pickle(Xfile)
    except:
        X = pd.read_excel(Xfile)
    try:
        X = X.set_index('BD')
    except:
        print 'X Index is already "bd" or no such column exist'
    
    phenotypeDFfile = command_args.path_to_phenotype_DF
    try:
        phenotypeDF = pd.read_pickle(phenotypeDFfile)
    except:
        phenotypeDF = pd.read_excel(phenotypeDFfile)
    try:
        phenotypeDF = phenotypeDF.set_index('BD')
    except:
        print 'phenotypeDF index is already "bd" or no such column exist'
        
    
    phen_list = command_args.phen_list
    print ('phenotypes to correct for are:', phen_list)
    
    if command_args.checkAssociationsOnPNP530Only:
        with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP530', 'rb') as fp:
            PNP530 = pickle.load(fp)
        sampleSubset = PNP530
    
    # loop over all phenotypes:
    tested_phen_list = []
    for phen in phen_list:
        tested_phen_list.append(phen)
        print ('correcting ' + Xname + ' for ' + phen)
        phenDF = pd.DataFrame(phenotypeDF[phen])
        print phenDF.head()
        
        # check which features to correct:
        featToCorrect = check_features_to_correct(X, Xname, phenDF, phen, sampleSubset, command_args)
        
        # correct features:
        X_Corr, Xname = correct_feature_by_reg(X, Xname, phen, phenDF, featToCorrect, command_args)
        f3 = Xdir + '/' + Xname + '%s.dat' %str(command_args.FDR).replace('.','')
        f3=f3.replace(' ','')
        X_Corr.to_pickle(f3)
        X = X_Corr
        
        # check if corrected matrix is still correlated with all phenotypes tested so far:
        for tested_phen in tested_phen_list:
            print 'checking which features are still associated with %s:' % tested_phen
            tested_phenDF = pd.DataFrame(phenotypeDF[tested_phen])
            stillcorrelated = check_features_to_correct(X_Corr, Xname, tested_phenDF, tested_phen, PNP530, command_args)

    # save parameters to a text file:
#     print ('saving parameters to file')
#     argFile=command_args.output_dir+'argFile.txt'
#     new_file=open(argFile,mode="w")
#     new_file.write("model: "+command_args.model)
#     new_file.write("predictor_params: "+str(command_args.predictor_params)+"\n")
#     new_file.write("-path_to_X: "+command_args.path_to_X+"\n")
#     new_file.write("path_to_Y: "+command_args.path_to_Y+"\n")
#     new_file.write("k_folds: "+str(command_args.k_folds)+"\n")
#     new_file.write("parameter_search_type: "+command_args.parameter_search_type+"\n")
#     new_file.write("n_random: "+str(command_args.n_random)+"\n")   
#     new_file.close()
    
if __name__ == "__main__":
    sethandlers()
    main()
