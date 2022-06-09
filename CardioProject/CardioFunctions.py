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
# from SufficientStatistics.SufficientStatistics import *
from ShaniBA.MyFunctionsShani import *
import math
# from myplots import roundup, rounddown, find_decimal_fold
from skbio.stats.distance import mantel
from scipy.spatial.distance import braycurtis, pdist
from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import *
from ShaniBA.TCR_feature_generation.publicSeqAnalysis import *
from ShaniBA.TCR_feature_generation.SubsamplingFunctions import *
from ShaniBA.TCR_microbiome_interactions.TCR_microbiome_interactions_functions import *
from ShaniBA.PhenotypicData.PhenotypeGenerationFunctions import *


import os
# from Utils import cacheOnDisk
from SegalQueue.qp import qp, fakeqp
from addloglevels import sethandlers

MyPath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF'

#----------------------------------------------------------------------------
'''
the following function is used to calculate fisher association between sequence presence/absence and cohorts

usage example:
datasetFolder1='%s/TCR_real_data' %MyPath
datasetName1='PNP530'
minNshared1=5
datasetFolder2='%s/TCR_real_data/CardioSamples' %MyPath
datasetName2='Cardio126'
minNshared2=2


### get relevant PNP sequences:
f1='%s/sharingAnalysis/sharingMatrix_%s_minNshared%s_RA_onlyProductiveTrue__percShared5_OLtrimmed_binary' %(datasetFolder1,
                                                                                          datasetName1,minNshared1)                 
tcrDF_binary1=pd.read_pickle(f1)

tcrDF_binary1['cohort']=datasetName1

print 'tcrDF_binary1 shape is %s_%s' %(tcrDF_binary1.shape[0],tcrDF_binary1.shape[1])

f2='%s/sharingAnalysis/sharingMatrix_%s_minNshared%s_RA_onlyProductiveTrue__percShared5_OLtrimmed_binary' %(datasetFolder2,
                                                                                          datasetName2,minNshared2)                 
tcrDF_binary2=pd.read_pickle(f2)
tcrDF_binary2['cohort']=datasetName2

print 'tcrDF_binary2 shape is %s_%s' %(tcrDF_binary2.shape[0],tcrDF_binary2.shape[1])

twoCohortsDF=pd.concat([tcrDF_binary1,tcrDF_binary2])
twoCohortsDF=twoCohortsDF.fillna(0)
print twoCohortsDF.shape

nPerm=99

minSeq=20
maxSeq=minSeq+20
df=pd.DataFrame()

for n, seq in enumerate(twoCohortsDF.columns.values[minSeq:maxSeq]):
    print minSeq+n
    df=calc_permFisher_twoCohorts(seq,twoCohortsDF,datasetName1,datasetName2,nPerm)
print 'done'   


'''



def calc_permFisher_twoCohorts(minSeq, maxSeq, twoCohortsDF, datasetName1, datasetName2, nPerm, folderToSave):
    df = pd.DataFrame()
    
    for n, seq in enumerate(twoCohortsDF.columns.values[minSeq:maxSeq]):
        print minSeq + n
 
    
    
        # calculate real fisher OR:
        seqDF = twoCohortsDF[[seq, 'cohort']]
        tab = pd.crosstab(seqDF[seq], seqDF.cohort)
    #     print tab
        
        abs1 = tab.loc[0, datasetName1]
        abs2 = tab.loc[0, datasetName2]
        pres1 = tab.loc[1, datasetName1]
        pres2 = tab.loc[1, datasetName2]
        
        OR, p = fisher_exact(tab, alternative='two-sided')
        print 'real abs(log10OR=%s)' % abs(np.log10(OR))
       
        moreSigCountFisher = 1
        totalCountFisher = 1
        
    
        for n in range(nPerm):
            shuffled = seqDF.apply(lambda x: x.sample(frac=1).values)
            tab_shuff = pd.crosstab(shuffled[seq], shuffled.cohort)
            OR_shuff, p_shuff = fisher_exact(tab_shuff, alternative='two-sided')    
            print n, abs(np.log10(OR_shuff))
    
            if abs(np.log10(OR_shuff)) >= abs(np.log10(OR)):
                            moreSigCountFisher = moreSigCountFisher + 1
            totalCountFisher = totalCountFisher + 1
    
        permP_fisher = float(moreSigCountFisher) / totalCountFisher
        print 'permP=%s' % permP_fisher
        
        df.loc[seq, 'realAbsLog10OR'] = abs(np.log10(OR))
        df.loc[seq, 'realP'] = p
        df.loc[seq, 'permP'] = permP_fisher
        df.loc[seq, 'nPerm'] = nPerm
        df.loc[seq, 'abs_%s' % datasetName1] = abs1
        df.loc[seq, 'pres_%s' % datasetName1] = pres1
        df.loc[seq, 'abs_%s' % datasetName2] = abs2
        df.loc[seq, 'pres_%s' % datasetName2] = pres2   
        
    
    f1 = '%s/permFisher_%s_%s_%s' % (folderToSave, nPerm, minSeq, maxSeq)
    df.to_pickle(f1)
    
    return df


#------------------------------------------------------------------
def logReg_withWrapper(X, Xname, yDF, yname, resultFolder, yMinCol=None, yMaxCol=None, maxiter=35):
    
    resultFileName = 'resultDFlogReg_%s_%s_ycol%s_ycol%s_mi%s.xlsx' % (Xname, yname, yMinCol, yMaxCol, maxiter)
    filenames = [f for f in listdir(resultFolder) if isfile(join(resultFolder, f))]
    if resultFileName not in filenames:
    
        import statsmodels.formula.api as sm
        resultsDF = pd.DataFrame()
        if yMinCol is None:
            yMinCol = 0
        if yMaxCol is None:
            yMaxCol = len(yDF.columns)
            
        xsamples = X.index.astype(str).tolist()  # get list of samples to use
    
        for n, target in enumerate(yDF.columns.values[yMinCol:yMaxCol]):
        #     if n<10:
                print n, target
                y = yDF[target]
    
                # (1) process y
    
        #         print 'original y shape is %s' %y.shape[0]
                y = y.dropna(how='any')
        #         print 'y shape after na removal=%s' %y.shape[0]
                y = y.loc[xsamples]
        #         print 'y shape after filtering for x samples=%s' %y.shape[0]
        #         print y.shape
    
                # (2) get seq pres/abse info and fisher test:
                cohortDF = pd.merge(pd.DataFrame(X.isCardio), pd.DataFrame(y), how='left', left_index=True, right_index=True)
        #         print cohortDF.head()
                tab = pd.crosstab(cohortDF[target], cohortDF.isCardio)
                abs0 = tab.loc[0, 0]
                abs1 = tab.loc[0, 1]
                pres0 = tab.loc[1, 0]
                pres1 = tab.loc[1, 1]
                OR, p = fisher_exact(tab, alternative='two-sided')
    
    
                # (3) run model:
    
                try:
                    model = sm.Logit(y, X)
                    result = model.fit(maxiter=maxiter)
    
                    # (4) summarize p-values
    
                    if n == 0:
                        resultDF = pd.DataFrame(result.pvalues).T
                    else:
                        resultDF = pd.concat([resultDF, pd.DataFrame(result.pvalues).T], ignore_index=True)
                    resultDF.loc[n, 'target'] = target
                    resultDF.loc[n, 'p_fisher'] = p
                    resultDF.loc[n, 'abs_PNP'] = abs0
                    resultDF.loc[n, 'pres_PNP'] = pres0
                    resultDF.loc[n, 'abs_Cardio'] = abs1
                    resultDF.loc[n, 'pres_Cardio'] = pres1
    
                except:
                    print ' couldnt fit model'
                    if n == 0:
                        resultDF = pd.DataFrame()
                    else:
                        resultDF.loc[n, 'target'] = target
                        resultDF.loc[n, 'p_fisher'] = p
                        resultDF.loc[n, 'abs_PNP'] = abs0
                        resultDF.loc[n, 'pres_PNP'] = pres0
                        resultDF.loc[n, 'abs_Cardio'] = abs1
                        resultDF.loc[n, 'pres_Cardio'] = pres1
                                                
    
                
                
                
    #     resultDF['isCardioBest']=np.where(resultDF.isCardio==resultDF.min(axis=1),1,0)
    #     resultDF['bestP']=resultDF.min(axis=1)
        
        resultLength = len(resultDF)
        nNans = len(resultDF[resultDF['isCardio'].isnull()])
        nSig = len(resultDF[resultDF['isCardio'] <= 0.05])
        percNan = 100 * float(nNans) / resultLength
        percSig = 100 * float(nSig) / resultLength
        print 'percent null results is %s' % round(percNan, 2)
        print 'percent sig results is %s' % round(percSig, 2)
        
                
        f1 = '%s/resultDFlogReg_%s_%s_ycol%s_ycol%s_mi%s.xlsx' % (resultFolder, Xname, yname, yMinCol, yMaxCol, maxiter)
        resultDF.to_excel(f1)
    else:
        print 'this result df already exist in folder'
        resultDF = pd.DataFrame()
    
    return resultDF    
    
#---------------------------------------------------
'''
the following function takes predictor variable matrix and target matrix and tries to fit a logistic regression or linear regression
method for the data. then it calculates the p-values of all coeeficients.
it adds fisher test/ks+ttest p_values

this function doesn't enable prediction!

input:
X-predicting variable matrix. index should be sample banes
Xname- string. 
yDF- a dataframe including one or more targets, each in a seperate column
yname-string
resultFolder-the folder to save the results
modelType - 'ols' for linear regression and 'logit' for logistic regression
yMinCol- the index of the first column of yDF to use. int, optional (used mainly for parallizing the function)
yMaxCol- the index of the next-to-last column of yDF to use. int, optional (used mainly for parallizing the function)
maxiter-int,optional. defines the maximal iteration for the model

usage example:
datasetFolder='%s/TCR_real_data/PNP530Cardio126Combined' %MyPath
datasetName='PNP530Cardio126'

f1='%s/TCR_real_data/PNP530Cardio126Combined/featureSummaryDFs/PNP530Cardio126_allFeatures' %MyPath
PNP530Cardio126_allFeatures=pd.read_pickle(f1)
print PNP530Cardio126_allFeatures.shape



f2='%s/TCR_real_data/PNP530Cardio126Combined/Phenotypes/PNP530Cardio126_phen.xlsx' %MyPath
PNP530Cardio126_phen=pd.read_excel(f2)
PNP530Cardio126_phen.head()

# Xtypes:
allNum=['BD', 'Age', 'BMI', 'Creatinine', 'isCardio', 'nTemplates', 'Female', 'Smoking_Past', 'Smoking_Yes',
       'PCR_Plate1', 'PCR_Plate10', 'PCR_Plate2', 'PCR_Plate3',
       'PCR_Plate4', 'PCR_Plate5', 'PCR_Plate6', 'PCR_Plate7',
       'PCR_Plate8', 'PCR_Plate9']
small=['BD', 'Age', 'BMI', 'Creatinine', 'isCardio', 'Female']
smallNoCardio=['BD', 'Age', 'BMI', 'Creatinine', 'Female']
allNumNoPlate=['BD', 'Age', 'BMI', 'Creatinine', 'isCardio', 'nTemplates', 'Female', 'Smoking_Past', 'Smoking_Yes']

#DEFINE WHICH FEATURES TO USE:

Xcols=allNumNoPlate
Xname='allNumNoPlate'

X=PNP530Cardio126_phen[Xcols].set_index('BD')
#drop all samples that don't have all feature info. check how many samples were lost
print 'X shape before na removal=%s_%s' %(X.shape[0],X.shape[1])
X=X.dropna(how='any')
print 'X shape before na removal=%s_%s' %(X.shape[0],X.shape[1])
xsamples=X.index.astype(str).tolist() # get list of samples to use


yDF=PNP530Cardio126_allFeatures
yname='PNP530Cardio126_allFeatures'
resultFolder='%s/TCR_real_data/PNP530Cardio126Combined/FeaturePredictions/%s/%s_%s_DFs' %(MyPath,modelType,Xname,yname)
if not isdir(resultFolder):
    makedirs(resultFolder)
yMinCol=None
yMaxCol=None
maxiter=1000
modelType='ols'

df=linReg_withWrapper(X,Xname,yDF,yname,resultFolder,modelType,yMinCol,yMaxCol,maxiter)

'''





def linearOrlogistic_reg_withWrapper(X, Xname, yDF, yname, resultFolder, modelType, yMinCol=None, yMaxCol=None, maxiter=35):
    
    resultFileName = 'resultDFlinReg_%s_%s_ycol%s_ycol%s_mi%s.xlsx' % (Xname, yname, yMinCol, yMaxCol, maxiter)
    filenames = [f for f in listdir(resultFolder) if isfile(join(resultFolder, f))]
    if resultFileName not in filenames:
    
        import statsmodels.api as sm
        resultsDF = pd.DataFrame()
        if yMinCol is None:
            yMinCol = 0
        if yMaxCol is None:
            yMaxCol = len(yDF.columns)
            
#         print 'editing sample names:'
#         X=editSampleNames(X)
#         yDF=editSampleNames(yDF)
            
        xsamples = X.index.astype(str).tolist()  # get list of samples to use
    
        for n, target in enumerate(yDF.columns.values[yMinCol:yMaxCol]):
            print n, target
            y = yDF[target]

            # (1) process y

            y = y.dropna(how='any')
            y = y.loc[xsamples]
#             print 'y shape after filtering for x samples=%s' %y.shape[0]
#             print y.shape

            # (2) get seq pres/abse info and fisher test:
            cohortDF = pd.merge(pd.DataFrame(X.isCardio), pd.DataFrame(y), how='left', left_index=True, right_index=True)
            if modelType == 'logit':
                print 'modelType is logit,calculating fisher test...'
                tab = pd.crosstab(cohortDF[target], cohortDF.isCardio)
                abs0 = tab.loc[0, 0]
                abs1 = tab.loc[0, 1]
                pres0 = tab.loc[1, 0]
                pres1 = tab.loc[1, 1]
                OR, p = fisher_exact(tab, alternative='two-sided')

            elif modelType == 'ols':
                print 'modelType is ols,calculating KS and t tests...'
                data = {}
                for name, group in cohortDF.groupby('isCardio'):
                        data[name] = list(group[target])
                dataPNP = data[0]
                dataPNP = [x for x in dataPNP if not np.isnan(x)]
                dataCardio = data[1]
                dataCardio = [x for x in dataCardio if not np.isnan(x)]
                meanPNP = np.mean(dataPNP)
                meanCardio = np.mean(dataCardio)

                s_t, p_t = stats.ttest_ind(dataPNP, dataCardio)
                s_k, p_k = stats.ks_2samp(dataPNP, dataCardio)


            # (3) run model:

            try:
                X = sm.add_constant(X)
                if modelType == 'ols': 
                    model = sm.OLS(y, X, missing='drop')
                elif modelType == 'logit':
                    model = sm.Logit(y, X, missing='drop')
                result = model.fit()
                
                # (4) summarize p-values
                if n == 0:
                    resultDF = pd.DataFrame(result.pvalues).T
                else:
                    resultDF = pd.concat([resultDF, pd.DataFrame(result.pvalues).T], ignore_index=True)
                resultDF.loc[n, 'target'] = target
                if modelType == 'logit':
                    resultDF.loc[n, 'p_fisher'] = p
                    resultDF.loc[n, 'abs_PNP'] = abs0
                    resultDF.loc[n, 'pres_PNP'] = pres0
                    resultDF.loc[n, 'abs_Cardio'] = abs1
                    resultDF.loc[n, 'pres_Cardio'] = pres1
                    resultDF.loc[n, 'prsquared'] = result.prsquared
                    resultDF.loc[n, 'adj_prsquared'] = 1-((result.llf-len(X.columns))/result.llnull)
                elif modelType == 'ols':
                    resultDF.loc[n, 'p_t'] = p_t
                    resultDF.loc[n, 'p_k'] = p_k
                    resultDF.loc[n, 'meanPNP'] = meanPNP
                    resultDF.loc[n, 'meanCardio'] = meanCardio
                    resultDF.loc[n, 'rsquared_adj'] = result.rsquared_adj
                    resultDF.loc[n, 'f_pvalue'] = result.f_pvalue


            except:
                print ' couldnt fit model'
                if n == 0:
                    resultDF = pd.DataFrame()
                else:
                    resultDF.loc[n, 'target'] = target
                    if modelType == 'logit':
                        resultDF.loc[n, 'p_fisher'] = p
                        resultDF.loc[n, 'abs_PNP'] = abs0
                        resultDF.loc[n, 'pres_PNP'] = pres0
                        resultDF.loc[n, 'abs_Cardio'] = abs1
                        resultDF.loc[n, 'pres_Cardio'] = pres1 
                    elif modelType == 'ols':
                        resultDF.loc[n, 'p_t'] = p_t
                        resultDF.loc[n, 'p_k'] = p_k
                        resultDF.loc[n, 'meanPNP'] = meanPNP
                        resultDF.loc[n, 'meanCardio'] = meanCardio


#         resultDF['isCardioBest']=np.where(resultDF.isCardio==resultDF.min(axis=1),1,0)
#         resultDF['bestP']=resultDF.min(axis=1)

        print 'adding info to the full df and saving..'
        
        resultLength = len(resultDF)
        nNans = len(resultDF[resultDF['isCardio'].isnull()])
        nSig = len(resultDF[resultDF['isCardio'] <= 0.05])
        percNan = 100 * float(nNans) / resultLength
        percSig = 100 * float(nSig) / resultLength
        print 'percent null results is %s' % round(percNan, 2)
        print 'percent sig results is %s' % round(percSig, 2)
        
                
        f1 = '%s/resultDF_regType%s_%s_%s_ycol%s_ycol%s_mi%s.xlsx' % (resultFolder, modelType, Xname, yname, yMinCol, yMaxCol, maxiter)
        resultDF.to_excel(f1)
    else:
        print 'this result df already exist in folder'
        resultDF = pd.DataFrame()
    
    return resultDF    
    
        
#--------------------------------------------------------------------------------
'''
the following function predicts binary variable based on phenotype features that is added one variable at a time from a TCRvariableDF. 
this enables to examine the effect of each variable - does it contribute to the prediction (improves the prsquared and/or get its own
significant p-value) or not
the functio also calculates ttest+ks test or Fisher test to check the association between the variable and the binary target regardless of
other parameters

input:
phenotypeDF=sampleXphenotype df
phenotypeDFname=string
TCRvariableDF= TCR feature DFs or TCR sequences
TCRvariableName=String
targetDF = actually a series: the index is sample names and the values are the binary target
targetName=String
resultFolder=folder to save the summarizin DF. when min and max columns are used and the job is paralleled, the seperate resulting DFs
to be concatenated could be found in this folder
modelType= 'ols' or 'logit'. currently the function is built for 'logit' only!
TCRdfMinCol - int or None
TCRdfMaxCol- int or None
maxiter= int, default is 35

usage example;
datasetFolder='%s/TCR_real_data/PNP530Cardio126Combined' %MyPath
datasetName='PNP530Cardio126'

#get feature df:
f1='%s/TCR_real_data/PNP530Cardio126Combined/featureSummaryDFs/PNP530Cardio126_allFeatures' %MyPath
PNP530Cardio126_allFeatures=pd.read_pickle(f1)
print PNP530Cardio126_allFeatures.shape

#get phenotype df:
f2='%s/TCR_real_data/PNP530Cardio126Combined/Phenotypes/PNP530Cardio126_phen.xlsx' %MyPath
PNP530Cardio126_phen=pd.read_excel(f2)
PNP530Cardio126_phen.head()

# Xtypes:
allNum=['BD', 'Age', 'BMI', 'Creatinine', 'isCardio', 'nTemplates', 'Female', 'Smoking_Past', 'Smoking_Yes',
       'PCR_Plate1', 'PCR_Plate10', 'PCR_Plate2', 'PCR_Plate3',
       'PCR_Plate4', 'PCR_Plate5', 'PCR_Plate6', 'PCR_Plate7',
       'PCR_Plate8', 'PCR_Plate9']
small=['BD', 'Age', 'BMI', 'Creatinine', 'isCardio', 'Female','Smoking_Past', 'Smoking_Yes']
smallNoCardio=['BD', 'Age', 'BMI', 'Creatinine', 'Female']
allNumNoPlate=['BD', 'Age', 'BMI', 'Creatinine', 'isCardio', 'nTemplates', 'Female', 'Smoking_Past', 'Smoking_Yes']

#DEFINE WHICH FEATURES TO USE:

Xcols=small
phenotypeDFname='small'
phenotypeDF=PNP530Cardio126_phen[Xcols].set_index('BD')
TCRvariableDF=PNP530Cardio126_allFeatures
TCRvariableName='PNP530Cardio126_allFeatures'
targetDF=phenotypeDF['isCardio']
phenotypeDF=phenotypeDF.drop('isCardio',axis=1)
targetName='isCardio'
modelType='logit'
# print 'target original:'
# print targetDF.head(10)
# print targetDF.value_counts(dropna=False)



resultFolder='%s/TCR_real_data/PNP530Cardio126Combined/%s_Predictions/%s/%s_%s_DFs' %(MyPath,targetName,modelType,
                                                                                      phenotypeDFname,TCRvariableName)
if not isdir(resultFolder):
    makedirs(resultFolder)
    
    
TCRdfMinCol=None
TCRdfMaxCol=None
maxiter=1000



df=predictBinaryPhenotype_logReg_TCRfeatures(phenotypeDF, phenotypeDFname, TCRvariableDF, TCRvariableName, targetDF, 
                                               targetName, resultFolder, modelType, TCRdfMinCol, TCRdfMaxCol, 
                                               maxiter)





'''

def predictBinaryPhenotype_logReg_TCRfeatures(phenotypeDF, phenotypeDFname, TCRvariableDF, TCRvariableName, targetDF, 
                                               targetName, resultFolder, modelType, TCRdfMinCol=None, TCRdfMaxCol=None, 
                                               maxiter=35):
    
    resultFileName = '%s_%s_%s_%s_ycol%s_ycol%s_mi%s.xlsx' % (targetName,modelType,phenotypeDFname, TCRvariableName, TCRdfMinCol, TCRdfMaxCol, maxiter)
    filenames = [f for f in listdir(resultFolder) if isfile(join(resultFolder, f))]
    if resultFileName not in filenames:
    
        import statsmodels.api as sm
        resultsDF = pd.DataFrame()
        if TCRdfMinCol is None:
            TCRdfMinCol = 0
        if TCRdfMaxCol is None:
            TCRdfMaxCol = len(TCRvariableDF.columns)
            
        print 'editing sample names...'
        phenotypeDF=editSampleNames(phenotypeDF)
        TCRvariableDF=editSampleNames(TCRvariableDF)
        targetDF=editSampleNames(targetDF)  
        
        #predict based on phenotypes only:
        print 'predicting based on phenotype only...'
        y=targetDF
        y=y.rename(targetName)
        X=phenotypeDF
        X=X.dropna(how='any')         
        xsamples = X.index.astype(str).tolist()  # get list of samples to use
        y = y.dropna(how='any')
        y = y.loc[xsamples]
        
        print 'y value counts:'
        print y.value_counts(dropna=False)
        if modelType == 'ols': 
            model = sm.OLS(y, X, missing='drop')
        elif modelType == 'logit':
            model = sm.Logit(y, X, missing='drop')
        result = model.fit(maxiter=maxiter)
        print 'phenotype p-values:'
        print pd.DataFrame(result.pvalues).T
        phenotype_prsquared=result.prsquared
        print 'model prsquared is %s' %result.prsquared
        print 'model manual calculated prsquared=%s' %(1-(result.llf/result.llnull))
        phenotype_adj_prsquared=1-((result.llf-len(X.columns))/result.llnull)
        print 'model manual calculated adjusted prsquared=%s' %phenotype_adj_prsquared
                      
        print 'start looping over variables...'
        for n, variable in enumerate(TCRvariableDF.columns.values[TCRdfMinCol:TCRdfMaxCol]):

            print n, variable
            
            #'restart y':
            y=targetDF
            y=y.rename(targetName)
            print 'y name is %s' %targetName
            print y.head()
            print y.value_counts(dropna=False)
            
            print 'merging variable to phenotypeDF to yieldX:'
            X=pd.merge(phenotypeDF,pd.DataFrame(TCRvariableDF[variable]),how='outer',left_index=True,right_index=True)
            X=X.dropna(how='any')
#             print 'X shape after merging is %s_%s' %(X.shape[0],X.shape[1])
#             print X.head()
            
            xsamples = X.index.astype(str).tolist()  # get list of samples to use


            # (1) process y

            y = y.dropna(how='any')
            try:
                y = y.loc[xsamples]
                print 'y shape after filtering for x samples=%s' %y.shape[0]
                print y.shape
            except:
                print 'y shape before filtering for x samples=%s' %y.shape[0]
                print 'couldnt filter for x samples'
                continue

            # (2) get seq pres/abse info and fisher test:
            cohortDF = pd.merge(pd.DataFrame(X[variable]), pd.DataFrame(y), how='left', left_index=True, right_index=True)
            if len(cohortDF[variable].value_counts()) <3:
                try:
                    tab = pd.crosstab(cohortDF[variable], cohortDF[targetName])
#                     print tab
                    abs0 = tab.loc[0, 0]
                    abs1 = tab.loc[0, 1]
                    pres0 = tab.loc[1, 0]
                    pres1 = tab.loc[1, 1]
                    OR, p = fisher_exact(tab, alternative='two-sided')
                    print 'variable is binary,calculating fisher test...'
                except:
                    print 'couldnt execute fisher test'
                    abs0=np.nan
                    abs1=np.nan
                    pres0=np.nan
                    pres1=np.nan
                    p=np.nan

            else:
                try:
                    data = {}
                    for name, group in cohortDF.groupby(targetName):
                            data[name] = list(group[variable])
                    dataPNP = data[0]
                    dataPNP = [x for x in dataPNP if not np.isnan(x)]
                    dataCardio = data[1]
                    dataCardio = [x for x in dataCardio if not np.isnan(x)]
                    meanPNP = np.mean(dataPNP)
                    meanCardio = np.mean(dataCardio)

                    s_t, p_t = stats.ttest_ind(dataPNP, dataCardio)
                    s_k, p_k = stats.ks_2samp(dataPNP, dataCardio)

                    print 'variable is not binary,calculating KS and t tests...'
                except:
                    print 'couldnt execute t and ks test...'
                    meanPNP=np.nan
                    meanCardio=np.nan
                    p_t=np.nan
                    p_k=np.nan
    


            # (3) run model:

            try:
#             X = sm.add_constant(X)
                print 'y value counts:'
                print y.value_counts(dropna=False)
#                 print y.shape
#                 print 'variable value count:'
#                 print X[variable].value_counts(dropna=False)
#                 print X.head()
#                 print X.shape
                if modelType == 'ols': 
                    model = sm.OLS(y, X, missing='drop')
                elif modelType == 'logit':
                    model = sm.Logit(y, X, missing='drop')
                result = model.fit(maxiter=maxiter)

                # (4) summarize p-values
                if n == 0:
                    resultDF = pd.DataFrame(result.pvalues).T.rename(columns={variable:'variable_p'})
                else:
                    resultDF = pd.concat([resultDF, pd.DataFrame(result.pvalues).T.rename(columns={variable:'variable_p'})], ignore_index=True)
                resultDF.loc[n, 'variable'] = variable
                if modelType == 'logit':
                    resultDF.loc[n, 'prsquared'] = result.prsquared
                    resultDF.loc[n, 'adj_prsquared'] = 1-((result.llf-len(X.columns))/result.llnull)
                elif modelType == 'ols':
                    resultDF.loc[n, 'rsquared_adj'] = result.rsquared_adj
                    resultDF.loc[n, 'f_pvalue'] = result.f_pvalue
                if len(cohortDF[variable].value_counts()) <3:
                    resultDF.loc[n, 'p_fisher'] = p
                    resultDF.loc[n, 'abs_PNP'] = abs0
                    resultDF.loc[n, 'pres_PNP'] = pres0
                    resultDF.loc[n, 'abs_Cardio'] = abs1
                    resultDF.loc[n, 'pres_Cardio'] = pres1
                else:     
                    resultDF.loc[n, 'p_t'] = p_t
                    resultDF.loc[n, 'p_k'] = p_k
                    resultDF.loc[n, 'meanPNP'] = meanPNP
                    resultDF.loc[n, 'meanCardio'] = meanCardio


            except:
                print ' couldnt fit model'
                if n == 0:
                    resultDF = pd.DataFrame()
                else:
                    resultDF.loc[n, 'variable'] = variable
                    if len(cohortDF[variable].value_counts()) <3:
                        resultDF.loc[n, 'p_fisher'] = p
                        resultDF.loc[n, 'abs_PNP'] = abs0
                        resultDF.loc[n, 'pres_PNP'] = pres0
                        resultDF.loc[n, 'abs_Cardio'] = abs1
                        resultDF.loc[n, 'pres_Cardio'] = pres1
                    else:     
                        resultDF.loc[n, 'p_t'] = p_t
                        resultDF.loc[n, 'p_k'] = p_k
                        resultDF.loc[n, 'meanPNP'] = meanPNP
                        resultDF.loc[n, 'meanCardio'] = meanCardio
    
        resultDF['varContribute']=np.where(resultDF['adj_prsquared']>phenotype_adj_prsquared,1,0)
#         resultDF['isVariableBest']=np.where(resultDF.variable_p==resultDF.min(axis=1),1,0)
#         resultDF['bestP']=resultDF.min(axis=1)

        print 'adding info to the full df and saving..'
        
#         resultLength = len(resultDF)
#         nNans = len(resultDF[resultDF['isCardio'].isnull()])
#         nSig = len(resultDF[resultDF['isCardio'] <= 0.05])
#         percNan = 100 * float(nNans) / resultLength
#         percSig = 100 * float(nSig) / resultLength
#         print 'percent null results is %s' % round(percNan, 2)
#         print 'percent sig results is %s' % round(percSig, 2)
        
                
        f1 = '%s/resultDF_regType%s_%s_%s_ycol%s_ycol%s_mi%s.xlsx' % (resultFolder, modelType, phenotypeDFname, TCRvariableName, TCRdfMinCol, TCRdfMaxCol, maxiter)
        resultDF.to_excel(f1)
    else:
        print 'this result df already exist in folder'
        resultDF = pd.DataFrame()
    
    return resultDF  


#-------------------------------------------------------------------------------------------------------------
'''
the following function divide cardio samples into groups according to a specific disease-phenotype (Dphenotype) and 
compares features, gene usage and phenotype between groups.

input:
Dphenotype - string
ssCardio - True/False - whether to subsample cardio samples
divideByDphenotype - True/False
ss= int - number of sequences to subsample/to use
repeat - int
phenotypeDF- existing DF or None. if None, will generate phen new dummies file that includes all PNP530 and Cardio126 samples
compare features - True/False whether to compare features
comparePhenotypes- True/False whether to compare phenotypes



'''
def compareGroups_DphenotypeBinary(Dphenotype,ssCardio,divideByDphenotype,ss=None,repeat=None,phenotypeDF=None,
                           compareFeature=True,comparePhenotypes=True):
      
    partialDatasetFolder='TCR_real_data/CardioSamples' 
    datasetName='Cardio126'       
    datasetFolder='%s/%s' %(MyPath,partialDatasetFolder)
    
    # (1) subsample Cardio:
    if ssCardio: 
        print '*********subsampling Cardio cohort with %s sequences per sample:***************'
        nTemplates=ss
        repeat=repeat 
        datasetName_forSS='Cardio126' 
        fullSamplesFolder='%s/SamplesForAnalysis_corrected' %datasetFolder 
        data_folder_full='TCR_real_data/CardioSamples' #change if necessary
        TakeSameSamples=True
        subsampling_and_featureExtraction(fullSamplesFolder,nTemplates,repeat,datasetName_forSS,data_folder_full,
                                          TakeSameSamples=TakeSameSamples)

    #(2): get phenotypeDF if necessary:
    if phenotypeDF is None:
        f1='%s/phenotypicData/Cardio126_phen_new_dummies.xlsx' %datasetFolder
        phen_new_dummies=pd.read_excel(f1)
        f2='%s/phenotypicData/Medical_database_27-06-2018.xlsx' %datasetFolder
        origData=pd.read_excel(f2)
        origDataCols=['RegistrationCode','Hypertension','Dyslipidemia','Microvascular Complications','LDL','HDL',
                      'Triglycerides','HbA1C','PLT','Hemoglobin','Initial CPK','Maximal CPK','LDH','AST','Glucose']
        phenotypeDF=pd.merge(phen_new_dummies,origData[origDataCols],how='left',left_on='RegistrationCode', 
                             right_on='RegistrationCode')
        phenotypeDF=phenotypeDF.replace(9999,np.nan)
        phenotypeDF['PreviousPCI_binary']=np.where(phenotypeDF['Previous PCI'].isnull(),np.nan,
                np.where(phenotypeDF['Previous PCI']>0,1,0))
        f3='%s/phenotypicData/Cardio126_phen_new_dummies_extended.xlsx' %datasetFolder
        phenotypeDF.to_excel(f3)  
        
    #(3) if ss is not None, update folders and names to use
    if ss is not None:
        partialDatasetFolder='TCR_real_data/CardioSamples/Cardio126_SubSampled%sdata_rep%s' %(ss,repeat)
        datasetName='Cardio126_SubSampled%sdata_rep%s' %(ss,repeat)
        
    datasetFolder='%s/%s' %(MyPath,partialDatasetFolder)  
    
    #(4) divide samples according to Dphenotype
    DphenotypeGroupListFolder='%s/phenotypicData/groupSampleLists' %datasetFolder
    if not isdir(DphenotypeGroupListFolder):
        makedirs(DphenotypeGroupListFolder)
        
    categoryCount=phenotypeDF[Dphenotype].value_counts()
    print 'number of categories is %s' %len(categoryCount)
    
    sampleListDict={}
    for categ in categoryCount.index:
        groupDF=phenotypeDF[phenotypeDF[Dphenotype]==categ]
        sampleList=groupDF['BD'].tolist()
        if type(categ)==float:
            categ=int(categ)
        sampleListDict[categ]=sampleList
        print '%s_%s' %(Dphenotype,categ)
        print sampleList
        print len(sampleList)
        
        with open('%s/%s_%s_list' %(DphenotypeGroupListFolder,Dphenotype,categ),'wb') as fp:
            pickle.dump(sampleList,fp)
            
     #(5) compre features
    if compareFeature:
        print '************comparing TCRfeatures between cohorts (not matched!):*************'
        data_folder1=partialDatasetFolder
        data_folder2=partialDatasetFolder
#         datasetName1=datasetName
#         datasetName2=datasetName
        if ss is not None:
            datasetName1='Cardio126_ss%s_rep%s' %(ss,repeat)
            datasetName2='Cardio126_ss%s_rep%s' %(ss,repeat)
        else:
            datasetName1='Cardio126'
            datasetName2='Cardio126'
        TakeSameSamples=False
        filteringList1=sampleListDict.values()[0]
        filteringList2=sampleListDict.values()[1]
        filteringList1Name='_'.join([Dphenotype,str(sampleListDict.keys()[0])])
        filteringList2Name='_'.join([Dphenotype,str(sampleListDict.keys()[1])])
        plotType='bar'

        compare_features_between_datasets(data_folder1, datasetName1,  data_folder2, datasetName2,TakeSameSamples, 
                                      filteringList1,filteringList2,filteringList1Name,filteringList2Name)
        plot_gene_usage_comparison(data_folder1, datasetName1,  data_folder2, datasetName2,plotType,TakeSameSamples, 
                                          filteringList1,filteringList2,filteringList1Name,filteringList2Name)
    
    
    
    # (6) compare phenotypes:
    if comparePhenotypes:
        numericalphenotypes=['Age','BMI','eGFR by CKD-EPI','WBC','LDL','HDL','Triglycerides','CRP','HbA1C','PLT','Hemoglobin',
        'Initial CPK','Maximal CPK','LDH','AST','Glucose','Initial Troponin','Maximal Troponin']
        categoricalphenotypes=['Gender', 'Smoking Status','Hypertension','Dyslipidemia','Microvascular Complications',
                               'Glucose Disorder','Previous CABG','PreviousPCImapped','Admission Diagnosis','LVEF','Known CAD']
        
        phenotypeDF1=phenotypeDF.set_index('BD')
        phenotypeDF2=phenotypeDF.set_index('BD')
        sampleList1=sampleListDict.values()[0]
        sampleList2=sampleListDict.values()[1]
        datasetName1=datasetName
        datasetName2=datasetName
        sampleListName1='_'.join([Dphenotype,str(sampleListDict.keys()[0])])
        sampleListName2='_'.join([Dphenotype,str(sampleListDict.keys()[1])])
        nBins=20

        fig1=compare_phenotypes(numericalphenotypes,categoricalphenotypes,nBins,phenotypeDF1,sampleList1,sampleListName1,datasetName1,
                               phenotypeDF2,sampleList2,sampleListName2,datasetName2,widthUnit=3,heightUnit=5)
        folderToSave='%s/phenotypicData/phenotypeComparisons' %datasetFolder
        if not isdir(folderToSave):
            makedirs(folderToSave)
        fig1.savefig('%s/%s_%s' %(folderToSave,datasetName1,datasetName2),dpi=300)


    print 'end of function!!!' 
    
    
#----------------------------------------------------------------------------------------
'''
the following function generates different TCRdf for each group of Dphenotype in the cardio cohort
input:
Dphenotype - string
ssCardio-True/False
percShared=int 
ss-int
repeat-int
phenotypeDF-existing DF or None

usage example:
Dphenotype='Known CAD'
ssCardio=False
percShared=10

genTCRdfForDphenotypeGroup(Dphenotype,ssCardio,percShared)

'''
def genTCRdfForDphenotypeGroup(Dphenotype,ssCardio,percShared,ss=None,repeat=None,phenotypeDF=None):
    
    partialDatasetFolder='TCR_real_data/CardioSamples' 
    datasetName='Cardio126'       
    datasetFolder='%s/%s' %(MyPath,partialDatasetFolder)

    
    #(1): get phenotypeDF if necessary:
    if phenotypeDF is None:
        f1='%s/phenotypicData/Cardio126_phen_new_dummies.xlsx' %datasetFolder
        phen_new_dummies=pd.read_excel(f1)
        f2='%s/phenotypicData/Medical_database_27-06-2018.xlsx' %datasetFolder
        origData=pd.read_excel(f2)
        origDataCols=['RegistrationCode','Hypertension','Dyslipidemia','Microvascular Complications','LDL','HDL',
                      'Triglycerides','HbA1C','PLT','Hemoglobin','Initial CPK','Maximal CPK','LDH','AST','Glucose']
        phenotypeDF=pd.merge(phen_new_dummies,origData[origDataCols],how='left',left_on='RegistrationCode', 
                             right_on='RegistrationCode')
        phenotypeDF=phenotypeDF.replace(9999,np.nan)
        phenotypeDF['PreviousPCI_binary']=np.where(phenotypeDF['Previous PCI'].isnull(),np.nan,
                np.where(phenotypeDF['Previous PCI']>0,1,0))
        f3='%s/phenotypicData/Cardio126_phen_new_dummies_extended.xlsx' %datasetFolder
        phenotypeDF.to_excel(f3)  
    
    #(2) if ss is not None, update folders and names to use
    if ss is not None:
        partialDatasetFolder='TCR_real_data/CardioSamples/Cardio126_SubSampled%sdata_rep%s' %(ss,repeat)
        datasetName='Cardio126_SubSampled%sdata_rep%s' %(ss,repeat)
        datasetFolder='%s/%s' %(MyPath,partialDatasetFolder)
    
    #(3) divide samples according to Dphenotype
    DphenotypeGroupListFolder='%s/phenotypicData/groupSampleLists' %datasetFolder
    if not isdir(DphenotypeGroupListFolder):
        makedirs(DphenotypeGroupListFolder)
        
    categoryCount=phenotypeDF[Dphenotype].value_counts()
    print 'number of categories is %s' %len(categoryCount)
    
    sampleListDict={}
    for categ in categoryCount.index:
        groupDF=phenotypeDF[phenotypeDF[Dphenotype]==categ]
        sampleList=groupDF['BD'].tolist()
        if type(categ)==float:
            categ=int(categ)
        sampleListDict[categ]=sampleList
        print '%s_%s' %(Dphenotype,categ)
        print sampleList
        print len(sampleList)
        
        with open('%s/%s_%s_list' %(DphenotypeGroupListFolder,Dphenotype,categ),'wb') as fp:
            pickle.dump(sampleList,fp)
            
        DFtype='TCR'
        genDF=True #False=generate new df
        toBinary=True
        mbLevel='g'
        useShortName=True
        datasetFolder=datasetFolder
        datasetName=datasetName
        minVal=None    
        minSharedT=None 
        percShared=percShared 
        removeOutliers=True
        normData=False
        logTransform=True
        extractUniqueAA=True 
        minNshared=2
        onlyProductive=True
        mbDataFolder='AllSeqProjects'
        SampleList=sampleList
        SampleListName='%s_%s' %(Dphenotype,categ)
        libPrepMethod=None

        genTCRorMBdfWithManipulations(DFtype, genDF, toBinary, removeOutliers, normData, logTransform,
                             minVal, minSharedT, percShared,
                             mbLevel, useShortName, datasetFolder, datasetName,extractUniqueAA,
                            minNshared, onlyProductive, mbDataFolder, SampleList,
                              SampleListName)
        

#--------------------------------------------------------------------
'''
the following functions are used together to compare sharing rates between PNP and Cardio cohorts
which are defined by age range, and by the AllUniqueWithCounts_path used (usually use the one 
from the PNP530Cardio126 combined folder, but can also use those from matched or unmatched
ss cohort

run the last function to execute. it calls all other functions

USAGE EXAMPLE:
min_age=50
max_age=60
nToSample=20
nSeqs=7500
isRandom=False

nSampledShared=2

sharing_rate_comparison_PNP_Cardio(AllUniqueWithCounts_path=None,min_age=min_age,
                max_age=max_age,nToSample=nToSample,nSeqs=nSeqs,
                nSampledShared=nSampledShared,isRandom=isRandom,
                random_state_PNP=1,random_state_Cardio=1)
'''
        
        
def get_data(AllUniqueWithCounts_path):
    
    ### this function gets all data required for sharing comparisons between PNP and Cardio
    ### (Age, gender, isCardio, AllUniqueWithCounts)
    
    print 'age data'
    Age=pd.read_excel('%s/TCR_real_data/PNP530Cardio126Combined/Phenotypes/\
PNP530Cardio126_Age.xlsx' %MyPath).set_index('BD')
    print Age.shape
    print Age.head()

    print 'gender data'
    Gender=pd.read_excel('%s/TCR_real_data/PNP530Cardio126Combined/Phenotypes/\
PNP530Cardio126_Gender_Male.xlsx' %MyPath).set_index('BD')
    print Gender.shape
    print Gender.head()
    
    print 'isCardio data'
    isCardio=pd.read_pickle(PRED_RESULTS_DIR+'TargetDFs/isCardio.dat')

    print 'merging all'
    AgeGender=pd.merge(Age,Gender,how='inner',left_index=True,right_index=True).dropna()
    AgeGenderIsCardio=pd.merge(AgeGender,isCardio,how='inner',left_index=True,right_index=True).dropna()
    print AgeGenderIsCardio.shape
    print AgeGenderIsCardio.head()
    
    print 'AllUniqueWithCounts'
    if AllUniqueWithCounts_path is None:
        AllUniqueWithCounts_path='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/sharingAnalysis/\
AllUniqueWithCounts'
    AllUniqueWithCounts=pd.read_pickle(AllUniqueWithCounts_path)

    print AllUniqueWithCounts.shape
    print AllUniqueWithCounts.head()
    
    print 'AllUniqueWithCounts_prod'
    AllUniqueWithCounts_prod=AllUniqueWithCounts[AllUniqueWithCounts['prod_stat']==1]
    print AllUniqueWithCounts_prod.shape
    print AllUniqueWithCounts_prod.head()
    
    return AllUniqueWithCounts_prod,AgeGenderIsCardio


def calc_seq_per_sample_stats(AllUniqueWithCounts_prod,nToSample):
    
    #this function prints the number of sequences in the 10t,20th and 30th sample in the df
    ## it plots the distribution of number of sequences per sample
    ## and returns a series with number of seqeunces per sample, ordered from highest to lowest
    
    
    AllUniqueWithCounts_prod_count=AllUniqueWithCounts_prod.groupby('Sample').count()['isPublic']
    
    print AllUniqueWithCounts_prod_count.describe()
    AllUniqueWithCounts_prod_count=AllUniqueWithCounts_prod_count.sort_values(ascending=False)
    
    for n in [10,20,30,nToSample]:
        try:
            print ('n seqs in %sth sample= ' %n,AllUniqueWithCounts_prod_count.iloc[n-1])
        except:
            print ('there are only %s samples' %len(AllUniqueWithCounts_prod_count))
                
    fig,ax=plt.subplots()
    ax.hist(AllUniqueWithCounts_prod_count,bins=25)
    plt.show()
    
    return AllUniqueWithCounts_prod_count


def get_samples_for_age_decade(min_age,max_age,AgeGenderIsCardio,AllUniqueWithCounts_prod,nToSample=20,nSeqs=None,isRandom=True,random_state_PNP=1,random_state_Cardio=1):
    
    ### this function takes the age range, number of samples and number of sequences to subsample.
    ### first it gets all PNP and Cardio samples in the age range.
    ### than it checks how many sequences are in each of those samples
    ### than it selects samples - if there are enough samples with the desired number of sequences, 
    ### it selects samples randomly, or the samples with the highest number of sequences (according to)
    ### the isRandom parameter)
    ### if there are not enough samples, it won't be able to return dataframes and an error will occure
    
    ###get all avaiable PNP and Cardio samples for the age range
    #PNP:
    df=AgeGenderIsCardio

    print 'PNP samples'
    PNP_males_ageDecade=df[(df['Age']>min_age) & (df['Age']<=max_age) & (df['Gender_Male']==1) & (df['isCardio']==0)]
       
    L=drop_relatives(PNP_males_ageDecade,onlyBloodRels=False)
    print ('relatives in list: ',L)
    if len(L)>0:
        print ('dropping the following samples because they are relatives: ', L)
        PNP_males_ageDecade=PNP_males_ageDecade.drop(L)
        
    print 'removing samples that are not in AllUniqueWithCounts_prod:'
    PNP_males_ageDecade=PNP_males_ageDecade[PNP_males_ageDecade.index.isin(AllUniqueWithCounts_prod.Sample.unique())]
    
    print PNP_males_ageDecade.shape
    print PNP_males_ageDecade.head()
    
    print 'getting stats on PNP samples:'
    AllUniqueWithCounts_prod_count_PNP=\
calc_seq_per_sample_stats(AllUniqueWithCounts_prod[AllUniqueWithCounts_prod['Sample'].isin(PNP_males_ageDecade.index)],
                         nToSample)

    #Cardio
    print 'cardio samples'
    Cardio_males_ageDecade=df[(df['Age']>min_age) & (df['Age']<=max_age) & (df['Gender_Male']==1) & (df['isCardio']==1)]
       
    print 'removing samples that are not in AllUniqueWithCounts_prod:'
    Cardio_males_ageDecade=Cardio_males_ageDecade[Cardio_males_ageDecade.index.isin(AllUniqueWithCounts_prod.Sample.unique())]
    
    print Cardio_males_ageDecade.shape
    print Cardio_males_ageDecade.head()
    
    print 'getting stats on Cardio samples:'
    AllUniqueWithCounts_prod_count_Cardio=\
calc_seq_per_sample_stats(AllUniqueWithCounts_prod[AllUniqueWithCounts_prod['Sample'].isin(Cardio_males_ageDecade.index)],
                         nToSample)
    
    
    ###check that there are enough samples with the desired nSeqs, and select samples:
    
    
    
    #PNP
    if nSeqs is None:
        if isRandom:
            print 'randomly selecting %s samples from PNP samples' %nToSample
            PNP_males_ageDecade_s=PNP_males_ageDecade.sample(n=nToSample,random_state=random_state_PNP)
        else:
            print 'selecting %s samples with the highest number of sequences' %nToSample
            PNP_males_ageDecade_s=PNP_males_ageDecade.loc[ AllUniqueWithCounts_prod_count_PNP[:nToSample].index,:] 
    elif AllUniqueWithCounts_prod_count_PNP.iloc[nToSample-1]<nSeqs:
        print ('nSeq= ', nSeqs)
        print ('n seqs in %sth sample: ' %(nToSample-1), AllUniqueWithCounts_prod_count_PNP.iloc[nToSample-1])
        print 'not enough samples with the desired number of seqs in PNP'
    elif AllUniqueWithCounts_prod_count_PNP.iloc[nToSample-1]==nSeqs:
        print 'number of samples with the desired number of seqs is equal to the number of desired samples (PNP)'
        PNP_males_ageDecade_s=PNP_males_ageDecade
    else:
        print 'There are more samples with the desired number of seqs than the number of desired samples (PNP)'
        if isRandom:
            print 'randomly selecting %s samples from PNP samples' %nToSample
            PNP_males_ageDecade_s=PNP_males_ageDecade.sample(n=nToSample,random_state=random_state_PNP)
        else:
            print 'selecting %s samples with the highest number of sequences' %nToSample
            PNP_males_ageDecade_s=PNP_males_ageDecade.loc[ AllUniqueWithCounts_prod_count_PNP[:nToSample].index,:]
    print PNP_males_ageDecade_s.shape
    print PNP_males_ageDecade_s.head()
    PNP_males_ageDecade_samples=PNP_males_ageDecade_s.index.tolist()
    
    #Cardio:
    if nSeqs is None:
        if isRandom:
            print 'randomly selecting %s samples from Cardio samples' %nToSample
            Cardio_males_ageDecade_s=Cardio_males_ageDecade.sample(n=nToSample,random_state=random_state_Cardio)
        else:
            print 'selecting %s samples with the highest number of sequences' %nToSample
            Cardio_males_ageDecade_s=Cardio_males_ageDecade.loc[AllUniqueWithCounts_prod_count_Cardio[:nToSample].index,:]
    elif AllUniqueWithCounts_prod_count_Cardio.iloc[nToSample-1]<nSeqs:
        print 'not enough samples with the desired number of seqs in Cardio'
    elif AllUniqueWithCounts_prod_count_Cardio.iloc[nToSample-1]==nSeqs:
        print 'number of samples with the desired number of seqs is equal to the number of desired samples (Cardio)'
        Cardio_males_ageDecade_s=Cardio_males_ageDecade
    else:
        print 'There are more samples with the desired number of seqs than the number of desired samples (Cardio)'
        if isRandom:
            print 'randomly selecting %s samples from Cardio samples' %nToSample
            Cardio_males_ageDecade_s=Cardio_males_ageDecade.sample(n=nToSample,random_state=random_state_Cardio)
        else:
            print 'selecting %s samples with the highest number of sequences' %nToSample
            Cardio_males_ageDecade_s=Cardio_males_ageDecade.loc[AllUniqueWithCounts_prod_count_Cardio[:nToSample].index,:]
 
    print Cardio_males_ageDecade_s.shape
    print Cardio_males_ageDecade_s.head()
    Cardio_males_ageDecade_samples=Cardio_males_ageDecade_s.index.tolist()
        
    ###comparing age distributions:    
    print 'comparing age distribution'
    dataList=[('Healthy' ,PNP_males_ageDecade_s.Age.tolist()),('Patients' ,Cardio_males_ageDecade_s.Age.tolist())]
    title=''

    fig,ax=plt.subplots()
    ax,ks_p_cohort1_cohort2,t_p_cohort1_cohort2,p_Anov,filename=plotHistComprison(dataList,ax,title,showLegend=True,nBins=5,
                            toAnnotate=True, colorList=['Grey','darkred'],alpha=None, plotType='hist')

    plt.show()
    print ('n samples in PNP sample list: ', len(PNP_males_ageDecade_samples))
    print ('n samples in Cardio sample list: ', len(Cardio_males_ageDecade_samples))
    
    return fig,PNP_males_ageDecade_samples,Cardio_males_ageDecade_samples
    
    
    
def subsample_sequences_for_samples(AllUniqueWithCounts_prod,PNP_males_ageDecade_samples,Cardio_males_ageDecade_samples,nSeqs=None):
    
    ### this function selects randomly the desired number of sequences from each sample in the PNP and Cardio dfs,
    ### if there is a sample with not enough sequences, it will fail
    print ('n samples in PNP sample list: ', len(PNP_males_ageDecade_samples))
    print ('n samples in Cardio sample list: ', len(Cardio_males_ageDecade_samples))
    
    PNP_data=AllUniqueWithCounts_prod[AllUniqueWithCounts_prod['Sample'].isin(PNP_males_ageDecade_samples)]
    Cardio_data=AllUniqueWithCounts_prod[AllUniqueWithCounts_prod['Sample'].isin(Cardio_males_ageDecade_samples)]
    
    print ('PNP_data.Sample.nunique(): ',PNP_data.Sample.nunique())
    print ('Cardio_data.Sample.nunique(): ',Cardio_data.Sample.nunique())
    
    print 'generating df containing selected PNP and Cardio samples...'
    df=pd.concat([PNP_data,Cardio_data])
    print df.shape
    print df.iloc[:5,:5]
    
    if nSeqs is None:
        print 'nSeqs is None, no subsampling, all sequences are used'
        subsampled_df=df
    else:
        subsampled_df_list=[]
        print 'subsampling sequences...'
        for n,sample in enumerate(df['Sample'].unique().tolist()):
            if n%10==0: print n
            df2=df[df['Sample']==sample]
            df_s=df2.sample(n=nSeqs,random_state=1)
            subsampled_df_list.append(df_s)

        print 'concatenating subsampled dfs'
        subsampled_df=pd.concat(subsampled_df_list)
        print subsampled_df.shape
        print subsampled_df.head()

    return subsampled_df
    


def calc_sharing_rate_per_sample(subsampled_df,PNP_males_ageDecade_samples,Cardio_males_ageDecade_samples,nseqs,nSampledShared):
    subsampled_df=subsampled_df.reset_index()
    
    dataList=[]
    for item in [('Healthy',PNP_males_ageDecade_samples),('Patients',Cardio_males_ageDecade_samples)]:
        dataset_name=item[0]
        sample_list=item[1]
        
        print ('******dataset: '+ dataset_name+'********')
        
        subsampled_df_dataset=subsampled_df[subsampled_df.Sample.isin(sample_list)]
        sampleBySeq_ss=subsampled_df_dataset.pivot(index='Sample',columns='index',values='prod_stat').fillna(0)
    
        print sampleBySeq_ss.shape
        print sampleBySeq_ss.iloc[:5,:5]
        print sampleBySeq_ss.sum().min()
        print sampleBySeq_ss.sum().max()
    
        public_seqs=sampleBySeq_ss.sum()[sampleBySeq_ss.sum() >= nSampledShared].index
        sampleBySeq_ss_public=sampleBySeq_ss.loc[:,public_seqs]
        print ('sampleBySeq_ss.shape:' ,sampleBySeq_ss.shape)
        print ('sampleBySeq_ss_public.shape:' ,sampleBySeq_ss_public.shape)
    
        seqSum=sampleBySeq_ss.sum(axis=1).rename('seqSum')
        publicSeqSum=sampleBySeq_ss_public.sum(axis=1).rename('publicSeqSum')
        merged=pd.merge(pd.DataFrame(seqSum),pd.DataFrame(publicSeqSum),how='inner',left_index=True,right_index=True)
        print ('merged.shape:' ,merged.shape)
        merged['perc_public']=100*merged['publicSeqSum'].astype(float) / merged['seqSum']
        print 'merged:'
        print merged.head(10)
    
        perc_public_list=merged['perc_public'].tolist()
        print ('perc_public_list[:10]: ',perc_public_list[:10])
    
        dataList.append(tuple([item[0],perc_public_list]))
    
    fig,ax=plt.subplots()
    title=''
    
    s_MW,p_MW=mannwhitneyu(dataList[0][1],dataList[1][1])
    ax,ks_p_cohort1_cohort2,t_p_cohort1_cohort2,p_Anov,filename=plotHistComprison(dataList,ax,
                        title,showLegend=True,nBins=20,toAnnotate=True,alpha=None,plotType='kde',
                        colorList=['grey','darkred'])
    ax.text(0.02,0.05,'p-value=%s (MW)' %p_MW, transform=ax.transAxes,ha='left',va='top')
    plt.show()
    
    return fig,ax,s_MW,p_MW

def sharing_rate_comparison_PNP_Cardio(AllUniqueWithCounts_path=None,min_age=50,max_age=60,nToSample=20,
                                       nSeqs=2000,isRandom=True,nSampledShared=2,random_state_PNP=1,random_state_Cardio=1):
    
    ## get data:
    
    AllUniqueWithCounts_prod,AgeGenderIsCardio=get_data(AllUniqueWithCounts_path)
        
    ## select samples
    fig,PNP_males_ageDecade_samples,Cardio_males_ageDecade_samples=get_samples_for_age_decade(min_age,max_age,AgeGenderIsCardio,
                    AllUniqueWithCounts_prod,nToSample=nToSample,
                    nSeqs=nSeqs,isRandom=isRandom,random_state_PNP=random_state_PNP,random_state_Cardio=random_state_Cardio)
    
    ##subsample sequences:
    subsampled_df = subsample_sequences_for_samples(AllUniqueWithCounts_prod,PNP_males_ageDecade_samples,
                                                    Cardio_males_ageDecade_samples,nSeqs=nSeqs)
    
    #compare sharing rates:
    fig,ax,s_MW,p_MW = calc_sharing_rate_per_sample(subsampled_df,PNP_males_ageDecade_samples,Cardio_males_ageDecade_samples,
                                                    nSeqs,nSampledShared)
    
    #save figure:
    
    
    return
    
    
#-----------------------------------------------------------------


