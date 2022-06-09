from os import listdir,mkdir
from os.path import isfile, join, isdir,exists
import pandas as pd
import numpy as np
from scipy import stats
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot,draw_correlation_scatter
from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
from Bio.SeqUtils import GC
import seaborn as sns
import random
from scipy.stats import pearsonr
from skbio.diversity.alpha import shannon, simpson, berger_parker_d

from pop_organize import get_sample_data, get_sample_with_dfs
from SufficientStatistics import *
from MyFunctionsShani import *
import math
from myplots import roundup, rounddown, find_decimal_fold
from skbio.stats.distance import mantel
from scipy.spatial.distance import braycurtis, pdist


from queue.qp import qp,fakeqp
from addloglevels import sethandlers
import logging 
from Utils import cacheOnDisk
import os
from skbio.stats.distance import mantel
from Feature_phenotype_functions import common_processing_feature_phenotype_matrices


#------------------------------------------------------------------------------------------------------------
basePath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/NewMantelFromEclipse2'


#------------------------------------------------------------------------------------------------------

def plot_corr_for_distMats(featureDistMat,PhenotypeDistMat,FeatureName,PhenotypeName,minPhenotypeValue):
    
    ##x=feature matrix
    ##y=phenotype matrix

    ##x and y should be preprocessed to be at the same shape!
    
    
    # generate lists of feature and phenotype distances for each sample combination:
    indList=[]
    xValueList=[]
    yValueList=[]

    for i,row in enumerate(featureDistMat.index):
        for j,column in enumerate(featureDistMat.columns.values):
            if row!=column:
                ind=row+'_'+column
                indList.append(ind)
                xValue=featureDistMat.loc[row,column]
                xValueList.append(xValue)
                yValue=PhenotypeDistMat.loc[row,column]
                yValueList.append(yValue)

    x=xValueList
    y=yValueList

    fig1=plt.figure(figsize=(8,8))

    ymean=np.mean(y)                                                                                                                      

    plt.scatter(x,y, alpha=0.1)
    plt.xlabel(FeatureName,fontsize=14)
    plt.ylabel(PhenotypeName,fontsize=14)
    plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)),c='blue',linewidth=1)
    plt.title('Distance matrix correlation',fontsize=16)
    
    from scipy.stats import pearsonr
    r,p = pearsonr(x,y)
    
    nSamples=len(featureDistMat)


    plt.annotate("r=%.4f p=%.6f, nSamples=%s" %(r,p,nSamples),  xy=(0.02, 0.96), xycoords='axes fraction', fontsize=14,
        horizontalalignment='left', verticalalignment='top', fontweight='bold')
    
    if minPhenotypeValue is not None:
        plt.ylim(minPhenotypeValue,np.max(y)*1.1)
#         plt.margins(0.2)
    
#     file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/DistMat_correlation_plots/'
#     fig1.savefig(file1,dpi=200)


    return fig1

#------------------------------------------------------------------------------------------------------


def mantelTest_feature_phenotype(feature_dist_file, phenotype_dist_file,feature_name,phenotype_name, method,n_permut,alternative,removeSameUser):
    
    
    
    MantelFolder='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/NewMantelTestResults_newFeatures'
    MantelName='MantelDF_%s_%s_%s_%s_%s' %(feature_name,phenotype_name, method,n_permut,alternative)
    MantelFile='%s/%s' %(MantelFolder,MantelName)
    existingMantels=[f for f in listdir( MantelFolder) if isfile(join(MantelFolder, f))]
    
    if MantelName not in existingMantels:
    
        #load distance matrix files:
        print 'loading and processing distance matrix files...'
        feature_dist_mat=pd.read_pickle(feature_dist_file)
        phenotype_dist_mat=pd.read_pickle(phenotype_dist_file)

        x=feature_dist_mat
        y=phenotype_dist_mat

        x,y=common_processing_feature_phenotype_matrices(x,y,removeSameUser)
        
        print 'running mantel test...'
        corr_coeff,p_value,N = mantel(x, y, method=method, permutations=n_permut, alternative=alternative, strict=True, lookup=None)

        mantelDF=pd.DataFrame()
        mantelDF.loc[0,'feature_name']=feature_name
        mantelDF.loc[0,'phenotype_name']=phenotype_name
        mantelDF.loc[0,'method']=method
        mantelDF.loc[0,'n_permut']=n_permut
        mantelDF.loc[0,'alternative']=alternative
        mantelDF.loc[0,'corr_coeff']=corr_coeff
        mantelDF.loc[0,'p_value']=p_value
        mantelDF.loc[0,'feature_name']=feature_name
        mantelDF.loc[0,'N']=N
        
        mantelDF.to_pickle(MantelFile)
        print p_value
        
        #generate correlation plot if desired:
        if p_value<0.05 or corr_coeff>0.15:
            fig1=plot_corr_for_distMats(x,y,feature_name,phenotype_name,minPhenotypeValue=None)
            file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/DistMat_correlation_plots_newFeatures/%s_%s' %(feature_name,phenotype_name)
            fig1.savefig(file1,dpi=200)
        
        
        print 'DONE!'
    
    else:
        print 'This Mantel Test result already exist in folder'
    
    
#-----------------------------------------------------------------------------------------------------

PhenotypeDistFolder='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/newPhenotypes'
phenotypeFiles = [f for f in listdir(PhenotypeDistFolder) if isfile(join(PhenotypeDistFolder, f))]
# phenotypeFiles = ['distMat_PNP434_Albumin_euclidean','distMat_PNP434_Hips_euclidean']
# phenotypeFiles = ['distMat_PNP434_Lymphocytes %_euclidean','distMat_PNP434_Age_euclidean','distMat_PNP434_Neutrophils %_euclidean',\
# 'distMat_PNP434_Chloride - blood_euclidean','distMat_PNP434_Total Bilirubin_euclidean','distMat_PNP434_Cholesterol, total - blood_euclidean',\
# 'distMat_PNP434_WHR_euclidean','distMat_PNP434_Total Protein_euclidean','distMat_PNP434_Triglycerides_euclidean']
FeatureDistFolder='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/NewFeatures'
featureFiles = [f for f in listdir(FeatureDistFolder) if isfile(join(FeatureDistFolder, f))]

#-------------------------------------------------------------------------------------------------------------


## send jobs to the queue (code adopted from Elad, see lab wiki for more details:
sethandlers()
os.chdir(basePath)
## the inputs for qp are: jobname, q=machine list, *** add max_r to prevent exploding the cluster!!***
with qp('Mantel_job',  q = ['himemint.q','himem7.q'], mem_def="10G", trds_def=2, deleteCSHwithnoerr = True, tryrerun = False, max_u=120) as q:
    q.startpermanentrun()
    wait_on =[]
    
##now define a loop that divide the job and send each part seperatly:

method='spearman'
n_permut=99999
alternative='greater'


for i,fileNameF in enumerate(featureFiles):
#     if i<2:
        feature_dist_file='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/NewFeatures/%s' %fileNameF
        feature_name=fileNameF.replace('distMat_PNP434_','')
#         fileNameF.split('_')[-2]+'_'+fileNameF.split('_')[-1]
#         print i,feature_name

        for j,fileNameP in enumerate(phenotypeFiles):
#             if j<4:
                phenotype_dist_file='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/newPhenotypes/%s' %fileNameP
                phenotype_name=fileNameP.split('_')[-2]+'_'+fileNameP.split('_')[-1]
                print i, j, feature_name, phenotype_name
               

#                 mantelTest_feature_phenotype(feature_dist_file, phenotype_dist_file,feature_name,phenotype_name, method,n_permut,alternative)

                wait_on.append(q.method(mantelTest_feature_phenotype,kwargs={'feature_dist_file':feature_dist_file,'phenotype_dist_file':phenotype_dist_file,
                                                                     'feature_name':feature_name,'phenotype_name':phenotype_name,
                                                                     'method':method, 'n_permut':n_permut, 'alternative':alternative,'removeSameUser':True}))
       
q.wait(wait_on)
