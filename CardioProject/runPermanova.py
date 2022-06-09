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
from Feature_phenotype_functions import common_processing_feature_phenotype_matrices
import os
from Utils import cacheOnDisk
from queue.qp import qp,fakeqp
from addloglevels import sethandlers


import time
cdate=str(time.strftime("%d%m%Y"))
cdate

#----------------------------------------------------------------------------------------------
basePath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/NewPERMANOVAFromEclipse'

#----------------------------------------------
file6='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/fullXgroupbyBD_only434'
fullXgroupbyBD_only434=pd.read_pickle(file6)

#-------------------------------------------------------------------------------------------
# @cacheOnDisk(basePath=basePath, filename='PERMANOVA_%(i)s_%(j)s_%(nPerm)s', force=True)
def compare_feature_distance_for_binary_phenotypes(feature_dist_file,feature_name,phenotype,nPerm,removeSameUser,i,j,filteringList):
    
    from skbio.stats.distance import permanova
    from skbio.stats.distance import DistanceMatrix
    
    permanovaFolder='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/permAnovaResults_FemalesOnly_newfeatures'
    permanovaName='%s_%s_%s_%s' %(feature_name,phenotype,nPerm,removeSameUser)
    permanovaFile='%s/%s' %(permanovaFolder,permanovaName)
    existingPermanovas=[f for f in listdir( permanovaFolder) if isfile(join(permanovaFolder, f))]
    
    if permanovaName not in existingPermanovas:


        #(1)load and process distance matrix files:
        print 'loading and processing distance matrix files...'
        feature_dist_mat=pd.read_pickle(feature_dist_file)
        
    
        x=feature_dist_mat
        y=pd.DataFrame(fullXgroupbyBD_only434[phenotype])
        
        x,y=common_processing_feature_phenotype_matrices(x,y,True)
        
        if filteringList is not None:
            for sample in x.index:
                if sample not in filteringList:
                    x=x.drop(sample,axis=0)
                    x=x.drop(sample,axis=1)
            for sample in y.index:
                if sample not in filteringList:
                    y=y.drop(sample,axis=0)
    #                 y=y.drop(sample,axis=1)
            
            print ' filtered matrices to include only samples in filteringList:'
            print 'new x array shape is %s_%s' %(x.shape[0],x.shape[1])
            print 'new y array shape is %s_%s' %(y.shape[0],y.shape[1])
    
        #(2)
        
        DM=DistanceMatrix(x,ids=x.index)
        results=permanova(DM, y, column=phenotype, permutations=nPerm)
    #     print results
        p=results['p-value']
        s=results['test statistic']
    #     print p
        
        permAnovaResDF=pd.DataFrame()
        permAnovaResDF.loc[0,'featureName']=feature_name
        permAnovaResDF.loc[0,'phenotype']=phenotype
        permAnovaResDF.loc[0,'nPerm']=nPerm
        permAnovaResDF.loc[0,'removeSameUser']=removeSameUser
        permAnovaResDF.loc[0,'s']=s
        permAnovaResDF.loc[0,'p']=p
    
#         file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/permAnovaResults_FemalesOnly_newfeatures/%s_%s_%s_%s' %(feature_name,phenotype,nPerm,removeSameUser)
        permAnovaResDF.to_pickle(permanovaFile)
        
        print p
        print 'done permanova'
    else:
        print 'This permanova test result already exist in folder'
    
    #------------------------------------------
    
FeatureDistFolder='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/NewFeatures'
featureFiles = [f for f in listdir(FeatureDistFolder) if isfile(join(FeatureDistFolder, f))]    
    

## send jobs to the queue (code adopted from Elad, see lab wiki for more details:
sethandlers()
os.chdir(basePath)
## the inputs for qp are: jobname, q=machine list, *** add max_r to prevent exploding the cluster!!***
with qp('PERMANOVA_job',  q = ['himemint.q','himem7.q'], mem_def="10G", trds_def=2, deleteCSHwithnoerr = True, tryrerun = False, max_u=120) as q:
    q.startpermanentrun()
    wait_on =[]
    
##now define a loop that divide the job and send each part seperatly:

nPerm=9999
removeSameUser=True

with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/PNPfemales','rb') as fp:
    PNPfemales=pickle.load(fp)
filteringList=PNPfemales

FeatureDistFolder='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/NewFeatures'
featureFiles = [f for f in listdir(FeatureDistFolder) if isfile(join(FeatureDistFolder, f))]

# featureFilesToDrop=['distMat_PNP434_sharingMatrixMoreThan1RAlog2scale_euclidean','distMat_PNP434_sharingMatrixMoreThan2_euclidean','distMat_PNP434_sharingMatrixMoreThan2RA_braycurtis',
#                     'distMat_PNP434_sharingMatrixMoreThan5_euclidean','distMat_PNP434_sharingMatrixMoreThan5RA_braycurtis']
# for toDrop in featureFilesToDrop:
#     featureFiles.remove(toDrop)
    
#get binary variable list:
file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/binaryVariableList.txt'

with open(file1,'rb') as fp:
    binaryVariableList=pickle.load(fp)
    
phenotypeList=binaryVariableList
phenotypeList.remove('Is pregnant')
phenotypeList.remove('Home delivery')
phenotypeList.remove('Gender')


# for i,Ffile in enumerate(featureFiles):
# #     if n<2:
#     feature_dist_file='%s/%s' %(FeatureDistFolder,Ffile)
    
    
for i,fileNameF in enumerate(featureFiles):
#     if i<2:
        feature_dist_file='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/NewFeatures/%s' %fileNameF
        feature_name=fileNameF.replace('distMat_PNP434_','')    
        
        for j,phenotype in enumerate(phenotypeList):
#             if j<2:
                print i,j,feature_name,phenotype
                wait_on.append(q.method(compare_feature_distance_for_binary_phenotypes,kwargs={'feature_dist_file':feature_dist_file,'feature_name':feature_name, 'phenotype':phenotype,
                                                                                       'nPerm':nPerm,'removeSameUser':removeSameUser,'i':i,'j':j,
                                                                                       'filteringList':filteringList}))

       
q.wait(wait_on)
        
        