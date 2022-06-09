from os import listdir,mkdir,makedirs
from os.path import isfile, join, isdir,exists
import pandas as pd
import numpy as np
from scipy import stats
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from ShaniBA.myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot,draw_correlation_scatter
from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
from Bio.SeqUtils import GC
import seaborn as sns
import random
from scipy.stats import pearsonr,fisher_exact
from skbio.diversity.alpha import shannon, simpson, berger_parker_d

from ShaniBA.pop_organize import get_sample_data, get_sample_with_dfs
# from SufficientStatistics import *
from ShaniBA.MyFunctionsShani import *
import math
from skbio.stats.distance import mantel
from scipy.spatial.distance import braycurtis, pdist
# from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import * 
from ShaniBA.TCR_microbiome_interactions.TCR_microbiome_interactions_functions import *
from ShaniBA.SampleLists.SampleFileFunctions import *
from ShaniBA.CardioProject.CardioFunctions import *


import os
from Utils import cacheOnDisk
from SegalQueue.qp import qp,fakeqp
from addloglevels import sethandlers

MyPath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF'

datasetFolder='%s/TCR_real_data/PNP530Cardio126Combined' %MyPath
datasetName='PNP530Cardio126'

f1='%s/sharingAnalysis/sharingMatrix_PNP530Cardio126_minNshared5_RA_onlyProductiveTrue__minSharedT5_OLtrimmed_binary' %datasetFolder
combTCRdfmin5=pd.read_pickle(f1)
print combTCRdfmin5.shape

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
Xname=phenotypeDFname
phenotypeDF=PNP530Cardio126_phen[Xcols].set_index('BD')
TCRvariableDF=combTCRdfmin5
TCRvariableName='combTCRdfmin5'
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

filenames=[f for f in listdir(resultFolder) if isfile(join(resultFolder, f))]
print filenames
    
    
TCRdfMinCol=0
TCRdfMaxCol=500
maximal=len(combTCRdfmin5.columns)
# maximal=1500
maxiter=1000

sethandlers()
os.chdir(resultFolder)

# send jobs to the cluster:
with qp(jobname='logReg', q=['himemint.q', 'himem7.q'], trds_def=1, mem_def='2G') as q:  # consider changing SegalQueue list!!
    q.startpermanentrun()
    tk = {}  # a dictionary storing arguments and results
    while TCRdfMinCol<maximal:
        if TCRdfMaxCol>maximal:
            TCRdfMaxCol=maximal
        resultFileName = 'resultDF_regType%s_%s_%s_ycol%s_ycol%s_mi%s.xlsx'   % (modelType, phenotypeDFname, TCRvariableName, TCRdfMinCol, TCRdfMaxCol, maxiter)
        print resultFileName
        print TCRdfMinCol, TCRdfMaxCol, maximal 
        if resultFileName not in filenames:
            print TCRdfMinCol
            tk[TCRdfMinCol] = q.method(predictBinaryPhenotype_logReg_TCRfeatures,(phenotypeDF, phenotypeDFname, TCRvariableDF, TCRvariableName, targetDF, 
                                               targetName, resultFolder, modelType, TCRdfMinCol, TCRdfMaxCol, 
                                               maxiter))
            print 'blabla'
            
        else:
            print 'this resultDF already exist'
#             tk[TCRdfMinCol]=pd.DataFrame()
        TCRdfMinCol=TCRdfMinCol+500
        TCRdfMaxCol=TCRdfMaxCol+500
    tk = {k:q.waitforresult(v) for k, v in tk.iteritems()}

print 'done'


