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
# from SufficientStatistics.SufficientStatistics import *
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
print twoCohortsDF.head()

f1='%s/TCR_real_data/CardioSamples/PNP_Cardio_seq_comparison/FisherSeqComparison_9999PERMUTATIONS_PNP530_Cardio126.xlsx' %MyPath
FisherPerms9999=pd.read_excel(f1)

seqListTo126700=FisherPerms9999[FisherPerms9999['permP']<0.0002].index.tolist()
print 'n sequences to 9999 perms=%s' %len(seqListTo126700)
columns=seqListTo126700+['cohort']

twoCohortsDFto126700=twoCohortsDF[columns]


nPerm=126700

folderToSave='%s/TCR_real_data/CardioSamples/PNP_Cardio_seq_comparison/fisherSeqComparisonDFs_9999' %MyPath
if not isdir(folderToSave):
    makedirs(folderToSave)
sethandlers()
os.chdir(folderToSave)

minSeq=0
maxSeq=minSeq+1
maximal=len(twoCohortsDFto126700.columns)-1



print minSeq, maxSeq, maximal
 
 
 # send jobs to the cluster:
with qp(jobname='quickPerm', q=['himemint.q', 'himem7.q'], trds_def=1, mem_def='2G') as q:  # consider changing queue list!!
     q.startpermanentrun()
     tk = {}  # a dictionary storing arguments and results
     while minSeq<maximal:
         if maxSeq>maximal:
             maxSeq=maximal
         print minSeq
         tk[minSeq] = q.method(calc_permFisher_twoCohorts,
                               (minSeq,maxSeq,twoCohortsDFto126700,datasetName1,datasetName2,nPerm,folderToSave))
         print 'blabla'
         minSeq=minSeq+1
         maxSeq=maxSeq+1
     tk = {k:q.waitforresult(v) for k, v in tk.iteritems()}
print 'done'


