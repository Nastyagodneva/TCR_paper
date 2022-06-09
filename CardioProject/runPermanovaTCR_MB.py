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

import os
from Utils import cacheOnDisk
from queue.qp import qp,fakeqp
from addloglevels import sethandlers

MyPath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF'


### define phenotypeDF:
datasetFolder='%s/TCR_real_data/SubSampled15000data_rep2' %MyPath
datasetName='PNP515_ss15000_rep2'
file2='%s/sharingAnalysis/sharingMatrix_PNP515_ss15000_rep2_minNshared2_RA_onlyProductiveTrue__percShared20_OLtrimmed_binary' %datasetFolder
TCRdf=pd.read_pickle(file2)
TCRdf.head()

TCRdfSorted=TCRdf.copy()
TCRdfSorted.loc['sum',:]=TCRdfSorted.sum()
TCRdfSorted=TCRdfSorted.sort_values(by='sum',axis=1,ascending=False)
TCRdfSorted.tail()


sequenceList=['CASSLAGSYEQYF','CASSLSGSSYNEQFF','CASSDRDTGELFF','CASSLGQGAYEQYF','CASSLYNEQFF','CASSLAGGTDTQYF','CASSLEETQYF',
              'CASSLGSSYEQYF','CASSPGETQYF','CASSELAGGQETQYF' ,'CSARLAGGQETQYF','CASSLGGQETQYF','CASSLDRNTEAFF']
print sequenceList
phenotypeDF=pd.DataFrame(TCRdfSorted.loc[:,sequenceList])
phenotypeDF=phenotypeDF.drop('sum',axis=0)
print phenotypeDF.shape


#define distMats to use:
folder='%s/MicrobiomeDataTables/distanceMatrices/forPermanova' %MyPath
distMats = [f for f in listdir(folder) if isfile(join(folder, f))]

nPerm=99999
removeSameUser=True
permanovaDFfolder='%s/TCR_mb_results/seqTCRpermanova/dfs2' %datasetFolder
if not isdir(permanovaDFfolder):
    makedirs(permanovaDFfolder)
    
sethandlers()
os.chdir('%s/TCR_mb_results/seqTCRpermanova' %datasetFolder)

with qp(jobname='permanova', q=['himemint.q', 'himem7.q'], trds_def=1, mem_def='2G') as q:  # consider changing queue list!!
    q.startpermanentrun()   
    count=0
    for n, distMat in enumerate(distMats):
        tk = {}  # a dictionary storing arguments and results
        if distMat!='MPA_g_standardParams_capped0_0001_percShared5_OLtrimmed_binary_jaccard_distMat':
            print  n, distMat
            file1='%s/%s' %(folder,distMat)
            feature_distMat=pd.read_pickle(file1)
            feature_name=distMat
            
            for k,phenotypeColumn in enumerate(phenotypeDF):
                print count,k,phenotypeColumn
                
                TCR_Mb_pair=(distMat,phenotypeColumn)         
                tk[TCR_Mb_pair] = q.method(permanova_withWrapper,(feature_distMat,feature_name,phenotypeDF,phenotypeColumn, 
                                             nPerm,removeSameUser,permanovaDFfolder))
                count=count+1
    tk = {k:q.waitforresult(v) for k, v in tk.iteritems()}
print 'done'


