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


def serial_merging(filenames,DFfolder, toSaveFolder,df_min,df_max,columns_to_merge_list,columns_to_rename):
    filenames_partial=filenames[df_min:df_max]
    print 'total number of dfs to merge is %s' %len(filenames_partial)
    f1=filenames_partial[0]
    file1='%s/%s' %(DFfolder,f1)
    df1=pd.read_pickle(file1)
    merged=df1
    nPer=f1.split('DF')[0]
    for col in columns_to_rename:
        merged=merged.rename(columns={col:col+'_'+str(nPer)})
    
    for n,f in enumerate(filenames_partial[1:]):
        print n,f
        file1='%s/%s' %(DFfolder,f)
        df=pd.read_pickle(file1)
        merged=pd.merge(merged,df,how='inner',left_on=columns_to_merge_list,
                        right_on=columns_to_merge_list)
        nPer=f.split('DF')[0]
        for col in columns_to_rename:
            merged=merged.rename(columns={col:col+'_'+str(nPer)})
        
        mergedFile='%s/merged_%s_%s' %(toSaveFolder,df_min,df_max)
        merged.to_pickle(mergedFile)
    return merged


datasetFolder='%s/TCR_real_data/SubSampled15000data_rep2' %MyPath
NsharedSamplesForSpecies=None
NsharedSamplesForSeqs=None
topNspecies=250
topNseqs=250


permFisherTestFolder='%s/TCR_mb_results/permFisherTest_%s%s%s%s' %(datasetFolder,NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies, topNseqs)
permFisherDFsFolder= '%s/permDFs' %permFisherTestFolder
DFfolder=permFisherDFsFolder
toSaveFolder='%s/mergedDFs' %permFisherTestFolder
if not isdir(toSaveFolder):
    makedirs(toSaveFolder)
filenames = [f for f in listdir(permFisherDFsFolder) if isfile(join(permFisherDFsFolder, f))]
nPermUsed=len(filenames)
columns_to_merge_list=['species','seq']
columns_to_rename=['OR','log10OR']

sethandlers()
os.chdir(toSaveFolder)

#send jobs to the cluster:
df_min=0
df_max=100

with qp(jobname='FisherDiv', q = ['himemint.q','himem7.q'],trds_def=1,mem_def='2G') as q: #consider changing queue list!!
    q.startpermanentrun()
    tk = {} #a dictionary storing arguments and results
    while df_min<9999:
        print df_min
        tk[df_min]=q.method(serial_merging, (filenames,DFfolder, toSaveFolder,df_min,df_max,columns_to_merge_list,columns_to_rename))
        df_min=df_min+100
        df_max=df_max+100
    tk = {k:q.waitforresult(v) for k,v in tk.iteritems()}
print 'done'

