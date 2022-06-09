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
from queue.qp import qp,fakeqp
from addloglevels import sethandlers

MyPath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF'


datasetFolder='%s/TCR_real_data/SubSampled15000data_rep2' %MyPath
datasetName='PNP515_ss15000_rep2'
permNum=9999
minNshared=2
libPrepMethod=None
corrTest=True
FisherTest=True
mergeResults=False


FolderToSave='%s/TCR_mb_results/quickPermutDFs' %datasetFolder
if not isdir(FolderToSave):
    makedirs(FolderToSave)
sethandlers()
os.chdir(FolderToSave)
nPairs=279
minPair=0
maxPair=10


# send jobs to the cluster:
with qp(jobname='quickPerm', q=['himemint.q', 'himem7.q'], trds_def=1, mem_def='2G') as q:  # consider changing queue list!!
    q.startpermanentrun()
    tk = {}  # a dictionary storing arguments and results
    while minPair<nPairs:
        if maxPair>279:
            maxPair=279
        PairNumberRange=(minPair,maxPair)
        print minPair
        tk[minPair] = q.method(fast_perm_association_test_for_SeqSpeciesPairs,(datasetFolder,datasetName,permNum,minNshared,
                                        libPrepMethod,corrTest,FisherTest,mergeResults,PairNumberRange))
        minPair=minPair+10
        maxPair=maxPair+10
    tk = {k:q.waitforresult(v) for k, v in tk.iteritems()}
print 'done'


