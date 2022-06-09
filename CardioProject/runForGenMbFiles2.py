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


#load PNP515 sample list:
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP515','rb') as fp:
    PNP515=pickle.load(fp)
print len(PNP515)
print PNP515[:5]


datasetFolder='%s/TCR_real_data/SubSampled15000data_rep1' %MyPath
datasetName='PNP515_ss15000_rep1'

gen_mb=True
mbLevel='s'
mbDataFolder='AllSeqProjects'
SampleList=PNP515
SampleListName='PNP515'
filterGenotek=False
filterMinimalReads=False
filterlibPrepMethod=None

# (NsharedSamplesForSpecies,NsharedSamplesForSeqs,topNspecies,topNseqs):
thresholdList=[(None,None,250,250)]
outlierSTD=None
corrTest='spearman'
completecorrelation=True
correlationForHeatMap=False

usePermCorr=None
runFisher=None
runCorr=None
usePermFisher=None
nPermCorr=None
nPermFisher=None
firstPermCorr=None
firstPermFisher=None
runRealCorr=None
runRealFisher=None
mergeResults=None


merged_fisher_correlation_withIdentity=TCR_Mb_interactions_for_dataset(datasetFolder, datasetName, gen_mb, mbLevel, mbDataFolder, SampleList,
                                    SampleListName, filterGenotek, filterMinimalReads,filterlibPrepMethod,
                                    thresholdList, outlierSTD, corrTest, completecorrelation,
                                    correlationForHeatMap,usePermCorr,runCorr,runFisher,usePermFisher,nPermCorr,nPermFisher,
                                    firstPermCorr,firstPermFisher,runRealCorr,runRealFisher,mergeResults)
       

