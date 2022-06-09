from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
from scipy import stats
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot, draw_correlation_scatter
from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
from Bio.SeqUtils import GC
import seaborn as sns
import random
from scipy.stats import pearsonr
from queue.qp import qp, fakeqp
from addloglevels import sethandlers
import logging 
from Utils import cacheOnDisk
import os

from pop_organize import get_sample_data, get_sample_with_dfs
from SufficientStatistics import *
from MyFunctionsShani import *


basePath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/descriptiveStatsSamplesForAnalysis\GeneUsageCount'

@cacheOnDisk(basePath=basePath, filename='%(param)s_%(sample_name)s_variableCounts', force=True)
def variableCountForSample(sample_name, param):
    
    sample_df1 = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/SamplesForAnalysis/%s.tsv" % sample_name)
    countDF1=pd.DataFrame(sample_df1[param].value_counts(normalize=True))
    countDF1=countDF1.rename(columns={'%s' %param: '%s' %sample_name})
    return countDF1.T
   
    


 # send jobs to the queue (code adopted from Elad, see lab wiki for more details:
sethandlers()
os.chdir(basePath)
## the inputs for qp are: jobname, q=machine list, *** add max_r to prevent exploding the cluster!!***
with qp('recombEventVariableCount', q=['himem7.q', 'himemint.q'], mem_def="10G", trds_def=2, deleteCSHwithnoerr=True, tryrerun=False, max_u=120) as q:
    q.startpermanentrun()
    wait_on = []
    
# #now define a loop that divide the job and send each part seperatly:
    df_file_names,samples_with_df=get_sample_with_dfs()
    sample_list = samples_with_df
    param_list=['dFamilyName','vGeneName','n1Insertion','n2Insertion']
    count = 1
    
    for sample_name in sample_list:
            for param in param_list:
                print count, sample_name, param
                wait_on.append(q.method(variableCountForSample, kwargs={'sample_name':sample_name,'param':param}))
                count += 1
                            # #q.method takes the desired function with its arguments and send it to the queue.
    q.wait(wait_on) 
   
    
    