## necessary imports:
from os import listdir
from os.path import isfile, join
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
from pop_organize import get_sample_data, get_sample_with_dfs
from queue.qp import qp,fakeqp
from addloglevels import sethandlers
import logging 
from Utils import cacheOnDisk
import os

from SufficientStatistics import *

basePath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/Ins1PredictionAllSamples' ## define the path to 
n_samples=587
generate_dfs=False

@cacheOnDisk(basePath=basePath, filename='calc_predictionPower_sampleByAllOthers_%(min_sample)s_%(max_sample)s', force=True)
def ins1Prediction_AllSample_whenExist(min_sample, max_sample):
    pickleFile='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins1_dfs/Alln1SeqsDFsmall'
    Allns1SeqsDFsmall=pd.read_pickle(pickleFile)
    
    import time
    cdate=str(time.strftime("%d%m%Y"))
    cdate
    
    print min_sample, max_sample
    if max_sample>587:
        max_sample=587
    n=1
    df_file_names,samples_with_df=get_sample_with_dfs()
    
    threshLength=20
    sample_list=samples_with_df[min_sample:max_sample]
    train_fraction_list=[0.8]
      
    samplesList=[]
    rList=[]
    pList=[]
       
    
    for sample_name in sample_list:
        r,p=predictInsoutSampleBasedOnAllOthers_whenExist(sample_name,Allns1SeqsDFsmall,1)
        samplesList.append(sample_name)
        rList.append(r)
        pList.append(p)
        
        
    ins1predictionStrengthDF_AllSamples=pd.DataFrame({'Sample': samplesList,
                                       'expected to observed frequency correlation r':rList,
                                      'expected to observed frequency correlation p':pList})
    
    return ins1predictionStrengthDF_AllSamples



## send jobs to the queue (code adopted from Elad, see lab wiki for more details:
sethandlers()
os.chdir(basePath)
## the inputs for qp are: jobname, q=machine list, *** add max_r to prevent exploding the cluster!!***
with qp('CalcPredictionIns1_job',  q = ['himem7.q','himemint.q'], mem_def="10G", trds_def=2, deleteCSHwithnoerr = True, tryrerun = False, max_u=120) as q:
    q.startpermanentrun()
    wait_on =[]
    
##now define a loop that divide the job and send each part seperatly:
    min_sample=0
    max_sample=360
    while min_sample<361:                                     
        print min_sample
        wait_on.append(q.method(ins1Prediction_AllSample_whenExist,kwargs={'min_sample':min_sample,'max_sample':max_sample}))
            ##q.method takes the desired function with its arguments and send it to the queue.
        min_sample=min_sample+3
        max_sample=max_sample+3
    q.wait(wait_on)