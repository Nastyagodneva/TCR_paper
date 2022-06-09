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




basePath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/Ins2PredictionAllSamples' ## define the path to 
n_samples=587
generate_dfs=False

@cacheOnDisk(basePath=basePath, filename='calc_predictionPower_allSamples_ins2_%(min_sample)s_%(max_sample)s', force=True)
def calc_predictionPower_allSamples(min_sample, max_sample):
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
    sampleLengthList=[]
    TrainFracsList=[]
    rExpToObsCorrelList=[]
    pExpToObsCorrelList=[]
    #figExpToObsCorrelList=[]
    #axExpToObsCorrelList=[]       
    
    #fig1, ((ax1, ax3,ax5,ax7),(ax2,ax4,ax6,ax8))= plt.subplots(nrows=2,ncols=4,figsize=(13,5),sharex=True,sharey=True)
    #fig1.suptitle('Ins1 sequence frequency - correlation between expected and observed values', fontsize=18)
    #axList=[ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8]
    #count=0
    
    
    for sample_name in sample_list:
        for train_fraction in train_fraction_list:
            print sample_name
            #print train_fraction
            axB=None
        figExpToObsCorrel,axExpToObsCorrel,rExpToObsCorrel,pExpToObsCorrel,sample_length=calc_plot_save_modelParamPrediction(
            sample_name,2, threshLength,train_fraction,axB)
        samplesList.append(sample_name)
        sampleLengthList.append(sample_length)
        TrainFracsList.append(train_fraction)
        rExpToObsCorrelList.append(rExpToObsCorrel)
        pExpToObsCorrelList.append(pExpToObsCorrel)
        #figExpToObsCorrelList.append(figExpToObsCorrel)
        #axExpToObsCorrelList.append(axExpToObsCorrel)
        #count+=1

    inspredictionStrengthDF=pd.DataFrame({'Sample': samplesList, '# nonProd unique seqeunces':sampleLengthList,
                                      'train set fraction':TrainFracsList, 
                                       'expected to observed frequency correlation r':rExpToObsCorrelList,
                                      'expected to observed frequency correlation p':pExpToObsCorrelList})
    
    return inspredictionStrengthDF

    
## send jobs to the queue (code adopted from Elad, see lab wiki for more details:
sethandlers()
os.chdir(basePath)
## the inputs for qp are: jobname, q=machine list, *** add max_r to prevent exploding the cluster!!***
with qp('CalcPredictionIns2_job',  q = ['himemint.q','himem7.q'], mem_def="10G", trds_def=2, deleteCSHwithnoerr = True, tryrerun = False, max_u=120) as q:
    q.startpermanentrun()
    wait_on =[]
    
##now define a loop that divide the job and send each part seperatly:
    min_sample=380
    max_sample=400
    while min_sample<381:                                     
        print min_sample
        wait_on.append(q.method(calc_predictionPower_allSamples,kwargs={'min_sample':min_sample,'max_sample':max_sample}))
            ##q.method takes the desired function with its arguments and send it to the queue.
        min_sample=min_sample+20
        max_sample=max_sample+20
    q.wait(wait_on)
