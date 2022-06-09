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


basePath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/Ins1PredictionAllSamples'

@cacheOnDisk(basePath=basePath, filename='ins1predOpt_%(sample_name)s_%(basicUnit)s_%(addPrior)s_%(insType)s_%(capping)s_%(Trim)s', force=True)
def ins1Prediction_AllSample_optimization(sample_name, basicUnit, insType, addPrior, capping,Trim):
    # load ins1 lists for seq and comb cases - *** must use read_pickle and not read_excel, as the file is too long for excel,
    # and also the index is reset while saving to excel
    
    print 'loading ins table...'
    
    if insType == 1 and basicUnit == 'comb':
        file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins1_dfs/AllinsSeqsDFsmall_ins1_comb'
        AllnsSeqsDFsmall = pd.read_pickle(file1)
    elif insType == 1 and basicUnit == 'seq':
        file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins1_dfs/Alln1SeqsDFsmall'
        AllnsSeqsDFsmall = pd.read_pickle(file1)
    elif insType == 2 and basicUnit == 'comb':
        file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins2_dfs/AllinsSeqsDFsmall_ins2_comb'  
        AllnsSeqsDFsmall = pd.read_pickle(file1)
    elif insType == 2 and basicUnit == 'seq':
        file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins2_dfs/AllinsSeqsDFsmall_ins2_seq'
        AllnsSeqsDFsmall = pd.read_pickle(file1)
        
    axB = None
    threshLength = 20
  
    print 'checking if expObsFreqDFopt already exist for the sample...'
    filename='expObsFreqDFopt_%s_%s_%s_%s_%s' % (sample_name, insType, basicUnit, addPrior, capping)
    dfs_folder='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins_opt'
    filenames = [f for f in listdir(dfs_folder) if isfile(join(dfs_folder, f))]
    
    print 'calculating predictions and correlations...'
    if filename in filenames:
        rExpToObsCorrel, pExpToObsCorrel, perc_predicted_to_appear = predictInsoutSampleBasedOnAllOthers_whenExist(sample_name, AllnsSeqsDFsmall, basicUnit, insType, axB, addPrior,capping, Trim)
    else:
        rExpToObsCorrel, pExpToObsCorrel, perc_predicted_to_appear = predictInsoutSampleBasedOnAllOthers(sample_name, AllnsSeqsDFsmall, basicUnit, insType, axB, addPrior,capping, Trim)
   
    
    ins1predictionStrengthDF = pd.DataFrame()
    ins1predictionStrengthDF.loc[1,'Sample']= sample_name
    ins1predictionStrengthDF.loc[1,'basicUnit']= basicUnit
    ins1predictionStrengthDF.loc[1,'addPrior']= addPrior
    ins1predictionStrengthDF.loc[1,'insType']= insType
    ins1predictionStrengthDF.loc[1,'capping']= capping
    ins1predictionStrengthDF.loc[1,'Trim']= Trim
    ins1predictionStrengthDF.loc[1,'expected to observed frequency correlation r']= rExpToObsCorrel
    ins1predictionStrengthDF.loc[1,'expected to observed frequency correlation p']= pExpToObsCorrel
    ins1predictionStrengthDF.loc[1,'perc_predicted_to_appear']= perc_predicted_to_appear
    
                    
                
    return ins1predictionStrengthDF
        
      
    
    
# # send jobs to the queue (code adopted from Elad, see lab wiki for more details:
sethandlers()
os.chdir(basePath)
## the inputs for qp are: jobname, q=machine list, *** add max_r to prevent exploding the cluster!!***
with qp('ins1OptJob', q=['himem7.q', 'himemint.q'], mem_def="10G", trds_def=2, deleteCSHwithnoerr=True, tryrerun=False, max_u=120) as q:
    q.startpermanentrun()
    wait_on = []
    
# #now define a loop that divide the job and send each part seperatly:
    insType_List = [1, 2]
    sample_list = ['HIP13505', 'HIP14071', 'HIP13518', 'HIP01091']
    basicUnit_list = ['seq', 'comb']
    prior_list = ['weak constant', 'weak differ', 'strong constant', 'strong differ']
    capping_list = [None, 'weak', 'strong']
    Trim_list=['without','drop','trimming']
    
    count = 1
    for insType in insType_List:
        for sample_name in sample_list:
            for basicUnit in basicUnit_list:
                for addPrior in prior_list:
                    for capping in capping_list:
                        for Trim in Trim_list:
                            print count, insType, sample_name, basicUnit, addPrior,capping,Trim
                            wait_on.append(q.method(ins1Prediction_AllSample_optimization, kwargs={'sample_name':sample_name,
                                                                                              'basicUnit':basicUnit,
                                                                                              'insType':insType,
                                                                                              'addPrior':addPrior,
                                                                                              'capping':capping,
                                                                                              'Trim':Trim}))
                            count += 1
                            # #q.method takes the desired function with its arguments and send it to the queue.
    q.wait(wait_on) 
