## this script includes functions to check predicted and observed frequencies of the while CDR3 combination
## based on 0.5:0.5 train:test set division of the data (of all 587 samples)

## this functions are based on the 'Sufficient Statistics 5 - find real sequence frequencies' jupyter notebook

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

#-----------------------------------------------------------------------------------


#internal functions:

#calculate the sum of logs of conditional probs for a specific combination:

def calc_cond_prob(indParVal, depParVal,condProbDF):
    
    CondProb=condProbDF.loc[indParVal,depParVal]
    logCondProb=np.log10(CondProb)
    return logCondProb

def calc_all_condProbs(seqParamsValDict,setType):
    indPar_list=['jGeneName','vGeneName','dFamilyName','dFamilyName','jGeneName']
    depPar_list=['dFamilyName','vDeletion','d5Deletion','d3Deletion','jDeletion']
    condProb_Train_df_list=[condProb_jGeneName_dFamilyName_Train_293,condProb_vGeneName_vDeletion_Train_293,
                           condProb_dFamilyName_d5Deletion_Train_293,condProb_dFamilyName_d3Deletion_Train_293,
                           condProb_jGeneName_jDeletion_Train_293]
    condProb_Test_df_list=[condProb_jGeneName_dFamilyName_Test_294,condProb_vGeneName_vDeletion_Test_294,
                           condProb_dFamilyName_d5Deletion_Test_294,condProb_dFamilyName_d3Deletion_Test_294,
                           condProb_jGeneName_jDeletion_Test_294]

    totalLogCondPro=0   
    
    for n in range(len(indPar_list)):
        indPar=indPar_list[n]
        depPar=depPar_list[n]
        #print indPar,depPar
        indParVal=seqParamsValDict[indPar]
        depParVal=seqParamsValDict[depPar]
        if setType=='Train':
            condProbDF=condProb_Train_df_list[n]
        else:
            condProbDF=condProb_Test_df_list[n]
        logCondPro=calc_cond_prob(indParVal, depParVal,condProbDF)
        totalLogCondPro=totalLogCondPro+logCondPro


    return totalLogCondPro
    

def load_all_train_probs_files():
    
    
    global vGeneTrainProbs_294,jGeneTrainProbs_294,condProb_jGeneName_dFamilyName_Train_293, \
            condProb_vGeneName_vDeletion_Train_293,condProb_jGeneName_jDeletion_Train_293, \
            condProb_dFamilyName_d5Deletion_Train_293,condProb_dFamilyName_d3Deletion_Train_293, \
            nt1FreqsDictDict_Train_ins1, lengthCount_Train_ins1, dinucNormDFDict_Train_ins1,nt1FreqsDictDict_Train_ins2,  \
            lengthCount_Train_ins2, dinucNormDFDict_Train_ins2
    
    
    #V:
    file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/vGeneTrainProbs_294'
    vGeneTrainProbs_294=pd.read_pickle(file1)

    #J:
    file2='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/jGeneTrainProbs_294'
    jGeneTrainProbs_294=pd.read_pickle(file2)

    #D|J:
    file3='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/condProb_jGeneName_dFamilyName_Train_293'
    condProb_jGeneName_dFamilyName_Train_293=pd.read_pickle(file3)

    #delV|V:
    file4='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/condProb_vGeneName_vDeletion_Train_293'
    condProb_vGeneName_vDeletion_Train_293=pd.read_pickle(file4)

    #delD5|D:
    file5='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/condProb_dFamilyName_d5Deletion_Train_293'
    condProb_dFamilyName_d5Deletion_Train_293=pd.read_pickle(file5)

    #delD3|D:
    file6='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/condProb_dFamilyName_d3Deletion_Train_293'
    condProb_dFamilyName_d3Deletion_Train_293=pd.read_pickle(file6)

    #delJ|J:
    file7='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/condProb_jGeneName_jDeletion_Train_293'
    condProb_jGeneName_jDeletion_Train_293=pd.read_pickle(file7)
    
    #ins1, ins2 parameters:

    file11='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/nt1FreqsDictDict_train_ins1'
    nt1FreqsDictDict_Train_ins1=pd.read_pickle(file11)

    file12='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/nt1FreqsDictDict_train_ins2'
    nt1FreqsDictDict_Train_ins2=pd.read_pickle(file12)

    file13='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/lengthCount_train_ins1'
    lengthCount_Train_ins1=pd.read_pickle(file13)

    file14='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/lengthCount_train_ins2'
    lengthCount_Train_ins2=pd.read_pickle(file14)

    file15='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/dinucNormDFDict_train_ins1'
    dinucNormDFDict_Train_ins1=pd.read_pickle(file15)

    file16='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/dinucNormDFDict_train_ins2'
    dinucNormDFDict_Train_ins2=pd.read_pickle(file16)
    
def load_all_test_probs_files():
    
    global vGeneTestProbs_293,jGeneTestProbs_293,condProb_jGeneName_dFamilyName_Test_294, \
    condProb_vGeneName_vDeletion_Test_294,     condProb_jGeneName_jDeletion_Test_294, \
    condProb_dFamilyName_d5Deletion_Test_294,condProb_dFamilyName_d3Deletion_Test_294, \
    nt1FreqsDictDict_Test_ins1, lengthCount_Test_ins1, dinucNormDFDict_Test_ins1,nt1FreqsDictDict_Test_ins2, \
    lengthCount_Test_ins2,     dinucNormDFDict_Test_ins2
    
    #V:
    file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/vGeneTestProbs_293'
    vGeneTestProbs_293=pd.read_pickle(file1)

    #J:
    file2='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/jGeneTestProbs_293'
    jGeneTestProbs_293=pd.read_pickle(file2)

    #D|J:
    file3='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/condProb_jGeneName_dFamilyName_Test_294'
    condProb_jGeneName_dFamilyName_Test_294=pd.read_pickle(file3)

    #delV|V:
    file4='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/condProb_vGeneName_vDeletion_Test_294'
    condProb_vGeneName_vDeletion_Test_294=pd.read_pickle(file4)

    #delD5|D:
    file5='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/condProb_dFamilyName_d5Deletion_Test_294'
    condProb_dFamilyName_d5Deletion_Test_294=pd.read_pickle(file5)

    #delD3|D:
    file6='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/condProb_dFamilyName_d3Deletion_Test_294'
    condProb_dFamilyName_d3Deletion_Test_294=pd.read_pickle(file6)

    #delJ|J:
    file7='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/condProb_jGeneName_jDeletion_Test_294'
    condProb_jGeneName_jDeletion_Test_294=pd.read_pickle(file7)
    
    #ins1, ins2 parameters:

    file11='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/nt1FreqsDictDict_test_ins1'
    nt1FreqsDictDict_Test_ins1=pd.read_pickle(file11)

    file12='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/nt1FreqsDictDict_test_ins2'
    nt1FreqsDictDict_Test_ins2=pd.read_pickle(file12)

    file13='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/lengthCount_test_ins1'
    lengthCount_Test_ins1=pd.read_pickle(file13)

    file14='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/lengthCount_test_ins2'
    lengthCount_Test_ins2=pd.read_pickle(file14)

    file15='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/dinucNormDFDict_test_ins1'
    dinucNormDFDict_Test_ins1=pd.read_pickle(file15)

    file16='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/dinucNormDFDict_test_ins2'
    dinucNormDFDict_Test_ins2=pd.read_pickle(file16)

# a function calculating the probability of a specific combination based on either the train or test data set:
#when generalizing this function. need to take the read_pickles out

def calc_logProb_per_seq(seqParamsValDict,setType):
    totalLogProb=0
       
    load_all_train_probs_files()
    load_all_test_probs_files()
       
    
    #calculate sum of logs of conditional probs and add to the total log prob:
    totalLogCondPro=calc_all_condProbs(seqParamsValDict,setType)
    totalLogProb=+totalLogCondPro
    
    #calc log probs for v and j and add to the total log prob:
    if setType=='Train':
        vGeneProbsDF=vGeneTrainProbs_294
        jGeneProbsDF=jGeneTrainProbs_294
    else: 
        vGeneProbsDF=vGeneTestProbs_293
        jGeneProbsDF=jGeneTestProbs_293          
       
    vGeneName=seqParamsValDict['vGeneName']
    vProb=vGeneProbsDF[vGeneName]
    #print 'vprob=%s' %vProb
    logVProb=np.log10(vProb)
    
    jGeneName=seqParamsValDict['jGeneName']
    jProb=jGeneProbsDF[jGeneName]
    #print 'jprob=%s' %jProb
    logJProb=np.log10(jProb)
   
    totalLogProb=totalLogProb+logVProb+logJProb
    
    #calc log probs for ins1 and ins2 and add to the total log prob:
    if setType=='Train':
        n1InsSeq=seqParamsValDict['n1InsSeq']
        logProbins1Seq=calc_insSeq_prob(n1InsSeq, 20, nt1FreqsDictDict_Train_ins1, dinucNormDFDict_Train_ins1, lengthCount_Train_ins1)
        n2InsSeq=seqParamsValDict['n2InsSeq']
        logProbins2Seq=calc_insSeq_prob(n2InsSeq, 20, nt1FreqsDictDict_Train_ins2, dinucNormDFDict_Train_ins2, lengthCount_Train_ins2)
    else:
        n1InsSeq=seqParamsValDict['n1InsSeq']
        logProbins1Seq=calc_insSeq_prob(n1InsSeq, 20, nt1FreqsDictDict_Test_ins1, dinucNormDFDict_Test_ins1, lengthCount_Test_ins1)
        n2InsSeq=seqParamsValDict['n2InsSeq']
        logProbins2Seq=calc_insSeq_prob(n2InsSeq, 20, nt1FreqsDictDict_Test_ins2, dinucNormDFDict_Test_ins2, lengthCount_Test_ins2)
        
        
        
    totalLogProb=totalLogProb+logProbins1Seq+logProbins2Seq
    
    #print 'ins1 prob=%s' %10**logProbins1Seq
    #print 'ins2 prob=%s' %10**logProbins2Seq
    
    return totalLogProb
    
    #-------------------------------------------------------------------------------------------------


basePath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred' ## define the path to 

@cacheOnDisk(basePath=basePath, filename='compare_cdr3_prediction_%(min_comb)s_%(max_comb)s', force=True)
def compare_cdr3_prediction(min_comb,max_comb):

    file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/WholeCDR3SeqPred/params_summary_tables/params_summary_TestSet_grouped_replicates'
    df=pd. read_pickle(file1)
        
    n_combinations_testSet=8352567
    n_combs_inDF=len(df.index)
    
    print min_comb, max_comb
    if max_comb>n_combs_inDF:
        max_comb=n_combs_inDF
    
    
    n_list=[]
    totalLogProb_Train_list=[]
    totalLogProb_Test_list=[]
    logObsFreqInTestSet_list=[]
   
    for n in range(min_comb,max_comb):
        print n
        seqParamsVals=df.index[n]
        paramsNames=df.index.names
        seqParamsValDict=dict(zip(paramsNames, seqParamsVals))
        
        totalLogProb_Test=calc_logProb_per_seq(seqParamsValDict,'Test')
        totalLogProb_Train=calc_logProb_per_seq(seqParamsValDict,'Train')
        obsfreq=float(df['count'][n])/n_combinations_testSet
        logObsFreqInTestSet=np.log10(obsfreq)
                
        
        #exp_n_comb_appears_in_sample=n_combinations_sample*10**totalLogProb
        #exp_n_comb_appears_in_testSet=n_combinations_testSet*10**totalLogProb
        #real_n_in_testSet=df['count'][n]
        
        
        n_list.append(n)
        totalLogProb_Train_list.append(totalLogProb_Train)
        totalLogProb_Test_list.append(totalLogProb_Test)
        logObsFreqInTestSet_list.append(logObsFreqInTestSet)

    cdr3_prediction_compare_df=pd.DataFrame({'n':n_list, 'cdr3LogProb_byTrainSet':totalLogProb_Train_list,
                                        'cdr3LogProb_byTestSet':totalLogProb_Test_list,
                                        'cdr3LogObsFreq': logObsFreqInTestSet_list})
        
        
        
    return cdr3_prediction_compare_df


## send jobs to the queue (code adopted from Elad, see lab wiki for more details:
sethandlers()
os.chdir(basePath)
## the inputs for qp are: jobname, q=machine list, *** add max_r to prevent exploding the cluster!!***
with qp('compare_cdr3_prediction_job',  q = ['himem7.q','himemint.q'], mem_def="10G", trds_def=2, deleteCSHwithnoerr = True, tryrerun = False, max_u=120) as q:
    q.startpermanentrun()
    wait_on =[]
    
##now define a loop that divide the job and send each part seperatly:
    min_comb=0
    max_comb=2000
    while min_comb<171100:                                     
        print min_comb
        wait_on.append(q.method(compare_cdr3_prediction,kwargs={'min_comb':min_comb,'max_comb':max_comb}))
            ##q.method takes the desired function with its arguments and send it to the queue.
        min_comb=min_comb+2000
        max_comb=max_comb+2000
    q.wait(wait_on)






