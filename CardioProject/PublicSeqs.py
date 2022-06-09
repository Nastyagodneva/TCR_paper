from os import listdir
from os.path import isfile, join
# from Utils import Load, Write
import pandas as pd
import numpy as np
#from scipy import stats
# import math
#import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib import cm
# import plotly.plotly as py
from myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot
#from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
#from Bio.SeqUtils import GC
from collections import Counter
from pop_organize import get_sample_data, get_sample_with_dfs
from queue.qp import qp,fakeqp
from addloglevels import sethandlers
import logging 
from Utils import cacheOnDisk
import os

#----------------------------------------------------------------------------
'''
def count_aa(df,total_unique_aa, sample_name):
    print('adding unique aas from sample %s...' %sample_name)
    list_aa=list(df['aminoAcid'])
    list_aa_new=[i for i in list_aa if isinstance(i, str)]
    sample_unique_aa=list(set(list_aa_new))
    print sample_unique_aa[:5]
    print total_unique_aa[:5]
    total_unique_aa=total_unique_aa+sample_unique_aa
    print total_unique_aa[:5]
    print len(total_unique_aa)
    return total_unique_aa

#-----------------------------------------------------------------------------


def aa_counter(total_unique_aa):
    print'counting aa sequence occurence number...'
    aa_seq_counter_dic=Counter(total_unique_aa)
    print'finished counting aa sequence occurence number...'
    return aa_seq_counter_dic 


#-----------------------------------------------------------------------------
################################
#DEFINE SOME PREFERENCES HERE: #
################################

n_samples=573
generate_unique_aa=False ## define if new list of unique aa's are needed
generate_counters=False ## define whether to re-generate counters (recommended if the lists of unique aa's were re-generated)
generate_grouped_counter_dfs=False
show_counter_dfs=True
plot_n_samples=False ##plotting should be defined more carefully (subplots, pdfs...)
generate_dfs=False
df_file_names,samples_with_df=get_sample_with_dfs() ## get list of samples with df's


## if new lists if unique aa's are needed, they are generated with a loop *** consider dividing the task
## to batches of 10 samples and sent to the queue in parallel.
if generate_unique_aa:
    total_unique_aa_prod=[]
    total_unique_aa_non_prod=[]
    for sample in samples_with_df[:n_samples]:
        sample_df, sample_df_prod, sample_df_non_prod=get_sample_data(sample, generate_dfs)
        total_unique_aa_prod=count_aa(sample_df_prod,total_unique_aa_prod,sample)
        total_unique_aa_non_prod=count_aa(sample_df_non_prod,total_unique_aa_non_prod,sample)
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/total_unique_aa_prod-n_samples-%s' %n_samples,"wb") as f1:
        pickle.dump(total_unique_aa_prod, f1)
    f1. close()
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/total_unique_aa_non_prod-n_samples-%s' %n_samples,"wb") as f2:
        pickle.dump(total_unique_aa_non_prod, f2)
    f2. close()

## if no new unique_aa's lists are needed, load the existing lists from pickles:
else:    
    print 'loading unique_aa lists from pickles...'
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/total_unique_aa_prod-n_samples-%s' %n_samples, "rb") as f1:
        total_unique_aa_prod=pickle.load(f1)
    f1. close()
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/total_unique_aa_non_prod-n_samples-%s' %n_samples,"rb") as f2:
        total_unique_aa_non_prod=pickle.load(f2)
    f2. close()
    print 'Finished loading unique_aa lists from pickles'


if generate_counters:
    ## if new counters are needed, they are calculated here *** consider dividing the task
    ## to batches- non trivial
    print 'generating new counters...'
    prod_aa_seq_counter_dic=aa_counter(total_unique_aa_prod)
    prod_aa_seq_counter_df=pd.DataFrame(prod_aa_seq_counter_dic.items(), columns=['aa_sequence', 'n_samples'])
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/prod_aa_seq_counter_df-n_samples-%s' %n_samples, "wb") as f1:
        pickle.dump(prod_aa_seq_counter_df, f1)
    f1.close()
    non_prod_aa_seq_counter_dic=aa_counter(total_unique_aa_non_prod)
    non_prod_aa_seq_counter_df=pd.DataFrame(non_prod_aa_seq_counter_dic.items(), columns=['aa_sequence', 'n_samples'])
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/non_prod_aa_seq_counter_df-n_samples-%s' %n_samples, "wb") as f2:
        pickle.dump(non_prod_aa_seq_counter_df, f2)
    f2.close()
    print 'Finished generating new counters...'
else:
    ## if counters exist, they are loaded here:
    print 'loading counters...'
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/prod_aa_seq_counter_df-n_samples-%s' %n_samples, "rb") as f1:
        prod_aa_seq_counter_df=pickle.load(f1)
    f1.close()
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/non_prod_aa_seq_counter_df-n_samples-%s' %n_samples, "rb") as f2:
        non_prod_aa_seq_counter_df=pickle.load(f2)
    f2.close()
    print 'fINISHED loading counters...'
    


if generate_grouped_counter_dfs:
    print 'generating grouped counter dfs...'
    prod_aa_seq_counter_grouped=prod_aa_seq_counter_df.groupby(['n_samples']).count() 
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/prod_aa_seq_counter_grouped-n_samples-%s' %n_samples, "wb") as f1:
        pickle.dump(prod_aa_seq_counter_grouped, f1)
    f1.close()
    non_prod_aa_seq_counter_grouped=non_prod_aa_seq_counter_df.groupby(['n_samples']).count() 
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/non_prod_aa_seq_counter_grouped-n_samples-%s' %n_samples, "wb") as f2:
        pickle.dump(non_prod_aa_seq_counter_grouped, f2)
    f2.close()
    print 'finished generating grouped counter dfs...'
    
    
## showing the grouped caller dfs by number of samples per sequence (only top and bottom 20 sequences are shown:   
if show_counter_dfs:
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/prod_aa_seq_counter_grouped-n_samples-%s' %n_samples, "rb") as f1:
        prod_aa_seq_counter_grouped=pickle.load(f1)
    f1.close()
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/non_prod_aa_seq_counter_grouped-n_samples-%s' %n_samples, "rb") as f2:
        non_prod_aa_seq_counter_grouped=pickle.load(f2)
    f2.close()
       
    print 'prod counter df grouped:'
    print prod_aa_seq_counter_grouped[:20]
    print prod_aa_seq_counter_grouped[-20:]
    
    print 'non_prod counter df grouped:'
    print non_prod_aa_seq_counter_grouped[:20]
    print non_prod_aa_seq_counter_grouped[-20:]
    

##plotting should be defined more carefully (subplots, pdfs...)    
if plot_n_samples:
    prod_n_samples=list(prod_aa_seq_counter_grouped.index)
    prod_public_counter=list(prod_aa_seq_counter_grouped['aa_sequence'])
    plt.bar(prod_n_samples,prod_public_counter)
    plt.yscale('log')
    plt.ylim(0.5,10**8)
    plt.xlabel('Number of samples')
    plt.ylabel('Number of Sequences')
    plt.title('Shared samples per aa sequence distribution')
    plt.show()
    non_prod_n_samples=list(non_prod_aa_seq_counter_grouped.index)
    non_prod_public_counter=list(non_prod_aa_seq_counter_grouped['aa_sequence'])
    plt.bar(non_prod_n_samples,non_prod_public_counter)
    plt.yscale('log')
    plt.ylim(0.5,10**8)
    plt.xlabel('Number of samples')
    plt.ylabel('Number of Sequences')
    plt.title('Shared samples per aa sequence distribution')
    plt.show()
    
 '''   
#-----------------------------------------------------------------------------------------------
def add_public(df,counter_df):
    sample_df_wityPublic=pd.merge(df, counter_df, how='inner', left_on='aminoAcid', right_index=True,
         suffixes=('_x', '_y'), copy=True, indicator=True)
    sample_df_wityPublic['shareStatus']=np.where(sample_df_wityPublic['n_samples']>1, 1, 0)
    sample_df_wityPublic_grouped=sample_df_wityPublic.groupby('aminoAcid').mean()
    columns_to_keep=['count (reads)', 'frequencyCount (%)', 'cdr3Length', 'vDeletion', 'n1Insertion', 'd5Deletion', 'd3Deletion', 'n2Insertion', 'jDeletion', 'vIndex', 'n1Index',
                     'dIndex','n2Index', 'jIndex', 'estimatedNumberGenomes', 'fractionNucleated', 'n_samples', 'shareStatus']
    sample_df_wityPublic_grouped=sample_df_wityPublic_grouped[columns_to_keep]
    
    return sample_df_wityPublic, sample_df_wityPublic_grouped

#-----------------------------------------------------------------------------------------------------

def calculate_public_res_df(result_df,function_list, sample_name, sample_df_prod_wityPublic_grouped, sample_df_non_prod_wityPublic_grouped):
    print 'calculating public results....'
    ind=0
    df_list=[sample_df_prod_wityPublic_grouped, sample_df_non_prod_wityPublic_grouped] ## df list is generated here and doesnt come 
                                                 ## as input as it contains the specfic dfs
                                                 ##for each sample   
    
    for j, func in enumerate(function_list):
        for i, df in enumerate(df_list):
            result=func(df)
            result_df.loc[sample_name, result_df.columns[ind]]=result
            ind+=1
    print 'finished calculating results....'
    return result_df

#---------------------------------------------------------------------------------------------------------------
def generate_res_DF(sample_names, function_list, df_n):
    print 'generating result dataframe...'
    result_df=pd.DataFrame({'Sample': sample_names})
    result_df=result_df.set_index('Sample')
    function_names=[f.__name__ for f in function_list]
    for function in function_names:
        for i in range(df_n):
            col_title=function+'_df_%s' %str(i)
            result_df[col_title]=np.NaN
    print 'finished generating result dataframe'
    return result_df
#---------------------------------------------------------------------------------------------------------
## public population functions:


def perc_public(df):
    n_public=len(df[df['n_samples']>1])
    perc_public=(float(n_public)/len(df['n_samples']))*100
    return perc_public

def public10perc(df):
    n_public10=len(df[df['n_samples']>57])
    perc_public10=(float(n_public10)/len(df['n_samples']))*100
    return perc_public10

def public50perc(df):
    n_public50=len(df[df['n_samples']>286])
    perc_public50=(float(n_public50)/len(df['n_samples']))*100
    return perc_public50

def public95perc(df):
    n_public95=len(df[df['n_samples']>544])
    perc_public95=(float(n_public95)/len(df['n_samples']))*100
    return perc_public95

def meanSharedSamples(df):
    mean_shared=df['n_samples'].mean()
    return mean_shared

def cdr3PriToPub(df):
    mean_cdr3_pri=df[df['n_samples']==1]['cdr3Length'].mean()
    mean_cdr3_pub=df[df['n_samples']>1]['cdr3Length'].mean()
    cdr3PriToPub=float(mean_cdr3_pri)/mean_cdr3_pub
    return cdr3PriToPub

def cdr3PriToPub95(df):
    mean_cdr3_pri=df[df['n_samples']==1]['cdr3Length'].mean()
    mean_cdr3_pub95=df[df['n_samples']>544]['cdr3Length'].mean()
    cdr3PriToPub95=float(mean_cdr3_pri)/mean_cdr3_pub95
    return cdr3PriToPub95

#------------------------------------------------------------------------------------------------------

'''
for calculating public sequence statistics per sample:
>Merge the n_sample (number of samples per sequence) from the counter df into the prod and non-prod df of each 
sample. Define for each sequence whether it is private or public and save df ('sample_df_prod_wityPublic',
'sample_df_non_prod_wityPublic' (Note! In this dfs, the raws are nt sequences, while the share status is defined 
by the aa sequence. Therefore, duplication are possible, pending on the analysis done!
 
>Group the df according to aa sequences
'''
basePath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/public_analysis' ## define the path to 
n_samples=573
generate_dfs=False

@cacheOnDisk(basePath=basePath, filename='public_analysis_results_df_%(min_sample)s_%(max_sample)s', force=True)
def calculate_public_stats(min_sample, max_sample):   
    print min_sample, max_sample
    if max_sample>573:
        max_sample=573
    n=1
    df_n=2
    save_pickles=True
    df_file_names,samples_with_df=get_sample_with_dfs() ## get list of samples with df's
    current_samples=samples_with_df[min_sample:max_sample]
    public_func_list=[perc_public, public10perc, public50perc, public95perc, meanSharedSamples, cdr3PriToPub, cdr3PriToPub95]
    public_result_df=generate_res_DF(current_samples, public_func_list, df_n)
    
    print 'loading counters...'
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/prod_aa_seq_counter_df-n_samples-%s' %n_samples, "rb") as f1:
        prod_aa_seq_counter_df=pickle.load(f1)
    f1.close()
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/non_prod_aa_seq_counter_df-n_samples-%s' %n_samples, "rb") as f2:
        non_prod_aa_seq_counter_df=pickle.load(f2)
    f2.close()
    prod_aa_seq_counter_df.set_index('aa_sequence', inplace=True)
    non_prod_aa_seq_counter_df.set_index('aa_sequence', inplace=True)
    print 'fINISHED loading counters...'
    
    for d in range(len(current_samples)): ##***change here for more samples!***
        sample_name=current_samples[d]
        print n
        print sample_name
    ## extract prod and non-prod dfs: 
        sample_df, sample_df_prod, sample_df_non_prod=get_sample_data(sample_name, generate_dfs)
        sample_df_prod_wityPublic, sample_df_prod_wityPublic_grouped=add_public(sample_df_prod,prod_aa_seq_counter_df)
        sample_df_non_prod_wityPublic, sample_df_non_prod_wityPublic_grouped=add_public(sample_df_non_prod,
        non_prod_aa_seq_counter_df)
        if save_pickles:
            with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_prod_wityPublic_%s' %sample_name,'wb') as f1:
                pickle.dump(sample_df_prod_wityPublic,f1)
            f1.close()
            with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_prod_wityPublic_grouped_%s' %sample_name,'wb') as f2:
                pickle.dump(sample_df_prod_wityPublic_grouped,f2)
            f2.close()
            with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_non_prod_wityPublic_%s' %sample_name,'wb') as f3:
                pickle.dump(sample_df_non_prod_wityPublic,f3)
            f3.close()
            with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_non_prod_wityPublic_grouped_%s' %sample_name,'wb') as f4:
                pickle.dump(sample_df_non_prod_wityPublic_grouped,f4)
            f4.close()
        public_result_df=calculate_public_res_df(public_result_df,public_func_list, sample_name, sample_df_prod_wityPublic_grouped, sample_df_non_prod_wityPublic_grouped)
        n=n+1
    return public_result_df

## send jobs to the queue (code adopted from Elad, see lab wiki for more details:
sethandlers()
os.chdir(basePath)
## the inputs for qp are: jobname, q=machine list, *** add max_r to prevent exploding the cluster!!***
with qp('public_stats_job',  q = ['himem7.q','16g.q'], mem_def="10G", trds_def=2, deleteCSHwithnoerr = True, tryrerun = False, max_u=120) as q:
    q.startpermanentrun()
    wait_on =[]
    
##now define a loop that divide the job and send each part seperatly:
    min_sample=0
    max_sample=10 
    while min_sample<573:                                     
        print min_sample
        wait_on.append(q.method(calculate_public_stats,kwargs={'min_sample':min_sample,'max_sample':max_sample}))
            ##q.method takes the desired function with its arguments and send it to the queue.
        min_sample=min_sample+10
        max_sample=max_sample+10
    q.wait(wait_on)


##with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/public_res_df_%s_samples' %n_samples, "wb") as f5:
##    pickle.dump(public_result_df, f5)    
##f5.close() 

    



