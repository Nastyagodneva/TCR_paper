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
from SegalQueue.qp import qp,fakeqp
from addloglevels import sethandlers
import logging 
from Utils import cacheOnDisk
import os

#----------------------------------------








basePath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/public_analysis' ## define the path to 
n_samples=573
generate_dfs=False

@cacheOnDisk(basePath=basePath, filename='PairSharing_%(min_sample)s_%(max_sample)s', force=True)
def PairSharingCalc(min_sample, max_sample):   
    print min_sample, max_sample
    if max_sample>573:
        max_sample=573
    n=1
    df_file_names,samples_with_df=get_sample_with_dfs()
    unique_aa_prod_dict1={}
    unique_aa_non_prod_dict1={}
    for d in samples_with_df[min_sample:max_sample]: ##***change here for more samples!***
        sample_name=d
        print n
        print sample_name
        sample_df, sample_df_prod, sample_df_non_prod=get_sample_data(sample_name, False)
        unique_aa_prod_dict1[sample_name]=sample_df_prod['aminoAcid'].unique()
        unique_aa_non_prod_dict1[sample_name]=sample_df_non_prod['aminoAcid'].unique()
        n+=1
    return unique_aa_prod_dict1,unique_aa_non_prod_dict1


## send jobs to the queue (code adopted from Elad, see lab wiki for more details:
sethandlers()
os.chdir(basePath)
## the inputs for qp are: jobname, q=machine list, *** add max_r to prevent exploding the cluster!!***
with qp('PairSharing_job',  q = ['himem7.q','16g.q'], mem_def="10G", trds_def=2, deleteCSHwithnoerr = True, tryrerun = False, max_u=120) as q:
    q.startpermanentrun()
    wait_on =[]
    
##now define a loop that divide the job and send each part seperatly:
    min_sample=0
    max_sample=20 
    while min_sample<573:                                     
        print min_sample
        wait_on.append(q.method(PairSharingCalc,kwargs={'min_sample':min_sample,'max_sample':max_sample}))
            ##q.method takes the desired function with its arguments and send it to the queue.
        min_sample=min_sample+20
        max_sample=max_sample+20
    q.wait(wait_on)


