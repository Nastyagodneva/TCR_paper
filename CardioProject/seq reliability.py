'''
this function takes all the sequences in a dataframe and com[are them to each other. it looks for matches
between sequences with x allowed mismatches the function divides the sequences into bulks of 500
sequences each and send them to the cluster. 
*** next time send fewer sequences in the first bulks and more sequences in the last bulks as the first 
are much longer, maybe also divide to smaller bulks. 
need to marge the resulting dfs


'''



import os
#from os import listdir
#from os.path import isfile, join
# from Utils import Load, Write
import pandas as pd
import numpy as np
#from scipy import stats
# import math
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib import cm
# import plotly.plotly as py
from myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot
#from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
#from Bio.SeqUtils import GC
from pop_calcs import *
from pop_organize import *
from queue.qp import qp,fakeqp
from addloglevels import sethandlers
import logging 
from Utils import cacheOnDisk
from bm_preproc import BoyerMoore



#----------------------------------------------------------------------

## the following function compare two sequences: p and t, and returns the number of mismatches between them, up to a maximal number of
## mismatches - x.
def naive_mm_seq_rel(p,t,x):
##  x - number of mismatches allowed
    print 'start sequence matching'
    if len(t)<len(p): ## if len(t)< len(P), the sequences of p and t are exchanged
        p2=p
        t2=t
        p=t2
        t=p2
        print 'replaced t and p'
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        mm_count=0
        print mm_count
        match = True
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                mm_count+=1
            if mm_count>=x: ## if the number of mismatches exceeds x, break.
                match = False
                break
    return mm_count



#----------------------------------------------------------
def boyer_moore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences

#-----------------------------------------------------------------------------------------
def approximate_match (p,t,n):
    #n-max number of mismatches allowed, 
    
    #make the sequences reciprocal:
    if len(t)<len(p): ## if len(t)< len(P), the sequences of p and t are exchanged
        p2=p
        t2=t
        p=t2
        t=p2
        print 'replaced t and p'
    n_mis=[]
    segment_length=round(len(p)/(n+1))
    if segment_length>=2:  
        for i in range(n+1): # loop over partitions
            start=int(i*segment_length) #define start and end site for each partition
            end=int(min((i+1)*segment_length,len(p)))
            p_bm=BoyerMoore(p[start:end], alphabet='ACGT')
            matches=boyer_moore(p[start:end], p_bm,t) 

            for m in matches: #verify match outside the partition, taking into account
                              #the allowed number of mismatches n
                if m<start or m-start+len(p) >len(t):
                    continue

                mismatches=0
                for j in range(0,start):
                    if not p[j] == t[m-start+j]:
                        mismatches+=1
                        if mismatches>n:
                            break

                for j in range(end, len(p)):
                    if not p[j] == t[m-start+j]:
                        mismatches+=1
                        if mismatches>n:
                            break

                n_mis.append(mismatches)
        if len(n_mis)>0:
            min_mis=min(n_mis)
        else:
            min_mis=1000
    else:
        min_mis=1000
   
    return min_mis

#--------------------------------------------------------------------------------





'''
from now on, the script is intended to divide the sample_df into parts, and run the sequence comparison process over each part
seperately in the queue. I adopted a code from Elad.
basepath=the path in which the helper files will be saved, and the resulting match dataframes will be saved. 
*** for each run: define sample_name (if changed, need to change to force=True in the df_len cacheOnDisk for 1 run.
                  

'''

basePath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/seqRel' ## define the path to 

## count the number of rows in the dataframe tested:
@cacheOnDisk(basePath=basePath,force=True)
def df_len(sample_name):
    print ('counting number of rows in sample %s' %sample_name)
    sample_df, sample_df_prod, sample_df_non_prod=get_sample_data(sample_name, True)
    df_len=len(sample_df.index)
    print('finished counting number of rows in sample %s' %sample_name)
    return df_len

## cacheOnDisk saves the result of the underlying function into the file 'filename' in the path 'basepPath' (force=True if
## overwrite is desired.
## the following function iterates over all the sequences in a df and compare each to all the subsequent sequences.
## it defines sequences p and t for the naive_mm_seq_rel function and calls it. 
## then it gets the number of mm between the sequences, and document all matches with less then the maximal number of mm in a 
## dataframe called seq_match_df.

@cacheOnDisk(basePath=basePath, filename='BMmatch_df_%(min_seq)s_%(max_seq)s', force=True)
def generate_match_df_new(min_seq,max_seq,x):
    print 'start generate_match_df function...'
    sample_df, sample_df_prod, sample_df_non_prod=get_sample_data('HIP00640', False)
    print 'finished generating sample dataframes...'
    columns_to_keep = ['nucleotide', 'count (reads)']
    seq_rel_df = sample_df[columns_to_keep]
    seq_match_df=pd.DataFrame({'#seq1':[],'seq1':[], 'seq1 reads':[], '#seq2':[],'seq2':[], 'seq2 reads':[], '#mm':[]})
    print 'finished generating preliminary match df'
    n=0
    if max_seq>len(seq_rel_df):
        max_seq=len(seq_rel_df)
    for seq1 in range(min_seq, max_seq):
        print seq1
        for seq2 in range(seq1+1,len(seq_rel_df.index)):
            p=seq_rel_df.loc[seq1,'nucleotide']
            t=seq_rel_df.loc[seq2,'nucleotide']
            min_mis=approximate_match (p,t,x)
            if min_mis<=x:
                seq_match_df.loc[n,'#mm']=min_mis
                seq_match_df.loc[n,'#seq1']=seq1
                seq_match_df.loc[n,'seq1']=p
                seq_match_df.loc[n,'seq1 reads']=seq_rel_df.loc[seq1,'count (reads)']
                seq_match_df.loc[n,'#seq2']=seq2
                seq_match_df.loc[n,'seq2']=t
                seq_match_df.loc[n,'seq2 reads']=seq_rel_df.loc[seq2,'count (reads)']
            n=n+1
    return seq_match_df

## send jobs to the queue (code adopted from Elad, see lab wiki for more details:
sethandlers()
os.chdir(basePath)
## the inputs for qp are: jobname, q=machine list, *** add max_r to prevent exploding the cluster!!***
with qp('BMmatch_df_job',  q = ['himem7.q','16g.q'], mem_def="10G", trds_def=2, deleteCSHwithnoerr = True, tryrerun = True, max_u=120) as q:
    q.startpermanentrun()
    wait_on =[]
    
##now define a loop that divide the job and send each part seperately:
## consider making a list of integer for the min/max_seq so the first jobs will be shorter than the last ones. 
    min_seq=0
    max_seq=100 ##change back to 500
    while min_seq<1000:                                      ## change back to: df_len(sample_name):
        print min_seq
        wait_on.append(q.method(generate_match_df_new,kwargs={'min_seq':min_seq,'max_seq':max_seq,'x':3}))
            ##q.method takes the desired function which its arguments and send it to the queue.
        min_seq=min_seq+100
        max_seq=max_seq+100
    q.wait(wait_on)

    
#---------------------------------------------------------------------------------
'''
def gen_combined_match_df():
##get the df file names as a list (note!! to take only the files relevant to the specific samples!!***
    filenames = [f for f in listdir('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/seqRel/generate_match_df') if isfile(join('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/seqRel/generate_match_df', f))]
    #filenames = [datafile for datafile in filenames if datafile.startswith ('HIP') and datafile.endswith('.csv')]
    #df_names=[re.sub('.csv', '', datafile) for datafile in filenames]
    print len(filenames)

##generate a list of all dfs extracted:    
    match_df_list=[]
    for df in filenames:
        with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/seqRel/generate_match_df/%s' %df, 'rb') as f:
            match_df=pickle.load(f)
        f.close()
        match_df_list.append(match_df)
    print len(match_df_list)
    
    match_df_combined=pd.concat(match_df_list)
    print len(match_df_combined)
    return match_df_combined

#---------------------------------------------------------------------------------------------------------------------

match_df_combined=gen_combined_match_df()

## save the combined dataframe as pickle
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/seqRel/generate_match_df/match_df_combined', 'wb') as f:
        pickle.dump(match_df_combined,f)
f.close()
        
## generate a dataframe which is sub-df of the combined matched df which includes only matches in which at least one of the 
## sequences has only 2 reads.
err_df=match_df_combined[(match_df_combined['seq1 reads']<3) | (match_df_combined['seq2 reads']<3)]
print len(err_df)

## look for sequences in which the more prevalent sequence in the pair of matched sequences has more than 100 reads:
err_df_final=err_df[(err_df['seq1 reads']>100) | (err_df['seq2 reads']>100)]
print len(err_df_final)

## ***decide the final parameters to cut, and decide what to do with the eeoneous sequences (
'''
