'''
this file contains 'calculating' functions for the 'PopulationAnalysis.py' module:

'''



from os import listdir
from os.path import isfile, join
# from Utils import Load, Write
import pandas as pd
import numpy as np
from scipy import stats
# import math
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib import cm
# import plotly.plotly as py
from myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot
from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
from Bio.SeqUtils import GC
from skbio.diversity.alpha import shannon, simpson, berger_parker_d
import time



#---------------------------------------------------------------------------------------------      

## this function generates a basline df for the results:
## rows are sample names
## columns are function names multiplied by df numbers, for example, if there are two dfs (prod and non prod)
## each function will repeat in two columns:
## normally 'df_0' maens 'prod' and 'df_1' means 'non_prod'


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

#-----------------------------------------------------------------------------------------------------------------
## the following functions calculate general population parameters:

def perc_prod(full_df):
    print 'calculating percent productive...'
    sample_prod = full_df.groupby(['sequenceStatus'])[['frequencyCount (%)']].sum()
    perc_prod=list(sample_prod['frequencyCount (%)'])[0]
    return perc_prod

def unique_seq_n(df): 
    n_seqs=len(list(df['nucleotide']))
    return n_seqs

def unique_aa_n(df):
    list_aa=list(df['aminoAcid'])
    list_aa_new=[i for i in list_aa if isinstance(i, str)]
    n_seq_aa=len(set(list_aa_new))
    return n_seq_aa

def max_nt_per_aa(df):
    nt_per_aa = df.groupby(['aminoAcid'])[['nucleotide']].count()
    max_nt_per_aa=max(nt_per_aa['nucleotide'])
    return max_nt_per_aa

def mean_nt_per_aa(df):
    nt_per_aa = df.groupby(['aminoAcid'])[['nucleotide']].count()
    mean_nt_per_aa=round(np.mean(nt_per_aa['nucleotide']),3)
    return mean_nt_per_aa

def norm_uniqe_nt_sequences(df): 
    repeats=10
    samp_size=100000
    reads=list(df['count (reads)'])
    df=df.set_index('nucleotide')
    seqs=[str(i) for i in list(df.index)]
    seq_popped=[]
    for i in range(0,len(seqs)):
        for j in range(0,reads[i]):
            seq_popped.append(seqs[i])        
    seq_n_list=[]
    for t in range(repeats):
        rand_seq=np.random.choice(seq_popped, samp_size, replace=False)
        seq_n=len(set(rand_seq))
        seq_n_list.append(seq_n)
        mean_seq_n=np.mean(seq_n_list)
    return mean_seq_n         

def norm_uniqe_aa_sequences(df): 
    repeats=10
    samp_size=5000
    reads=list(df['count (reads)'])
    list_aa=list(df['aminoAcid'])  
    seq_aa_popped=[]
    for i in range(0,len(list_aa)):
        for j in range(0,reads[i]):
            seq_aa_popped.append(list_aa[i])
    seq_aa_popped_new=[i for i in seq_aa_popped if isinstance(i, str)]
    seq_n_list=[]
    for t in range(repeats):
        rand_seq=np.random.choice(seq_aa_popped_new, samp_size, replace=False)
        seq_n=len(set(rand_seq))
        seq_n_list.append(seq_n)
        mean_seq_n=np.mean(seq_n_list)
    return mean_seq_n         

def gc_content(df):
    seqs=list(df['nucleotide'])
    gc_values = [GC(seq) for seq in seqs]
    mean_gc=np.mean(gc_values)
    return mean_gc

#----------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------
##clonality functions:

def gen_clonality_list_nt(df):
    clonality_list=list(df['frequencyCount (%)'])
    clonality_total=sum(clonality_list)
    normed_clonality_list=[i/clonality_total for i in clonality_list]
    sorted_clonality_list_nt=sorted(normed_clonality_list, reverse=True)
    return sorted_clonality_list_nt

def top_clonal_nt(df):
    sorted_clon_list_nt=gen_clonality_list_nt(df)
    top_clone_freq=sorted_clon_list_nt[0]
    return  top_clone_freq

def top_1000clons_nt(df):
    sorted_clon_list_nt=gen_clonality_list_nt(df)
    top1000fract=0
    for i in range(1000):
        top1000fract+=sorted_clon_list_nt[i]
    return top1000fract

def Percentile1_clone_nt(df):
    sorted_clon_list_nt=gen_clonality_list_nt(df)
    perc1=np.percentile(sorted_clon_list_nt,1)
    return perc1

def median_clone_nt(df):
    sorted_clon_list_nt=gen_clonality_list_nt(df)
    median_clone=np.percentile(sorted_clon_list_nt,50)
    return median_clone

def Percentile999_clone_nt(df):
    sorted_clon_list_nt=gen_clonality_list_nt(df)
    perc999=np.percentile(sorted_clon_list_nt,99.9)
    return perc999

def mean_clonal_nt(df):
    sorted_clon_list_nt=gen_clonality_list_nt(df)
    mean_clonal=np.mean(sorted_clon_list_nt)
    return mean_clonal

def std_clonal_nt(df):
    sorted_clon_list_nt=gen_clonality_list_nt(df)
    std_clonal=np.std(sorted_clon_list_nt)
    return std_clonal

def gen_clonality_list_aa(df):
    aa_freqs = df.groupby(['aminoAcid'])[['frequencyCount (%)']].sum()
    aa_freq_list=list(aa_freqs['frequencyCount (%)'])  
    aa_freq_total=sum(aa_freq_list)
    normed_clonality_list_aa=[i/aa_freq_total for i in aa_freq_list]
    sorted_clonality_list_aa=sorted(normed_clonality_list_aa, reverse=True)
    return sorted_clonality_list_aa

def top_clonal_aa(df):
    sorted_clon_list_aa=gen_clonality_list_aa(df)
    top_clone_freq=sorted_clon_list_aa[0]
    return  top_clone_freq

def top_1000clons_aa(df):
    sorted_clon_list_aa=gen_clonality_list_aa(df)
    if len(sorted_clon_list_aa)>1000:
        top1000fract=0
        for i in range(1000):
            top1000fract+=sorted_clon_list_aa[i]
    else:
        top1000fract=1.0
    return top1000fract

def Percentile1_clone_aa(df):
    sorted_clon_list_aa=gen_clonality_list_aa(df)
    perc1=np.percentile(sorted_clon_list_aa,1)
    return perc1

def median_clone_aa(df):
    sorted_clon_list_aa=gen_clonality_list_aa(df)
    median_clone=np.percentile(sorted_clon_list_aa,50)
    return median_clone

def Percentile999_clone_aa(df):
    sorted_clon_list_aa=gen_clonality_list_aa(df)
    perc999=np.percentile(sorted_clon_list_aa,99.9)
    return perc999

def mean_clonal_aa(df):
    sorted_clon_list_aa=gen_clonality_list_aa(df)
    mean_clonal=np.mean(sorted_clon_list_aa)
    return mean_clonal

def std_clonal_aa(df):
    sorted_clon_list_aa=gen_clonality_list_aa(df)
    std_clonal=np.std(sorted_clon_list_aa)
    return std_clonal
      

#------------------------------------------------------------------------------------------------------------
## diversity functions:

def shannon_div_nt(df):
    read_counts=list(df['count (reads)'])
    H=shannon(read_counts, base=2)
    return H

def simpson_div_nt(df):
    read_counts=list(df['count (reads)'])
    D=simpson(read_counts)
    return D

def bpi_div_nt(df):
    read_counts=list(df['count (reads)'])
    D=berger_parker_d(read_counts)
    return D

def shannon_div_aa(df):
    aa_reads = df.groupby(['aminoAcid'])[['count (reads)']].sum()
    read_counts=list(aa_reads['count (reads)']) 
    H=shannon(read_counts, base=2)
    return H

def simpson_div_aa(df):
    aa_reads = df.groupby(['aminoAcid'])[['count (reads)']].sum()
    read_counts=list(aa_reads['count (reads)']) 
    D=simpson(read_counts)
    return D

def bpi_div_aa(df):
    aa_reads = df.groupby(['aminoAcid'])[['count (reads)']].sum()
    read_counts=list(aa_reads['count (reads)']) 
    D=berger_parker_d(read_counts)
    return D
    
    #--------------------------------------------------------------------------------------------------------
##length functions: 


def mean_cdr3(df):
    mean_cdr3=np.mean(df['cdr3Length'])
    return mean_cdr3

def sv_cdr3(df):
    sv_cdr3=np.std(df['cdr3Length'])/np.mean(df['cdr3Length'])
    return sv_cdr3

def mean_vDeletion(df):
    mean_vDeletion=np.mean(df['vDeletion'])
    return mean_vDeletion

def sv_vDeletion(df):
    sv_vDeletion=np.std(df['vDeletion'])/np.mean(df['vDeletion'])
    return sv_vDeletion

def mean_n1Insertion(df):
    mean_n1Insertion=np.mean(df['n1Insertion'])
    return mean_n1Insertion

def sv_n1Insertion(df):
    sv_n1Insertion=np.std(df['n1Insertion'])/np.mean(df['n1Insertion'])
    return sv_n1Insertion

def mean_d5Deletion(df):
    mean_d5Deletion=np.mean(df['d5Deletion'])
    return mean_d5Deletion

def sv_d5Deletion(df):
    sv_d5Deletion=np.std(df['d5Deletion'])/np.mean(df['d5Deletion'])
    return sv_d5Deletion

def mean_d3Deletion(df):
    mean_d3Deletion=np.mean(df['d3Deletion'])
    return mean_d3Deletion

def sv_d3Deletion(df):
    sv_d3Deletion=np.std(df['d3Deletion'])/np.mean(df['d3Deletion'])
    return sv_d3Deletion

def mean_n2Insertion(df):
    mean_n2Insertion=np.mean(df['n2Insertion'])
    return mean_n2Insertion

def sv_n2Insertion(df):
    sv_n2Insertion=np.std(df['n2Insertion'])/np.mean(df['n2Insertion'])
    return sv_n2Insertion

def mean_jDeletion(df):
    mean_jDeletion=np.mean(df['jDeletion'])
    return mean_jDeletion

def sv_jDeletion(df):
    sv_jDeletion=np.std(df['jDeletion'])/np.mean(df['jDeletion'])
    return sv_jDeletion

#----------------------------------------------------------------------------------------------------------
    

'''
the following function calls different sub-functions to calculate  features of
the sample population and generates a result dataframe:
***note that the result df is filled column by column, regardless of the column titles, so need to make
sure that the column order is compatible with the column titles as generated by the 'generate_res_DF' 
function

input:
***'result_df' - the evolving result df. in the first call to this function, it contains
all samples as its rows, and the number of columns equal the number of functionsXnumber 
of dfs. the column content is defined as 0 in the begining. after each call to the function, 
the row of the 'calling' sample is filled with result, and the df is re-entered as input in 
the next call. 

***function list, sample name and sample dfs

     
'''

def calculate_res_df(result_df,function_list, sample_name, sample_df_prod, sample_df_non_prod):
    print 'calculating results....'
    ind=0
    df_list=[sample_df_prod, sample_df_non_prod] ## df list is generated here and doesnt come 
                                                 ## as input as it contains the specfic dfs
                                                 ##for each sample   
    
    for j, func in enumerate(function_list):
        for i, df in enumerate(df_list):
            result=func(df)
            result_df.loc[sample_name, result_df.columns[ind]]=result
            ind+=1
    print 'finished calculating results....'
    return result_df

#----------------------------------------------------------------------------------------



'''
the following function plots the histograms of the data collected in the result_df
the image contains subplots in number equal to the number of functions. all dfs with the same function
are plotted in the same subplot. 
the 'title' input should be a string indicating the name of the result df (general/clonality/
length/geneUsage/public)

'''

def plot_res_df(result_df, function_list, title, yscale, hspace):
    print 'plotting results...'
    result_df.dropna(inplace=True)
    function_names=[f.__name__ for f in function_list]
    fig=plt.figure(figsize=(6,11))
    plt.suptitle('%s population features' %title, fontsize=16)
    n_plots=len(function_list)
    print n_plots
    for p in range(n_plots):
        print p
        function_name=function_names[p]
        prod_col=2*p
        non_prod_col=2*p+1
        ax= plt.subplot2grid((n_plots,1), (p,0)) 
        plot=plot_population_view(ax, p,  result_df, function_name, yscale, prod_col, non_prod_col)
    plt.subplots_adjust(left=0.14,bottom=0.08, right=0.88, top=0.92, wspace=0.24,hspace=hspace)
    return fig

def plot_population_view(ax,p, result_df,function_name, yscale, prod_col,non_prod_col):
      
    prod_l=list(result_df[result_df.columns[prod_col]])
    if non_prod_col!=None:
        non_prod_l=list(result_df[result_df.columns[non_prod_col]])
        plot=ax.hist((prod_l,non_prod_l), bins=50, color=('blue', 'red'), label=('Productive','Non-Productive'), alpha=0.7)
        ks_p, t_p=stat_tests(prod_l,non_prod_l)
        ax.annotate('KS_p_value=%s\nt-test_p_value=%s' %(ks_p, t_p), xy=(0.95, 0.95), xycoords='axes fraction', fontsize=8,
        horizontalalignment='right', verticalalignment='top', fontweight='bold')
        
    else:
        plot=ax.hist(prod_l, color='blue', bins=50)
    ax.set_title(str(function_name), fontsize=8)
    ax.set_ylabel('Frequency', fontsize=7)
    ax.tick_params(labelsize=6)
    ax.set_yscale(yscale)
    if p==0:
        ax.legend(loc='upper center', fontsize=6)
    return plot


def stat_tests(prod_l,non_prod_l):
    ks_s,ks_p=stats.ks_2samp(prod_l,non_prod_l)
    t_s,t_p=stats.ttest_ind(prod_l,non_prod_l)
    if ks_p<=10**-4:
        ks_p='<10^-4'
    else:
        ks_p=str(round(ks_p,4))
    if t_p<=10**-4:
        t_p='<10^-4'
    else:
        t_p=str(round(t_p,4))
    return ks_p, t_p



#---------------------------------------------------------------------------------------------------------------------------
## the following function is an alternative for the plot_res_df function, which 
##which should be used for plotting result dfs in which the function was executed
##over the sample_df and not sample_df_prod and sample_df_non_prod

def plot_1_res_df(result_df, function_list, title,yscale):
    print 'plotting results...'
    result_df.dropna(inplace=True)
    function_names=[f.__name__ for f in function_list]
    fig=plt.figure(figsize=(10,8))
    plt.suptitle('%s population features' %title, fontsize=16)
    n_plots=len(function_list)
    print n_plots
    for p in range(n_plots):
        print p
        function_name=function_names[p]
        df_col=p
        ax= plt.subplot2grid((n_plots,1), (p,0)) 
        plot=plot_population_view(ax, p,  result_df, function_name, yscale, df_col, non_prod_col=None)
    plt.subplots_adjust(left=0.09,bottom=0.11, right=0.95, top=0.89, wspace=0.24,hspace=0.50)
    return fig

#-----------------------------------------------------------------------------------------

