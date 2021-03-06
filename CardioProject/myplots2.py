'''
this code contains general functions to generate plots. or relevant data, for sample view generations
(for files generated by the immunoseq kit
'''




import pandas as pd
import numpy as np
# from scipy.stats import mannwhitneyu
import math
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# import collections
# from itertools import chain
from Utils import Load, Write
from venn import venn
import cPickle as pickle



#-----------------------------------------------------------------------------------------------
def roundup(x, fold):
## this function rounds x up to the closest ten/hundred/thoudand etc, according to 'fold' definition
    return int(math.ceil(x / float(fold))) * fold  

#-------------------------------------------------------------------------------------------

def rounddown(x, fold):
## this function rounds x down to the closest ten/hundred/thoudand etc, according to 'fold' definition
    return int(math.floor(x / float(fold))) * fold  

#-------------------------------------------------------------------------------------------------
'''
the following function calculates the values corresponding to the p percentile in each of n lists, and
then a roundup to the maximal value of all these values.
the result can be used as a cutoff in graphs
data-list of lists. make sure that if there is only one data list, it is still incapsulated with []
cut_percentile - the desired percentile
round -the fold to which to roundup
'''

def percentile_cut_off(data,cut_percentile, round_fold):
    if data is None:
        raise Exception("No data!")

    perc_cuts=[]
    for d in data:
        perc_cut=np.percentile(d,cut_percentile)
        perc_cuts.append(perc_cut)
    perc_cut_max=max(perc_cuts)
    perc_cut_round=round(perc_cut_max, round_fold)+10**-round_fold
    
    return perc_cut_round

#----------------------------------------------------------------------------------------

'''
the following function gets a float number as input, and return the position of the first digit after the
decimal point which is not 0 - useful for round clculations
'''

def find_decimal_fold(num):
    s_num=str(num)
    fraction=s_num.split('.')[1]
    print fraction
    fold=0
    for i in range(len(fraction)):
        if fraction[i]!='0':
            fold=i+1
            break
    return fold

#-------------------------------------------------------------------------------------------------


'''
goal: to generate two lists of sequence frequencies, corresponsing to two different samples/replicates. 
input: two dfs generated from standard Analyzer TSV files.
SeqType can be either 'nt' or 'aa'
output should be used to generate scatter plot.
note! when the seqeunce exist only in one sample, the value in the other sample is 1e-6
'''   
def seq_freq_calc(df1, df2, SeqType, sample_name):
    global freqs11, freqs2, max_freq, min_freq
    print ('replicate analysis for sample %s...' %sample_name)  
    if SeqType=='nt':
        df1=df1.set_index('nucleotide')
        df2=df2.set_index('nucleotide')
    elif SeqType=='aa':
        df1=df1.set_index('aminoAcid')
        df2=df2.set_index('aminoAcid')
    else:
        raise NameError('please indicate SeqType')
        
        
    seq1=[str(i) for i in list(df1.index)] #generate list of sequences from df1
    seq2=[str(i) for i in list(df2.index)] #generate list of sequences from df2
    total_seq_unique=list(set(seq1+seq2))  #generate combined set of sequences from the two lists
## comment from here to avoid timewaste, when value are already stored in pickles:     
    freqs1=[] #list of frequencies for each sequence in sample 1
    freqs2=[] #list of frequencies for each sequence in sample 1
           
    for seq in total_seq_unique:           
        if seq in seq1 and seq in seq2:
            freqs1.append(df1.loc[seq,'frequencyCount (%)'])
            freqs2.append(df2.loc[seq,'frequencyCount (%)'])
        elif seq in seq1 and seq not in seq2:
            freqs1.append(df1.loc[seq,'frequencyCount (%)'])
            freqs2.append(1e-6)
        elif seq not in seq1 and seq in seq2:
            freqs1.append(1e-6)
            freqs2.append(df2.loc[seq,'frequencyCount (%)'])
      
## since it is time consuming to generate this list, save them as pickle and next time load the pickles:
    pickle.dump(freqs1, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freqs1_%(name)_%(type).p' %{'name': sample_name, 'type': SeqType},"wb" ))
    pickle.dump(freqs2, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freqs2_%(name)_%(type).p' %{'name': sample_name, 'type': SeqType},"wb"))
    freqs1= pickle.load( open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freqs1_%(name)_%(type).p' %{'name': sample_name, 'type': SeqType},"rb" ))
    freqs2= pickle.load( open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freqs2_%(name)_%(type).p' %{'name': sample_name, 'type': SeqType},"rb"))
    print ('finished replicate analysis for sample %s...' %sample_name)
    
    max_freq=max(freqs1, freqs2)
    min_freq=min(freqs1, freqs2)
    print max_freq, min_freq
#     plt.figure(figsize=(8, 8))
#     plt.suptitle('Sequences Frequency Correlation Between a1 and b1 Replicates')

    
#----------------------------------------------------------------------------------


'''
this function generates a scatter plot of sequence frequencies from two samples/frequencies
data1 and data2 are corresponding lists of frequencies of each sequence in the two samples. 
name1 and name2 are the names of the datasets compared
note that the order is critical here, as these lists carry data for the *same* sequences.
'''
      
def replicate_scatter(data1, data2, name1, name2, ax):    
    if ax is None:
        ax = plt.gca()
    scatter1=ax.scatter(data1, data2, alpha=0.1)
    dist = int(math.log10(abs(min_freq))) #calculates the number od zeros after the decimal point in the min_freq value
    ticks=np.logspace(dist-2,0,-(dist-3))#generates the x,y ticks based on the min_freq. the lower tick is used only for margin purposes and the second tick is used as '0'
    labels=[]  # generates the x,y ticklabels based on the ticks
    for i in range(0,len(ticks)):
        if i==0:
            labels.append('')
        elif i==1:
            labels.append('0')
        else:
            labels.append(str(ticks[i]))
    
    ax.set_xscale('log')
    ax.set_yscale('log')       
    ax.set_xticks(ticks,minor=False)
    ax.set_yticks(ticks,minor=False)
    ax.set_xticklabels(labels,fontsize=8)
    ax.set_yticklabels(labels,fontsize=8)
    ax.set_xlim(5*10**(dist-2),1)
    ax.set_ylim(5*10**(dist-2),1)
    ax.set_title('Sequences Frequency Correlation',fontsize=12)
    ax.set_xlabel('%s (Sequence Frequency)' %name1,fontsize=9)
    ax.set_ylabel('%s (Sequence Frequency)' %name2,fontsize=9)
    ax.margins(0.2)
    return scatter1

#-----------------------------------------------------------------------------------------------
'''
the following functions were copied from the OptSeqAnalysis script - need to generelaized!!
'''

#--------------------------------------------------------------------------------------------------
'''
the following function calculated the data need for rarefactions plots, based on df data for specific samples
input:
sample_name (need to modify when the graph shows farefaction curves for different samples)
df_list-meta list, in which each item is a df of a sample or sub-sample
column_name-the column in the df to extract the relevant quantitative data. e.g. - 'count (reads)' / 
'count (templates/reads)' etc.

'''

def rarefaction_calc(sample_name, df_list, column_name):
    
    print ('rarefaction analysis for  %s' %sample_name)

    n=0
    all_data=[]
    all_seqs=[]
    all_seq_popped=[]
    popped_list_lengths=[]
    '''
    the following loop generate three lists:
    (1) all_seqs list=contains lists of all uniqe sequences in each df contains within the df_list.
    (2) all_data list=contains lists of reads/templates number for each sequence in the all_seq_list
    (3) all_seq_popped list=contains list of popped squences lists, in which each sequence appears the number
    of times indicated by the number of reads/ template
    ***** note that each of the  lists is mata-list containing n lists of the data*******
    '''
    for n in range(len(df_list)):
        data=list(df_list[n][column_name])
        df_list[n]=df_list[n].set_index('nucleotide')
        seq=[str(i) for i in list(df_list[n].index)]
        all_data.append(data)
        all_seqs.append(seq)
        seq_popped=[]
        for i in range(0,len(seq)):
            for j in range(0,data[i]):
                seq_popped.append(seq[i])        
        all_seq_popped.append(seq_popped)
        popped_list_lengths.append(len(seq_popped))       
        n=n+1
     
     
    max_samp_size=min(popped_list_lengths)
    print('the maximal sample size will be: %s' %max_samp_size)
    samp_size=range(0,max_samp_size,10000)
        
    n=0
    n_samp_seq_all=[]
    last_values=[]
    for popped in all_seq_popped:
        n_samp_seq=[]
        for size in samp_size:
            rand_seq=np.random.choice(popped, size, replace=False)
            seq_n=len(set(rand_seq))
            n_samp_seq.append(seq_n)
 ##since it is time consuming to generate theses lists, save them as pickles and next time load the pickles
## (comment/uncomment each of the next couple of lines:):          
        pickle.dump(n_samp_seq, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/n_samp_seq_%s_%s.p' %(n, sample_name),"wb" ))
        n_samp_seq_all.append(pickle.load( open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/n_samp_seq_%s_%s.p' %(n, sample_name),"rb" )))
        last_values.append(n_samp_seq[-1])
        n=n+1
    
    ymax_rare=max(last_values)
    print('ymax for rarefaction plot is %s' %ymax_rare)
    print ('finished rarefaction analysis for sample %s...' %sample_name)
   
    return samp_size,max_samp_size, n_samp_seq_all, ymax_rare      

#------------------------------------------------------------------------------------------------------------




   
def rarefaction_plot(ax, sample_name, data, names, colors, samp_size, max_samp_size,ymax_rare):
    '''
    this function generates rarefaction graphs for one or more datasets at the same time.
    the datasets can be generated by the function rarefaction_calc. 
    *** the data is a meta-list that can contains several lists of rarefaction values, but they must be generated
    using the same samp_size list***
    input:
    sample_name
    data - see above (the udea is that in each sample size, the list of all sequences in the sample 
    is randomly samples, and the number of unique clones sampled is indicated).
    names corresponding names for the lists in 'data' rarefaction graphs.
    ***need to insert at least name1*** 
    colors: corresponding names for the lists in 'data' rarefaction graphs.
    ***need to insert at least color1***
    samp_size: list of increasing sample sizes, corresponding to the values in the data sets
    max_samp_size: the length of the shorter list among the dataset
    y_max - the maximal value in all data sets together
    ax - the number of subplot in the original image   
    '''
    
    
    print 'plotting rarefaction graphs'
    if ax is None:
        ax = plt.gca()
    curves=[]
    for i in range(len(data)):
        curve=ax.plot(samp_size,data[i],label=names[i],color=colors[i])
        curves.append(curve)
  
    ax.set_title('Rarefaction Plots', fontsize=12)
    ax.set_xlabel('Sample Size',fontsize=9)
    ax.set_ylabel('Number of unique sequences',fontsize=9)
    xticks=range(0,roundup(max_samp_size,100000)+100000,100000)
    yticks=range(0,roundup(ymax_rare,10000)+10000,10000)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(xticks, rotation=45, fontsize=8)
    ax.set_yticklabels(yticks, fontsize=8)
    ax.legend(loc='best')
    ax.set_xlim(0, roundup(max_samp_size,10000))
    ax.set_ylim(0, roundup(ymax_rare,10000))
    print 'Finished plotting rarefaction graphs'
    return curves