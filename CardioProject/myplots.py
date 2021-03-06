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
'''
the following function roundup a number to the next round fold (for example if the number has 4 digits, it will be 
rounded to the next 1000s, and so on:

'''
def adjusted_roundup(num,use1LevelBelow=False):
    if num>=1:
        s_num=str(num)
        fraction=s_num.split('.')[0]
        fold=len(fraction)
        if use1LevelBelow:
            fold=fold-1
        roundup_fold=10**(fold-1)
    if num<1:
        s_num=str(num)
        fraction=s_num.split('.')[1]
        fold=len(fraction)
        roundup_fold=10**-(fold-1)
        
    return int(math.ceil(float(num) / roundup_fold)) * roundup_fold

def adjusted_rounddown(num,use1LevelBelow=False):
    if num>1:
        s_num=str(num)
        fraction=s_num.split('.')[0]
        fold=len(fraction)
        if use1LevelBelow:
            fold=fold-1
        roundup_fold=10**(fold-1)
    if num<1:
        s_num=str(num)
        fraction=s_num.split('.')[1]
        fold=len(fraction)
        roundup_fold=10**-(fold-1)
        
    return int(math.floor(float(num) / roundup_fold)) * roundup_fold

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
*sample_name (need to modify when the graph shows farefaction curves for different samples)
*df_list-meta list, in which each item is a df of a sample or sub-sample
*column_name-the column in the df to extract the relevant quantitative data. e.g. - 'count (reads)' / 
'count (templates/reads)' etc.
*sampling_interval-the interval size for sampling the number of unique sequences

'''

def rarefaction_calc(sample_name, df_list, column_name,sampling_interval):
    
    print ('rarefaction analysis for sample %s' %sample_name)

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
     
     
    max_samp_size=max(popped_list_lengths)
    print('the maximal sample size will be: %s' %max_samp_size)
    samp_size=range(0,max_samp_size,sampling_interval)
        
    n=0
    n_samp_seq_all=[]
    last_values=[]
    for popped in all_seq_popped:
        n_samp_seq=[]
        for size in samp_size:
            if size<=len(popped): ## if the sample size is larger than the amount of
                                  ## sequences in the sample, do not employ random sampling
                rand_seq=np.random.choice(popped, size, replace=False)
                seq_n=len(set(rand_seq))
                n_samp_seq.append(seq_n)
            else:
                n_samp_seq.append(np.nan)
 ##since it is time consuming to generate theses lists, save them as pickles and next time load the pickles
## (comment/uncomment each of the next couple of lines:):          
        pickle.dump(n_samp_seq, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/n_samp_seq_%s_%s.p' %(n, sample_name),"wb" ))
        n_samp_seq_all.append(pickle.load( open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/n_samp_seq_%s_%s.p' %(n, sample_name),"rb" )))
        last_values.append(n_samp_seq[-1])
        
        n=n+1
    
    ymax_rare=max(last_values)
    print('last values list is %s' %last_values)
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
    fold=np.log10(max_samp_size)
    fold_floor=int(math.floor(fold))
    xticks=range(0,roundup(max_samp_size,10**fold_floor)+10**fold_floor,10**fold_floor)
    yticks=range(0,roundup(ymax_rare,10000)+10000,10000)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(xticks, rotation=45, fontsize=8)
    ax.set_yticklabels(yticks, fontsize=8)
    ax.legend(loc='best')
    ax.set_xlim(0, roundup(max_samp_size,10000))
    #ax.set_ylim(0, roundup(ymax_rare,10000))
    print 'Finished plotting rarefaction graphs'
    return curves



#----------------------------------------------------------------------------------------------------------------------

'''
This function calculate the correlations among all numeric variables in a df, and
the significance of the correlation.
the significance is computed in 3 different ways:
1. the inherent p-value of the scipy.stats.pearsonr function.
2. corrected p-value which equals the p-value divided by the total number of tests. 
3. permutation analysis, in which the data is permutated n times, the correlations
are calculated each time and the 2.5 and 97.5 percentile of these r's is calculated,
then the real r is compared to this confidence interval to see if it is outside
of it (=significant)

input:
df - the reault df of which correlations should be calculated,
n-per: the number of permutations to execute (100 is recommended)

'''


def calc_sig_corr_all2all(df, n_per):
    from scipy.stats import pearsonr

    ## stage 1: calculate real correlations and generate dataframe
    print 'calculating real correlations...'
    corr_list=[]
    for n1, column1 in enumerate(list(df.columns.values)):  
        for n2, column2 in enumerate(list(df.columns.values)[n1+2:]):
            if (df[column1].dtype == np.float64 or df[column1].dtype == np.int64)&(df[column2].dtype == np.float64 or df[column2].dtype == np.int64) :
                nc1=np.isnan(df[column1])
                nc2=np.isnan(df[column2])
                n=nc1|nc2
                newx=list(df[column1][~n])
                newy=list(df[column2][~n])
                r,p = pearsonr(newx,newy)
            else:
                r=np.nan
                p=np.nan
            corr_list.append({'column1':column1,'column2':column2, 'real_r':r, 'real_p':p})
    print 'finished calculating real correlations'
    res_corr_df=pd.DataFrame(corr_list)
    res_corr_df['abs_r'] = res_corr_df['real_r'].abs()
    res_corr_df.sort_index(by=['real_p','abs_r','column1','column2'], ascending=[True,False,True,True],inplace=1)
    res_corr_df.drop('abs_r', axis=1,inplace=1)


## stage 2: permutate result df and calculate correlations over permutated data:
    print 'calculating suffle correlations...'
    for i in range(n_per):
        print i
        print 'start shuffling df'
        shuffle=df.apply(np.random.permutation)
        shuf_corr_list=[]
        #column_list=[]
        for n1, column1 in enumerate(list(shuffle.columns.values)):
            for n2, column2 in enumerate(list(shuffle.columns.values)[n1+2:]):
                if (shuffle[column1].dtype == np.float64 or shuffle[column1].dtype == np.int64)&(shuffle[column2].dtype == np.float64 or shuffle[column2].dtype == np.int64):
                    nc1=np.isnan(shuffle[column1])
                    nc2=np.isnan(shuffle[column2])
                    n=nc1|nc2
                    newx=list(shuffle[column1][~n])
                    newy=list(shuffle[column2][~n])
                    r,p = pearsonr(newx,newy)
                    shuf_corr_list.append(r)
                else:
                    shuf_corr_list.append(np.nan)
            
        res_corr_df.loc[:,('r_shuf_%s'%i)] = shuf_corr_list
    
    print 'finished calculating suffle correlations'
    
    ## stage 3: calculate confidence interval for shuffled r's:
    print 'calculating percentile values...'
    #print shuffle_r_df.columns.values
    col_for_percentile=[col for col in res_corr_df.columns.values if col.startswith('r_shuf_')]
    #shuffle_r_df.loc[:,'r_mean']=np.mean(shuffle_r_df[col_for_percentile])
    #res_corr_df['avg'] = res_corr_df[col_for_percentile].mean(axis=1)
    res_corr_df['r_perc_2_5'] = res_corr_df[col_for_percentile].quantile(q=0.025,axis=1)
    res_corr_df['r_perc_97_5'] = res_corr_df[col_for_percentile].quantile(q=0.975,axis=1)
    print 'finished calculating percentile values'
    
    ## stage 4: check significance of the correlations:
    n_tests=sigma(len(df.columns.values))
    res_correct_p=0.05/n_tests
    res_corr_df['real_p_sig']=np.where(res_corr_df['real_p']<0.05,1,0)
    res_corr_df['real_p_sig_corrected']=np.where(res_corr_df['real_p']<res_correct_p,1,0)
    res_corr_df['r_outof_CI']=np.where((res_corr_df['real_r']<res_corr_df['r_perc_2_5']) | (res_corr_df['real_r']>res_corr_df['r_perc_97_5']) ,1,0)
    col_to_keep=['column1','column2','real_p','real_r','r_perc_2_5','r_perc_97_5','real_p_sig','real_p_sig_corrected','r_outof_CI']
    res_corr_df=res_corr_df[col_to_keep]
    res_corr_df.set_index('column1',inplace=1)
    res_corr_df=res_corr_df[~np.isnan(res_corr_df['real_r'])]
    only_sig_corr_df=res_corr_df[(res_corr_df['real_p_sig']==1) & (res_corr_df['r_outof_CI']==1)]
    
    return res_corr_df,only_sig_corr_df


#-------------------------------------------------------------------------------

def sigma(x):
    sigma=0
    for n in range(x):
        sigma=+n
    return sigma

#---------------------------------------------------------------------------------

def calc_sig_corr(df, corr_col, n_per):
    from scipy.stats import pearsonr

    ## stage 1: calculate real correlations and generate dataframe
    print 'calculating real correlations'
    corr_list=[]
    n_tests=len(df.columns.values)
    for column in list(df.columns.values):
        if df[column].dtype == np.float64 or df[column].dtype == np.int64 :
            n1=np.isnan(df[column])
            n2=np.isnan(df[corr_col])
            n=n1|n2
            newx=list(df[column][~n])
            newy=list(df[corr_col][~n])
            r,p = pearsonr(newx,newy)
        else:
            r=np.nan
            p=np.nan
        corr_list.append({'column':column,'real_r':r, 'real_p':p})
    res_corr_df=pd.DataFrame(corr_list)
    res_corr_df.sort(columns='real_p', inplace=1)


    ## stage 2: permutate result df and calculate correlations over permutated data:
    for i in range(n_per):
        print i
        print 'start shuffling df'
        shuffle=df.apply(np.random.permutation)
        shuf_corr_list=[]
        #column_list=[]
        for column in list(shuffle.columns.values):
            #print ('start calculating r for column %s' %column)
            #column_list.append(column)
            if shuffle[column].dtype == np.float64 or shuffle[column].dtype == np.int64 :
                n1=np.isnan(shuffle[column])
                n2=np.isnan(shuffle[corr_col])
                n=n1|n2
                newx=list(shuffle[column][~n])
                newy=list(shuffle[corr_col][~n])
                r,p = pearsonr(newx,newy)
                shuf_corr_list.append(r)
            else:
                shuf_corr_list.append(np.nan)
        #shuffle_r_df.loc[:,('column_shuf_%s'%i)] = column_list
        res_corr_df.loc[:,('r_shuf_%s'%i)] = shuf_corr_list
    
    ## stage 3: calculate confidence interval for shuffled r's:
    #print shuffle_r_df.columns.values
    col_for_percentile=[col for col in res_corr_df.columns.values if col.startswith('r_shuf_')]
    #shuffle_r_df.loc[:,'r_mean']=np.mean(shuffle_r_df[col_for_percentile])
    #res_corr_df['avg'] = res_corr_df[col_for_percentile].mean(axis=1)
    res_corr_df['r_perc_2_5'] = res_corr_df[col_for_percentile].quantile(q=0.025,axis=1)
    res_corr_df['r_perc_97_5'] = res_corr_df[col_for_percentile].quantile(q=0.975,axis=1)
    
    ## stage 4: check significance of the correlations:
    res_correct_p=0.05/n_tests
    res_corr_df['real_p_sig']=np.where(res_corr_df['real_p']<0.05,1,0)
    res_corr_df['real_p_sig_corrected']=np.where(res_corr_df['real_p']<res_correct_p,1,0)
    res_corr_df['r_outof_CI']=np.where((res_corr_df['real_r']<res_corr_df['r_perc_2_5']) | (res_corr_df['real_r']>res_corr_df['r_perc_97_5']) ,1,0)
    col_to_keep=['column','real_p','real_r','r_perc_2_5','r_perc_97_5','real_p_sig','real_p_sig_corrected','r_outof_CI']
    res_corr_df=res_corr_df[col_to_keep]
    res_corr_df.set_index('column',inplace=1)
    only_sig_corr_df=res_corr_df[(res_corr_df['real_p_sig']==1) & (res_corr_df['r_outof_CI']==1)]
    
    return res_corr_df,only_sig_corr_df


#-------------------------------------------------------------------------------------

def draw_correlation_scatter(x, y, figsize = (3, 3), xticks=None, yticks=None,\
                             xlim = None, ylim = None, r = None, ms=4, logd = False,\
                             xlab = None, ylab = None, filename = None, title = None,
                             color = "#a0a0a0", grid = True, dpi = 800, xticklabels = None, 
                             contour = False, **figkwargs):
    from scipy.stats import pearsonr,spearmanr
    fig = plt.figure(figsize = figsize, dpi = dpi)
    axB = fig.add_subplot(111)
    if contour:
        print "Contour plot are experimental here"
        #axB.hist2d(x,y,bins = 40,norm=LogNorm())
    else:
        axB.plot(x, y, 'o', color = color, ms=ms, **figkwargs)
    if logd:
        axB.set_xscale('log',basex=2)
        axB.set_yscale('log',basey=2)
    if xticks is not None:
        axB.set_xticks(xticks)
        axB.set_xticklabels(xticks)
    if xticklabels is not None:
        axB.set_xticklabels(xticklabels)
    if yticks is not None:
        axB.set_yticks(yticks)
        axB.set_yticklabels(yticks)
    if xlim is not None:
        axB.set_xlim(xlim)
    if ylim is not None:
        axB.set_ylim(ylim)
    if r is not None: 
        if r == 'pearson':
            n=np.isnan(x)
            newx=list(x[~n])
            newy=list(y[~n])
            r,p = pearsonr(newx,newy)
            axB.text(0.01,0.99,"r=%.4f p=%.4f" %(r,p), transform=axB.transAxes, verticalalignment = 'top', ha = 'left',fontsize=14,color='red')
        if r == 'spearman':
            n=np.isnan(x)
            newx=list(x[~n])
            newy=list(y[~n])
            r,p = spearmanr(newx,newy)
            axB.text(0.01,0.99,"r=%.4f p=%.4f" %(r,p), transform=axB.transAxes, verticalalignment = 'top', ha = 'left',fontsize=14,color='red')
    if xlab is not None:
        axB.set_xlabel(xlab)
    if ylab is not None:
        axB.set_ylabel(ylab)
    if title is not None:
        plt.title(title, fontsize=16)
    if grid:
        axB.grid()
    if filename is not None:
        fig.savefig(filename, bbox_inches='tight', dpi = dpi)
    axB.margins(0.1, 0.1)
    #axB.set_xmargin(0.2); axB.autoscale_view()
    return fig, axB
