'''
new sample view generator- general and documented
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


#---------------------------------------------------------------------------------------
 ## generates dataframes from the relevant tsv file(s), and then generates "sub-dfs" for only productive and
 ## only non-productive seqeunces. 
 
  
def get_sample_data(sample_name): 
    
    print 'getting sample data...'

    sample_df=pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/%s.tsv" %sample_name)  
    sample_df_prod = sample_df[sample_df['sequenceStatus'] == 'In']
    sample_df_non_prod = sample_df[sample_df['sequenceStatus'] != 'In']
 
 ## save dfs to pickles:   
    pickle.dump(sample_df, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_%s' %sample_name, "wb"))
    pickle.dump(sample_df_prod, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_prod_%s' %sample_name, "wb"))
    pickle.dump(sample_df_non_prod, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_non_prod_%s' %sample_name, "wb"))


##uncomment to generate df and/or csv files:  
# sample_df_file=Write('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/DFfiles/%s_data.df' %sample_name, sample_df)
# sample_df.to_csv('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/CSVfiles/%s_dat.csv' %sample_name)
    print 'finished getting sample data...'
    return sample_df, sample_df_prod, sample_df_non_prod

#-----------------------------------------------------------------------------------------------------------


    
def prod_dist(ax, sample_df, sample_df_prod, sample_df_non_prod, sample_name):
##  generate pie charts showing the proportions of productive and non productive reads
## returns the total number of reads, productive reads and non productive reads, and the pie figure data

    print 'working on the productive pie chart...'
    sample_prod = sample_df.groupby(['sequenceStatus'])[['frequencyCount (%)']].sum()
    sample_prod_list=list(sample_prod['frequencyCount (%)'])
    sample_prod_list=[round(i,2) for i in sample_prod_list]
    labels=[str(i)+'%' for i in sample_prod_list]  
    
    legend=sample_prod.index.tolist()
    sum_reads = sample_df['count (reads)'].sum()
    sum_reads_prod = sample_df_prod['count (reads)'].sum()
    sum_reads_non_prod = sample_df_non_prod['count (reads)'].sum()
    
    
    pie=ax.pie(sample_prod_list, shadow=True, startangle=90,labels=labels, colors=('black','red','pink'))
    ax.set_title('Total Reads:  %(total)s \n Productive Reads:  %(prod)s \n Non-Productive Reads: %(non_prod)s ' %{'total': sum_reads, 'prod': sum_reads_prod, 'non_prod': sum_reads_non_prod}, fontsize=7)
    ax.legend(legend,loc='best',fontsize=7)
    ax.set_aspect('equal')
    
   
    return sum_reads, sum_reads_prod, sum_reads_non_prod, pie
    print 'finish pie graph'

  

#----------------------------------------------------------------------------------------------------
## this function calculates data for generating a histogram of clonality frequencies.
## clonality= the fraction of a specific sequence out of the total number of reads or templates. 
## this function also calculated cut off values for seperating high and low resolution graphs, and 
##maximal x value for the graphs.  
## *** note that the clonality lists are declared as global. can be changed so they would be returned***


def seq_freq_calc(sample_df, sample_df_prod, sample_df_non_prod, sample_name):
    
    print 'sequence frequency calculations...'
    global read_freq_prod, read_freq_non_prod, read_freq_aa, sum_prod, sum_non_prod
##generating a list of frequencies for each sequence in the df. the frequencies are calculated as percent of
## total reads in the relevant 'group' (productive or non-productive sequences, while aa reads are identical
##to productive reads:
    read_freq_prod = list((sample_df_prod['count (reads)']/sample_df_prod['count (reads)'].sum())*100)
    read_freq_non_prod = list((sample_df_non_prod['count (reads)']/sample_df_non_prod['count (reads)'].sum())*100)        
    read_freq_aa_df = sample_df_prod.groupby(['aminoAcid'])[['count (reads)']].sum()
    read_freq_aa = list((read_freq_aa_df['count (reads)']/sample_df_prod['count (reads)'].sum())*100)
## the following print output validates that each list sum up to 100%
    print ('total x prod=%s' %sum(read_freq_prod))
    print ('total x non prod=%s' %sum(read_freq_non_prod))
    print ('total x aa=%s' %sum(read_freq_aa))
## calculate some optional x cutoff values:
## (1) the x_max for the low resolution hist (the 'ceiling' of the largest value in all the frequency lists): 
    x_max_prod=max(read_freq_prod)
    x_max_non_prod=max(read_freq_non_prod)
    x_max_aa=max(read_freq_aa)
    x_max=max(x_max_prod,x_max_non_prod, x_max_aa)
    x_max=round(x_max,1)+0.1
    print x_max
                               
  
##calculate the x_cut between the high and low resolution histogrtams (using my function: percentile_cut_off):
    x_cut=percentile_cut_off([read_freq_prod, read_freq_non_prod, read_freq_aa],cut_percentile=99.5, round_fold=2)
    if x_cut<0.01:
        x_cut=0.01
    print x_cut
##calculate the x value to end the cumulative function(using my function: percentile_cut_off):
    x_max_cumul=percentile_cut_off([read_freq_prod, read_freq_non_prod, read_freq_aa],cut_percentile=98.5, round_fold=2)
    if x_max_cumul<0.01:
        x_max_cumul=0.01
    print x_max_cumul
    
##calculate the number of unique productive and non productive sequences
    sum_prod=len(sample_df_prod['count (reads)'])
    sum_non_prod=len(sample_df_non_prod['count (reads)'])
    print'finish sequence frequency calculations...'
    return x_max, x_cut, x_max_cumul

#--------------------------------------------------------------------------------
## this function plots the clonality histogram based on the calculations in the seq_freq_calc
## I use to call this function twice - for low and high resolution graphs.  

def seq_freq_hists(ax, resolution, n_bins, x_min, x_max, yscale, y_max=None):
    '''
    this function uses the lists of sequences frequencies generated by the seq_freq_calc function, 
    and the cutoff values, to generate histograms showing the distribution of those frequencies in the sample
    input:
    ax=axes parameters from main function
    resolution-'high', 'low' or 'None' [string]
    n_bins= number of bins to use in the hist [int]
    x_min, x_max= minimum and maximal values of the x-axis shown (not calculated!) [int/float]
    yscale= 'linear' or 'log' [string]
    y_max - maximum value of y-axis to show (not calculated!) [int/float]
    '''
    
    global weights1, weights2, weights3
    print 'generating clonality distribution...'
    
## plotting the histograms:  
    x_eval=percentile_cut_off([read_freq_prod, read_freq_non_prod, read_freq_aa],cut_percentile=99.999999, round_fold=1)
    if resolution=='high':
        if x_eval>9:
            n_bins=n_bins*1.5
        elif x_eval<5 and x_eval>=3:
            n_bins=n_bins/2
        elif x_eval<3 and x_eval>=1:
            n_bins=n_bins/4
        elif x_eval<1:
            n_bins=n_bins/8
    print ('n_bins=%s' %n_bins)
    
    if ax is None:
        ax = plt.gca()
## caculating the weights to normalize each value in the read frequency lists, in order to acheive normalized
##histogram which is not affected by the total number of unique sequences. the frequencies are normalized to 
## the number of unique sequences.      
    weights1 = np.ones_like(read_freq_prod)/len(read_freq_prod) ## actually 1/len(list)
    weights2 = np.ones_like(read_freq_non_prod)/len(read_freq_non_prod)
    weights3 = np.ones_like(read_freq_aa)/len(read_freq_aa)
 ## plot the three normalized histograms on the same plot
 ##the histogram calculation is based on the number of bins, and the weights   
    x,bins,p=ax.hist((read_freq_prod, read_freq_non_prod,read_freq_aa) , n_bins, weights=(weights1, weights2, weights3), label=('Productive','Non Productive', 'Amino Acids') , color=('black','red','green'))
    ax.set_xlim(x_min,x_max)
#     ax.set_ylim(0, y_max)
    ax.set_yscale(yscale)
    if resolution=='None':
        ax.set_title('Sequence Abundance Distribution', fontsize=10)
    else:
        ax.set_title('Sequence Abundance Distribution \n %s resolution' %resolution.capitalize(), fontsize=10)
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlabel('Sequence Clonality (% of reads)', fontsize=9)
    ax.set_ylabel('% of Unique Sequences')
    ax.margins(0.02,0.02)
    print 'finish generating clonality distribution...'
    return x

        
def seq_freq_cumul(ax, n_bins_cumul, x_min, x_max, y_min, y_max):   
    '''
    this function uses the lists of sequences frequencies generated by the seq_freq_calc function, 
    and the cutoff values, to generate cumulative histogram of these frequencies.
    input:
    ax=axes parameters from main function
    n_bins_cumul= number of bins to calculate the hist [int]
    x_min, x_max, ymin, ymax= minimum and maximal values of the x and y--axes shown (not calculated!) [int/float]   
    '''   
    
    print 'generating cumulative clonality distribution...'
## calculating bins data for the cumulative plot:    
    prod_hist, prod_e= np.histogram(read_freq_prod, bins=5000, weights=weights1)
    non_prod_hist, non_prod_e= np.histogram(read_freq_non_prod,bins=5000, weights=weights2)
    aa_hist, aa_e= np.histogram(read_freq_aa,bins=5000, weights=weights3)
#     
# ##now generate the cumulative distribution from this data:
    prod_hist=list(prod_hist)
    non_prod_hist=list(non_prod_hist)
    aa_hist=list(aa_hist)  
    prod_hist.insert(0,0.0) ## add the value '0' to the begining of each list, so th ecumulative function 
                            ## begins at 0:
    non_prod_hist.insert(0,0.0)
    aa_hist.insert(0,0.0)
    prod_e=np.insert(prod_e,0, 0.0) ## add the value '0' to the begining of each list, so th ecumulative function 
                            ## begins at 0:
    non_prod_e=np.insert(non_prod_e,0, 0.0)
    aa_e=np.insert(aa_e,0, 0.0)
    
    cumulative_prod = np.cumsum(prod_hist)
    cumulative_non_prod = np.cumsum(non_prod_hist)
    cumulative_aa = np.cumsum(aa_hist)
    print('the first x values in the prod_cumul list are: %s' %prod_e[0:20])
    print('the first y values in the prod_cumul list are: %s' %cumulative_prod[0:20])
    
    print('the first x values in the non_prod_cumul list are: %s' %non_prod_e[0:20])
    print('the first y values in the non_prod_cumul list are: %s' %cumulative_non_prod[0:20])
    
 ## calcculate the KS-test2:
    
    
    KS2_Prod_NonProd, p_Prod_NonProd=stats.ks_2samp(read_freq_prod, read_freq_non_prod)
    KS2_Prod_AA, p_Prod_AA=stats.ks_2samp(read_freq_prod, read_freq_aa)
    

# ##plotting the cumulative distributions:               
    if ax is None:
        ax = plt.gca()
    plot1=ax.plot(prod_e[:-1], cumulative_prod , label='Productive', alpha=0.5, color='black')
    plot2=ax.plot(non_prod_e[:-1], cumulative_non_prod, label='Non Productive', alpha=0.5,color='red')
    plot3=ax.plot(aa_e[:-1], cumulative_aa , label='Amino Acids' , alpha=0.5, color='green')
#     
    ax.text(0.005, 0.4, ' Prod vs. Non-prod, p=%s \n Prod vs. AA, p=%s \n (KS test)' %(p_Prod_NonProd, p_Prod_AA), bbox=dict(facecolor='white', alpha=0.5))
    ax.set_xlim(0,x_max_cumul)
    ax.set_title('Cumulative Sequence Frequencies')
#     ax.set_ylim(0,1.2)
    ax.legend(loc='lower right', fontsize=9)
    ax.set_xlabel('Sequence Clonality (% of reads)', fontsize=9)
    ax.set_ylabel('cumulative Percent \n of Unique Sequences')
    ax.margins(0.02,0.02)
    print 'Finished generating cumulative clonality distribution...'
    return plot1, plot2, plot3

#-----------------------------------------------------------------
'''
the following function plots a bar graph for a specific categorial sample variable (a column in the df, such
as 'cdr3Length' 'jFamilyName' etc.
input:
ax: axes parameters from the calling function.
var: the name of the variable (exact column name)
count=if the function is called within a for loop, the count should be zero before the loop starts, and
incremented after each call to the loop. this argument is used only to indicate that the legend should
appear only in the first subplot of a figure, and can be removed from the argument list if the last lines 
are deleted.
sample_df_prod, sample_df_non_prod- dfs of productive and non-productive sequences of a sample, should be 
generated by the 'get_sample_data' function for each sample, or loaded from pickles. this part of the function 
should be generalizaed to enable comparison of different samples, regardless of their nature or number.
sample name: string of sample name without file extension.
fontsizelabels: fontsize for the xtick labels - should be small for variables with many categories and 
vice versa. 
'''    
    

        
def plot_var_bar(ax,var, count, sample_df_prod, sample_df_non_prod, sample_name,fontsizelabels):        
        print '%s' % var
        print count
        print ax
## calculate how many sequences (measured by the total number of their reads) are in each category of 
## the variable, for example: how many sequences have cdr in lengths of 39, 42,45,48 bps etc. 
        sum_var_prod = sample_df_prod.groupby(['%s' % var])[['count (reads)']].sum()
        sum_var_non_prod = sample_df_non_prod.groupby(['%s' % var])[['count (reads)']].sum()
## normalize the number of occurences to fractions by dividing by the total sum of reads:
        sum_var_prod['weighted reads']=sum_var_prod['count (reads)']/(sum_var_prod['count (reads)'].sum())
        sum_var_non_prod['weighted reads']=sum_var_non_prod['count (reads)']/(sum_var_non_prod['count (reads)'].sum())

##generate a combined list of all the categories found in the prod and non prod data.  
## the list is generated a little different if the categories are strings or integrers.        
        varlist_prod=sum_var_prod.index.tolist() 
        varlist_non_prod=sum_var_non_prod.index.tolist()
        print type(varlist_prod[0])
        if isinstance(varlist_prod[0], str):
            print 'string'
            varlist_total=list(set(varlist_prod+varlist_non_prod))
        else:
            print 'not string'
            max_varlist=max(varlist_prod+varlist_non_prod)
            print max_varlist
            varlist_total=range(max_varlist+1)
## generate lists of normalized fraction of reads per category in either prod or non-prod seqeunces
## based on the existing dfs. 
        weighted_reads_prod=[]
        weighted_reads_non_prod=[]
        for i in varlist_total:
            if i in varlist_prod and i in varlist_non_prod:
                weighted_reads_prod.append(sum_var_prod.loc[i,'weighted reads'])
                weighted_reads_non_prod.append(sum_var_non_prod.loc[i,'weighted reads'])
            elif i in varlist_prod and i not in varlist_non_prod:
                weighted_reads_prod.append(sum_var_prod.loc[i,'weighted reads'])
                weighted_reads_non_prod.append(0.0)
            elif i not in varlist_prod and i in varlist_non_prod:
                weighted_reads_prod.append(0.0)
                weighted_reads_non_prod.append(sum_var_non_prod.loc[i,'weighted reads'])
            else:
                weighted_reads_prod.append(0.0)
                weighted_reads_non_prod.append(0.0)
                
        print('total weighted reads in prod=%s' %sum(weighted_reads_prod))
        print('total weighted reads in non prod=%s' %sum(weighted_reads_non_prod))
        print weighted_reads_prod[0:10]
        print weighted_reads_non_prod[0:10]
 
 ## for categories that are not strings, calculate the ks and t- tests to compare productive and
 ##non--productive data:
        if not isinstance(varlist_prod[0], str):
            ks_s,ks_p=stats.ks_2samp(sample_df_prod[var], sample_df_non_prod[var])
            print('KS_p_value=%s' %ks_p)
            t_s,t_p=stats.ttest_ind(sample_df_prod[var], sample_df_non_prod[var])
            print('ttest_p_value=%s' %t_p)
            if ks_p<=10**-6:
                text_ks='KS_p_value<10^-6'
            else:
                text_ks='KS_p_value=%s' %round(ks_p,6)
            if t_p<=10**-6:
                text_t='ttest_p_value<10^-6'
            else:
                text_t='ttest_p_value=%s' %round(t_p,6)
            ax.annotate(text_ks+'\n'+text_t, xy=(0.95, 0.95), xycoords='axes fraction', fontsize=8,
                horizontalalignment='right', verticalalignment='top', fontweight='bold')
            
 ## calculate parameters for the plot, including bar positions, tick positions and labels.           
        space = 0.3
        width = (1 - space) / 2
        if isinstance(varlist_prod[0], str):
            pos_prod = [i - 0.5 * width for i in range(len(varlist_total))]
            pos_non_prod = [i + 0.5 * width for i in range(len(varlist_total))]
            print varlist_total
            print len(varlist_total)
            xtickspos=range(len(varlist_total))
            xtickslab=varlist_total
            rot=90
            
        else:
            pos_prod = [i - 0.5 * width for i in varlist_total]
            pos_non_prod = [i + 0.5 * width for i in varlist_total]
            xtickspos=[i for i in range(len(varlist_total)) if i%3==0]
            xtickslab=xtickspos
            rot=0
           
        print 'showing bar graph'
        bar1=ax.bar(pos_prod,weighted_reads_prod, align='center', label='Productive', width=width, color='black', edgecolor='None')
        bar2=ax.bar(pos_non_prod, weighted_reads_non_prod, align='center', label='Non-Productive', width=width, color='red', edgecolor='None')
        ax.set_xticks(xtickspos)
        #ax.set_yticks(yticks)
        ax.set_xticklabels(xtickslab, rotation=rot, fontsize=fontsizelabels)
        #ax.set_yticklabels(yticks, fontsize=8)
        ax.set_title('%s (bps)' % var, fontsize=9)
        ax.set_ylabel("Frequency (%)", fontsize=8)
        print 'now set yscale'
        ax.set_yscale('log')
        if not isinstance(varlist_prod[0], str):
            ax.set_xlim(-1.5*width,max_varlist+1)
        
        if count==0:
            ax.legend(loc='upper left', fontsize=8) 
        return bar1,bar2

#-----------------------------------------------------------------------------------------------
## this function generates the figure and general outline for categorial parameter histograms and call the
## plot_var_bar function for generating each plot.
 

def categorial_par_hist (sample_df_prod, sample_df_non_prod, sample_name):
    
## replace all 'unresolved' values in the dfs and generate new parameters for combine v and j usage (genes 
## and families. 
    print 'working on length parameters histograms...'
    sample_df_prod=sample_df_prod.replace('unresolved', np.nan)
    sample_df_non_prod=sample_df_non_prod.replace('unresolved', np.nan)
    sample_df_prod['v-j']=sample_df_prod['vGeneName']+'_'+sample_df_prod['jGeneName']
    sample_df_non_prod['v-j']=sample_df_non_prod['vGeneName']+'_'+sample_df_non_prod['jGeneName']
    sample_df_prod['v-j-family']=sample_df_prod['vFamilyName']+'_'+sample_df_prod['jFamilyName']
    sample_df_non_prod['v-j-family']=sample_df_non_prod['vFamilyName']+'_'+sample_df_non_prod['jFamilyName']
    
    var_list = ['cdr3Length', 'vDeletion', 'n1Insertion', 'd5Deletion', 'd3Deletion', 'n2Insertion', 'jDeletion','vGeneName','vFamilyName', 'jGeneName', 'jFamilyName','dGeneName', 'v-j','v-j-family']
    
    fig2=plt.figure(figsize=(8, 10))    
    ax1 = plt.subplot2grid((7,1), (0,0)) 
    ax2 = plt.subplot2grid((7,1), (1,0))
    ax3 = plt.subplot2grid((7,1), (2,0))
    ax4 = plt.subplot2grid((7,1), (3,0))
    ax5 = plt.subplot2grid((7,1), (4,0))
    ax6 = plt.subplot2grid((7,1), (5,0))
    ax7 = plt.subplot2grid((7,1), (6,0))
    print 'start working on length parameters...'
    
#     ((ax1), (ax2),(ax3), (ax4),(ax5),(ax6),(ax7)) = plt.subplots(7, 1, sharey=False, figsize=(8,10))
    ax_list1=[ax1,ax2,ax3,ax4,ax5,ax6,ax7]
    plt.suptitle('CDR3 Length Parameters - Sample %s' % sample_name, fontsize=12, fontweight='bold')
    
    count=0
    for var in var_list[0:7]:
        plot_var_bar(ax_list1[count],var, count, sample_df_prod, sample_df_non_prod, sample_name,8)
        count=count+1
    plt.subplots_adjust(left=0.12,bottom=0.04, right=0.97, top=0.90, wspace=0.28,hspace=0.60)
    print 'finished working on length parameters'
  
    print 'start working on Gene Usage Distributions...'
    fig3=plt.figure(figsize=(8, 10))    
    ax1 = plt.subplot2grid((5,1), (0,0)) 
    ax2 = plt.subplot2grid((5,1), (1,0))
    ax3 = plt.subplot2grid((5,1), (2,0))
    ax4 = plt.subplot2grid((5,1), (3,0))
    ax5 = plt.subplot2grid((5,1), (4,0))  
    ax_list1=[ax1,ax2,ax3,ax4,ax5,ax6,ax7]
    plt.suptitle('Gene Usage Distributions - Sample %s' % sample_name, fontsize=12, fontweight='bold')
    count=0
    for var in var_list[7:12]:
        print ax_list1[count]
        plot_var_bar(ax_list1[count],var, count, sample_df_prod, sample_df_non_prod, sample_name,8)
        count=count+1
    plt.subplots_adjust(left=0.12,bottom=0.10, right=0.97, top=0.90, wspace=0.28,hspace=0.99)
    
    fig4=plt.figure(figsize=(11, 6))    
    ax1 = plt.subplot2grid((2,1), (0,0)) 
    ax2 = plt.subplot2grid((2,1), (1,0))
    ax_list1=[ax1,ax2,ax3,ax4,ax5,ax6,ax7]
    plt.suptitle('Combined Gene Usage Distributions - Sample %s' % sample_name, fontsize=12, fontweight='bold')
    count=0
    for var in var_list[12:14]:
        print ax_list1[count]
        plot_var_bar(ax_list1[count],var, count, sample_df_prod, sample_df_non_prod, sample_name,5)
        count=count+1
    plt.subplots_adjust(left=0.10,bottom=0.13, right=0.97, top=0.90, wspace=0.28,hspace=0.97)
    
    
    
    
    print 'finished working on Gene Usage Distributions...'
    return fig2, fig3, fig4


def gc_content(ax, sample_df, sample_df_prod, sample_name):
    
## this function calculates the gc content of all seqeunces in the dfs (unique sequneces, regardless of their
## number of reads/templates) and plots an histogram:

    print 'start working on gc graph...'     
    seq_prod=list(sample_df_prod['nucleotide'])
    gc_values_prod = sorted(GC(seqp) for seqp in seq_prod)
    seq_non_prod=list(sample_df_non_prod['nucleotide'])
    gc_values_non_prod = sorted(GC(seqn) for seqn in seq_non_prod)
    weights_prod=(np.ones_like(seq_prod, dtype=np.float))/len(seq_prod)
    weights_non_prod=(np.ones_like(seq_non_prod, dtype=np.float))/len(seq_non_prod)
    hist1=ax.hist([gc_values_prod,gc_values_non_prod], weights=[weights_prod, weights_non_prod], label=['Productive','Non-Productive'], color=['black','red'])
    ax.set_title("sample name: %s\n%i unique productive sequences(GC%% %0.1f to %0.1f)\n%i unique non productive sequences(GC%% %0.1f to %0.1f)" \
            % (sample_name, len(gc_values_prod),min(gc_values_prod),max(gc_values_prod),len(gc_values_non_prod),min(gc_values_non_prod),max(gc_values_non_prod)))
    ax.set_xlabel("GC%")
    ax.set_ylabel("Frequency")
    ax.legend()
    return hist1
    print 'finish working on gc graph...'     
#---------------------------------------------------------------------------------------------------------------------------------
## this function calculates the usage frequency of each amino acid in all seqeunces in the dfs (unique
## sequneces, regardless of their number of reads/templates) and plots an histogram:
    
def aa_frequencies (ax, sample_df, sample_df_prod, sample_name): 
    print 'start working on aa graph...'                                                  
    aa_prod=list(sample_df_prod['aminoAcid'])
    aa_non_prod=list(sample_df_non_prod['aminoAcid'])
    aa_list=["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    
## generate lists of aminoacid counters:     
    aa_prod_count=[0]*len(aa_list)
    for aap in aa_prod:
        aa_ind=0
        for aa in aa_list:
            aa_prod_count[aa_ind]+=aap.count(aa)
            aa_ind+=1
    print aa_prod_count
    aa_non_prod_count=[0]*len(aa_list)
    for aap in aa_non_prod:
        if isinstance(aap, basestring):
            aa_ind=0
            for aa in aa_list:
                aa_non_prod_count[aa_ind]+=aap.count(aa)
                aa_ind+=1
    print aa_non_prod_count
## normalize the counts to get fractions:
    total_aa_prod=sum(aa_prod_count)
    norm_aa_prod_count=[float(i)/total_aa_prod for i in aa_prod_count]
    total_aa_non_prod=sum(aa_non_prod_count)
    norm_aa_non_prod_count=[float(i)/total_aa_non_prod for i in aa_non_prod_count]
    
## calculate parameters for plotting the histograms, including bar positions, tick positions and labels:    
    space = 0.3
    width = (1 - space) / 2
    pos_prod = [i - 0.5 * width for i in range(len(aa_list))]
    pos_non_prod = [i + 0.5 * width for i in range(len(aa_list))]
    xtickspos=range(len(aa_list))
    xtickslab=aa_list
##plot the hitograms: 
    bar1=ax.bar(pos_prod,norm_aa_prod_count, align='center', label='Productive', width=width, color='black', edgecolor='None')
    bar2=ax.bar(pos_non_prod, norm_aa_non_prod_count, align='center', label='Non-Productive', width=width, color='red', edgecolor='None')
    ax.set_xticks(xtickspos)
    ax.set_xticklabels(xtickslab, fontsize=8)
    ax.set_title('AminoAcid', fontsize=14)
    ax.set_ylabel("Frequency (%)", fontsize=8)
    ax.legend(loc='upper left', fontsize=8) 
    return bar1,bar2
    print 'finish working on aa graph...'     


'''
main function:
'''    
# # uncomment the following rows for using more than one sample:

onlyfiles = [f for f in listdir("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/CSVfiles") if isfile(join("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/CSVfiles", f))]
onlyfiles = [datafile for datafile in onlyfiles if datafile.startswith ('HIP') and datafile.endswith('.csv')]
print onlyfiles
countfiles=0
## insert 'get headers' here
for datafile in onlyfiles[0:6]:
    print datafile
    countfiles=countfiles+1
    global sample_name
    sample_name=re.sub('.csv', '', datafile)
    print sample_name
    print countfiles

# # uncomment from here to specify single sample:
# sample_name='HIP00110' ## use instead of reading the file names from the directory
## from here on, relevant to any number of samples:    
    
    global create_pdf, figlist
    
    figlist=[] ## generate figlist for pdf creation
    create_pdf=True
    
    
## if need to generate dfs:
#    sample_df, sample_df_prod, sample_df_non_prod = get_sample_data(sample_name)## generate dataframe from
                                                                                ## the sample TSV file, and also
                                                                                ## sub-dfs for productive and non
                                                                                ## productive sequences alone.                              ## if dfs should be loaded:
    sample_df= pickle.load( open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_%s' %sample_name,"rb" ))
    sample_df_prod= pickle.load( open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_prod_%s' %sample_name,"rb" ))
    sample_df_non_prod= pickle.load( open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_non_prod_%s' %sample_name,"rb" ))
    # get_headers(sample_name, sample_df) #this function is currently not in use!
    
    #################
    ## first figure##
    #################
    fig1=plt.figure(figsize=(10,8))
      
    ## the following lines define the location and size of each subplot in the figure, the genral format is 3 lines
    ## and two columns (3,2). then the specific location within this grid is defined (now 1 is 0, 2 is 1 and so on)
    ## then the width of the subplot is defined in terms of how many columns it spans, when the default is one. 
    ax1 = plt.subplot2grid((3,2), (0,0), colspan=1) 
    ax2 = plt.subplot2grid((3,2), (0,1), colspan=1)
    ax3 = plt.subplot2grid((3,2), (1, 0), colspan=2)
    ax4 = plt.subplot2grid((3,2), (2, 0))
    ax5 = plt.subplot2grid((3,2), (2, 1))
      
    #  ((ax1, ax2), (ax3, ax4), (ax5,ax6)) = plt.subplots(3, 2, sharey=False, figsize=(10,8))
    plt.suptitle ('Sample %s - Sample Overview (Fig1)' %sample_name, fontsize=16)
      
    sum_reads, sum_reads_prod, sum_reads_non_prod, prod_pie_fig=prod_dist(ax5, sample_df, sample_df_prod, sample_df_non_prod, sample_name)
    x_max, x_cut, x_max_cumul = seq_freq_calc(sample_df, sample_df_prod, sample_df_non_prod, sample_name)
    seq_freq_hists(ax=ax1, resolution='high', n_bins=4000, x_min=0, x_max=x_cut, yscale='log')
    seq_freq_hists(ax=ax2, resolution='low', n_bins=30, x_min=0, x_max=x_max, yscale='log')
    seq_freq_cumul(ax=ax3, n_bins_cumul=5000, x_min=0, x_max=x_max_cumul, y_min=0, y_max=1.2)      
    samp_size,max_samp_size, n_samp_seq_all, ymax_rare=rarefaction_calc(sample_name, [sample_df_prod, sample_df_non_prod],column_name='count (reads)')   
    rarefaction_plot(ax=ax4, sample_name=sample_name, data=n_samp_seq_all, names=['Productive', 'Non-Productive'], colors=['black', 'red'], samp_size=samp_size, max_samp_size=max_samp_size,ymax_rare=ymax_rare)   
    plt.subplots_adjust(left=0.09,bottom=0.11, right=0.95, top=0.89, wspace=0.24,hspace=0.43)
    figlist.append(fig1)

##############################################################
## 2nd, 3rd and 4th figures (length parameters, gene usage): ##
##############################################################    
    
    fig2, fig3, fig4=categorial_par_hist(sample_df_prod, sample_df_non_prod, sample_name)
    figlist.append(fig2)
    figlist.append(fig3)
    figlist.append(fig4)
    


    #################
    ## figure 5th:##
    #################

    fig5=plt.figure(figsize=(10,8))
    ax1 = plt.subplot2grid((2,1), (0,0), colspan=1) 
    ax2 = plt.subplot2grid((2,1), (1,0), colspan=1)
    
    gc_content(ax1, sample_df, sample_df_prod, sample_name)
    aa_frequencies (ax2, sample_df, sample_df_prod, sample_name)
    plt.subplots_adjust(left=0.09,bottom=0.11, right=0.95, top=0.89, wspace=0.24,hspace=0.43)
    figlist.append(fig5)
   

    if create_pdf:
        with PdfPages('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Images/Sample_View_%s.pdf' %sample_name) as pdf:
            for fig in figlist:
                pdf.savefig(fig)
        pdf.close
    else:
        plt.show()


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## olf cunctions and call for the functions:        


def get_headers(sample_name, sample_df): 
## this function is currently not in use!!
## uncomment to use existing df file:    
# df=Load ('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/DFfiles/%s_data_updated.df' %sample_name)  
    global head_list, head_dic
    head_list = list(sample_df.columns.values)
    head_dic = {}
    for col in head_list:
        col_id = head_list.index('%s' % col)
        head_dic[col] = col_id


        
## call for old functions:        
    # seq_hists(sample_df, sample_df_prod, sample_name)
    # gene_usage_hist (sample_df, sample_df_prod, 'HIP00110')
    # aa_frequencies (sample_df, sample_df_prod, 'HIP00110')
        
    