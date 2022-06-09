from os import listdir
from os.path import isfile, join, isdir, exists
import pandas as pd
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot,draw_correlation_scatter
from matplotlib.ticker import FormatStrFormatter

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
# import cPickle as pickle
from Bio.SeqUtils import GC
import seaborn as sns
import random
import math



print('done1')

def concat_summarizing_dfs(dfs_folder):
    filenames = [f for f in listdir(dfs_folder) if isfile(join(dfs_folder, f))]
    print ('number of dfs in directory: %s' %len(filenames))
    
    df_list=[]
    for file_name in filenames:
        with open('%s/%s' %(dfs_folder, file_name), 'rb') as f:
            current_df=pd.read_pickle(f)
        f.close()
        df_list.append(current_df)
    print ('the length of df list is %s' %len(df_list))
           
    df_all=pd.concat(df_list)
           
    return df_all

def concat_summarizing_dfs_excel(dfs_folder):
    filenames = [f for f in listdir(dfs_folder) if isfile(join(dfs_folder, f))]
    print ('number of dfs in directory: %s' %len(filenames))
    
    df_list=[]
    for file_name in filenames:
        try:
            with open('%s/%s' %(dfs_folder, file_name), 'rb') as f:
                current_df=pd.read_excel(f)
            f.close()
            df_list.append(current_df)
        except:
            with open('%s/%s' %(dfs_folder, file_name), 'rb') as f:
                current_df=pd.read_pickle(f)
            f.close()
            df_list.append(current_df)
            
    print ('the length of df list is %s' %len(df_list))
           
    df_all=pd.concat(df_list)
           
    return df_all


def MyPearsonr(x,y):
    import numpy as np
    from scipy.stats import pearsonr
    x=list(x)
    y=list(y)
    nx=np.isnan(x)
    ny=np.isnan(y)
    n=nx+ny
    newx=[x[i] for i in range(len(x)) if not n[i]]
    newy=[y[i] for i in range(len(y)) if not n[i]]
    r,p = pearsonr(newx,newy)
    
    return r,p

def MySpearmanr(x,y):
    import numpy as np
    from scipy.stats import spearmanr
    x=list(x)
    y=list(y)
    nx=np.isnan(x)
    ny=np.isnan(y)
    n=nx+ny
    newx=[x[i] for i in range(len(x)) if not n[i]]
    newy=[y[i] for i in range(len(y)) if not n[i]]
    r,p = spearmanr(newx,newy)
    
    return r,p


def plot_bestFitLine(x,y,ax,color):
    if ax==None:
        fig,ax=plt.subplots()
    nx=np.isnan(x)
    ny=np.isnan(y)
    n=nx+ny
    newx=list(x[~n])
    newy=list(y[~n])
    ax.plot(np.unique(newx), np.poly1d(np.polyfit(newx, newy, 1))(np.unique(newx)),color=color)
#-----------------------------------------------------------------

#the following function takes df, a categorial (category) column within it and a continuous varaible 
#within it (PredictionR), plot boxplot and calculate ttest/ANOVA to see the category effect on the 
#continous variable

def check_category_effect_on_PredictionR(df,category,PredictionR,ax):
    fig = plt.figure()
    if ax==None:
        fig,ax=plt.subplots()
            
    groups = df.groupby(category)
    
    from scipy.stats import ttest_ind, f_oneway

    data_list=[]
    name_list=[]
    ticks_locations=range(1,len(groups)+1)
    for name,group in groups:
        if name!='Unknown':
            data=group[PredictionR]
            data_list.append(data)
            name_list.append(name)
    ticks_locations=range(1,len(name_list)+1)
    ax.boxplot(data_list)
    ax.set_xticks(ticks_locations)
    ax.set_xticklabels(name_list,rotation=90)
    ax.set_xlabel(category)
    ax.set_ylabel(PredictionR)

    #get locations for p-value text
    ylim=ax.get_ylim()
    ypos=ylim[1]
    xlim=ax.get_xlim()
    xpos=xlim[0]
    
    ax.set_title('%s effect on %s' %(category,PredictionR),fontsize=14)


#if the number of groups=2, conduct t-test, if more then 2 conduct 1-way ANOVA. in any case, print the p-value on the graph:
    if len(data_list)==2:
        s,p=ttest_ind(data_list[0],data_list[1])
        ax.text(xpos,ypos,"T test p-value=%.5f" %p,  verticalalignment = 'top', ha = 'left',fontsize=14,color='red')
    elif len(data_list)>2:
        s,p=f_oneway(*data_list)
        ax.text(xpos,ypos,"1way ANOVA p-value=%.5f" %p,  verticalalignment = 'top', ha = 'left',fontsize=14,color='red')
    else:
        print 'data_list is shorter than 2'
        s=0
        p=1
        
        
    return fig,ax,s,p


#--------------------------------------

#the following function takes df, and two continuous  varaibles (Tag, PredictionR) within it, plot their scatter plot
#and calculate the pearson correlation between them:

def check_correlation_PredictionToTag(df,Tag,PredictionR,ax):
    
    if ax==None:
        fig,ax=plt.subplots()
    
    
    x=df[Tag]
    y=df[PredictionR]
   
    ymean=df[PredictionR].mean()                                                                                                                      
        
    ax.scatter(x,y)
    ax.set_xlabel(Tag)
    ax.set_xlabel(Tag,fontsize=12)
    ax.set_ylabel(PredictionR,fontsize=12)
    ylim=ax.get_ylim()
    ypos=ylim[1]
    xlim=ax.get_xlim()
    xpos=xlim[0]
    
    
    # fit with np.polyfit
    nx=np.isnan(x)
    ny=np.isnan(y)
    n=nx+ny
    newx=list(x[~n])
    newy=list(y[~n])
    ax.plot(np.unique(newx), np.poly1d(np.polyfit(newx, newy, 1))(np.unique(newx)),c='black')
    ax.plot(x=x,c='black', linewidth=4)
    ax.set_title('Correl_%s_to_%s' %(Tag,PredictionR),fontsize=14)
                                                                                                                        
                                                                                                                       
    from scipy.stats import pearsonr
    r,p = pearsonr(newx,newy)
    
    ax.text(xpos,ypos,"r=%.4f p=%.6f" %(r,p),  verticalalignment = 'top', ha = 'left',fontsize=14,color='red')
    #transform=ax.transAxes,                                                                                                                 
                                                                                                                        
    return fig,ax,r,p

#-----------------------------------------------------------------------------

# this function plot a scatter plot, color the dots according to the  specified category, and calculate weather the ratios between x and y (transformed, 
# see below) is significantly different between the different groups
# the ratio calculation is - transform each variable to z scores, multiply each zscore by 10 and add 100, calculate x_var transformed z_score/y_var transformed zscore



def plot_scatter_with_category_colors(df,x_var,y_var,category,ax=None,figsize=(3,3)):
    df=df.dropna() #get rid of all na values in df
    from scipy.stats.mstats import zscore
    from scipy.stats import ttest_ind, f_oneway

    #calculate transformed zscores and their ratio:
    df['%s_zscore' %x_var]=zscore(df[x_var])*10+100
    df['%s_zscore' %y_var]=zscore(df[y_var])*10+100
    df['%s_%s_zscore ratio' %(x_var, y_var)]=df['%s_zscore' %x_var]/df['%s_zscore' %y_var]
    correl_ratios_normed_column='%s_%s_zscore ratio' %(x_var, y_var)

    fig = plt.figure()
    if ax is None:
        ax = fig.add_subplot(111)

    groups = df.groupby(category) #divide data into groups according to the chosen category

    color_list = plt.cm.Set1(np.linspace(0, 1, 20)) #choose cmap for coloring the different groups
    count=0


    #plot each group separately, divide the transformed z score ratios into groups
    data_list=[]
    name_list=[]
    for name,group in groups:
        ax.plot(group[x_var], group[y_var], marker='o', linestyle='', ms=6, label=name,c=color_list[count])
        count+=1
        if name!='Unknown':
            data=group[correl_ratios_normed_column]
            data_list.append(data)
            name_list.append(name)

    #ax.set_xlabel(x_var) #uncomment if only one subplot
    #ax.set_ylabel(y_var) # uncomment if only one subplot

    ax.legend(loc=1)


    # fit with np.polyfit
    x=df[x_var]
    y=df[y_var]
    nx=np.isnan(x)
    ny=np.isnan(y)
    n=nx+ny
    newx=list(x[~n])
    newy=list(y[~n])
    ax.plot(np.unique(newx), np.poly1d(np.polyfit(newx, newy, 1))(np.unique(newx)),c='black')
    ax.plot(x=x,c='black', linewidth=4)

    #get locations for the p-value text
    ylim=ax.get_ylim()
    ypos=ylim[1]
    xlim=ax.get_xlim()
    xpos=xlim[0]

    ax.set_title('%s' %category)


    #if the number of groups=2, conduct t-test, if more then 2 conduct 1-way ANOVA. in any case, print the p-value on the graph:
    if len(data_list)==2:
        s,p=ttest_ind(data_list[0],data_list[1])
        ax.text(xpos,ypos,"T test p-value=%.5f" %p,  verticalalignment = 'top', ha = 'left',fontsize=14,color='red')

    elif len(data_list)>2:
        s,p=f_oneway(*data_list)
        ax.text(xpos,ypos,"1way ANOVA p-value=%.5f" %p,  verticalalignment = 'top', ha = 'left',fontsize=14,color='red')
    
    return fig,ax,s,p

#------------------------------------------------------------------

'''
This function generates scatter plot comparing the seq frequencies of two samples. sequences can be nucleotides (seqType-'nt;) 
or aminoacid (seqType='aa'). the scatter plot can be colored using a specified color (color='blue', 'red' etc) or by density 
(color='density'). sequences that appear in one sample but not in the other are indicated with seq freq that is lower than the 
minimum real frequency in the datasets. 
pearson correlation r and p are also calculated, and the percent of sequences overlapping between sequences. 

this function was copied to 'MyFunctionsShani.py'

input:
*sample_name1, sample_name2 - names of samples to compare
*seqType: 'nt' or 'aa'
*color= 'density' for density map, any other legitimate color for uniform color
*ax
'''

def sample_dup_analysis(sample_name1, sample_name2,seqType,color,ax=None):
    from scipy.stats import gaussian_kde
    fig = plt.figure()
    if ax==None:
        fig,ax=plt.subplots()

    df1=pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/%s.tsv" %sample_name1)
    df2=pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/%s.tsv" %sample_name2)
    column1='%s frequencyCount (perc)' %sample_name1
    column2='%s frequencyCount (perc)' %sample_name2

    if seqType=='nt':
        df1nt=df1[['nucleotide','frequencyCount (%)']]
        df1nt=df1nt.set_index('nucleotide')
        df1nt=df1nt.rename(columns={'frequencyCount (%)':column1})
        df2nt=df2[['nucleotide','frequencyCount (%)']]
        df2nt=df2nt.set_index('nucleotide')
        df2nt=df2nt.rename(columns={'frequencyCount (%)':column2})
        df1df2=pd.merge(df1nt, df2nt, how='outer', left_index=True, right_index=True)
        
    elif seqType=='aa':
        print 'aa'
        df1=pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/%s.tsv" %sample_name1)
        df1aa=df1[['aminoAcid','frequencyCount (%)']]
        df1aagrouped=pd.DataFrame(df1.groupby('aminoAcid')['frequencyCount (%)'].sum())
        df1aagrouped=df1aagrouped.rename(columns={'frequencyCount (%)':column1})
        df2=pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/%s.tsv" %sample_name2)
        df2aa=df2[['aminoAcid','frequencyCount (%)']]
        df2aagrouped=pd.DataFrame(df2.groupby('aminoAcid')['frequencyCount (%)'].sum())
        df2aagrouped=df2aagrouped.rename(columns={'frequencyCount (%)':column2})

        df1df2=pd.merge(df1aagrouped, df2aagrouped, how='outer', left_index=True, right_index=True)

   
    data1=df1df2[column1]
    data2=df1df2[column2]
    
    min_freq=min(data1.min(),data2.min())
    dist = int(math.log10(abs(min_freq)))-1 #calculates the number od zeros after the decimal point in the min_freq value

    n_total_seqs=len(df1df2)
    
    n_overlap=len(df1df2[(df1df2[column1].notnull())&(df1df2[column2].notnull())])
    perc_overlap=round(float(n_overlap)*100/n_total_seqs,2)

   
    data1=data1.fillna(10**dist)
    data2=data2.fillna(10**dist)

    if ax==None:
        fig,ax=plt.subplots()

    if color=='density':
        # Calculate the point density
        xy = np.vstack([data1,data2])
        z = gaussian_kde(xy)(xy)
        cmap=plt.cm.rainbow
        c=z
        bestfitColor='grey'
    else:
        c=color
        cmap=plt.cm.rainbow
        bestfitColor=color
        

    ax.scatter(data1, data2, c=c, cmap=cmap,  alpha=0.2)


    ticks=np.logspace(dist-1,0,-(dist-2))#generates the x,y ticks based on the min_freq. the lower tick is used only for margin purposes and the second tick is used as '0'
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
    ax.set_xlim(5*10**(dist-1),1)
    ax.set_ylim(5*10**(dist-1),1)
#     ax.set_title('NT Sequences Frequency Correlation',fontsize=12)
    ax.set_xlabel(column1 ,fontsize=9)
    ax.set_ylabel(column2,fontsize=9)
    ax.margins(0.2,0.2)
   
#     ax.set_aspect('equal', adjustable='box')
#     cax = fig.add_axes([0.95, 0.2, 0.02, 0.6])
#     ax.add_colormap(cax, cmap=cmap)
    
    
#    
    r,p=MyPearsonr(data1,data2)

    #get locations for p-value text
    ylim=ax.get_ylim()
    ypos=ylim[1]
    xlim=ax.get_xlim()
    xpos=xlim[0]
    print xpos, ypos
    ax.text(xpos,ypos,"r=%.3f p=%.2f, perc overlap=%.2f" % (r, p,perc_overlap),  verticalalignment = 'top', ha = 'left',fontsize=8,color='red')
    

    #plot bestfit and identity curves:
    
    ax.plot(np.unique(data1), np.poly1d(np.polyfit(data1, data2, 1))(np.unique(data1)),c=bestfitColor)
    ax.plot(data1,data1,color='black')
    
    
#     filename='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/QC/dup_comparison_%s_%s_%s_%s' %(sample_name1,sample_name2,seqType,color)
#     fig.savefig(filename=filename, dpi=200)
    print sample_name1, sample_name2
    
    
    ## generate summarizing df:
    summaryDF=pd.DataFrame()
    summaryDF.loc[1,'Sample1']=sample_name1
    summaryDF.loc[1,'Sample2']=sample_name2
    summaryDF.loc[1,'seqType']=seqType
    summaryDF.loc[1,'r']=r
    summaryDF.loc[1,'p']=p
    summaryDF.loc[1,'perc_overlap']=perc_overlap
    
    df_file='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/QC/dup analysis dfs/dup_comparison_%s_%s_%s' %(sample_name1,sample_name2,seqType)
    summaryDF.to_pickle(df_file)
    
    
    return r,p, ax,fig


#--------------------------------------------------------------------------------
'''
the following function get as input 'dataList': a list of tuples, each tuple is build of (dataName,data), and plot an histogram 
of all datasets, comparing their means, calculate t and ks test (or Anova) and print out the p-values




'''
def plotHistComprison(dataList,ax,title,showLegend=True,nBins=20,toAnnotate=True,
                      colorList=['lightskyblue','mediumorchid'],alpha=None,
                      text_kws={'fontsize':'large','fontweight':'bold','color':'red',
                                'horizontalalignment':'left','verticalalignment':'top'},
                      plotType='hist'):
    
    if alpha is None:
        alpha=1
    ks_p_cohort1_cohort2=np.nan
    t_p_cohort1_cohort2=np.nan
    mean2=np.nan
    
    dataNameList=[]
    allData=[]
    weightList=[]
    meanList=[]
    meanText=''
    for data in dataList:
        dataNameList.append(data[0])
        allData.append(data[1])
        weight=(np.ones_like(data[1]).astype(float))/len(data[1])
        weightList.append(weight)
        mean1=round(np.mean(data[1]),4)
        meanList.append(mean1)
        meanText='%s%s mean=%s\n' %(meanText, data[0], mean1)
        std1=round(np.std(data[1]),4)     
    colorList2=colorList[:len(allData)]
    if plotType=='hist':

        plot=ax.hist(allData, bins=nBins, color=colorList2, weights=weightList,
                         label=dataNameList, alpha=0.7)
    else:
        for n in range(len(allData)):
            print 'generating kde plot'
            plot=sns.kdeplot(np.array(allData[n]), shade=False, color=colorList2[n],ax=ax,label=dataNameList[n],lw=3)
    if len(allData)==2:
        p_Anov=np.nan
        ks_s_cohort1_cohort2, ks_p_cohort1_cohort2=stats.ks_2samp(*allData)
        t_s_cohort1_cohort2, t_p_cohort1_cohort2=stats.ttest_ind(*allData)
        if toAnnotate:
            ax.annotate('ttest_p=%s\n%s mean=%s\n%s mean=%s' %(round(t_p_cohort1_cohort2,6),dataNameList[0],meanList[0],dataNameList[1],meanList[1]),
            xy=(0.02, 0.96), xycoords='axes fraction', **text_kws)
#         ax.annotate('KS_p=%s\nttest_p=%s\n%s mean=%s\n%s mean=%s' %(round(ks_p_cohort1_cohort2,6), 
#         round(t_p_cohort1_cohort2,6),dataNameList[0],meanList[0],dataNameList[1],meanList[1]),
#         xy=(0.96, 0.95), xycoords='axes fraction', fontsize=10, horizontalalignment='right', verticalalignment='top', fontweight='bold')
        
    elif  len(allData)>2:
        ks_p_cohort1_cohort2=np.nan
        t_p_cohort1_cohort2=np.nan
        from scipy.stats import f_oneway
        F_anov,p_Anov=f_oneway(*allData)
        if toAnnotate:
            ax.annotate('p_Anov=%s\n%s' %(round(p_Anov,4),meanText),xy=(0.02, 0.96), xycoords='axes fraction', 
            **text_kws)
    else:
        p_Anov=np.nan
        ks_p_cohort1_cohort2=np.nan
        t_p_cohort1_cohort2=np.nan
        if toAnnotate:
            ax.annotate(meanText,xy=(0.02, 0.96), xycoords='axes fraction', **text_kws)
              
    if title is not None:
            ax.set_title(title, fontsize=16,fontweight='bold')
    if showLegend:
        ax.legend(bbox_to_anchor=(0.05, 1.25), loc='upper left', borderaxespad=0.,fontsize=16)   
    else:
        ax.legend([],[])    
    try:    
        filename='%s_distComparison' %'_'.join(dataNameList)  
    except:
        dataNameListStr=[str(x) for x in dataNameList]
        filename='%s_distComparison' %'_'.join(dataNameListStr)           
            
    return ax,ks_p_cohort1_cohort2,t_p_cohort1_cohort2,p_Anov,filename,meanText


#----------------------------------------------------------------------------------------
# from tunneltoamazondb import getengine

#________________________________

'''
the following function takes us input two lists of data in the same length (can include nans), and their names, and generate a scatter 
plot of their correlation, including trend line, r and p of pearson/spearman correlation

'''
def plot_corr(data1,data2,data1name,data2name,ax,title=None,corrType='pearson',
              toAnnotate=True,plotTrendLine=True,scatter_kws={'alpha':0.4,'s':100},
               text_kws={'fontsize':'x-large'},linecolor='black'):

    
    #clean data: remove nans and outliers
    assert len(data1)==len(data2)
    
    newx=[]
    newy=[]
    for n in range(len(data1)):
        if not np.isnan(data1[n]) and not np.isnan(data2[n]):
            newx.append(data1[n])
            newy.append(data2[n])
    print len(newx)
    print len(newy)          
        
    ymean=np.mean(newy)
    nsamples=len(newx)

    ax.scatter(newx,newy,**scatter_kws)
    ax.set_xlabel(data1name,fontsize='x-large')
    ax.set_ylabel(data2name,fontsize='x-large')
    ax.tick_params(labelsize='large')
    if plotTrendLine:
        ax.plot(np.unique(newx), np.poly1d(np.polyfit(newx, newy, 1))(np.unique(newx)),c=linecolor,linewidth=1,label=data1name+'-'+data2name)
    ax.locator_params(axis='x', nbins=6,rotation=90)
    
    if title is not None:
        ax.set_title(title,fontsize='large')

    from scipy.stats import pearsonr,spearmanr
    if corrType=='spearman':
        r,p = spearmanr(newx,newy)
    else:
        r,p = pearsonr(newx,newy)
    if toAnnotate:
        text=ax.annotate("%s_r=%.2f\np=%.4f" %(corrType,r,p),  xy=(0.02, 0.96), xycoords='axes fraction',
            horizontalalignment='left', verticalalignment='top',**text_kws)
#             fontweight='bold')
    else:
        text=None
        
    handles, labels = ax.get_legend_handles_labels()

   
    return ax, nsamples,r,p,text, handles, labels

#---------------------------------------------------
'''
the following function scans all subfolder in the folder dir_name
and deletes all files within them that their names contain the string
string_in_filename_todelete
'''

def delete_specific_files_from_dir(dir_name,string_in_filename_todelete,print_deleted=True):
    import os
    count=0

    prediction_dirs=listdir(dir_name)
    print ('number of dirs=',len(prediction_dirs))

    for n,pred_dir in enumerate(prediction_dirs):
        if n%5==0: print n
        if isdir(dir_name+pred_dir+'/'):
            for f in listdir(dir_name+pred_dir+'/'):
                if string_in_filename_todelete in f:
                    os.remove(dir_name+pred_dir+'/'+f)
                    count+=1
                    if print_deleted:
                        print 'deleted: '+dir_name+pred_dir+'/'+f
                    
    return count

