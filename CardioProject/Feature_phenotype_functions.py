# from os import listdir, mkdir, makedirs
# from os.path import isfile, join, isdir, exists
# import pandas as pd
# import numpy as np
# from scipy import stats
# from scipy.stats import pearsonr,fisher_exact,ttest_ind,mannwhitneyu, sem
# from scipy.spatial.distance import braycurtis, pdist, squareform
# import re
# import math
# import cPickle as pickle
# import random


# import matplotlib.pyplot as plt
# from matplotlib.ticker import FormatStrFormatter
# from matplotlib.backends.backend_pdf import PdfPages
# import seaborn as sns
#
# from skbio.diversity.alpha import shannon, simpson, berger_parker_d
# from skbio.stats.distance import mantel
# from Bio.SeqUtils import GC
#
# from ShaniBA.MyFunctionsShani import *
# from ShaniBA.myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot, draw_correlation_scatter
from ShaniBA.TCR_microbiome_interactions.TCR_microbiome_interactions_functions import *

from sklearn import metrics, preprocessing


#------------------------------------------------------
def editSampleNamesList(l):
    newNames = [sample.split('_')[0] for sample in l]
    newNames = [sample.split('a')[0] for sample in newNames]
    newNames = [sample.split('b')[0] for sample in newNames]
    newNames = [sample.split('.')[0] for sample in newNames]
    return newNames


def editSampleNames(df):
    df.index = editSampleNamesList(df.index)
#     for n,sample in enumerate(df.index):
#             if '_' in sample:
#                 NewName=sample.split('_')[0]
#             else:
#                 NewName=sample
#             if 'b' in NewName:
#                 NewName=NewName.split('b')[0]
#             df.rename(index={sample:NewName},inplace=True)
#             df.rename(columns={sample:NewName},inplace=True)
    
    return df
    
#----------------------------------

def filter_outliers(df, outlierSTD, columnList=None, trim=None):
    
    df2 = pd.DataFrame()
    if columnList is None:
        columnList = df.columns.values
    for column in columnList:
        columnMean = df[column].mean()
        columnSTD = df[column].std()
        upperThreshold = columnMean + outlierSTD * columnSTD
        bottomThreshold = columnMean - outlierSTD * columnSTD
        
        if trim:
            bottomValue = bottomThreshold
            upperValue = upperThreshold
        else:
            bottomValue = np.nan
            upperValue = np.nan
#         print column,columnMean,columnSTD,bottomThreshold,upperThreshold
        df2[column] = np.where(df[column] < bottomThreshold, bottomValue, df[column])
        df2[column] = np.where(df[column] > upperThreshold, upperValue, df[column])
        df2.index = df.index
    return df2
#---------------------------------------------------------------------------------------
'''
This function takes a sample df and 
1. remove non numeric columns/rows
2.remove empty columns/rows
3.edit sample names
4.order columns/rows by method:
'indexAB' = alphabetically
'sum'=rows and columns sums
'means'rows and column menas
'mins'-rows and columns mins
5. if desired (removeSameUser=True): remove second samples from the same user
***note that those samples to remove are defined manually within the function. need to update if necessary***


'''
def process_sample_matrix(df, removeSameUser, orderMethod):
   
#     print df.head()
    print 'original df shape is %s_%s' % (df.shape[0], df.shape[1])
        
    # remove non numeric columns and empty rows
    for column in df:  # remove nonnumeric columns
        if df[column].dtype == 'object':
            df = df.drop(column, axis=1)
#             print '%s column was dropped from x' %column
    df = df.dropna(axis=(0, 1), how='all')  # remove empty rows
#     print df.head()
 
    # edit sample names:
    df = editSampleNames(df)
    
    # order matrix:
    if orderMethod == 'indexAB':
        print 'ordering matrix alphabetically...'
        df = df.sort_index(axis=0)
        if df.shape[1] > 1:
            df = df.sort_index(axis=1)
    elif orderMethod == 'sum':
        print 'ordering matrix by row and column sums...'
        df = df.reindex_axis(df.sum(axis=0).sort_values().index, axis=0)
        if df.shape[1] > 1:
            df = df.reindex_axis(df.sum(axis=1).sort_values().index, axis=1)
        
    elif orderMethod == 'mean':
        print 'ordering matrix by row and column means...'
        df = df.reindex_axis(df.mean(axis=0).sort_values().index, axis=0)
        if df.shape[1] > 1:
            df = df.reindex_axis(df.mean(axis=1).sort_values().index, axis=1)
        
    elif orderMethod == 'min':
        print 'ordering matrix by row and column mins...'
        df = df.reindex_axis(df.min(axis=0).sort_values().index, axis=0)
        if df.shape[1] > 1:
            df = df.reindex_axis(df.min(axis=1).sort_values().index, axis=1)
        
    else:
        print 'matrix was not re-ordered...'
    
#     print df.head()
    
    
    
    if removeSameUser:  # Remove second sample from the same user:
        print 'removing second sample from the same user...'
        removedSampleList = ['BD714', 'BD838']  # UPDATE LIST WITH MORE SAMPLES COMING FROM AN ALREADY EXIST USER
        for removedSample in removedSampleList:
            if removedSample in df.index:
                df = df.drop(removedSample, axis=0)
#                 print 'removed from df rows'
                if df.shape[1] > 1:
                    df = df.drop(removedSample, axis=1)
#                 print 'removed from df columns'
                
                
    print 'new df shape is %s_%s' % (df.shape[0], df.shape[1])
#     print df.head()
    
            
    return df


#-----------------------------------
'''
This function takes feature and phenotype matrices and:
1. make sure both have the same sample names
2. make sure both have the same sample order
3. remove non-numeric columns
4. remove empty rows/columns
5. if desired (removeSameUser=True): remove second samples from the same user
***note that those samples to remove are defined manually within the function. need to update if necessary***


'''
def common_processing_feature_phenotype_matrices(featureDF, phenotypeDF, removeSameUser):
    x = featureDF
    y = phenotypeDF

    print 'original x array shape is %s_%s' % (x.shape[0], x.shape[1])
    print 'original y array shape is %s_%s' % (y.shape[0], y.shape[1])

    
    # remove non numeric columns and empty rows
    for column in x:  # remove nonnumeric columns
        if x[column].dtype == 'object':
            x = x.drop(column, axis=1)
#             print '%s column was dropped from x' %column
    x = x.dropna(axis=(0, 1), how='all')  # remove empty rows
    
    for column in y:  # remove nonnumeric columns
        if y[column].dtype == 'object':
            y = y.drop(column, axis=1)
#             print '%s column was dropped from y' %column
    y = y.dropna(axis=(0, 1), how='all')  # remove empty rows
    
    
    # make sure samples in the two matrices are identical:
    for sample in x.index:
        if sample not in y.index:
            x = x.drop(sample, axis=0)
            x = x.drop(sample, axis=1)
    for sample in y.index:
        if sample not in x.index:
            y = y.drop(sample, axis=0)
            if y.shape[1] != 1:
                y = y.drop(sample, axis=1)

    
    # make sure both matrices have the same order:
    x = x.sort_index(axis=0)
    x = x.sort_index(axis=1)
    y = y.sort_index(axis=0)
    y = y.sort_index(axis=1)
    
    
    
    if removeSameUser:  # Remove second sample from the same user:
        print 'removing second sample from the same user...'
        removedSampleList = ['BD714', 'BD838']  # UPDATE LIST WITH MORE SAMPLES COMING FROM AN ALREADY EXIST USER
        for removedSample in removedSampleList:
            if removedSample in x.index:
                x = x.drop(removedSample, axis=0)
#                 print 'removed from x rows'
                x = x.drop(removedSample, axis=1)
#                 print 'removed from x columns'
                y = y.drop(removedSample, axis=0)
#                 print 'removed from y rows'
                if len(y.columns.values) > 1:
                    print len(y.columns.values)
                    y = y.drop(removedSample, axis=1)
                    print 'removed from y columns'
#                 print 'sample %s was removed' %removedSample
                
    print 'new x array shape is %s_%s' % (x.shape[0], x.shape[1])
    print 'new y array shape is %s_%s' % (y.shape[0], y.shape[1])
#     print x.head()
#     print y.head()
            
    return x, y

#--------------------------------------------------------------------------------------------------
'''
this function was copied from notebook: Feature-Phenotype interactions in the PNP cohort

'''

file6 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/fullXgroupbyBD_only434'
fullXgroupbyBD_only434 = pd.read_pickle(file6)


def compare_feature_distance_for_binary_phenotypes(feature_dist_file, phenotype, nPerm, removeSameUser, filteringList):
    
    from skbio.stats.distance import permanova
    from skbio.stats.distance import DistanceMatrix
    featureName = feature_dist_file.split('_')[-2] + '_' + feature_dist_file.split('_')[-1]
    print featureName

    # (1)load and process distance matrix files:
    print 'loading and processing distance matrix files...'
    feature_dist_mat = pd.read_pickle(feature_dist_file)
    

    x = feature_dist_mat
    y = pd.DataFrame(fullXgroupbyBD_only434[phenotype])
    
    x, y = common_processing_feature_phenotype_matrices(x, y, True)
    
    if filteringList is not None:
        for sample in x.index:
            if sample not in filteringList:
                x = x.drop(sample, axis=0)
                x = x.drop(sample, axis=1)
        for sample in y.index:
            if sample not in filteringList:
                y = y.drop(sample, axis=0)
                y = y.drop(sample, axis=1)
        
        print ' filtered matrices to include only samples in filteringList:'
        print 'new x array shape is %s_%s' % (x.shape[0], x.shape[1])
        print 'new y array shape is %s_%s' % (y.shape[0], y.shape[1])

    # (2)
    
    DM = DistanceMatrix(x, ids=x.index)
    results = permanova(DM, y, column=phenotype, permutations=nPerm)
#     print results
    p = results['p-value']
    s = results['test statistic']
#     print p
    
    permAnovaResDF = pd.DataFrame()
    permAnovaResDF.loc[0, 'featureName'] = featureName
    permAnovaResDF.loc[0, 'phenotype'] = phenotype
    permAnovaResDF.loc[0, 'nPerm'] = nPerm
    permAnovaResDF.loc[0, 'removeSameUser'] = removeSameUser
    permAnovaResDF.loc[0, 's'] = s
    permAnovaResDF.loc[0, 'p'] = p

    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/permAnovaResults/%s_%s_%s_%s' % (featureName, phenotype, nPerm, removeSameUser)
    permAnovaResDF.to_pickle(file1)
    
    print p
    print 'done'
    
#---------------------------------------------------------------------------------------------------------------------------------
'''
the following function can be used to plot interesting results from the preceding functions. an example to use this function:

#################
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/PNPfemales','rb') as fp:
    PNPfemales=pickle.load(fp)
filteringList=PNPfemales


for n in sigResults.index:
    print n
    featureName=sigResults.loc[n,'featureName']
    phenotype=sigResults.loc[n,'phenotype']
    if phenotype!='Gender':
        print featureName
        print phenotype
        p=sigResults.loc[n,'p']
        nSamplesForVariable=sigResults.loc[n,'nSamples for variable']

        feature_dist_file='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/Features/distMat_PNP434_%s' %featureName
        phenotype_dist_file='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/BinaryPhenotypes/distMat_PNP434_%s_euclidean' %phenotype

        fig2=plot_distances_for_each_group(feature_dist_file,phenotype_dist_file,filteringList)
        fig2.suptitle('Distance comparison between sample pairs with identical and opposite phenotype\nFemales only! Phenotype=%s,p-value=%s,samples positive for phenotype=%s' %(phenotype,p,nSamplesForVariable))     

        figfile='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/permANOVAplotsFemales/%s_%s_permanovaPlot' %(featureName,phenotype)
        fig2.savefig(figfile,dpi=200)
##########

'''
    
def plot_distances_for_each_group(feature_dist_file, phenotype_dist_file, filteringList):
    featureName = feature_dist_file.split('_')[-2] + '_' + feature_dist_file.split('_')[-1]
    phenotypeName = phenotype_dist_file.split('_')[-2] + '_' + feature_dist_file.split('_')[-1]
    
    # (1)load and process distance matrix files:
    print 'loading and processing distance matrix files...'
    x = pd.read_pickle(feature_dist_file)
    y = pd.read_pickle(phenotype_dist_file)

    x, y = common_processing_feature_phenotype_matrices(x, y, True)
    
    if filteringList is not None:
        print 'filtering samples using filteringList...'
        for sample in x.index:
            if sample not in filteringList:
                x = x.drop(sample, axis=0)
                x = x.drop(sample, axis=1)
        for sample in y.index:
            if sample not in filteringList:
                y = y.drop(sample, axis=0)
    
    # #(2) comparing distances between samples positive for the phenotype and samples negative for the phenotype:

    positive = x[y == 1]
    print 'positive mean distance=%s' % positive.mean().mean()

    negative = x[y == 0]
    print 'negative mean distance=%s' % negative.mean().mean()

    negative = negative.values.flatten().tolist()
    positive = positive.values.flatten().tolist()

    newNegative = [i for i in negative if not np.isnan(i)]
    newPositive = [i for i in positive if not np.isnan(i)]
    
    perc5Negative = np.percentile(newNegative, 5)
    perc95Negative = np.percentile(newNegative, 95)
    perc5Positive = np.percentile(newPositive, 5)
    perc95Positive = np.percentile(newPositive, 95)
    minPlot = np.min([perc5Negative, perc5Positive]) * 0.9
    maxPlot = np.max([perc95Negative, perc95Positive]) * 1.1
    

    ks_s, ks_p = stats.ks_2samp(newNegative, newPositive)
    print 'ks_s=%s ks_p=%s' % (ks_s, ks_p)
    t_s, t_p = stats.ttest_ind(newNegative, newPositive)
    print 't_s=%s t_p=%s' % (t_s, t_p)
    
    # #(3)plotting:
    
    fig2 = plt.figure(figsize=(10, 10))


    plt.boxplot((newNegative, newPositive), labels=('Same category', 'Different category'), whis=[5, 95])
    plt.ylabel(featureName + ' distance', fontsize=12)
    plt.ylim(minPlot, maxPlot)
    plt.legend()

    return fig2

#----------------------------------------------------------------------------------------------
'''
general functions to calculate diversity correlations between a sharing sequence matrix and a microbiome matrix
copied from notebook: Feature-Phenotype interactions in the PNP cohort-2

to use this function, use the following code:
#load TCR file ***NOTE THAT IT IS THE RA file!***:
print 'loading TCR file'
file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PublicSeqAnalysis/sharingMatrix_moreThan1_434Samples_RA'
TCRfile=pd.read_pickle(file1)

#load MB file ***NOTE THAT IT IS THE RA file!***:

print 'loading MB file'
file2='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/SampleSpeciesDFgroupedByBD_only434'
microbiomeFile=pd.read_pickle(file2)

TCRfile_binary,TCRfile,microbiomeFile_binary,microbiomeFile=calc_corr_between_TCR_and_microbiome_preprocessing(TCRfile,microbiomeFile)

file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PublicSeqAnalysis/organized_TCR_and_mb_files_for434samples/TCRfile_binary'
TCRfile_binary.to_pickle(file1)

file2='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PublicSeqAnalysis/organized_TCR_and_mb_files_for434samples/TCRfile'
TCRfile.to_pickle(file2)

file3='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PublicSeqAnalysis/organized_TCR_and_mb_files_for434samples/microbiomeFile_binary'
microbiomeFile_binary.to_pickle(file3)

file4='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PublicSeqAnalysis/organized_TCR_and_mb_files_for434samples/microbiomeFile'
microbiomeFile.to_pickle(file4)

TCRcutoff=50
mbCutoff=50


calc_corr_between_TCR_and_microbiome(TCRfile_binary,TCRfile,microbiomeFile_binary,microbiomeFile,TCRcutoff,mbCutoff)
'''
def reject_outliers(data, m):
    outlier_ind = abs(data - np.mean(data)) < m * np.std(data)
    return outlier_ind


def plot_corr_diversity(divDF_seqs, divDF_mb, measure, ax, stdToReject):

    x = divDF_seqs[measure]
    y = divDF_mb[measure]
    
    print 'checking TCR and MB df order...'
    print x.head()
    print y.head()
    
    # clean data: remove nans and outliers:
    nx = np.isnan(x)
    ny = np.isnan(y)
    n = nx + ny
    newx = x[~n]
    newy = y[~n]
    
    
    if stdToReject is not None:
        nx_outliers = reject_outliers(newx, m=stdToReject)
        ny_outliers = reject_outliers(newy, m=stdToReject)
        n_outliers = nx_outliers + ny_outliers
        finalx = newx[~n_outliers]
        finaly = newy[~n_outliers]
    else:
        finalx = newx
        finaly = newy
        
        
    ymean = np.mean(finaly)
    nsamples = len(finalx)

    ax.scatter(finalx, finaly, alpha=0.4)
    ax.set_xlabel('TCR sequences')
    ax.set_ylabel('Microbiome Species')
    ax.plot(np.unique(finalx), np.poly1d(np.polyfit(finalx, finaly, 1))(np.unique(finalx)), c='blue', linewidth=1)
    ax.set_title('%s' % measure, fontsize=16)

    from scipy.stats import pearsonr
    r, p = pearsonr(finalx, finaly)

    ax.annotate("r=%.4f p=%.6f,n=%s" % (r, p, nsamples), xy=(0.02, 0.96), xycoords='axes fraction', fontsize=11,
        horizontalalignment='left', verticalalignment='top')

    # if minPhenotypeValue is not None:
    #     plt.ylim(minPhenotypeValue,np.max(y)*1.1)
    #         plt.margins(0.2)

    #     file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/DistMat_correlation_plots/'
    #     fig1.savefig(file1,dpi=200)
    
    return nsamples, r, p



def calc_corr_between_TCR_and_microbiome_preprocessing(TCRfile, microbiomeFile):   
    # (1)process sample names:
    
    
    for dfFile in [TCRfile, microbiomeFile]:
        print 'processing sample names...'
        for n, sample in enumerate(dfFile.index):
            print n
            if '_' in sample:
                NewName = sample.split('_')[0]
            else:
                NewName = sample
            if 'b' in NewName:
                NewName = NewName.split('b')[0]
            dfFile.rename(index={sample:NewName}, inplace=True)
            dfFile.rename(columns={sample:NewName}, inplace=True)
    
    # (2)transform RA files to binary:
    
    print 'now converting TCR counts to binary indications...'          
    TCRfile_binary = pd.DataFrame()
    for column in TCRfile.columns.values:
        TCRfile_binary[column] = np.where(TCRfile[column] > 0, 1, 0)
        TCRfile_binary.index = TCRfile.index

    print 'now converting mb counts to binary indications...'          
    microbiomeFile_binary = pd.DataFrame()
    for column in microbiomeFile.columns.values:
        microbiomeFile_binary[column] = np.where(microbiomeFile[column] > 0, 1, 0)
        microbiomeFile_binary.index = microbiomeFile.index
        
    return TCRfile_binary, TCRfile, microbiomeFile_binary, microbiomeFile


def calc_corr_between_TCR_and_microbiome(TCRfile_binary, TCRfile, microbiomeFile_binary, microbiomeFile, TCRcutoff, mbCutoff):
            
    # (3)truncate files to include only sequences/species shared by more than a cutoff number of samples:
    if 'FD' in TCRfile_binary.columns.values:
            TCRfile_binary = TCRfile_binary.drop('FD', axis=1)
            print 'FD column was dropped from TCR file'
    if 'FD' in microbiomeFile_binary.columns.values:
            microbiomeFile_binary = microbiomeFile_binary.drop('FD', axis=1)
            print 'FD column was dropped from MB file'
    
    if TCRcutoff is not None:
        print 'truncating TCR File to include only sequences shared by more than %s' % TCRcutoff
        nSeqsBefore = len(TCRfile_binary.columns.values)
        print 'number of sequences before truncation=%s' % nSeqsBefore
        columnList = []
        for n, column in enumerate(TCRfile_binary.columns.values):
            if n % 1000 == 0:
                print n
            nSamples = TCRfile_binary[column].sum()
            if nSamples > TCRcutoff:
                columnList.append(column)
        
        print 'number of sequences after truncation=%s' % len(columnList)

        TCRfile_binary_truncated = TCRfile_binary[columnList]
        TCRfile_RA_truncated = TCRfile[columnList]
    else:
        print 'TCR files doesnt need truncation...'
        TCRfile_binary_truncated = TCRfile_binary
        TCRfile_RA_truncated = TCRfile
        
        
    if mbCutoff is not None:
        print 'truncating MB File to include only sequences shared by more than %s' % mbCutoff
        nSeqsBefore = len(microbiomeFile_binary.columns.values)
        print 'number of sequences before truncation=%s' % nSeqsBefore
        columnList = []
        for n, column in enumerate(microbiomeFile_binary.columns.values):
            if n % 1000 == 0:
                print n
            nSamples = microbiomeFile_binary[column].sum()
            if nSamples > mbCutoff:
                columnList.append(column)
        if 'FD' in columnList:
            columnList.remove('FD')
        print 'number of sequences after truncation=%s' % len(columnList)

        microbiomeFile_binary_truncated = microbiomeFile_binary[columnList]
        microbiomeFile_RA_truncated = microbiomeFile[columnList]
    else:
        print 'MB files doesnt need truncation...'
        microbiomeFile_binary_truncated = microbiomeFile_binary
        microbiomeFile_RA_truncated = microbiomeFile
        
        
    # (4)calculate diversity measures for TCRs:
    
    print 'calculating diversity measures for TCR...'
    
    dfList = [TCRfile_binary_truncated, TCRfile_RA_truncated]
    dfName = 'TCR_moreThan%s' % TCRcutoff
    isRAList = [False, True]

    divDF_seqs = pd.DataFrame(index=TCRfile_binary_truncated.index)

    for n, df in enumerate(dfList):
        print 

#         if 'FD' in df.columns.values:
#             df=df.drop('FD',axis=1)

        isRA = isRAList[n]
        if isRA:
            RA = 'RA'
            df = df.round(5) * 100000
            df = df.astype(int, errors='ignore')      
        else:
            RA = 'binary'

        for sample in df.index:
    #         print sample
            divDF_seqs.loc[sample, 'shannon_%s' % RA] = shannon(df.loc[sample, :], base=2)
            divDF_seqs.loc[sample, 'simpson_%s' % RA] = simpson(df.loc[sample, :])
            divDF_seqs.loc[sample, 'berger_parker_d_%s' % RA] = berger_parker_d(df.loc[sample, :])
            if isRA:
                divDF_seqs.loc[sample, 'maxFreq_%s' % RA] = np.max(df.loc[sample, :])
                divDF_seqs.loc[sample, 'meanFreq_%s' % RA] = np.mean(df.loc[sample, :])
            else:
                divDF_seqs.loc[sample, 'nUnique'] = np.sum(df.loc[sample, :])
                
    print 'calculating diversity measures for MB...'
    
    dfList = [microbiomeFile_binary_truncated, microbiomeFile_RA_truncated]
    dfName = 'microbiome_moreThan%s' % TCRcutoff
    isRAList = [False, True]

    divDF_mb = pd.DataFrame(index=microbiomeFile_binary_truncated.index)

    for n, df in enumerate(dfList):
        print 'df number=%s' % n

#         if 'FD' in df.columns.values:
#             df=df.drop('FD',axis=1)

        isRA = isRAList[n]
        if isRA:
            RA = 'RA'
            df = df.round(5) * 100000
            df = df.astype(int, errors='ignore')      
        else:
            RA = 'binary'

        for sample in df.index:
    #         print sample
            divDF_mb.loc[sample, 'shannon_%s' % RA] = shannon(df.loc[sample, :], base=2)
            divDF_mb.loc[sample, 'simpson_%s' % RA] = simpson(df.loc[sample, :])
            divDF_mb.loc[sample, 'berger_parker_d_%s' % RA] = berger_parker_d(df.loc[sample, :])
            if isRA:
                divDF_mb.loc[sample, 'maxFreq_%s' % RA] = np.max(df.loc[sample, :])
                divDF_mb.loc[sample, 'meanFreq_%s' % RA] = np.mean(df.loc[sample, :])
            else:
                divDF_mb.loc[sample, 'nUnique'] = np.sum(df.loc[sample, :])
                
    # plotting correlation scatters between TCR and microbiome, with different stdToReject:
    
    print'plotting correlation scatter for stdToRejectList=...'
    stdToRejectList = [None, 0.25, 0.5]
    
    for stdToReject in stdToRejectList:
        print stdToReject

        fig1 = plt.figure(figsize=(12, 12))
        fig1.suptitle('Correlations between TCR and microbiome diversities\nRemoved outliers=%s' % stdToReject,
                     fontsize=18)
        sumDF = pd.DataFrame()



        for n, measure in enumerate(divDF_seqs.columns.values):
            print n, measure
            ax = fig1.add_subplot(3, 3, n + 1)
            nsamples, r, p = plot_corr_diversity(measure, ax, stdToReject)
            sumDF.loc[n, 'TCRcutoff'] = TCRcutoff
            sumDF.loc[n, 'mbCutoff'] = mbCutoff
            sumDF.loc[n, 'stdToReject'] = stdToReject
            sumDF.loc[n, 'measure'] = measure
            sumDF.loc[n, 'nsamples'] = nsamples
            sumDF.loc[n, 'r'] = r
            sumDF.loc[n, 'p'] = p
            
        file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/Diversity analysis/sumDFs/DF_sharedSeqMoreThan%s_mbMoreThan%s_rejectMoreThan%s' % (TCRcutoff, mbCutoff, stdToReject)        
        sumDF.to_pickle(file1)  
            
            

        fig1.subplots_adjust(left=0.09, right=0.98, top=0.9, bottom=0.02, wspace=0.25, hspace=0.30)
      
        stdToRejectNameList = str(stdToReject).split('.')
        if len(stdToRejectNameList) == 1:
            stdToRejectName = stdToRejectNameList[0]
        else:
            stdToRejectName = stdToRejectNameList[0] + stdToRejectNameList[1]
      
        file2 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/Diversity analysis/sharedSeqMoreThan%s_mbMoreThan%s_rejectMoreThan%s' % (TCRcutoff, mbCutoff, stdToRejectName)
        fig1.savefig(file2, dpi=200) 


#--------------------------------------------------------------------------------------------------------

'''
this function generates distance matrix for feature/phenotype matrices
copied from notebook: Distance Matrices preperation and mantel tests

matrix_name - free text to describe the feature/phenotype matrix (String)
matrix_df - matrix df. in the former version of the function I used the matrix_file as input and the function downloaded the matrix from the file
distance_measure = 'euclidean' for most cases, 'braycurtis' for 'ecological' matrices such as microbiome
species relative abundance or public TCR relative abundance
do_binary - True if necessary to convert the original matrix into a binary form

an example for using the function:
matrix_name='MicrobiomeSpeciesRA'
matrix_file='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/SampleSpeciesDFgroupedByBD_only434'
matrix_df=pd.read_pickle(matrix_file)
distance_measure='braycurtis'
do_binary = False

compute_distance_matrix_general(matrix_name,matrix_df, distance_measure, do_binary = False)
'''


def compute_distance_matrix_general(matrix_name, matrix_df, distance_measure, do_binary):
    
    from scipy.spatial.distance import braycurtis, pdist
    
#     #load and process matrix file: ### I updated the function and now it gets as input the df itself and doesn't load it###
#     samps_by_features=pd.read_pickle(matrix_file)
    samps_by_features = matrix_df
    
    
    print 'original df shape is %s_%s' % (samps_by_features.shape[0], samps_by_features.shape[1])
    for column in samps_by_features:  # remove nonnumeric columns
        if samps_by_features[column].dtype == 'object':
            samps_by_features = samps_by_features.drop(column, axis=1)
            print '%s column was dropped from df' % column
    samps_by_features = samps_by_features.dropna(axis=(0, 1), how='all')
    print 'df shape after processing is %s_%s' % (samps_by_features.shape[0], samps_by_features.shape[1])
    
    if do_binary:
        isBinary = 'binary'
    else:
        isBinary = ''
    
    
    # calculate the distance matrix:
    print 'calculating distance matrix...'
    A = samps_by_features.T.copy()
    temp_res = pd.DataFrame(1, index=A.columns, columns=A.columns, dtype=np.float64)
    A = A.values
    if do_binary:
        A[A > 0] = 1
    for i in range(A.shape[1]):
        a = A[:, i]
        vals = []
        if i % 50 == 0:
            print i
        for j in range(A.shape[1]):
            b = A[:, j]
            val = None
            val = pdist(np.stack((a, b), axis=1).T, distance_measure)
            vals.append(val)
        temp_res.iloc[i, :] = vals
        
    # edit the distance matrix:
    # (1) edit sample names to match those in the phenotype files
    # (2) remove rows and columns with nan values:
    print 'edit distance matrix...'
    for n, sample in enumerate(temp_res.index):
        if '_' in sample:
            NewName = sample.split('_')[0]
        else:
            NewName = sample
        if 'b' in NewName:
            NewName = NewName.split('b')[0]
        temp_res.rename(index={sample:NewName}, inplace=True)
        temp_res.rename(columns={sample:NewName}, inplace=True)
        
#     print 'original distMat shape is %s_%s' %(temp_res.shape[0],temp_res.shape[1])
#     bool_row_idx = temp_res.isnull().any(axis=0)
#     row_idx_list=temp_res[bool_row_idx].index
#     bool_col_idx = temp_res.isnull().any(axis=1)
#     col_idx_list=temp_res[bool_col_idx].index
#     for row_idx in row_idx_list:
#         for col_idx in col_idx_list:
#             print '%s row was dropped from distMat' %row_idx
#             temp_res=temp_res.drop(row_idx,axis=0)
#             print '%s column was dropped from distMat' %col_idx
#             temp_res=temp_res.drop(col_idx,axis=1)
    print 'distMat shape after processing is %s_%s' % (temp_res.shape[0], temp_res.shape[1])
    
    # save resulting distance matrix to file:
    print 'saving to file...'
    file2 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/distMat_PNP434_%s%s_%s' % (matrix_name,
                                                          isBinary, distance_measure)   
    temp_res.to_pickle(file2)
        
    return temp_res

#-----------------------------------------------------------------------------
'''
define function to run matel test for different distance matrix combinations
copied from notebook: Distance Matrices preperation and mantel tests
this function had its own .py file in eclipse

an example for using this function"
method='spearman'
n_permut=9999
alternative='two-sided'

fileNameP='distMat_PNP434_Age_euclidean'
fileNameF='distMat_PNP434_sharingMatrixMoreThan1RA_braycurtis'

feature_dist_file='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/Features/%s' %fileNameF
feature_name=fileNameF.split('_')[-2]+'_'+fileNameF.split('_')[-1]

phenotype_dist_file='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/newPhenotypes/%s' %fileNameP
phenotype_name=fileNameP.split('_')[-2]+'_'+fileNameP.split('_')[-1]
          
mantelTest_feature_phenotype(feature_dist_file, phenotype_dist_file,feature_name,phenotype_name, method,n_permut,alternative)

'''


def mantelTest_feature_phenotype(distMat1, distMat2, distMat1name, distMat2name, MantelFolder, method, n_permut, alternative, removeSameUser):
    
    from skbio.stats.distance import mantel
    
    
    # (1)load and process distance matrix files:
    x = distMat1
    y = distMat2
    
    print 'original x array shape is %s_%s' % (x.shape[0], x.shape[1])
    print 'original y array shape is %s_%s' % (y.shape[0], y.shape[1])
    print 'loading and processing distance matrix files...'
    

    print 'verifying same samples in compared distance matrices:'

    for sample in x.index:
        if sample not in y.index:
            x = x.drop(sample, axis=0)
            x = x.drop(sample, axis=1)
    for sample in y.index:
        if sample not in x.index:
            y = y.drop(sample, axis=0)
            y = y.drop(sample, axis=1)
            
    for column in x:  # remove nonnumeric columns
        if x[column].dtype == 'object':
            x = x.drop(column, axis=1)
            print '%s column was dropped from x' % column
    x = x.dropna(axis=(0, 1), how='all')  # remove empty rows
    
    for column in y:  # remove nonnumeric columns
        if y[column].dtype == 'object':
            y = y.drop(column, axis=1)
            print '%s column was dropped from y' % column
    y = y.dropna(axis=(0, 1), how='all')  # remove empty rows

    
    # make sure both matrices have the same order:
    x = x.sort_index(axis=0)
    x = x.sort_index(axis=1)
    y = y.sort_index(axis=0)
    y = y.sort_index(axis=1)

    
    if removeSameUser:  # Remove second sample from the same user:
        from ShaniBA.TCR_microbiome_interactions.TCR_microbiome_interactions_functions import filterSamePerson
        print 'filtering same person samples'
        x = filterSamePerson(df=x, axisList=[0, 1])
        y = filterSamePerson(df=y, axisList=[0, 1])
        
                
    print 'new x array shape is %s_%s' % (x.shape[0], x.shape[1])
    print 'new y array shape is %s_%s' % (y.shape[0], y.shape[1])
    print x.iloc[:5, :5]
    print y.iloc[:5, :5]
    print 'make sure matrices dimensions and heads are identical!'
             
    
    # (2) run mantel:
    print 'running mantel test...'
    corr_coeff, p_value, N = mantel(x, y, method=method, permutations=n_permut, alternative=alternative, strict=True, lookup=None)

    mantelDF = pd.DataFrame()
    mantelDF.loc[0, 'distMat1name'] = distMat1name
    mantelDF.loc[0, 'distMat2name'] = distMat2name
    mantelDF.loc[0, 'method'] = method
    mantelDF.loc[0, 'n_permut'] = n_permut
    mantelDF.loc[0, 'alternative'] = alternative
    mantelDF.loc[0, 'corr_coeff'] = corr_coeff
    mantelDF.loc[0, 'p_value'] = p_value
    mantelDF.loc[0, 'N'] = N

    print corr_coeff, p_value
    
    print 'saving mantel results'
    
#         mantelDF.to_pickle(MantelFile)
    
#         #generate correlation plot if desired:
#         if p_value<0.05 or corr_coeff>0.15:
#             fig1=plot_corr_for_distMats(x,y,feature_name,phenotype_name)
#             file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/DistMat_correlation_plots/%s_%s' %(feature_name,phenotype_name)
#             fig1.savefig(file1,dpi=200)
        
    
    
    print 'DONE!'
           
    return mantelDF

#---------------------------------------------------------------------

'''
copied from notebook: Distance Matrices preperation and mantel tests
The function: plot correlation scatter from two distance matrices:

x=feature matrix
y=phenotype matrix

x and y should be preprocessed to be at the same shape!

'''


def plot_corr_for_distMats(featureDistMat, PhenotypeDistMat, FeatureName, PhenotypeName, minPhenotypeValue=None):
    
    
    
    # edit matrices:
    print 'processing distance matrix tables...'
    for sample in featureDistMat.index:
            if sample not in PhenotypeDistMat.index:
                featureDistMat = featureDistMat.drop(sample, axis=0)
                featureDistMat = featureDistMat.drop(sample, axis=1)
    for sample in PhenotypeDistMat.index:
        if sample not in featureDistMat.index:
            PhenotypeDistMat = PhenotypeDistMat.drop(sample, axis=0)
            PhenotypeDistMat = PhenotypeDistMat.drop(sample, axis=1)
            
    for column in featureDistMat:  # remove nonnumeric columns
        if featureDistMat[column].dtype == 'object':
            featureDistMat = featureDistMat.drop(column, axis=1)
            print '%s column was dropped from x' % column
    featureDistMat = featureDistMat.dropna(axis=(0, 1), how='all')  # remove empty rows
    
    for column in PhenotypeDistMat:  # remove nonnumeric columns
        if PhenotypeDistMat[column].dtype == 'object':
            PhenotypeDistMat = PhenotypeDistMat.drop(column, axis=1)
            print '%s column was dropped from PhenotypeDistMat' % column
    PhenotypeDistMat = PhenotypeDistMat.dropna(axis=(0, 1), how='all')  # remove empty rows

    
    # make sure both matrices have the same order:
    featureDistMat = featureDistMat.sort_index(axis=0)
    featureDistMat = featureDistMat.sort_index(axis=1)
    PhenotypeDistMat = PhenotypeDistMat.sort_index(axis=0)
    PhenotypeDistMat = PhenotypeDistMat.sort_index(axis=1)
    
    # generate lists of feature and phenotype distances for each sample combination:
    indList = []
    xValueList = []
    yValueList = []

    for i, row in enumerate(featureDistMat.index):
        for j, column in enumerate(featureDistMat.columns.values):
            if row != column:
#                 print row,column
                ind = row + '_' + column
                indList.append(ind)
                xValue = featureDistMat.loc[row, column]
                xValueList.append(xValue)
                yValue = PhenotypeDistMat.loc[row, column]
                yValueList.append(yValue)

    x = xValueList
    y = yValueList

    fig1 = plt.figure(figsize=(8, 8))

    ymean = np.mean(y)                                                                                                                      

    plt.scatter(x, y, alpha=0.1)
    plt.xlabel(FeatureName, fontsize=14)
    plt.ylabel(PhenotypeName, fontsize=14)
    plt.plot(np.unique(x), np.poly1d(np.polyfit(x, y, 1))(np.unique(x)), c='blue', linewidth=1)
    plt.title('Distance matrix correlation', fontsize=16)
    
    from scipy.stats import pearsonr
    r, p = MyPearsonr(x, y)
    
    nSamples = len(featureDistMat)


    plt.annotate("r=%.4f p=%.6f, nSamples=%s" % (r, p, nSamples), xy=(0.02, 0.96), xycoords='axes fraction', fontsize=14,
        horizontalalignment='left', verticalalignment='top', fontweight='bold')
    
    if minPhenotypeValue is not None:
        plt.ylim(minPhenotypeValue, np.max(y) * 1.1)
#         plt.margins(0.2)
    
#     file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/DistMat_correlation_plots/'
#     fig1.savefig(file1,dpi=200)


    return fig1
#-------------------------------------------------------------------------
'''
this function takes a distance matrix and order all pairs in ascending order of distnaces.
an example for using this function:

distMatNameList=['50_euclidean','50RA_braycurtis','1Binary_euclidean','1RA_braycurtis',
                '10_euclidean','10RA_braycurtis']

for n, distMatName in enumerate(distMatNameList):
    distMatFile='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/Features/distMat_PNP434_sharingMatrixMoreThan%s' %distMatName  
    distMat=pd.read_pickle(distMatFile)
    
    print n, distMatName
    order_sample_pairs_by_distance(distMat,distMatName)
'''



def order_sample_pairs_by_distance(distMat, distMatName):
    
    distanceVector = pd.DataFrame()
    nSamples = distMat.shape[0]
    count = 0
    
    for i in range(nSamples):
        if i % 10 == 0:
            print i
        for j in range(i + 1, nSamples):
    #         print i,j
            sample1 = distMat.index[i]
            sample2 = distMat.columns.values[j]
            distance = distMat.iloc[i, j]
    #         print sample1,sample2,distance
            distanceVector.loc[count, 'sample1'] = sample1
            distanceVector.loc[count, 'sample2'] = sample2
            distanceVector.loc[count, 'distance'] = distance
            count = count + 1
    distanceVector = distanceVector.sort_values(by='distance')
    
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Distance Matrices/samplePairOrderByDistance/distOrder_%s' % distMatName
    distanceVector.to_pickle(file1)
    
    return distanceVector

#---------------------------------------------------------------------------------

def add_corrected_pValues(resultsDF, pValueColumn, nTests, FDR):
    
    from statsmodels.sandbox.stats.multicomp import fdrcorrection0
    
    # Bonferroni correction:
    bon_correct_p = 0.05 / nTests
    resultsDF['%s_Sig_Bonferroni' % pValueColumn] = np.where(resultsDF[pValueColumn] < bon_correct_p, 1, 0)
    
    
    # FDR correction:
    
    pList = resultsDF[pValueColumn]
    if FDR is None:
        FDR = round(1. / len(pList), 6)    
    b, p = fdrcorrection0(pList, alpha=FDR)
    sigList = [1 if x else 0 for x in b]
    resultsDF['%s_sigFDR=%s' % (pValueColumn, FDR)] = sigList
    resultsDF['%s_corrPval_DR%s' % (pValueColumn, FDR)] = p

    return resultsDF

#---------------------------------------------------------------------------------
'''
USE THE FOLLOWING FUNCTION AS GROUPBY-AGG FUNCTION
this function takes data, removes outliers and calculate data mean accordng to the following rules:
1. if std==0, or data length is less than nMinSamples, regular mean is calculated.
2. if the resulting clean Mean is Nan, regular mean is calculated.
3. in other cases, only datapoints lying between mean-nSTD*std and mean+nSTD*std are being used. 

usage example:
nSTD=10
nMinSamples=3
df2=Ndf1.groupby('ConnectionID').agg(lambda x: noOutlierMean(x,nSTD,nMinSamples))
'''

def noOutlierMean(x, nSTD, nMinSamples):
    mean = np.mean(x)
    
    if (len(x) > nMinSamples) & (np.std(x) != 0):
        cleanMean = np.mean([i for i in x if (i > np.mean(x) - nSTD * np.std(x)) & (i < np.mean(x) + nSTD * np.std(x))])
        if cleanMean.dtype == float and np.isnan(cleanMean):
            cleanMean = mean
    else:
        cleanMean = mean
    
    
                          
    return cleanMean
        
#-----------------------------------------------------------------------------------------
'''
the following function takes a df with phenotype columns and sample rows and filter out outlier samples for each phenotype (column).
the filtering is done by assigning 'np.nan' to the cell instead of a value

inputs:
df=phenotype dataframe
nSTD= number of std to filter out samples
nMinSamples- number of minimal samples required for outlier filtering (for example, it doesnt make sense to filter DISTANCE when there
are 2 samples only...
'''


def filter_phenotypiclly_outlier_samples(df, nSTD, nMinSamples):
    for column in df.columns.values:
        if ('int' in str(df[column].dtype)) or ('float' in str(df[column].dtype)):
            mean = df[column].mean()
            std = df[column].std()
            
            if (std > 0) and (df[column].count() > nMinSamples):
                lowLim = mean - nSTD * std
                highLim = mean + nSTD * std
                df[column] = np.where((df[column] > lowLim) & (df[column] < highLim), df[column], np.nan)
            
    return df



#----------------------------------------------------------------------------------------
''' usage example for the two functions together:

file1='%s/MicrobiomeDataTables/FilteredAndMergedOnBD/MPA_s_AllSeqProjects_SampleListPNP515_filterGenotekTrue_filterMinimalRead9000000_libMethNone_noOutlierMeannSTD5nMinSamples3_fOutlierSampleFalse_fSamePersonFalse' %MyPath

df=pd.read_pickle(file1)
dfName='MbdfOld'

df_distMat_square,df_distMat_square_binary=genDistMatsForAll(df,dfName)
'''



     
def genDistMat(df, metric,generateSquareform=True):

    # generate distance series:
    print 'generating condensed distance matrix using %s' % metric
    distMat_condensed = pd.Series(pdist(df, metric))

    # add sample names:
    sample1List = []
    sample2List = []

    count = 0
    for i in range(len(df.index)):
        for j in range(i + 1, len(df.index)):
            sample1 = df.index[i]
            sample2 = df.index[j]
            sample1List.append(sample1)
            sample2List.append(sample2)
            count = count + 1
    #     print count
    df_condensed_org = pd.DataFrame({'sample1':sample1List, 'sample2':sample2List, 'dist':distMat_condensed})
    print 'top similar pairs:'
    print df_condensed_org.sort_values(by='dist').head(10)
    

    # GENERATE SQUAREFORM DISTANCE MATRIX AND SAVE:
    if generateSquareform:
        print 'generating squareform distance matrix using %s' % metric
        distMat_square = pd.DataFrame(squareform(pdist(df, metric)), columns=df.index, index=df.index)
        #     print TCRdf_RA_distMat_square.iloc[:5,:5]
    else:
        distMat_square=None
    
    return df_condensed_org, distMat_square     
          


def genDistMatsForAll(df, dfName, dataset_folder, condensedFile,):
    TCR_distMatFolder = '%s/distanceMatrices' % dataset_folder
    MB_distMatFolder = '%s/MicrobiomeDataTables/distanceMatrices' % MyPath
    
    print 'generating distance matrix for RA data, df=%s' % dfName
    metric1 = 'braycurtis'
    df_condensed_org, df_distMat_square = genDistMat(df, metric=metric1)
    
    print 'saving RA distMat CONDENSED'
   
    if 'Mbdf' in dfName:
        folder = MB_distMatFolder   
    else:
        folder = TCR_distMatFolder
    
    print 'saving RA distMat'
    RAdistMatFile = '%s/%s_%s_distMat' % (folder, dfName, metric1)
    RAdistMatFileCondensed = '%s/%s_%s_distMat_CONDENSED' % (folder, dfName, metric1)
    df_distMat_square.to_pickle(RAdistMatFile)
    df_condensed_org.to_pickle(RAdistMatFileCondensed)
    
    if 'norm' not in dfName:
        print 'converting data to binary'
        minVal = df.min().min()
        print 'the min value of this df is %s' % minVal
        df2 = df.copy()
        df2 = (df2 > minVal).astype(int)
        print df2.iloc[:3, :3]

        print 'saving binary df...'
        binaryDFfile = '%s/%s_binary' % (folder, dfName)
        df2.to_pickle(binaryDFfile)

        print 'generating binary distance metric, df=%s' % dfName
        metric2 = 'jaccard'

        df_condensed_org_binary, df_distMat_square_binary = genDistMat(df2, metric=metric2)
        print 'saving binary distMat'
        binaryDistMatFile = '%s/%s_binary_%s_distMat' % (folder, dfName, metric2)
        binaryDistMatFileCondensed = '%s/%s_binary_%s_distMat_CONDENSED' % (folder, dfName, metric2)
        df_distMat_square_binary.to_pickle(binaryDistMatFile)
        df_condensed_org_binary.to_pickle(binaryDistMatFileCondensed)
    else:
        df_distMat_square_binary = pd.DataFrame()
        print 'No binary df for normed df!'
        
    print 'DONE!'
    
    return df_distMat_square, df_distMat_square_binary
    
#-----------------------------------------------------------------------------------------------

def permanova_withWrapper(feature_distMat, feature_name, phenotypeDF, phenotypeColumn, nPerm, removeSameUser, permanovaDFfolder,
                          i=None, j=None):
    
    from skbio.stats.distance import permanova
    from skbio.stats.distance import DistanceMatrix
    
    permanovaName = '%s_%s_%s_%s' % (feature_name, phenotypeColumn, nPerm, removeSameUser)
    permanovaFile = '%s/%s' % (permanovaDFfolder, permanovaName)
    existingPermanovas = [f for f in listdir(permanovaDFfolder) if isfile(join(permanovaDFfolder, f))]
    permAnovaResDF = pd.DataFrame()
    
    if permanovaName not in existingPermanovas:


        # (1)load and process distance matrix files:
        print 'load and process distMat and phenotype...'
        x = feature_distMat
        y = pd.DataFrame(phenotypeDF[phenotypeColumn])
        
        # compare indices and order:
        # leave only common samples in each df (X and seq1data)
        
        x = x.sort_index()
        xSamples = [str(i) for i in x.index.tolist() if i in y.index]
        x = x.loc[xSamples, xSamples]
#         x = x.dropna(how='any')
        
        ySamples = [str(i) for i in y.index.tolist() if i in x.index]
        y = y.loc[ySamples, :]
        y = y.sort_index()
        
        
        print 'feature distMat shape is %s_%s' % (x.shape[0], x.shape[1])
        print '5th sample in x is %s and 30th is %s' % (x.index[5], x.index[30])
        print x.iloc[:4, :4]
        print 'phenotype shape is %s_%s' % (y.shape[0], y.shape[1])
        print '5th sample in y is %s and 30th is %s' % (y.index[5], y.index[30])
        print y.iloc[:4, :4]
    
        # (2)
        print 'conducting permanova...'
        
        DM = DistanceMatrix(x, ids=x.index)
        results = permanova(DM, y, column=phenotypeColumn, permutations=nPerm)
    #     print results
        p = results['p-value']
        s = results['test statistic']
        print p
        
        
        permAnovaResDF.loc[0, 'featureName'] = feature_name
        permAnovaResDF.loc[0, 'phenotype'] = phenotypeColumn
        permAnovaResDF.loc[0, 'nPerm'] = nPerm
        permAnovaResDF.loc[0, 'removeSameUser'] = removeSameUser
        permAnovaResDF.loc[0, 's'] = s
        permAnovaResDF.loc[0, 'p'] = p
    
        permAnovaResDF.to_pickle(permanovaFile)
        
        print p
#         print 'done permanova'
    else:
        print 'This permanova test result already exist in folder'
                                           
                                           
    return permAnovaResDF


#--------------------------------------------------------------------------------
'''
new mantel function! June2018

usage example:
see in notebook:  AgeGenderBMIcreatenineSmoking_effect_on_TCRfeatures_cardio+PNP




'''

def mantel_withWrapper(feature_distMat, feature_name, phenotype_distMat, phenotype_name, mantelDFfolder,
                       nPerm, method='spearman', alternative='greater', removeSameUser=True, i=None, j=None):
    
    from skbio.stats.distance import mantel
    from skbio.stats.distance import DistanceMatrix
    
    mantelName = '%s_%s_%s_%s_%s_%s_%s' % (feature_name, phenotype_name, nPerm, method, alternative, removeSameUser, method)
    mantelFile = '%s/%s' % (mantelDFfolder, mantelName)
    existingmantels = [f for f in listdir(mantelDFfolder) if isfile(join(mantelDFfolder, f))]
    mantelResDF = pd.DataFrame()
    
    if mantelName not in existingmantels:


        # (1)load and process distance matrix files:
        print 'load and process distMat and phenotype...'
        x = feature_distMat
        x = x.dropna(axis=(0, 1), thresh=2)  # remove empty rows
        y = phenotype_distMat
        y = y.dropna(axis=(0, 1), thresh=2)  # remove empty rows
        
        # compare indices and order:
        # leave only common samples in each df (X and seq1data)
        
        xSamples = [str(i) for i in x.index.tolist() if i in y.index]
        x = x.loc[xSamples, xSamples]
        x = x.sort_index()
        x = x.sort_index(axis=1)
        
        ySamples = [str(i) for i in y.index.tolist() if i in x.index]
        y = y.loc[ySamples, ySamples]
        y = y.sort_index()
        y = y.sort_index(axis=1)
        
        
        print 'feature distMat shape is %s_%s' % (x.shape[0], x.shape[1])
        print '5th sample in x is %s and 30th is %s' % (x.index[5], x.index[30])
        print x.iloc[:4, :4]
        print 'phenotype distMat shape is %s_%s' % (y.shape[0], y.shape[1])
        print '5th sample in y is %s and 30th is %s' % (y.index[5], y.index[30])
        print y.iloc[:4, :4]
    
        # (2)
        print 'conducting mantel...'
        
        DMx = DistanceMatrix(x, ids=x.index)
        DMy = DistanceMatrix(y, ids=y.index)
        
#         print 'feature distMat shape is now %s_%s' %(DMx.shape[0], DMx.shape[1])
# #         print '100th sample in x is %s and 200th is %s' %(DMx.index[100],DMx.index[200])
#         print DMx.iloc[:4,:4]
#         print 'phenotype distMat shape is %s_%s' %(DMy.shape[0], DMy.shape[1])
# #         print '100th sample in y is %s and 200th is %s' %(DMy.index[100],DMy.index[200])
#         print DMy.iloc[:4,:4]    
        
        
        
        results = mantel(DMx, DMy, method, permutations=nPerm, alternative=alternative)
        corrCoef = results[0]
        p = results[1]
        print p
        
        
        mantelResDF.loc[0, 'featureName'] = feature_name
        mantelResDF.loc[0, 'phenotype_name'] = phenotype_name
        mantelResDF.loc[0, 'nPerm'] = nPerm
        mantelResDF.loc[0, 'method'] = method
        mantelResDF.loc[0, 'alternative'] = alternative
        mantelResDF.loc[0, 'corrCoef'] = corrCoef
        mantelResDF.loc[0, 'p'] = p
    
        mantelResDF.to_pickle(mantelFile)
        
        print 'done mantel'
    else:
        print 'This mantel test result already exist in folder'
                                           
                                           
    return mantelResDF

#---------------------------------------------------------------------------------------
def gen_dummies(df, toDummyColList):
    print 'original number of columns is %s' % len(df.columns)
    print 'number of toDummy columns is %s' % len(toDummyColList)
    for n, col in enumerate(toDummyColList):
        dummy = pd.get_dummies(df[col], prefix=col).iloc[:, 1:]
        if n == 0:
            df_new = pd.merge(df, dummy, how='left', left_index=True, right_index=True)
        else:
            df_new = pd.merge(df_new, dummy, how='left', left_index=True, right_index=True)
    print 'number of columns after dummy generation is %s' % len(df_new.columns)

    return df_new

#-----------------------------------------------------------------------------------------

def PCAfunc(TCRdf, n_comp, isSparse, ax, groupingByDF=None, groupbyName=None):
    from sklearn.decomposition import PCA
    from sklearn.decomposition import TruncatedSVD
    print 'TCRdf shape is %s_%s' % (TCRdf.shape[0], TCRdf.shape[1])
    if isSparse:
        pca = TruncatedSVD(n_comp, random_state=42)
    else:
        pca = PCA(n_comp)
    pca.fit(TCRdf.T)
    print 'PCA explained variance:'
    print pca.explained_variance_[:10]
    expVar = (pca.explained_variance_)[:10]
    print expVar
    PCAdf = pd.DataFrame(pca.components_).T
    PCAdf.index = TCRdf.index
    
    for column in PCAdf.columns:
        PCAdf = PCAdf.rename(columns={column:'PC' + str(column + 1)})
    PCAdf['BDindex'] = PCAdf.index.str.replace('BD', '').astype(int)
    if groupingByDF is None:
        PCAdf['isCardio'] = np.where(PCAdf['BDindex'] > 949, 1, 0)
        groupbyName = 'isCardio'
    else:
        PCAdf = pd.merge(PCAdf, pd.DataFrame(groupingByDF), how='left', left_index=True, right_index=True)
    for column in PCAdf.columns:
        if 'PC' in column:
            data = {}
            for name, group in PCAdf.groupby(groupbyName):
                data[name] = list(group[column])
            s_t, p_t = stats.ttest_ind(data[0], data[1])
            s_k, p_k = stats.ks_2samp(data[0], data[1])
            print 'p_value ttest for %s is %s' % (column, round(p_t, 6))
            print 'p_value kstest for %s is %s' % (column, round(p_k, 6))
            if column == 'PC1':
                p_ttest_PC1 = p_t
            if column == 'PC2':
                p_ttest_PC2 = p_t
    print PCAdf.head()
          
    ax.scatter(PCAdf.PC1, PCAdf.PC2, s=100, alpha=0.3)
    for n, ind in enumerate(PCAdf.index):
        BDind = int(ind.replace('BD', ''))
        if BDind > 949:
            xpoint = PCAdf.loc[ind, 'PC1']
            ypoint = PCAdf.loc[ind, 'PC2']
            ax.scatter(xpoint, ypoint, color='red', s=100, alpha=0.3)
#             plt.annotate(ind,xy=(xpoint,ypoint),xytext=(xpoint*1.05,ypoint*1.05),fontsize=8)
    plt.show()
    
    return PCAdf, ax, p_ttest_PC1, p_ttest_PC2

## the following is a new PCA function:
def calc_PCA(df,isSparse=False,n_comp=10,sample_list_list=None,color_list=None,pca_n1_toplot=0,pca_n2_toplot=1,toScale=False,
             toPlot=True, toAnnotate=False,fig=None,ax=None,calculateSeperation=True):
    
    ## sample_list_list = list of tuples, each tuple is composed of a string which is the name of the sample list, and the sample list
    from sklearn.decomposition import PCA
    from sklearn.decomposition import TruncatedSVD
    
    ### scale data:
    if toScale:
        print 'scaling data'
        df=pd.DataFrame(index=df.index, columns=df.columns, data=preprocessing.scale(df,copy=False))
        print df.iloc[:5,:5]
    
        
    ###calcualte pca:
    print 'df shape is %s_%s' %(df.shape[0],df.shape[1])
    print 'calculating PCA...'
    if isSparse:
        pca= TruncatedSVD(n_comp, random_state=42)
    else:
        pca = PCA(n_comp)
    pca.fit(df.T)
    expVar=(pca.explained_variance_)[:10]
    print ('PCA explained variance: ',expVar) 
    
    PCAdf=pd.DataFrame(pca.components_).T
    PCAdf.index=df.index
                        
                        
    if toPlot:
        print 'plotting...'
        if ax is None:
            fig,ax=plt.subplots()
        for n,item in enumerate(sample_list_list):
            ax.scatter(PCAdf.loc[item[1], pca_n1_toplot], PCAdf.loc[item[1], pca_n2_toplot],
                       color=color_list[n], alpha=1, s=100, label=item[0])

        ax.set_title('PCA', fontsize=30)
        ax.set_xlabel('PC%s (exp. var.=%0.2f)' %(pca_n1_toplot, pca.scores().proportion_explained[pca_n1_toplot]), fontsize=20)
        ax.set_ylabel('PC%s (exp. var.=%0.2f)' %(pca_n2_toplot, pca.scores().proportion_explained[pca_n2_toplot]), fontsize=20)
    #     ax.set_xticks(fontsize=15)
    #     ax.set_yticks(fontsize=15)
        ax.legend(bbox_to_anchor=(1.01, 0.95),loc='upper left', fontsize=12)
        if toAnnotate:
            print 'annotating'
            for item in sample_list_list:
                for sample in item[1]:
                    x=PCAdf.loc[sample,pca_n1_toplot]
                    y=PCAdf.loc[sample,pca_n2_toplot]
                    print (sample,x,y)
                    ax.annotate(sample,(x,y),xycoords='data',fontsize='medium')
    
                        
    ### calculate pca seperation:
    if calculateSeperation:
        print 'calculate PC seperation...'
        if len(sample_list_list)==2:
            for n in range(10):
                try:
                    item1_data=PCAdf.loc[sample_list_list[0][1], n].tolist()
                    item2_data=PCAdf.loc[sample_list_list[1][1], n].tolist()
                    MW_s1, p1 = mannwhitneyu(item1_data,item2_data)
                    print ('%s mean (PC%s):' %(sample_list_list[0][0],n), np.mean(item1_data))
                    print ('%s mean (PC%s):' %(sample_list_list[1][0],n), np.mean(item2_data))
                    print MW_s1, p1
                    print ('proportion explained: ',pca.scores().proportion_explained[n] )
                except:
                    print ' couldnt caclculate MW'
        else:
            print 'number of groups is %s' %len(sample_list_list)
        
    
    
    
    return PCAdf,fig,ax
#--------------------------------------------------------------------------------------------
def calc_catPhens_TCRfeatures_associations_ttest(featureDF, featureDFname, phenDF, phenDFname, FDR=0.1, testType='ttest'):
    feature_phen_ttest = pd.DataFrame()
    
    phenDF = phenDF.loc[featureDF.index, :]
    
    if phenDFname is not None:
        phenDF = pd.DataFrame(phenDF[phenDFname])
    
#     print ('featureDF shape is', featureDF.shape)
#     print ('fisrt samples in featureDF are:', featureDF.index[:4])
#     print ('phenDF shape is', phenDF.shape)
#     print ('fisrt samples in phenDF are:', phenDF.index[:4])
    
    for n, feature in enumerate(featureDF.columns.tolist()):
        featureDF[feature] = featureDF[feature].fillna(0)
        for k, phen in enumerate(phenDF.columns):
            df = pd.DataFrame()
            phen_feat_df = pd.merge(pd.DataFrame(featureDF[feature]), pd.DataFrame(phenDF[phen]), how='inner', left_index=True, right_index=True)
            data = {}
            lend = {}
            names = []
            for name, group in phen_feat_df.groupby(phen):
    #             print name
                data[name] = group[feature]
                lend[name] = len(group[feature])
                names.append(name)

            if testType == 'ttest':
                s, p = stats.ttest_ind(data[names[0]], data[names[1]])
            else:
                s, p = stats.mannwhitneyu(data[names[0]], data[names[1]])   
            df.loc[0, 'Feature'] = feature
            df.loc[0, 'Phenotype'] = phen
            df.loc[0, 'testType'] = testType
            df.loc[0, 'statistics'] = s
            df.loc[0, 'p'] = p
            df.loc[0, 'n0'] = lend[names[0]]
            df.loc[0, 'n1'] = lend[names[1]]

            feature_phen_ttest = pd.concat([feature_phen_ttest, df])
        
    nTests = len(feature_phen_ttest) 
    feature_phen_ttest = add_corrected_pValues(feature_phen_ttest, 'p', nTests, FDR)
    feature_phen_ttest = feature_phen_ttest.sort_values(by='p')
    
    return feature_phen_ttest



def calc_numPhens_TCRfeatures_associations_correlation(featureDF, featureDFname, phenDF, phenDFname, testType='pearson', FDR=0.1):
    feature_phen_corr = pd.DataFrame()
    
    phenDF = phenDF.loc[featureDF.index, :]
    
#     print ('featureDF shape is', featureDF.shape)
#     print ('fisrt samples in featureDF are:', featureDF.index[:4])
#     print ('phenDF shape is', phenDF.shape)
#     print ('fisrt samples in phenDF are:', phenDF.index[:4])
    
    for n, feature in enumerate(featureDF.columns):
#         featureDF[feature]=featureDF[feature].fillna(0)
        for k, phen in enumerate(phenDF.columns):

            df = pd.DataFrame()
#             print 'featureDF:'
#             print featureDF[feature].head()
#             print 'phenF:'
#             print phenDF[phen].head()
            
            if testType == 'pearson':
                r, p = MyPearsonr(featureDF[feature], phenDF[phen])
            else:
                r, p = MySpearmanr(featureDF[feature], phenDF[phen])
            df.loc[0, 'Feature'] = feature
            df.loc[0, 'Phenotype'] = phen
            df.loc[0, 'r'] = r
            df.loc[0, 'p'] = p
            df.loc[0, 'testType'] = testType

            feature_phen_corr = pd.concat([feature_phen_corr, df])
        
    nTests = len(feature_phen_corr) 
    feature_phen_corr = add_corrected_pValues(feature_phen_corr, 'p', nTests, FDR)
    feature_phen_corr = feature_phen_corr.sort_values(by='p')
    
    return feature_phen_corr

def calc_PCoA(df,metric,sample_list_list,color_list,pcoa_n1_toplot=0,pcoa_n2_toplot=1,toScale=False,
             toAnnotate=False,fig=None,ax=None,calculateSeperation=True):
    
    ## sample_list_list = list of tuples, each tuple is composed of a string which is the name of the sample list, and the sample list
    from skbio.stats.ordination import PCoA
    from skbio.stats.distance import DistanceMatrix
    
    ### scale data:
    if toScale:
        print 'scaling data'
        df=pd.DataFrame(index=df.index, columns=df.columns, data=preprocessing.scale(df,copy=False))
        print df
    
    ###generate distance matrix:
    print 'generating distance matrix...'
    df_condensed_org, distMat_square=genDistMat(df,metric)
    distMat_square_distmat=DistanceMatrix(distMat_square.values,distMat_square.index)
    
    ###calcualte PCoA:
    print 'calculate PCoA'
    pcoa = PCoA(distMat_square_distmat)
    pcoa_df = pd.DataFrame(pcoa.scores().site, index=distMat_square_distmat.ids, columns=range(pcoa.scores().site.shape[1]))
    print pcoa_df
    ###plot PCoA:
    print 'plotting'
   
    if ax is None:
        fig,ax=plt.subplots()
    for n,item in enumerate(sample_list_list):
        ax.scatter(pcoa_df.loc[item[1], pcoa_n1_toplot], pcoa_df.loc[item[1], pcoa_n2_toplot],
                   color=color_list[n], alpha=1, s=100, label=item[0])
   
    ax.set_title('PCoA', fontsize=30)
    ax.set_xlabel('PC%s (exp. var.=%0.2f)' %(pcoa_n1_toplot, pcoa.scores().proportion_explained[pcoa_n1_toplot]), fontsize=20)
    ax.set_ylabel('PC%s (exp. var.=%0.2f)' %(pcoa_n2_toplot, pcoa.scores().proportion_explained[pcoa_n2_toplot]), fontsize=20)
#     ax.set_xticks(fontsize=15)
#     ax.set_yticks(fontsize=15)
    ax.legend(bbox_to_anchor=(1.01, 0.95),loc='upper left', fontsize=12)
    
    ### calculate PCoA seperation:
    if calculateSeperation:
        if len(sample_list_list)==2:
            for n in range(10):
                try:
                    item1_data=pcoa_df.loc[sample_list_list[0][1], n].tolist()
                    item2_data=pcoa_df.loc[sample_list_list[1][1], n].tolist()
                    MW_s1, p1 = mannwhitneyu(item1_data,item2_data)
                    print ('%s mean (PCo%s):' %(sample_list_list[0][0],n), np.mean(item1_data))
                    print ('%s mean (PCo%s):' %(sample_list_list[1][0],n), np.mean(item2_data))
                    print MW_s1, p1
                    print ('proportion explained: ',pcoa.scores().proportion_explained[n] )
                except:
                    print ' couldnt caclculate MW'
        else:
            print 'number of groups is %s' %len(sample_list_list)
        
    if toAnnotate:
        print 'annotating'
        for item in sample_list_list:
            for sample in item[1]:
                x=pcoa_df.loc[sample,pcoa_n1_toplot]
                y=pcoa_df.loc[sample,pcoa_n2_toplot]
                print (sample,x,y)
                ax.annotate(sample,(x,y),xycoords='data',fontsize='medium')
    
    
    return pcoa_df,fig,ax,df_condensed_org

########################################################
def compare_TCR_between_groups(sample_list1,groupName1,sample_list2,groupName2,featureDF,TCRdf_binary):
    
    
    ### this function takes TCRdf_binary and featureDF and compare the sequences and features between
    ### two groups whose samples are included in these dfs (for example, PNP530 and Cardio)
    ####### the function was written originally for Ravid's data - check that it generalize well!#############
    
    
    #compare all TCR features:
    print 'comparing features among sample groups'
    print 'number of samples=%s' %len(featureDF.columns)
    featureDF1=featureDF.loc[sample_list1,:]
    featureDF2=featureDF.loc[sample_list2,:]
    
    feature_comparison_df=pd.DataFrame()
    for n,feature in enumerate(featureDF.columns):
        print n,feature
        data1=featureDF1[feature].tolist()
        data2=featureDF2[feature].tolist()
        mean1=np.mean(data1)
        mean2=np.mean(data2)
        try:
            MW_s, MW_p=mannwhitneyu(data1,data2)
        except:
            MW_s=999999; MW_p=1
        feature_comparison_df.loc[feature,'mean_%s' %groupName1]=mean1
        feature_comparison_df.loc[feature,'mean_%s' %groupName2]=mean2
        feature_comparison_df.loc[feature,'MW_s']=MW_s
        feature_comparison_df.loc[feature,'MW_p']=MW_p
        
    feature_comparison_df=add_corrected_pValues(feature_comparison_df,pValueColumn='MW_p',nTests=feature_comparison_df.shape[0],FDR=0.01)
    feature_comparison_df=feature_comparison_df.sort_values(by='MW_p')
    
    #compare TCR sequences between groups:
    
    print 'comparing TCR seqs among sample groups'
    print 'number of seqs=%s' %len(TCRdf_binary.columns)
    
    TCRdf_binary['group']=np.where(TCRdf_binary.index.isin(sample_list1),0,
                          np.where(TCRdf_binary.index.isin(sample_list2),1,np.nan))
    
    TCR_comparison_df=pd.DataFrame()
    for n,seq in enumerate(TCRdf_binary.columns[:-1]):
        print n,seq
        tab = pd.crosstab(TCRdf_binary[seq],TCRdf_binary['group'], dropna=False).fillna(0)
        seqAbs_group0=tab.iloc[0,0]
        seqAbs_group1=tab.iloc[0,1]
        seqPres_group0=tab.iloc[1,0]
        seqPres_group1=tab.iloc[1,1]
        try:
            fisher_OR, fisher_p = fisher_exact(tab, alternative='two-sided')
        except:
            print ('couldnt execute fisher test for seq %s' %seq)
            print (tab)
            fisher_OR= 999999; fisher_p=1
                
        TCR_comparison_df.loc[seq,'seqAbs_group0']=seqAbs_group0
        TCR_comparison_df.loc[seq,'seqAbs_group1']=seqAbs_group1
        TCR_comparison_df.loc[seq,'seqPres_group0']=seqPres_group0
        TCR_comparison_df.loc[seq,'seqPres_group1']=seqPres_group1
        TCR_comparison_df.loc[seq,'fisher_OR']=fisher_OR
        TCR_comparison_df.loc[seq,'fisher_p']=fisher_p

    TCR_comparison_df=add_corrected_pValues(TCR_comparison_df,pValueColumn='fisher_p',nTests=feature_comparison_df.shape[0],FDR=0.01)
    TCR_comparison_df=TCR_comparison_df.sort_values(by='fisher_p')
    
        
    return feature_comparison_df,TCR_comparison_df
        



