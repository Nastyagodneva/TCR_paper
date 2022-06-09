from os import listdir, mkdir, makedirs, remove
from os.path import isfile, join, isdir, exists
import pandas as pd
import numpy as np
from scipy import stats
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from ShaniBA.myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot, draw_correlation_scatter
from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
from Bio.SeqUtils import GC
import seaborn as sns
import random
from scipy.stats import pearsonr
from skbio.diversity.alpha import shannon, simpson, berger_parker_d

from ShaniBA.pop_organize import get_sample_data, get_sample_with_dfs
from ShaniBA.SufficientStatistics import *
from ShaniBA.MyFunctionsShani import *
import math
from skbio.stats.distance import mantel
from scipy.spatial.distance import braycurtis, pdist
# from tunneltoamazondb import getengine
from pandas import concat, read_csv, Series, read_sql
from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import * 


from ShaniBA.TCR_microbiome_interactions.TCR_microbiome_interactions_functions import *
from ShaniBA.TCR_feature_generation.TCR_feature_generation_functions import *

MyPath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF'

#------------------------------------------------------------------------------
'''
the following function takes a dataset name and folder, and generates a sharing matrix according to the minSharedList indicated.
during the process, a few other outputs are generated:
    df with all unique aa sequences is generated for each sample
    a summarizing df with all aa sequences and their counts 
    some statistics (currently just printed out and not saved

the sharing matrix is saved in the 'sharingAnalysis' folder in the dataset folder
NOTE: if the pivot table is too big, it is not saved and an exception is printed.
NOTE2: for the full PNP cohort, 'corrected' samples should be used


input: 
datasetName: string
dataset_folder - path, make sure it is corresponding to the name!!
minNsharedList - list of integers. if only 1 number is desired, include it in a list was well.
extractUniqueAA - true/false. if dfs were already extracted, change to flase to save time
etSharingStatistics - true/false

usage example:
datasetName='PNP515_ss12500_rep1'
dataset_folder='%s/TCR_real_data/SubSampled12500 data' %MyPath #change as necessary
minNsharedList=[2,50]

gen_sharingMatrix_perDataset(datasetName,dataset_folder,minNsharedList,extractUniqueAA=True,getSharingStatistics=True)
'''


def gen_sharingMatrix_perDataset(datasetName, dataset_folder, minNsharedList, sampleList=None, sampleListName=None, extractUniqueAA=True, getSharingStatistics=True, onlyProductive=True):
    
    # (1) define files and folders:
    try:
        sample_folder = '%s/SamplesForAnalysis_corrected' % dataset_folder
    except:
        sample_folder = '%s/SamplesForAnalysis' % dataset_folder
    if sampleList is None:
        sharingDFfolder = '%s/sharingDFs' % dataset_folder
        sharingAnalysisFolder = '%s/sharingAnalysis' % dataset_folder
    else:
        print 'only samples included in %s sample list are included in the analysis' % sampleListName
        sharingDFfolder = '%s/%s_sharingDFs' % (dataset_folder, sampleListName)
        sharingAnalysisFolder = '%s/%s_sharingAnalysis' % (dataset_folder, sampleListName)
     
    if not isdir(sharingDFfolder):
        makedirs(sharingDFfolder)
    if not isdir(sharingAnalysisFolder):
        makedirs(sharingAnalysisFolder)
    
    # (2) #get sample files:
    sharingFiles = [f for f in listdir(sharingDFfolder) if isfile(join(sharingDFfolder, f))]
    samples = [f for f in listdir(sample_folder) if isfile(join(sample_folder, f))]
    
    if sampleList is not None:
        print 'getting all samples from the sampleList...'
        # get 'clean name' samples to enable comparison
        cleanSamples = editSampleNamesList(samples)
        cleanSampleList = editSampleNamesList(sampleList)

        filteredSamples = []
        for n, sample in enumerate(samples):
            cleanSample = cleanSamples[n]
            if cleanSample in cleanSampleList:
                filteredSamples.append(sample)
        samples = filteredSamples

    
    # (3)extract unique aa sequences from all samples:
    if extractUniqueAA:
        print 'extracting unique aa sequenecs from all samples...(n=%s)' % len(samples)
        for n, sample in enumerate(samples):
        #     if n<5:
                print n, sample
                uniqueAAfile = '%s/%s' % (sharingDFfolder, sample)
                if uniqueAAfile not in sharingFiles:
                    sample_df = pd.read_table('%s/%s' % (sample_folder, sample))
                    sample_df = sample_df.rename(columns={'count (templates/reads)':'count (templates)', 'count (reads)':'count (templates)'})
                    sample_df['prod_stat'] = np.where(sample_df['sequenceStatus'] == 'In', 1, 0)
                    f = {'frequencyCount (%)':'sum', 'prod_stat':'mean'}
                    uniqAA = sample_df.groupby('aminoAcid').agg(f)
        #             uniqAA=uniqAA.drop(('frequencyCount (%)','count'),axis=1)
                    newSampleName=sample.split('_')[0]
#                     newSampleName=newSampleName.split('a')[0]
                    newSampleName=newSampleName.split('b')[0]
                    newSampleName=newSampleName.split('.')[0]
                    uniqAA['Sample'] = newSampleName
                    uniqAA.to_pickle(uniqueAAfile)
        #             print uniqAA.head()
                else:
                    print('found uniqueAAdf for this sample...')
    else:
        print 'skipping unique aa sequences extraction...'
                
    # (4) generate df with all sequences and nShared counts:
    print 'generating a df with all sequences and number of samples shared...'
    # concatenate all uniqueAA lists from all samples:
    AllUnique = concat_summarizing_dfs(sharingDFfolder)
    print AllUnique.head()
        
    # count number of repeats:
    nSharedCount = AllUnique.index.value_counts()
    # merge with former df:
    AllUniqueWithCounts = pd.merge(AllUnique, pd.DataFrame(nSharedCount), how='left', left_index=True, right_index=True)
    # edit df:
    AllUniqueWithCounts = AllUniqueWithCounts.rename(columns={'aminoAcid':'nShared'})
    AllUniqueWithCounts['isPublic'] = np.where(AllUniqueWithCounts['nShared'] > 1, 1, 0)
    # save AllUniqueWithCounts
    file1 = '%s/AllUniqueWithCounts' % sharingAnalysisFolder
    AllUniqueWithCounts.to_pickle(file1)
    
    if onlyProductive:
        print 'USING ONLY PRODUCTIVE SEQUENCES TO GENERATE SHARING MATRIX!!!'
        AllUniqueWithCounts = AllUniqueWithCounts[AllUniqueWithCounts['prod_stat'] == 1]
        print 'filtering out strange sequences:...'
        AllUniqueWithCounts['seqLen'] = AllUniqueWithCounts.index.str.len()
        AllUniqueWithCounts[~((AllUniqueWithCounts['seqLen'] < 3) & (AllUniqueWithCounts['nShared'] < 3))]
    
    
    # (5) get some sharing sequence statistics:
    if getSharingStatistics:
        print 'maximal number of shared samples per sequence for dataset %s (filtered by sample List %s) is %s' % (datasetName, sampleListName, AllUniqueWithCounts['nShared'].max())
        noDups = AllUniqueWithCounts.reset_index().drop_duplicates(subset='index')
        print 'percPublic=%s' % (noDups['isPublic'].mean() * 100)
        print 'mean number of shared samples per sequence is %s' % noDups['nShared'].mean()
        print 'mean percPublic per sample=%s' % AllUniqueWithCounts.groupby('Sample')['isPublic'].apply(lambda x: np.mean(x) * 100).mean()
        print 'mean nShared per sample=%s' % AllUniqueWithCounts.groupby('Sample')['nShared'].apply(lambda x: np.mean(x)).mean()
    
    # get only sequences shared by more than minNshared:
    print 'total number of unique aa seqeunces in the dataset is %s' % AllUniqueWithCounts.index.nunique()
    for minNshared in minNsharedList:
        print 'generating sharing matrix for minNshared=%s...' % minNshared
        sharedSeqs = AllUniqueWithCounts[AllUniqueWithCounts['nShared'] >= minNshared]
        print 'number of unique sequences shared by %s samples or more is %s' % (minNshared, sharedSeqs.index.nunique())
        sharedSeqsClean = sharedSeqs.reset_index().rename(columns={'index':'Sequence'})
        
        # generate pivot tables:
        sharingMatrix_RA = sharedSeqsClean.pivot(index='Sample', columns='Sequence', values='frequencyCount (%)')
        sharingMatrix_RA = sharingMatrix_RA.fillna(0)
        sharingMatrix_RA = sharingMatrix_RA = editSampleNames(sharingMatrix_RA)
        print 'isPublic status is %s and should be 1' % sharedSeqsClean['isPublic'].mean()
        print 'sharingMatrix_RA shape is %s_%s' % (sharingMatrix_RA.shape[0], sharingMatrix_RA.shape[1])
        
        # divide all values by 100 in order to get real RA values and not percents
        sharingMatrix_RA = sharingMatrix_RA.div(100)
        
        RA_file = '%s/sharingMatrix_%s_minNshared%s_RA_onlyProductive%s' % (sharingAnalysisFolder, datasetName, minNshared, onlyProductive)

        try:
            sharingMatrix_RA.to_pickle(RA_file)
        except SystemError:
            print '*****df is too big and was not saved****'
    
    return sharingMatrix_RA, RA_file 

#--------------------------------------------------------------------------------------------

'''
the following function generate pie diagram showing the distribution of different TCR identities in a list.

input:
*ax
*df - dataframe containing the sequences as index and at least the 'column' columns which includes TCR seq identity annotations according to
one of the databases (usually 'Epitope species_VDJDB'/'Pathology_McPAS', can be also 'Category_McPAS'. consider generating a column that summarizes identities from different databases
MAKE SURE THIS DF is depleted of whitespaces.
*nSharedThreshold - None or int>1, indicate that minimal number of samples shared by a sequence required to consider this sequence in the 
analysis. if this input is not none, there should be a column named 'nShared' which contains integers reflecting number of shared samples per
sequence
*useMore - True/False. relevant only if nSharedThreshold is not None, and indicates whether to consider all sequence shared by more than 
nSharedThreshold (True) or shared by less (True)
*column - string. name of column holding the identity information
*dropna - True/False - whether or not to show unrecognized identities (usually >95% of the sequences are un-identified so the use True)
*size = 



'''
def plot_identity_pie_plot(ax, df, dfName, nSharedThreshold, useMore, column, dropna,ThresoldToOther):
    
    # (1) DEFINE PARAMETERS:
    if nSharedThreshold is None:
        df2 = df.copy()
        if dropna:
            nSeq = len(df2[df2[column].notnull()])
        else:
            nSeq = len(df2[column])
        title = '%s (%s sequences)\nTCR identity distribution (%s)' % (dfName, nSeq, column)
    else:
        if useMore:
            df2 = df[df['nShared'] >= nSharedThreshold]
            sharedBy = 'more'
        else:
            df2 = df[df['nShared'] <= nSharedThreshold]
            sharedBy = 'less'
        if dropna:
            nSeq = len(df2[df2[column].notnull()])
        else:
            nSeq = len(df2[column])
        title = 'TCR identity distribution (%s)\nSequences shared by %s or %s samples\n(%s sequences)' % (column, nSharedThreshold, sharedBy, nSeq)
    
    # (2) CALCULATE NORMALIZED NUMBER OF IDENTITITES AND PLOT
    valueCounts = pd.DataFrame(df2[column].value_counts(dropna=dropna, normalize=True))
    if ThresoldToOther is None:
        valueCounts['count2'] = np.where(valueCounts[column] > valueCounts[column].min() * 2,
                                             valueCounts.index, 'Other')
    else:
        valueCounts['count2'] = np.where(valueCounts[column] > ThresoldToOther,
                                             valueCounts.index, 'Other')
    valueCounts['count2'] = valueCounts['count2'].astype(str)
    valueCountsgrouped = valueCounts.groupby('count2').sum()
    # print valueCountsgrouped.head()
    valueCountsgrouped = valueCountsgrouped.sort_values(by=column, ascending=False)
    cmap = plt.cm.get_cmap('Spectral')
    colors = [cmap(x * 20) for x in range(len(valueCountsgrouped))]

    patches, texts, autotexts=ax.pie(valueCountsgrouped[column], labels=valueCountsgrouped.index, autopct='%.2f', colors=colors, 
           textprops={'fontsize': 18})
    for a in autotexts:
        a.set_size('xx-large')
    for t in texts:
        t.set_size('xx-large')
#     valueCountsgrouped[column].plot.pie(ax=ax,autopct='%.2f')
    ax.set_ylabel('')
    ax.set_title(title, fontsize=20)
#     ax.legend(fontsize=16)
    
    # (3) CALCULATE NON-NORMALIZED NUMBER OF IDENTITES TO ENABLE RETURNING COUNTS:
    valueCounts2 = pd.DataFrame(df2[column].value_counts(dropna=dropna, normalize=False))
    totalCounts=valueCounts2[column].sum()
    
    valueCounts2['count2'] = np.where(valueCounts2[column] > ThresoldToOther*totalCounts, valueCounts2.index, 'Other')
    print 'ThresholdToOther is %s, total counts=%s, final threshold for chi square is %s' %(ThresoldToOther,totalCounts,
                                                                                            ThresoldToOther*totalCounts)
#     if ThresoldToOther is None:
#         valueCounts2['count2'] = np.where(valueCounts2[column] > valueCounts2[column].min() * 2,
#                                              valueCounts2.index, 'Other')
#     else:
#         valueCounts2['count2'] = np.where(valueCounts2[column] > ThresoldToOther,
#                                              valueCounts2.index, 'Other')
    valueCounts2['count2'] = valueCounts2['count2'].astype(str)
    valueCountsgrouped2 = valueCounts2.groupby('count2').sum()

    
    return ax, valueCountsgrouped2

#-----------------------------------------------------------------------------------------------------

'''
the following function is an helper function for 'comparePublicSeqIdentities'
'''
def addSeqIdentitiesAndPlotPie(publicDF, publicDFName, publicAnalysisFolder,identsProcessedDF, ax, nSharedThreshold,
                       useMore, column, dropna,ThresholdToOther):
    
#     publicDF_WithIdents=pd.merge(publicDF,identsDF,how='left',left_index=True,right_index=True)
    publicDF_WithIdents_processed = pd.merge(publicDF, identsProcessedDF, how='left', left_index=True, right_index=True)
    f1 = '%s/%swithIdents.xlsx' % (publicAnalysisFolder, publicDFName)
#     publicDF_WithIdents.to_excel(f1)
    f2 = '%s//%swithIdents_processed.xlsx' % (publicAnalysisFolder, publicDFName)
    publicDF_WithIdents_processed.to_excel(f2)
    df = publicDF_WithIdents_processed
    dfName = publicDFName
    
    ax, valueCountsgrouped2 = plot_identity_pie_plot(ax, df, dfName, nSharedThreshold, useMore, column, dropna,
                                                     ThresholdToOther)
    
    return ax, valueCountsgrouped2, publicDF_WithIdents_processed

'''
the following funciton takes a dataset List and generate the following information:
1. distributions of public sequence per samples (and comparison between datasets)
2. venn diagram of public sequences - which are common to all datasets and which only for specific dataset
3. pie charts of public sequence identities and comparison between datasets (also an excel file of chi square p-values for all 
comparisons
4.excel file with all public sequences and their identity if exist

input:
datasetList: a list of tuples. each tuple contains: (0) datasetName (string)(1) datasetFolder (string) (2) TCRdfPercShared (int or None
use None for all public sequenes (shared by two samples or more) (3) sampleListName 
make sure there is a TCRdf with the indicated percShared for each dataset!
*identityColumnForPie='Epitope species_VDJDB'/'Pathology_McPAS' - the column in the CDR3 identity table that should be used for 
the identity pie figure.
dropnaFromIdentitiesPie=True/False
nSharedMin- minimal number of shared samples per sequence in order to be included in the analysis FOR THE VENN DIAGRAM AND IDENTITY ANALYSIS!
(note that the public rate per sample distribution is based on the percShared defined within the dataset info)


NOTE! in the pie square, only annotations that appear in more than 2% of the samples are plotted and counted for chi square test (others are aggregated to 'others)
'''

def comparePublicSeqIdentities(datasetList, identityColumnForPie, dropnaFromIdentitiesPie=True, nSharedMin=2,
                               identityDF=None, identityDFname=None,ThresholdToOther=None,RemoveSamplesWithRelative=False):
   
    from scipy.stats import chi2_contingency
    
    # (1) loop 1: (a) plot distribution of public sequence number per sample in each dataset (b) generate df of all public 
    # sequences in each dataset.   
    uniqueSeqDFdict = {}
    percPublicDataList = []
    datasetNameList = []
    for n, dataset in enumerate(datasetList):            
        print 'step 1: collecting public sequence rates per sample in each dataset...'
        if dataset[3] is not None:
            datasetNameList.append('_'.join([dataset[0], dataset[3]]))
        else:
            datasetNameList.append(dataset[0])
        print n
        print 'dataset Name=%s' % dataset[0]
        print 'dataset file=%s' % dataset[1]
        print 'dataset sampleList=%s' % dataset[3]
        
        # load relevant TCRdf:
        TCRdfPercShared = dataset[2]
        print 'TCRdf percent shared is %s' % TCRdfPercShared
        if dataset[3] is not None:
            if TCRdfPercShared is not None:
                TCRdfFiles = '%s/%s_sharingAnalysis/sharingMatrix_%s_minNshared2_RA_onlyProductiveTrue__\
percShared%s_OLtrimmed_binary' % (dataset[1], dataset[3], dataset[0], dataset[2])
            else:
                TCRdfFiles = '%s/%s_sharingAnalysis/sharingMatrix_%s_minNshared2_RA_onlyProductiveTrue' % (dataset[1],
                                    dataset[3], dataset[0])
            TCRdf = pd.read_pickle(TCRdfFiles)
        else:
            if TCRdfPercShared is not None:
                try:
                    TCRdfFiles = '%s/sharingAnalysis/sharingMatrix_%s_minNshared2_RA_onlyProductiveTrue__\
percShared%s_OLtrimmed_binary' % (dataset[1], dataset[0], dataset[2])
                    TCRdf = pd.read_pickle(TCRdfFiles)
                except:
                    try:
                        TCRdfFiles = '%s/%s_sharingAnalysis/sharingMatrix_%s_minNshared2_RA_onlyProductiveTrue__\
percShared%s_OLtrimmed_binary' % (dataset[1], dataset[0], dataset[0], dataset[2])
                        TCRdf = pd.read_pickle(TCRdfFiles)
                    except:
                        print 'couldnt load TCRdf file'
            else:
                try:
                    TCRdfFiles = '%s/sharingAnalysis/sharingMatrix_%s_minNshared2_RA_onlyProductiveTrue' % (dataset[1], dataset[0])
                    TCRdf = pd.read_pickle(TCRdfFiles)
                except:
                    try:
                        TCRdfFiles = '%s/%s_sharingAnalysis/sharingMatrix_%s_minNshared2_RA_onlyProductiveTrue' % (dataset[1], dataset[0], dataset[0])
                        TCRdf = pd.read_pickle(TCRdfFiles)
                    except:
                        print 'couldnt load TCRdf file'

        TCRdf = editSampleNames(TCRdf)
        if RemoveSamplesWithRelative:
            TCRdf=removeRelatives(TCRdf)
                        
        # generate df with perc public per sample
        # loading number of unique seqeunce per sample:
        f3 = '%s/featureSummaryDFs/%s_allFeatures' % (dataset[1], dataset[0])
        featuresDF = pd.read_pickle(f3)
        featuresDF = editSampleNames(featuresDF)
        nUnique = pd.DataFrame(featuresDF['NT count_1'])
        nPublicPerSampleDF = pd.merge(pd.DataFrame(TCRdf.sum(axis=1)), nUnique, how='left', left_index=True, right_index=True)
        
        # calculate percentage of public sequences out of all unique and plot distributions of this statistic in both cohorts:
        nPublicPerSampleDF['percPublic'] = 100 * nPublicPerSampleDF[0].astype('float') / nPublicPerSampleDF['NT count_1']
        percPublicDataList.append(list(nPublicPerSampleDF['percPublic']))
        
        print 'step 2: generating merged public DF...'
        # get uniqueAA df per dataset:
        # load uniqueAAdf:
        if dataset[3] is not None:
            uniqueAAfile = '%s/%s_sharingAnalysis/AllUniqueWithCounts' % (dataset[1], dataset[3])
            AllUniqueWithCounts = pd.read_pickle(uniqueAAfile)
            print 'sharing analysis folder is called *%s_sharing analysis*' % dataset[3]
        else:
            try:
                uniqueAAfile = '%s/sharingAnalysis/AllUniqueWithCounts' % dataset[1]
                AllUniqueWithCounts = pd.read_pickle(uniqueAAfile)
                print 'sharing analysis folder is called *sharing analysis*'
            except:
                try:
                    uniqueAAfile = '%s/%s_sharingAnalysis/AllUniqueWithCounts' % (dataset[1], dataset[0])
                    AllUniqueWithCounts = pd.read_pickle(uniqueAAfile)
                    print 'sharing analysis folder is called *%s_sharing analysis*' % dataset[0]
                except:
                    print 'coutldnt load uniqueAA df'
        # get only sequences shared by more than nSharedMin times, remove unnecessary columns and then duplicates
        AllUniqueWithCounts = AllUniqueWithCounts[AllUniqueWithCounts['nShared'] >= nSharedMin]  # take only shared by more than nShared
        UniqueSequencesWithCounts = AllUniqueWithCounts.drop(['Sample', 'frequencyCount (%)', 'prod_stat', 'isPublic'], axis=1)
        UniqueSequencesWithCounts = UniqueSequencesWithCounts[~UniqueSequencesWithCounts.index.duplicated(keep='first')]   
        if dataset[3] is not None:
            uniqueSeqDFdict['_'.join([dataset[0], dataset[3]])] = UniqueSequencesWithCounts
        else:
            uniqueSeqDFdict[dataset[0]] = UniqueSequencesWithCounts
        print len(UniqueSequencesWithCounts)
        
        if n == 0:
            mergedPublicDF = UniqueSequencesWithCounts
        else:
            mergedPublicDF = pd.merge(mergedPublicDF, UniqueSequencesWithCounts, how='outer', left_index=True, right_index=True)
        if dataset[3] is not None:
            mergedPublicDF = mergedPublicDF.rename(columns={'nShared':'nShared_%s_%s' % (dataset[0], dataset[3])}) 
        else:
            mergedPublicDF = mergedPublicDF.rename(columns={'nShared':'nShared_%s' % dataset[0]})
    print 'merged public list length length=%s' % len(mergedPublicDF)
    
    print 'step 3: plotting public rates per samples for all datasets...'    
    # plot perc public per sample distribution:
    if dataset[3] is not None:
        publicAnalysisFolder = '%s/publicAnalysis/%s' % (datasetList[0][1], datasetList[0][3])
    else:
        publicAnalysisFolder = '%s/publicAnalysis' % datasetList[0][1]
    if not isdir(publicAnalysisFolder):
        makedirs(publicAnalysisFolder)
        
    folderToSave = publicAnalysisFolder
    fig1, ax = plt.subplots(figsize=(8, 6))
    title = 'Sharing rates per sample - comparison between cohorts\n(Shared \
among %s percent of samples or more)' % datasetList[0][2]
    showLegend = True
    nBins = 20
    dataList = zip(datasetNameList, percPublicDataList)    
    ax, ks_p_cohort1_cohort2, t_p_cohort1_cohort2, p_Anov, filename = plotHistComprison(dataList,ax, title, showLegend, nBins)

#     ax.set_title(filename)
    filename = filename + '_percShared%s' % datasetList[0][2]
    if RemoveSamplesWithRelative:
        filename=filename + 'relativedRemoved'
    plt.show()
    figFile = '%s/%s' % (publicAnalysisFolder, filename)
    fig1.savefig(figFile, dpi=200, bbox_inches="tight")

    print 'step 4: find  sequences that are public in all datasets, add identities and plot pie diagram...'
    # general definitions for pie diagrams super-figure:
    ncols = len(datasetList) + 1
    fig4, axes4 = plt.subplots(nrows=1, ncols=ncols, figsize=(10 * ncols, 9))
    if 'popped' in identityDFname:
        publicIdentsFigFile = '%s/publicIdentsFig_%s_%s_sharedBy%s_dropna%s_thOther%s_popped' % (publicAnalysisFolder, '_'.join(datasetNameList),
                                                                  identityColumnForPie, nSharedMin,dropnaFromIdentitiesPie,ThresholdToOther)
    else:
        publicIdentsFigFile = '%s/publicIdentsFig_%s_%s_sharedBy%s_dropna%s_thOther%s' % (publicAnalysisFolder, '_'.join(datasetNameList),
                                                                  identityColumnForPie, nSharedMin,dropnaFromIdentitiesPie,ThresholdToOther)
    publicIdentsFigFile=publicIdentsFigFile.replace('.','-')
    
    valueCountDict = {}
    nSharedThreshold = None
    useMore = True
    column = identityColumnForPie
    dropna = dropnaFromIdentitiesPie
    
    # load identity files:
#     identsFile='%s/TCR CDR3 sequence databases/CDR3identityTable_06082018.xlsx' %MyPath
    if identityDF is None:
        identsProcessedFile = '%s/TCR CDR3 sequence databases/CDR3identityTable_06082018_processed.xlsx' % MyPath
#         identsDF=pd.read_excel(identsFile)
        identsProcessedDF = pd.read_excel(identsProcessedFile)
    else:
        identsProcessedDF = identityDF
    
    # add identity columns to sequences shared by all groups, save df to file, plot pie and count values for each
    # category:
    publicDFInAll = mergedPublicDF[mergedPublicDF.notnull().all(axis=1)]
    publicDF = publicDFInAll
    publicDFName = 'publicAll'
    print 'number of sequence shared by %s or more samples in %s: %s (%s perc)' % (nSharedMin, publicDFName,
                                len(publicDF), 100 * float(len(publicDF)) / len(mergedPublicDF))
    ax4 = axes4.flatten()[0]
    ax4, valueCountsgrouped_all, publicDF_WithIdents_processed = addSeqIdentitiesAndPlotPie(publicDF, publicDFName,
                        publicAnalysisFolder, identsProcessedDF, ax4, nSharedThreshold, useMore, column, dropna,
                        ThresholdToOther)
    valueCountDict[publicDFName] = valueCountsgrouped_all
    
    # for each dataset, find sequences that are public only in it, add identities, save df to files, plot pie
    # and count values for each category:
    print 'step 5: find sequences that are public only in one dataset, add identities, save files and plot pie diagram...'
    count = 0
    chi_df = pd.DataFrame()
    for n, dataset in enumerate(datasetList): 
        if dataset[3] is not None:
            col = 'nShared_%s_%s' % (dataset[0], dataset[3])
        else:
            col = 'nShared_%s' % dataset[0]
        publicDFsingle = mergedPublicDF[(mergedPublicDF[col].notnull()) & (mergedPublicDF.notnull().sum(axis=1) == 1)]
        publicDF = publicDFsingle
        if dataset[3] is not None:
            publicDFName = 'publicOnlyIn%s_%s' % (dataset[0], dataset[3])
        else:
            publicDFName = 'publicOnlyIn%s' % dataset[0]
        print 'number of sequence shared by %s or more samples in %s: %s (%s perc)' % (nSharedMin, publicDFName,
                                len(publicDF), 100 * float(len(publicDF)) / len(mergedPublicDF))
        ax4 = axes4.flatten()[n + 1]
        ax4, valueCountsgrouped2, publicDF_WithIdents_processed = addSeqIdentitiesAndPlotPie(publicDF, publicDFName,
                publicAnalysisFolder, identsProcessedDF, ax4, nSharedThreshold, useMore, column, dropna,ThresholdToOther)
        
        print 'calculating chi square test comparing %s to public sequence common to all datasets' % publicDFName      
#         print 'valueCountsgrouped_all:'
#         print valueCountsgrouped_all
#         print 'valueCountsgrouped2:'
#         print valueCountsgrouped2
        contTable = pd.merge(pd.DataFrame(valueCountsgrouped_all), pd.DataFrame(valueCountsgrouped2),
                           how='outer', left_index=True, right_index=True)
        contTable = contTable.fillna(0)
        print contTable
        chi, p, dof, expctd = chi2_contingency(contTable.T)
        ax4.annotate('p_chi square tests\n(vs. public in all)=%s' % round(p, 6), xy=(0.96, 0), xycoords='axes fraction',
                fontsize=18, horizontalalignment='right', verticalalignment='top', fontweight='bold')
        
        for k, v in valueCountDict.items():
            contTable2 = pd.merge(pd.DataFrame(valueCountsgrouped2), pd.DataFrame(v), how='outer',
                                left_index=True, right_index=True)
            contTable2 = contTable2.fillna(0)
            chi2, p2, dof2, expctd2 = chi2_contingency(contTable2.T)
            chi_df.loc[count, 'dataset1Name'] = publicDFName
            chi_df.loc[count, 'dataset2Name'] = k
            chi_df.loc[count, 'p_chi'] = p2
            count = count + 1
        valueCountDict[publicDFName] = valueCountsgrouped2
    
    
    print 'step 6: finish pie diagram figure and chi_df generation...'
    
    if 'popped' in identityDFname:
        chi_df_File = '%s/chi_df_%s_%s_sharedBy%s_dropna%s_thOther%s_popped.xlsx' % (publicAnalysisFolder, '_'.join(datasetNameList),
                                                                  identityColumnForPie, nSharedMin,dropnaFromIdentitiesPie,ThresholdToOther)  
    else: 
        chi_df_File = '%s/chi_df_%s_%s_sharedBy%s_dropna%s_thOther%s.xlsx' % (publicAnalysisFolder, '_'.join(datasetNameList),
                                                                  identityColumnForPie, nSharedMin, dropnaFromIdentitiesPie,ThresholdToOther)
    chi_df_File=chi_df_File.replace('.','-')
    chi_df_File=chi_df_File.replace('-xlsx','.xlsx')
    chi_df.to_excel(chi_df_File)
    
    fig4.subplots_adjust(left=0.09, right=0.9, top=0.9, bottom=0.2, wspace=0.5, hspace=0.25)
    fig4.savefig(publicIdentsFigFile, dpi=200, bbox_inches="tight")
    plt.show()
    
    print 'step 7: plotting venn diagrams...'
    # plot venn diagram :   
    if len(datasetList) == 2:
        print 'plotting venn diagram for two groups:'
        from matplotlib_venn import venn2
        fig3, ax3 = plt.subplots(figsize=(8, 6))
        if datasetList[0][3] is not None:
            label1 = '_'.join([datasetList[0][0], datasetList[0][3]])
            label2 = '_'.join([datasetList[1][0], datasetList[1][3]])
            venn = venn2([set(uniqueSeqDFdict[label1].index.tolist()),
                        set(uniqueSeqDFdict[label2].index.tolist())],
                       set_labels=(label1, label2), ax=ax3)
            VennPlotFile = '%s/%s%s_%s%s_PublicSeq_VennPlot_nSharedMin%s' % (publicAnalysisFolder, datasetList[0][0], datasetList[0][3],
                datasetList[1][0], datasetList[1][3], nSharedMin)
        else:
            label1 = datasetList[0][0]
            label2 = datasetList[1][0]
            venn = venn2([set(uniqueSeqDFdict[label1].index.tolist()),
                        set(uniqueSeqDFdict[label2].index.tolist())],
                       set_labels=(datasetList[0][0], datasetList[1][0]), ax=ax3)
            VennPlotFile = '%s/%s_%s_PublicSeq_VennPlot_nSharedMin%s' % (publicAnalysisFolder, datasetList[0][0],
                datasetList[1][0], nSharedMin)
        ax3.set_title('Number of Public sequences in each cohort')
        fig3.savefig(VennPlotFile, dpi=200, bbox_inches="tight")
        plt.show()
        
    elif len(datasetList) == 3:
        print 'plotting venn diagram for two groups:'
        from matplotlib_venn import venn3
        fig3, ax3 = plt.subplots(figsize=(8, 6))
        if datasetList[0][3] is not None:
            label1 = '_'.join([datasetList[0][0], datasetList[0][3]])
            label2 = '_'.join([datasetList[1][0], datasetList[1][3]])
            label3 = '_'.join([datasetList[1][0], datasetList[2][3]])
            venn = venn3([set(uniqueSeqDFdict[label1].index.tolist()),
                        set(uniqueSeqDFdict[label2].index.tolist()),
                       set(uniqueSeqDFdict[label3].index.tolist())],
                       set_labels=(label1, label2, label3), ax=ax3)
            VennPlotFile = '%s/%s%s_%s%s_%s%s_PublicSeq_VennPlot_nSharedMin%s' % (publicAnalysisFolder, datasetList[0][0], datasetList[0][3],
                    datasetList[1][0], datasetList[1][3], datasetList[2][0], datasetList[2][3], nSharedMin)    
        else:
            label1 = datasetList[0][0]
            label2 = datasetList[1][0]
            label3 = datasetList[1][0]
            venn = venn3([set(uniqueSeqDFdict[label1].index.tolist()),
                        set(uniqueSeqDFdict[label2].index.tolist()),
                       set(uniqueSeqDFdict[label3].index.tolist())],
                       set_labels=(datasetList[0][0], datasetList[1][0], datasetList[2][0]), ax=ax3)
            VennPlotFile = '%s/%s_%s_%s_PublicSeq_VennPlot_nSharedMin%s' % (datasetList[0][0], datasetList[1][0],
                datasetList[2][0], publicAnalysisFolder, nSharedMin)
        ax3.set_title('Number of Public sequences in each cohort')
        fig3.savefig(VennPlotFile, dpi=200)
        plt.show()
    else:
        print 'number of datasets is %s and therefore venn plots were not drown' % len(datasetList)
        
    print 'end of function!!'
    
    
def removeRelatives(df,alsoCols=False):
    df2=df.copy()
    samplesToRemoveFile='%s/Sample files/SamplesWithRelativeInPNP530.xlsx' %MyPath
    samplesToRemoveDF=pd.read_excel(samplesToRemoveFile)
    samplesToRemove=samplesToRemoveDF['BD'].tolist()
    newIndex=[x for x in df2.index if x not in samplesToRemove]
    print 'removing samples that have relatives in PNP530'
    print 'n samples before removal=%s' %len(df.index)
    df2=df2.loc[newIndex,:]
    if alsoCols:
        df2=df2.loc[:,newIndex]
    print 'n samples after removal=%s' %len(df2.index)
    
    return df2
        
        
  
   
