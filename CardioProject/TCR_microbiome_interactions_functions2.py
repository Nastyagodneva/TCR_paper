import sys
sys.path.insert(0, '/net/mraid08/export/genie/workspace/Microbiome/ShaniBA')

from ShaniBA.myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot
from os import listdir, mkdir, makedirs
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
from scipy.stats import pearsonr, fisher_exact, spearmanr
from skbio.diversity.alpha import shannon, simpson, berger_parker_d
# 
# from pop_organize import get_sample_data, get_sample_with_dfs
# from SufficientStatistics import *
from ShaniBA.MyFunctionsShani import *
import math
# from myplots import roundup, rounddown, find_decimal_fold
from skbio.stats.distance import mantel
from scipy.spatial.distance import braycurtis, pdist
from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import *
from ShaniBA.TCR_feature_generation.publicSeqAnalysis import *
from ShaniBA.TCR_microbiome_interactions.TCR_microbiome_interactions_functions import *


import os
from Utils import cacheOnDisk
from queue.qp import qp, fakeqp
from addloglevels import sethandlers

MyPath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF'
# load PNP515 sample list:
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP515', 'rb') as fp:
    PNP515 = pickle.load(fp)
print 'done'

'''
the following funciton generates Mb/TCR table and process it:

USAGE EXAMPLE:

# parameters to change:
#general parameters:
DFtype='TCR'
genDF=False #False=generate new df
toBinary=True

#basic prameters for generating df:
mbLevel='g'
useShortName=True

datasetFolder='%s/TCR_real_data/SubSampled9000data_rep2' %MyPath
datasetName='PNP515_ss9000_rep2'

#parameters for manipulations on df 
minVal='dfMinVal'    #minVal can be None,0, float, or 'dfMinVal' or dfMinVal2:
minSharedT=None #minimal number of samples shared by seq/species required to leave sample in the database (int or None)
percShared=5 #minimal number of samples shared by seq/species required to leave sample in the database (int [ percentage]
                #or None)
removeOutliers=True
normData=True
logTransform=True


#constant parameters:
extractUniqueAA=True # use True when this is the first time to analyze this dataset, otherwise, use False
minNshared=2
onlyProductive=True

mbDataFolder='AllSeqProjects'
#load PNP515 sample list:
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP515','rb') as fp:
    PNP515=pickle.load(fp)
SampleList=PNP515
SampleListName='PNP515'
libPrepMethod=None
filterGenotek=True
filterMinimalReads=9000000
filterlibPrepMethod=libPrepMethod
groupFunction='noOutlierMean'
nSTD=5
nMinSamples=3
ignoreNotSameYear=True
removeSamePerson=False

df=genTCRorMBdfWithManipulations(DFtype,genDF,toBinary,removeOutliers,normData,logTransform,
                                 minVal,minSharedT,percShared,
                                 mbLevel,useShortName,datasetFolder,datasetName,extractUniqueAA,
                                minNshared,onlyProductive,mbDataFolder,SampleList,SampleListName,filterlibPrepMethod,
                                 filterGenotek, filterMinimalReads, nSTD,nMinSamples,ignoreNotSameYear,removeSamePerson)
'''

def genTCRorMBdfWithManipulations(DFtype, genDF, toBinary, removeOutliers, normData, logTransform,
                                 minVal, minSharedT, percShared,
                                 mbLevel, useShortName, datasetFolder, datasetName,extractUniqueAA,
                                minNshared, onlyProductive, mbDataFolder, SampleList,
                                  SampleListName, filterlibPrepMethod=None, filterGenotek=True,filterMinimalReads=True,
                                  groupFunction='noOutlierMean', nSTD=5, nMinSamples=3, ignoreNotSameYear=True, removeSamePerson=False):


     # # (1)load/generate df:

    if DFtype == 'Mb':
        if genDF:
            print 'generating MbDF...'
            df, filename = genFilterMergeMbData(mbLevel, mbDataFolder, SampleList, SampleListName, filterGenotek, filterMinimalReads,
                                    filterlibPrepMethod, groupFunction, nSTD, nMinSamples, ignoreNotSameYear, useShortName)
        else:
            folderToSave = '%s/MicrobiomeDataTables/FilteredAndMergedOnBD' % MyPath
            try:
                filename = '%s/MPA_%s_standardParams' % (folderToSave, mbLevel)
                df = pd.read_pickle(filename)
                print 'loading existing MbDF: level=%s, standard params...' % mbLevel
            except:
                filename = '%s/MPA_%s_%s_%s_fGenotek%s_fMinimalRead%s_libMeth%s_igAnotherYear%s_%snSTD%snMinSamples%s_fSamePerson%s' % (folderToSave,
                    mbLevel, mbDataFolder, SampleListName, filterGenotek, filterMinimalReads, filterlibPrepMethod, ignoreNotSameYear,
                    groupFunction, nSTD, nMinSamples, removeSamePerson)
                try:
                    df = pd.read_pickle(filename)
                    print 'loading existing MbDF: level=%s, long name...' % mbLevel
                except:
                    print filename
                    print 'file was not found...'
    elif DFtype == 'TCR':
        if genDF:
            minNsharedList = [int(x) for x in str(minNshared)]
            print 'preparing sharing Matrix'
            df, filename = gen_sharingMatrix_perDataset(datasetName, datasetFolder, minNsharedList, SampleList,SampleListName,extractUniqueAA=extractUniqueAA,
                                                  getSharingStatistics=True,onlyProductive=onlyProductive)
        else:
            print 'loading existing TCR: datasetName=%s, minNshared=%s' % (datasetName, minNshared)
            try:
                filename = '%s/sharingAnalysis/sharingMatrix_%s_minNshared%s_RA_onlyProductiveTrue' % (datasetFolder, datasetName, minNshared)
                print filename
                df = pd.read_pickle(filename)
                print 'load completed'
            except:
                try:
                    filename = '%s/%s_sharingAnalysis/sharingMatrix_%s_minNshared%s_RA_onlyProductiveTrue' % (datasetFolder, SampleListName,datasetName, minNshared)
                    print filename
                    df = pd.read_pickle(filename)
                    print 'load completed'
                except:
                    print 'TCR file was not found'

        print 'checking if RA values are percentage values or real RA values:'
        currMinVal = df[df > 0].min().min()
        print 'the minimal value of the DF is %s' % currMinVal

        if currMinVal > 0.001:
            print 'the minimal value is %s, which is higher than 0.001. this means that RA values are percentage. dividing by 100\
    to get real RA values...' % currMinVal
            df = df.div(100)
        else:
            print 'RA values are real RA values, no division is required'

    else:
        print 'DFtype is not known'

    print 'Number of analyzed units (seqs/species/genus) is %s' % len(df.columns)
    print 'filename is %s' % filename
    print 'df shape is %s_%s' % (df.shape[0], df.shape[1])

    origDF = df.copy()
    
    # # (2) cap value to minVal
    
    # ## define minimal val:
    # minVal can be None,0, float, or 'dfMinVal' 
    print 'setting minVal'
    if minVal == 'dfMinVal':
        minVal2 = df[df > 0].min().min()
        print 'minVal is set to the existing df minimal value: %s' % minVal2
    elif minVal == 'dfMinVal2':
        minVal2 = (df[df > 0].min().min()) * 2
        print 'minVal is set to the existing df minimal value: %s' % minVal2
    else:
        minVal2 = minVal
        print 'minVal remains %s as defined by user' % minVal2
    
    # ## cap all values to minimal val:
    if minVal is not None:
        print 'capping to minVal=%s' % minVal2
        print 'df sum before capping is %s' % df.sum().sum()
        df = df.mask(df < minVal2)
        df = df.fillna(minVal2)
        df.head()
        print 'df sum after capping is %s' % df.sum().sum()
        name1 = 'capped%s' % minVal
    else:
        print 'no capping was applied'
        minVal2 = 0
        name1 = ''
        
    # # (3)Leave only seq/species appear in specific number of samples  (use '2' to leave only shared seqs/species
    # define minSharedT=minimal number of shared samples by each seq/species to leave in df, or-
    # percShared = leave only seqs/species shared by specific percent of samples


    if (minSharedT is not None) or (percShared is not None):
        print 'Number of analyzed units is %s' % len(df.columns)
        if minSharedT is not None:
            print 'leaving only units shared by %s samples or more' % minSharedT
            checkShared = df.apply(lambda x: True if len([i for i in x if i > minVal2]) >= minSharedT else False)
            name2 = 'minSharedT%s' % minSharedT
        elif percShared is not None:
            minSharedT = round(float(len(df)) * percShared / 100)
            print 'leaving only units shared by more than %s percent of the samples (%s samples)' % (percShared, minSharedT)
            checkShared = df.apply(lambda x: True if len([i for i in x if i > minVal2]) >= minSharedT else False)
            name2 = 'percShared%s' % percShared
        shared = checkShared[checkShared]
        df = df[shared.index]
        print 'number of analyzed units is now %s' % len(df.columns)
    else:
        print 'no units were removed due to low sharing level, number of units remains %s' % len(df.columns)
        name2 = ''
        
    # # (4) Check outliers and trim to the maximal allowed value (not for binary!)
    # # trim outliers
    if removeOutliers:
        print 'removing outliers, outlierSTD=5, trim=True...'
        df2 = df.copy()
        df = filter_outliers(df, outlierSTD=5, columnList=None, trim=True)
        name3 = 'OLtrimmed'
        
           # ## plot max and min values distribution to follow correct removal   
        print 'plotting values distribution before and after trimming...'
        fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 2))
        ax1 = axes[0]
        x1 = df2.min().tolist()
        y1 = df.min().tolist()
        ax1.hist([x1, y1], color=['blue', 'red'], bins=20, label=['before trimming', 'after trimming'], alpha=0.5)
        print 'now the second'
        ax1.set_yscale('log')
        ax1.legend()
        ax1.set_title('min RA value distribution before\nand after trimming')
        ax1.set_xlabel('RA values')
        ax1.set_ylabel('freq')
        ax1.xaxis.set_major_locator(plt.MaxNLocator(4))
    
        ax2 = axes[1]
        x2 = df2.max().tolist()
        y2 = df.max().tolist()
        ax2.hist([x2, y2], color=['blue', 'red'], bins=20, label=['before trimming', 'after trimming'], alpha=0.5)
        print 'now the second'
        ax2.set_yscale('log')
        ax2.legend()
        ax2.set_title('max RA value distribution before\nand after trimming')
        ax2.set_xlabel('RA values')
        ax2.set_ylabel('freq')
        ax2.xaxis.set_major_locator(plt.MaxNLocator(4))
        plt.show()
    else:
        name3 = ''
        
        
 
    
    
    # # (5) Binarize data (use minimal value as zero)
    if toBinary:
        print 'binarizing df...'
        dfBinary = df.mask(df <= minVal2)
        dfBinary = dfBinary.fillna(0)
        dfBinary = dfBinary.mask(dfBinary != 0)
        dfBinary = dfBinary.fillna(1)
        print 'existing values in df are now:'
        print pd.unique(dfBinary.values.ravel('K'))
    
    # # (6) Norm data (not for binary!)
    if normData:
        print 'normalizing data so each sample  sum will be 1...'
        normeddf = df.div(df.sum(axis=1), axis=0)
        print df.sum(axis=1).head()
        print normeddf.sum(axis=1).head()
    
    # #  (7)log transformation (not for binary!)
    if logTransform:
        if 0 in df.values:
            print 'values contain zeros, log transformation was not applied'
        else:
            print 'transforming df to log scale'
            logdf = df.apply(lambda x: np.log10(x))
            print 'df descriptives:'
            logdf.describe()
    
    
    
    # # (8) generate some summarizing parameters on the df and save
    # ## generate final filename
    finalFileName = '%s_%s_%s_%s' % (filename, name1, name2, name3)
    finalFileName = finalFileName.replace('.', '_')
    print 'final file name is:'
    print finalFileName    
    name = finalFileName.split('/')[-1]
    
    
    # ## (9) plot mean sample value distribution for each df:
    dfList = [origDF, df]
    dfListNames = ['origDF', name]
    try:
        dfList.append(dfBinary)
        dfListNames.append('dfBinary')
    except:
        dfList.append(None)
        dfListNames.append(None)

    try:
        dfList.append(normeddf)
        dfListNames.append('normeddf')
    except:
        dfList.append(None)
        dfListNames.append(None)

    try:
        dfList.append(logdf)
        dfListNames.append('logdf')
    except:
        dfList.append(None)
        dfListNames.append(None)
        
    print 'plotting mean sample distributions for all generated dfs...'
    dfListNew = [x for x in dfList if x is not None]
    dfListNamesNew = [x for x in dfListNames if x is not None]

    llen = len(dfListNew)
    fig, axes = plt.subplots(nrows=1, ncols=llen, figsize=(4 * llen, 2))
    fig.suptitle ('RA value distributions', fontweight='bold')
    for n, dfx in enumerate(dfListNew):
        print n, dfx.shape
        if len(dfx.columns.values)<100000:
            ax = axes[n]
            x = dfx.mean(axis=1).tolist()
            ax.hist(x)
            ax.set_yscale('log')
            ax.set_title(dfListNamesNew[n])
        #     ax.set_xlabel('RA values')
            if n == 0:
                ax.set_ylabel('freq')
            if n == 1:
                title = dfListNamesNew[n]
                firstpart, secondpart = title[:len(title) / 2], title[len(title) / 2:]
                ax.set_title('%s\n%s' % (firstpart, secondpart))
            ax.xaxis.set_major_locator(plt.MaxNLocator(3))
        else:
            print '%s is too big, plot was not generated' %dfListNamesNew[n]

    fig.subplots_adjust(top=0.7, wspace=0.4)

    plt.show()
    
    print 'plotting RA value distributions for all generated dfs IN A RANDOM SAMPLE'
    import random
    dfListNew = [x for x in dfList if x is not None]
    dfListNamesNew = [x for x in dfListNames if x is not None]

    llen = len(dfListNew)
    fig, axes = plt.subplots(nrows=1, ncols=llen, figsize=(4 * llen, 2))

    sample = random.sample(df.index, 1)

    fig.suptitle ('RA value distributions in sample %s' % sample[0], fontweight='bold')
    for n, dfx in enumerate(dfListNew):
        ax = axes[n]
        x = dfx.loc[sample, :].values.flatten().tolist()
        ax.hist(x, bins=20)
        ax.set_yscale('log')
        ax.set_title(dfListNamesNew[n])
    #     ax.set_xlabel('RA values')
        if n == 0:
            ax.set_ylabel('freq')
        if n == 1:
            title = dfListNamesNew[n]
            firstpart, secondpart = title[:len(title) / 2], title[len(title) / 2:]
            ax.set_title('%s\n%s' % (firstpart, secondpart))
        ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    fig.subplots_adjust(top=0.62, wspace=0.4)

    plt.show()
    
    
    print 'saving df to file...'
    df.to_pickle(finalFileName)

    try:
        dfBinary.to_pickle('%s_binary' % finalFileName)
        print 'saving binary df to file'
    except:
        pass

    try:
        normeddf.to_pickle('%s_normed' % finalFileName)
        print 'saving normed df to file'
    except:
        pass

    try:
        logdf.to_pickle('%s_log' % finalFileName)
        print 'saving log df to file'
    except:
        pass


    print 'ALL is done!!!'
    
    return df



#----------------------------------------------------------------------

'''
the following function calculates diversity and clonality features for MPA table

USAGE EXAMPLE:

MPAtableName='MPA_g_Filtered'
MBclonDivFeatures=calcMBclonDivFeaturesPerMPAtable(MPAtableName)



'''
def calcMBclonDivFeaturesPerMPAtable(MPAtableName,removeOutliers=True):
    
    #load relevant MPA table
    print 'loading relevant MPA table:'
    MPAtableFile='%s/MicrobiomeDataTables/MbFeatureTables/%s' %(MyPath,MPAtableName)
    MPAtable=pd.read_pickle(MPAtableFile)
    print MPAtable.iloc[:5,:5]
    
    # generate new dataframe with index identical to the MPAtable index:
    print 'generating MBclonDivFeatures...'
    MBclonDivFeatures = pd.DataFrame(index=MPAtable.index.copy())
    # calculate features and add to the df:
    MBclonDivFeatures['nTaxa']=MPAtable[MPAtable>0.0001].count(axis=1)
    MBclonDivFeatures['meanRA']=MPAtable.mean(axis=1)
    MBclonDivFeatures['meanRAno0']=MPAtable[MPAtable>0.0001].mean(axis=1)
    MBclonDivFeatures['meadianRAno0']=MPAtable[MPAtable>0.0001].median(axis=1)
    MBclonDivFeatures['stdRA']=MPAtable.std(axis=1)
    MBclonDivFeatures['stdRAno0']=MPAtable[MPAtable>0.0001].std(axis=1)
    MBclonDivFeatures['maxRA']=MPAtable.max(axis=1)
    MBclonDivFeatures['max/meanRA']=MBclonDivFeatures['maxRA']/MBclonDivFeatures['meanRA']
    MBclonDivFeatures['max/meanRAno0']=MBclonDivFeatures['maxRA']/MBclonDivFeatures['meanRAno0']
    MBclonDivFeatures['sumTop10']=MPAtable.apply(lambda x: np.sum(np.sort(x)[-10:]), axis=1)

    MBclonDivFeatures['shannon']=MPAtable.copy().multiply(1000).round(0).apply(lambda x: shannon(x.astype('int64'), base=2),axis=1)
    MBclonDivFeatures['simpson']=MPAtable.copy().multiply(1000).round(0).apply(lambda x: simpson(x.astype('int64')),axis=1)
    MBclonDivFeatures['berger_parker_d']=MPAtable.copy().multiply(1000).round(0).apply(lambda x: berger_parker_d(x.astype('int64')),axis=1)

    print 'MBclonDivFeatures table shape is %s_%s' %(MBclonDivFeatures.shape[0],MBclonDivFeatures.shape[1])
    
    #removing outliers:
    if removeOutliers:
        print 'removing outliers, outlierSTD=10, trim=True...'
        MBclonDivFeatures = filter_outliers(MBclonDivFeatures, outlierSTD=10, columnList=None, trim=True)
    
    
    
    #calc feature correlation
    print 'calculating feature correlation...'
    f, ax = plt.subplots(figsize=(10, 8))
    corr = MBclonDivFeatures.corr()
    sns.heatmap(corr, mask=np.zeros_like(corr, dtype=np.bool), cmap=sns.diverging_palette(220, 10, as_cmap=True),
                square=True, ax=ax)
    ax.set_title('%s_clonDivFeatures_CorrelationMap' %MPAtableName,fontsize=18)
    plt.show()
    
    #save table and figure:
    file1= '%s/MicrobiomeDataTables/MbFeatureTables/%s_MBclonDivFeatures_removeOL%s' %(MyPath,MPAtableName,removeOutliers)
    MBclonDivFeatures.to_pickle(file1)

    corrFigFile='%s/MicrobiomeDataTables/MbFeatureTables/featuresFigs/%s_MBclonDivFeaturesCorr_removeOL%s_Fig.png' %(MyPath,MPAtableName,
                                                                                                                    removeOutliers)
    f.savefig(corrFigFile,dpi=200)
    
    return MBclonDivFeatures
    
    

