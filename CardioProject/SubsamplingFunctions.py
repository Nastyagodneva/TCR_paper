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
from scipy.spatial.distance import braycurtis, pdist, euclidean
# from tunneltoamazondb import getengine
from pandas import concat, read_csv, Series, read_sql
from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import *
from ShaniBA.PhenotypicData.PhenotypeGenerationFunctions import *

from sklearn import metrics, preprocessing
from sklearn.svm import LinearSVC
from sklearn.feature_selection import SelectFromModel, SelectKBest, chi2, mutual_info_classif, f_classif
from sklearn.naive_bayes import GaussianNB


from ShaniBA.TCR_microbiome_interactions.TCR_microbiome_interactions_functions import *
from ShaniBA.TCR_feature_generation.TCR_feature_generation_functions import *

MyPath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF'

#------------------------------------------------------------------------------
'''
downsampling function from Yael
'''

def downsample_vector_v2(counts, dsn):
    """
    downsample given vector to exactly dsn counts
    """

    assert counts.sum() >= dsn  # debugging statement
    if counts.sum() == dsn:
        return counts

    useful = counts.nonzero()[0]
    l = np.array([0] * int(counts.sum()))
    cumsum = np.cumsum(counts[useful])
    cumsum = np.insert(cumsum, 0, 0).astype(int)
    for jind in range(len(useful)):
        l[cumsum[jind]:cumsum[jind + 1]] = useful[jind]

    tmp = np.random.choice(l, int(dsn), replace=False)
    tmpu, tmpc = np.unique(tmp, return_counts=True)
    counts_ds = np.zeros(counts.shape)
    counts_ds[tmpu] = tmpc

    return counts_ds


#--------------------------------------------------
'''
the following function subsamples a df to take a specified number of templates

input:
sample_df
nTemplates=number of templates to sample

usage example:
fullSamplesFolder='%s/TCR_real_data/SamplesForAnalysis_corrected' %MyPath
nTemplates=18000
nSampled=0
subsampledSamplesFolder='%s/TCR_real_data/SubSampled%s data' %(MyPath,nTemplates)
if not isdir(subsampledSamplesFolder):
        makedirs(subsampledSamplesFolder) 
FullFiles = [f for f in listdir(fullSamplesFolder) if isfile(join(fullSamplesFolder, f))]
FullFiles=[f.strip('.tsv') for f in FullFiles]
print 'number of samples in folder is %s' %len(FullFiles)

# #there is a mistake in the column names in the batch released on 25/2/18:
# #get right column names:
# f='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/SamplesForAnalysis/BD438.tsv'
# BD438=pd.read_table(f)
# right_column_names=BD438.columns.values
newColumnList=None


for n,f in enumerate(FullFiles):
#     if n>408:
        print n,f

        #correct column names:
        sample_df=pd.read_table('%s/%s.tsv' %(fullSamplesFolder,f))
        if newColumnList is not None: #make sure column names are correct:
            sample_df=sample_df.iloc[:,:44]
            sample_df.columns=newColumnList
        sample_df=sample_df.rename(columns={'count (templates/reads)':'count (templates)'})


        nTemplatesSample=sample_df['count (templates)'].sum()
        if nTemplatesSample>=nTemplates:
            subsampled15full=subsampling(sample_df,nTemplates)
            file1='%s/%s_%s.tsv' %(subsampledSamplesFolder,f,nTemplates)
            subsampled15full.to_csv(file1,sep='\t')
            nSampled=nSampled+1
        else:
            print 'nTemplates in sample is %s and therefore this sample is not subsampled' %nTemplatesSample

with open('%s/nSampled' %subsampledSamplesFolder,'wb') as fp:
    pickle.dump(nSampled,fp)

'''
def subsampling(sample_df, nTemplates):
    
    # generating a popped nucleotide sequenes:
    print 'generating a popped nucleotide sequenes:...'
    inds = sample_df.index
    counts = sample_df['count (templates)']
    l = np.array([0] * int(counts.sum()))
    print 'n Unique sequences in sample is %s' % len(inds)
    print 'nTemplates in sample is %s' % sample_df['count (templates)'].sum()
    print 'l list length is %s and should be equal to nTemplates' % len(l)
    
    cumsum = list(np.cumsum(counts))
    print 'cumsum list length is %s and should be equal to n Unique' % len(cumsum)
    cumsum = np.insert(cumsum, 0, 0).astype(int)
    
    for jind in range(len(inds)):  # loop over all numbers from zero to the length of the vector
        l[cumsum[jind]:cumsum[jind + 1]] = int(inds[jind])  # replace zeros in the l vector with the relevant index
                                                     # repeat this value n times according to the number of templates per this index
    indCountdf = pd.DataFrame(l)  # generate a df contains all sequences' indices and their counts
    # generate a popped nucleotide df:
    popped = pd.merge(indCountdf, pd.DataFrame(sample_df['nucleotide']), how='left', left_on=0, right_index=True)
    popped = popped.drop(0, axis=1)
    print 'length of popped nucleotide dataframe is %s and should be equal to nTemplates' % len(popped)
    
    # subsample:
    print 'subsampling...'
    subsampled = popped.sample(n=nTemplates, replace=False)
    print 'subsampled df length is %s and should be equal to %s' % (len(subsampled), nTemplates)
    
    
    subsampledGrouped = pd.DataFrame(subsampled['nucleotide'].value_counts())
    subsampledGrouped.rename(columns={'nucleotide':'count (templates)'})
    subsampledFull = pd.merge(subsampledGrouped, sample_df, how='left', left_index=True, right_on='nucleotide')
    subsampledFull = subsampledFull.rename(columns={'count (templates)':'oldTemplateCount', 'nucleotide_x':'count (templates)',
                                                   'frequencyCount (%)':'frequencyCount (%) old'})

    subsampledFull = subsampledFull.drop('nucleotide_y', axis=1)
    subsampledFull['frequencyCount (%)'] = subsampledFull['count (templates)'] * 100 / subsampledFull['count (templates)'].sum()
    subsampledFull.head()
    print 'templates sum in subsampled file is %s and should be equal to %s' % (subsampledFull['count (templates)'].sum(), nTemplates)
    
    return subsampledFull

#--------------------------------------------------------------
'''
the following function is a mega functions that subsample, extract features, summarize features and compare to original datasets:

input:
nTemplates - number of templates to sample #always change!
repeat - repeat number for this nTemplate sampling - to enable multisampling of the same number of templates #always change!
datasetName= - a string. for example 'PNP515' #change if necessary
fullSamplesFolder - the folder of the original samples, that should be sub-sampled #change if necessary
data_folder_full='TCR_real_data' - the high level folder
TakeSameSamples = True/False - use true when comparing subsampled dataset to its original set

usage example:
nTemplates=12500 
repeat=2 
datasetName='PNP515' 
fullSamplesFolder='%s/TCR_real_data/SamplesForAnalysis_corrected' %MyPath 
data_folder_full='TCR_real_data' #change if necessary
TakeSameSamples=True
'''


def subsampling_and_featureExtraction(fullSamplesFolder, nTemplates, repeat, datasetName, data_folder_full, TakeSameSamples):
    
    data_folder = '%s/%s_SubSampled%sdata_rep%s' % (data_folder_full, datasetName, nTemplates, repeat)
    if not isdir(data_folder):
            makedirs(data_folder) 
    # (1) SUBSAMPLING:

    print 'step 1: subsampling (long)'

    subsampledSamplesFolder = '%s/%s/SamplesForAnalysis_corrected' % (MyPath, data_folder)
    if not isdir(subsampledSamplesFolder):
            makedirs(subsampledSamplesFolder) 
    FullFiles = [f for f in listdir(fullSamplesFolder) if isfile(join(fullSamplesFolder, f))]
    FullFiles = [f.strip('.tsv') for f in FullFiles]
    print 'number of samples in folder is %s' % len(FullFiles)

    for n, f in enumerate(FullFiles):
    #     if n>408:
            print n, f
            sample_df = pd.read_table('%s/%s.tsv' % (fullSamplesFolder, f))
            nTemplatesSample = sample_df['count (templates)'].sum()
            if nTemplatesSample >= nTemplates:
                subsampled15full = subsampling(sample_df, nTemplates)
                file1 = '%s/%s_%s.tsv' % (subsampledSamplesFolder, f, nTemplates)
                subsampled15full.to_csv(file1, sep='\t')
            else:
                print 'nTemplates in sample is %s and therefore this sample is not subsampled' % nTemplatesSample

    # (2) feature extraction:
    print 'step 2: feature extraction (long)'

    filenames = [f for f in listdir(subsampledSamplesFolder) if isfile(join(subsampledSamplesFolder, f))]
    filenames = [f.strip('.tsv') for f in filenames]
    filenames = [f.strip('.xlsx') for f in filenames]
    print 'number of samples for feature extraction is %s' % len(filenames)

    newColumnList = None
    for n, sample_name in enumerate(filenames): 
    #     if n>91: 
            print  n, sample_name
            if 'nSampled' not in sample_name:
                gen_descriptive_stats(sample_name, data_folder, newColumnList)
                gen_geneUsageCount(sample_name, data_folder, newColumnList)


    # (3) generate feature summary DFs:

    print 'generating seperate feature data dfs: (long)'

    seqTypeList = ['Total', 'Prod', 'nonProd']
    newDatasetName = '%s_ss%s_rep%s' % (datasetName, nTemplates, repeat) 
    FeatureFolder = '%s/%s/descriptiveStatsSamplesForAnalysis/Total' % (MyPath, data_folder)
    FeatureGroups = listdir(FeatureFolder)
    for seqType in seqTypeList:
        print seqType
        gen_featureSummaryDF_forSeqType(seqType, data_folder, newDatasetName, FeatureGroups)

    print 'generating merged feature data df: (short)'
    gen_merged_feature_df_for_dataset(data_folder, newDatasetName)


    # (4) compare features between subsampled dataset and original one:
    print 'comparing features between ss and original datasets... (short)'

    data_folder1 = data_folder
    data_folder2 = data_folder_full
    datasetName1 = '%s_ss%s_rep%s' % (datasetName, nTemplates, repeat)
    datasetName2 = datasetName



    compare_features_between_datasets(data_folder1, datasetName1, data_folder2, datasetName2, TakeSameSamples)
    plot_gene_usage_comparison(data_folder1, datasetName1, data_folder2, datasetName2, TakeSameSamples)

    print 'done'
    

#-------------------------------------------------------------------------------------------------------------------------

########################################################################################################################
### the following functions serve to generate PNP and diseased cohorts which are matched in terms of subsampling    ####
### and relevant phenotypes:                                                                                        ####
### function 1: subsampling, generation of common TCRdf to use with mantel and permanova test, and conduct the      ####
### mantel/permanova test to check which phenotypes are correlated with TCRdf structure                             ####
### function 2: generate matched cohorts based on important phenotypes as evident from function 1, and compare      ####
### phenotypes in the resulting cohorts to validate matching.                                                       ####
### function 2:                                                                                                     ####
### generate matched cohorts and compare phenotype                                                                  ####
### function 3:                                                                                                     ####
### generate a folder with all samples, compare features, seperate by PCA and find unique sequences                 ####
########################################################################################################################

 
'''
 Function 1:
 ubsampling, generation of common TCRdf to use with mantel and permanova test, and conduct the      ####
### mantel/permanova test to check which phenotypes are correlated with TCRdf structure
 
 
 usage example:
ss=12500 (int. how many sequences to subsample)
repeat=1 (int. ordinal of the repeat of the same process)
ssPNP=False #True/False. whether or not to subsample (use True when subsampling haven't occured before)
ssCardio=False #True/False. whether or not to subsample (use True when subsampling haven't occured before)
genTCRdfPNP=False #True/False. 
genTCRdfCardio=False #True/False. 


testPhenotypeAffectsOnsubsampledCohorts(ss,repeat,ssPNP,ssCardio,genTCRdfPNP,genTCRdfCardio)
 
 
 '''  
    
def testPhenotypeAffectsOnsubsampledCohorts(ss, repeat, ssPNP, ssCardio, genTCRdfPNP, genTCRdfCardio):

    
    # load sample lists:
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP530', 'rb') as fp:
        PNP530 = pickle.load(fp)
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/Cardio126', 'rb') as fp:
        Cardio126 = pickle.load(fp)
    
    
    
    # (1) subsample PNP and Cardio:
    if ssPNP: 
        print '********subsampling PNP530 cohort with %s sequences per sample:***********'
        nTemplates = ss
        repeat = repeat 
        datasetName = 'PNP530' 
        fullSamplesFolder = '%s/TCR_real_data/SamplesForAnalysis_corrected' % MyPath 
        data_folder_full = 'TCR_real_data'  # change if necessary
        TakeSameSamples = True
        subsampling_and_featureExtraction(fullSamplesFolder, nTemplates, repeat, datasetName, data_folder_full,
                                          TakeSameSamples=TakeSameSamples)
    if ssCardio: 
        print '*********subsampling Cardio cohort with %s sequences per sample:***************'
        nTemplates = ss
        repeat = repeat 
        datasetName = 'Cardio126' 
        fullSamplesFolder = '%s/TCR_real_data/CardioSamples/SamplesForAnalysis_corrected' % MyPath 
        data_folder_full = 'TCR_real_data/CardioSamples'  # change if necessary
        TakeSameSamples = True
        subsampling_and_featureExtraction(fullSamplesFolder, nTemplates, repeat, datasetName, data_folder_full,
                                          TakeSameSamples=TakeSameSamples)
        
    # compare TCRfeatures between ss cohorts:
    if ss is None:
        partialDatasetFolderPNP = 'TCR_real_data' 
        datasetNamePNP = 'PNP530'
        partialDatasetFolderCardio = 'TCR_real_data/CardioSamples' 
        datasetNameCardio = 'Cardio126'
    else:
        partialDatasetFolderPNP = 'TCR_real_data/PNP530_SubSampled%sdata_rep%s' % (ss, repeat)
        datasetNamePNP = 'PNP530_ss%s_rep%s' % (ss, repeat)
        partialDatasetFolderCardio = 'TCR_real_data/CardioSamples/Cardio126_SubSampled%sdata_rep%s' % (ss, repeat) 
        datasetNameCardio = 'Cardio126_ss%s_rep%s' % (ss, repeat)
        
    datasetFolderPNP = '%s/%s' % (MyPath, partialDatasetFolderPNP)
    datasetFolderCardio = '%s/%s' % (MyPath, partialDatasetFolderCardio)
    
    
    print '************comparing TCRfeatures between cohorts (not matched!):*************'
    data_folder1 = partialDatasetFolderPNP
    data_folder2 = partialDatasetFolderCardio
    datasetName1 = datasetNamePNP
    datasetName2 = datasetNameCardio
    TakeSameSamples = False
    filteringList1 = None
    filteringList2 = None
    filteringList1Name = None
    filteringList2Name = None
    plotType = 'bar'
    
    compare_features_between_datasets(data_folder1, datasetName1, data_folder2, datasetName2, TakeSameSamples,
                                  filteringList1, filteringList2, filteringList1Name, filteringList2Name)
    plot_gene_usage_comparison(data_folder1, datasetName1, data_folder2, datasetName2, plotType, TakeSameSamples,
                                      filteringList1, filteringList2, filteringList1Name, filteringList2Name)
    
    
    # (2) produce TCRdfs to use with permanova/mantel:
    if genTCRdfPNP:
        print '*****generating 10perc TCRdfs******'
        DFtype = 'TCR'
        genDF = True  # False=generate new df
        toBinary = True
        mbLevel = 'g'
        useShortName = True
        datasetFolder = '%s/%s' % (MyPath, data_folder1)
        datasetName = datasetName1
        minVal = None  # minVal can be None,0, float, or 'dfMinVal' or dfMinVal2:
        minSharedT = None  # minimal number of samples shared by seq/species required to leave sample in the database (int or None)
        percShared = 10  # minimal number of samples shared by seq/species required to leave sample in the database (int [ percentage]
                        # or None)
        removeOutliers = True
        normData = True
        logTransform = True
        extractUniqueAA = True  # use True when this is the first time to analyze this dataset, otherwise, use False
        minNshared = 2
        onlyProductive = True
        mbDataFolder = 'AllSeqProjects'
        PNP530Cardio126list = PNP530 + Cardio126
        SampleList = PNP530
        SampleListName = 'PNP530'
        libPrepMethod = None
        filterGenotek = True
        filterMinimalReads = 9000000
        filterlibPrepMethod = libPrepMethod
        groupFunction = 'noOutlierMean'
        nSTD = 5
        nMinSamples = 3
        ignoreNotSameYear = True
        removeSamePerson = False
    
        df = genTCRorMBdfWithManipulations(DFtype, genDF, toBinary, removeOutliers, normData, logTransform,
                                         minVal, minSharedT, percShared,
                                         mbLevel, useShortName, datasetFolder, datasetName, extractUniqueAA,
                                        minNshared, onlyProductive, mbDataFolder, SampleList,
                                          SampleListName, filterlibPrepMethod, filterGenotek,
                                          groupFunction, nSTD, nMinSamples, ignoreNotSameYear, removeSamePerson)
    if genTCRdfCardio:
        print '***********generating 10perc TCRdfs for Cardio************'
        datasetFolder = '%s/%s' % (MyPath, data_folder2)
        datasetName = datasetName2
        SampleList = Cardio126
        SampleListName = 'Cardio126'
    
        df = genTCRorMBdfWithManipulations(DFtype, genDF, toBinary, removeOutliers, normData, logTransform,
                                         minVal, minSharedT, percShared,
                                         mbLevel, useShortName, datasetFolder, datasetName, extractUniqueAA,
                                        minNshared, onlyProductive, mbDataFolder, SampleList,
                                          SampleListName, filterlibPrepMethod, filterGenotek,
                                          groupFunction, nSTD, nMinSamples, ignoreNotSameYear, removeSamePerson)
        
    # (3) run permanova and mantel to test which phenotype are associated with TCRdf structure
    # load PNP phenotypeDF:
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/newPhenotypesPNPAllInfo_withDummies.xlsx'
    PNP530_phen_new_dummies = pd.read_excel(f1).set_index('BD')
    
    # load Cardio phenotypeDF:
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/Cardio126phenAllInfo_withDummies.xlsx'
    Cardio126_phen_new_dummies = pd.read_excel(f1).set_index('BD')
    
    numericalphenotypesPNP = ['Age', 'BMI', 'WBC', 'eGFR by CKD-EPI', 'HbA1C', 'Glucose', 'HDL', 'Hemoglobin', 'AST',
                            'Total Cholesterol']
    numericalphenotypesCardio = ['Age', 'BMI', 'eGFR by CKD-EPI', 'WBC', 'LDL', 'HDL', 'Triglycerides', 'CRP', 'HbA1C', 'PLT',
                               'Hemoglobin', 'Initial CPK', 'Maximal CPK', 'LDH', 'AST', 'Glucose', 'Initial Troponin',
                               'Maximal Troponin', 'Total Cholesterol', 'GRACE Score', 'Admission Systolic BP',
                               'Admission Diastolic BP', 'Admission Pulse Rate', 'Admission Saturation Value',
                               'Admission Temperature']
    categoricalphenotypesPNP = ['Gender_Male', 'Smoking Status_Past', 'Smoking Status_Yes']
    categoricalphenotypesCardio = ['Gender_Male', 'Smoking Status_Past', 'Smoking Status_Yes', 'LVEFmapped', 'PreviousPCImapped_1',
            'PreviousPCImapped_2+', 'Hypertension', 'Dyslipidemia',
           'Microvascular Complications', 'Previous CABG', 'Known CAD', 'PCI_binary', 'Religion', 'PVD',
           'Glucose Disorder_DM2', 'Glucose Disorder_No', 'Glucose Disorder_PreDM', 'Admission Diagnosis_NSTEMI',
           'Admission Diagnosis_STEMI',
           'Admission Statins', 'Chief Complaint_anginal pain', 'Chief Complaint_atypical pain', 'Chief Complaint_dyspnea']
  
    
    
    numMetric = 'euclidean'
    catMetric = 'jaccard'
    
    # generate distMats:
    featureDFfile = '%s/TCR_real_data/PNP530Cardio126Combined/featureSummaryDFs/PNP530Cardio126_filteredByPNP530Cardio126\
_allFeatures_noCorrelated_noConsts_filledna' % MyPath
    TCRfeatureDF2 = pd.read_pickle(featureDFfile)
    
    PNP530featureDF = TCRfeatureDF2.loc[PNP530, :]
    Cardio126featureDF = TCRfeatureDF2.loc[Cardio126, :]
    datasetList = ['PNP530', 'Cardio126']
    for dataset in datasetList: 
        if 'PNP' in dataset:  # DEFINE RELEVANT DATASET FOLDER, NAME AND PHENOTYPE DF FILE
            print 'generating distmats for PNP'
            datasetFolder = '%s/%s' % (MyPath, data_folder1)
            datasetName = datasetName1
            phenotypeDF = PNP530_phen_new_dummies
            categoricalphenotypes = categoricalphenotypesPNP
            numericalphenotypes = numericalphenotypesPNP
            SampleListName = 'PNP530'
            featureDF = PNP530featureDF
        elif 'Cardio' in dataset:  # DEFINE RELEVANT DATASET FOLDER, NAME AND PHENOTYPE DF FILE
            print 'generating distmats for Cardio'
            datasetFolder = '%s/%s' % (MyPath, data_folder2)
            datasetName = datasetName2
            phenotypeDF = Cardio126_phen_new_dummies
            categoricalphenotypes = categoricalphenotypesCardio
            numericalphenotypes = numericalphenotypesCardio
            SampleListName = 'Cardio126'
            featureDF = Cardio126featureDF
            
        print '****************generating phenotype distMats*******************'
        distMatFolder = '%s/phenDistMats' % datasetFolder
        if not isdir(distMatFolder):
            makedirs(distMatFolder)
        for n, phenotype in enumerate(numericalphenotypes):
            metric = 'euclidean'
            print 'generating distMat for phenotype %s' % phenotype
    
            df = pd.DataFrame(phenotypeDF[phenotype])
            df_condensed_org, distMat_square = genDistMat(df, metric)
            
            if ss is None:
                distMatFileSquare = '%s/%s_%s_%s_distMat' % (distMatFolder, dataset, phenotype, metric)
                distMatFileCondensed = '%s/%s_%s_%s_distMat_CONDENSED' % (distMatFolder, dataset, phenotype, metric)
            else:
                distMatFileSquare = '%s/%sss%s_rep%s_%s_%s_distMat' % (distMatFolder, dataset, ss, repeat, phenotype, metric)
                distMatFileCondensed = '%s/%sss%s_rep%s_%s_%s_distMat_CONDENSED' % (distMatFolder, dataset, ss, repeat, phenotype, metric)                 
            distMat_square.to_pickle(distMatFileSquare)
            df_condensed_org.to_pickle(distMatFileCondensed)
            
        print '****************generating feature distMats*******************'
        distMatFolder = '%s/featureDistMats' % datasetFolder
        if not isdir(distMatFolder):
            makedirs(distMatFolder)
        featureList = ['top10clonal_nt_T', 'frequencyCount (%)_max_T', 'frequencyCount (%)_mean_T', 'shannon_nt_T', 'cdr3Length_max_T',
                       'n1Insertion_mean_T', 'n2Insertion_mean_T']
        # 'normSeqNums_per2000_NT_T',
        for n, feature in enumerate(featureList):
            metric = 'euclidean'
            print 'generating distMat for feature %s' % feature
    
            df = pd.DataFrame(featureDF[feature])
            df = df[df[feature].notnull()]
            df_condensed_org, distMat_square = genDistMat(df, metric)
            
            if ss is None:
                distMatFileSquare = '%s/%s_%s_%s_distMat' % (distMatFolder, dataset, feature, metric)
                distMatFileCondensed = '%s/%s_%s_%s_distMat_CONDENSED' % (distMatFolder, dataset, feature, metric)
            else:
                distMatFileSquare = '%s/%sss%s_rep%s_%s_%s_distMat' % (distMatFolder, dataset, ss, repeat, feature, metric)
                distMatFileCondensed = '%s/%sss%s_rep%s_%s_%s_distMat_CONDENSED' % (distMatFolder, dataset, ss, repeat, feature, metric)                 
            distMat_square.to_pickle(distMatFileSquare)
            df_condensed_org.to_pickle(distMatFileCondensed)
        
        
        # ## generating TCRdf for TCRdf
        if ss is None:
            if dataset == 'PNP':
                TCRdfName = 'sharingMatrix_%s_minNshared5_RA_onlyProductiveTrue__percShared10_OLtrimmed' % dataset
            else:
                TCRdfName = 'sharingMatrix_%s_minNshared2_RA_onlyProductiveTrue__percShared10_OLtrimmed' % dataset                
        else:  
            TCRdfName = 'sharingMatrix_%s_ss%s_rep%s_minNshared2_RA_onlyProductiveTrue__percShared10_OLtrimmed' % (dataset, ss, repeat)
        
        
        print '******generating distMat for TCRdfs:**********'       
        distMatFolder = '%s/distanceMatrices' % datasetFolder
        if not isdir(distMatFolder):
            makedirs(distMatFolder)
        TCRdfnewName = TCRdfName.replace('sharingMatrix_', 'TCRdf')
        TCRdfnewName = TCRdfnewName.replace(datasetName, '')
        TCRdfnewName = TCRdfnewName.replace('_', '')
        print (TCRdfnewName)
        print 'generating distance matrix for RA data:'
        print 'sampleListName=%s' % SampleListName
        try:
            file1 = '%s/%s_sharingAnalysis/%s' % (datasetFolder, SampleListName, TCRdfName)
            print 'file name=%s' % file1
            TCRdf = pd.read_pickle(file1)
        except:
            try:
                file1 = '%s/sharingAnalysis/%s' % (datasetFolder, TCRdfName)
                print 'still trying...'
                print 'file name=%s' % file1
                TCRdf = pd.read_pickle(file1)
            except:
                print 'couldnt load TCRdf'
                
        metric1 = 'braycurtis'
        print 'TCRdf sum is %s' % TCRdf.sum().sum()
        distMatFileSquare = '%s/%s_%s_distMat' % (distMatFolder, TCRdfnewName, metric1)
        distMatFileCondensed = '%s/%s_%s_distMat_CONDENSED' % (distMatFolder, TCRdfnewName, metric1)
        df_condensed_org, distMat_square = genDistMat(TCRdf, metric1)
        print 'distMat_square sum is %s' % distMat_square.sum().sum()
        distMat_square.to_pickle(distMatFileSquare)
        df_condensed_org.to_pickle(distMatFileCondensed)  
    
        print 'generating distance matrix for binary data:'
        TCRdfName_binary = TCRdfName + '_binary'
        TCRdfnewName_binary = TCRdfnewName + '_binary'
        
        try:
            file2 = '%s/%s_sharingAnalysis/%s' % (datasetFolder, SampleListName, TCRdfName_binary)
            print 'file name=%s' % file2
            TCRdf_binary = pd.read_pickle(file2)
        except:
            try:
#                 print 'still trying...'
#                 print 'sampleListName=%s' %SampleListName
                file2 = '%s/sharingAnalysis/%s' % (datasetFolder, TCRdfName_binary)
                print 'still trying...'
                print 'file name=%s' % file2
                TCRdf_binary = pd.read_pickle(file2)
            except:
                print 'couldnt load TCRdf'
        
        metric2 = 'jaccard'
        distMatFileSquare_binary = '%s/%s_%s_distMat' % (distMatFolder, TCRdfnewName_binary, metric2)
        distMatFileCondensed_binary = '%s/%s_%s_distMat_CONDENSED' % (distMatFolder, TCRdfnewName_binary, metric2)
        df_condensed_org_binary, distMat_square_binary = genDistMat(TCRdf_binary, metric2)
        distMat_square_binary.to_pickle(distMatFileSquare_binary)
        df_condensed_org_binary.to_pickle(distMatFileCondensed_binary) 
    
        print '******done with distance matrices******'
    
        permanovaDFfolder = '%s/TCR_phenotype_relations/permanovaDFs' % datasetFolder
        if not isdir(permanovaDFfolder):
            makedirs(permanovaDFfolder)
                                                                                                            
        TCRdistMatFolder = '%s/distanceMatrices' % datasetFolder
#         TCRdistMatNameList=[f for f in listdir(TCRdistMatFolder) if isfile(join(TCRdistMatFolder,f))]
#         TCRdistMatNameList=[f for f in TCRdistMatNameList if '_CONDENSED' not in f]
#         print 'number of distMats to analyze in folder is %s' %len(TCRdistMatNameList)
        TCRdistMatNameList = ['TCRdfminNshared2RAonlyProductiveTruepercShared10OLtrimmed_binary_jaccard_distMat',
                            'TCRdfminNshared2RAonlyProductiveTruepercShared10OLtrimmed_braycurtis_distMat']
        for phenotype in categoricalphenotypes:
            print phenotype
            for TCRdistMatName in TCRdistMatNameList:
                nPerm = 99
                removeSameUser = True
                print TCRdistMatName
                TCRdistMat = pd.read_pickle(join(TCRdistMatFolder, TCRdistMatName))
    
                print 'TCRdistMat head:'
                print TCRdistMat.iloc[:5, :5]
                newName = TCRdistMatName.replace('minNshared2RAonlyProductiveTrue', '')
                phenotypeDF = phenotypeDF[phenotypeDF[phenotype].notnull()]
                print 'phenotype head:'
                print phenotypeDF[phenotype].head()
                print 'conducting permanova'
                permAnovaResDF = permanova_withWrapper(TCRdistMat, newName, phenotypeDF, phenotype,
                                             nPerm, removeSameUser, permanovaDFfolder)
                try:
                    if permAnovaResDF.loc[0, 'p'] < 0.05:
                        print 'repeating permanova for phenotype %s with 9999 permutations' % phenotype
                        nPerm = 9999
                        permAnovaResDF = permanova_withWrapper(TCRdistMat, newName, phenotypeDF, phenotype,
                                                 nPerm, removeSameUser, permanovaDFfolder)
                except:
                    pass
                
            for feature in featureList:
                nPerm = 99
                removeSameUser = True
                print feature
                metric = 'euclidean'
                distMatFolder = '%s/featureDistMats' % datasetFolder
                distMatFileSquare = '%s/%s_%s_%s_distMat' % (distMatFolder, dataset, feature, metric)
                featuredistMat = pd.read_pickle(distMatFileSquare)
    
                print 'featuredistMat head:'
                print featuredistMat.iloc[:5, :5]
                print 'conducting permanova'
                permAnovaResDF = permanova_withWrapper(featuredistMat, feature, featureDF, feature,
                                             nPerm, removeSameUser, permanovaDFfolder)
                try:
                    if permAnovaResDF.loc[0, 'p'] < 0.05:
                        print 'repeating permanova for feature %s with 9999 permutations' % phenotype
                        nPerm = 9999
                        permAnovaResDF = permanova_withWrapper(featuredistMat, feature, featureDF, feature,
                                             nPerm, removeSameUser, permanovaDFfolder)
                except:
                    pass
                
    
        print '****get PERMANOVA RESULTS FOR %s*****' % dataset                                                                                            
        dfs_folder = permanovaDFfolder
        PermanovaResults1 = concat_summarizing_dfs(dfs_folder)
        nTests = len(PermanovaResults1)
        PermanovaResults1 = PermanovaResults1.sort_values(by='p')
        PermanovaResults1 = add_corrected_pValues(PermanovaResults1, pValueColumn='p', nTests=nTests, FDR=0.1)
        PermanovaResults1['cohort'] = dataset
        print PermanovaResults1
        permanovaResultFile = '%s/TCR_phenotype_relations/phenotype_TCRdf_permanovaResults.xlsx' % datasetFolder
        PermanovaResults1.to_excel(permanovaResultFile)   
                                                                                                            
        print '******performing mantel test now, still with dataset %s******' % dataset
        mantelDFfolder = '%s/TCR_phenotype_relations/mantelDFs' % datasetFolder
        if not isdir(mantelDFfolder):
            makedirs(mantelDFfolder)
        for phenotype in numericalphenotypes:
            nPerm = 99
            method = 'spearman'
            alternative = 'greater'
            removeSameUser = True
                                                                                                            
            phenotypedistMatFolder = '%s/phenDistMats' % datasetFolder
            if ss is None:
                phenotypedistMatName = '%s_%s_euclidean_distMat' % (dataset, phenotype)
            else:
                phenotypedistMatName = '%sss%s_%s_rep%s_euclidean_distMat' % (dataset, ss, repeat, phenotype)
            phenotype_distMat = pd.read_pickle('%s/%s' % (phenotypedistMatFolder, phenotypedistMatName))
    
            TCRdistMatFolder = '%s/distanceMatrices' % datasetFolder
            TCRdistMatNameList = [f for f in listdir(TCRdistMatFolder) if isfile(join(TCRdistMatFolder, f))]
            TCRdistMatNameList = [f for f in TCRdistMatNameList if '_CONDENSED' not in f]
    
            print 'number of TCR distMats to analyze in folder is %s' % len(TCRdistMatNameList)
            for TCRdistMatName in TCRdistMatNameList:
    
                print 'conducting mantel test for %s in %s' % (phenotype, dataset) 
                TCR_distMat = pd.read_pickle(join(TCRdistMatFolder, TCRdistMatName))
                newName = TCRdistMatName.replace('minNshared2RAonlyProductiveTrue', '')
                feature_distMat = TCR_distMat
                feature_name = newName
                phenotype_name = phenotypedistMatName
                mantelResDF = mantel_withWrapper(feature_distMat, feature_name, phenotype_distMat, phenotype_name, mantelDFfolder,
                                       nPerm, method, alternative, removeSameUser, i=None, j=None)                                                                                            
                try:
                    if mantelResDF.loc[0, 'p'] < 0.05:
                        print 'repeating mantel for phenotype %s with 9999 permutations' % phenotype
                        nPerm = 9999
                        mantel_withWrapper(feature_distMat, feature_name, phenotype_distMat, phenotype_name, mantelDFfolder,
                                           nPerm, method, alternative, removeSameUser, i=None, j=None)
                except:
                    pass
    
        print 'getting results for mantel tests:'    
        mantelResults1 = concat_summarizing_dfs(mantelDFfolder)
        nTests = len(mantelResults1)
        mantelResults1 = mantelResults1.sort_values(by='p')
        mantelResults1 = add_corrected_pValues(mantelResults1, pValueColumn='p', nTests=nTests, FDR=0.1)
        f1 = '%s/TCR_phenotype_relations/phenotype_TCRdf_mantelResults.xlsx' % datasetFolder
        mantelResults1.to_excel(f1)
        print mantelResults1


print 'end of function!!!' 

#---------------------------------------------------------------------------------------------------
'''
function 2:
generate matched cohorts and compare phenotype 
the function 'findClosestSample' is an helper function to the main function 'gen_matched_cohorts_by_phenotypes' 

input to main function:
*cohort1phenotypeDF - phenotypeDF of the cohort to be matched.
*cohort1name - string
*cohort2phenotypeDF - phenotypeDF of the cohort from which samples are matched to cohort1
*cohort2name - string 
*phenotypeList - henotypeList is a list of tuples, in each tuple the first item is the phenotype name (string)the second is the conditionType-'identity' or 'delta',
the third is the delta (int) or None and the forth is isStopingCrit (True/False)
the order of the items in the list is according to their importance!
*columnsToUse: list of strings. phenotype columns to use for comparison - ***THIS COLUMNS SHOULD BE INCLUDED IN BOTH cohort1phenotypeDF AND 
cohort2phenotypeDF. generated by list comprehension from phenotypeList
*restrictedPhenotypes: list of strings. phenotypes with stopping indiction. generated by list comprehension from phenotypeList
*folderToSaveLists: string. path to save the list of samples (cohort1 and cohort2 matched samples)
*random_state: int. defailt=None. to keep df shuffling consistent. this number appears in the sampleList names.
imputeMissing=True/False. default=True


usage example:
f3='%s/TCR_real_data/NewPhenotypicData/PNP530ss12500_phen_new_dummies.xlsx' %MyPath
PNP530ss12500_phen_new_dummies=pd.read_excel(f3).set_index('BD')

f4='%s/TCR_real_data/CardioSamples/phenotypicData/Cardio126ss12500_phen_new_dummies.xlsx' %MyPath
Cardio126ss12500_phen_new_dummies=pd.read_excel(f4).set_index('BD')

phenotypeList=[('Gender_Male','identity',None,True),('Age','delta',9,True),('Smoking Status_Past','identity',None,True),('Smoking Status_Yes','identity',None,True),
               ('BMI','delta',4,False),('eGFR_CKD-EPI_new','delta',None,False),('WBC','delta',3,False)]
columnsToUse=[x[0] for x in phenotypeList]
restrictedPhenotypes=[x[0] for x in phenotypeList if x[3]==True]
print restrictedPhenotypes

cohort1phenotypeDF=Cardio126ss12500_phen_new_dummies
cohort1name='Cardio126ss12500'
cohort2phenotypeDF=PNP530ss12500_phen_new_dummies
cohort2name='PNP530ss12500'

folderToSaveLists='%s/TCR_real_data/PNP530Cardio126Combined/MatchedSamples/ss12500rep1' %MyPath
random_state=5
               
matchedCardio12500samples,matchedPNP12500samples=gen_matched_cohorts_by_phenotypes(cohort1phenotypeDF,cohort1name,cohort2phenotypeDF,cohort2name, 
                                      phenotypeList, columnsToUse,restrictedPhenotypes,folderToSaveLists,random_state)
                
#phenotypeList is a list of tuples, in each tuple the first item is the phenotype name (string)the second is the conditionType-'identity' or 'delta',
#the third is the delta (int) or None and the forth is isStopingCrit (True/False)

#matching samples from cohort2 to cohort1
#p




'''
def findClosestSample(potCohort2, cohort1organizedDF, sample):
#     sampleDF=pd.DataFrame(cohort1organizedDF)
    allSamples = pd.concat([potCohort2.copy(), pd.DataFrame(cohort1organizedDF.loc[sample, :]).T])
    print 'allSamples:'
    print allSamples
    colsToScale = []
    for col in allSamples.columns:
        allSamples[col] = allSamples[col].fillna(allSamples[col].median())
        if len(allSamples[col].value_counts()) > 2:
            colsToScale.append(col)
    allSamples[colsToScale] = preprocessing.scale(allSamples[colsToScale])
    # then, calculate euclidean distance between each potential sample and the sample-tobe-matched and 
    # choose the closetst sample:
    print 'now choosing the closest sample:'
    v2 = allSamples.iloc[-1, :]
    print 'v2: %s' % v2
    distDict = {}
    for n, s in enumerate(allSamples.index[:-1]):
        print n, s
        v1 = allSamples.loc[s, :]
        print 'v1: %s' % v1
        dist = euclidean(v1, v2)
        distDict[s] = dist
    sampleTaken = min(distDict, key=distDict.get)
    sampleTakenDF = pd.DataFrame(potCohort2.loc[sampleTaken, :])
    print 'sample %s was found to be the closet to sample %s, so it is selected' % (sampleTaken, sample)

    return sampleTaken, sampleTakenDF


def gen_matched_cohorts_by_phenotypes(cohort1phenotypeDF, cohort1name, cohort2phenotypeDF, cohort2name,
                                      phenotypeList, columnsToUse, restrictedPhenotypes, folderToSaveLists, random_state=None, imputeMissing=True,
                                      numericalPhenotypes=None):
                                      
    if not isdir(folderToSaveLists):
        makedirs(folderToSaveLists)
    
    
    # (1) generate matched cohort:
    print 'number of samples in cohort 1 is %s' % len(cohort1phenotypeDF)
    print 'number of samples in cohort 2 is %s' % len(cohort2phenotypeDF)
    nPhenotypes = len(phenotypeList)
       
    cohort1organizedDF = cohort1phenotypeDF.copy()[columnsToUse]
    cohort1organizedDF = cohort1organizedDF.sample(frac=1, random_state=random_state)
    cohort2organizedDF = cohort2phenotypeDF.copy()[columnsToUse]
    if imputeMissing:
        print 'imputing missing values for cohort %s' % cohort2name
        for col in cohort2organizedDF.columns:
            cohort2organizedDF[col] = cohort2organizedDF[col].fillna(cohort2organizedDF[col].median())
            
    # loop over all samples in cohort1 to find matched samples in cohort2:
    for n, sample in enumerate(cohort1organizedDF.index):
        sampleTaken = None
        sampleTakenDF = None
        cohort2organizedDF = cohort2organizedDF.sample(frac=1, random_state=random_state)
        print 'now the number of samples in cohort 2 is %s' % len(cohort2organizedDF)
        potCohort2 = cohort2organizedDF
        print n, sample
        print pd.DataFrame(cohort1organizedDF.loc[sample, :])
        
        sort = True
        for i, phen in enumerate(phenotypeList):
                
            phenotype = phen[0]
            conditionType = phen[1]
            maxDelta = phen[2]
            isStopping = phen[3]
            phenotypeValue = cohort1organizedDF.loc[sample, phenotype]
            print phenotype, conditionType, maxDelta, isStopping
            
            if (np.isnan(phenotypeValue)) and (i != nPhenotypes - 1):
                    print 'phenotype value is nan, skip this phenotype'
                    continue
            elif (np.isnan(phenotypeValue)) and (i == nPhenotypes - 1):
                print 'phenotype value is nan, and this is the last phenotype, thus taking the closest sample:'
                sampleTaken, sampleTakenDF = findClosestSample(potCohort2, cohort1organizedDF, sample)
                break

            if conditionType == 'identity':
                same = potCohort2[potCohort2[phenotype] == phenotypeValue]
                print 'len same=%s' % len(same)
                if len(same) < 1:
                    if isStopping:
                        print 'no more %s samples in %s cohort' % (phenotypeValue, cohort2name)
                        break
                    else:
                        potCohort2 = potCohort2
                else:
                    potCohort2 = same
                    
                    
            elif conditionType == 'delta':
                if maxDelta is None:
                        maxDelta = abs(phenotypeValue - potCohort2[phenotype]).max()
                same = potCohort2[(abs(potCohort2[phenotype] - phenotypeValue)) <= maxDelta]      
#                 minPhenDif=abs(phenotypeValue-potCohort2[phenotype]).min()
#                 if minPhenDif>maxDelta:
                if len(same) < 1:
                    if isStopping:
                        print 'no more samples with similar %s' % phenotype
                        break
                    else:
                        potCohort2 = potCohort2
                else: 
                    print 'len same=%s' % len(same)
                    potCohort2 = same
            else:
                print 'no conditionalType was specified!'
                
            print 'len potCohort=%s' % len(potCohort2)
            
            if (len(potCohort2) >= 2) and (i != nPhenotypes - 1):  # if there are more than one potential sample and it isn't the last phenotype:
                print 'len potCohort>2, continue to the next phenotype'
                continue
            elif len(potCohort2) < 2:  # if there is only one sample, take it and 
                sampleTakenDF = pd.DataFrame(potCohort2.iloc[0, :])
                sampleTaken = sampleTakenDF.columns[0]
                print 'len potCohort<2, so after phenotype %s, sample %s is taken' % (phenotype, sampleTaken)
                break
            else:  # if there are more than 1 potential sample and it is the last phenotype:
                 # take the closest sample:
                # first scale the potential samples toegther with the sample to be matched:
                print 'len potCohort>2 but this is the last phenotype'
                sampleTaken, sampleTakenDF = findClosestSample(potCohort2, cohort1organizedDF, sample)
                break

        if  sampleTaken is not None:   
            if n == 0:
                matchedCohort1sampleDF = pd.DataFrame(cohort1organizedDF.loc[sample, :])
                matchedCohort2sampleDF = sampleTakenDF
            else:
                matchedCohort2sampleDF = pd.merge(matchedCohort2sampleDF, sampleTakenDF, how='outer', left_index=True, right_index=True)
                matchedCohort1sampleDF = pd.merge(matchedCohort1sampleDF, pd.DataFrame(cohort1organizedDF.loc[sample, :]), how='outer', left_index=True, right_index=True)
            cohort2organizedDF = cohort2organizedDF.drop(sampleTaken)

  
    print 'done matching'
    print 'length of cohort1 matched sample DF is %s' % len(matchedCohort1sampleDF.columns)
#     print matchedCohort1sampleDF.iloc[:4,:4]
    print 'length of cohort2 matched sample DF is %s' % len(matchedCohort2sampleDF.columns)
#     print matchedCohort2sampleDF.iloc[:4,:4]

    matchedCohort1samples = matchedCohort1sampleDF.T.index.tolist()
    matchedCohort2samples = matchedCohort2sampleDF.T.index.tolist()
    
    print 'length of cohort1 matched sample list is %s' % len(matchedCohort1samples)
#     print matchedCohort1samples[:5]
    print 'length of cohort2 matched sample list is %s' % len(matchedCohort2samples)
#     print matchedCohort2samples[:5]
    
    # save matched Lists:
    with open('%s/%s_samples_matchedTo_%s_n%s_resBy%s_rs%s' % (folderToSaveLists, cohort1name, cohort2name, len(matchedCohort1samples),
                                                              ''.join(restrictedPhenotypes).replace('_', ''), random_state), 'wb') as fp1: 
        pickle.dump(matchedCohort1samples, fp1)
    with open('%s/%s_samples_matchedTo_%s_n%s_resBy%s_rs%s' % (folderToSaveLists, cohort2name, cohort1name, len(matchedCohort2samples),
                                                              ''.join(restrictedPhenotypes).replace('_', ''), random_state), 'wb')as fp2:
        pickle.dump(matchedCohort2samples, fp2)
        
    # (2) compare phenotypes
    if numericalPhenotypes is None:
        numericalphenotypes = ['Age', 'BMI', 'HbA1C', 'Total Cholesterol', 'Glucose', 'WBC']
    categoricalphenotypes = ['Gender', 'Smoking Status']
    
#     print 'cohort1phenotypeDF:'
#     print cohort1phenotypeDF.iloc[:4,:4]
#     print 'cohort2phenotypeDF:'
#     print cohort2phenotypeDF.iloc[:4,:4]
    
    
    
    phenotypeDF1 = cohort1phenotypeDF.loc[matchedCohort1samples, :]
    phenotypeDF2 = cohort2phenotypeDF.loc[matchedCohort2samples, :]
    sampleList1 = None
    sampleList2 = None
    sampleListName1 = None
    sampleListName2 = None
    datasetName1 = '%s_matched_phenotypeDF' % cohort1name   
    datasetName2 = '%s_matched_phenotypeDF' % cohort2name
    nBins = 20

    fig1 = compare_phenotypes(numericalphenotypes, categoricalphenotypes, nBins, phenotypeDF1, sampleList1, sampleListName1, datasetName1,
                           phenotypeDF2, sampleList2, sampleListName2, datasetName2)
                       
    return matchedCohort1samples, matchedCohort2samples


#----------------------------------------------------------------------------------
'''
function 3: when the matched cohorts are good:
7. generate a folder with all samples
8. compare TCR features
9. generate TCRdf
10. compare sharing rates between PNP and cardio
11. try to seperate by PCA on 10perc binary TCR
12. find unique sequences per cohort and their identity


input:
*ss - int, number of subsampled sequences (usually 5000/9000/12500/15000)
*repeat - int, number of repeat for the same ss number (usually 1)
*PNPmatchedSampleList-list of sample names, all PNP samples that were matched to the Cardio samples - output of function 2
*CardioMatchedSampleList-list of sample names, all cardio samples that were matched to the PNP samples - output of function 2
*copySamples - True/False - should or shouldn't copy all data files for the matched samples to a dedicated folder
*compareFeatures - True/False - should ot shouldnt compare TCR features between cohorts
*compareSharing = True/False - should or shouldn't compare the sepetate TCRdfs for PNP and cardio cohorts (all common for nSharedMin or more samples)
*genSepTCRdf - True/False - should or shouldn't generate the seperate TCRdf (no need when they already exist)
*nSharedMin - the minimum number of shared samples need to be defined as shared in the SEPERATE TCRdfs. affects the annotation analysis
*delFolderFirst - True/False, default=False - should or shouldn't delete the data file folder (the one generated under 'copySamples')
*genComTCRdf - True/False, default=True. should or shouldn't generate commonTCRdf - all public sequences in both cohorts TOGETHER, which are shared 
by at least percSharedCom percent of samples or more.
*n_comp - int,default - None. how many PCs should be calculated. not so important as we use only the first 2 for the figure. however, all others
can be inspected manually.
*percSharedCom - int. minimal percent of samples that should be shared by each sequence in order to be included in the commonTCRdf used for the
PCA analysis.
*identityColumnForPie='Epitope species_VDJDB'/'Pathology_McPAS'/'combined annotation_list_clean' - the column in the CDR3 identity table that should be used for the identity pie figure.
*dropnaFromIdentitiesPie - whether or not to include unidentified sequences in the identity analysis
*ThresholdToOther - float, usually 0.01-0.05. defiens what is the fraction of specific annotation out ofall annotations in order 
to be includes in the pie diagram and the chi square calculation. 
*identityDF - None or specific DF. if None, the regular CDR3 identity file is used. if not, the specified df is used (change to use
the 'popped' files instead)
identityDFname - string. use 'popped' for popped files.  name of the identity (annotation) df used
*RemoveSamplesWithRelative: True/ False/ done. if 'done': the function should be given sample list generated by matching procedure 
with PNP file from which relatives were removed. in this case all results will be generated in a slightly different folder. if True, will use the 
normal lists and folder, but in the sharing rates distribution calculations, the relatives will be removed. otherwise, no effect for samples 
from relatives. 

*******Note! sharing rate distribution comparison and PCA analysis are controled by the percSharedCom variable. if it is None, then all sequences
shared by 2 or more samples will be considered. the venn disgrams uses all sequences shared by 2 or more samples (can't be modified).
and the annotation analysis is controled by the nSahredMin variable, which indicates the minimum number of samples to use for sharing definition 
in order to use all sequences including provate, use nSharedMin=1 *************
usage example:
FeatureMeanSummary_all_dataSets=compareTCRforMatchedCohorts(ss,repeat,PNPmatchedSampleList,CardioMatchedSampleList,
                            copySamples,compareFeatures,compareSharing, genSepTCRdf,nSharedMin,delFolderFirst,
                            genComTCRdf,n_comp,percSharedCom,identityColumnForPie)
'''  
    


def compareTCRforMatchedCohorts(ss, repeat, PNPmatchedSampleList, CardioMatchedSampleList,
                                copySamples, compareFeatures, compareSharing, genSepTCRdf, nSharedMin,
                                delFolderFirst=False, genOrigComTCRdf=False, genComTCRdf=True, n_comp=None, percSharedCom=10,
                                identityColumnForPie='Pathology_McPAS', dropnaFromIdentitiesPie=True, ThresholdToOther=None,
                                identityDF=None, identityDFname=None, RemoveSamplesWithRelative=False):
    
    from shutil import copyfile, rmtree
    
    if RemoveSamplesWithRelative=='done':
        datasetFolder = '%s/TCR_real_data/PNP530Cardio126Combined/MatchedSamples/ss%srep%s_relativesRemoved' % (MyPath, ss, repeat)
    else:
        datasetFolder = '%s/TCR_real_data/PNP530Cardio126Combined/MatchedSamples/ss%srep%s' % (MyPath, ss, repeat)
    FeatureMeanSummary_all_dataSets = []
    
    if not isdir(datasetFolder):
        makedirs(datasetFolder)
    
    # (1) generate a folder with all matched samples:
    
    matchedSampleFolder = '%s/SamplesForAnalysis_corrected' % datasetFolder
    if copySamples:
        if delFolderFirst:
            rmtree(matchedSampleFolder)
        if not isdir(matchedSampleFolder):
            makedirs(matchedSampleFolder)
        # look for data files of all PNP samples that were matched to cardio sample and copy them to the common folder:
        PNPsourceFolder = '%s/TCR_real_data/PNP530_SubSampled%sdata_rep%s/SamplesForAnalysis_corrected' % (MyPath, ss, repeat)
        CardiosourceFolder = '%s/TCR_real_data/CardioSamples/Cardio126_SubSampled%sdata_rep%s/SamplesForAnalysis_corrected' % (MyPath, ss, repeat)

        folderTuple = [('PNP', PNPsourceFolder, PNPmatchedSampleList), ('Cardio', CardiosourceFolder, CardioMatchedSampleList)]
        for folder in folderTuple:
            print 'now copying files from %s folder' % folder[0]
            files = [f for f in listdir(folder[1]) if isfile(join(folder[1], f))]
            print 'number of %s files = %s' % (folder[0], len(files))
            count = 0                           
            for f in files:
                    f1 = f.replace('b', '')
                    f1 = f1.replace('a', '')
                    f1 = f1.split('_')[0]
                    f1 = f1.split('.')[0]
        #             print f1
                    if f1 in folder[2]:
                        src = '%s/%s' % (folder[1], f)
                        dst = '%s/%s' % (matchedSampleFolder, f)
                        copyfile(src, dst)
                        count = count + 1
                    else:
    #                     print 'file %s is not in the %s matched samples list' %(f,folder[0])
                        pass

            print '%s files were copied from %s folder' % (count, folder[0])
        print 'in total, %s files were copied' % len([f for f in listdir(matchedSampleFolder) if isfile(join(matchedSampleFolder, f))])
        print 'done copying sample data files'
    
    # (2) compare features:
    data_folder1 = 'TCR_real_data/CardioSamples/Cardio126_SubSampled%sdata_rep%s' % (ss, repeat)
    data_folder2 = 'TCR_real_data/PNP530_SubSampled%sdata_rep%s' % (ss, repeat)
    if compareFeatures:
        data_folder1 = 'TCR_real_data/CardioSamples/Cardio126_SubSampled%sdata_rep%s' % (ss, repeat)
        data_folder2 = 'TCR_real_data/PNP530_SubSampled%sdata_rep%s' % (ss, repeat)
        datasetName1 = 'Cardio126_ss%s_rep%s' % (ss, repeat)
        datasetName2 = 'PNP530_ss%s_rep%s' % (ss, repeat)

        TakeSameSamples = False
        filteringList1 = CardioMatchedSampleList
        filteringList2 = PNPmatchedSampleList
        filteringList1Name = 'CardioMatchedss%srep%s' % (ss, repeat)
        filteringList2Name = 'PNPMatchedss%srep%s' % (ss, repeat)
        MatchedFolderToSave = '%s/featureSummaryDFs' % datasetFolder

        print 'comparing features between datasets:'
        FeatureMeanSummary_all_dataSets = compare_features_between_datasets(data_folder1, datasetName1, data_folder2, datasetName2, TakeSameSamples,
                                          filteringList1, filteringList2, filteringList1Name, filteringList2Name, MatchedFolderToSave)

        print 'comparing geneset usage between datasets:'
        plotType = 'bar'
        plot_gene_usage_comparison(data_folder1, datasetName1, data_folder2, datasetName2, plotType, TakeSameSamples,
                                   filteringList1, filteringList2, filteringList1Name, filteringList2Name)

        # check significance of the results for different types of results:
        FeatComp = FeatureMeanSummary_all_dataSets.copy()
        FeatComp_T = FeatComp[FeatComp.index.str.contains('_T')]

        # # generate the geneUsageDF:
        regexList = [('vFamily', 'V.._T'), ('vGene', 'V..-.._T'), ('jGene', 'J..-.._T'), ('dFamily', 'D.._T'),
                  ('VJ', 'V.._J.._T'), ('DJ', 'D.._J..-.._T')]
        df_list = []
        for regex in regexList:
            regex_item = re.compile(regex[1])
            indices = [n for n in FeatComp.index if re.match(regex_item, n)]
            df = FeatComp.loc[indices, :]
            df_list.append(df)

        geneUsageResults_T = pd.concat(df_list)
        print len(geneUsageResults_T)

        nonGeneUsageResults_T = FeatComp.loc[[x for x in FeatComp_T.index if x not in geneUsageResults_T.index], :]
        print len(nonGeneUsageResults_T)

        # # check sig with FDR for each result df seperately:
        resultDFlist = [('all', FeatComp), ('onlyTotal', FeatComp_T), ('onlyGeneUsage_T', geneUsageResults_T), ('onlyFeatures_T',
                nonGeneUsageResults_T)]
        resultFolder = '%s/featureSummaryDFs/comparisonResults' % datasetFolder
        if not isdir(resultFolder):
            makedirs(resultFolder)
        count = 0
        for n, resultDF in enumerate(resultDFlist):
            print '******** %s *********' % resultDF[0]
            nTests = len(resultDF[1])
            print nTests
            FDR = 0.1

            pValueColumnList = ['ks_p', 't_p']
            for pValueColumn in pValueColumnList:
                print pValueColumn
                resultDF_FDR = add_corrected_pValues(resultDF[1], pValueColumn, nTests, FDR)
                resultDF_FDR = resultDF_FDR.rename(columns={'Sig by bonferroni corrected pVal':'%s_Sig by bonferroni corrected pVal'\
    % pValueColumn, 'sig. by FDR=0.1':'%s_sig. by FDR=0.1' % pValueColumn})
                sigResults = resultDF_FDR[resultDF_FDR['%s_sig. by FDR=0.1' % pValueColumn] == 1]
                print 'number of sig results with FDR=0.1 is %s' % len(sigResults)
                print 'sig results:'
                print sigResults
                f1 = '%s/sigResults_%s_%s.xlsx' % (resultFolder, resultDF[0], pValueColumn)
                sigResults.to_excel(f1)
                for col in sigResults.columns:
                    newCol = '%s_%s' % (col, resultDF[0])
                    sigResults = sigResults.rename(columns={col:newCol})
                print 'results were saved to file'
                if count == 0:
                    allSigResults = sigResults
                else:
                    allSigResults = pd.merge(allSigResults, sigResults, how='outer', left_index=True, right_index=True)
                count = count + 1
        for col in allSigResults.columns:
            if '_y' in col:
                allSigResults = allSigResults.drop(col, axis=1)
            if '_x' in col:
                newCol = col.replace('_x', '')
                allSigResults = allSigResults.rename(columns={col:newCol})
        f1 = '%s/allSigResults.xlsx' % resultFolder
        allSigResults.to_excel(f1)
                
            
                
    
    # (3) generate TCRdf, calculate sharing rates and seperate by PCA:
    # general definition, define anyway so they will be accessible: 
    print '*****generating %s perc TCRdfs******' % percSharedCom
    
    publicAnalysisFolder = '%s/publicAnalysis' % datasetFolder
    if not isdir(publicAnalysisFolder):
        makedirs(publicAnalysisFolder)
    
    DFtype = 'TCR'
    genDF = genOrigComTCRdf  # False=generate new df
    toBinary = True
    mbLevel = 'g'
    useShortName = True
    

    datasetName = 'MatchedSamples_ss%srep%s' % (ss, repeat)
    minVal = None  # minVal can be None,0, float, or 'dfMinVal' or dfMinVal2:
    minSharedT = None  # minimal number of samples shared by seq/species required to leave sample in the database (int or None)
    percShared = percSharedCom  # minimal number of samples shared by seq/species required to leave sample in the database (int [ percentage]
                    # or None)
    removeOutliers = True
    normData = True
    logTransform = True
    extractUniqueAA = True  # use True when this is the first time to analyze this dataset, otherwise, use False
    minNshared = 2
    onlyProductive = True
    mbDataFolder = 'AllSeqProjects'

    SampleList = [f for f in listdir(matchedSampleFolder) if isfile(join(matchedSampleFolder, f))]
    SampleListName = 'MatchedSamples_ss%srep%s' % (ss, repeat)
    libPrepMethod = None
    filterGenotek = True
    filterMinimalReads = 9000000
    filterlibPrepMethod = libPrepMethod
    groupFunction = 'noOutlierMean'
    nSTD = 5
    nMinSamples = 3
    ignoreNotSameYear = True
    removeSamePerson = False
    
    
    if genComTCRdf:
        print 'generating %sperc TCRdfs' % percSharedCom
        TCRdf0 = genTCRorMBdfWithManipulations(DFtype, genDF, toBinary, removeOutliers, normData, logTransform,
                                         minVal, minSharedT, percShared,
                                         mbLevel, useShortName, datasetFolder, datasetName, extractUniqueAA,
                                        minNshared, onlyProductive, mbDataFolder, SampleList,
                                          SampleListName, filterlibPrepMethod, filterGenotek,
                                          groupFunction, nSTD, nMinSamples, ignoreNotSameYear, removeSamePerson)
    
    print 'loading %sperc TCRdfs' % percSharedCom
    if percSharedCom is not None:
        TCRdfFileName = 'sharingMatrix_MatchedSamples_ss%srep%s_minNshared2_RA_onlyProductiveTrue__percShared%s_OLtrimmed_binary' % (ss, repeat, percSharedCom)
        TCRdfShortName = 'TCRdf_percShared%s_binary' % percSharedCom
    else:
        TCRdfFileName = 'sharingMatrix_MatchedSamples_ss%srep%s_minNshared2_RA_onlyProductiveTrue___OLtrimmed_binary' % (ss, repeat)
        TCRdfShortName = 'TCRdf_percShared%s_binary' % percSharedCom
        
    f1 = '%s/%s_sharingAnalysis/%s' % (datasetFolder, SampleListName, TCRdfFileName)
    TCRdf = pd.read_pickle(f1)
    print 'TCRdf shape is %s_%s' % (TCRdf.shape[0], TCRdf.shape[1])
    print 'TCRdf head:'
    print TCRdf.iloc[:4, :4]
    
    # #get PCAdf:
    # # **generate PCs as features:
    if n_comp is None:
        n_comp = 5
    try:
        TCRdf = TCRdf.set_index('BD')
    except:
        print 'TCRdf already has BD column as index'
    print 'generating PCAdf with %s PCs based on TCRdf:' % n_comp
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    print 'percSharedCom=%s' % percSharedCom
#     print 'percShared=%s' % percShared
    PCAdf, ax2, p_ttest_PC1, p_ttest_PC2 = PCAfunc(TCRdf, n_comp, isSparse=True, ax=ax2)
    ax2.set_title('PC1 - PC2 TCRdf %sperc-shared sequences' % percSharedCom)
    ax2.set_xlabel('PC1')
    ax2.set_ylabel('PC2')
    ax2.annotate('p_ttest_PC1=%s\np_ttest_PC2=%s' % (round(p_ttest_PC1, 4), round(p_ttest_PC2, 4)), xy=(0.96, 0.95), xycoords='axes fraction',
                    fontsize=12, horizontalalignment='right', verticalalignment='top', fontweight='bold')
    print 'PCAdf shape is %s_%s' % (PCAdf.shape[0], PCAdf.shape[1])
#     print 'PCAdf HEAD:'
    PCAplotFile = '%s/PC1PC2_%spercShared' % (publicAnalysisFolder, percSharedCom)
    fig2.savefig(PCAplotFile, dpi=200, bbox_inches="tight")
#     print PCAdf.iloc[:4,:4]
    
    # (4) compare public sequences between cohorts:
    # # number of public sequence per sample distributions, total numbers of shared sequences, identity of shared sequences
    if compareSharing:
        TCRdfInfoTuple = [('PNP', PNPmatchedSampleList, 'PNPmatchedSampleList'), ('Cardio', CardioMatchedSampleList, 'CardioMatchedSampleList')]
        TCRdfDict = {}
        uniqueSeqDict = {}
#         minSharedT = 2
        if percSharedCom is not None:
            percShared = percSharedCom 
            minSharedT = None
        else:
            percShared = None 
            minSharedT = 2
        if genSepTCRdf:
            # loop over PNP and Cardio matched samples (seperately) and generate TCRdf of shared sequences for each, seperately
            for TCRdfInfo in TCRdfInfoTuple:
                removeOutliers = True
                normData = False
                logTransform = False
    #             extractUniqueAA=False

                SampleList = TCRdfInfo[1]
                SampleListName = TCRdfInfo[2]

                print 'generating TCRdf only for %s samples from the matched samples:' % TCRdfInfo[0]
                TCRdfDict['%s' % TCRdfInfo[0]] = genTCRorMBdfWithManipulations(DFtype, genDF, toBinary, removeOutliers, normData, logTransform,
                                                     minVal, minSharedT, percShared,
                                                     mbLevel, useShortName, datasetFolder, datasetName, extractUniqueAA,
                                                    minNshared, onlyProductive, mbDataFolder, SampleList,
                                                      SampleListName, filterlibPrepMethod, filterGenotek,
                                                      groupFunction, nSTD, nMinSamples, ignoreNotSameYear, removeSamePerson)

        # #load binary DFs and compare number of shared sequences per sample distribution between PNP and cardio:
        # # +load uniqueAA file for each
        for TCRdfInfo in TCRdfInfoTuple:     
            print 'loading seperate minSharedT2 TCRdfs for %s' % TCRdfInfo[0]
            if percSharedCom is not None:
                TCRdfFileName = 'sharingMatrix_%s_minNshared2_RA_onlyProductiveTrue__percShared%s_OLtrimmed_binary' %(datasetName,percSharedCom)
            else:
                TCRdfFileName = 'sharingMatrix_%s_minNshared2_RA_onlyProductiveTrue__minSharedT2_OLtrimmed_binary' % datasetName
            f2 = '%s/%s_sharingAnalysis/%s' % (datasetFolder, TCRdfInfo[2], TCRdfFileName)
            TCRdfDict['%s' % TCRdfInfo[0]] = pd.read_pickle(f2)
            print '%s TCRdf shape is %s_%s' % (TCRdfInfo[0], TCRdfDict['%s' % TCRdfInfo[0]].shape[0], TCRdfDict['%s' % TCRdfInfo[0]].shape[1])
            print '%s TCRdf head:' % TCRdfInfo[0]
            print TCRdfDict['%s' % TCRdfInfo[0]].iloc[:4, :4]
            
            print 'loading number of unique sequences per sample:'
            if TCRdfInfo[0] == 'PNP':
                f3 = '%s/%s/featureSummaryDFs/PNP530_ss%s_rep%s_allFeatures' % (MyPath, data_folder2, ss, repeat)
            else:
                f3 = '%s/%s/featureSummaryDFs/Cardio126_ss%s_rep%s_allFeatures' % (MyPath, data_folder1, ss, repeat)
            featuresDF = pd.read_pickle(f3)
            featuresDF = editSampleNames(featuresDF)
            nUnique = pd.DataFrame(featuresDF['NT count_1'])
            
            print 'loading uniqueSeqWithCounts....'
            f4 = '%s/%s_sharingAnalysis/AllUniqueWithCounts' % (datasetFolder, TCRdfInfo[2])
            AllUniqueWithCounts = pd.read_pickle(f4)
            AllUniqueWithCounts = AllUniqueWithCounts[AllUniqueWithCounts['nShared'] >= nSharedMin]  # take only shared by more than nShared
            UniqueSequencesWithCounts = AllUniqueWithCounts.drop(['Sample', 'frequencyCount (%)', 'prod_stat', 'isPublic'], axis=1)
            UniqueSequencesWithCounts = UniqueSequencesWithCounts[~UniqueSequencesWithCounts.index.duplicated(keep='first')]   
            uniqueSeqDict['%s' % TCRdfInfo[0]] = UniqueSequencesWithCounts
            
            if TCRdfInfo[0] == 'PNP':
                PNP_nShared = pd.merge(pd.DataFrame(TCRdfDict['PNP'].sum(axis=1)), nUnique, how='left', left_index=True, right_index=True)
            else:
                Cardio_nShared = pd.merge(pd.DataFrame(TCRdfDict['Cardio'].sum(axis=1)), nUnique, how='left', left_index=True, right_index=True)
       
        # calculate percentage of public sequences out of all unique and plot distributions of this statistic in both cohorts:
        PNP_nShared['percPublic'] = 100 * PNP_nShared[0] / PNP_nShared['NT count_1']
        Cardio_nShared['percPublic'] = 100 * Cardio_nShared[0] / Cardio_nShared['NT count_1']
         
        data1 = Cardio_nShared['percPublic']
        data2 = PNP_nShared['percPublic']
        data1 = data1[data1.notnull()].tolist()
        data2 = data2[data2.notnull()].tolist()
        data1Name = 'Cardio'
        data2Name = 'PNP'
        folderToSave = publicAnalysisFolder
        fig1, ax = plt.subplots(figsize=(8, 6))
        if percSharedCom is not None:
            title = 'Sharing rates per sample - comparison between cohorts\n(Shared \
among %s percent of samples or more)' % percSharedCom
        else:
            title='Sharing rates per sample - comparison between cohorts\n(Shared \
among 2 samples or more)'
        showLegend = True
        nBins = 20
        dataList = [(data1Name, data1), (data2Name, data2)]
#          - a list of tuples, each tuple is build of (dataName,data)
        
        ax, ks_p_cohort1_cohort2, t_p_cohort1_cohort2, p_Anov, filename = plotHistComprison(dataList, ax, title, showLegend, nBins)
#         plotHistComprison(dataList, ax, title, showLegend=True, nBins=20, toAnnotate=True, alpha=None, plotType='hist')
        filename=filename+'_percShared%s' %percSharedCom
        if RemoveSamplesWithRelative:
            filename=filename+'_relativesRemoved' 
        
        print 'figure filename is %s' % filename
        fig1.savefig('%s/%s' % (publicAnalysisFolder, filename), dpi=200, bbox_inches="tight")
        plt.show()
       
        # merge tables and get statistics : how much are public in one cohort but not the other:
        mergedPublicList = pd.merge(uniqueSeqDict['PNP'], uniqueSeqDict['Cardio'], how='outer', left_index=True, right_index=True)
        for col in mergedPublicList.columns:
            newCol = col.replace('_x', 'PNP')
            newCol = newCol.replace('_y', 'Cardio')
            mergedPublicList = mergedPublicList.rename(columns={col:newCol})
#         print mergedPublicList.head(20)
        
        publicInBoth = mergedPublicList[(mergedPublicList['nSharedPNP'].notnull()) & (mergedPublicList['nSharedCardio'].notnull())]
        publicOnlyInPNP = mergedPublicList[(mergedPublicList['nSharedPNP'].notnull()) & (mergedPublicList['nSharedCardio'].isnull())]
        publicOnlyInCardio = mergedPublicList[(mergedPublicList['nSharedPNP'].isnull()) & (mergedPublicList['nSharedCardio'].notnull())]
        
        print 'number of sequence shared by %s or more samples in both cohorts: %s (%s perc)' % (nSharedMin,
                                                    len(publicInBoth), 100 * float(len(publicInBoth)) / len(mergedPublicList))
        print 'number of sequence shared by %s or more samples only in the PNP cohort: %s (%s perc)' % (nSharedMin,
                                                            len(publicOnlyInPNP), 100 * float(len(publicOnlyInPNP)) / len(mergedPublicList))
        print 'number of sequence shared by %s or more samples only in the Cardio cohort: %s (%s perc)' % (nSharedMin,
                                                            len(publicOnlyInCardio), 100 * float(len(publicOnlyInCardio)) / len(mergedPublicList))
        print 'total=%s' % (len(publicInBoth) + len(publicOnlyInPNP) + len(publicOnlyInCardio))
        print 'merged public list length length=%s' % len(mergedPublicList)
        
        # plot venn diagram:
        from matplotlib_venn import venn2
        fig3, ax3 = plt.subplots(figsize=(8, 6))
        venn = venn2([set(uniqueSeqDict['PNP'].index.tolist()), set(uniqueSeqDict['Cardio'].index.tolist())],
                   set_labels=('PNP', 'Cardio'), ax=ax3)
        ax3.set_title('Number of Public sequences in each cohort')
        VennPlotFile = '%s/PublicSeq_VennPlot_nSharedMin%s' % (publicAnalysisFolder, nSharedMin)
        fig3.savefig(VennPlotFile, dpi=200, bbox_inches="tight")
        plt.show()
        
        # add sequence identities:
        print 'analyzing public sequence identities...'
        if identityDF is None:
            identsProcessedFile = '%s/TCR CDR3 sequence databases/CDR3identityTable_06082018_processed.xlsx' % MyPath
    #         identsDF=pd.read_excel(identsFile)
            identsProcessedDF = pd.read_excel(identsProcessedFile)
        else:
            identsProcessedDF = identityDF
        
        
        publicDFlist = [('AllDB', identsProcessedDF), ('Both', publicInBoth), ('onlyPNP', publicOnlyInPNP), ('onlyCardio', publicOnlyInCardio)]
        fig4, axes = plt.subplots(nrows=2, ncols=2, figsize=(24, 22))
        valueCountDict = {}
        for n, publicDF in enumerate(publicDFlist):
            print n, publicDF[0]
            if n != 0:
#                 publicDFwithIdents = pd.merge(publicDF[1], identsDF, how='left', left_index=True, right_index=True)
#                 f1 = '%s/publicDFwithIdents_%s.xlsx' % (publicAnalysisFolder, publicDF[0])
#                 publicDFwithIdents.to_excel(f1)
                publicDFwithIdents_processed = pd.merge(publicDF[1], identsProcessedDF, how='left', left_index=True, right_index=True)
                df = publicDFwithIdents_processed
            else:
                df = publicDF[1]
            
            dfName = publicDF[0]
            ax = axes.flatten()[n]
            nSharedThreshold = None
            useMore = True
            column = identityColumnForPie
            dropna = dropnaFromIdentitiesPie          
            size = None
            if ThresholdToOther is None:
                ThresholdToOther = 0.05
            ax, valueCountsgrouped2 = plot_identity_pie_plot(ax, df, dfName, nSharedThreshold, useMore, column, dropna, ThresholdToOther)

            valueCountDict[publicDF[0]] = valueCountsgrouped2
#             print valueCountsgrouped2
        if 'popped' in identityDFname:
            publicIdentsFigFile = '%s/publicIdentsFig_%s_sharedBy%s_dropna%s_thOther%s_popped' % (publicAnalysisFolder,
                                                                      identityColumnForPie, nSharedMin, dropnaFromIdentitiesPie, ThresholdToOther)
        else:
            publicIdentsFigFile = '%s/publicIdentsFig_%s_sharedBy%s_dropna%s_thOther%s' % (publicAnalysisFolder,
                                                                      identityColumnForPie, nSharedMin, dropnaFromIdentitiesPie, ThresholdToOther)
        publicIdentsFigFile = publicIdentsFigFile.replace('.', '-')
                
        # generate contingency table and calculate chi test:
        onlyPNPvalueCount = pd.DataFrame(valueCountDict['onlyPNP']).rename(columns={identityColumnForPie: 'onlyPNP'})
        onlyCardiovalueCount = pd.DataFrame(valueCountDict['onlyCardio']).rename(columns={identityColumnForPie: 'onlyCardio'})
        BothvalueCount = pd.DataFrame(valueCountDict['Both']).rename(columns={identityColumnForPie: 'Both'})
        
        
        valueCountComb = pd.merge(onlyPNPvalueCount, onlyCardiovalueCount, how='outer',
                               left_index=True, right_index=True)
        valueCountComb = pd.merge(valueCountComb, BothvalueCount, how='outer',
                               left_index=True, right_index=True)
        valueCountComb2 = valueCountComb.T.fillna(0)
        
        print valueCountComb2
        
#         print valueCountComb2
        from scipy.stats import chi2_contingency
        comparisonList = [('onlyPNP', 'onlyCardio'), ('onlyPNP', 'Both'), ('onlyCardio', 'Both')]
        pDict = {}
        for item in comparisonList:
            contTable = valueCountComb2.loc[[item[0], item[1]], :]
            print contTable
            for col in contTable.columns:
                if contTable[col].sum() < 1:
                    contTable = contTable.drop(col, axis=1)
            chi, p, dof, expctd = chi2_contingency(contTable)
            pDict['_'.join([item[0], item[1]])] = round(p, 4)
                    
        ax.annotate('p_chi square tests:\nonly PNP vs. only Cardio=%s\nonly PNP vs. Both=%s\nonly Cardio vs. Both=%s' % (pDict['onlyPNP_onlyCardio'],
                    pDict['onlyPNP_Both'], pDict['onlyCardio_Both']), xy=(0.96, 0), xycoords='axes fraction',
                    fontsize=18, horizontalalignment='right', verticalalignment='top', fontweight='bold')
        fig4.subplots_adjust(left=0.09, right=0.9, top=0.9, bottom=0.2, wspace=0.5, hspace=0.25)
        fig4.savefig(publicIdentsFigFile, dpi=200, bbox_inches="tight")
        plt.show()
        
        print 'end of function!'
    
    return FeatureMeanSummary_all_dataSets
    
    

    













