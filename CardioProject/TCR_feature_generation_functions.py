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

MyPath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF'
#-----------------------------------------------------------------------------


# # define functions to calculate statistics
# ## non-gene usage functions:

# this function gets as input the tuple (df type, df data)
# df type can be 'Total', 'Prod' or 'nonProd', and the df data is the dataframe itself. 


def gen_LengthFeaturesAndMore(df, data_folder, sample_name):
    
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/LengthFeaturesAndMore/%s' % (data_folder, df[0], sample_name) 
    dfs_folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/LengthFeaturesAndMore' % (data_folder, df[0])
    if not isdir(dfs_folder):
        makedirs(dfs_folder) 
    files = [f for f in listdir(dfs_folder) if isfile(join(dfs_folder, f))]

    
    if sample_name not in files:
        print('calculating length features...')
        DF = df[1]  # DF is only the data
        DF = DF.drop(['nucleotide', 'aminoAcid'], axis=1)  # drop those columns as they interfere with the mean,max,std calculations
        LengthFeaturesAndMore = pd.DataFrame()  # generate empty dataframe
        LengthFeaturesAndMore.loc[0, 'Sample'] = sample_name 
        LengthFeaturesAndMore.loc[0, 'df type'] = df[0]
        for column in DF.columns.values:  # loop over each column and calculate its mean, std and max, store information in the dataframe
            LengthFeaturesAndMore.loc[0, '%s_mean' % column] = DF[column].mean()
            LengthFeaturesAndMore.loc[0, '%s_std' % column] = DF[column].std()
            LengthFeaturesAndMore.loc[0, '%s_max' % column] = DF[column].max()
    #     if df[0]=='Total':
    #         LengthFeaturesAndMore.loc[0,'perc_prod']=DF['prod_stat'].mean()
    #     else:
    #         LengthFeaturesAndMore.loc[0,'perc_prod']=np.nan

        # save dataframe for each df type in each sample that contains length features:
        LengthFeaturesAndMore.to_pickle(file1)
    else:
        print('found length features for this sample...')
        

# calculate the mean %gc in all NT sequences in a specific df:
# this function is used within gen_generalFeatures
def gc_content(df):
    seqs = list(df['nucleotide'])
    gc_values = [GC(seq) for seq in seqs]
    mean_gc = np.mean(gc_values)
    return mean_gc


# this function gets as input the tuple (df type, df data)
# df type can be 'Total', 'Prod' or 'nonProd', and the df data is the dataframe itself. 

def gen_generalFeatures(df, data_folder, sample_name):
    
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/GeneralFeatures/%s' % (data_folder, df[0], sample_name)   
    dfs_folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/GeneralFeatures' % (data_folder, df[0])
    if not isdir(dfs_folder):
            makedirs(dfs_folder)    
    files = [f for f in listdir(dfs_folder) if isfile(join(dfs_folder, f))]
    
    if sample_name not in files:
        print('calculating general features...')

        generalFeatures = pd.DataFrame()  # generate empty dataframe
        generalFeatures.loc[0, 'Sample'] = sample_name
        generalFeatures.loc[0, 'df type'] = df[0]
        # calculate general features:
        generalFeatures.loc[0, 'NT count'] = len(df[1])
        n_aa_in_total = len(df[1].groupby('aminoAcid'))
        generalFeatures.loc[0, 'AA count'] = n_aa_in_total
        nt_per_aa = df[1].groupby(['aminoAcid'])[['nucleotide']].count()  # for each aa sequence, count how many nt sequence generate the same one
        max_nt_per_aa = max(nt_per_aa['nucleotide'])
        mean_nt_per_aa = round(np.mean(nt_per_aa['nucleotide']), 3)
        generalFeatures.loc[0, 'max_nt_per_aa'] = max_nt_per_aa
        generalFeatures.loc[0, 'mean_nt_per_aa'] = mean_nt_per_aa
        generalFeatures.loc[0, 'gc_content'] = gc_content(df[1])


        # save dataframe for each df type in each sample that contains general  features:
        generalFeatures.to_pickle(file1)
    
    else:
        print('found General Features for this sample...')
        


# this function counts the number of unique sequences per defined number of templates, using random sampling without replacement:
def norm_uniqe_nt_sequences(df,samp_size_nt=None): 
    repeats = 10
    if samp_size_nt is None:
        samp_size_nt = 2000  # this number was selected based on the PNP cohort data, can be changed but need to change the column name in the df!
    reads = list(df['count (templates)'])
    df = df.set_index('nucleotide')
    seqs = [str(i) for i in list(df.index)]
    seq_popped = []
    for i in range(0, len(seqs)):
        for j in range(0, reads[i]):
            seq_popped.append(seqs[i])        
    seq_n_list = []
    n_templates_NT = len(seq_popped)
    if len(seq_popped) > samp_size_nt:
        for t in range(repeats):  # as the calculation has a random component, repeat this calculation several times and average
            rand_seq = np.random.choice(seq_popped, samp_size_nt, replace=False)
            seq_n = len(set(rand_seq))
            seq_n_list.append(seq_n)
            mean_seq_n = np.mean(seq_n_list)
    else:
        mean_seq_n = np.NaN
        
    return mean_seq_n 



# this function counts the number of unique sequences per defined number of templates, using random sampling without replacement:
def norm_uniqe_aa_sequences(df,samp_size_aa=None): 
    repeats = 10
    if samp_size_aa is None:
        samp_size_aa = 200  # this number was selected based on the PNP cohort data, can be changed but need to change the column name in the df!
    reads = list(df['count (templates)'])
    list_aa = list(df['aminoAcid'])  
    seq_aa_popped = []
    for i in range(0, len(list_aa)):
        for j in range(0, reads[i]):
            seq_aa_popped.append(list_aa[i])
    seq_aa_popped_new = [i for i in seq_aa_popped if isinstance(i, str)]
    seq_n_list = []
    n_templates_aa = len(seq_aa_popped_new)
    if len(seq_aa_popped_new) > samp_size_aa:
        for t in range(repeats):  # as the calculation has a random component, repeat this calculation several times and average
            rand_seq = np.random.choice(seq_aa_popped_new, samp_size_aa, replace=False)
            seq_n = len(set(rand_seq))
            seq_n_list.append(seq_n)
        mean_seq_n = np.mean(seq_n_list)
    else:
        mean_seq_n = np.NaN
    return mean_seq_n         



def gen_normSeqNums(df, data_folder, sample_name,samp_size_nt=None,samp_size_aa=None):
    
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/normSeqNums/%s' % (data_folder, df[0], sample_name)   
    dfs_folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/normSeqNums' % (data_folder, df[0])
    if not isdir(dfs_folder):
            makedirs(dfs_folder)  
    files = [f for f in listdir(dfs_folder) if isfile(join(dfs_folder, f))]
    
    if sample_name not in files:
        print('calculating normSeqNums features...')

        normSeqNums = pd.DataFrame()  # generate empty dataframe
        normSeqNums.loc[0, 'Sample'] = sample_name
        normSeqNums.loc[0, 'df type'] = df[0]

        DF = df[1]

        # calculate general features:
        normSeqNums.loc[0, 'normSeqNums_per%s_NT' %samp_size_nt] = norm_uniqe_nt_sequences(DF,samp_size_nt)
        normSeqNums.loc[0, 'normSeqNums_per%s_AA' %samp_size_aa] = norm_uniqe_aa_sequences(DF, samp_size_aa)

        # save dataframe for each df type in each sample that contains general  features:

        normSeqNums.to_pickle(file1)
    else:
        print('found normSeqNums for this sample...')
        

# this function gets as input the tuple (df type, df data)
# df type can be 'Total', 'Prod' or 'nonProd', and the df data is the dataframe itself. 

def gen_clonalityFeatures(df, data_folder, sample_name):
    
    
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/ClonalityFeatures/%s' % (data_folder, df[0], sample_name)
    dfs_folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/ClonalityFeatures' % (data_folder, df[0])
    if not isdir(dfs_folder):
            makedirs(dfs_folder)   
    files = [f for f in listdir(dfs_folder) if isfile(join(dfs_folder, f))]
    
    if sample_name not in files:
        print('calculating clonality features...')

        DF = df[1]  # DF is only the data
        DF_aa = pd.DataFrame(DF.groupby('aminoAcid').sum()['frequencyCount (%)'])
        DF['clonality_nt'] = DF['frequencyCount (%)'] / DF['frequencyCount (%)'].sum()
        DF_aa['clonality_aa'] = DF_aa['frequencyCount (%)'] / DF_aa['frequencyCount (%)'].sum()
        clonalityList_nt = list(DF['clonality_nt'].sort_values(ascending=False))
        clonalityList_aa = list(DF_aa['clonality_aa'].sort_values(ascending=False))

        if len(clonalityList_nt) < 10:
            top10clonal_nt = np.nan
            top1000clonal_nt = np.nan
        elif len(clonalityList_nt) < 1000:
            top10clonal_nt = np.sum(clonalityList_nt[:10])
            top1000clonal_nt = np.nan
        else:
            top10clonal_nt = np.sum(clonalityList_nt[:10])
            top1000clonal_nt = np.sum(clonalityList_nt[:1000])


        if len(clonalityList_aa) < 10:
            top10clonal_aa = np.nan
            top1000clonal_aa = np.nan
        elif len(clonalityList_aa) < 1000:
            top10clonal_aa = np.sum(clonalityList_aa[:10])
            top1000clonal_aa = np.nan
        else:
            top10clonal_aa = np.sum(clonalityList_aa[:10])
            top1000clonal_aa = np.sum(clonalityList_aa[:1000])

#         topClonal_nt = np.max(clonalityList_nt)
#         meanClonal_nt = np.mean(clonalityList_nt)
#         medianClonal_nt = np.median(clonalityList_nt)
#         stdClonal_nt = np.std(clonalityList_nt)
#         percentile1_nt = np.percentile(clonalityList_nt, 1)
#         percentile999_nt = np.percentile(clonalityList_nt, 99.9)

#         topClonal_aa = np.max(clonalityList_aa)
#         meanClonal_aa = np.mean(clonalityList_aa)
#         medianClonal_aa = np.median(clonalityList_aa)
#         stdClonal_aa = np.std(clonalityList_aa)
#         percentile1_aa = np.percentile(clonalityList_aa, 1)
#         percentile999_aa = np.percentile(clonalityList_aa, 99.9)

        clonFeatures = pd.DataFrame()  # generate empty dataframe
        clonFeatures.loc[0, 'Sample'] = sample_name
        clonFeatures.loc[0, 'df type'] = df[0]
#         clonFeatures.loc[0, 'topClonal_nt'] = topClonal_nt
#         clonFeatures.loc[0, 'meanClonal_nt'] = meanClonal_nt
#         clonFeatures.loc[0, 'medianClonal_nt'] = medianClonal_nt
#         clonFeatures.loc[0, 'stdClonal_nt'] = stdClonal_nt
        clonFeatures.loc[0, 'top10clonal_nt'] = top10clonal_nt
        clonFeatures.loc[0, 'top1000clonal_nt'] = top1000clonal_nt
#         clonFeatures.loc[0, 'percentile1_nt'] = percentile1_nt
#         clonFeatures.loc[0, 'percentile999_nt'] = percentile999_nt
#         clonFeatures.loc[0, 'topClonal_aa'] = topClonal_aa
#         clonFeatures.loc[0, 'meanClonal_aa'] = meanClonal_aa
#         clonFeatures.loc[0, 'medianClonal_aa'] = medianClonal_aa
#         clonFeatures.loc[0, 'stdClonal_aa'] = stdClonal_aa
        clonFeatures.loc[0, 'top10clonal_aa'] = top10clonal_aa
        clonFeatures.loc[0, 'top1000clonal_aa'] = top1000clonal_aa
#         clonFeatures.loc[0, 'percentile1_aa'] = percentile1_aa
#         clonFeatures.loc[0, 'percentile999_aa'] = percentile999_aa

        # save dataframe for each df type in each sample that contains general  features:
        clonFeatures.to_pickle(file1)
        
    else:
        print('found clonality features for this sample...')

    
  # this function gets as input the tuple (df type, df data)
# df type can be 'Total', 'Prod' or 'nonProd', and the df data is the dataframe itself. 

def gen_diversityFeatures(df, data_folder, sample_name):
    
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/DiversityFeatures/%s' % (data_folder, df[0], sample_name)
    dfs_folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/DiversityFeatures' % (data_folder, df[0])
    if not isdir(dfs_folder):
            makedirs(dfs_folder)
    files = [f for f in listdir(dfs_folder) if isfile(join(dfs_folder, f))]
    
    if sample_name not in files:
        print('calculating diversity features...')

        DF = df[1]  # DF is only the data
        DF_aa = pd.DataFrame(DF.groupby('aminoAcid').sum()['count (templates)'])

        shannon_nt = shannon(DF['count (templates)'], base=2)
        simpson_nt = simpson(DF['count (templates)'])
        berger_nt = berger_parker_d(DF['count (templates)'])
        shannon_aa = shannon(DF_aa['count (templates)'], base=2)
        simpson_aa = simpson(DF_aa['count (templates)'])
        berger_aa = berger_parker_d(DF_aa['count (templates)'])

        diversityFeatures = pd.DataFrame()  # generate empty dataframe
        diversityFeatures.loc[0, 'Sample'] = sample_name
        diversityFeatures.loc[0, 'df type'] = df[0]
        diversityFeatures.loc[0, 'shannon_nt'] = shannon_nt
        diversityFeatures.loc[0, 'simpson_nt'] = simpson_nt
        diversityFeatures.loc[0, 'berger_nt'] = berger_nt
        diversityFeatures.loc[0, 'shannon_aa'] = shannon_aa
        diversityFeatures.loc[0, 'simpson_aa'] = simpson_aa
        diversityFeatures.loc[0, 'berger_aa'] = berger_aa

        # save dataframe for each df type in each sample that contains general  features:
        diversityFeatures.to_pickle(file1)
    else:
        print('found diversityFeatures for this sample...')
        
        
# ## gene usage functions:

# # this function counts independent gene usage **without unresolved!!**:


def count_geneUsage(df, data_folder, sample_name):
    print('counting gene usage features...')
    
    param_list = ['vGeneName', 'vFamilyName', 'dFamilyName', 'jGeneName']
    DF = df[1]
    
    for param in param_list:
        file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/%s/%s' % (data_folder, df[0], param, sample_name)
        dfs_folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/%s/' % (data_folder, df[0], param)
        if not isdir(dfs_folder):
            makedirs(dfs_folder)
        files = [f for f in listdir(dfs_folder) if isfile(join(dfs_folder, f))]
    
        if sample_name not in files:
            print('calculating gene usage for param %s' % param)
            GeneUsage = pd.DataFrame(DF[param].value_counts(normalize=True))
            GeneUsage = GeneUsage.rename(columns={'%s' % param: '%s' % sample_name})
            GeneUsage = GeneUsage.T
            GeneUsage.loc[sample_name, 'df type'] = df[0]
            GeneUsage.to_pickle(file1)
        else:
            print('found gene usage count for %s' % param)

            
# # this function counts dependent (combined) gene usage **without unresolved!!**:
def count_dependent_geneUsage(df, data_folder, sample_name):
    print('counting dependent gene usage features...')
    
    dep_param_list = ['V-J family combination', 'D-J gene combination']
    DF = df[1]
    DF['V-J family combination'] = DF['vFamilyName'] + '_' + DF['jFamilyName']
    DF['D-J gene combination'] = DF['dFamilyName'] + '_' + DF['jGeneName']
    
    for dep_param in dep_param_list:
        file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/%s/%s' % (data_folder, df[0], dep_param, sample_name)        
        dfs_folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/%s/' % (data_folder, df[0], dep_param)
        
        if not isdir(dfs_folder):
            makedirs(dfs_folder) 
        files = [f for f in listdir(dfs_folder) if isfile(join(dfs_folder, f))]
    
        if sample_name not in files:
            print('calculating dep gene usage for dep param %s' % dep_param)

            depGeneUsage = pd.DataFrame(DF[dep_param].value_counts(normalize=True))
            depGeneUsage = depGeneUsage.rename(columns={'%s' % dep_param: '%s' % sample_name})
            depGeneUsage = depGeneUsage.T
            depGeneUsage.loc[sample_name, 'df type'] = df[0]
            depGeneUsage.to_pickle(file1)
        else:
            print('found dep gene usage count for %s' % dep_param)  

#----------------------------------------------------------------------------------------------------
'''
the following functions are higher level functions that uses the former functions.

usage example:

data_folder='TCR_real_data'
dfs_folder='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/SamplesForAnalysis_corrected' %data_folder
filenames = [f for f in listdir(dfs_folder) if isfile(join(dfs_folder, f))]
filenames=[f.strip('.tsv') for f in filenames]
print len(filenames)

for n,sample_name in enumerate(filenames):    
        print  n,sample_name
        gen_descriptive_stats(sample_name)
        gen_geneUsageCount(sample_name)
'''




# # define highest function to call all feature calculating functions:

def gen_descriptive_stats(sample_name, data_folder, newColumnList,samp_size_nt=None,samp_size_aa=None):
    # (1) read sample data, add indications for productive, get only interesting
    # columns and generate dfs for total, only productive, only non productive

    file1 = "/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/SamplesForAnalysis_corrected/%s.tsv" % (data_folder, sample_name)
    if isfile(file1):
        sample_df = pd.read_table(file1)
    else:
        file2 = "/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/SamplesForAnalysis_corrected/%s.xlsx" % (data_folder, sample_name)
        sample_df = pd.read_excel(file2)
    if newColumnList is not None:  # make sure column names are correct:
        sample_df = sample_df.iloc[:, :44]
        sample_df.columns = newColumnList
    sample_df = sample_df.rename(columns={'count (templates/reads)':'count (templates)', 'count (reads)':'count (templates)'})
    
#     print sample_df['sequenceStatus'].head()
    sample_df['prod_stat'] = np.where(sample_df['sequenceStatus'].astype(str) == 'In', 1, 0)
    interesting_columns = ['nucleotide', 'aminoAcid', 'count (templates)', 'frequencyCount (%)', 'cdr3Length',
                         'vDeletion', 'n1Insertion', 'd5Deletion', 'd3Deletion',
                         'n2Insertion', 'jDeletion', 'prod_stat']    

    sample_df = sample_df[interesting_columns]
    prod_df = sample_df[sample_df['prod_stat'] == 1]
    nonprod_df = sample_df[sample_df['prod_stat'] == 0]
    
    df_dict = [('Total', sample_df),
     ('Prod', prod_df),
     ('nonProd', nonprod_df)]
 
 # (2) call selected functions to calculate statistics for each df
    for df in df_dict:
        print df[0]
        gen_LengthFeaturesAndMore(df, data_folder, sample_name)
        gen_generalFeatures(df, data_folder, sample_name)
        gen_normSeqNums(df, data_folder, sample_name,samp_size_nt,samp_size_aa)
        gen_clonalityFeatures(df, data_folder, sample_name)
        gen_diversityFeatures(df, data_folder, sample_name)
        

# geneUsage count has its own function as the column needed are different

def gen_geneUsageCount(sample_name, data_folder, newColumnList):
    # (1) read sample data, add indications for productive, get only interesting
    # columns and generate dfs for total, only productive, only non productive
    file1 = "/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/SamplesForAnalysis_corrected/%s.tsv" % (data_folder, sample_name)
    if isfile(file1):
        sample_df = pd.read_table(file1)
    else:
        file2 = "/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/SamplesForAnalysis_corrected/%s.xlsx" % (data_folder, sample_name)
        sample_df = pd.read_excel(file2)
    
    if newColumnList is not None:  # make sure column names are correct:
        sample_df = sample_df.iloc[:, :44]
        sample_df.columns = newColumnList
    sample_df = sample_df.rename(columns={'count (templates/reads)':'count (templates)', 'count (reads)':'count (templates)'})
    sample_df['prod_stat'] = np.where(sample_df['sequenceStatus'] == 'In', 1, 0)
    
    interesting_columns = ['vGeneName', 'vFamilyName', 'dFamilyName', 'jGeneName', 'jFamilyName', 'prod_stat']       

    sample_df = sample_df[interesting_columns]
    for column in sample_df.columns.values:  # remove 'unresolved' counts and calculate frequencies without it. remove
                                            # 'TCRB' indication
        if sample_df[column].dtype != 'int64':
            sample_df = sample_df[sample_df[column].str.contains('unresolved') == False]
            sample_df[column] = sample_df[column].str.replace('TCRB', '')
    
    prod_df = sample_df[sample_df['prod_stat'] == 1]
    nonprod_df = sample_df[sample_df['prod_stat'] == 0]
    
    df_dict = [('Total', sample_df),
     ('Prod', prod_df),
     ('nonProd', nonprod_df)]
 
 # (2) call selected functions to calculate statistics for each df
    for df in df_dict:
        count_geneUsage(df, data_folder, sample_name)
        count_dependent_geneUsage(df, data_folder, sample_name)
        
        
        
#-----------------------------------------------------------------------------------

### generate feature DFs###

'''the following functions should be used in sequence in order to generate seperate and merged feature summary df for
a cohort.
!!!MAKE sure that the datasetName matches the data_folder!!!!

example usage:

(1) seperate dfs for each seqType:
seqTypeList = ['Total', 'Prod','nonProd']
datasetName='PNP490'
data_folder='TCR_real_data'
dfs_folder='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/Total' %data_folder
FeatureGroups=listdir(dfs_folder)
for seqType in seqTypeList:
    print seqType
    gen_featureSummaryDF_forSeqType(seqType,data_folder,datasetName,FeatureGroups)
    
the resulting df will be in folder '\featureSummaryDFs'
    
(2) merge all features to one df:
datasetName='PNP490'
data_folder='TCR_real_data'
gen_merged_feature_df_for_dataset(data_folder,datasetName)

'''
        

def gen_featureSummaryDF_forSeqType(seqType, data_folder, datasetName, FeatureGroups):
    
    if seqType != 'Total':
        if 'SharingFeatures' in FeatureGroups:
            FeatureGroups.remove('SharingFeatures')
    AllSamplesDFLIST = []
    for group in FeatureGroups:
        print group
        dfs_folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/descriptiveStatsSamplesForAnalysis/%s/%s' % (data_folder, seqType, group)  # take the folder to concatentate sample dfs
        # remove unnecessary files if exist:
        filesToRemove = ['BD139b_05ug', 'BD202b_05ug', 'BD207b_05ug']
        for f in filesToRemove:
            if isfile('%s/%s' % (dfs_folder, f)):
                remove('%s/%s' % (dfs_folder, f))
                print 'file %s was removed from folder %s' % (f, dfs_folder)
        
        
        
        AllSamplesDF = concat_summarizing_dfs(dfs_folder)
        AllSamplesDFLIST.append(AllSamplesDF)
    
    for n, df in enumerate(AllSamplesDFLIST):
        if df.index[0] == 0:
            df = df.set_index('Sample')
            print df.columns.values[0]
        if n == 0:
            AllFeatures = df
        if n > 0:
            AllFeatures = pd.merge(AllFeatures, df, how='inner', left_index=True, right_index=True)

    print AllFeatures[['df type_x', 'df type_y']].iloc[0]  # check that all df type indications are the same
    print('n columns in merged df is %s' % len(AllFeatures.columns.values))
    columns_to_drop1 = [column for column in AllFeatures.columns.values if 'df type' in column]
    columns_to_drop2 = ['prod_stat_std', 'prod_stat_max']
    AllFeatures = AllFeatures.drop(columns_to_drop1, axis=1)
    AllFeatures = AllFeatures.drop(columns_to_drop2, axis=1)
    print('n columns in merged df after drop is %s' % len(AllFeatures.columns.values))
    
    # save to pickle and excel:
    folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/featureSummaryDFs' % data_folder
    if not isdir(folder):
            makedirs(folder)
    filePickle = '%s/%s_%s' % (folder, datasetName, seqType)
    filexlsx = '%s/%s_%s.xlsx' % (folder, datasetName, seqType)
    AllFeatures.to_pickle(filePickle)
    AllFeatures.to_excel(filexlsx)
    
    return AllFeatures
    
    

def gen_merged_feature_df_for_dataset(data_folder, datasetName):
    
    ##### compare column lists in all generated feature dfs:
    
    # total:
    Total_file = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/featureSummaryDFs/%s_Total' % (data_folder, datasetName)
    TotalFeatureDF = pd.read_pickle(Total_file)
    TotalColumns = TotalFeatureDF.columns.values
    print 'Total col number is %s' % len(TotalFeatureDF)
    
    
    # prod:
    Prod_file = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/featureSummaryDFs/%s_Prod' % (data_folder, datasetName)
    ProdFeatureDF = pd.read_pickle(Prod_file)
    ProdColumns = ProdFeatureDF.columns.values
    print 'Prod col number is %s' % len(ProdFeatureDF)
    
    # nonProd:
    nonProd_file = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/featureSummaryDFs/%s_nonProd' % (data_folder, datasetName)
    nonProdFeatureDF = pd.read_pickle(nonProd_file)
    nonProdColumns = nonProdFeatureDF.columns.values
    print 'nonProd col number is %s' % len(nonProdFeatureDF)
    
    # ## compare columns:
    
    TotalNotInProd = [column for column in TotalColumns if column not in ProdColumns ]
    print 'the following columns appear in total but not in prod:'
    print TotalNotInProd
    
    TotalNotInnonProd = [column for column in TotalColumns if column not in nonProdColumns ]
    print 'the following columns appear in total but not in nonProd:'
    print TotalNotInnonProd
    
    ProdNotInNonProd = [column for column in ProdColumns if column not in nonProdColumns ]
    print 'the following columns appear in prod but not in nonProd:'
    print ProdNotInNonProd
    
    nonProdNotInProd = [column for column in nonProdColumns if column not in ProdColumns]
    print 'the following columns appear in nonProd but not in Prod:'
    print nonProdNotInProd
    
    
    ##### merge the three feature summary DFs into one:
    TotalFeatureDF = TotalFeatureDF.add_suffix('_T')
    ProdFeatureDF = ProdFeatureDF.add_suffix('_1')
    nonProdFeatureDF = nonProdFeatureDF.add_suffix('_0')

    cohort_allFeatures = pd.merge(TotalFeatureDF, ProdFeatureDF, how='inner', left_index=True, right_index=True)
    cohort_allFeatures = pd.merge(cohort_allFeatures, nonProdFeatureDF, how='inner', left_index=True, right_index=True)
    print 'number of columns in merged feature file for cohort %s is %s' % (datasetName, len(cohort_allFeatures.columns.values))
    
    ##### save to pickle and excel:
    
    folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/featureSummaryDFs' % data_folder
    if not isdir(folder):
        makedirs(folder)
    filePickle = '%s/%s_allFeatures' % (folder, datasetName)
    filexlsx = '%s/%s_allFeatures.xlsx' % (folder, datasetName)
    cohort_allFeatures.to_pickle(filePickle)
    cohort_allFeatures.to_excel(filexlsx)
    
    ##### check which column have nan's? 
    expectedCount = len(cohort_allFeatures)
    print 'expected count is %s' % expectedCount
    print 'the following columns dont have the expected count:'
    for column in cohort_allFeatures.columns.values:
        if cohort_allFeatures[column].count() != expectedCount:
            print column
            print cohort_allFeatures[column].count()
        
    return cohort_allFeatures



#------------------------------------------------------------------------------------------------------

### dataset comparison functions####

'''
the first 3 functions are used in combination to generate an excel file summarizing all feature means in the two compared datasets, and ttest and ks test p-values
(need to be corrected for multiple comparisons), and a plot of the most interesting feature comparison.

usage example:
data_folder1='TCR_real_data'
data_folder2='TCR_real_data/CardioSamples'
datasetName1='PNP490'
datasetName2='cardio_plate8'
compare_features_between_datasets(data_folder1,datasetName1,data_folder2,datasetName2)

the forth function is used to generate the same excel file, and a plot of gene usage comparison

usage example:
data_folder1='TCR_real_data'
data_folder2='TCR_real_data/CardioSamples'
datasetName1='PNP490'
datasetName2='cardio_plate8'


plot_gene_usage_comparison(data_folder1,data_folder2,datasetName1,datasetName2)



'''
def plot_feature_comparison(feature, ax, datasetName1, cohort1_allFeatures, datasetName2, cohort2_allFeatures,filteringList1Name=None,
                            filteringList2Name=None):
    cohort1data = cohort1_allFeatures[feature]
    cohort1data = list(cohort1data[~np.isnan(cohort1data)])
    cohort1weights = np.ones_like(cohort1data, dtype=np.float) / len(cohort1data)
    cohort2data = cohort2_allFeatures[feature]
    cohort2data = list(cohort2data[~np.isnan(cohort2data)])
    cohort2weights = np.ones_like(cohort2data, dtype=np.float) / len(cohort2data)  
    
    FeatureMeans_dataset1 = np.mean(cohort1data)
    FeatureMeans_dataset2 = np.mean(cohort2data)
    

    fig = plt.figure()
    if ax == None:
        fig, ax = plt.subplots()

    plot = ax.hist((cohort1data, cohort2data), bins=50, color=('black', 'green'), weights=[cohort1weights, cohort2weights],
                 label=('%s_%s' %(datasetName1,filteringList1Name), '%s_%s' %(datasetName2,filteringList2Name)), alpha=0.7)
    
    ks_p_cohort1_cohort2, t_p_cohort1_cohort2 = compare_feature_statistics(feature, datasetName1, cohort1_allFeatures,
                                                                            datasetName2, cohort2_allFeatures)


    ax.annotate('KS_p=%s\nttest_p=%s\n%s_%s mean=%s\n%s_%s mean=%s' % (round(ks_p_cohort1_cohort2, 6), round(t_p_cohort1_cohort2, 6),
                    datasetName1, filteringList1Name, round(FeatureMeans_dataset1, 3), datasetName2, filteringList2Name, round(FeatureMeans_dataset2, 3)),
                    xy=(0.96, 0.95), xycoords='axes fraction', fontsize=9, horizontalalignment='right', verticalalignment='top', fontweight='bold')

    ax.set_title(feature, fontsize=16, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=9)
    
    return fig, ax


def compare_feature_statistics(feature, datasetName1, cohort1_allFeatures, datasetName2, cohort2_allFeatures):
    cohort1data = cohort1_allFeatures[feature]
    cohort1data = list(cohort1data[~np.isnan(cohort1data)])
    cohort1weights = np.ones_like(cohort1data) / len(cohort1data)
    cohort2data = cohort2_allFeatures[feature]
    cohort2data = list(cohort2data[~np.isnan(cohort2data)])
    cohort2weights = np.ones_like(cohort2data) / len(cohort2data)
    
    try:
        ks_s_cohort1_cohort2, ks_p_cohort1_cohort2 = stats.ks_2samp(cohort1data, cohort2data)
    except ValueError:
        print 'couldnt execute ks test'
        ks_s_cohort1_cohort2 = np.nan
        ks_p_cohort1_cohort2 = np.nan
    
    try:
        t_s_cohort1_cohort2, t_p_cohort1_cohort2 = stats.ttest_ind(cohort1data, cohort2data)
    except ValueError:
        print 'couldnt execute t test'
        t_s_cohort1_cohort2 = np.nan
        t_p_cohort1_cohort2 = np.nan
    
    return ks_p_cohort1_cohort2, t_p_cohort1_cohort2

def compare_features_between_datasets(data_folder1, datasetName1,  data_folder2, datasetName2,TakeSameSamples=False, 
                                      filteringList1=None,filteringList2=None,filteringList1Name=None,
                                      filteringList2Name=None,MatchedFolderToSave=None,samp_size_nt=None,samp_size_aa=None):
    
    ##### generate excel file with mean values for the compared datasets:
    print 'generating excel file with mean values for compared dataset (find in -featureSummaryDFs folder-)'
    folder1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/featureSummaryDFs' % data_folder1
    filePickle1 = '%s/%s_allFeatures' % (folder1, datasetName1)
    cohort1_allFeatures = pd.read_pickle(filePickle1)
    cohort1_allFeatures = editSampleNames(cohort1_allFeatures)
    
    if samp_size_nt is None:
        samp_size_nt = 2000
    if samp_size_aa is None:
        samp_size_aa = 200
    
    if filteringList1 is not None:
        cohort1_allFeatures=cohort1_allFeatures.loc[filteringList1,:]
        print 'filtering cohort1 to include %s samples' %len(filteringList1)
    else:
        filteringList1Name=''
    
    
    folder2 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/featureSummaryDFs' % data_folder2
    filePickle2 = '%s/%s_allFeatures' % (folder2, datasetName2)
    cohort2_allFeatures = pd.read_pickle(filePickle2)
    cohort2_allFeatures = editSampleNames(cohort2_allFeatures)
    
    if filteringList2 is not None:
        cohort2_allFeatures=cohort2_allFeatures.loc[filteringList2,:]
        print 'filtering cohort2 to include %s samples' %len(filteringList2)
    else:
        filteringList2Name=''
    
    
    # take identical samples from each dataset:
    if TakeSameSamples:
        print 'taking only samples common to the two cohorts:'
        cohort1_sampleList = list(cohort1_allFeatures.index)
        cohort2_sampleList = list(cohort2_allFeatures.index)
#         cohort1_sampleList=editSampleNamesList(list(cohort1_allFeatures.index))
#         cohort2_sampleList=editSampleNamesList(list(cohort2_allFeatures.index))
        CommonSamples = list(set(cohort1_sampleList).intersection(set(cohort2_sampleList)))
        print 'number of samples in cohort 1 is %s' % len(cohort1_sampleList)
        print 'number of samples in cohort 2 is %s' % len(cohort2_sampleList)
        print 'number of common samples is %s' % len(CommonSamples)
        cohort1_allFeatures = cohort1_allFeatures.loc[CommonSamples, :]  # take only the common samples
        cohort2_allFeatures = cohort2_allFeatures.loc[CommonSamples, :]  # take only the common samples
        print ' cohort1 df length is now %s' % len(cohort1_allFeatures)
        print ' cohort2 df length is now %s' % len(cohort2_allFeatures)
    
    
    FeatureMeans_dataset1 = pd.DataFrame(cohort1_allFeatures.mean())
    FeatureMeans_dataset2 = pd.DataFrame(cohort2_allFeatures.mean())
    
    FeatureMeanSummary_all_dataSets = pd.merge(FeatureMeans_dataset1, FeatureMeans_dataset2, how='outer', left_index=True, right_index=True)
    
    
    FeatureMeanSummary_all_dataSets = FeatureMeanSummary_all_dataSets.rename(columns={'0_x':'%s_%s' %(datasetName1,filteringList1Name),
                                                                               '0_y':'%s_%s' %(datasetName2,filteringList2Name)})
    # compare datasets and plot selected features:
    print 'comparing datasets and plot selected features:...'
    fig1 = plt.figure(figsize=(12, 20))
    fig1.suptitle('Feature Comparison - %s vs. %s' % ('%s_%s' %(datasetName1,filteringList1Name), '%s_%s' %(datasetName2,filteringList2Name)), fontsize=20)
    
    
    GeneralFeaturesToCompare = ['prod_stat_mean_T', 'max_nt_per_aa_T', 'mean_nt_per_aa_T',
                          'normSeqNums_per%s_AA_T' %samp_size_aa, 'normSeqNums_per%s_NT_T' %samp_size_nt]
    LengthFeaturesToCompare = ['cdr3Length_mean_T', 'n1Insertion_mean_T', 'n2Insertion_mean_T', 'vDeletion_mean_T', 'd3Deletion_mean_T', 'd5Deletion_mean_T',
                            'jDeletion_mean_T']
    ClonalityFeaturesToCompare = ['top10clonal_nt_T','top1000clonal_nt_T','frequencyCount (%)_mean_T',
                               'frequencyCount (%)_max_T', 'shannon_nt_T']
    
    
    GeneralCount = 0
    LengthCount = 0
    ClonalityCount = 0
    for n, feature in enumerate(FeatureMeanSummary_all_dataSets.index):
        if (feature in FeatureMeans_dataset1.index) and (feature in FeatureMeans_dataset2.index):
            
            
            ks_p_cohort1_cohort2, t_p_cohort1_cohort2 = compare_feature_statistics(feature, datasetName1, cohort1_allFeatures,
                                                                                datasetName2, cohort2_allFeatures)
            FeatureMeanSummary_all_dataSets.loc[feature, 'ks_p'] = ks_p_cohort1_cohort2
            FeatureMeanSummary_all_dataSets.loc[feature, 't_p'] = t_p_cohort1_cohort2
            
            FeatureMean_dataset1 = FeatureMeanSummary_all_dataSets.loc[feature, '%s_%s' %(datasetName1,filteringList1Name)]
            FeatureMean_dataset2 = FeatureMeanSummary_all_dataSets.loc[feature, '%s_%s' %(datasetName2,filteringList2Name)]

            
            if feature in GeneralFeaturesToCompare:
                ax = fig1.add_subplot(7, 3, 3 * GeneralCount + 1)
                plot_feature_comparison(feature, ax, datasetName1, cohort1_allFeatures, datasetName2,
                                        cohort2_allFeatures,filteringList1Name,filteringList2Name)
                if GeneralCount == 4:
                    ax.legend(bbox_to_anchor=(0.05, -2.85), loc='lower left', borderaxespad=0., fontsize=18)
                GeneralCount = GeneralCount + 1
            if feature in LengthFeaturesToCompare:
                ax = fig1.add_subplot(7, 3, 3 * LengthCount + 2)
                plot_feature_comparison(feature, ax, datasetName1, cohort1_allFeatures, datasetName2,
                                        cohort2_allFeatures,filteringList1Name,filteringList2Name)
                LengthCount = LengthCount + 1
            if feature in ClonalityFeaturesToCompare:
                ax = fig1.add_subplot(7, 3, 3 * ClonalityCount + 3)
                plot_feature_comparison(feature, ax, datasetName1, cohort1_allFeatures, datasetName2,
                                        cohort2_allFeatures,filteringList1Name,filteringList2Name)
                ClonalityCount = ClonalityCount + 1
    
    #ADD FDRs to results:
    FDR=0.1
    nTests_ks=len(FeatureMeanSummary_all_dataSets[FeatureMeanSummary_all_dataSets['ks_p'].notnull()])
    FeatureMeanSummary_all_dataSets=add_corrected_pValues(FeatureMeanSummary_all_dataSets,pValueColumn='ks_p',nTests=nTests_ks,FDR=FDR)
    FeatureMeanSummary_all_dataSets=FeatureMeanSummary_all_dataSets.rename(columns={'Sig by bonferroni corrected pVal':'Sig by bonferroni corrected pVal_ks',
                                                                            'sig. by FDR=%s' %FDR : 'sig. by FDR=%s_ks' %FDR })
    nTests_t=len(FeatureMeanSummary_all_dataSets[FeatureMeanSummary_all_dataSets['t_p'].notnull()])
    FeatureMeanSummary_all_dataSets=add_corrected_pValues(FeatureMeanSummary_all_dataSets,pValueColumn='t_p',nTests=nTests_ks,FDR=0.1)
    FeatureMeanSummary_all_dataSets=FeatureMeanSummary_all_dataSets.rename(columns={'Sig by bonferroni corrected pVal':'Sig by bonferroni corrected pVal_t',
                                                                            'sig. by FDR=%s' %FDR : 'sig. by FDR=%s_t' %FDR })
    
    
    fig1.subplots_adjust(left=0.09, right=0.98, top=0.95, bottom=0.02, wspace=0.22, hspace=0.28)
        
    
    ##### saving excel file and plot:
    if MatchedFolderToSave is not None:
        folder=MatchedFolderToSave
    else:
        folder = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/featureSummaryDFs' 
    if not isdir(folder):
        makedirs(folder)
    
    filePickle = '%s/FeatureMeanSummary_%s%s_%s%s' % (folder, datasetName1, filteringList1Name, datasetName2,filteringList2Name)
    filexlsx = '%s/FeatureMeanSummary_%s%s_%s%s.xlsx' % (folder, datasetName1, filteringList1Name, datasetName2,filteringList2Name)
    FeatureMeanSummary_all_dataSets.to_pickle(filePickle)
    FeatureMeanSummary_all_dataSets.to_excel(filexlsx)
    
    filename = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/realAnalysis/Feature_Comparison_%s%s_%s%s.png' % (datasetName1, filteringList1Name, datasetName2,filteringList2Name)
    fig1.savefig(filename, bbox_inches='tight', dpi=200) 
    print 'plot file can be found in -realAnalysis folder'   

    plt.show()
    return FeatureMeanSummary_all_dataSets
    

def plot_gene_usage_comparison(data_folder1, datasetName1,  data_folder2, datasetName2,plotType='bar',TakeSameSamples=False, 
                                      filteringList1=None,filteringList2=None,filteringList1Name=None,filteringList2Name=None):
    
    ##### generate excel file with mean values for the compared datasets:
    print 'generating excel file with mean values for compared dataset (find in -featureSummaryDFs folder-)'
    folder1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/featureSummaryDFs' % data_folder1
    filePickle1 = '%s/%s_allFeatures' % (folder1, datasetName1)
    cohort1_allFeatures = pd.read_pickle(filePickle1)
    cohort1_allFeatures = editSampleNames(cohort1_allFeatures)
    
    if filteringList1 is not None:
        cohort1_allFeatures=cohort1_allFeatures.loc[filteringList1,:]
        print 'filtering cohort1 to include %s samples' %len(filteringList1)
    else:
        filteringList1Name=''
    
    
    folder2 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/featureSummaryDFs' % data_folder2
    filePickle2 = '%s/%s_allFeatures' % (folder2, datasetName2)
    cohort2_allFeatures = pd.read_pickle(filePickle2)
    cohort2_allFeatures = editSampleNames(cohort2_allFeatures)
    
    if filteringList2 is not None:
        cohort2_allFeatures=cohort2_allFeatures.loc[filteringList2,:]
        print 'filtering cohort2 to include %s samples' %len(filteringList2)
    else:
        filteringList2Name=''
    
    if TakeSameSamples:
        print 'taking only samples common to the two cohorts:'
        cohort1_sampleList = list(cohort1_allFeatures.index)
        cohort2_sampleList = list(cohort2_allFeatures.index)
#         cohort1_sampleList=editSampleNamesList(list(cohort1_allFeatures.index))
#         cohort2_sampleList=editSampleNamesList(list(cohort2_allFeatures.index))
        CommonSamples = list(set(cohort1_sampleList).intersection(set(cohort2_sampleList)))
        print 'number of samples in cohort 1 is %s' % len(cohort1_sampleList)
        print 'number of samples in cohort 2 is %s' % len(cohort2_sampleList)
        print 'number of common samples is %s' % len(CommonSamples)
        cohort1_allFeatures = cohort1_allFeatures.loc[CommonSamples, :]  # take only the common samples
        cohort2_allFeatures = cohort2_allFeatures.loc[CommonSamples, :]  # take only the common samples
        print ' cohort1 df length is now %s' % len(cohort1_allFeatures)
        print ' cohort2 df length is now %s' % len(cohort2_allFeatures)
        
    FeatureMeans_dataset1 = pd.DataFrame(cohort1_allFeatures.mean())
    FeatureMeans_dataset2 = pd.DataFrame(cohort2_allFeatures.mean())

    
    FeatureMeanSummary_all_dataSets = pd.merge(FeatureMeans_dataset1, FeatureMeans_dataset2, how='outer', left_index=True, right_index=True)   
    FeatureMeanSummary_all_dataSets = FeatureMeanSummary_all_dataSets.rename(columns={'0_x':'%s_%s' % (datasetName1,filteringList1Name),
                                                                               '0_y':'%s_%s' % (datasetName2,filteringList2Name)})
#     print 'FeatureMeanSummary:'
#     print FeatureMeanSummary_all_dataSets.head()
    ###### define sub-dfs for each gene usage type
    regex_vFamily = re.compile('V.._T')
    indices_vFamily = [n for n in FeatureMeanSummary_all_dataSets.index if re.match(regex_vFamily, n) ]
    vFamily = FeatureMeanSummary_all_dataSets.loc[indices_vFamily, ['%s_%s' %(datasetName1,filteringList1Name), '%s_%s' %(datasetName2,filteringList2Name)]]

    regex_vGene = re.compile('V..-.._T')
    indices_vGene = [n for n in FeatureMeanSummary_all_dataSets.index if re.match(regex_vGene, n) ]
    vGene = FeatureMeanSummary_all_dataSets.loc[indices_vGene, ['%s_%s' %(datasetName1,filteringList1Name), '%s_%s' %(datasetName2,filteringList2Name)]]

    regex_jGene = re.compile('J..-.._T')
    indices_jGene = [n for n in FeatureMeanSummary_all_dataSets.index if re.match(regex_jGene, n) ]
    jGene = FeatureMeanSummary_all_dataSets.loc[indices_jGene, ['%s_%s' %(datasetName1,filteringList1Name), '%s_%s' %(datasetName2,filteringList2Name)]]

    regex_dFamily = re.compile('D.._T')
    indices_dFamily = [n for n in FeatureMeanSummary_all_dataSets.index if re.match(regex_dFamily, n) ]
    dFamily = FeatureMeanSummary_all_dataSets.loc[indices_dFamily, ['%s_%s' %(datasetName1,filteringList1Name), '%s_%s' %(datasetName2,filteringList2Name)]]

    regex_VJ = re.compile('V.._J.._T')
    indices_VJ = [n for n in FeatureMeanSummary_all_dataSets.index if re.match(regex_VJ, n) ]
    VJ = FeatureMeanSummary_all_dataSets.loc[indices_VJ, ['%s_%s' %(datasetName1,filteringList1Name), '%s_%s' %(datasetName2,filteringList2Name)]]

    regex_DJ = re.compile('D.._J..-.._T')
    indices_DJ = [n for n in FeatureMeanSummary_all_dataSets.index if re.match(regex_DJ, n) ]
    DJ = FeatureMeanSummary_all_dataSets.loc[indices_DJ, ['%s_%s' %(datasetName1,filteringList1Name), '%s_%s' %(datasetName2,filteringList2Name)]]

    
    
    #### plot gene usage:
    print 'plotting gene usage comparison...'
    usageFeatureList = [vFamily, vGene, jGene, dFamily, VJ, DJ]
    usageFeatureNameList = ['vFamily', 'vGene', 'jGene', 'dFamily', 'VJ', 'DJ']
    colspanList = [3, 3, 2, 1, 3, 3]

    fig1 = plt.figure(figsize=(10, 20))
    fig1.suptitle('Gene Usage Comparison - %s_%s vs. %s_%s' % (datasetName1,filteringList1Name,datasetName2,filteringList2Name), fontsize=20)

    for n, feature in enumerate(usageFeatureList):
#         print feature
        colspan = colspanList[n]
        ax1 = plt.subplot2grid((6, 3), (n, 0), colspan=colspan)
    #     ax= fig1.add_subplot(6,1,n+1)
#         if 'box' in plotType:
#             print 'plotting box plot'
#             feature.plot.box(ax=ax1, label=('%s' % datasetName1, '%s' % datasetName2))
#         else:
        feature.plot.bar(ax=ax1, color=('black', 'green'), label=('%s_%s' % (datasetName1,filteringList1Name), '%s_%s' % (datasetName2,filteringList2Name)), alpha=0.7)
        
        ax1.set_title(usageFeatureNameList[n], fontsize=16, fontweight='bold')
        try:
            if n == 0:
                ax1.legend(bbox_to_anchor=(1.01, 0.95), loc='upper left', borderaxespad=0., fontsize=16)
            else:
                ax1.legend().set_visible(False)
        except:
            pass

    fig1.subplots_adjust(left=0.09, right=0.98, top=0.94, bottom=0.01, wspace=0.22, hspace=0.54)

    filename = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/realAnalysis/geneUsage_Comparison_%s%s_%s%s' % (datasetName1, filteringList1Name, datasetName2,filteringList2Name)
    fig1.savefig(filename, bbox_inches='tight', dpi=200)  
    print 'plot file can be found in -realAnalysis folder'

    plt.show()
    
    
    
#-----------------------------------------------------------------------------------------------------------------
'''
the following function compares datasets for main phenotypes: age, gender, hemoglobin, BMI, CRP, eversmoked
it was copied from the notebook 'Generate features DF'

input:
sampleList1 and 2 - to define for which samples to take the phenotypic data from the phtnotypeDF.
datasetName1 and 2

THIS FUNCTION APPLIES ONLY FOR COMPARING SAMPLES INCLUDED IN THE SAME PHENOTYPIC DATA FILE!
NEED TO ADAPT FOR OTHER SITUATIONS.

usage example:
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP515','rb') as fp:
    PNP515=pickle.load(fp)
    
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP434samples','rb') as fp:
    PNP434samples=pickle.load(fp)
    
file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/fullXgroupbyBD'
fullXgroupbyBD=pd.read_pickle(file1)


sampleList1=PNP434samples
sampleList2=PNP515
phenotypeDF=fullXgroupbyBD
datasetName1='PNPsamples_434'
datasetName2='PNP515'

fig1=compare_mainPhenotypes_between_datasets(sampleList1,datasetName1,sampleList2,datasetName2,phenotypeDF)


'''
    
    
def roundup2(a, digits=0):
    n = 10 ** -digits
    return round(math.ceil(a / n) * n, digits)


def compare_mainPhenotypes_between_datasets(sampleList1, datasetName1, sampleList2, datasetName2, phenotypeDF,filteringList1Name='',
                                            filteringList2Name=''):
    
    numericalFeatures = ['Age', 'Hemoglobin', 'BMI', 'CRP (WIDE RANGE)']
    binaryFeatures = ['Gender', 'Ever smoked']

    fig1 = plt.figure(figsize=(6, 18))
    fig1.suptitle('Main Phenotype Distributions', fontsize=22)

    # plotting numerical phenotypes:
       
    for n, feature in enumerate(numericalFeatures):
        print n, feature
        data1 = phenotypeDF.loc[sampleList1, feature]
        data1 = data1[data1.notnull()]
        data1 = list(data1)
        weights1 = np.ones_like(data1) / len(data1)
        mean1 = round(np.mean(data1), 2)
        std1 = round(np.std(data1), 2)
    #     print len(data1)
        data2 = phenotypeDF.loc[sampleList2, feature]
        data2 = data2[data2.notnull()]
        data2 = list(data2)
        weights2 = np.ones_like(data2) / len(data2)
        mean2 = round(np.mean(data2), 2)
        std2 = round(np.std(data2), 2)
        Alldata = data1 + data2
        Allweights = np.ones_like(Alldata) / len(Alldata)

        ax = fig1.add_subplot(5, 1, n + 1)
        plot = ax.hist((data1, data2), bins=50, color=('black', 'green'), weights=[weights1, weights2],
                     label=('%s_%s' %(datasetName1,filteringList1Name), '%s_%s' %(datasetName2,filteringList2Name)), alpha=0.7)


        ks_s_cohort1_cohort2, ks_p_cohort1_cohort2 = stats.ks_2samp(data1, data2)
        t_s_cohort1_cohort2, t_p_cohort1_cohort2 = stats.ttest_ind(data1, data2)

        ax.annotate('KS_p=%s\nttest_p=%s\n%s mean=%s\n%s mean=%s' % (round(ks_p_cohort1_cohort2, 6), round(t_p_cohort1_cohort2, 6),
                        datasetName1, round(mean1, 3), datasetName2, round(mean2, 3)),
                        xy=(0.96, 0.95), xycoords='axes fraction', fontsize=9, horizontalalignment='right', verticalalignment='top', fontweight='bold')

        ax.set_title(feature, fontsize=16, fontweight='bold')
    #     ax.set_ylabel('Frequency',fontsize=9)
        if n == 0:
            ax.legend(bbox_to_anchor=(1.01, 0.95), loc='upper left', borderaxespad=0., fontsize=16)
        else:
            ax.legend().set_visible(False)
    
    
    # plotting binary phenotypes:
    
    for n, Bfeature in enumerate(binaryFeatures):
        print n, Bfeature
        if Bfeature == 'Gender':
            ticklabels = ['Male', 'Female']
        elif Bfeature == 'Ever smoked':
            ticklabels = ['No', 'Yes']

        ax = fig1.add_subplot(5, 2, n + 9)   
        a1 = phenotypeDF.loc[sampleList1, Bfeature].value_counts(normalize=True)
        a2 = phenotypeDF.loc[sampleList2, Bfeature].value_counts(normalize=True)
        combined = list(a1) + list(a2)
        Ymax = roundup2(np.max(combined), 1)
        ax.bar([-0.15, 0.85], a1, align='center', width=0.3, tick_label=ticklabels, color='black')
        ax.bar([0.15, 1.15], a2, align='center', width=0.3, tick_label=ticklabels, color='green')
        ax.set_xlabel(Bfeature, fontsize=16, fontweight='bold')
        ax.set_ylim(0, Ymax)
        
        
    fig1.text(0, 0.5, "Frequency", ha='left', fontsize=18, rotation=90)
    fig1.subplots_adjust(left=0.13, right=0.98, top=0.92, bottom=0.02, wspace=0.22, hspace=0.25)
    filename = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/realAnalysis/MainPhenotypecomparison_%s_%s__%s_%s' \
    % (datasetName1,filteringList1Name,datasetName2,filteringList2Name)
    fig1.savefig(filename, bbox_inches='tight', dpi=200) 
    
    print 'plot file can be found in -real analysis- folder'
    
    plt.show()
    return fig1


#-----------------------------------------------------------------------------------------------------
'''
the following function plots a cumulative distribution of number of templates per samples in a dataset, and annotate the exact
number of samples for percentiles 1,5,15

this function was copied from notebook 'PopulationAnalysis_new_version', step 2.1

input:
nTemplates- a df generated by the nTemplates calulcation in the notebook. see usage example
datasetName-for saving the figure
data_folder - for saving the figure

usage example:
file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/nTemplatesSummary_515_08042018' 
nTemplates=pd.read_pickle(file1)

nTemplates=nTemplates
datasetName='PNP515'
data_folder='TCR_real_data'

fig=plot_nTemplates_percentile(nTemplates, datasetName,data_folder)
plt.show()

'''

def plot_nTemplates_percentile(nTemplates, datasetName, data_folder):
    templates = list(nTemplates['n_templates_total_nt'])
    ymax = roundup(np.max(templates), 10000)
    fig, ax = plt.subplots(figsize=(12, 9))
    perc = [0, 1, 2, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 95, 99, 100]
    temps = []
    for n in perc:
        print n, np.percentile(templates, n)
        temps.append(np.percentile(templates, n))
    ax.plot(perc, temps)
    percToShow = [1, 5, 15]
    tempsToShow = []
    for p in percToShow:
        ind = perc.index(p)
        t = temps[ind]
        ax.scatter(p, t, color='r', s=20)
        ax.annotate('%sperc-%s templates' % (p, t), xy=(p, t),
                             xytext=(p * 1.05, t * 1.05), fontsize=11, color='r')
    ax.set_xlabel('percent of samples', fontsize=20)
    ax.set_ylabel('number of templates', fontsize=20)
    ax.set_title('Cumulative distribution of #templates per sample - %s' % datasetName, fontsize=26)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, ymax)
    
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/%s/realAnalysis/CumDistNtemplates_%s' % (data_folder, datasetName)
    fig.savefig(file1, dpi=200)
    
    return fig


       
    

#----------------------------------------

def complete_TCRfeatureExtractionForDataset(datasetFolder,datasetName, ss=None,
                    repeat=None,samp_size_nt=None,samp_size_aa=None):

    SampleFolder='%s/SamplesForAnalysis_corrected' %datasetFolder
    partial_datasetFolder=datasetFolder.replace('%s/' %MyPath,'')
    print partial_datasetFolder
    
    # (1) feature extraction:
    print 'step 1: feature extraction (long)'

    filenames = [f for f in listdir(SampleFolder) if isfile(join(SampleFolder, f))]
    filenames=[f.strip('.tsv') for f in filenames]
    filenames=[f.strip('.xlsx') for f in filenames]
    print 'number of samples for feature extraction is %s' %len(filenames)

    newColumnList=None
    for n,sample_name in enumerate(filenames): 
    #     if n>91: 
            print  n,sample_name
            if 'nSampled' not in sample_name:
                gen_descriptive_stats(sample_name,partial_datasetFolder,newColumnList,samp_size_nt,samp_size_aa)
                gen_geneUsageCount(sample_name,partial_datasetFolder,newColumnList)


    # (3) generate feature summary DFs:

    print 'generating seperate feature data dfs: (long)'

    seqTypeList = ['Total', 'Prod','nonProd']
    if ss is not None:
        newDatasetName='%s_ss%s_rep%s' %(datasetName,ss,repeat)
    else:
        newDatasetName=datasetName
    
    FeatureFolder='%s/descriptiveStatsSamplesForAnalysis/Total' %(datasetFolder)
    FeatureGroups=listdir(FeatureFolder)
    for seqType in seqTypeList:
        print seqType
        gen_featureSummaryDF_forSeqType(seqType,partial_datasetFolder,newDatasetName,FeatureGroups)

    print 'generating merged feature data df: (short)'
    gen_merged_feature_df_for_dataset(partial_datasetFolder,newDatasetName)


#-----------------------------------------------------------------------------------------------

def remove_correlated_features(TCRfeatureDF,corrThreshold):
    # # ** remove columns with too many nans
    print 'TCR feature DF shape after loading and before any removal  is %s_%s' %(TCRfeatureDF.shape[0],TCRfeatureDF.shape[1])
    print 'removing feature columns with less than 5 samples that are not nan...'
    nanColumns=TCRfeatureDF.loc[:, TCRfeatureDF.notnull().sum() <=5].columns.tolist()
    print nanColumns
    print 'number of such columns is %s' %len(TCRfeatureDF.loc[:, TCRfeatureDF.notnull().sum() <=5].columns)
    TCRfeatureDF=TCRfeatureDF.loc[:, TCRfeatureDF.loc[:, TCRfeatureDF.notnull().sum() > 5].columns.tolist()]
    print 'TCRfeatureDF shape after removing all columns that dont have at least 5 valid values is %s_%s'\
    %(TCRfeatureDF.shape[0],  TCRfeatureDF.shape[1])

    #generating list of columns to drop with correlation higher or equal to threshold                                                                                                        
    print 'TCRfeatureDF2 shape before dropping highly correlated columns=%s_%s' % (TCRfeatureDF.shape[0],                                                                             TCRfeatureDF.shape[1])

    TCRfeatureDF2=TCRfeatureDF.copy()
    print 'filtering out columns with correlation larger than %s' % corrThreshold
    colToDrop = []
    for i in range(len(TCRfeatureDF2.columns)):
        for j in range(i + 1, len(TCRfeatureDF2.columns)):
            col1 = TCRfeatureDF2.columns.values[i]
            col2 = TCRfeatureDF2.columns.values[j]
            try:
                r, p = MyPearsonr(TCRfeatureDF2[col1], TCRfeatureDF2[col2])
                if r >= corrThreshold:
                    print col1, col2, r
                    dropCol=col2
                    if ('_T' not in col1)& ('_T' in col2):
                        dropCol=col1
                    else:
                        if ('_0' in col1)& ('_0' not in col2):
                            dropCol=col1
                        else:
                            if ('_std' in col1)& ('_std' not in col2):
                                dropCol=col1
                            else:
                                if ('-' in col1)& ('-' not in col2):
                                    dropCol=col1
    #                 print col1, col2, r, len(TCRfeatureDF2[TCRfeatureDF2[col1].notnull()]),len(TCRfeatureDF2[TCRfeatureDF2[col2].notnull()])
                    colToDrop.append(dropCol)
            except:
                print 'couldnt calc correlation for %s and %s' %(col1,col2)
                continue
                
    TCRfeatureDF3 = TCRfeatureDF2.drop(colToDrop, axis=1)
    print 'TCRfeatureDF3 shape after dropping highly correlated columns=%s_%s' % (TCRfeatureDF3.shape[0], TCRfeatureDF3.shape[1])
    # TCRfeatureDF2Name = '%s_fCorr%s' % (TCRfeatureDF2Name, corrThreshold)
    
    colToDrop2=list(set(colToDrop))
    colToDrop3=colToDrop2+nanColumns
    print 'colToDrop list length with repeats: %s' %len(colToDrop)
    print 'colToDrop list length without repeats: %s' %len(colToDrop2)
    print 'colToDrop list length without repeats and with nan columns: %s' %len(colToDrop3)
    print colToDrop3
    
    return TCRfeatureDF3,colToDrop3



def explore_feature_table(featureDF):
    df=pd.DataFrame()
    for n,col in enumerate(featureDF.columns):
        nans=len(featureDF[featureDF[col].isnull()])
        n9999=len(featureDF[featureDF[col]==9999])
        VCnorm=featureDF[col].value_counts(normalize=True).sort_values(ascending=False)
#         print VCnorm
        nCategs=len(VCnorm)
        topValueFrac=VCnorm.iloc[0]
        topValue=VCnorm.index[0]
    #     print col,nans,n9999,nCategs
        df.loc[n,'number of Nans']=nans
        df.loc[n,'number of 9999s']=n9999
        df.loc[n,'number of categories']=nCategs
        df.loc[n,'columns']=col
        df.loc[n,'topValue']=topValue
        df.loc[n,'topValueFrac']=topValueFrac
        df.loc[n,'type']=featureDF[col].dtype
    df['total nans']=df['number of Nans']+df['number of 9999s']
    df['Type']=df['number of categories'].apply(lambda x: 'constant' if x<2 else ('binary' if x==2 else ('categorial' if (x>2)&\
    (x<6) else 'Continuous')))
    
    return df

#------------------------------------------------------------------------------------------------------

def gen_featureDF_for_MatchedCohorts(ss,repeat,sampleList,sampleListName):
    print 'getting TCRfeature DF...'
    if ss is None:
        PNPfeatureFile='%s/TCR_real_data/featureSummaryDFs/PNP530_allFeatures' %MyPath
    else:
        PNPfeatureFile='%s/TCR_real_data/PNP530_SubSampled%sdata_rep%s/featureSummaryDFs/\
PNP530_ss%s_rep%s_allFeatures' %(MyPath,ss,repeat,ss,repeat)
    PNPfeatureDF=pd.read_pickle(PNPfeatureFile)
    if ss is None:
        CardiofeatureFile='%s/TCR_real_data/CardioSamples/featureSummaryDFs/Cardio126_allFeatures' %MyPath
    else:
        CardiofeatureFile='%s/TCR_real_data/CardioSamples/Cardio126_SubSampled%sdata_rep%s/featureSummaryDFs/\
Cardio126_ss%s_rep%s_allFeatures' %(MyPath,ss,repeat,ss,repeat)
    CardiofeatureDF=pd.read_pickle(CardiofeatureFile)
    
    if ss is None:
        datasetFolder='%s/TCR_real_data/PNP530Cardio126Combined' %MyPath
        datasetName='PNP530Cardio126'
    else:
        datasetFolder='%s/TCR_real_data/PNP530Cardio126Combined/MatchedSamples/ss%srep%s/' %(MyPath,ss,repeat)
        datasetName='MatchedSamples_ss%srep%s' %(ss,repeat)

    TCRfeatureDF=pd.concat([PNPfeatureDF,CardiofeatureDF])
    TCRfeatureDF=editSampleNames(TCRfeatureDF)
    # TCRfeatureDF.head()
    
    if sampleListName is None:
        featureFileXls='%s/featureSummaryDFs/%s_allFeatures.xlsx' %(datasetFolder,datasetName)
        featureFilePickle='%s/featureSummaryDFs/%s_allFeatures' %(datasetFolder,datasetName)
    else:
        featureFileXls='%s/featureSummaryDFs/%s_%s_allFeatures.xlsx' %(datasetFolder,sampleListName,datasetName)
        featureFilePickle='%s/featureSummaryDFs/%s_%s_allFeatures' %(datasetFolder,sampleListName,datasetName)
        TCRfeatureDF_matched=TCRfeatureDF.loc[sampleList,:]
    print 'TCR feature DF shape is %s_%s' %(TCRfeatureDF_matched.shape[0],TCRfeatureDF_matched.shape[1])
    TCRfeatureDF_matched.to_excel(featureFileXls)
    TCRfeatureDF_matched.to_pickle(featureFilePickle)
        
    return TCRfeatureDF_matched

    
    
#----------------------------------------------------------------------------
'''
the following function takes a feature DF and clean it - delete features with mroe than 95% nans, delete features which are completely
correlated with another features, explore featuredf after this cleaning, removes constnat features and fillna values for all gene
usage features

'''
def process_featureDF(featureDF,corrThreshold=1,removeCorr=True,exploreDF=True,deleteConsts=True,fillnans=True):
    featureDFnameAddition=''
    if removeCorr:
        print 'removing all columns correlated to other column by r=%s or more' %corrThreshold
        featureDF,colToDrop=remove_correlated_features(featureDF,corrThreshold) #remove nans and highly correlated
        featureDFnameAddition=featureDFnameAddition+'_noCorr'+str(corrThreshold)
        featureDFnameAddition=featureDFnameAddition.replace('.','-')
    if fillnans:
        print 'filling nans with 0...'
        featureDFnameAddition=featureDFnameAddition+'_nanFilled'
        colsTofill=[x for x in featureDF.columns.values if 'norm' not in x]
        print 'number of columns to fill nans=%s' %len(colsTofill)
        print 'total numvber of columns=%s' %len(featureDF.columns)
        featureDF[colsTofill]=featureDF[colsTofill].fillna(0)
    if exploreDF:
        print 'exploring feature df...'
        analysisDF=explore_feature_table(featureDF) #explore featureDF
        print featureDF.shape
    if deleteConsts:
        print 'deleting constant variables...'
        featureDFnameAddition=featureDFnameAddition+'_noConsts'
        constantdf=analysisDF[analysisDF['Type']=='constant']
        constantVars=constantdf['columns'].tolist()
        print 'constant variables are: %s' %constantVars
        print 'featureDF shape before constant features removal is %s_%s' %(featureDF.shape[0],featureDF.shape[1])
        featureDF=featureDF.drop(constantVars,axis=1)
        print 'Final featureDF shape  is %s_%s' %(featureDF.shape[0],featureDF.shape[1])
        print 'end of function!'
        
        return featureDF,analysisDF, featureDFnameAddition
    


#-----------------------------------------------------------------------------------
'''
the following function takes a dataset and sample List and generate a df that summarizes 'TCR annotation loads' for all samples in the 
sample list 
***** note that if the annotation file is updated/changed, it needs to be updated within the function!************

input: datasetfolder- string
datasetName-string
sampleList- list of sample names (e.g. ['BD100', 'BD200'])
sampleListName-string:

output:
annotationFeatureDF- all calculated annotation feature load without removal of fullu corrected and na filling
annotationFeatureDF_processed - after processing
'''

def gen_annotation_features_for_dataset(datasetFolder,datasetName,sampleList,sampleListName):
    #get unique aa list for dataset:
    f1='%s/sharingAnalysis/AllUniqueWithCounts' %datasetFolder
    AllUniqueWithCounts=pd.read_pickle(f1)
    print ('AllUniqueWithCounts shape is',AllUniqueWithCounts.shape)
    
    # get annoation file:
    f2='%s/TCR CDR3 sequence databases/combined annotation_list_clean_popped.xlsx' %MyPath
    annotations=pd.read_excel(f2)
    print ('annotations shape is',annotations.shape)
    
    #merge unique aa list with annotation list and save
    AllUniqueWithCountsAndAnnotations=pd.merge(AllUniqueWithCounts,annotations,how='left',left_index=True,right_index=True)
    print ('AllUniqueWithCountsAndAnnotations:',AllUniqueWithCountsAndAnnotations.shape)
    
    print 'saving AllUniqueWithCountsAndAnnotations'
    f3='%s/sharingAnalysis/AllUniqueWithCountsAndAnnotations' %datasetFolder
    AllUniqueWithCountsAndAnnotations.to_pickle(f3)
    
    annotationFolder='%s/annotationFeatures' %datasetFolder
    if not isdir(annotationFolder):
        makedirs(annotationFolder)

    annotationFeatureDF=pd.DataFrame(index=sampleList)
    AllUniqueWithCountsAndAnnotations['Sample_clean']=AllUniqueWithCountsAndAnnotations['Sample'].str.replace('.tsv','')
    
    print 'calculating annotation load for each sample:'
    for sample in sampleList:
        sampleDF=AllUniqueWithCountsAndAnnotations[AllUniqueWithCountsAndAnnotations['Sample_clean']==sample]
        print sample, len(sampleDF)
        sampleAnnotationDF=sampleDF.groupby('combined annotation_list_clean').agg({'Sample_clean': 'count',
                                                        'frequencyCount (%)': 'sum'})
        sampleAnnotationDF=sampleAnnotationDF.rename(columns={'Sample_clean':'seq_count',
                                                                                 sampleAnnotationDF.columns[1]:'cum_freq(%)'})
        sampleAnnotationDF.index=sampleAnnotationDF.index.rename('Annotation')
        
        #calculate total annotation load and relative anotation load for each annotation:
        totalAnnotateSeqs=sampleAnnotationDF['seq_count'].sum()
        totalAnnotatefreqs=sampleAnnotationDF['cum_freq(%)'].sum()
        sampleAnnotationDF['rel_seq_count']=(sampleAnnotationDF['seq_count'].astype('float'))*100/totalAnnotateSeqs.round(2)
        sampleAnnotationDF['rel_cum_freq']=(sampleAnnotationDF['cum_freq(%)'].astype('float'))*100/totalAnnotatefreqs.round(2)
        
        #summarize annotation features for all samples:
        print sampleAnnotationDF.sort_values(by='cum_freq(%)',ascending=False)
        for annot in sampleAnnotationDF.index:
            annotationFeatureDF.loc[sample,'%s_seq_count' %annot]=sampleAnnotationDF.loc[annot,'seq_count']
            annotationFeatureDF.loc[sample,'%s_cum_freq(perc)' %annot]=sampleAnnotationDF.loc[annot,'cum_freq(%)']
            annotationFeatureDF.loc[sample,'%s_rel_seq_count' %annot]=sampleAnnotationDF.loc[annot,'rel_seq_count']
            annotationFeatureDF.loc[sample,'%s_rel_cum_freq(perc)' %annot]=sampleAnnotationDF.loc[annot,'rel_cum_freq']
        annotationFeatureDF.loc[sample,'totalAnnotateSeqs']=totalAnnotateSeqs
        annotationFeatureDF.loc[sample,'totalAnnotatefreqs']=totalAnnotatefreqs 
        
          
    f4='%s/%s_%s_annotationFeatures.xlsx' %(annotationFolder,datasetName,sampleListName)
    annotationFeatureDF.to_excel(f4)
    
    #explore df and consider cleaning it using the designated function:
    print 'exploring features, removing constants adter filling nas with 0s'
    annotationFeatureDF_processed,analysisDF,featureDFnameAddition=process_featureDF(annotationFeatureDF,corrThreshold=0.99,
                                        removeCorr=True,exploreDF=True,deleteConsts=True,fillnans=True)
    
    f5='%s/%s_%s_annotationFeatures_ANALYSIS.xlsx' %(annotationFolder,datasetName,sampleListName)
    analysisDF.to_excel(f5)
    
    f6='%s/%s_%s_annotationFeatures' %(annotationFolder,datasetName,sampleListName) + featureDFnameAddition +'.xlsx'
    annotationFeatureDF_processed.to_excel(f6)
    
    print 'end of annotation function!'
    
    return annotationFeatureDF,annotationFeatureDF_processed
    
    
      
#----------------------------------------------------------------------------------------------------------------
'''
the following function merges feature DFs, PCAdf and TCRdf into one X matrix
'''


def gen_TCRfeatureX(datasetFolder,datasetName,incSingleSeqs=False,percShared=10,withRels=True):
    
    #get classical TCR feature df:
    featureFilePickle2='%s/featureSummaryDFs/%s_filteredBy%s_allFeatures_noCorrelated_noConsts_filledna' %(datasetFolder,datasetName,datasetName)
    TCRfeatureDF2=pd.read_pickle(featureFilePickle2)
    print 'TCRfeatureDF shape is %s_%s' %(TCRfeatureDF2.shape[0],TCRfeatureDF2.shape[1])
    print TCRfeatureDF2.iloc[:4,:4]
    
    #get annot_features
    if withRels:
        annotFeatFile='%s/annotationFeatures/%s_%s_annotationFeatures_withRels_noCorr0-99_nanFilled_noConsts.xlsx' %(datasetFolder,datasetName,datasetName)
        nameAddition='_withRels'
    else:
        annotFeatFile='%s/annotationFeatures/%s_%s_annotationFeatures_noCorr0-99_nanFilled.xlsx' %(datasetFolder,datasetName,
                                                                                                   datasetName)
        nameAddition=''
        
    TCRannotDF=pd.read_excel(annotFeatFile)
    print 'TCRannotDF shape is %s_%s' %(TCRannotDF.shape[0],TCRannotDF.shape[1])
    print TCRannotDF.iloc[:4,:4]
    
    #get TCRdf:
    TCRdfFile='%s/sharingAnalysis/sharingMatrix_%s_minNshared2_RA_onlyProductiveTrue__percShared%s_OLtrimmed_binary' %(datasetFolder,datasetName,percShared)
    TCRdf=pd.read_pickle(TCRdfFile)
    print 'TCRdf shape is %s_%s' %(TCRdf.shape[0],TCRdf.shape[1])
    print TCRdf.iloc[:4,:4]
    
    #GET PCAdf:
    PCAdfFile='%s/sharingAnalysis/PCAdf_%spercShared' %(datasetFolder,percShared)
    PCAdf=pd.read_pickle(PCAdfFile)
    print 'PCAdf shape is %s_%s' %(PCAdf.shape[0],PCAdf.shape[1])
    print PCAdf.iloc[:4,:4]
    
    #merge all:
    TCRfeatureX=TCRfeatureDF2
    TCRfeatureX=pd.merge(TCRfeatureX,TCRannotDF,how='left',left_index=True,right_index=True)
    TCRfeatureX=pd.merge(TCRfeatureX,PCAdf,how='left',left_index=True,right_index=True)
    
    if incSingleSeqs:
        TCRfeatureX=pd.merge(TCRfeatureX,TCRdf,how='left',left_index=True,right_index=True)
        nameAddition=nameAddition+'withTCRdf'
        
    print 'TCRfeatureX shape is %s_%s' %(TCRfeatureX.shape[0],TCRfeatureX.shape[1])
    print TCRfeatureX.iloc[:4,:4]
    
    featureFolder='%s/TCR_real_data/Predictions/featureDFs' %MyPath
    if not isdir(featureFolder):
        makedirs(featureFolder)
        
        
    f1='%s/allTCRfeatures_percShared%s_%s.dat' %(featureFolder,percShared,nameAddition)
    TCRfeatureX.to_pickle(f1)
    
    print 'saved feature df'
    
    return TCRfeatureX
    

    
    
    
      
      

    

  
    
    
