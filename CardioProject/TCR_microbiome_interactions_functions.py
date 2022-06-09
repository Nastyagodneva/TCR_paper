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
from scipy.stats import pearsonr, fisher_exact,spearmanr
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
from ShaniBA.TCR_microbiome_interactions.TCR_microbiome_interactions_functions2 import *


import os
# from Utils import cacheOnDisk
from queue.qp import qp, fakeqp
from addloglevels import sethandlers

MyPath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF'
print 'done'

#---------------------------------------------------------------------------------------------------------------
'''
the following functions enable the calculations of Fisher tests for each combination of seqsXspecies using the designated seq and species
files

in order to run them, run:
find_associations_binaryTCR_binaryMBspecies(NsharedSamplesForSeqs,topNspecies)

where:
NsharedSamplesForSeqs= is the number of minimum samples shared for each sequence (minimal value is 50!)
topNspecies= the number of species to use (the species that appear in the highest number of samples will be chosen



'''

def fisher_TCRseq_mbSpecies(speciesList, sequenceList, species_sequences_df, resultFolder):
    from scipy.stats import fisher_exact

    countTests = 0
    for i, species in enumerate(speciesList):
#         if i<10:
            # check if at least 10 of the samples don't have the species, otherwise, skip species

            nSamplesSpecies = species_sequences_df[species].sum()
#           if (nSamplesSpecies < len(species_sequences_df)-10) and (nSamplesSpecies > 0):
            print i, species
            groups = species_sequences_df.groupby(species)
            for j, seq in enumerate(sequenceList):
                nSamplesSeq = species_sequences_df[seq].sum()
#                 if j % 200 == 0:
#                     print 'seq number %s' % j
                df4 = pd.DataFrame(index=['seq_abs', 'seq_pres'], columns=['spec_abs', 'spec_pres'])
                for name, group in groups:
                    if name == 0:
                        seq_pres = group[seq].sum()
                        df4.loc['seq_pres', 'spec_abs'] = seq_pres
                        df4.loc['seq_abs', 'spec_abs'] = len(group) - seq_pres
                    if name == 1:
                        seq_pres = group[seq].sum()
                        df4.loc['seq_pres', 'spec_pres'] = seq_pres
                        df4.loc['seq_abs', 'spec_pres'] = len(group) - seq_pres
                df4 = df4.fillna(0)
                if not df4.isnull().values.any():
                    countTests = countTests + 1
                    OR, p = fisher_exact(df4, alternative='two-sided')
#                     print 'p=%s' %p
                    if p < 0.05:
                        nSamplesSeq = species_sequences_df[seq].sum()
                        speciesNegSeqNeg = df4.loc['seq_abs', 'spec_abs']
                        speciesNegSeqPos = df4.loc['seq_pres', 'spec_abs']
                        speciesPosSeqNeg = df4.loc['seq_abs', 'spec_pres']
                        speciesPosSeqPos = df4.loc['seq_pres', 'spec_pres']
                        result = pd.DataFrame()
                        result.loc[0, 'species'] = species
                        result.loc[0, 'seq'] = seq
                        result.loc[0, 'seq #'] = j
                        result.loc[0, 'speciesNegSeqNeg'] = speciesNegSeqNeg
                        result.loc[0, 'speciesNegSeqPos'] = speciesNegSeqPos
                        result.loc[0, 'speciesPosSeqNeg'] = speciesPosSeqNeg
                        result.loc[0, 'speciesPosSeqPos'] = speciesPosSeqPos
                        result.loc[0, 'OR'] = OR
                        result.loc[0, 'p'] = p
#                         print result
                        file1 = '%s/%s_seq%s' % (resultFolder, species, j)
                        result.to_pickle(file1)



                else:
                    print i, j
                    print 'contingent table contains null values:'
                    print df4
#         else:
#             print 'species %s in absent in less then 10 samples' % species

#     print countTests
    return countTests, result

def find_associations_binaryTCR_binaryMBspecies(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, datasetFolder, NsharedSamplesForSpecies, NsharedSamplesForSeqs,
                                          topNspecies, topNseqs):
    
    print 'preparing combined seq and species file...'
    species_sequences_df, speciesList, sequenceList = prepare_MbDF_and_TCRdf_for_analysis(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, NsharedSamplesForSpecies, NsharedSamplesForSeqs,
                                          topNspecies, topNseqs)
    # convert species_sequences_df to binary:
    print 'converting species_sequences_df to binary...'
    for column in species_sequences_df.columns.values:
        species_sequences_df[column] = np.where(species_sequences_df[column] > 0, 1, 0)
    
    # calculating fisher test:
    print 'executing fisher tests:...'
    resultFolder = '%s/TCR_mb_results/FisherResults_%s_%s_NsharedSpecies%s_Nsharedseqs%s_topSpecies%s_topSeqs%s.xlsx' % (datasetFolder, MbDFName,
                                    TCRdfName, NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies, topNseqs)
    if not isdir(resultFolder):
        makedirs(resultFolder)
    resultDFfolder = '%s/singleDFs' % resultFolder
    if not isdir(resultDFfolder):
        makedirs(resultDFfolder)
    
    countTests, result = fisher_TCRseq_mbSpecies(speciesList, sequenceList, species_sequences_df, resultDFfolder)
    with open('%s/countTests_fisher_%s_%s_%s_%s' % (datasetFolder, MbDFName, TCRdfName, NsharedSamplesForSeqs, topNspecies), 'wb') as fp:
        pickle.dump(countTests, fp)
    
    # generate resultSummary file:
    print 'generating a result summary file:...'
    FisherResults = concat_summarizing_dfs(resultDFfolder)
    FisherResults_full = add_corrected_pValues(FisherResults, pValueColumn='p', nTests=countTests, FDR=None)
    FisherResults_full = add_corrected_pValues(FisherResults, pValueColumn='p', nTests=countTests, FDR=0.01)
    FisherResults_full = FisherResults_full.sort_values(by='p')
    
    file1 = '%s/FisherResults.xlsx' % resultFolder
    FisherResults_full.to_excel(file1)
    print 'done'
    

    return FisherResults_full
    
#--------------------------------------------------------

'''
the following functions enable the calculations of correlation correlation tests for each combination of seqsXspecies using 
the designated seq and species RA files

in order to run them, run:
find_associations_TCR_RA_MBspecies_RA(MbDF_RA, MbDFName,TCRdf_RA, TCRdfName, NsharedSamplesForSpecies,NsharedSamplesForSeqs,
                                          topNspecies, topNseqs,outlierSTD,resultFolder)

where:
MbDF_RA=species DF with log relative abundances (df)
MbDFName - species DF name (string)
TCRdf_RA=seqs DF with log relative abundances (df)
TCRdfName - seq DF name (string)
NsharedSamplesForSpecies= is the number of minimum samples shared for each species. USE NONE if filtering by topNspecies
NsharedSamplesForSeqs= is the number of minimum samples shared for each sequence (minimal value is 50!) USE NONE if filtering by topNseqs
topNspecies= the number of species to use (the species that appear in the highest number of samples will be chosen. USE NONE if filtering by
 NsharedSamplesForSpecies
topNseqs= the number of seqs to use (the seqs that appear in the highest number of samples will be chosen. USE NONE if filtering by
 NsharedSamplesForSeqs
outlierSTD = number of STDs below/above the mean to filter outliers (chosse 'None' for not filtering)
resultFolder=folder to save results, an example (merge the two lines!):
'%s/TCR_real_data/PhenotypicData/Microbiome-TCR interactions/CorrelResults_%s_%s_NsharedSpecies%s_Nsharedseqs%s_topSpecies%s
_topSeqs%s_outlierSTD%s' %(MyPath,MbDFName,TCRdfName,NsharedSamplesForSpecies,NsharedSamplesForSeqs,
                                          topNspecies, topNseqs,outlierSTD)
completecorrelation=True/False (execute correlation, save df only if p<0.05, save many parameters)
correlationForHeatMap=True/False (execute correlation, save all results df, save only sample names and r)



'''

def correlation_TCRseq_mbSpecies(corrTest, speciesList, sequenceList, species_sequences_df, resultFolder, outlierSTD):
    

    countTests = 0
    for i, species in enumerate(speciesList):
        print i, species
        if outlierSTD is not None:  # remove outliers in species abundances
            species_sequences_df = filter_outliers(df=species_sequences_df, outlierSTD=outlierSTD, columnList=[species])
            print 'after outlier removal by species abundances, number of samples is %s' % len(species_sequences_df)
        for j, seq in enumerate(sequenceList):
            if outlierSTD is not None:  # remove outliers in seq abundances:
                species_sequences_df = filter_outliers(df=species_sequences_df, outlierSTD=outlierSTD, columnList=[seq])
#                 print 'after outlier removal by seq abundances, number of samples is %s' %len(species_sequences_df)
            nSamplesSeq = species_sequences_df[seq].sum()
#             if j % 200 == 0:
#                 print 'seq number %s' % j
            species_sequences_df4 = species_sequences_df[(species_sequences_df[species] != 0) | (species_sequences_df[seq] != 0)]  
            NsamplesWithSpeciesOrSeq = len(species_sequences_df4)
            if NsamplesWithSpeciesOrSeq > len(species_sequences_df) * 0.1:  # take only pairs for which at least 10% of samples are not double negative
                countTests = countTests + 1
                x = species_sequences_df[seq]
                y = species_sequences_df[species]
                if corrTest == 'spearman':
                    r, p = MySpearmanr(x, y)
                else:
                    r, p = MyPearsonr(x, y)
                
                if p < 0.05:
                    nSamplesSeq = species_sequences_df[seq].sum()
                    result = pd.DataFrame()
                    result.loc[0, 'species'] = species
                    result.loc[0, 'seq'] = seq
                    result.loc[0, 'seq #'] = j
                    result.loc[0, 'NsamplesWithSpeciesOrSeq'] = NsamplesWithSpeciesOrSeq
                    result.loc[0, 'r'] = r
                    result.loc[0, 'p'] = p
#                         print result
                    file1 = '%s/%s_seq%s' % (resultFolder, species, j)
                    result.to_pickle(file1)
            else:
                print i, j
                print 'not enough samples containing the species or the seqs'

#     print countTests
    return countTests
'''
below is another version of this function, used for calculating a df of r values which will then be used for generating a pivot table,
a heat map and a cluster map. this function is called by 'find_associations_TCR_RA_MBspecies_RA_forHeatMap' function.
see comment for former function for more information.
this function was copied from 'species X TCR correlation heatmap.ipynb'

'''

def correlation_TCRseq_mbSpecies_forHeatMap(corrTest, speciesList, sequenceList, species_sequences_df, resultFolder, outlierSTD):
    

    countTests = 0
    for i, species in enumerate(speciesList):
        print i, species
        if outlierSTD is not None:  # remove outliers in species abundances
            species_sequences_df = filter_outliers(df=species_sequences_df, outlierSTD=outlierSTD, columnList=[species])
            print 'after outlier removal by species abundances, number of samples=%s' % len(species_sequences_df)
        for j, seq in enumerate(sequenceList):
            if outlierSTD is not None:  # remove outliers in species abundances
                species_sequences_df = filter_outliers(df=species_sequences_df, outlierSTD=outlierSTD, columnList=[seq])
              #  print 'after outlier removal by seq abundances, number of samples is %s' %len(species_sequences_df)
            
            nSamplesSeq = species_sequences_df[seq].sum()
#             if j % 200 == 0:
#                 print 'seq number %s' % j
            countTests = countTests + 1               
            x = species_sequences_df[seq]
            y = species_sequences_df[species]
            if corrTest == 'spearman':
                r, p = MySpearmanr(x, y)
            else:
                r, p = MyPearsonr(x, y)
                
            result = pd.DataFrame()
            result.loc[0, 'species'] = species
            result.loc[0, 'seq'] = seq
            result.loc[0, 'r'] = r
#                         print result
            file1 = '%s/%s_seq%s' % (resultFolder, species, j)
            result.to_pickle(file1)
            
#     print countTests
    return countTests

'''
The following function calles 'Correlation_TCRseq_mbSpecies' and together they calculate correlation correlations for all possible
seq X species pairs. see 'Correlation_TCRseq_mbSpecies' function documentation for requested input explanations.
'''

def find_associations_TCR_RA_MBspecies_RA(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, datasetFolder, NsharedSamplesForSpecies, NsharedSamplesForSeqs,
                                          topNspecies, topNseqs, outlierSTD, corrTest, completecorrelation, correlationForHeatMap):
    
    # get combined df of seqs and species according to the requested thresholds:
    print 'preparing combined seq and species file...'
    species_sequences_df, speciesList, sequenceList = prepare_MbDF_and_TCRdf_for_analysis(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, NsharedSamplesForSpecies, NsharedSamplesForSeqs,
                                          topNspecies, topNseqs)
    
    # calculating correlation test:
    resultFolder = '%s/TCR_mb_results/CorrResults_%s_%s_%s_Thresholds(%s%s%s%s)_olSTD%s'\
    % (datasetFolder, corrTest, MbDFName, TCRdfName, NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies,
       topNseqs, outlierSTD)
    if not isdir(resultFolder):
        makedirs(resultFolder)

    corrResults_full = pd.DataFrame()
    corrResultsHeatmap = pd.DataFrame()
    
    if completecorrelation:
        allDFsFolder = '%s/singleDFs' % resultFolder
        if not isdir(allDFsFolder):
            makedirs(allDFsFolder)
        print 'executing complete correlation correlation tests:...'
        countTests = correlation_TCRseq_mbSpecies(corrTest, speciesList, sequenceList, species_sequences_df, allDFsFolder, outlierSTD)  
        print 'generating a result summary file:...'
        corrResults = concat_summarizing_dfs(allDFsFolder)
        corrResults_full = add_corrected_pValues(corrResults, pValueColumn='p', nTests=countTests, FDR=None)
        corrResults_full = add_corrected_pValues(corrResults_full, pValueColumn='p', nTests=countTests, FDR=0.01)
        corrResults_full = corrResults_full.sort_values(by='p')
        
        print 'saving files:'
        with open('%s/countTests_completecorrelation' % (resultFolder), 'wb') as fp:
            pickle.dump(countTests, fp) 
        file1 = '%s/corrResults_%s.xlsx' % (resultFolder, corrTest)
        corrResults_full.to_excel(file1)
        print 'done'
    
    if correlationForHeatMap:
        heatMapDFFolder = '%s/singleHeatMapResultDF' % resultFolder
        if not isdir(heatMapDFFolder):
            makedirs(heatMapDFFolder)
        print 'executing correlation test for heatmap...'
        countTests = correlation_TCRseq_mbSpecies_forHeatMap(corrTest, speciesList, sequenceList, species_sequences_df, heatMapDFFolder, outlierSTD) 
        print 'generating a result summary file:...'
        corrResultsHeatmap = concat_summarizing_dfs(heatMapDFFolder)
        corrResultsHeatmap = corrResultsHeatmap.sort_values(by='r')
    
        file2 = '%s/corrResultsHeatmap_%s.xlsx' % (resultFolder, corrTest)
        corrResultsHeatmap.to_excel(file2)
        print 'done'
    return corrResults_full

#------------------------------------------------------
'''
the following functions is called by the binary or RA TCR X MB analysis function, and prepare the df for analysis:
'''



def prepare_MbDF_and_TCRdf_for_analysis(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, NsharedSamplesForSpecies, NsharedSamplesForSeqs,
                                          topNspecies, topNseqs):
    
    # prepare seq file:
    print 'preparing seq file:...'   

    # #convert SEQ file to binary:
    TCRdf_binary = pd.DataFrame()
    for column in TCRdf_RA.columns.values:
        TCRdf_binary[column] = np.where(TCRdf_RA[column] > 0, 1, 0)
    TCRdf_binary = TCRdf_binary.set_index(TCRdf_RA.index)
    
    
    # ## filter seq binary file according to minimal number of samples shared OR number of top seqs to use:
    if NsharedSamplesForSeqs is not None:
        # #filter seq binary file to include only species shared by more than NsharedSamplesForSeqs:
        print 'filtering sequence df according to minimal number of shared samples by seqs=%s' % NsharedSamplesForSeqs
        seqCount = TCRdf_binary.sum()
        if NsharedSamplesForSeqs < seqCount.min():
            raise NameError('NsharedSamplesForSeqs must be %s or higher!' % seqCount.min())
        elif (NsharedSamplesForSeqs > seqCount.min()) and (NsharedSamplesForSeqs < len(TCRdf_binary)):    
            selectedSeqs = list(seqCount[seqCount > NsharedSamplesForSeqs].index)  
            TCRdf_RA = TCRdf_RA[selectedSeqs]
        else:
            raise NameError('NsharedSamplesForSeqs must be lower than %s!' % len(TCRdf_binary))
    elif topNseqs is not None:
        # #filter SEQ binary file to include only sequences shared by more than NsharedSamplesForSeqs:
        print 'filtering sequence df to include only top %s sequences ' % topNseqs
        seqCount = TCRdf_binary.sum().sort_values(ascending=False)
        topSeqs = list(seqCount.index[:topNseqs])
        TCRdf_RA = TCRdf_RA.loc[:, topSeqs]
    else:  
        print 'filtering to leave only sequences that appear in more than 10% percent of the samples:' 
        NsharedSamplesForSeqs=round(len(TCRdf_binary)*0.1,0)
        print '10 percent of samples is %s' %NsharedSamplesForSeqs
        seqCount = TCRdf_binary.sum()
        selectedSeqs = list(seqCount[seqCount > NsharedSamplesForSeqs].index)  
        TCRdf_binary = TCRdf_binary[selectedSeqs]
        TCRdf_RA = TCRdf_RA[selectedSeqs]
        
    TCRdf_RA = editSampleNames(TCRdf_RA)
        
    # prepare mb species files:
    print 'preparing mb file:...'   

    # #convert mb file to binary:
    MbDF_binary = pd.DataFrame()
    for column in MbDF_RA.columns.values:
        MbDF_binary[column] = np.where(MbDF_RA[column] > 0, 1, 0)
    MbDF_binary = MbDF_binary.set_index(MbDF_RA.index)
    
    
    # ## filter species binary file according to minimal number of samples shared OR number of top species to use:
    if NsharedSamplesForSpecies is not None:
        # #filter species binary file to include only species shared by more than NsharedSamplesForSpecies:
        print 'filtering species df according to minimal number of shared samples by species=%s' % NsharedSamplesForSpecies
        speciesCount = MbDF_binary.sum()
        if NsharedSamplesForSpecies < speciesCount.min():
            raise NameError('NsharedSamplesForSpecies must be %s or higher!' % speciesCount.min())
        elif (NsharedSamplesForSpecies > speciesCount.min()) and (NsharedSamplesForSpecies < len(MbDF_binary)):    
            selectedSpecies = list(speciesCount[speciesCount > NsharedSamplesForSpecies].index)  
            MbDF_RA = MbDF_RA[selectedSpecies]
        else:
            raise NameError('NsharedSamplesForSpecies must be lower than %s!' % len(MbDF_binary))
    elif topNspecies is not None:
        # #filter species binary file to include only topNspecies:
        print 'filtering species df to include only top %s species ' % topNspecies
        speciesCount = MbDF_binary.sum().sort_values(ascending=False)
        topSpecies = list(speciesCount.index[:topNspecies])
        MbDF_RA = MbDF_RA.loc[:, topSpecies]
    else: 
        print 'filtering according to minimal relative abundance=0.0001'
        print 'converting MbDF_RA to binary according to minimal relative abundace of 0.0001...'
        MbDF_binary = pd.DataFrame()
        print 'number of columns to convert is %s' %len(MbDF_RA.columns.values)
        for n,column in enumerate(MbDF_RA.columns.values):
            if n%200==0:
                print n
            MbDF_binary[column] = np.where(MbDF_RA[column] > 0.0001, 1, 0)
            if MbDF_binary[column].sum()==0:
                print 'species %s is being dropped...' %column
                MbDF_binary=MbDF_binary.drop(column,axis=1)
        MbDF_binary = MbDF_binary.set_index(MbDF_RA.index)
        MbDF_RA = MbDF_RA[MbDF_binary.columns.values]  

 
    # generate a merged file containing all sequences and species:
    print 'generating merged file:...'
    sequenceList = TCRdf_RA.columns.values
    speciesList = MbDF_RA.columns.values

    species_sequences_df = pd.merge(TCRdf_RA, MbDF_RA, how='inner',
                             left_index=True, right_index=True)
    print 'validating columns in merged file:'
    print 'sequenceList length is %s' % len(sequenceList)
    print 'SpeciesList length is %s' % len(speciesList)  
    print 'total length of merged file columns =%s' % len(species_sequences_df.columns.values)
    print 'validating rows in merged file:'
    print 'sequence matrix row number is %s' % len(TCRdf_RA)
    print 'species matrix row number is %s' % len(MbDF_RA)
    print 'merged matrix row number is %s' % len(species_sequences_df)
    print ' potential number of tests is %s' % str(len(sequenceList) * len(speciesList))
    
    # #remove same person samples:
    species_sequences_df = filterSamePerson(species_sequences_df, [0])
    print 'same person samples were removed, now matrix row number is %s' % len(species_sequences_df)
    
    return species_sequences_df, speciesList, sequenceList



'''
this function was copied from 'species X TCR correlation heatmap.ipynb'
This function generates an heatmap and a cluster map for Ratble -  a df with r values for each seqXspecies pairs:

input:
*Rtable= a df with the following columns: 'species','seq','r' (generated by the function find_associations_TCR_RA_MBspecies_RA using
*RtableName=string, used to generate the figure names
correlationForHeatMap=True) 
*figFolder=folder to save the heatmap and clustermap, recommended to use the Rtable folder
*figsizeX= fig width
*figsizeY = fig height
*fontsizeX = xtickslabels font size 
*fontsizeY= ytickslabels font size 
*fontSizeLabel = x and ylabel fontsize
*sortMethod = 'max'/'mean' - how to order the rows and columns of the heatmap. 
'''
def gen_heatmap_and_clustermap(Rtable, RtableName, figFolder,indexCol,columnsCol,valueCol,metric='euclidean',
                               method='average',sortMethod='max',cmapName='hsv'):
                               
#                                figsizeX, figsizeY, fontsizeX, fontsizeY, fontSizeLabel, sortMethod):
    
           
                  
    figsizeX=float(Rtable[indexCol].nunique())/8
    figsizeY=float(Rtable[columnsCol].nunique())/8
    
    print 'Rtable shape is %s_%s' % (Rtable.shape[0], Rtable.shape[1])
    print 'number of unique %s=%s' % (Rtable[indexCol].nunique(),indexCol)
    print 'number of unique %s=%s' % (Rtable[columnsCol].nunique(),columnsCol)  
        
    # generate pivot table:
    print 'generating pivot table...'
    pivot = Rtable.pivot(index=indexCol, columns=columnsCol, values=valueCol)
    print pivot.shape
    if sortMethod == 'max':
        pivot[sortMethod] = pivot.max(axis=1)
        pivot.loc[sortMethod, :] = pivot.max(axis=0)
    elif sortMethod == 'mean':
        pivot[sortMethod] = pivot.mean(axis=1)
        pivot.loc[sortMethod, :] = pivot.mean(axis=0)
    pivot = pivot.sort_values(by=sortMethod, axis=1)
    pivot = pivot.sort_values(by=sortMethod, axis=0)
    pivot = pivot.drop(sortMethod, axis=1)
    pivot = pivot.drop(sortMethod, axis=0)
    
    
    if indexCol=='species':
        # generate genus color list:
        uniqueSpecies = pivot.index.unique()
        genus = [x.split('_')[0] for x in uniqueSpecies ]
        genus2 = list(set(genus))
        genus2n = len(genus2)
        if genus2n < 10:
            fact = 0.1
        elif genus2n < 20:
            fact = 0.05
        elif genus2n < 50:
            fact = 0.02
        else:
            fact = 0.01 
        species_genus_df = pd.DataFrame({'species':uniqueSpecies, 'genus':genus})
        genus2color = []
        for n, g in enumerate(genus2):
            color = plt.cm.get_cmap(cmapName)(fact * n)
            genus2color.append(color)
        genus2colorDF = pd.DataFrame({'genus':genus2, 'color':genus2color})
        speciesGenusColorDF = pd.merge(species_genus_df, genus2colorDF, how='left', left_on='genus', right_on='genus')
        cluster_colors = list(speciesGenusColorDF['color'])
    
    # generate heatmap:
    print 'generate heatmap...'
    heatMapFig, ax = plt.subplots(figsize=(figsizeX, figsizeY))  # Sample figsize in inches
    sns.heatmap(pivot, ax=ax)
#     ax.yaxis.label.set_size(fontSizeLabel)
#     ax.tick_params(axis='y', labelsize=fontsizeY)
#     ax.xaxis.label.set_size(fontSizeLabel)
#     ax.tick_params(axis='x', labelsize=fontsizeX)
    heatMapName = '%s_heatMap_%sBy%s_sortBy%s' % (RtableName, indexCol,columnsCol, sortMethod)
    heatMapFile='%s/%s' %(figFolder,heatMapName)
    heatMapFig.savefig(heatMapFile, dpi=200)
    
    # generate clustermap
    print 'generate clustermap...'
    if indexCol=='species': 
        cg = sns.clustermap(pivot, figsize=(figsizeX, figsizeY), col_cluster=True, row_colors=cluster_colors,metric=metric,
                            method=method)
    else:
        cg = sns.clustermap(pivot, figsize=(figsizeX, figsizeY), col_cluster=True,metric=metric,method=method)
    
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels())
    clusterMapName = '%s_clusterMap_%sBy%s_%s_%s' % (RtableName, indexCol,columnsCol,metric,method)
    clusterMapFile='%s/%s' %(figFolder,clusterMapName)
    cg.savefig(clusterMapFile, dpi=200)
    
    print 'done'
    
    return heatMapFig, cg
    
    #-------------------------------------------------------------------------------------------------------

'''
the following function plots the relationship between a specific sequence and a specific species.
make sure to download the sequences and the species files!

an example for usage:

file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PublicSeqAnalysis/TCRdf'
TCRdf=pd.read_pickle(file1)
file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/SampleSpeciesDFgroupedByBD_only434'
SampleSpeciesDFgroupedByBD_only434_RA=pd.read_pickle(file1)
seq='CSARQGDTEAFF'
species='Ruminococcus_albus'
check_seq_species_relations_sns(seq,species,color_binary=True)

'''

def check_seq_species_relations_sns(datasetFolder, datasetName, seq, species, color_binary,outlierSTD=None): 
    
    print 'loading most updated Mb file...'
    f1 = '%s/MicrobiomeDataTables/FilteredAndMergedOnBD/MPA_s_AllSeqProjects_SampleListPNP515_filterGenotekTrue_filterMinimalRead9000000_meannSTDNonenMinSamplesNone_filterOutlierSampleFalse_filterSamePersonFalse' % MyPath
    groupedFilteredMB = pd.read_pickle(f1)
       
        
    # (2) prepare MB and TCR files for analysis:
    print 'preparing mb and TCR files for analysis:...'
    
    # prepare MB file:
    # load BD_FD converter, to make sure we remove BD_FD columns from the MB file:
    f1 = '%s/Sample files/updatedBDandFDlists/final_BD_FD_converter_mergedOnDate_RCfolderAllSeqProjects_31052018' % MyPath
    BD_FD = pd.read_pickle(f1)
    
    try:
        MbDF_RA = groupedFilteredMB.set_index('BD')  #### very important!!
    except KeyError:
        print 'couldnt reset index to BD,probably BD is already the index'
        MbDF_RA = groupedFilteredMB
        

    # remove unnecessary columns:
    for column in BD_FD.columns.values:
        if column in MbDF_RA.columns.values:
            print 'dropping %s column' % column
            MbDF_RA = MbDF_RA.drop(column, axis=1)
            
    # prepare TCR file:
    minNshared = 50  # add to function variables if necessary
    file1 = '%s/sharingAnalysis/sharingMatrix_%s_minNshared%s_RA' % (datasetFolder, datasetName, minNshared)
    TCRdf_RA = pd.read_pickle(file1)
 
    seqData = pd.DataFrame(TCRdf_RA[seq])
    seqData = seqData.rename(columns={0:seq})
    seqData = editSampleNames(seqData)
    speciesData = pd.DataFrame(MbDF_RA[species])
    merged = pd.merge(seqData, speciesData, how='inner', left_index=True, right_index=True)
    print 'n rows after merging=%s' % len(merged)
    merged = filterSamePerson(merged,[0])
    print 'n rows after same person removal=%s' % len(merged)
        
    x = merged[seq]
    y = merged[species]
    r, p = MyPearsonr(x, y)
    print 'correlation:'
    print r, p
    if outlierSTD is not None:
        merged = filter_outliers(df=merged, outlierSTD=outlierSTD, columnList=[seq])
        merged = filter_outliers(df=merged, outlierSTD=outlierSTD, columnList=[species])
    print 'n rows after outlier removal=%s' % len(merged)
    x = merged[seq]
    y = merged[species]
    r, p = MyPearsonr(x, y)
    print 'correlation after outlier removal'
    print r, p
    
    xmax = round(np.max(x) * 1.05, 4)
    ymax = round(np.max(y) * 1.05, 4)
    from scipy import stats
    g = sns.JointGrid(x, y)
    g = g.plot_joint(plt.scatter, color="0.5", edgecolor='white', s=100, alpha=0.3)
    g = g.plot_joint(sns.regplot)
    for ind in range(0, len(merged)):
        xpoint = merged[seq][ind]
        ypoint = merged[species][ind]
        if (xpoint > np.percentile(merged[seq], 99)) | (ypoint > np.percentile(merged[species], 99)):
            plt.annotate(merged.index[ind], xy=(xpoint, ypoint),
                     xytext=(xpoint * 1.05, ypoint * 1.05), fontsize=8)
        if color_binary:
            if xpoint > 0 and ypoint == 0:
                plt.scatter(xpoint, ypoint, color='red', edgecolor='red', s=100, alpha=0.3)
            if ypoint > 0 and xpoint == 0:
                plt.scatter(xpoint, ypoint, color='blue', edgecolor='blue', s=100, alpha=0.3)
    plt.annotate('n=%s' % len(x), xy=(0.16, 0.83), xycoords='axes fraction', fontsize=16)
    plt.xlim(0, xmax)
    plt.ylim(0, ymax)
    plt.xticks(rotation=90)
    g = g.plot_marginals(sns.distplot, kde=False, color="0.5")
    g = g.annotate(stats.spearmanr, fontsize=16)
    
    
    plt.show()
    
    return g

#----------------------------------------------------------------
'''
this function takes specific sequence and species from the sharingMatrix and the SampleSpeciesDF
and plot the correlation between them

not that in the script calling the function, there should be additional information such as:
fig, axes=plt.subplots(nrows=2, ncols=3,figsize=(20,10))
...
plt.subplots_adjust(left=0.02, bottom=None, right=0.99, top=0.84,
                wspace=0.2, hspace=0.25)

'''

def plot_seq_species_scatter(ax,corrTest,datasetFolder, datasetName, seq, species, color_binary,onlyNextera,outlierSTD=None):

    if onlyNextera:
        print 'loading most updated Mb file...filtered by libPrepMethod'
        f1 = '%s/MicrobiomeDataTables/FilteredAndMergedOnBD/\
MPA_s_AllSeqProjects_SampleListPNP515_filterGenotekTrue_filterMinimalRead9000000_libMethNextera_meannSTDNonenMinSamplesNone_fOutlierSampleFalse_fSamePersonFalse' % MyPath
    else:
        print 'loading most updated Mb file...NOT filtered by libPrepMethod'
        f1 = '%s/MicrobiomeDataTables/FilteredAndMergedOnBD/\
MPA_s_AllSeqProjects_SampleListPNP515_filterGenotekTrue_filterMinimalRead9000000_libMethNone_meannSTDNonenMinSamplesNone_fOutlierSampleFalse_fSamePersonFalse_RA' % MyPath

    
    
    groupedFilteredMB = pd.read_pickle(f1)
       
        
    # (2) prepare MB and TCR files for analysis:
#     print 'preparing mb and TCR files for analysis:...'
    
    # prepare MB file:
    # load BD_FD converter, to make sure we remove BD_FD columns from the MB file:
    f1 = '%s/Sample files/updatedBDandFDlists/final_BD_FD_converter_mergedOnYear_RCfolderAllSeqProjects_30042018' % MyPath
    BD_FD = pd.read_pickle(f1)
    
    try:
        MbDF_RA = groupedFilteredMB.set_index('BD')  #### very important!!
    except KeyError:
#         print 'couldnt reset index to BD,probably BD is already the index'
        MbDF_RA = groupedFilteredMB
        

    # remove unnecessary columns:
    for column in BD_FD.columns.values:
        if column in MbDF_RA.columns.values:
#             print 'dropping %s column' % column
            MbDF_RA = MbDF_RA.drop(column, axis=1)
            
    # prepare TCR file:
    minNshared = 50  # add to function variables if necessary
    file1 = '%s/sharingAnalysis/sharingMatrix_%s_minNshared%s_RA' % (datasetFolder, datasetName, minNshared)
    TCRdf_RA = pd.read_pickle(file1)
    
    
    
    seqData = pd.DataFrame(TCRdf_RA[seq])
    seqData = seqData.rename(columns={0:seq})
    seqData = editSampleNames(seqData)
    speciesData = pd.DataFrame(MbDF_RA[species])
    merged = pd.merge(seqData, speciesData, how='inner', left_index=True, right_index=True)
#     print 'n rows after merging=%s' % len(merged)
    merged = filterSamePerson(merged,[0])
#     print 'n rows after same person removal=%s' % len(merged)
    if outlierSTD is not None:
        merged = filter_outliers(df=merged, outlierSTD=outlierSTD, columnList=[seq])
        merged = filter_outliers(df=merged, outlierSTD=outlierSTD, columnList=[species])
#     merged2 = merged[(merged[species] != 0) | (merged[seq] != 0)]
    x = merged[seq]
    y = merged[species]
    xy=merged[[seq,species]]
    usedSamples=len(xy.dropna())
    if corrTest=='pearson':
        r, p = MyPearsonr(x, y)
    else:
        r, p = MySpearmanr(x, y)
    ax.scatter(x, y, color='black', alpha=0.3)
    xmax = round(np.max(x) * 1.05, 4)
    ymax = round(np.max(y) * 1.05, 4)

    for ind in range(0, len(merged)):
            xpoint = merged[seq][ind]
            ypoint = merged[species][ind]
#             if (xpoint>np.percentile(merged2[seq],99))|(ypoint>np.percentile(merged2[species],99)):
#                 ax.annotate(merged2.index[ind],xy=(xpoint,ypoint),
#                          xytext=(xpoint*1.05,ypoint*1.05),fontsize=8)
            if color_binary:
                if xpoint > 0 and ypoint == 0:
                    ax.scatter(xpoint, ypoint, color='red', edgecolor='red', s=100, alpha=0.3)
                if ypoint > 0 and xpoint == 0:
                    ax.scatter(xpoint, ypoint, color='blue', edgecolor='blue', s=100, alpha=0.3)
    ax.annotate('%sr= %s,p= %s\nn=%s' % (corrTest,round(r, 2), round(p, 6), usedSamples), xy=(0.16, 0.8), xycoords='axes fraction', fontsize=10)
    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    ax.set_xlabel(seq,fontsize=12)
    ax.set_ylabel(species,fontsize=12)
    plot_bestFitLine(x, y, ax, color='black')
    
    return ax

#----------------------------------------------------------------------
'''
ths following function can be used to generate a dataframe containing all sample information
binary counts and RA) for input sequence and species
'''

### this function was copied to 'TCR_microbiome_interactions_functions.py' - use and update there!###


def get_seq_species_data(seq, species):
    # load binary seq file:
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PublicSeqAnalysis/sharingMatrix_moreThan50_434Samples'
    sharingMatrix_moreThan50_434Samples = pd.read_pickle(file1)
    # get the sequence binary counts:
    seqBinaryData = pd.DataFrame(sharingMatrix_moreThan50_434Samples[seq])
    # get the sequence RA data:
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PublicSeqAnalysis/sharingMatrix_moreThan50_434Samples_RA'
    sharingMatrix_moreThan50_434Samples_RA = pd.read_pickle(file1)
    # get the sequence RA data:
    seqRAData = pd.DataFrame(sharingMatrix_moreThan50_434Samples_RA[seq])
    # merge and edit names:
    seqData = pd.merge(seqBinaryData, seqRAData, how='inner', left_index=True, right_index=True)
    seqData = editSampleNames(seqData)
    # get the mb RA data:
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/SampleSpeciesDFgroupedByBD_only434'
    SampleSpeciesDFgroupedByBD_only434 = pd.read_pickle(file1)
    # get the species RA data:
    speciesRAData = pd.DataFrame(SampleSpeciesDFgroupedByBD_only434[species])
    # transform to binary:
    SpeciesBinaryData = pd.DataFrame(np.where(speciesRAData > 0, 1, 0)).set_index(speciesRAData.index)
    # combine species Data:
    SpeciesData = pd.merge(SpeciesBinaryData, speciesRAData, how='inner', left_index=True, right_index=True)
    # combineAll:
    seqSpeciesData = pd.merge(seqData, SpeciesData, how='inner', left_index=True, right_index=True)
    for column in seqSpeciesData.columns.values:
        newColumn = str(column)
        newColumn = newColumn.replace('_x', '_binary')
        newColumn = newColumn.replace('_y', '_RA')
        newColumn = newColumn.replace(species, species + '_RA')
        newColumn = newColumn.replace('0', species + '_binary')
        seqSpeciesData = seqSpeciesData.rename(columns={column:newColumn})
    print len(seqSpeciesData)
    print seqSpeciesData.head()
    
    return seqSpeciesData

#--------------------------------------------------------

'''
THIS FUNCTION WAS DEVELOPED IN NOTEBOOK 'microbiomeDataGeneration'
the following function extracts MPA data for the required taxonomial level, organized it, and save RA and binary versions of the data
*level input - 'd'/'k'/'p'/'c'/'o'/'f'/'g'/'s'
(domain/kingdom/phylum/class/order/family/genus/species)
*dataFolder='Metabolon2'/'AllSeqProjects'

usage example:
MPAlevel,MPAlevelBinary=gen_mb_df('s','Metabolon2')

'''

def gen_mb_df(level, dataFolder):
    
    # load species info:
    print 'loading MPASpid file from DFOut folder in %s' % dataFolder
    if dataFolder == 'Metabolon2':
        f1 = '/net/mraid08/export/jafar/Microbiome/Analyses/Metabolon2/DFOut/MPA.dat'
    else:
        f1 = '/net/mraid08/export/jafar/Microbiome/Analyses/AllSeqProjects/DFOut/MPASpid.dat'
    
    
    
    MPASpid = pd.read_pickle(f1)
    print 'MPASpid table length (number of all taxa) is %s' % len(MPASpid)
    
    # take only data for the desired tax level ()
    MPASpid = MPASpid.reset_index()
    MPAlevel = MPASpid[MPASpid['TaxLevel'] == level]
    MPAlevel = MPAlevel.set_index('Tax')
    MPAlevel = MPAlevel.drop('TaxLevel', axis=1)
    print 'MPA data table length for level %s is %s' % (level, len(MPAlevel))
    
    # correct column names:
    print 'correcting sample names (columns)...'
    for column in MPAlevel.columns.values:
        newColumn = column.split('_')[0]
        MPAlevel = MPAlevel.rename(columns={column:newColumn})
        
    # transpose table:
    print 'transposing table, now rows are samples and columns are taxa...'
    MPAlevel = MPAlevel.T
    MPAlevel.index.rename('FD', inplace=True)
    #
    # generate binary version:
    print 'generating binary version of the table...'
    MPAlevelBinary = pd.DataFrame(index=MPAlevel.index)
    # MPAlevelBinary.index=MPAlevel.index
    for column in MPAlevel.columns.values:
        MPAlevelBinary[column] = np.where(MPAlevel[column] > 0, 1, 0)
        
    # fillna in RA table:
    MPAlevel = MPAlevel.fillna(0)
    
    print 'RA table shape is %s_%s' % (MPAlevel.shape[0], MPAlevel.shape[1])
    print 'binary table shape is %s_%s' % (MPAlevelBinary.shape[0], MPAlevel.shape[1])
    print 'should be identical'
    
    # save tables:

    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/MicrobiomeDataTables/MPA_%s_%s_RA' % (level, dataFolder)
    MPAlevel.to_pickle(f1)
    f3 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/MicrobiomeDataTables/MPA_%s_%s_binary' % (level, dataFolder)
    MPAlevelBinary.to_pickle(f3)
    
    return MPAlevel, MPAlevelBinary


'''

THIS function was copied from 'Filtering sample lists.ipynb'

this function takes an FD-based data table (can use also other identifiers such as SPID, connection)
merge it with the BD-FD conversion table, and then filter the data based on specific sample list, only
swab samples and/or minimal reads per sample.

the required input:
dataDF: FD-based table, with samples as rows and phenotypes as columns (for example - microbiome species 
relative abundances) ####make sure the column to merge on is not the index!! ###
dataDFname - a string, such as MPA_s
dataMergeOn- the name of the column in the dataDF table to merge on ###don't forget to reset_index if necessary
BDFDMergeOn- the name of the column in the BD_FD table to merge on
Sample list - None, or python list of BD numbers in the format: ['BD1','BD2'...]
SampleListName-None or the list name
filterGenotek - True/False. True will lead to filtering out of all samples with isGenotek==1. samples with
no indication at all won't be filtered out.
filterMinimalReads - None, or a number. recommended value is 9000000
filterlibPrepMethod - None, 'TrueSeq' or 'Nextera'
folderToSave=None or full path to the required folder

USAGE EXAMPLE:
# SPECIES LEVEL, PNP434samples, no additional filtering:

f1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/MicrobiomeDataTables/notFiltered/MPA_s_Metabolon2_RA'
MPAlevel=pd.read_pickle(f1)

#LOAD pnp434 sample list:
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP434samples','rb') as fp:
    PNP434samples=pickle.load(fp)
    

readCountFolder='Metabolon2'
dataDF=MPAlevel.reset_index()
dataDFname='MPA_s_Metabolon2'
dataMergeOn='FD'
BDFDMergeOn='FD'
SampleList=PNP434samples
SampleListName='PNP434'
filterGenotek=False
filterMinimalReads=None
filterlibPrepMethod='TrueSeq'
folderToSave='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/MicrobiomeDataTables/Filtered'



dataWithMeta=filter_data_on_FDproperties(dataDF,dataDFname,dataMergeOn,BDFDMergeOn,SampleList=SampleList,
                                SampleListName=SampleListName,filterGenotek=filterGenotek, filterMinimalReads=filterMinimalReads,
                               filterlibPrepMethod=filterlibPrepMethod,folderToSave=folderToSave)

'''
def filter_data_on_FDproperties(dataDF, dataDFname, dataMergeOn, BDFDMergeOn, SampleList=None,
                                SampleListName=None, filterGenotek=False, filterMinimalReads=9000000,
                               filterlibPrepMethod=None, folderToSave=None, BD_FD=None):
    
    # ## merge with bd-fd table:
    # load bd-fd table (can change the BD-FD converter version by changing the date in the file name)
    
    if BD_FD is None:
        print 'using BD_FD converter based on AllSeqProjects folder, date 31052018, merged on date'
        f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/final_BD_FD_converter_mergedOnDate_RCfolderAllSeqProjects_31052018'          
        BD_FD = pd.read_pickle(f1)
    print 'row number in BD_FD table is %s' % len(BD_FD)
    
    # merging:
    print 'merging dataDF with FD-BD table:'

    dataWithMeta = pd.merge(dataDF, BD_FD, how='left', left_on=dataMergeOn, right_on=BDFDMergeOn)
    print 'row number in dataDF before merging was %s' % len(dataDF)
    print 'row number after merging is %s' % len(dataWithMeta)
    
    # ## if required - filter on specific BDs, if not-filter on samples that have any BD:
    if SampleList is None:
        print 'filtering out all FD samples without BDs'
        dataWithMeta = dataWithMeta[dataWithMeta['BD'].notnull()]
        print 'row number after filtering is %s' % len(dataWithMeta)
    else:
        print 'filtering out all FD that are not matched to any BD in the sample list' 
        dataWithMeta = dataWithMeta[dataWithMeta['BD'].isin(SampleList)]
        print 'row number after filtering is %s' % len(dataWithMeta)
        
    # ## if required filter on isGenotek
    if filterGenotek:
        print 'filtering out genotek samples'
        dataWithMeta = dataWithMeta[dataWithMeta['isGenotek'] != 1]
        print 'row number after filtering is %s' % len(dataWithMeta)
        
    # ## if required, filter on minimal number of reads
    if filterMinimalReads is not None:
        print 'filtering out samples with read number after HGF=%s' % filterMinimalReads
        dataWithMeta = dataWithMeta[dataWithMeta['PostHGF'] > filterMinimalReads]
        print 'row number after filtering is %s' % len(dataWithMeta)
        
    # ## if required, filter on libPrepMethod
    if filterlibPrepMethod is not None:
        print 'using only samples prepared with %s' % filterlibPrepMethod
        dataWithMeta = dataWithMeta[dataWithMeta['libPrepMethod'] == filterlibPrepMethod]
        print 'row number after filtering is %s' % len(dataWithMeta)
    
        
    # ## indicating final number of BDs and FDs in the table:
    nBDs = dataWithMeta['BD'].nunique()
    nFDs = dataWithMeta['FD'].nunique()

    print 'final number of unique BDs is %s and unique FDs is %s' % (nBDs, nFDs)
    
    # ## save to file:
    fileName = '%s_%s_fGenotek%s_fMinimalRead%s_libMeth%s' % (dataDFname, SampleListName, filterGenotek,
                                                                filterMinimalReads, filterlibPrepMethod)
    if folderToSave is None:
        folderToSave = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/FilteredDataTables'

    print 'saving filtered data table to folder %s' % folderToSave   
    if not isdir(folderToSave):
            makedirs(folderToSave)
    f1 = '%s/%s' % (folderToSave, fileName)    
    dataWithMeta.to_pickle(f1)
    
    return dataWithMeta, fileName

#--------------------------------------------------------------------------------------
    
'''
the following functions takes DF tables that include 'BD' column and groups the data based on this column


input:
*DFtoGroup - sampleXphenotype dataframe, 'BD' column should exist
*DFtoGroupName - a string
*groupFunction='noOutlierMean'/'mean'/'mostFrequent'/'binarySelection'
*filterOutlierSamples=True/False
*filterSamePerson=True/False
*folderToSave=None/full path to folder

*nSTD - must be None if mergeFunction is not 'noOutlierMean' this defines the number of STD from the mean
to define low/high threshold for non-outlier samples. this apply to the BD merging procedure and not the 
outlier removal which is currently constant at 3. 
*nMinSamples - must be None if mergeFunction is not 'noOutlierMean'.this defines the minimal number of 
samples to apply samples filtering. this apply to the BD merging procedure and not the 
outlier removal which is currently constant at 3. 

USAGE EXAMPLE:
f1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/MicrobiomeDataTables/Filtered/MPA_s_SampleListNone_filterGenotekFalse_filterMinimalReadNone'
MPAnotfiltered=pd.read_pickle(f1)


DFtoGroup=MPAnotfiltered
DFtoGroupName='MPAnotfiltered'
groupFunction='noOutlierMean'
filterOutlierSamples=True
filterSamePerson=False
folderToSave='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/MicrobiomeDataTables/FilteredAndMergedOnBD'
nSTD=3,
nMinSamples=3

groupedDF=mergeOnBD(DFtoGroup,DFtoGroupName,groupFunction,filterOutlierSamples,filterSamePerson,folderToSave,nSTD=nSTD,
              nMinSamples=nMinSamples)
'''

def mergeOnBD(DFtoGroup, DFtoGroupName, groupFunction, filterOutlierSamples, filterSamePerson, nSTD=None,
              nMinSamples=None,ignoreNotSameYear=True):
    print ' n unique BD is %s' %DFtoGroup['BD'].nunique()
    if ignoreNotSameYear:
        DFtoGroup=DFtoGroup[DFtoGroup['SameYear']==1]
        print ' n unique BD after removing samples not matched at least by year is %s' %DFtoGroup['BD'].nunique()
    DFtoGroup=DFtoGroup.drop('SameYear',axis=1)
    
    ## remove unecessary columns:
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/final_BD_FD_converter_mergedOnDate_RCfolderAllSeqProjects_31052018'          
    BD_FD = pd.read_pickle(f1)
    for column in DFtoGroup.columns.values:
        if column in BD_FD.columns.values:
            if (column !='BD') and (column !='FD'):
                print 'dropping column %s from matrix' %column
                DFtoGroup=DFtoGroup.drop(column,axis=1)
    
    # ## group on BD, using the selected mergeFunction:
    print 'df length before grouping is %s' % len(DFtoGroup)
    print 'grouping on BD using %s' % groupFunction
    if groupFunction == 'noOutlierMean':
        groupedDF = DFtoGroup.groupby('BD').agg(lambda x: noOutlierMean(x, nSTD, nMinSamples))
    elif groupFunction == 'mean':
        groupedDF = DFtoGroup.groupby('BD').mean()
    elif groupFunction == 'mostFrequent':   
        import collections
        groupedDF = DFtoGroup.groupby('BD').agg(lambda x: collections.Counter(x).most_common()[0][0])
    elif groupFunction == 'binarySelection':
        groupedDF = DFtoGroup.groupby('BD').agg(lambda x: 1 if 1 in list(x) else 0)
    print 'grouped table length after grouping by BD and outlier removal=%s' % len(groupedDF)
    
    # ## filter same person:
    if filterSamePerson:     
        # get again UserIDs for BD samples:
        print 'loading again BD_FD table'
        f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/final_BD_FD_converter_mergedOnDate_31052018'
        BD_FD = pd.read_pickle(f1)
        BD_FD_toUse = BD_FD[['BD', 'UserID']]
        BD_FD_toUse = BD_FD_toUse.drop_duplicates()
        print 'number of BD-UserID pairs is %s' % len(BD_FD_toUse)

        print 'merging UserID info to BDgrouped table:'
        groupedDF = pd.merge(groupedDF, BD_FD_toUse, how='left', left_index=True, right_on='BD')

        # counting how many User IDs have more than one BD:
        BDrepeats = pd.DataFrame(groupedDF['UserID'].value_counts())
        nRepeatingBDs = len(BDrepeats[BDrepeats['UserID'] > 1])
        print '%s UserIDs have more than one BD' % nRepeatingBDs
        
        # dropping same person samples (leaving the last one)
        print 'groupedDF length before dropping same person samples is %s' % len(groupedDF)
        groupedDF = groupedDF.drop_duplicates(subset='UserID', keep='last')
        print 'groupedDF length after dropping same person samples is %s' % len(groupedDF)
        groupedDF = groupedDF.set_index('BD')
        
    # ## if required, filter outlier samples:
    if filterOutlierSamples:
        print 'filtering outlier samples...'
        groupedDF = filter_phenotypiclly_outlier_samples(groupedDF, nSTD=3, nMinSamples=3)
        print 'groupedDF length after dropping outlier samples is %s' % len(groupedDF)
        
    # ## save groupedDF:
    # ## save to file:
    fileName3 = '%s_%snSTD%snMinSamples%s_fSamePerson%s' % (DFtoGroupName, groupFunction,
                                                                nSTD, nMinSamples, filterSamePerson)   
    
    return groupedDF, fileName3    
 
 
'''
 the following function generates MbDF matrix, which is actually MPA matrix agter filtration and merging on BD number
 
 USAGE EXAMPLE:
 
 #load PNP515 sample list:
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP515','rb') as fp:
    PNP515=pickle.load(fp)
    
#DEFINE PARAMETERS:
mbLevel='s'
mbDataFolder='AllSeqProjects'
SampleList=PNP515
SampleListName='PNP515'
filterGenotek=True
filterMinimalReads=9000000
filterlibPrepMethod=None
groupFunction='noOutlierMean'
nSTD=5
nMinSamples=3
ignoreNotSameYear=True (don't use samples that are matched only based on UserID but their storage year is different

groupedFilteredMB=genFilterMergeMbData(mbLevel,mbDataFolder,SampleList,SampleListName,filterGenotek,filterMinimalReads,filterlibPrepMethod,
                        groupFunction,nSTD,nMinSamples,ignoreNotSameYear)
 
 '''
    
    
def genFilterMergeMbData(mbLevel,mbDataFolder,SampleList,SampleListName,filterGenotek,filterMinimalReads,filterlibPrepMethod,
                        groupFunction,nSTD,nMinSamples,ignoreNotSameYear,useShortName):
    
    print 'generating mb file for level %s based on folder %s' % (mbLevel, mbDataFolder)
    MPAlevel, MPAlevelBinary = gen_mb_df(mbLevel, mbDataFolder)
    dataDF = MPAlevel.reset_index()
    
    
    print 'filtering data...'
    dataDFname = 'MPA_%s_%s' % (mbLevel, mbDataFolder)
    folderToSave = '%s/MicrobiomeDataTables/Filtered' % MyPath
    filteredMb, fileName = filter_data_on_FDproperties(dataDF, dataDFname, dataMergeOn='FD', BDFDMergeOn='FD',
                                SampleList=SampleList, SampleListName=SampleListName,
                                filterGenotek=filterGenotek, filterMinimalReads=filterMinimalReads,
                                filterlibPrepMethod=filterlibPrepMethod, folderToSave=folderToSave, BD_FD=None)
    
    print 'grouping on BD...'
    folderToSave = '%s/MicrobiomeDataTables/FilteredAndMergedOnBD' % MyPath
    fileName2=fileName+'_igAnotherYear%s' %ignoreNotSameYear
    groupedFilteredMB,filename3 = mergeOnBD(filteredMb, fileName2, groupFunction=groupFunction, filterOutlierSamples=False,
                        filterSamePerson=False, nSTD=nSTD, nMinSamples=nMinSamples,
                        ignoreNotSameYear=ignoreNotSameYear)
    
       
    print 'saving mergedOnBDs data table to folder %s' % folderToSave   
    if not isdir(folderToSave):
            makedirs(folderToSave)
    if useShortName:
        f1 = '%s/MPA_%s_standardParams' % (folderToSave, mbLevel)  
    else:
        f1=  '%s/%s' % (folderToSave, filename3)
    groupedFilteredMB.to_pickle(f1)
    
    print 'Done'
        
    return groupedFilteredMB,f1




#-------------------------------------------------------------------------------------------------------------------------------------
    
def filterSamePerson(df, axisList):
    # get again UserIDs for BD samples, it is not important to load new file as the userID information is full in older files as well
    print 'loading BD_FD table'
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/BDfile_31052018.xlsx'
    BD_FD = pd.read_excel(f1)
    BD_FD = BD_FD.rename(columns={'DnaID':'BD'})
    BD_FD_toUse = BD_FD[['BD', 'UserID']]
    BD_FD_toUse = BD_FD_toUse.drop_duplicates()
    print 'number of BD-UserID pairs is %s' % len(BD_FD_toUse)
    
    print 'merging UserID info to BDgrouped table:'
    df = pd.merge(df, BD_FD_toUse, how='left', left_index=True, right_on='BD')
    
    # counting how many User IDs have more than one BD:
    BDrepeats = pd.DataFrame(df['UserID'].value_counts())
    nRepeatingBDs = len(BDrepeats[BDrepeats['UserID'] > 1])
    print '%s UserIDs have more than one BD' % nRepeatingBDs
    
    # dropping same person samples (leaving the last one) NOTE- the way the function is written now, it drops the first occurence of a userID.
    # therefore if there are more than 2 repeats of the same user, more than one will be left in df
    print 'groupedDF length before dropping same person samples is %s' % len(df)
    duplicated = df[df.duplicated(subset='UserID', keep='first')]
#     print duplicated
    UserIDtoRemove = list(duplicated['UserID'])
    BDsToRemove = list(duplicated['BD'])
    print 'removing UserIDs: %s' % UserIDtoRemove
    print 'removing BD: %s' % BDsToRemove
    df = df.set_index('BD')
    for droppedBD in BDsToRemove:
        for axis in axisList:
            df = df.drop(droppedBD, axis=axis)
    df = df.drop('UserID', axis=1)
    # df=df.drop_duplicates(subset='UserID',keep='last',axis=axis)
    print 'groupedDF shape after dropping same person samples is %s_%s' % (df.shape[0], df.shape[1])
    
    
    return df

#----------------------------------------------------------------------------------------------------

'''
the following function takes a dataset and check TCR-microbiome interactions within it. 
the main outputs are:
1. correlation correlation results for relative abundance data(all results with p<0.05)
2. fisher test results for binary data (p<0.05)
3. merged results with TCR identity. 

input:
datasetFolder=string
datasetName=string
gen_mb=True/False. if needed to generate new mb RA dataframe. if not, the data used is all FD 
samples corresponding to PNP515 BD samples, filtered for swabs only and for minimal number of 
reads (post HGF) of 9k.
mblevel - string, 's'=species, 'g'=genius etc. (None of gen_mb=False)
mbDataFolder=string (None of gen_mb=False)
sampleList=list of BD numbers to filter by (None of gen_mb=False)
SampleListName=string (None of gen_mb=False)
filterGenotek=True/False
filterMinimalReads=True/False
thresholdList=list of quadriple tuples, representing NsharedSamplesForSpecies,NsharedSamplesForSeqs,
topNspecies,topNseqs (the first couple or the second couple should be None but not both!)
outlierSTD=STD to use for outlier sample removal, normally=3
corrTest - 'spearman'  or 'pearson'
completecorrelation-True/False - whether to generate full result df for each seq-species pair with p<0.05
correlationForHeatMap=True/False whether to generate a short df for ALL seq-species pairs (used for heatmap/
clustermap generation)
usePermCorr - whether to use permutations or not in correlation test
nPermCorr- number of permutation to execute - note! this is an ordinal number. therefore if desire to generate permutations number 100-200, define 200 here and 
99 in firstPermCorr. if desired to skip permutations calculation, define nPermCorr lower then firstPermCorr
firstPermCorr - the number of first permutation to execute, minus 1.
runRealCorr- whether to skip the part of generating real r for each seq-species pair
runFisher - True/False
runCorr - True/False
mergeResults - True/False

NOTE!
*in order to use permutation analysis but skip the permutation process and go directly to permutation df
merging and p-value calculation - use nPermCorr=None, irstPermCorr=None (or the same for fisher).
*in order to skip the correlation and/or fisher test completely and go directly to result merging use runCorr=False 
and/or runFisher=False

USAGE EXAMPLE:
#load PNP515 sample list:
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP515','rb') as fp:
    PNP515=pickle.load(fp)
print len(PNP515)
print PNP515[:5]


datasetFolder='%s/TCR_real_data' %MyPath
datasetName='PNP515'

gen_mb=False
mbLevel='s'
mbDataFolder='AllSeqProjects'
SampleList=PNP515
SampleListName='PNP515'
filterGenotek=True
filterMinimalReads=9000000

# (NsharedSamplesForSpecies,NsharedSamplesForSeqs,topNspecies,topNseqs):
thresholdList=[(None,None,50,50)]
outlierSTD=3
corrTest='spearman'
completecorrelation=True
correlationForHeatMap=False

usePermCorr=True
nPermCorr=None
firstPermCorr=None
runRealCorr=False
runFisher=False
mergeResults=True

TCR_Mb_interactions_for_dataset(datasetFolder, datasetName, gen_mb, mbLevel, mbDataFolder, SampleList,
                                    SampleListName, filterGenotek, filterMinimalReads,
                                    thresholdList, outlierSTD, corrTest, completecorrelation,
                                    correlationForHeatMap,usePermCorr,nPermCorr,firstPermCorr,runRealCorr)



'''

def TCR_Mb_interactions_for_dataset(datasetFolder, datasetName, gen_mb, mbLevel, mbDataFolder, SampleList,
                                    SampleListName, filterGenotek, filterMinimalReads, filterlibPrepMethod,
                                    thresholdList, outlierSTD, corrTest, completecorrelation,
                                    correlationForHeatMap, usePermCorr, runCorr, runFisher, usePermFisher, nPermCorr, nPermFisher,
                                    firstPermCorr, firstPermFisher, runRealCorr, runRealFisher,mergeResults,minNshared = 50):
    
    nPermCorrUsed=0.0
    nPermFisherUsed=0.0
    # (1) if required, generate mb file for a specific level and filter it:
    if gen_mb:
        print 'generating mb file for level %s based on folder %s' % (mbLevel, mbDataFolder)
        MPAlevel, MPAlevelBinary = gen_mb_df(mbLevel, mbDataFolder)
        dataDF = MPAlevel.reset_index()
        dataDFname = 'MPA_%s_%s' % (mbLevel, mbDataFolder)
        folderToSave = '%s/MicrobiomeDataTables/Filtered' % MyPath
        print 'filtering data...'
        filteredMb, fileName = filter_data_on_FDproperties(dataDF, dataDFname, dataMergeOn='FD', BDFDMergeOn='FD',
                                    SampleList=SampleList, SampleListName=SampleListName,
                                    filterGenotek=filterGenotek, filterMinimalReads=filterMinimalReads,
                                    filterlibPrepMethod=filterlibPrepMethod, folderToSave=folderToSave, BD_FD=None)

        print 'grouping on BD...'
        folderToSave = '%s/MicrobiomeDataTables/FilteredAndMergedOnBD' % MyPath

        groupedFilteredMB = mergeOnBD(filteredMb, fileName, groupFunction='mean', filterOutlierSamples=False,
                            filterSamePerson=False, folderToSave=folderToSave, nSTD=None, nMinSamples=None)
    else:
        print 'loading most updated Mb file...'
        f1 = '%s/MicrobiomeDataTables/FilteredAndMergedOnBD/MPA_s_AllSeqProjects_SampleListPNP515_filterGenotekFalse_filterMinimalReadFalse_libMethNone_meannSTDNonenMinSamplesNone_fOutlierSampleFalse_fSamePersonFalse' % MyPath
        groupedFilteredMB = pd.read_pickle(f1)
        filterGenotek = True
        filterMinimalReads = 9000000
        
        
    # (2) prepare MB and TCR files for analysis:
    print 'preparing mb and TCR files for analysis:...'
    
    # prepare MB file:
    # load BD_FD converter, to make sure we remove BD_FD columns from the MB file:
    f1 = '%s/Sample files/updatedBDandFDlists/final_BD_FD_converter_mergedOnDate_RCfolderAllSeqProjects_31052018' % MyPath
    BD_FD = pd.read_pickle(f1)
    
    try:
        MbDF_RA = groupedFilteredMB.set_index('BD')  #### very important!!
    except KeyError:
        print 'couldnt reset index to BD,probably BD is already the index'
        MbDF_RA = groupedFilteredMB
        

    # remove unnecessary columns:
    for column in BD_FD.columns.values:
        if column in MbDF_RA.columns.values:
            print 'dropping %s column' % column
            MbDF_RA = MbDF_RA.drop(column, axis=1)
            
    # prepare TCR file:
    minNshared = minNshared #default is 50
    file1 = '%s/sharingAnalysis/sharingMatrix_%s_minNshared%s_RA' % (datasetFolder, datasetName, minNshared)
    TCRdf_RA = pd.read_pickle(file1)
   
    
    # (3) run tests and merge results - it is possible to run on different thresholds in a loop
    print 'running interaction tests...'
    
    MbDFName = 'Mb_Swabs%sMinReads%s' % (filterGenotek, filterMinimalReads)
    TCRdfName = 'TCRminNshared%s' % minNshared
    
    # load TCR identity table:
    file2 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR CDR3 sequence databases/CDR3identityTable_23042018'
    CDR3identityTable = pd.read_pickle(file2)
    CDR3identityTable.head()
    

    for t in thresholdList:
        NsharedSamplesForSpecies = t[0]
        NsharedSamplesForSeqs = t[1]
        topNspecies = t[2]
        topNseqs = t[3]
        print 'threshold parameters are:'
        print NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies, topNseqs
        
        print 'correlation test:...'
        permCorrTestFolder = '%s/TCR_mb_results/permCorrTest_%s%s%s%s_std%s' % (datasetFolder, NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies, topNseqs, outlierSTD)
        if not isdir(permCorrTestFolder):
            makedirs(permCorrTestFolder)
        if runCorr:
            if usePermCorr: 
                Final_permCorrResults_df, nPermCorrUsed = calc_permCorrTest(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies,
                          topNseqs, corrTest, nPermCorr, firstPermCorr, outlierSTD, permCorrTestFolder, runRealCorr)
                corrResults_full = Final_permCorrResults_df[Final_permCorrResults_df['real_p_Corr%sPerms' % nPermCorrUsed] <= 0.05]
            else:
                corrResults_full = find_associations_TCR_RA_MBspecies_RA(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, datasetFolder, NsharedSamplesForSpecies, NsharedSamplesForSeqs,
                                                                          topNspecies, topNseqs, outlierSTD, corrTest, completecorrelation, correlationForHeatMap)
        else:
            try:
                print 'correlation test was not executed...most updated permCorrResults file is loaded...'
                nPermCorrUsed = 9999
                file2 = '%s/Final_permCorrResults_%sPerms_df.xlsx' % (permCorrTestFolder, nPermCorrUsed)
                Final_permCorrResults_df = pd.read_excel(file2)
                corrResults_full = Final_permCorrResults_df[Final_permCorrResults_df['real_p_Corr%sPerms' % nPermCorrUsed] <= 0.05]
            except:
                print 'couldnt load former correlation result file...'                                        
        if runFisher:
            if usePermFisher:
                print 'Fisher test with permutations...'
                permFisherTestFolder = '%s/TCR_mb_results/permFisherTest_%s%s%s%s' % (datasetFolder, NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies, topNseqs)
                if not isdir(permFisherTestFolder):
                    makedirs(permFisherTestFolder)
                Final_permFisherResults_df, nPermFisherUsed = calc_permFisherTest(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies,
                      topNseqs, nPermFisher, firstPermFisher, permFisherTestFolder, runRealFisher)
                FisherResults_full = Final_permFisherResults_df[Final_permFisherResults_df['real_p_Fisher%sPerms' % nPermFisherUsed] <= 0.05]
            else: 
                print 'fisher test without permutations:...'
                FisherResults_full = find_associations_binaryTCR_binaryMBspecies(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, datasetFolder, NsharedSamplesForSpecies, NsharedSamplesForSeqs,
                                                  topNspecies, topNseqs)
        else:
            nPermFisherUsed = 9999
            resultFolder = '%s/TCR_mb_results/permFisherTest_%s%s%s%s' % (datasetFolder, NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies, topNseqs)
            file3 = '%s/Final_permFisherResults_%sPerms_df.xlsx' % (resultFolder, nPermFisherUsed)
            Final_permFisherResults_df = pd.read_excel(file3)
            FisherResults_full = Final_permFisherResults_df[Final_permFisherResults_df['real_p_Fisher%sPerms' % nPermFisherUsed] <= 0.05]
            
        if mergeResults:
            # # merge correlation and fisher results:
            merged_fisher_correlation = pd.merge(corrResults_full, FisherResults_full, how='inner',
                                          left_on=('species', 'seq'), right_on=('species', 'seq'))
            print 'Correlation result length is %s' % len(corrResults_full)
            print 'Fisher result length is %s' % len(FisherResults_full)
            print 'merged df length is %s' % len(merged_fisher_correlation)
            for column in merged_fisher_correlation:
                newColumn = column.replace('_x', '_correlation')
                newColumn = newColumn.replace('_y', '_Fisher')
                merged_fisher_correlation = merged_fisher_correlation.rename(columns={column:newColumn})
                if 'p_' in newColumn:
                    p_value_col = newColumn
            
            merged_fisher_correlation = merged_fisher_correlation.sort_values(by=p_value_col)
            
            # # add identity information to the table:
            merged_fisher_correlation_withIdentity = pd.merge(merged_fisher_correlation, CDR3identityTable, how='left',
                                               left_on='seq', right_index=True)
            print'CDR3 identity table length is %s' % len(CDR3identityTable)
            print 'merged result table with identity length is %s' % len(merged_fisher_correlation_withIdentity)
            
            merged_fisher_correlation_withIdentity.head()
            
            print 'saving result files:...'
            if usePermCorr:
                resultName = '%sPerm%s_t(%s%s%s%s)_std%s' \
                % (corrTest, usePermCorr, NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies, topNseqs, outlierSTD)
            else:
                resultName = '%sPerm%s_%sCorrPerms%FisherPerms_t(%s%s%s%s)_std%s' \
                % (corrTest, usePermCorr, nPermCorrUsed, nPermFisherUsed, NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies, topNseqs, outlierSTD)
    
            try:
                file4 = '%s/TCR_mb_results/mergedFisherCorr_%s.xlsx' % (datasetFolder, resultName)
                merged_fisher_correlation_withIdentity.to_excel(file4)
            except ValueError:
                print 'problames with saving file, therefore saved it in MyPath'
                file5 = '%s/mergedFisherCorr_%s.xlsx' % (MyPath, resultName)
                merged_fisher_correlation_withIdentity.to_excel(file5)
    else:
        print 'result merging was not executed'
        merged_fisher_correlation_withIdentity=pd.DataFrame()
            
        
        
    print 'done'
    return merged_fisher_correlation_withIdentity


#-------------------------------------------

'''
the following 3 functions are subfunctions of the next function: 'calc_permCorrTest'
'''


def corrTest_SeqSpecies_forPermutationTest(corrTest, speciesList, sequenceList, species_sequences_df, outlierSTD, isPerm):       
    corrResultsList = []
    countTests = 0
    for i, species in enumerate(speciesList):
        print i, species
        if outlierSTD is not None:  # remove outliers in species abundances
            species_sequences_df = filter_outliers(df=species_sequences_df, outlierSTD=outlierSTD, columnList=[species])
#             print 'after outlier removal by species abundances, number of samples=%s' % len(species_sequences_df)
        nSamplesSpecies = len(species_sequences_df[species_sequences_df[species] > 0])
        if nSamplesSpecies > 0:
            for j, seq in enumerate(sequenceList):
                if outlierSTD is not None:  # remove outliers in species abundances
                    species_sequences_df = filter_outliers(df=species_sequences_df, outlierSTD=outlierSTD, columnList=[seq])
                nSamplesSeq = len(species_sequences_df[species_sequences_df[seq] > 0])
                if nSamplesSeq > 0:
                    countTests = countTests + 1               
                    x = species_sequences_df[seq]
                    y = species_sequences_df[species]
                    if corrTest == 'spearman':
                        r, p = MySpearmanr(x, y)
                    else:
                        r, p = MyPearsonr(x, y)
        
                    result = pd.DataFrame()
                    result.loc[0, 'species'] = species
                    result.loc[0, 'seq'] = seq
                    result.loc[0, 'r'] = r
                    
                    if not isPerm:
                        nComSeqSpecies = len(species_sequences_df[(species_sequences_df[species] > 0) & (species_sequences_df[seq] > 0)])
                        result.loc[0, 'nSamplesSeq'] = nSamplesSeq
                        result.loc[0, 'nSamplesSpecies'] = nSamplesSpecies
                        result.loc[0, 'nComSeqSpecies'] = nComSeqSpecies
                    
                    corrResultsList.append(result)
                else:
                    print 'no samples with sequence %s' % seq
        else:
            print 'no samples with this species'

    corrResultsdf = pd.concat(corrResultsList)
                    
    return corrResultsdf
    
def permutateDFandCalcCorrelation(nPermCorr, species_sequences_df, corrTest, speciesList, sequenceList, outlierSTD, permDFsFolder):
    print nPermCorr
    print 'permutation number %s for seqSpecies def...' % nPermCorr
    perm_species_sequences_df = species_sequences_df.apply(lambda x: x.sample(frac=1).values)
#     print perm_species_sequences_df.shape
    print 'correlation test for permutated df...'
    perm_result_df = corrTest_SeqSpecies_forPermutationTest(corrTest, speciesList, sequenceList, perm_species_sequences_df,
                                                          outlierSTD, isPerm=True)
    print perm_result_df.shape
    permDFfile = '%s/perm%sDF' % (permDFsFolder, nPermCorr)
    perm_result_df.to_pickle(permDFfile)
    
    return perm_result_df
    
# the following function calculates how many of the values in a list are larger than the first object in the list
# good for calculating ampirical p-value from a list of ampirical statistics.
def perc_above(my_vals):
    t = abs(my_vals[0])
    data = [abs(x) for x in my_vals]
    count_vals = sum(i >= t for i in data)
    percentile_val = float(count_vals) / len(data)
    return percentile_val




'''
the following function execute correlation test (spearman/pearson) with permutations. see above its 3 sub-functions

input:
corrTest - 'spearman'/'pearson'
nPermCorr- number of permutation to execute - note! this is an ordinal number. therefore if desire to generate permutations number 100-200, define 200 here and 
99 in firstPermCorr. if desired to skip permutations calculation, define nPermCorr lower then firstPermCorr
firstPermCorr - the number of first permutation to execute, minus 1.
speciesList, sequenceList, species_sequences_df - outputs of the function ********
outlierSTD - to reject outliers, usually=3. write None if no outlier removal is desired
 permCorrTestFolder - the folder to hold results. define in the calling function as: '%s/TCR_mb_results/permCorrTest_%s%s%s%s' %(datasetFolder,
 NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies, topNseqs)
runRealCorr- whether to skip the part of generating real r for each seq-species pair

**** possible improvements: *** 
consider adding variables for trds_def and mem_def
define whether to use the cluster or not

Usage example:
nPermCorr=120
firstPermCorr=99
permCorrTestFolder='%s/TCR_mb_results/permCorrTest_%s%s%s%s' %(datasetFolder,NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies, topNseqs)

Final_permCorrResults_df=calc_permCorrTest(corrTest,nPermCorr,firstPermCorr,speciesList, sequenceList, species_sequences_df, 
                                           outlierSTD,permCorrTestFolder,runRealCorr=False)
'''

def calc_permCorrTest(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies,
                      topNseqs, corrTest, nPermCorr, firstPermCorr, outlierSTD, permCorrTestFolder, runRealCorr):
    
    # ## (1) get combined df of seqs and species according to the requested thresholds:
    print 'preparing combined seq and species file...'
    species_sequences_df, speciesList, sequenceList = prepare_MbDF_and_TCRdf_for_analysis(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, NsharedSamplesForSpecies, NsharedSamplesForSeqs,
                                          topNspecies, topNseqs)
    
        
    # ## (2) run correlation test for all seq-species pairs in the real data.
    # ##get all results in one df
    
    if runRealCorr:
        print 'running %s correlation test for all seq-species pairs in the real data' % corrTest
        real_result_df = corrTest_SeqSpecies_forPermutationTest(corrTest, speciesList, sequenceList, species_sequences_df,
                                                              outlierSTD, isPerm=False)
        real_result_df = real_result_df.rename(columns={'r':'r_real'})
        file1 = '%s/real_result_df' % (permCorrTestFolder)
        real_result_df.to_pickle(file1)
    else:
        print 'loading real R table...'
        file1 = '%s/real_result_df' % (permCorrTestFolder)
        real_result_df = pd.read_pickle(file1)
        
    
    # ## (3)permutate seq-species df nPermCorr times, and each time calculate full correlation table and save it - USE CLUSTER FOR CALCULATIONS!
    permCorrelationDFsFolder = '%s/permDFs' % permCorrTestFolder
    if nPermCorr > firstPermCorr:
        print 'permutations and correlation testing... running on the cluster!!! might be very long!!!'
        print 'total number of permutationsL %s' % nPermCorr
        print 'start after permutation number %s' % firstPermCorr
               
        
        if not isdir(permCorrelationDFsFolder):
            makedirs(permCorrelationDFsFolder)
    
        sethandlers()
        os.chdir(permCorrTestFolder)
    
        # send jobs to the cluster:
        with qp(jobname='permCorr', q=['himemint.q', 'himem7.q'], trds_def=1, mem_def='2G') as q:  # consider changing queue list!!
            q.startpermanentrun()
            tk = {}  # a dictionary storing arguments and results
            for i in range(nPermCorr + 1):
                if i > firstPermCorr:
                    print i
                    tk[i] = q.method(permutateDFandCalcCorrelation, (i, species_sequences_df, corrTest, speciesList, sequenceList, outlierSTD, permCorrelationDFsFolder))
            tk = {k:q.waitforresult(v) for k, v in tk.iteritems()}
        print 'done'
    else:
        print 'no new permutation tables are generated...'
    
    # ## (4) summarize all permutation results:
    print 'summarizing all permutation results...'
    
    Full_permCorrResults_df = real_result_df
    # merge all permDFs to the real r table:
    filenames = [f for f in listdir(permCorrelationDFsFolder) if isfile(join(permCorrelationDFsFolder, f))]
    nPermCorrUsed = len(filenames)
    print 'number of permutations used is %s' % nPermCorrUsed
    for n, f in enumerate(filenames):
    #     if n<5:
            print n, f
            file1 = '%s/%s' % (permCorrelationDFsFolder, f)
            df = pd.read_pickle(file1)
            Full_permCorrResults_df = pd.merge(Full_permCorrResults_df, df, how='inner',
                                            left_on=('species', 'seq'), right_on=('species', 'seq'))
            nPer = f.split('DF')[0]
            Full_permCorrResults_df = Full_permCorrResults_df.rename(columns={'r':'r_%s' % nPer})
            
    # calculate how many permutations have more extreme r value (this is the p-value, two-sided)
    print 'calculating real p values'
    columnsToApply = [column for column in Full_permCorrResults_df.columns.values if 'r_' in column]
    Full_permCorrResults_df['real_p_Corr%sPerms' % nPermCorrUsed] = Full_permCorrResults_df[columnsToApply].apply(perc_above, axis=1)
    
    nTests = len(speciesList) * len(sequenceList)

    Full_permCorrResults_df = add_corrected_pValues(Full_permCorrResults_df, pValueColumn='real_p_Corr%sPerms' % nPermCorrUsed, nTests=nTests, FDR=0.01)
    Full_permCorrResults_df = add_corrected_pValues(Full_permCorrResults_df, pValueColumn='real_p_Corr%sPerms' % nPermCorrUsed, nTests=nTests, FDR=0.1)
    Full_permCorrResults_df = add_corrected_pValues(Full_permCorrResults_df, pValueColumn='real_p_Corr%sPerms' % nPermCorrUsed, nTests=nTests, FDR=0.25)
    Full_permCorrResults_df = Full_permCorrResults_df.sort_values(by='r_real')
    
    columns_to_use = ['species', 'seq', 'r_real', 'nSamplesSeq', 'nSamplesSpecies',
       'nComSeqSpecies', 'real_p_Corr%sPerms' % nPermCorrUsed, 'Sig by bonferroni corrected pVal', 'sig. by FDR=0.01',
       'sig. by FDR=0.1', 'sig. by FDR=0.25']
    Final_permCorrResults_df = Full_permCorrResults_df[columns_to_use]
    file2 = '%s/Final_permCorrResults_%sPerms_df.xlsx' % (permCorrTestFolder, nPermCorrUsed)
    Final_permCorrResults_df.to_excel(file2)
    
    nSig = len(Final_permCorrResults_df[Final_permCorrResults_df['real_p_Corr%sPerms' % nPermCorrUsed] <= 0.05])
    percSig = float(nSig) * 100 / len(Final_permCorrResults_df)
    
    print 'number of pairs with p<0.05 is %s (%s percent)' % (nSig, percSig)
    print 'number of corrected significant results:'
    print Full_permCorrResults_df[['Sig by bonferroni corrected pVal', 'sig. by FDR=0.1', 'sig. by FDR=0.25']].sum()
    
    return Final_permCorrResults_df, nPermCorrUsed


#-----------------------------------------------------------------

'''
the following function are sub-function for running permFisher test:

'''
def fisher_SeqSpecies_forPermutationTest(speciesList, sequenceList, species_sequences_df, isPerm):
    FisherResultsList = []
    countTests = 0
    for i, species in enumerate(speciesList):
        print i, species
        nSamplesSpecies = species_sequences_df[species].sum()
        if (nSamplesSpecies > 0) and (nSamplesSpecies < len(species_sequences_df)):  # check that the species appear in at                                                                   #least one sample and not all sample
            groups = species_sequences_df.groupby(species)
            for j, seq in enumerate(sequenceList):
                nSamplesSeq = species_sequences_df[seq].sum()
                if nSamplesSeq > 0:  # check that the sequence appear in at least one sample:
                    df4 = pd.DataFrame(index=['seq_abs', 'seq_pres'], columns=['spec_abs', 'spec_pres'])
                    for name, group in groups:
                        if name == 0:
                            seq_pres = group[seq].sum()
                            df4.loc['seq_pres', 'spec_abs'] = seq_pres
                            df4.loc['seq_abs', 'spec_abs'] = len(group) - seq_pres
                        if name == 1:
                            seq_pres = group[seq].sum()
                            df4.loc['seq_pres', 'spec_pres'] = seq_pres
                            df4.loc['seq_abs', 'spec_pres'] = len(group) - seq_pres
                    df4 = df4.fillna(0)
    #                 if not df4.isnull().values.any():
    #                     countTests = countTests + 1
                    OR, p = fisher_exact(df4, alternative='two-sided')
                    result = pd.DataFrame()
                    result.loc[0, 'species'] = species
                    result.loc[0, 'seq'] = seq
                    result.loc[0, 'OR'] = OR
                    result.loc[0, 'log10OR'] = np.log10(OR)
                    
                    if not isPerm:
                        speciesNegSeqNeg = df4.loc['seq_abs', 'spec_abs']
                        speciesNegSeqPos = df4.loc['seq_pres', 'spec_abs']
                        speciesPosSeqNeg = df4.loc['seq_abs', 'spec_pres']
                        speciesPosSeqPos = df4.loc['seq_pres', 'spec_pres']
                        result.loc[0, 'speciesNegSeqNeg'] = speciesNegSeqNeg
                        result.loc[0, 'speciesNegSeqPos'] = speciesNegSeqPos
                        result.loc[0, 'speciesPosSeqNeg'] = speciesPosSeqNeg
                        result.loc[0, 'speciesPosSeqPos'] = speciesPosSeqPos
                    FisherResultsList.append(result)
                else:
                    print 'sequence %s (num %s) do not appear in any sample and therefore not used for calculations' % (seq, j)
        else:
            print 'species apears in %s samples and therefore not used for calculation' % nSamplesSpecies
    FisherResultsdf = pd.concat(FisherResultsList)
    return FisherResultsdf

def permutateDFandCalcFisher(nPermFisher, species_sequences_df, speciesList, sequenceList, permDFsFolder):
    print nPermFisher
    print 'permutation number %s for seqSpecies df...' % nPermFisher
    perm_species_sequences_df = species_sequences_df.apply(lambda x: x.sample(frac=1).values)
#     print perm_species_sequences_df.shape
    print 'Fisher test for permutated df...'
    perm_result_df = fisher_SeqSpecies_forPermutationTest(speciesList, sequenceList, perm_species_sequences_df, isPerm=True)
    print perm_result_df.shape
    permDFfile = '%s/perm%sDF' % (permDFsFolder, nPermFisher)
    perm_result_df.to_pickle(permDFfile)
    
    return perm_result_df


'''
the following function executed permFisherTest for seq-species dfs. it uses the former functions as subfunctions
see the documentation of the function 'calc_permCorrTest' for explanations on the input. 
'''   



def calc_permFisherTest(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies,
                      topNseqs, nPermFisher, firstPermFisher, permFisherTestFolder, runRealFisher):
       
    # ## (1) get combined df of seqs and species according to the requested thresholds:
    print 'preparing combined seq and species file...'
    species_sequences_df, speciesList, sequenceList = prepare_MbDF_and_TCRdf_for_analysis(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, NsharedSamplesForSpecies, NsharedSamplesForSeqs,
                                          topNspecies, topNseqs)
    
    # ## (2) convert species_sequences_df to binary:
    print 'converting species_sequences_df to binary...'
    for column in species_sequences_df.columns.values:
        species_sequences_df[column] = np.where(species_sequences_df[column] > 0, 1, 0)
    
        
    # ## (3) run Fisher test for all seq-species pairs in the real data.
    # ##get all results in one df
    
    if runRealFisher:
        print 'running Fisher test for all seq-species pairs in the real data'
        real_result_df = fisher_SeqSpecies_forPermutationTest(speciesList, sequenceList, species_sequences_df, isPerm=False)
        real_result_df = real_result_df.rename(columns={'log10OR':'log10OR_real', 'OR':'OR_real'})
        file1 = '%s/real_result_df' % (permFisherTestFolder)
        real_result_df.to_pickle(file1)
    else:
        print 'loading real R table...'
        file1 = '%s/real_result_df' % (permFisherTestFolder)
        real_result_df = pd.read_pickle(file1)
        
    
    # ## (4)permutate seq-species df nPerm times, and each time calculate full correlation table and save it - USE CLUSTER FOR CALCULATIONS!
    permFisherDFsFolder = '%s/permDFs' % permFisherTestFolder
    
    
    if nPermFisher > firstPermFisher:
        print 'permutations and Fisher testing... running on the cluster!!! might be very long!!!'
        print 'total number of permutationsL %s' % nPermFisher
        print 'start after permutation number %s' % firstPermFisher
               
        if not isdir(permFisherDFsFolder):
            makedirs(permFisherDFsFolder)
        
        sethandlers()
        os.chdir(permFisherDFsFolder)
    
        # send jobs to the cluster:
        with qp(jobname='permFisher', q=['himemint.q', 'himem7.q'], trds_def=1, mem_def='2G') as q:  # consider changing queue list!!
            q.startpermanentrun()
            tk = {}  # a dictionary storing arguments and results
            for i in range(nPermFisher + 1):
                if i > firstPermFisher:
                    print i
                    tk[i] = q.method(permutateDFandCalcFisher, (i, species_sequences_df, speciesList, sequenceList,
                                                              permFisherDFsFolder))
            tk = {k:q.waitforresult(v) for k, v in tk.iteritems()}
        print 'done'
    else:
        print 'no new permutation tables are generated...'
    
    # ## (5) summarize all permutation results:
    print 'summarizing all Fisher permutation results...'
    
    Full_permFisherResults_df = real_result_df
    # merge all permDFs to the real r table:
    filenames = [f for f in listdir(permFisherDFsFolder) if isfile(join(permFisherDFsFolder, f))]
    nPermUsed = len(filenames)
    print 'number of permutations used is %s' % nPermUsed
    for n, f in enumerate(filenames):
    #     if n<5:
            print n, f
            file1 = '%s/%s' % (permFisherDFsFolder, f)
            df = pd.read_pickle(file1)
            Full_permFisherResults_df = pd.merge(Full_permFisherResults_df, df, how='inner',
                                            left_on=('species', 'seq'), right_on=('species', 'seq'))
            nPer = f.split('DF')[0]
            Full_permFisherResults_df = Full_permFisherResults_df.rename(columns={'OR':'OR_%s' % nPer, 'log10OR':'log10OR_%s' % nPer})
            
    # calculate how many permutations have more extreme r value (this is the p-value, two-sided)
    print 'calculating real p values'
    columnsToApply = [column for column in Full_permFisherResults_df.columns.values if 'log10OR_' in column]
    Full_permFisherResults_df['real_p_Fisher%sPerms' % nPermUsed] = Full_permFisherResults_df[columnsToApply].apply(perc_above, axis=1)
    
    nTests = len(speciesList) * len(sequenceList)

    Full_permFisherResults_df = add_corrected_pValues(Full_permFisherResults_df, pValueColumn='real_p_Fisher%sPerms' % nPermUsed, nTests=nTests, FDR=0.01)
    Full_permFisherResults_df = add_corrected_pValues(Full_permFisherResults_df, pValueColumn='real_p_Fisher%sPerms' % nPermUsed, nTests=nTests, FDR=0.1)
    Full_permFisherResults_df = add_corrected_pValues(Full_permFisherResults_df, pValueColumn='real_p_Fisher%sPerms' % nPermUsed, nTests=nTests, FDR=0.25)
    Full_permFisherResults_df = Full_permFisherResults_df.sort_values(by='log10OR_real')
    
    columns_to_use = ['species', 'seq', 'log10OR_real', 'speciesNegSeqNeg', 'speciesNegSeqPos', 'speciesPosSeqNeg', 'speciesPosSeqPos',
       'real_p_Fisher%sPerms' % nPermUsed, 'Sig by bonferroni corrected pVal', 'sig. by FDR=0.01',
       'sig. by FDR=0.1']
    Final_permFisherResults_df = Full_permFisherResults_df[columns_to_use]
    file2 = '%s/Final_permFisherResults_%sPerms_df.xlsx' % (permFisherTestFolder, nPermUsed)
    Final_permFisherResults_df.to_excel(file2)
    
    nSig = len(Final_permFisherResults_df[Final_permFisherResults_df['real_p_Fisher%sPerms' % nPermUsed] <= 0.05])
    percSig = float(nSig) * 100 / len(Final_permFisherResults_df)
    
    print 'number of pairs with p<0.05 is %s (%s percent)' % (nSig, percSig)
    print 'number of corrected significant results:'
    print Full_permFisherResults_df[['Sig by bonferroni corrected pVal', 'sig. by FDR=0.01', 'sig. by FDR=0.1']].sum()
    
    return Final_permFisherResults_df, nPermUsed

#-------------------------------------------------------------------------------------
'''
the following function is an helper for the function 'calc_diversity_correlations_TCR_Mb' that follows

'''

def calc_diversity_features_for_public_matrix(df,threshold):
#     print 'calculate diversity features for df...'
    featureDF = pd.DataFrame(index=df.index)
    featureDF['nCom'] = df.apply(lambda x: len([i for i in x if i > 0]), axis=1)
    featureDF['nCom_stringest'] = df.apply(lambda x: len([i for i in x if i > threshold]), axis=1)
#     featureDF['meanRA0'] = df.apply(lambda x: np.mean(x), axis=1)
    featureDF['meanRA'] = df.apply(lambda x: np.mean([i for i in x if i > 0]), axis=1)
    featureDF['meanRA_stringent'] = df.apply(lambda x: np.mean([i for i in x if i > threshold]), axis=1)
    featureDF['meanLogRA'] = df.apply(lambda x: np.log10(np.mean([i for i in x if i > threshold])), axis=1)
    featureDF['maxRA'] = df.apply(np.max, axis=1)
    featureDF['maxRA_stringent'] = df.apply(lambda x: np.max([i for i in x if i > threshold]), axis=1)
    featureDF['max/meanRA'] = featureDF['maxRA']/featureDF['meanRA']
    featureDF['max/meanRA0'] = featureDF['maxRA']/featureDF['meanRA']
    featureDF['max/meanRA_stringent'] = featureDF['maxRA']/featureDF['meanRA_stringent']
    featureDF['simpson'] = df.apply(lambda x: simpson([int(i * 1000) for i in x if not np.isnan(i)]), axis=1)
#     featureDF['simpson_stringent'] = df.apply(lambda x: simpson([int(i * 1000) for i in x if  i > threshold]), axis=1)
    featureDF['berger_parker_d'] = df.apply(lambda x: berger_parker_d([int(i * 1000) for i in x if not np.isnan(i)]), axis=1)
#     featureDF['berger_parker_d_stringent'] = df.apply(lambda x: berger_parker_d([int(i * 1000) for i in x if i > threshold]), axis=1)
    featureDF['shannon'] = df.apply(lambda x: shannon([int(i * 1000) for i in x if not np.isnan(i)], base=2), axis=1)
#     featureDF['shannon_stringent'] = df.apply(lambda x: shannon([int(i * 1000) for i in x if i > threshold], base=2), axis=1)
#     print featureDF.head()
    
    print 'calculate features correlations...'
    for i in range(len(featureDF.columns)):
        for j in range(i + 1, len(featureDF.columns)):
            column1 = featureDF.columns[i]
            column2 = featureDF.columns[j]
            r, p = MyPearsonr(featureDF[column1], featureDF[column2])
            print column1, column2, r,p
            
    return featureDF

'''
the following function takes MbDF and TCRdf (shared species or shared sequences matrices
and calculates for each several diversity measures per sample - number of shared sequences, mean relative abundance
(RA), max RA, shannon, simpson and berger_parker.
than it calculates the correlation between TCR and MB for this parameters.

input:
datasetFolder- string
MbDF_RA - MPA table (samples Xspecies, can be other taxonomic level)
MbDFName - string. make sure it fits the MbDF
TCRdf_RA - sharing sequence matrix (sample X seqeunce) - the minNshared parameter indicates the number of samples
that should be shared in order to included in the df (typically 2 or 50)
TCRdfName - string. make sure it fits the TCRdf
thresholdList - list of quadriple tuples, representing NsharedSamplesForSpecies,NsharedSamplesForSeqs,
topNspecies,topNseqs (the first couple or the second couple should be None but not both!)
outlierSTD - number of STDs to define outliers in outliers removal (None=no outlier removal)
minNshared - number of samples to be shared by TCR sequences in order to be included in the TCRdf.MAKE SURE THIS
INPUT CORRESPOND TO THE TCRdf_RA and the TCRdfName used!!!

Usage example:
# top 50 and top 250, std=None/3

stdList=[None,3,5]
thresholdList=[(None,None,50,50),(None,None,250,250),(None,None,400,400)]

datasetFolder='%s/TCR_real_data/SubSampled15000data_rep2' %MyPath
datasetName='PNP515_ss15000_rep2'

f1 = '%s/MicrobiomeDataTables/FilteredAndMergedOnBD/MPA_s_AllSeqProjects_SampleListPNP515_filterGenotekTrue_filterMinimalRead9000000_meannSTDNonenMinSamplesNone_filterOutlierSampleFalse_filterSamePersonFalse' % MyPath
MbDF_RA = pd.read_pickle(f1)
filterGenotek = True
filterMinimalReads = 9000000
MbDFName = 'Mb_Swabs%sMinReads%s' % (filterGenotek, filterMinimalReads)

minNshared = 50  # add to function variables if necessary
file1 = '%s/sharingAnalysis/sharingMatrix_%s_minNshared%s_RA' % (datasetFolder, datasetName, minNshared)
TCRdf_RA = pd.read_pickle(file1)
TCRdfName = 'TCRminNshared%s' % minNshared

for outlierSTD in stdList:
    print 'outlierSTD=%s' %outlierSTD
    divCorrTCR_MB=calc_diversity_correlations_TCR_Mb(datasetFolder,MbDF_RA,MbDFName,TCRdf_RA,
                                                 TCRdfName,thresholdList,outlierSTD,minNshared)
'''


def calc_diversity_correlations_TCR_Mb(datasetFolder, MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, thresholdList,
                                      outlierSTD, minNshared):
    
    # (1) prepare MB file:
    # load BD_FD converter, to make sure we remove BD_FD columns from the MB file:
    f1 = '%s/Sample files/updatedBDandFDlists/final_BD_FD_converter_mergedOnDate_RCfolderAllSeqProjects_31052018' % MyPath
    BD_FD = pd.read_pickle(f1)

    try:
        MbDF_RA = MbDF_RA.set_index('BD')  #### very important!!
    except KeyError:
        print 'couldnt reset index to BD,probably BD is already the index'
        MbDF_RA = MbDF_RA


    # remove unnecessary columns:
    for column in BD_FD.columns.values:
        if column in MbDF_RA.columns.values:
            print 'dropping %s column' % column
            MbDF_RA = MbDF_RA.drop(column, axis=1)
                                       
    # ## (2) define cutoffs for TCR and Mb matrix:
    for t in thresholdList:
        NsharedSamplesForSpecies = t[0]
        NsharedSamplesForSeqs = t[1]
        topNspecies = t[2]
        topNseqs = t[3]
        print 'threshold parameters are:'
        print NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies, topNseqs
    
        # ## (3) prepare TCR and Mb tables according to the cutoffs:
        species_sequences_df, speciesList, sequenceList = prepare_MbDF_and_TCRdf_for_analysis(MbDF_RA, MbDFName, TCRdf_RA, TCRdfName, NsharedSamplesForSpecies, NsharedSamplesForSeqs,
                                              topNspecies, topNseqs)
        print 'species_sequences_df length is %s' % len(species_sequences_df)
        if outlierSTD is not None:
            print 'filtering outliers...'
            species_sequences_df_filtered = filter_outliers(df=species_sequences_df,
                                            outlierSTD=outlierSTD, columnList=species_sequences_df.columns.values)
        else:
            species_sequences_df_filtered = species_sequences_df
        print 'species_sequences_df length is %s' % len(species_sequences_df_filtered)

        # get seperate Mb table:
        species_df = species_sequences_df_filtered[speciesList]
        print 'species_df shape is %s_%s' % (species_df.shape[0], species_df.shape[1])
    #     print species_df.index[:5]
    #     print species_df.index[-5:]
    #     # print species_df.head()

        # get seperate TCR table:
        seq_df = species_sequences_df_filtered[sequenceList]
        print 'seq_df shape is %s_%s' % (seq_df.shape[0], seq_df.shape[1])
    #     print seq_df.index[:5]
    #     print seq_df.index[-5:]
    #     # print species_df.head()

        # ## (4) calculate features for both dfs:
        print 'calculating features for seqDF...'
        seqFeatureDF = calc_diversity_features_for_public_matrix(seq_df)
        print 'calculating features for speciesDF...'
        speciesFeatureDF = calc_diversity_features_for_public_matrix(species_df)

        # ## (5) calculate correlations between TCR and Mb:
        print 'calculating correlation between TCR and MB...'
        divCorrTCR_MB = pd.DataFrame()
        for n, column in enumerate(seqFeatureDF.columns):
            r, p = MyPearsonr(seqFeatureDF[column], speciesFeatureDF[column])
            print column, r, p
            divCorrTCR_MB.loc[n, 'column'] = column
            divCorrTCR_MB.loc[n, 'r'] = r
            divCorrTCR_MB.loc[n, 'p'] = p

        # ## save diversity correlation df....
        minP = round(divCorrTCR_MB['p'].min(), 3)
        print 'the best p-value is %s' % minP
        print 'saving file...'
        parameters = 't(%s%s%s%s)_std%s_minNshared%s' % (NsharedSamplesForSpecies, NsharedSamplesForSeqs, topNspecies,
                                         topNseqs, outlierSTD, minNshared)
        folder = '%s/TCR_mb_results/divCorrTCR_MB_dfs' % datasetFolder
        if not isdir(folder):
            makedirs(folder)
        file1 = '%s/divCorrTCR_MB_%s_minP%s.xlsx' % (folder, parameters, minP)
        divCorrTCR_MB.to_excel(file1)
                                       
    return divCorrTCR_MB
    

#-------------------------------------------------------------
'''
the following function takes seqxspecies r table and cut it to include only specific number of top sequences (topNseqs) and
species (topNspecies)

usage example:
# cut to top 100x100
Rtable=Rtable250250
topNseqs=251
topNspecies=260

Rtable100100,topNseqs,topNspecies=cutSeqSpeciesRtable(Rtable,topNseqs,topNspecies)

'''

def cutSeqSpeciesRtable(Rtable,topNseqs,topNspecies):
    if len(Rtable.drop_duplicates(subset='species'))<topNspecies:
        topNspecies=len(Rtable.drop_duplicates(subset='species'))
        print 'not enough species in table, topNspecies was set to %s' %topNspecies
    if len(Rtable.drop_duplicates(subset='seq'))<topNseqs:
        topNseqs=len(Rtable.drop_duplicates(subset='seq'))
        print 'not enough seqs in table, topNseqs was set to %s' %topNseqs
    
    topSpecies=list(Rtable.drop_duplicates(subset='species').sort_values(by='nSamplesSpecies', ascending=False)[:topNspecies]['species'])
    topSeqs=list(Rtable.drop_duplicates(subset='seq').sort_values(by='nSamplesSeq', ascending=False)[:topNseqs]['seq'])
    cutRtable=Rtable[(Rtable['species'].isin(topSpecies))&(Rtable['seq'].isin(topSeqs))]
    print len(cutRtable)
    return cutRtable,topNseqs,topNspecies


#-------------------------

'''
the following function takes speciesxseqeunces Rtable, calculate the distances between all species, or all sequences
(depending on the indeCol and columnsColl chosen), and sort the pairs according to ascending distances

usage example:
pd.set_option('display.width', 1000)
distMetricList=['euclidean','sqeuclidean','canberra']


Rtable=Rtable100100
indexCol='species'
columnsCol='seq'
valueCol='r_real'

for metric in distMetricList:
    print metric
    sortedDistMat=calc_dist_Rtable(Rtable,indexCol,columnsCol,valueCol,metric)
    print sortedDistMat.head(12)

'''


def calc_dist_Rtable(Rtable,indexCol,columnsCol,valueCol,distMetric):
    #GENERATE PIVOT TABLE:
    RtablePivot=Rtable.pivot(index=indexCol, columns=columnsCol, values=valueCol)
    
    #generate ditance matrix:
    Rtable_distMat=pd.Series(pdist(RtablePivot.fillna(0),distMetric))
    
    #generate a df with all pairs and their distance, ranked by increasing distance:
    sample1List=[]
    sample2List=[]

    count=0
    for i in range(len(RtablePivot.index)):
        for j in range(i+1,len(RtablePivot.index)):
            sample1=RtablePivot.index[i]
            sample2=RtablePivot.index[j]
            sample1List.append(sample1)
            sample2List.append(sample2)
            count=count+1
#     print count
    df=pd.DataFrame({'sample1':sample1List,'sample2':sample2List,'dist':Rtable_distMat})
    sortedDistMat=df.sort_values(by='dist')
    
    return sortedDistMat



'''
the following function takes sorted distance table between species or sequences based on their correlations with sequences/species as defined 
in the former function. the table should include a binary score column.
then, it compares the distances between the groups defined by the binary score column and plot histograms/density plot

input:
ax=axes object
sortedDistMat - df generated by the function calc_dist_Rtable
binaryScoreCol- string, the name of the scoring column. the scoring should be 0,1
plotType - 'hist' or 'density'



usage example:
binaryScoreCol='sameGenus'
sortedDistMat=sortedDistMat
plotType='density'
distColumn - the name of the column containing the distance information

fig,ax=plt.subplots()
ax=plot_distribution_for_binary_scoring(ax,sortedDistMat,binaryScoreCol,plotType)
plt.show()


### see notebook: developPairDistanceScoring step 6 for more complex usage of this function####
'''

def plot_distribution_for_binary_scoring(ax,sortedDistMat,distColumn,binaryScoreCol,plotType):
    data={}

    for name,group in sortedDistMat.groupby(binaryScoreCol):
        if plotType=='hist':
            data[name]=list(group[distColumn])
        else:
            data[name]=group[distColumn]
    data1=data[0]
    data2=data[1]
    weights1=np.ones_like(data1)/len(data1)
    weights2=np.ones_like(data2)/len(data2)
    label1='%s=0' %binaryScoreCol
    label2='%s=1' %binaryScoreCol
    
    if plotType=='hist':
        plot=ax.hist((data1,data2), bins=50, color=('blue', 'red'), weights=[weights1,weights2],
                         label=(label1,label2), alpha=0.5)
    else:
        plot=sns.kdeplot(data1, shade=True, color="b",ax=ax,label=label1)
        plot=sns.kdeplot(data2, shade=True, color="r",ax=ax,label=label2)

    mean1=round(np.mean(data1),2)
    mean2=round(np.mean(data2),2)
    normDif=float(abs(mean1-mean2))/np.mean(mean1+mean2)
#     Alldata=data1+data2
#     Allweights=np.ones_like(Alldata)/len(Alldata)


    ks_s_cohort1_cohort2, ks_p_cohort1_cohort2=stats.ks_2samp(data1,data2)
    t_s_cohort1_cohort2, t_p_cohort1_cohort2=stats.ttest_ind(data1,data2)

    ax.annotate('KS_p=%s\nttest_p=%s\n%s mean=%s\n%s mean=%s\nnormDif=%s' 
                %(round(ks_p_cohort1_cohort2,6), round(t_p_cohort1_cohort2,6),label1,round(mean1,3),label2,round(mean2,3),
                 round(normDif,3)),
                    xy=(0.96, 0.95), xycoords='axes fraction', fontsize=9, horizontalalignment='right', 
                verticalalignment='top', fontweight='bold')

    return ax,mean1,mean2,ks_p_cohort1_cohort2,t_p_cohort1_cohort2,normDif



#----------------------------------------------------------------
'''
The following function enables relatively efficient loading of MBDF/tcrdf - BINARY AND RA
usage example:
#define parameters:
datasetFolder='%s/TCR_real_data/SubSampled15000data_rep2' %MyPath
datasetName='PNP515_ss15000_rep2'
minNshared=2

#load files:
TCRdf_RA,TCRdf_binary=prepare_TCRdfRA_TCRdfbinary_for_distMat(datasetFolder,datasetName,minNshared)

'''
def prepare_TCRdfRA_TCRdfbinary_for_distMat(datasetFolder,datasetName,minNshared,onlyProductive=False):
       
    #load/generate TCRdf_RA file:
    if onlyProductive:
        print 'taking only productive'
        TCRdf_RAFile='%s/sharingAnalysis/sharingMatrix_%s_minNshared%s_RA_onlyProductive%s' %(datasetFolder,datasetName,
                                                                         minNshared,onlyProductive)    
    else:
        try:
            TCRdf_RAFile='%s/sharingAnalysis/sharingMatrix_%s_minNshared%s_RA' %(datasetFolder,datasetName,
                                                                         minNshared)
        except:
            TCRdf_RAFile='%s/sharingAnalysis/sharingMatrix_%s_minNshared%s_RA_onlyProductive%s' %(datasetFolder,datasetName,
                                                                         minNshared,onlyProductive)
    print TCRdf_RAFile
    if isfile(TCRdf_RAFile):
        print 'loading TCRdf_RA...'
        TCRdf_RA=pd.read_pickle(TCRdf_RAFile)
    else:
        minNsharedList=list(minNshared)
        print 'preparing sharing Matrix'
        TCRdf_RA=gen_sharingMatrix_perDataset(datasetName,datasetFolder,minNsharedList,extractUniqueAA=False,getSharingStatistics=True)
    print 'done...'
#     print TCRdf_RA.iloc[:5,:5]
    
    #load/generate TCRdf_binary file:
    TCRdfFile_binary='%s/sharingAnalysis/sharingMatrix_%s_minNshared%s_binary' %(datasetFolder,datasetName,
                                                                         minNshared)
    if isfile(TCRdfFile_binary):
        print 'loading TCRdf_binary...'
        TCRdf_binary=pd.read_pickle(TCRdfFile_binary)
    else:
        print 'converting TCRdf_RA to binary...'
        TCRdf_binary = pd.DataFrame()
        print 'number of columns to convert is %s' %len(TCRdf_binary.columns.values)
        for n,column in enumerate(TCRdf_RA.columns.values):
            if n%1000==0:
                print n
            TCRdf_binary[column] = np.where(TCRdf_RA[column] > 0, 1, 0)
        TCRdf_binary = TCRdf_binary.set_index(TCRdf_RA.index)
        TCRdf_binary.to_pickle(TCRdfFile_binary)

    print 'done'
#     print TCRdf_binary.iloc[:5,:5]
    
    #check same-person couples in dataset:
    print 'checking for same-person pairs within the dataset:'
    samePairs=[('BD617','BD838'),('BD705','BD714')]

    sampleList=TCRdf_RA.index
    print 'number of samples in the dataset - %s' %len(sampleList)

    for samplePair in samePairs:
        if (samplePair[0] in sampleList) and  (samplePair[1] in sampleList):
            print 'the sample pair %s and %s are in the dataset' %(samplePair[0],samplePair[1])
        elif samplePair[0] in sampleList:
            print 'only sample %s is in the dataset' %samplePair[0]
        elif samplePair[1] in sampleList:
            print 'only sample %s is in the dataset' %samplePair[1]
            
    return TCRdf_RA,TCRdf_binary



'''
USAGE EXAMPLE FOR THE FOLLOWING FUNCTION:
libPrepMethod='Nextera'

MbDF_RA_libNone,MbDF_binary_libNone=prepare_MbDFRA_MbDFbinary_for_distMat(libPrepMethod)

'''

def prepare_MbDFRA_MbDFbinary_for_distMat(mbLevel,libPrepMethod):
    
    #load/generate MbDF_RA file:
    MbFolder='%s/MicrobiomeDataTables/FilteredAndMergedOnBD' %MyPath
    MbDF_RAFile='%s/MPA_%s_AllSeqProjects_SampleListPNP515_filterGenotekTrue_filterMinimalRead9000000_\
libMeth%s_meannSTDNonenMinSamplesNone_fOutlierSampleFalse_fSamePersonFalse' %(MbFolder,mbLevel,libPrepMethod)
    
    if isfile(MbDF_RAFile):
        print 'loading MbDF_RA...'
        MbDF_RA=pd.read_pickle(MbDF_RAFile)
#         print MbDF_RA.iloc[:5,:5]
    else:
         #load PNP515 sample list:
        with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP515','rb') as fp:
            PNP515=pickle.load(fp)
        SampleList=PNP515
        SampleListName='PNP515'
        filterGenotek=True
        filterMinimalReads=9000000
        filterlibPrepMethod=libPrepMethod
        groupFunction='noOutlierMean'
        nSTD=5
        nMinSamples=3
        ignoreNotSameYear=True
        mbDataFolder='AllSeqProjects'
        
        
        print 'MbDF file doesnt exist, GENERATING MbDF '
        genFilterMergeMbData(mbLevel,mbDataFolder,SampleList,SampleListName,filterGenotek,filterMinimalReads,filterlibPrepMethod,
                        groupFunction,nSTD,nMinSamples,ignoreNotSameYear)   
        
    print 'done...'
   
    
    
   
    #load/generate MbDF_binary file:
    MbDFFile_binary='%s_RA' %MbDF_RAFile
    
    if isfile(MbDFFile_binary):
        print 'loading MbDF_binary...'
        MbDF_binary=pd.read_pickle(MbDFFile_binary)
    else:
        print 'converting MbDF_RA to binary...'
        MbDF_binary = pd.DataFrame()
        print 'number of columns to convert is %s' %len(MbDF_binary.columns.values)
        for n,column in enumerate(MbDF_RA.columns.values):
            if n%1000==0:
                print n
            MbDF_binary[column] = np.where(MbDF_RA[column] > 0, 1, 0)
        MbDF_binary = MbDF_binary.set_index(MbDF_RA.index)
        MbDF_binary.to_pickle(MbDFFile_binary)

    print 'done'
    print MbDF_binary.iloc[:5,:5]
    
    #filter to include only species that appear in more than one sample:
    print 'filtering species df according to minimal number of shared samples by species=2'
    speciesCount = MbDF_binary.sum()
    selectedSpecies = list(speciesCount[speciesCount > 1].index)
    print len(selectedSpecies)
    MbDF_RA = MbDF_RA[selectedSpecies]
    MbDF_binary=MbDF_binary[selectedSpecies]
    
    
    #check same-person couples in dataset:
   #check same-person couples in dataset:
    print 'checking for same-person pairs within the dataset:'
    samePairs=[('BD617','BD838'),('BD705','BD714')]

    sampleList=MbDF_RA.index
    print 'number of samples in the dataset - %s' %len(sampleList)

    for samplePair in samePairs:
        if (samplePair[0] in sampleList) and  (samplePair[1] in sampleList):
            print 'the sample pair %s and %s are in the dataset' %(samplePair[0],samplePair[1])
        elif samplePair[0] in sampleList:
            print 'only sample %s is in the dataset' %samplePair[0]
        elif samplePair[1] in sampleList:
            print 'only sample %s is in the dataset' %samplePair[1]
            
    return MbDF_RA,MbDF_binary


#-----------------------------------------------------
'''
the following function enables fast spearman correlation and/or fisher association for specific seq-species pairs
and yield sum of correlation tests, spearman tests and combined, with multiple comparison corrections

input:
datasetFolder
datasetName
permNum - int, number of desired permutations (usually add 1 to reach round number)
minNshared - definition of the TCR df used, default=2, change to 10 if using un-subsampled PNP515 dataset
libPrepMethod - definition of the MbDF used - None, 'Nextera' or 'TrueSeq'
corrTest - True/False - whetehr to execute the corrTest
FisherTest - True/False whetehr to execute the Fisher test
mergeResults - True/False
PairNumberRange - None or tuple of min and max sample numbers

*** required improvements - enable merging without executing the tests (loading), define more felxible 
definition of seq-species pairs, and of TCRdf and MbDF used***

USAGE EXAMPLE:
datasetFolder='%s/TCR_real_data/SubSampled15000data_rep2' %MyPath
datasetName='PNP515_ss15000_rep2'
permNum=2
minNshared=2
libPrepMethod=None
corrTest=True
FisherTest=True
mergeResults=True
PairNumberRange=None



permSumDf,permSumDf_Binary,mergedResults_withIdentity=fast_perm_association_test_for_SeqSpeciesPairs(datasetFolder,datasetName,permNum,
                                                    minNshared,libPrepMethod,corrTest,FisherTest,mergeResults)

'''

def fast_perm_association_test_for_SeqSpeciesPairs(datasetFolder,datasetName,permNum,minNshared=2,libPrepMethod=None,
                                                  corrTest=True,FisherTest=True,mergeResults=True,PairNumberRange=None):

    print PairNumberRange 

    # (1) generate top pair list:
    print 'generating top pairs list...'
    # interesting pairs:
    # from each test (spearman/fisher):
    # -take all pairs with sig. p_value after bonferroni correction
    # -take all pairs with sig. p_value after FDR correction with q=0.01
    # -if fewer than 100 pairs, complete to 100 with top p_values
    # -from the merged file, take top 100 by p_values multiplication
    # combine all and remove duplicates

    # top pairs from spearman:
    file1='%s/TCR_mb_results/CorrResults_spearman_Mb_SwabsTrueMinReads9000000_TCRminNshared2_Thresholds(NoneNoneNoneNone)\
_olSTDNone/corrResults_spearman.xlsx' %datasetFolder
        
    spearmanResults=pd.read_excel(file1)
    topCorrPairsDF=spearmanResults[(spearmanResults['Sig by bonferroni corrected pVal']==1) | (spearmanResults['sig. by FDR=0.01']==1)]
    if len(topCorrPairsDF)<100:
        topCorrPairsDF=spearmanResults.sort_values(by='p').head(100)
    topCorrPairs=list(zip(list(topCorrPairsDF['species']),list(topCorrPairsDF['seq'])))
    print 'number of topCorrPairs is %s' %len(topCorrPairs)

    # top pairs from Fisher:
    file2='%s/TCR_mb_results/FisherResults_Mb_SwabsTrueMinReads9000000_TCRminNshared2_NsharedSpeciesNone_NsharedseqsNone_topSpeciesNone_topSeqsNone/\
FisherResults.xlsx' %datasetFolder
    FisherResults=pd.read_excel(file2)
    topFisherPairsDF=FisherResults[(FisherResults['Sig by bonferroni corrected pVal']==1) | (FisherResults['sig. by FDR=0.01']==1)]
    if len(topFisherPairsDF)<100:
        topFisherPairsDF=spearmanResults.sort_values(by='p').head(100)
    topFisherPairs=list(zip(list(topFisherPairsDF['species']),list(topFisherPairsDF['seq'])))
    print 'number of topFisherPairs is %s' %len(topFisherPairs)

    #top merged pairs:
    file3='%s/TCR_mb_results/mergedFisherCorr_spearmanPermFalse_0.0CorrPerms0.000000isherPerms_t(NoneNoneNoneNone)_stdNone.xlsx' %datasetFolder
    merged=pd.read_excel(file3)
    merged['p_mult']=merged['p_correlation']*merged['p_Fisher']
    topMerged=merged.sort_values(by='p_mult').head(100)
    topMergedPairs=list(zip(list(topMerged['species']),list(topMerged['seq'])))
    print 'number of topMergedPairs is %s' %len(topMergedPairs)

    # merge all top pairs and drop duplicates:
    AllTopPairs=topCorrPairs+topFisherPairs+topMergedPairs
    AllTopPairsUnique=list(set(AllTopPairs))
    print 'TOTAL number  of UNIQUE top pairs for the analysis is %s' %len(AllTopPairsUnique)

    # (2) load relevant TCRdf, Mbdf:
    print 'loading TCRdf and MbDF files...'
    #load files:
    TCRdf_RA,TCRdf_binary=prepare_TCRdfRA_TCRdfbinary_for_distMat(datasetFolder,datasetName,minNshared)
    MbDF_RA,MbDF_binary=prepare_MbDFRA_MbDFbinary_for_distMat(libPrepMethod)

    
    print 'executing Corr test=%s' %corrTest
    print 'executing Fisher test=%s' %FisherTest
    
    TCRdf_RA=filterSamePerson(TCRdf_RA,[0])
    MbDF_RA=filterSamePerson(MbDF_RA,[0])
    nTestsCorr=0
    nTestsFisher=0
    permSumDf=pd.DataFrame()
    permSumDf_Binary=pd.DataFrame()
    mergedResults_withIdentity=pd.DataFrame()
    
    
    if PairNumberRange is None:
        minPair=0
        maxPair=len(AllTopPairsUnique)
        FolderToSave='%s/TCR_mb_results/' %datasetFolder
    else:
        minPair=PairNumberRange[0]
        maxPair=PairNumberRange[1]
        FolderToSave='%s/TCR_mb_results/quickPermutDFs' %datasetFolder
        if not isdir(FolderToSave):
            makedirs(FolderToSave)

    for n,pair in enumerate(AllTopPairsUnique[minPair:maxPair]):
        species=pair[0]
        seq=pair[1]
        print n,pair
        real_r=np.nan
        real_p=np.nan
        try:
            pairData=pd.merge(pd.DataFrame(TCRdf_RA[seq]),pd.DataFrame(MbDF_RA[species]),how='inner', left_index=True,right_index=True)

        except:
            print 'couldnt get data for this pair'
            continue
        if corrTest:
            real_r_corr,real_p_corr=MySpearmanr(pairData[seq],pairData[species])
            print 'real r=%s, real p=%s' %(real_r_corr,real_p_corr)
        if FisherTest:
            # generate contingency table for this pair:
            if pairData[seq].sum()==len(pairData):
                print 'seq appears in all samples, drop this pair'
                continue
            elif pairData[species].sum()==len(pairData):
                print 'species appears in all samples, drop this pair'
                continue

            #prepare contingency table and calculate fisher for real data:
            tab = pd.crosstab(pairData[species] > 0.0001,pairData[seq] > 0)
            try:
                real_OR_fisher, real_p_fisher=fisher_exact(tab)
                print 'real OR=%s, real p=%s' %(real_OR_fisher,real_p_fisher)
            except:
                print 'couldnt execute test for this pair'
                continue



        moreSigCountCorr=1
        moreSigCountFisher=1
        totalCountCorr=1
        totalCountFisher=1


        for i in range(permNum):
            shuffled = pairData.apply(lambda x: x.sample(frac=1).values)

            if corrTest:
                perm_r_corr,perm_p_corr=MySpearmanr(shuffled[seq],shuffled[species])
                if abs(perm_r_corr)>=abs(real_r_corr):
                    moreSigCountCorr=moreSigCountCorr+1
                totalCountCorr=totalCountCorr+1

            if FisherTest:           
                tab_shuffled = pd.crosstab(shuffled[species] > 0,shuffled[seq] > 0)
                perm_OR_fisher,perm_p_fisher=fisher_exact(tab_shuffled)     

                if abs(np.log10(perm_OR_fisher))>=abs(np.log10(real_OR_fisher)):
                    moreSigCountFisher=moreSigCountFisher+1
                totalCountFisher=totalCountFisher+1

        if corrTest:       
            permP_corr=float(moreSigCountCorr)/totalCountCorr
            print 'permP_corr=%s' %permP_corr
            permSumDf.loc[n,'species']=species
            permSumDf.loc[n,'seq']=seq
            permSumDf.loc[n,'real_r_corr']=real_r_corr
            permSumDf.loc[n,'real_p_corr']=real_p_corr
            permSumDf.loc[n,'permP_corr']=permP_corr
            permSumDf.loc[n,'permNum']=permNum

            nTestsCorr=nTestsCorr+1

        if FisherTest: 
            permP_fisher=float(moreSigCountFisher)/totalCountFisher
            print 'permP_fisher=%s' %permP_fisher
            permSumDf_Binary.loc[n,'species']=species
            permSumDf_Binary.loc[n,'seq']=seq
            permSumDf_Binary.loc[n,'real_OR_fisher']=real_OR_fisher
            permSumDf_Binary.loc[n,'real_p_fisher']=real_p_fisher
            permSumDf_Binary.loc[n,'permP_fisher']=permP_fisher
            permSumDf_Binary.loc[n,'permNum']=permNum

            nTestsFisher=nTestsFisher+1
            
            
            
    if corrTest:
        permSumDf=add_corrected_pValues(permSumDf, pValueColumn='permP_corr', nTests=nTestsCorr, FDR=0.01)
        permSumDf=add_corrected_pValues(permSumDf, pValueColumn='permP_corr', nTests=nTestsCorr, FDR=0.25)
        permSumDf=permSumDf.sort_values(by='permP_corr')

        print 'number of expected tests is %s' %len(AllTopPairsUnique)
        print 'number of corr tests executed is %s' %nTestsCorr
        print permSumDf.head()
        file4='%s/SpearmanCorr_top%spairsquickPerms%s_min%smax%s.xls' %(FolderToSave,len(AllTopPairsUnique),permNum,
                                                                       minPair,maxPair)
        permSumDf.to_excel(file4)
        
    if FisherTest:
        permSumDf_Binary=add_corrected_pValues(permSumDf_Binary, pValueColumn='permP_fisher', nTests=nTestsFisher, FDR=0.01)
        permSumDf_Binary=add_corrected_pValues(permSumDf_Binary, pValueColumn='permP_fisher', nTests=nTestsFisher, FDR=0.25)
        permSumDf_Binary=permSumDf_Binary.sort_values(by='permP_fisher')

        print 'number of expected tests is %s' %len(AllTopPairsUnique)
        print 'number of corr tests executed is %s' %nTestsFisher
        print permSumDf_Binary.head()
        file5='%s/Fisher_top%spairsquickPerms%s_min%smax%s.xls' %(FolderToSave,len(AllTopPairsUnique),permNum,
                                                                       minPair,maxPair)
        permSumDf_Binary.to_excel(file5)
        print 'results were save in folder %s' %FolderToSave
        
    if mergeResults:
        
        mergedResults=pd.merge(permSumDf,permSumDf_Binary,how='outer',
                            left_on=['species','seq','permNum'],right_on=['species','seq','permNum'])
        mergedResults['comb_p']=mergedResults['permP_fisher']*mergedResults['permP_corr']
        mergedResults=mergedResults.sort_values(by='comb_p')
        
        # # add identity information to the table:
        
        # load TCR identity table:
        file6 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR CDR3 sequence databases/CDR3identityTable_23042018'
        CDR3identityTable = pd.read_pickle(file6)
        
        mergedResults_withIdentity = pd.merge(mergedResults, CDR3identityTable, how='left',
                                           left_on='seq', right_index=True)
        
        file7='%s/mergedResults_top%spairsquickPerms%s_min%smax%s.xls' %(FolderToSave,len(AllTopPairsUnique),permNum,
                                                                       minPair,maxPair)
        permSumDf_Binary.to_excel(file7)     
        

    print 'done '
    
    return permSumDf,permSumDf_Binary,mergedResults_withIdentity

##############################################################

def analyze_tcr_Mb_relations(MBdf,TCRdf_binary,nShared_perSpecies,nShared_perTCR,binary_zero=0):
    
    #this function calculate the associations between TCR sequenceses and microbial species. the TCRdf should be binary and the MBdf should not be binary.
    # Fisher test is calculated for binarized species data and mann-whitney for not-binarized data
    #binary_zero - the value that will be regarded as zero for binarizing the df. for Mb
    # data it is usually 0.0001
    ######### this function was written for analyzing Ravid's sample - check if it really generalizes ########
    
    #filter species and generate binary df:
    sample_count=MBdf.count()
    species_filtered=MBdf.loc[:,sample_count[(sample_count>=nShared_perSpecies)&(sample_count<len(MBdf))].index].fillna(0)
    print ('species shape after filtering only for species that appear in more than %s: ' %nShared_perSpecies,species_filtered.shape)
    print species_filtered.iloc[:4,:4]
    species_filtered_binary=(species_filtered>binary_zero).astype(int)
    print 'binary:'
    print species_filtered_binary.iloc[:4,:4]
    
    #filter TCRs:
    cols_to_use=[col for col in TCRdf_binary if (TCRdf_binary[col].sum()>=nShared_perTCR) & (TCRdf_binary[col].sum()<len(TCRdf_binary))]
    TCRdf_binary_filtered=TCRdf_binary[cols_to_use]
    print ('TCRdf shape after filtering only for seqs that appear in more than %s: ' %nShared_perTCR,TCRdf_binary_filtered.shape)
    print TCRdf_binary_filtered.iloc[:4,:4]
    
    #run over all seqs and species and calculate MW and fisher test:
    summary_df=pd.DataFrame()
    count=0
    for i, seq in enumerate(TCRdf_binary_filtered.columns):
        if i%10==0: print i
    #     print i,seq
        for j, s in enumerate(species_filtered.columns):
            
            #calculate mann-whitney:
            merged1=pd.merge(pd.DataFrame(TCRdf_binary_filtered[seq]),pd.DataFrame(species_filtered[s]),
                                            how='inner',left_index=True,right_index=True)
#             print 'merged1'
#             print merged1.head()
            data={}
            for name,group in merged1.groupby(seq):
                data[name]=group[s].tolist()
            MW_s, MW_p=mannwhitneyu(data[0],data[1])          
            mean_speciesFreq_seqAbs=np.mean(data[0])
            mean_speciesFreq_seqPres=np.mean(data[1])
            
            #calculate Fisher:
            merged2=pd.merge(pd.DataFrame(TCRdf_binary_filtered[seq]),pd.DataFrame(species_filtered_binary[s]),
                                            how='inner',left_index=True,right_index=True)
#             print 'merged2'
#             print merged2.head()
            if (merged2[s].sum()<8) | (merged2[seq].sum()<8):
                tab = pd.crosstab(merged2[seq],merged2[s], dropna=False).fillna(0)
                seqAbs_speciesAbs=tab.iloc[0,0]
                seqAbs_speciesPres=tab.iloc[0,1]
                seqPres_speciesAbs=tab.iloc[1,0]
                seqPres_speciesPres=tab.iloc[1,1]
                
    #             print (tab)
                try:
                    fisher_OR, fisher_p = fisher_exact(tab, alternative='two-sided')
                except:
                    print ('couldnt execute fisher test for seq %s and species %s' %(seq,s))
                    print (tab)
                    fisher_OR= 999999; fisher_p=1
            else:
                fisher_OR= 999999; fisher_p=1
              
            #update summary_df:
            summary_df.loc[count,'seq']=seq
            summary_df.loc[count,'species']=s
            summary_df.loc[count,'mean_speciesFreq_seqAbs']=mean_speciesFreq_seqAbs
            summary_df.loc[count,'mean_speciesFreq_seqPres']=mean_speciesFreq_seqPres
            summary_df.loc[count,'MW_s']=MW_s
            summary_df.loc[count,'MW_p']=MW_p
            summary_df.loc[count,'MW_s']=MW_s
            summary_df.loc[count,'MW_p']=MW_p
            summary_df.loc[count,'seqAbs_speciesAbs']=seqAbs_speciesAbs
            summary_df.loc[count,'seqAbs_speciesPres']=seqAbs_speciesPres
            summary_df.loc[count,'seqPres_speciesAbs']=seqPres_speciesAbs
            summary_df.loc[count,'seqPres_speciesPres']=seqPres_speciesPres
            summary_df.loc[count,'fisher_OR']=fisher_OR
            summary_df.loc[count,'fisher_p']=fisher_p

            #update count:
            count=count+1

    summary_df=add_corrected_pValues(summary_df,pValueColumn='MW_p',nTests=summary_df.shape[0],FDR=0.01)
    
    summary_df=add_corrected_pValues(summary_df,pValueColumn='fisher_p',nTests=summary_df.shape[0],FDR=0.01)
    summary_df=summary_df.sort_values(by='MW_p')


    return summary_df

    
