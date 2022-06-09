from os import listdir, mkdir
from os.path import isfile, join, isdir, exists
import pandas as pd
import numpy as np
from scipy import stats
import re
import matplotlib.pyplot as plt
from ShaniBA.myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot, draw_correlation_scatter
from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
import seaborn as sns
import random
from scipy.stats import pearsonr
from ShaniBA.MyFunctionsShani import *
import math
from ShaniBA.myplots import roundup, rounddown, find_decimal_fold
from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import *

# from tunneltoamazondb import getengine
from pandas import concat, read_csv, Series, read_sql

MyPath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF'

#---------------------------------------------------------------------------
'''
the following function generates a table linking BD numbers to FD numbers, UserIDs, regNums and SPID

input:
mergoOn='Year' or 'Date'
readCountFolder='Metabolon2'/'AllSeqProjects'

note that all other information is loaded from the database, and also my file for correcting regNum based on Daphna's function 
is necessary ('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/CleanCorrectBloodDNASamples')

UPDATED- 16/3/18 with new procedure to match unmatched BD samples based on UserID only. 
18/3/18 - added the option to select the DFout folder to get read counts info 

USAGE EXAMPLE:
final_BD_FD_converter_YEAR=gen_BD_FD_conversion_file(mergeOn='Year',readCountFolder='AllSeqProjects')

'''

def gen_BD_FD_conversion_file(mergeOn, readCountFolder):
    import time
    cdate = str(time.strftime("%d%m%Y"))
    cdate
    
    
    # # extract blood DNA information:
    BDsamples = read_sql('select  DnaID,UserID,storageDT from Lab.hostdna', getengine())
    BDsamples = BDsamples.rename(columns={'UserID':'oldUserID'})
    BDsamples['DnaID_num'] = BDsamples['DnaID'].str.replace('BD', '').astype(np.float)
    emptyBDs = BDsamples[BDsamples['DnaID'].isnull()]
    BDsamples = BDsamples[BDsamples['DnaID'].notnull()] 
    
    print 'number of rows without BD=%s-those rows were removed' % len(emptyBDs)

    countingBDs = pd.DataFrame(BDsamples['DnaID_num'].value_counts())  # counting repeating BDs
    repeatingBDs = countingBDs[countingBDs['DnaID_num'] > 1]
    BDsamples = BDsamples.drop_duplicates(subset='DnaID')  # removing BD duplicates
    
    UserIDs_RegNum = read_sql('select UserID,RegistrationCode from pnp.users', getengine())
    
    # merge regNums to BD samples:
    BDsamples_withOldIDs=pd.merge(BDsamples,UserIDs_RegNum,how='left',left_on='oldUserID',right_on='UserID')
    BDsamples_withOldIDs=BDsamples_withOldIDs.drop('UserID', axis=1)
    BDsamples_withOldIDs=BDsamples_withOldIDs.rename(columns={'RegistrationCode':'old_RegistrationCode'})
    print 'number of rows after addition of regNums is %s' %len(BDsamples_withOldIDs)
    
    from ShaniBA.SampleLists.correctingRegNumsFromDaphna import getQCedBloodRegNumByRegNum
    BDsamples_corrected=BDsamples_withOldIDs
    BDsamples_corrected['corrected_regCode']=BDsamples_corrected['old_RegistrationCode'].apply(lambda x:getQCedBloodRegNumByRegNum(x))
    
    BDsamples_corrected['correction status (reg number)']=np.where(BDsamples_corrected['corrected_regCode']!=BDsamples_corrected['old_RegistrationCode'],'corrected','same')
    BDsamples_corrected['correction status (reg number)']=np.where(BDsamples_corrected['corrected_regCode'].isnull(),'unknown',BDsamples_corrected['correction status (reg number)'])
    BDsamples_corrected['corrected_regCode']=np.where(BDsamples_corrected['corrected_regCode'].isnull(),BDsamples_corrected['old_RegistrationCode'],BDsamples_corrected['corrected_regCode'])
    
    # add correct UserID:
    BDsamples_corrected_withUserID=pd.merge(BDsamples_corrected,UserIDs_RegNum,how='left',left_on='corrected_regCode',right_on='RegistrationCode')
    print len(BDsamples_corrected_withUserID)
    
    BDsamples_corrected_withUserID=BDsamples_corrected_withUserID[~((BDsamples_corrected_withUserID['UserID']==11343)&(BDsamples_corrected_withUserID['corrected_regCode']==541391))]
    print 'sample number is still %s' % len(BDsamples_corrected_withUserID)
    
    BDsamples_corrected_withUserID = BDsamples_corrected_withUserID[['DnaID', 'storageDT', 'UserID', 'corrected_regCode', 'correction status (reg number)']]
    
    BDwithCorrectUserIDandRegs=BDsamples_corrected_withUserID
    
#     # # extract blood DNA information:
#     BDsamples = read_sql('select  DnaID,UserID,storageDT from Lab.hostdna', getengine())
#     BDsamples = BDsamples.rename(columns={'UserID':'oldUserID'})
#     BDsamples['DnaID_num'] = BDsamples['DnaID'].str.replace('BD', '').astype(np.float)
#     emptyBDs = BDsamples[BDsamples['DnaID'].isnull()]
# 
#     print 'extracted list of  %s BD samples from DB' % len(BDsamples)
#     print 'last BD index is %s' % BDsamples['DnaID_num'].max()
# 
#     BDsamples = BDsamples[BDsamples['DnaID'].notnull()]  # remove rows with no BD
#     print 'number of rows without BD=%s-those rows were removed' % len(emptyBDs)
# 
#     countingBDs = pd.DataFrame(BDsamples['DnaID_num'].value_counts())  # counting repeating BDs
#     repeatingBDs = countingBDs[countingBDs['DnaID_num'] > 1]
#     print repeatingBDs
#     BDsamples = BDsamples.drop_duplicates(subset='DnaID')  # removing BD duplicates
#     print 'BD samples repeating more than once=%s, duplicates were removed' % list(repeatingBDs.index)
# 
#     print 'final number of rows should be the initial list minus empty BD rows and minus duplicate BD rows\nand equall to the last BD index'
#     print 'final number of rows is %s' % len(BDsamples)
#     
#     # # correct UserIDs using correct RegNums:
#     # open the file with corrected userIDs for each BD number based on Daphna's correction for registration number:
#     file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/CleanCorrectBloodDNASamples'
#     CleanCorrectBloodDNASamples = pd.read_pickle(file1)
#     BD_correctReg = CleanCorrectBloodDNASamples[['DnaID', 'original registration code', 'correct registration code', 'correction status (reg number)']]
# 
#     
#     # merge BD list with correct RegNum list:
#     BDwithCorrectReg = pd.merge(BDsamples, BD_correctReg, how='left', left_on='DnaID', right_on='DnaID')
#     BDwithCorrectReg['correct registration code'] = BDwithCorrectReg['correct registration code'].fillna('unknown')
#     BDwithCorrectReg['correction status (reg number)'] = BDwithCorrectReg['correction status (reg number)'].fillna('unknown')
#     print 'original BD sample list length is %s' % len(BDsamples)
#     print 'sample list length after merging is %s (SHOULD BE THE SAME AS ABOVE!)' % len(BDwithCorrectReg)
# 
#     # extract UserID-regNum from the database
#     UserIDs_RegNum = read_sql('select UserID,RegistrationCode from pnp.users', getengine())
# 
#     # merge new UserID to BD sample according to the corrected regNum
#     BDwithCorrectUserID = pd.merge(BDwithCorrectReg, UserIDs_RegNum, how='left',
#                                  left_on='correct registration code',
#                                 right_on='RegistrationCode')
#     # give old UserIDs and regNum to samples that don't have correct regNum:
#     noCorrectReg = len(BDwithCorrectUserID[BDwithCorrectUserID['correct registration code'] == 'unknown'])
#     print '%s samples dont have correct registration code and are given old UserID' % noCorrectReg
#     BDwithCorrectUserID['UserID'] = np.where(BDwithCorrectUserID['correct registration code'] == 'unknown',
#                                       BDwithCorrectUserID['oldUserID'], BDwithCorrectUserID['UserID'])
# 
#     BDwithCorrectUserID=BDwithCorrectUserID[~((BDwithCorrectUserID['UserID']==11343)&(BDwithCorrectUserID['correct registration code']==541391))]
#     print 'sample number is still %s' % len(BDwithCorrectUserID)
# 
#     BDwithCorrectUserID = BDwithCorrectUserID[['DnaID', 'storageDT', 'UserID', 'correct registration code', 'correction status (reg number)']]
# 
#     # add regNum according to userID and use it only for samples with unknown regNum:
#     BDwithCorrectUserIDandRegs = pd.merge(BDwithCorrectUserID, UserIDs_RegNum, how='left',
#                                        left_on='UserID', right_on='UserID')
#     BDwithCorrectUserIDandRegs['correct registration code'] = np.where(BDwithCorrectUserIDandRegs['correct registration code'] == 'unknown',
#                                       BDwithCorrectUserIDandRegs['RegistrationCode'], BDwithCorrectUserIDandRegs['correct registration code'])
# 
#     BDwithCorrectUserIDandRegs = BDwithCorrectUserIDandRegs.drop('RegistrationCode', axis=1)
    
    # # generate BD_year and date:
    BDwithCorrectUserIDandRegs['storageDT'] = BDwithCorrectUserIDandRegs['storageDT'].astype(str)
    BDwithCorrectUserIDandRegs['Blood_Year'] = BDwithCorrectUserIDandRegs['storageDT'].str.split('-').str[0]
    BDwithCorrectUserIDandRegs['Blood_Date'] = BDwithCorrectUserIDandRegs['storageDT'].str.split(' ').str[0]
    BDwithCorrectUserIDandRegs['UserID'] = BDwithCorrectUserIDandRegs['UserID'].astype(str)
    BDwithCorrectUserIDandRegs['UserID'] = BDwithCorrectUserIDandRegs['UserID'].str.split('.').str[0]
    BDwithCorrectUserIDandRegs['UserID_bYear'] = BDwithCorrectUserIDandRegs['UserID'].astype(str).str.cat(BDwithCorrectUserIDandRegs['Blood_Year'], sep='_')
    BDwithCorrectUserIDandRegs['UserID_bDate'] = BDwithCorrectUserIDandRegs['UserID'].astype(str).str.cat(BDwithCorrectUserIDandRegs['Blood_Date'], sep='_')

    # # extract FD table
    FDsamples = read_sql('select  DnaID,UserID,connectionID,storageDT,SPID,isGenotek,Blacklisted from Lab.dna', getengine())
    print 'extracted list of  %s FD samples from DB' % len(FDsamples)
    
    # # generate FD year and date:
    FDsamples['storageDT'] = FDsamples['storageDT'].astype(str)
    FDsamples['FD_Year'] = FDsamples['storageDT'].str.split('-').str[0]
    FDsamples['FD_Date'] = FDsamples['storageDT'].str.split(' ').str[0]
    FDsamples['UserID'] = FDsamples['UserID'].astype(str)
    FDsamples['UserID'] = FDsamples['UserID'].str.split('.').str[0]
    FDsamples['UserID_fYear'] = FDsamples['UserID'].astype(str).str.cat(FDsamples['FD_Year'], sep='_')
    FDsamples['UserID_fDate'] = FDsamples['UserID'].astype(str).str.cat(FDsamples['FD_Date'], sep='_')
    
    # ## add read number information to FD samples:
    # read count info can be taken from either AllSeqProjects (new) or Metabolon2 (old) folder:
    if readCountFolder == 'Metabolon2':
        f1 = '/net/mraid08/export/jafar/Microbiome/Analyses/Metabolon2/DFOut/ReadCountDF.dat'
    else:
        f1 = '/net/mraid08/export/jafar/Microbiome/Analyses/AllSeqProjects/DFOut/ReadCountSpidDF.dat'
    ReadCountDF = pd.read_pickle(f1)
    

    # check if there are FD samples repeating twice:
    maxRepeatsPerFD = ReadCountDF.reset_index()['FD'].value_counts().sort_values().max()
    print 'read count information for %s samples was loaded' % len(ReadCountDF)
    print 'maximal repeats per FD samples is %s - SHOULD BE ONE! IF NOT, NEED TO CHECK IT' % maxRepeatsPerFD
    # extract only important information to merge with FD table:
    ReadCountDFsmall = ReadCountDF.reset_index()[['FD', 'PostHGF', 'PostSubSamp']]
    # merge read counts info with FD table:
    FDsamples_withReads = pd.merge(FDsamples, ReadCountDFsmall, how='left', left_on='DnaID', right_on='FD')
    FDsamples_withReads = FDsamples_withReads.drop('FD', axis=1)
    
    ## add library preparation method:
    print 'adding library preparation method:'
    f1='%s/Sample files/all_fds_barcode_df.dat' %MyPath
    libPrep=pd.read_pickle(f1)
    FDsamples_withReads=pd.merge(FDsamples_withReads,libPrep,how='left', left_on='DnaID', right_index=True)
    FDsamples_withReads=FDsamples_withReads.rename(columns={'Barcode':'libPrepMethod'})
    
    #### try to merge unmerged samples based on SPID - CURRENTLY NOT IMPLEMENTED AS IT SEEMS TO GIVE WRONG INFO###
 
    
    
    # # save BD and FD files
    print 'saving BD and FD lists to excel files...'
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/BDfile_%s.xlsx' % cdate
    BDwithCorrectUserIDandRegs.to_excel(file1)

    file2 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/FDfile_RCfolder%s_%s.xlsx' % (readCountFolder, cdate) 
    FDsamples_withReads.to_excel(file2)
    
    
    ## MATCH BD and FD samples!
    print 'merging BD and FD samples...'
    
    BD_FD_converter = pd.merge(BDwithCorrectUserIDandRegs, FDsamples_withReads, how='left', left_on=['UserID', 'UserID_b%s' % mergeOn],
                        right_on=['UserID', 'UserID_f%s' % mergeOn])
    
    BD_FD_converter = BD_FD_converter.rename(columns={'DnaID_x':'BD', 'DnaID_y':'FD'})
    
    print 'after merging based on UserID and year, there are %s matches and %s BD samples without \
    a match' % (len(BD_FD_converter[BD_FD_converter['FD'].notnull()]),
          len(BD_FD_converter[BD_FD_converter['FD'].isnull()]))
    
    # ## merging unmatched BDs to FDs based on UserID only:
    # ## new strategy started on 16.3.18
    
    print 'merging unmatched samples based on UserID only...'
    BDwithNoFDs = BD_FD_converter[BD_FD_converter['FD'].isnull()]
    BDwithFDs = BD_FD_converter[BD_FD_converter['FD'].notnull()]


    converterLength = len(BD_FD_converter)
    noPostHGF = len(BD_FD_converter[BD_FD_converter['PostHGF'].notnull()])
    uniqueBDwithPostHGF = BD_FD_converter[BD_FD_converter['PostHGF'].notnull()]['BD'].nunique()
    print 'coverter length is %s' % converterLength
    print 'number of rows with noPostHGF info is %s' % noPostHGF
    print 'number of unique BDs with postHGF info is %s' % uniqueBDwithPostHGF
    
    BD_FD_converter['Comment'] = ''
    for BD in list(BDwithNoFDs['BD']):
        BDdf = BDwithNoFDs[BDwithNoFDs['BD'] == BD]
        UserID = list(BDdf.loc[:, 'UserID'])[0]
        n = BDdf.index
        BDtoMerge = BDwithCorrectUserIDandRegs[BDwithCorrectUserIDandRegs['DnaID'] == BD]
    #     print n.astype(int),BD, UserID
    #     print 'length of BDtoMerge is %s and should be 1' %len(BDtoMerge)
        correspondingFDdf = FDsamples_withReads[FDsamples_withReads['UserID'] == UserID]
        if len(correspondingFDdf) > 0:
            newMerge = pd.merge(BDtoMerge, correspondingFDdf, how='left',
                              left_on='UserID', right_on='UserID')
            newMerge['Comment'] = np.where(newMerge['UserID_b%s' % mergeOn]!=newMerge['UserID_f%s' % mergeOn],'Not the same %s!!' % mergeOn,'')
            newMerge = newMerge.rename(columns={'DnaID_x':'BD', 'DnaID_y':'FD'})
    #         print newMerge
    #         print 'number of corresponding FD is %s' %len(correspondingFDdf)
    #         print 'len of newMerge is %s' %len(newMerge)
    #         print 'len of converter before changes is %s' %len(BD_FD_converter)
            BD_FD_converter = BD_FD_converter.drop(n, axis=0)
            BD_FD_converter = pd.concat([BD_FD_converter, newMerge]) 
    #         print 'len of converter after changes is %s' %len(BD_FD_converter)
    #         print BD_FD_converter[BD_FD_converter['BD']==BD]
    
    converterLength = len(BD_FD_converter)
    noPostHGF = len(BD_FD_converter[BD_FD_converter['PostHGF'].isnull()])
    uniqueBDwithPostHGF = BD_FD_converter[BD_FD_converter['PostHGF'].notnull()]['BD'].nunique()
    print 'converter length is %s' % converterLength
    print 'number of rows with noPostHGF info after new merging is %s' % noPostHGF
    print 'number of unique BDs with postHGF info after new merging is %s' % uniqueBDwithPostHGF
    
    # ## remove rows with no FD number:
    print 'removing rows with no FD sample:'
    final_BD_FD_converter = BD_FD_converter[BD_FD_converter['FD'].notnull()]
    print 'n samples after removal=%s' % len(final_BD_FD_converter)
    
    
    ##
    print 'dropping duplicates'
    final_BD_FD_converter['SameYear'] = np.where(final_BD_FD_converter['Blood_Year']==final_BD_FD_converter['FD_Year'],1,0)
    final_BD_FD_converter['SameDate'] = np.where(final_BD_FD_converter['Blood_Date']==final_BD_FD_converter['FD_Date'],1,0)
    final_BD_FD_converter = final_BD_FD_converter.drop_duplicates(subset=['BD', 'Blood_Year', 'Blood_Date', 'FD', 'FD_Year', 'FD_Date', 'SPID', 'UserID',
                                            'corrected_regCode', 'correction status (reg number)', 'SameYear','SameDate', 'connectionID', 'isGenotek',
                                                 'Blacklisted', 'PostHGF', 'PostSubSamp','libPrepMethod']) 
    
    print 'n samples after removal=%s' % len(final_BD_FD_converter)
    
    ## check how many many-to-many pairs I have:
    final_BD_FD_converter['BDrepeats']= final_BD_FD_converter.groupby('BD')['BD'].transform('count')
    final_BD_FD_converter['FDrepeats']= final_BD_FD_converter.groupby('FD')['FD'].transform('count')
    final_BD_FD_converter['manyToMany']=np.where((final_BD_FD_converter['BDrepeats']>1)&(final_BD_FD_converter['FDrepeats']>1),1,0)
    countingBDs = pd.DataFrame(final_BD_FD_converter['BD'].value_counts())  # counting repeating BDs
    repeatingBDs = countingBDs[countingBDs['BD'] > 1]
    repeatingBDsDF=final_BD_FD_converter[final_BD_FD_converter['BD'].isin(repeatingBDs.index)]
#     print repeatingBDsDF
    countingFDsRepBDs=pd.DataFrame(repeatingBDsDF['FD'].value_counts())  # counting repeating BDs
    repeatingFDsWithinBDs = countingFDsRepBDs[countingFDsRepBDs['FD'] > 1]
    print ' n FDs that repeat severl times within the repeating BD: %s' %repeatingFDsWithinBDs.index.nunique()
    print 'n repeating BDs: %s' %repeatingBDs.index.nunique()
#     print 'number of repeating pairs of BD-FD is %s' %repeatingFDsWithinBDs['FD'].nunique()
#     print len(repeatingFDsWithinBDs)
#     print 'here are the repeating FDs-BDs' 
#     print final_BD_FD_converter[final_BD_FD_converter['FD'].isin(repeatingFDsWithinBDs.index)][['BD','FD','Blood_Date','FD_Date']]
    
    
    # ## save relevant columns only
    final_BD_FD_converter = final_BD_FD_converter[['BD', 'Blood_Year', 'Blood_Date', 'Comment', 'FD', 'FD_Year', 'FD_Date', 'SPID', 'UserID',
                                            'corrected_regCode', 'correction status (reg number)', 'UserID_bYear','SameYear','SameDate', 'connectionID', 'isGenotek',
                                                 'Blacklisted', 'PostHGF', 'PostSubSamp','libPrepMethod','BDrepeats','FDrepeats','manyToMany']]    
    print 'I checked multyBDperFD pairs in file final_BD_FD_converter_mergedOnYear_07032018 (see section 4 below) and saw that trying to resolve those pairs by using dates instead of years was not proved helpful. so i will leave these pairs as is.'
    # # add StudyType
    # get studyTypeID and UserID information:
    studyTypeSPIDs = read_sql('select  * from pnp.spid', getengine())
    
    studyTypeSPIDs['SPID_Year'] = studyTypeSPIDs['Start_date'].astype(str).str.split('-').str[0]
    studyTypeSPIDs['UserID'] = studyTypeSPIDs['UserID'].astype(str)
    studyTypeSPIDs['UserID'] = studyTypeSPIDs['UserID'].str.split('.').str[0]
    studyTypeSPIDs['UserID_SPIDYear'] = studyTypeSPIDs['UserID'].astype(str).str.cat(studyTypeSPIDs['SPID_Year'], sep='_')

    studyTypeSPIDs_small = studyTypeSPIDs[['StudyTypeID', 'UserID_SPIDYear']]
    
    studyTypelists = studyTypeSPIDs_small.groupby('UserID_SPIDYear').agg(lambda x: x.unique().tolist())
    
    final_BD_FD_converter = pd.merge(final_BD_FD_converter, studyTypelists, how='left',
                                          left_on='UserID_bYear', right_index=True)
    
    # # add BD index column:
    final_BD_FD_converter['BD_index'] = final_BD_FD_converter['BD'].str.replace('BD', '')
    final_BD_FD_converter['BD_index']=final_BD_FD_converter['BD_index'].astype('int64')
    final_BD_FD_converter = final_BD_FD_converter.sort_values(by='BD')
    
    # # save final file:
    print 'final number of samples in the list is %s' % len(final_BD_FD_converter)

    print 'saving final files...'
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/final_BD_FD_converter_mergedOn%s_RCfolder%s_%s.xlsx' % (mergeOn,
                                                                                                                                       readCountFolder, cdate)
    final_BD_FD_converter.to_excel(file1)

    file2 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/final_BD_FD_converter_mergedOn%s_RCfolder%s_%s' % (mergeOn,
                                                                                                                                  readCountFolder, cdate)
    final_BD_FD_converter.to_pickle(file2)

   
    print 'Done!!!'
    
    return final_BD_FD_converter
#----------------------------------------------------------------------------------------------
'''
THE FOLLOWING FUNCTION WAS NOT CHECKED!
make sure the merging columns are indeed columns and not the index (change by 'reset_index()' if necessary
'''
def noOutlierMean(data, nSTD, nMinSamples):  # ## helping function for the following function
    data = np.array(data)
    mean = np.mean(data)
    std = np.std(data)
    if (std == 0) or (len(data) < nMinSamples):
        cleanMean = mean
    else:
        lowLim = mean - nSTD * std
        highLim = mean + nSTD * std
        data2 = data[data > lowLim]
        data2 = data2[data2 < highLim]
        cleanMean = np.mean(data2)
        if np.isnan(cleanMean):
            cleanMean = mean
    return cleanMean

    
def merge_phenotype_and_BDtable(phenotypeDF, phenotypeMergeColumn, BDtableMergeColumn, groupByBD, nSTD, nMinSamples):
    # load BD_FD conversion table:
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/final_BD_FD_converter'   
    final_BD_FD_converter = pd.read_pickle(file1)
    
    mergedDF = pd.merge(final_BD_FD_converter, phenotypeDF, how='left',
                  left_on=BDtableMergeColumn, right_on=phenotypeMergeColumn)
    
    print 'phenotype table length=%s' % len(phenotypeDF)
    print 'BD_FD table length=%s' % len(final_BD_FD_converter)
    print 'merged table length=%s' % len(mergedDF)
    
    if groupByBD:
        mergedDF = mergedDF.groupby('BD').agg(lambda x: noOutlierMean(x, nSTD, nMinSamples))
        print 'merged table length after grouping by BD and outlier removal=%s' % len(mergedDF)
        
    print 'phenotype table length=%s' % len(phenotypeDF)
    print 'BD_FD table length=%s' % len(final_BD_FD_converter)
    
        
    return mergedDF


#-----------------------------------------------------------
'''
 the following function can be used to take an excel file containing one column with BD sample names in the format of
 *BD* and convert it to pythonic list object with sample names in the format of 'BD*'
'''

def gen_py_BD_list(file_fullPath):
    BDdf = pd.read_excel(file_fullPath)
    BDlist = list(BDdf.iloc[:, 0].astype(str))
    regex_BD = re.compile('BD')
    BDlist_clean = ['BD' + x.split('BD')[1] for x in BDlist if re.search(regex_BD, x)]
    
    return BDlist_clean


#------------------------------------------------------------------------------------------------------------
'''
the following function takes sample data files generated by 'adaptive' and correct them (correct table headers and save as csv):

usage example:
datasetFolder='%s/TCR_real_data/RavidSamples' %MyPath
correctSampleDataFiles(datasetFolder)


'''
def correctSampleDataFiles(datasetFolder):
    #get correct column headers:
    f='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/SamplesForAnalysis/BD438.tsv'
    BD438=pd.read_table(f)
    right_column_names=BD438.columns.values
    newColumnList=right_column_names[:-3]

    #get file names to be checked:
    fullSamplesFolder='%s/SamplesForAnalysis' %datasetFolder
    files = [f for f in listdir(fullSamplesFolder) if isfile(join(fullSamplesFolder, f))]
    print 'number of samples in folder is %s' %len(files)

    #define the folder to save the corrected samples:
    newSamplesFolder='%s/SamplesForAnalysis_corrected' %datasetFolder

    correctedFiles=[f for f in listdir(newSamplesFolder) if isfile(join(newSamplesFolder, f))]
    print 'number of samples in the corrected folder is %s' %len(correctedFiles)

    #looping over all samples and correct them:

    for n,f in enumerate(files):
    #     if n<10:
            print n,f
            if f not in correctedFiles:
                df=pd.read_table('%s/%s' %(fullSamplesFolder,f))
                columns=df.columns.values
                if set(columns)==set(right_column_names):
                    print 'correct column headers'
                    print df['sequenceStatus'].head()
                else:
                    print 'before correction:'
                    print df['sequenceStatus'].head()
                    df=df.iloc[:,:44]
                    df.columns=newColumnList
                    df=df.rename(columns={'count (templates/reads)':'count (templates)','count (reads)':'count (templates)'})
                    print 'after correction:'
                    print df['sequenceStatus'].head()

                file1='%s/%s' %(newSamplesFolder,f)
                df.to_csv(file1,sep='\t')
            else:
                print 'sample %s already exist in corrected file folder' %f

    #check that the correct number of samples exist in the corrected folder:
    correctedFiles=[f for f in listdir(newSamplesFolder) if isfile(join(newSamplesFolder, f))]
    print 'number of samples in the corrected folder is %s' %len(correctedFiles)


 #-------------------------------------------------------------------------------------------------------------
'''
 the following function takes phenotype files for two cohorts and generate matched sub-cohorts.
 the function is currently completely not modular and fits specific cohorts and phenotypes but should be used as a basis to 
 generate more modular function. 
'''
 
 
 
 
def find_matched_PNP530_Cardio126(maximalAgeDif,maximalBMIDif,maximalCreatDif):
    
    #load PNP530 phenotype file:
    f1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/PNP530ss9000_AgeGenderBMIcreatSmoking.xlsx'
    PNP530ss9000_AgeGenderBMIcreatSmoking=pd.read_excel(f1)
    
    #load cardio126 phenotype file:
    f1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/Cardio126ss9000_AgeGenderBMIcreatSmoking.xlsx'
    Cardio126ss9000_AgeGenderBMIcreatSmoking=pd.read_excel(f1)

    #define dfs for sampling and new dfs for matched samples:
    CardiomatchedSamples=pd.DataFrame()
    PNPmatchedSamples=pd.DataFrame()
    PNP=PNP530ss9000_AgeGenderBMIcreatSmoking.set_index('BD')
    Cardio=Cardio126ss9000_AgeGenderBMIcreatSmoking.set_index('BD')



    #loop over all samples in cardio list to find matched PNP samples 
    for n,sample in enumerate(Cardio.index):
    #     if n<3:
            print n,sample

            #get sample phenotypes:
            cardioSample=sample
            cardioAge=Cardio.loc[sample,'Age']
            cardioGender=Cardio.loc[sample,'Gender']
            cardioBMI=Cardio.loc[sample,'BMI']
            cardioCreatinine=Cardio.loc[sample,'Creatinine']
            cardioSmoking=Cardio.loc[sample,'Smoking']

            CardioSample=pd.DataFrame(Cardio.loc[sample,:])
            print CardioSample

            #search for similar sample in cardio and add to 'matched Samples' df

            ##same gender:
#             print 'PNP length is %s' %len(PNP)
            potentialPNP=PNP[PNP['Gender']==cardioGender]
            if len(potentialPNP)<1:
                print 'no more %s samples in PNP cohort' %cardioGender
                continue

            ## similar age:
            minAgeDif=abs(cardioAge-potentialPNP['Age']).min()
            print 'minAgeDif=%s' %minAgeDif
            if (maximalAgeDif is not None) & (minAgeDif>maximalAgeDif):
                print 'no more samples with similar age' 
                continue
            if np.isnan(minAgeDif):
                print 'no more samples with similar age' 
                continue

            potPNPsimAge=potentialPNP[(potentialPNP['Age']==cardioAge-minAgeDif)|(potentialPNP['Age']==cardioAge+minAgeDif)]
#             print potPNPsimAge

            ## similar BMI:
            minBMIdif=abs(cardioBMI-potPNPsimAge['BMI']).min()
            print 'minBMIdif=%s' %minBMIdif
            if (maximalBMIDif is not None) & (minBMIdif>maximalBMIDif):
                print 'no more samples with similar BMI'
                continue
            elif (maximalBMIDif is None) & (np.isnan(minBMIdif)):
                potPNPsimAgesimBMI=potPNPsimAge
            elif (maximalBMIDif is not None) & (np.isnan(minBMIdif)):
                print 'no more samples with similar age' 
                continue
            elif len(potPNPsimAge)>1:
                potPNPsimAgesimBMI=potPNPsimAge[(potPNPsimAge['BMI']==cardioBMI-minBMIdif)|(potPNPsimAge['BMI']==cardioBMI+minBMIdif)] #### CORRECT HERE!!!!
            else:
                potPNPsimAgesimBMI=potPNPsimAge

            ## similar Creatinine:
            minCreatdif=abs(cardioCreatinine-potPNPsimAgesimBMI['Creatinine']).min()
            print 'minCreatdif=%s' %minCreatdif
            if (maximalCreatDif is not None) & (minCreatdif>maximalCreatDif) :
                print 'no more samples with similar Creatinine' 
                continue 
            
            if np.isnan(minCreatdif):
                potPNPsimAgesimBMIsimCreat=potPNPsimAgesimBMI 

            if len(potPNPsimAgesimBMI)>1:

                potPNPsimAgesimBMIsimCreat=potPNPsimAgesimBMI[(potPNPsimAgesimBMI['Creatinine']==cardioCreatinine-minCreatdif)|(potPNPsimAgesimBMI['Creatinine']==cardioCreatinine+minCreatdif)] #### CORRECT HERE!!!!
            else:
                potPNPsimAgesimBMIsimCreat=potPNPsimAgesimBMI 

            print potPNPsimAgesimBMIsimCreat

            # 
            takeSample=pd.DataFrame(potPNPsimAgesimBMIsimCreat.iloc[0,:])
            sampleTaken=takeSample.columns[0]

            if n==0:
                PNPmatchedSamples=takeSample
                CardiomatchedSamples=CardioSample
            else:
                PNPmatchedSamples=pd.merge(PNPmatchedSamples, takeSample,how='outer', left_index=True,right_index=True)
                CardiomatchedSamples=pd.merge(CardiomatchedSamples, CardioSample,how='outer', left_index=True,right_index=True)

           
            print 'sample taken: %s' %sampleTaken

            print 'selected sample:'
            print takeSample

    #         print 'matched samples so far:'
    #         print PNPmatchedSamples

            #remove selected sample from PNP df
            PNP=PNP.drop(sampleTaken)


#             print 'PNP length is now %s' %len(PNP)


    print 'done matching'
    print 'n PNP samples = %s' %len(PNPmatchedSamples.columns)
    print 'n Cardio samples = %s' %len(CardiomatchedSamples.columns)

    PNPmatchedSamples=PNPmatchedSamples.T
    CardiomatchedSamples=CardiomatchedSamples.T

    # save matched Lists:
    f1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/PNPmatchedtoCardio_maxAgeDif%s_maxBMIDif%s.xlsx' %(maximalAgeDif,maximalBMIDif)
    PNPmatchedSamples.to_excel(f1)

    f2='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/CardiomatchedtoPNP_maxAgeDif%s_maxBMIDif%s.xlsx' %(maximalAgeDif,maximalBMIDif)
    CardiomatchedSamples.to_excel(f2)
    
    return PNPmatchedSamples,CardiomatchedSamples


        
#-------------------------------------------------------------------------------------------------------
def gen_sampleList_from_Folder(folder,toSaveName=None):
    filenames=[f for f in listdir(folder) if isfile(join(folder,f))]
    filenames=[f.split('.')[0] for f in filenames]
    sampleList=editSampleNamesList(filenames)
    
    if toSaveName is not None:
        with open('%s/Sample files/BD lists/%s' %(MyPath,toSaveName),'wb') as fp:
            pickle.dump(sampleList,fp)
    return sampleList


#------------------------------------------------------------------------------
'''
the following function takes a df whose index is BD sample list and check whether it contains relative samples
the output is a list of samples to drop (the relative samples)
need to define onlyBloodRels=True if would like to drop only relatives in blood

'''
def drop_relatives(df,onlyBloodRels=False):
    
    #df index should be BD numbers!
    
    f2='%s/Sample files/PNP530-relationships.xlsx' %MyPath
    PNP530_relationships=pd.read_excel(f2)
    samplesToDrop=[]
    
    if onlyBloodRels:
        RelativeTable=PNP530_relationships[(PNP530_relationships['is relative in PNP515?']==1) & (PNP530_relationships['BloodGroup'].notnull())].set_index('BD')
    else:
        RelativeTable=PNP530_relationships[PNP530_relationships['is relative in PNP515?']==1].set_index('BD')

    for sample in df.index:
    #     print (sample[0])
        if sample in RelativeTable.index.tolist() and sample not in samplesToDrop:
            sample2=RelativeTable.loc[sample,'Relative']
            if sample2 in df.index.tolist():
#                 print ('dropping' , sample, sample2)
                samplesToDrop.append(sample2)

#     print (samplesToDrop)
    return samplesToDrop
   
    
