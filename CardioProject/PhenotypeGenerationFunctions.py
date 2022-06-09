from os import listdir, mkdir
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
from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import * 
from ShaniBA.TCR_microbiome_interactions.TCR_microbiome_interactions_functions import *

# from SufficientStatistics import *
from ShaniBA.MyFunctionsShani import *
import math
from skbio.stats.distance import mantel
from scipy.spatial.distance import braycurtis, pdist
# from tunneltoamazondb import getengine
from pandas import concat, read_csv, Series, read_sql


#------------------------------------------------------



def extract_interestingQuestions_fromDB():
    quest = read_sql('select  UserID, QuestionID, Answer, Timestamp from pnp.questionnaire', getengine())
    print 'question list length is %s' % len(quest)
    # # extract question information:
    lookup = read_sql('select  QuestionID,Question from pnp.questions', getengine())
    print 'lookup list length is %s' % len(lookup)
    # merge:
    questWithLookup = pd.merge(quest, lookup, how='left', left_on='QuestionID', right_on='QuestionID')
    print len(questWithLookup)
    # define interesting questions:
    interestingQuestions = range(1005, 1039) + range(1083, 1092) + range(1680, 1689) + [1045]
    interestingQuestionsDF = pd.DataFrame(interestingQuestions)
    
    questWithLookup_Interesting = pd.merge(interestingQuestionsDF, questWithLookup, how='left', left_on=0, right_on='QuestionID')
    print 'list of interesting question length is %s' % len(questWithLookup_Interesting)
    
    # extract question year:
    for n in questWithLookup_Interesting.index:
        m = re.search('(?<=20)\w+', questWithLookup_Interesting.loc[n, 'Timestamp'])
        if m is not None:
            questWithLookup_Interesting.loc[n, 'year'] = '20' + m.group(0)
        
    questWithLookup_Interesting = questWithLookup_Interesting.drop(0, axis=1)
    questWithLookup_Interesting['UserID_qYear'] = questWithLookup_Interesting['UserID'].astype(str).str.cat\
    (questWithLookup_Interesting['year'], sep='_')
    
    # save file:
    print 'saving file...'
    file4 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/questWithLookup_Interesting.csv'
    questWithLookup_Interesting.to_csv(file4)
    
    print 'question divided by year:'
    print questWithLookup_Interesting['year'].value_counts()


    return questWithLookup_Interesting




#-------------------------------------------------------------------------

def getCorrectYear(x):  # helper function for extract_age_DB function
    try:
        y1 = str(x)[-2:]
        y2 = int(y1)
    except:
        try:
            y1 = str(int(float(x)))[-2:]
            y2 = int(y1)
        except:
            y2 = np.nan
    return y2


def extract_age_DB():
    # extract age information from DF:
    print 'extracting age info from DB'
    Age = read_sql('select  UserID, QuestionID, Answer,TimeStamp from pnp.questionnaire WHERE QuestionID=415 or QuestionID=517 ', getengine())
    print 'number of answers is %s' % len(Age)
    
    # extract answer's year:
    print 'extracting answers year'
    for n in Age.index:
        m = re.search('(?<=20)\w+', Age.loc[n, 'TimeStamp'])
        if m is not None:
            Age.loc[n, 'year'] = '20' + m.group(0)
    print 'years are:'
    print Age.year.unique().tolist()
    print Age.head()
       
    
    print 'extracting the correct age per User...'
    Age = Age.drop_duplicates(subset=['UserID', 'QuestionID'])
    Age2 = Age.pivot(index='UserID', columns='QuestionID', values='Answer')
    Age2['YOB'] = np.where(Age2[415].notnull(), Age2[415], Age2[517])
    Age2['YOB2'] = Age2['YOB'].apply(lambda x: getCorrectYear(x))
    Age2['Age'] = np.where(Age2['YOB2'] > 18, Age2.year - (Age2['YOB2'] + 1900), Age2.year - 2000 - Age2['YOB2'])
    
    
    print 'converting Users to BD samples...'
    
    Age3 = pd.DataFrame(Age2['Age'])
    
    # add BD:
    # load BD SAMPLE FILE:
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/BDfile_31052018.xlsx'          
    BDfile = pd.read_excel(f1)
    BDfile_small = BDfile[['DnaID', 'UserID']]
    BDfile_small = BDfile_small.drop_duplicates()
    BDfile_small = BDfile_small.rename(columns={'DnaID':'BD'})
    Age4 = pd.merge(BDfile_small, Age3, how='left', left_on='UserID', right_index=True)
    
    Age5 = Age4[Age4.Age.notnull()]
    
    print 'number of BD samples with age is %s' % len(Age5)

    print 'saving into excel file...'
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/PNPage.xlsx'          
    Age5.to_excel(f1)
    
    print 'done'
    
    return Age5

#-----------------------------------------------------------------------
def extract_Gender_DB():
    # extract Gender information from DF:
    print 'extracting Gender info from DB'
    GenderDF = read_sql('select  UserID, QuestionID, Answer from pnp.questionnaire WHERE QuestionID=414 or QuestionID=516 or QuestionID=2342 ', getengine())
    print 'number of answers is %s' % len(GenderDF)
       
    print 'extracting the correct answer per User...'
    GenderDF2 = GenderDF.drop_duplicates(subset=['UserID', 'QuestionID'])
    GenderDF3 = GenderDF2.pivot(index='UserID', columns='QuestionID', values='Answer')
    
    # CONVERT HEBREW WORDS TO ENGLISH:
    GenderDF3[GenderDF3 == u'\u05d6\u05db\u05e8'] = 'Male'
    GenderDF3[GenderDF3 == u'\u05e0\u05e7\u05d1\u05d4'] = 'Female'
    
    print 'check users with ambigiuos answers:...'
    GenderDF3['notNone'] = GenderDF3[[414, 516, 2342]].apply(lambda x: x.count(), axis=1)
    print  GenderDF3[GenderDF3['notNone'] > 1]
    
    print 'choose for each user the first value which is not None:'
    GenderDF3['Gender'] = GenderDF3.apply(lambda x: next(item for item in x if item is not None), axis=1)
    
    print 'converting Users to BD samples...'
    GenderDF4 = pd.DataFrame(GenderDF3['Gender'])

    # add BD:
    # load BD SAMPLE FILE:
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/BDfile_31052018.xlsx'          
    BDfile = pd.read_excel(f1)
    BDfile_small = BDfile[['DnaID', 'UserID']]
    BDfile_small = BDfile_small.drop_duplicates()
    BDfile_small = BDfile_small.rename(columns={'DnaID':'BD'})
    GenderDF5 = pd.merge(BDfile_small, GenderDF4, how='left', left_on='UserID', right_index=True)

    GenderDF6 = GenderDF5[GenderDF5.Gender.notnull()]

    print 'number of BD samples with Gender is %s' % len(GenderDF6)
    print 'number of unique BD samples is %s' % GenderDF6['BD'].nunique()    
    
    
    print 'saving into excel file...'
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/PNPgender.xlsx'          
    GenderDF6.to_excel(f1)
    
    print 'done'
    
    return GenderDF6
    
    
    
#--------------------------------------------------------------------
def extract_BMI_DB():
    
    print 'extracting BMI info from DB...'
    BMI = read_sql('select UserID, height, weight from user_connectionMeetingData', getengine())
    print len(BMI)
    
    print 'fixing nonsense height and weight data...'
    BMI2 = BMI.copy()
    BMI2.height = np.where((BMI2.height > 1) & (BMI2.height < 2), BMI2.height * 100, BMI2.height)  # all heights between 1 and 2 probably 
                                                                                        # reflects height in meters and not  cm
    BMI2.height = BMI2.height.mask(BMI2.height < 135)  # remove all heights below 140cm
    BMI2.height = BMI2.height.mask(BMI2.height > 210)  # remove all heights above 200cm
    BMI3 = BMI2.copy()
    BMI3.weight = BMI3.weight.mask(BMI3.weight < 38)  # remove all weights below 38kg
    BMI3.weight = BMI3.weight.mask(BMI3.weight > 200)  # remove all weights above 180kg
    
    print 'calculating BMI...'
    BMI4 = BMI3.copy()
    BMI4['BMI'] = BMI4['weight'] / (BMI4['height'] / 100) ** 2  # calculate BMI
    
    print 'adding BD info and converting to BD-based data'
# load BD SAMPLE FILE:
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/BDfile_31052018.xlsx'          
    BDfile = pd.read_excel(f1)
    BDfile_small = BDfile[['DnaID', 'UserID']]
    BDfile_small = BDfile_small.drop_duplicates()
    BDfile_small = BDfile_small.rename(columns={'DnaID':'BD'})
    BMI5 = pd.merge(BMI3, BDfile_small, how='left', left_on='UserID', right_on='UserID')
    print 'n samples after merging is %s' % len(BMI5)
    BMI5 = BMI5[BMI5.BD.notnull()]
    print 'n samples after nan removal is %s' % len(BMI5)
    
    from collections import Counter
    BMIgrouped = BMI5.groupby('BD').agg({'UserID':lambda x: Counter(x).most_common()[0][0],
                                      'height':lambda x: Counter(x).most_common()[0][0],
                                      'weight':lambda x: noOutlierMean(x, nSTD=3, nMinSamples=2)})
    print 'number of BD samples after grouping is %s' % len(BMIgrouped)
    BMIgrouped = BMIgrouped[BMIgrouped.weight.notnull() | BMIgrouped.height.notnull()]
    print 'number of BD samples with either weight or height info is %s' % len(BMIgrouped)
    
    # calculate BMI again, look for outliers
    BMIgrouped['BMI'] = BMIgrouped['weight'] / (BMIgrouped['height'] / 100) ** 2
    print 'number of BD samples with BMI info is %s' % len(BMIgrouped[BMIgrouped['BMI'].notnull()])
    
    print 'remove outliers...'
    BMIgroupedOLremoved = filter_outliers(df=BMIgrouped, outlierSTD=4, trim=False)
    
    print 'BMIgrouped describe:'
    print BMIgrouped.describe()
    
    print 'BMIgroupedOLremoved describe:'
    print BMIgroupedOLremoved.describe()
    
    print 'saving files...'
    # save with and without outlier removal
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/PNP_BMI.xlsx'          
    BMIgrouped.to_excel(f1)

    f2 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/PNP_BMI_OL4removed.xlsx'          
    BMIgroupedOLremoved.to_excel(f2)
    
    print 'done....'
    return BMIgrouped, BMIgroupedOLremoved

#---------------------------------------------------------------------------
def extract_singleQinfo_DB(QuestionID, QuestionName):
    print 'extracting %s info from DB' % QuestionName
    info = read_sql('select  UserID, QuestionID, Answer,TimeStamp from pnp.questionnaire WHERE QuestionID=%s' % QuestionID, getengine())
    print 'number of answers is %s' % len(info)
    
    #extract year of answer:
    info['Answer_Year'] = info['TimeStamp'].str.split('-').str[0]
    info['Answer_Year']= info['Answer_Year'].str.split('.').str[-1]
    info['Answer_Year']=info['Answer_Year'].astype(str)
    info['Answer_Year']=info['Answer_Year'].apply(lambda x: '20'+str(x) if len(str(x))==2 else str(x))
    info['Answer_Year']=info['Answer_Year'].apply(lambda x: (x.split('/')[2]).split(' ')[0] if '/' in x else str(x))
    print info['Answer_Year'].value_counts(dropna=False)
    
    info = info.drop_duplicates()
#     info2 = info[['UserID', 'Answer']].set_index('UserID')
    
    # add BD:
    # load BD SAMPLE FILE:
    print 'addd BD info:...'
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/BDfile_31052018.xlsx'          
    BDfile = pd.read_excel(f1)
    BDfile_small = BDfile[['DnaID', 'UserID','Blood_Year']]
    BDfile_small = BDfile_small.drop_duplicates()
    BDfile_small = BDfile_small.rename(columns={'DnaID':'BD'})
    # info3=pd.merge(BDfile_small,info2,how='left',left_on='UserID', right_index=True)
    info3 = pd.merge(info, BDfile_small, how='left', left_on=['UserID','Answer_Year'], right_on=['UserID','Blood_Year'])
    info3 = info3[info3.Answer.notnull()]
    info3 = info3[info3.BD.notnull()]
    print 'number of UserIDs with BD samples and %s info is %s' % (QuestionName, len(info3))
    uniqueBD = info3.BD.nunique()
    print 'number of Unique BDs with info is %s' % uniqueBD
    
    info3.Answer = info3.Answer.astype(float)
    if len(info3) > uniqueBD:
        print 'there are multiple BDs per User or multiple values per BD'
        info4 = info3.groupby('BD').agg({'Answer':lambda x: noOutlierMean(x, nSTD=3, nMinSamples=2)})
        print 'number of BD samples after grouping is %s' % len(info4)
    else:
        info4 = info3.set_index('BD')
        
    info5 = pd.DataFrame(info4['Answer'])
    info5 = info5.rename(columns={'Answer':QuestionName})
    
    info5[QuestionName].hist()
    plt.show()
    
    print 'remove outliers...'
    info5removed = filter_outliers(df=info5, outlierSTD=4, trim=False)
    
    print 'info5 describe:'
    print info5.describe()
    
  
    print 'info5removed describe:'
    print info5removed.describe()
    
    print 'saving files...'
    # save with and without outlier removal
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/PNP_%s.xlsx' % QuestionName          
    info5.to_excel(f1)

    f2 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/PNP_%s_OL4removed.xlsx' % QuestionName       
    info5removed.to_excel(f2)
    
    print 'done!'
    
    return info5, info5removed

    
#----------------------------------

def get_glucoseInfo_fromNastya():
    # load file from Nastya:
    file1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/PhenotypicData/X_pnp3_gdm_with_additional_connectons_PCA.dat'   
    NastyaTable = pd.read_pickle(file1)
    
    # GET Interesting continuous phenotypes from table:
    print 'getting interesting continuous phenotypes from table...'
    NastyaPhenotypeShortCont = NastyaTable.reset_index()[['ConnectionID', 'BMI', 'BodyFat', 'HbA1C%', 'Hips',
                                     'Max_Glucose', 'Median_Glucose', 'Min_Glucose', 'Waist', 'WakeupGlucose']]
    NastyaSmokingInfo = NastyaTable.reset_index()[['ConnectionID', 'Ever smoked', 'Currently smokes']]
    
    # merge on connectionID using noOutlierMean:
    NastyaPhenotypeShortCont = NastyaPhenotypeShortCont.astype('float64')
    NastyaPhenotypeShortContGrouped = NastyaPhenotypeShortCont.groupby('ConnectionID').agg(lambda x: noOutlierMean(x, nSTD=5, nMinSamples=3))
    NastyaSmokingInfoGrouped = NastyaSmokingInfo.groupby('ConnectionID').agg(lambda x: noOutlierMean(x, nSTD=5, nMinSamples=3))
    
    # merge on connectionID using most frequent
    from collections import Counter
    
    print 'grouping glucsoe info by connection ID....'
    NastyaSmokingInfoGrouped = NastyaSmokingInfo.groupby('ConnectionID').agg(lambda x: Counter(x).most_common()[0][0])                                                        
                                                                           
    print 'merging with BD info...'                                                                       
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/updatedBDandFDlists/final_BD_FD_converter_mergedOnDate_RCfolderAllSeqProjects_31052018.xlsx'          
    BDFDfile = pd.read_excel(f1)
    BDFDfile_small = BDFDfile[['BD', 'UserID', 'connectionID']]
    BDFDfile_small = BDFDfile_small.drop_duplicates()
    
    print 'grouping info by BD...'
    NastyaPhensBD = pd.merge(NastyaPhenotypeShortContGrouped, BDFDfile_small, how='inner', left_index=True, right_on='connectionID')
    NastyaPhensBDgrouped = NastyaPhensBD.groupby('BD').agg(lambda x: noOutlierMean(x, nSTD=5, nMinSamples=3))
    NastyaPhensBDgrouped = NastyaPhensBDgrouped.drop('connectionID', axis=1)
    
                                                                           
    NastyaSmokingInfoBD = pd.merge(NastyaSmokingInfoGrouped, BDFDfile_small, how='inner', left_index=True, right_on='connectionID')
    NastyaSmokingInfoBDgrouped = NastyaSmokingInfoBD.groupby('BD').agg(lambda x: Counter(x).most_common()[0][0])
    NastyaSmokingInfoBDgrouped = NastyaSmokingInfoBDgrouped.drop('connectionID', axis=1)
                                                                           
    print 'remove outliers...'
    NastyaPhensBDgroupedOLremoved = filter_outliers(df=NastyaPhensBDgrouped, outlierSTD=4, trim=False)

    print 'NastyaPhensBDgrouped describe:'
    print NastyaPhensBDgrouped.describe()
    
    print 'NastyaPhensBDgroupedOLremoved describe:'
    print NastyaPhensBDgroupedOLremoved.describe()
                                                                           
    print 'NastyaSmokingInfoBDgrouped describe:'
    print NastyaSmokingInfoBDgrouped.describe()

    print 'saving files...'
    # save with and without outlier removal
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/PNP_NastyaPhen.xlsx'          
    NastyaPhensBDgrouped.to_excel(f1)

    f2 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/PNP_NastyaPhens_OL4removed.xlsx'          
    NastyaPhensBDgroupedOLremoved.to_excel(f2)
                                                                           
    f3 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/PNP_NastyaSmokingInfo.xlsx'          
    NastyaSmokingInfoBDgrouped.to_excel(f3)
    
    print 'done!'
    
    return  NastyaPhensBDgrouped, NastyaPhensBDgroupedOLremoved, NastyaSmokingInfoBDgrouped
    
#--------------------------------------------------------------------------------------
'''
the following functions take 1 or 2 datasets and plots their phenotypes. 
ifnu there are two datasets, they are compared and statistics are presented (KS and ttest, means for numerical phenotypes, 
chi square p for categorial phenotypes.

input:
numericalphenotypes= list of strings - names of all continuous phenotypes to be compared
ategoricalphenotypes= list of strings - names of all categorical phenotypes to be compared
nBins - int, number of bins in histograms
phenotypeDF1,phenotypeDF2 = dataframe, containing all samples and all phenotypes to be compared (the second cab be None)
sampleList1sampleList2= list of samples to be taken out of the corresponding phenotypeDF. if None, all samples will be taken!
datasetName1,datasetName2 = string, name of datasets for titles and labels


Uasge example:

numericalphenotypes=['Age','BMI','Creatinine']
categoricalphenotypes=['Gender', 'Smoking']
phenotypeDF1=PNP530ss9000_AgeGenderBMIcreatSmoking
sampleList1=None
datasetName1='PNP530ss9000'
phenotypeDF2=Cardio126ss9000_AgeGenderBMIcreatSmoking
sampleList2=None
datasetName2='Cardio126ss9000'

nBins=20

fig1=compare_phenotypes(numericalphenotypes,categoricalphenotypes,nBins,phenotypeDF1,sampleList1,datasetName1,
                       phenotypeDF2,sampleList2,datasetName2)


'''
#### phenotypeDF1 and phenotypeDF2 index should be sample names
# ## second dataset is not obligatory. 


def roundup2(a, digits=0):
    n = 10 ** -digits
    return round(math.ceil(a / n) * n, digits)


def compare_phenotypes(numericalphenotypes, categoricalphenotypes, nBins, phenotypeDF1, sampleList1, sampleListName1,datasetName1,
                       phenotypeDF2=None, sampleList2=None, sampleListName2=None,datasetName2=None,widthUnit=4,heightUnit=4):

    nNum = len(numericalphenotypes)
    nCat = len(categoricalphenotypes)
    
    if sampleList1 is None:
        sampleList1 = phenotypeDF1.index.tolist()
    if (sampleList2 is None) & (phenotypeDF2 is not None):
        sampleList2 = phenotypeDF2.index.tolist()

    fig1 = plt.figure(figsize=(nCat * widthUnit, (nNum + 1) * heightUnit))
    if sampleList2 is not None:  
        fig1.suptitle('Main Phenotype Distribution comparison\n%s vs %s ' % (datasetName1, datasetName2))
    else:
        fig1.suptitle('Main Phenotype Distribution - %s' % datasetName1)

    # plotting numerical phenotypes:
       
    for n, phenotype in enumerate(numericalphenotypes):
        print n, phenotype
        data1 = phenotypeDF1.loc[sampleList1, phenotype]
        data1 = data1[data1.notnull()]
        data1 = list(data1)
        weights1 = np.ones_like(data1) / len(data1)
        mean1 = round(np.mean(data1), 2)
        std1 = round(np.std(data1), 2)
    #     print len(data1)
            
            
        if sampleList2 is not None:
            
            if sampleListName1 is not None:
                label1='_'.join([datasetName1,sampleListName1])
                label2='_'.join([datasetName2,sampleListName2])
            else:
                label1=datasetName1
                label2=datasetName2
            
            filename = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/realAnalysis/MainPhenotypecomparison_%s_%s' % (label1, label2)
        
            data2 = phenotypeDF2.loc[sampleList2, phenotype]
            data2 = data2[data2.notnull()]
            data2 = list(data2)
            weights2 = np.ones_like(data2) / len(data2)
            mean2 = round(np.mean(data2), 2)
            std2 = round(np.std(data2), 2)
            Alldata = data1 + data2
            Allweights = np.ones_like(Alldata) / len(Alldata)

            ax = fig1.add_subplot(nNum + 1, 1, n + 1)
            plot = ax.hist((data1, data2), bins=nBins, color=('black', 'green'), weights=[weights1, weights2],
                         label=(datasetName1, datasetName2), alpha=0.7)


            ks_s_cohort1_cohort2, ks_p_cohort1_cohort2 = stats.ks_2samp(data1, data2)
            t_s_cohort1_cohort2, t_p_cohort1_cohort2 = stats.ttest_ind(data1, data2)
            

            ax.annotate('KS_p=%s\nttest_p=%s\n%s mean=%s\n%s mean=%s' % (round(ks_p_cohort1_cohort2, 6), round(t_p_cohort1_cohort2, 6),
                                label1, round(mean1, 3), label2, round(mean2, 3)),
                                xy=(0.96, 0.95), xycoords='axes fraction', fontsize=9, horizontalalignment='right', verticalalignment='top', fontweight='bold') 
            ax.set_title(phenotype, fontsize=16, fontweight='bold')
        #     ax.set_ylabel('Frequency',fontsize=9)
            if n == 0:
                ax.legend(bbox_to_anchor=(1.01, 0.95), loc='upper left', borderaxespad=0., fontsize=16)
            else:
                ax.legend().set_visible(False)
                
        else:
            if sampleListName1 is not None:
                label1='_'.join([datasetName1,sampleListName1])
            else:
                label1=datasetName1
            filename = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/realAnalysis/MainPhenotypeDist_%s' % label1
                
            ax = fig1.add_subplot(nNum + 1, 1, n + 1)
            plot = ax.hist(data1, bins=nBins, color='black', weights=weights1,
                         label=label1, alpha=0.7)
            
            ax.annotate('mean=%s' % round(mean1, 3), xy=(0.96, 0.95), xycoords='axes fraction', fontsize=11,
                        horizontalalignment='right', verticalalignment='top', fontweight='bold')

            ax.set_title(phenotype, fontsize=16, fontweight='bold')
           
    
    
    # plotting binary phenotypes:
    
    for n, Cphenotype in enumerate(categoricalphenotypes):
        print n, Cphenotype

        ax = fig1.add_subplot(nNum + 1, nCat, nNum * nCat + n + 1)   
        a1norm = phenotypeDF1.loc[sampleList1, Cphenotype].value_counts(normalize=True)
        a1 = phenotypeDF1.loc[sampleList1, Cphenotype].value_counts()
        
        if sampleList2 is not None:
            if sampleListName1 is not None:
                label1='_'.join([datasetName1,sampleListName1])
                label2='_'.join([datasetName2,sampleListName2])
            else:
                label1=datasetName1
                label2=datasetName2
                
            from scipy.stats import chi2_contingency
            a2norm = phenotypeDF2.loc[sampleList2, Cphenotype].value_counts(normalize=True)
            a2 = phenotypeDF2.loc[sampleList2, Cphenotype].value_counts()
            
            mergedNorm = pd.merge(pd.DataFrame(a1norm), pd.DataFrame(a2norm), how='outer', left_index=True, right_index=True)
            merged = pd.merge(pd.DataFrame(a1), pd.DataFrame(a2), how='outer', left_index=True, right_index=True)
            merged=merged.fillna(0)
            try:
                chi, chi_p, dof, expected = chi2_contingency(merged)
            except: 
                print 'couldnt calculate chi test for this phenotype'
                chi=np.nan; chi_p=np.nan; dof=np.nan; expected=np.nan;
            mergedNorm.plot(kind='bar', ax=ax, label=(label1, label2), color=('black', 'green'))
            ax.legend('')
            
            ax.annotate('ChiSquare_p=%s' % round(chi_p, 6), xy=(0.96, 0.95), xycoords='axes fraction', fontsize=9, horizontalalignment='right', verticalalignment='top', fontweight='bold')
        else:
            a1.plot(kind='bar', ax=ax, color='black')
        
        ax.set_title(Cphenotype, fontsize=16, fontweight='bold')      
        
    fig1.text(0, 0.5, "Frequency", ha='left', fontsize=18, rotation=90)
    fig1.subplots_adjust(left=0.13, right=0.98, top=0.8, bottom=0.02, wspace=0.22, hspace=0.25)
   
    fig1.savefig(filename, bbox_inches='tight', dpi=200) 
    
    plt.show()
    
    
    print 'figure was saved in folder TCR_real_data/realAnalysis'
    print 'done'
    return fig1


#----------------------------------------------------------------------------------------------------------------


def extract_blood_and_urine_data_fromDB():
    import time
    cdate = str(time.strftime("%d%m%Y"))
    cdate

    print 'extracting questionaire from DB...'
    quest = read_sql('select  UserID, QuestionID, Answer, Timestamp from pnp.questionnaire', getengine())
    print 'question list length is %s' % len(quest)
    # # extract question information:
    lookup = read_sql('select  QuestionID,Question, QuestionType from pnp.questions', getengine())
    print 'lookup list length is %s' % len(lookup)
    # merge:
    questWithLookup = pd.merge(quest, lookup, how='left', left_on='QuestionID', right_on='QuestionID')
    print len(questWithLookup)
    # define interesting questions:
    # QuestionType=2: blood test measures
    # questionType=7 Urine measures
    # 1044-Hba1c
    # 1045-glucose(?)
    BloodTestQuestions = lookup[ lookup['QuestionType'] == 2]
    BloodTestQuestionsList = list(BloodTestQuestions['QuestionID'])
    UrineTestQuestions = lookup[ lookup['QuestionType'] == 7]
    UrineTestQuestionsList = list(UrineTestQuestions['QuestionID'])
    Hba1c = [1044]
    glucose = [1045]
    
    interestingQuestions = BloodTestQuestionsList + UrineTestQuestionsList + Hba1c + glucose
    interestingQuestionsDF = pd.DataFrame(interestingQuestions)
    
    questWithLookup_Interesting = pd.merge(interestingQuestionsDF, questWithLookup, how='left', left_on=0, right_on='QuestionID')
    print 'list of interesting question length is %s' % len(questWithLookup_Interesting)
    
    # extract question year:
    for n in questWithLookup_Interesting.index:
        TimeStamp = str(questWithLookup_Interesting.loc[n, 'Timestamp'])
#         print n,TimeStamp
        if TimeStamp is not None:
            p = re.compile('(?<=20)\w+')
            m = p.search(TimeStamp)
            if m is not None:
#                 print str(m.group(0))
                questWithLookup_Interesting.loc[n, 'year'] = '20' + str(m.group(0))
        
    questWithLookup_Interesting2 = questWithLookup_Interesting.drop(0, axis=1)
    questWithLookup_Interesting2['UserID'] = questWithLookup_Interesting2['UserID'].astype(str)
    print questWithLookup_Interesting2.head()
    questWithLookup_Interesting2['UserID'] = questWithLookup_Interesting2['UserID'].str.split('.').str[0]
    questWithLookup_Interesting2['UserID_qYear'] = questWithLookup_Interesting2['UserID'].astype(str).str.cat\
    (questWithLookup_Interesting2['year'], sep='_')
    
    # save file:
    print 'saving file...'
    file4 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/NewPhenotypicData/blood_and_urine_phenotypes_%s.csv' % cdate
    questWithLookup_Interesting2.to_csv(file4)
    
    print 'question divided by year:'
    print questWithLookup_Interesting2['year'].value_counts()


    return questWithLookup_Interesting2
    
    


    
