'''
this file contains functions involved in calculation of suffcient statistics
for the VDJ recombination modeling
'''
#-------------------------------------------------------------------------
# # necessary imports:
from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
from scipy import stats
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot, draw_correlation_scatter
from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
from Bio.SeqUtils import GC
import seaborn as sns
import random
from pop_organize import get_sample_data, get_sample_with_dfs


#-----------------------------------------------------------------------


# # generate a df of ins1 sequences and length for a specific sample_name:


def genIns1SeqDF(sample_name, basicUnit):
    
    columns_to_check = ['cdr3Length', 'vGeneName', 'dFamilyName', 'jGeneName', 'vDeletion', 'n1Insertion', 'd5Deletion', 'd3Deletion',
                 'n2Insertion', 'jDeletion', 'vIndex', 'n1Index', 'dIndex', 'n2Index', 'jIndex']
    
    sample_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/%s.tsv" % sample_name) 
    sample_df_non_prod = sample_df[sample_df['sequenceStatus'] != 'In']
    sample_df_non_prod_combinations = sample_df_non_prod.drop_duplicates(subset=columns_to_check)  # work in the combination level
    print ('n combinations in sample %s is %s' % (sample_name, len(sample_df_non_prod_combinations)))
    if basicUnit == 'seq':
        df = sample_df_non_prod
    else:
        df = sample_df_non_prod_combinations
    
    
    seqList = list(df['nucleotide'])
    n1IndexList = list(df['n1Index'])
    dIndexList = list(df['dIndex'])
    n1InsSeqList = []
    n1lengthList = []
    for n, seq in enumerate(seqList):
        n1Index = n1IndexList[n]
        dIndex = dIndexList[n]
        if n1Index == -1 or dIndex == -1:
            n1InsSeqList.append(None)
            n1lengthList.append(0)
        else:
            n1InsSeq = seq[n1Index:dIndex]
            n1InsSeqList.append(n1InsSeq)
            n1lengthList.append(len(n1InsSeq))
    ins1Seq_df = pd.DataFrame({'sequence': n1InsSeqList, 'length':n1lengthList})
    print ('n ins1 seqs in sample %s is %s' % (sample_name, len(ins1Seq_df)))
    return ins1Seq_df

def genIns2SeqDF(sample_name, basicUnit):
    from Bio.Seq import Seq
    columns_to_check = ['cdr3Length', 'vGeneName', 'dFamilyName', 'jGeneName', 'vDeletion', 'n1Insertion', 'd5Deletion', 'd3Deletion',
                 'n2Insertion', 'jDeletion', 'vIndex', 'n1Index', 'dIndex', 'n2Index', 'jIndex']
    
    sample_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/%s.tsv" % sample_name) 
    sample_df_non_prod = sample_df[sample_df['sequenceStatus'] != 'In']
    sample_df_non_prod_combinations = sample_df_non_prod.drop_duplicates(subset=columns_to_check)  # work in the combination level
    print ('n combinations in sample %s is %s' % (sample_name, len(sample_df_non_prod_combinations)))
    
    if basicUnit == 'seq':
        df = sample_df_non_prod
    else:
        df = sample_df_non_prod_combinations
    
    seqList = list(df['nucleotide'])
    n2IndexList = list(df['n2Index'])
    jIndexList = list(df['jIndex'])
    n2InsSeqList = []
    n2lengthList = []
    # seqListShort=seqList[:5]
    for n, seq in enumerate(seqList):
        n2Index = n2IndexList[n]
        jIndex = jIndexList[n]
        if n2Index == -1 or jIndex == -1:
            n2InsSeqList.append(None)
            n2lengthList.append(0)
        else:
            n2InsSeqRev = Seq(seq[n2Index:jIndex])
            n2InsSeq = n2InsSeqRev.reverse_complement()
            n2InsSeq = ''.join(n2InsSeq)
            n2InsSeqList.append(n2InsSeq)
            n2lengthList.append(len(n2InsSeq))
        
    n2seq_df = pd.DataFrame({'sequence': n2InsSeqList, 'length':n2lengthList})
    print ('n ins2 seqs in sample %s is %s' % (sample_name, len(n2seq_df)))
    return n2seq_df





# # generate ins1 sequence and length df from ALL sample:
def genInsSeqDF_allSamples(toSave, insType, basicUnit):
    df_file_names, samples_with_df = get_sample_with_dfs()
    
    insSeq_sample_df_list = []
    count = 1
    for sample_name in samples_with_df:
        print count
        if insType == 1:
            insSeq_sample_df = genIns1SeqDF(sample_name, basicUnit)
        if insType == 2:
            insSeq_sample_df = genIns2SeqDF(sample_name, basicUnit)
        insSeq_sample_df['Sample'] = sample_name
        insSeq_sample_df_list.append(insSeq_sample_df)
        # print len(insSeq_sample_df)
        count += 1
    
    AllinsSeqsDF = pd.concat(insSeq_sample_df_list)
    print len(AllinsSeqsDF)
    AllinsSeqsDFsmall = AllinsSeqsDF.set_index('Sample')
    AllinsSeqsDFsmall.index = AllinsSeqsDFsmall.index.str.replace("HIP", "")
    
    
    if toSave:
        with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins%s_dfs/AllinsSeqsDFsmall_ins%s_%s' % (insType,
                                                                                                                                                    insType, basicUnit), "wb") as f:
            pickle.dump(AllinsSeqsDFsmall, f)
        f.close()
        writer = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins%s_dfs/AllinsSeqsDFsmall_ins%s_%s.xlsx' % (insType, insType, basicUnit)
        AllinsSeqsDFsmall.to_excel(writer, sheet_name='Sheet1', na_rep='', float_format=None, columns=None, header=True, index=True, index_label=None, startrow=0, startcol=0, engine=None, merge_cells=True, encoding=None, inf_rep='inf', verbose=True)  # # saving the correct Reg table to pickles and excel:
    return AllinsSeqsDFsmall


def divide_data_into_TrainAndTest(dataDF, TestFraction, toShuffle=True):
    dataLength = len(dataDF)
    trainLength = int(round(TestFraction * dataLength, 0))
    
    if toShuffle:
        shuffledDF = dataDF.apply(np.random.permutation)  
        trainSet = shuffledDF[:trainLength]
        TestSet = shuffledDF[trainLength:]
    else:
        trainSet = dataDF[:trainLength]
        TestSet = dataDF[trainLength:]
        
    
    return trainSet, TestSet

def divide_data_into_TrainAndTest_AllSamples(AllinsSeqsDFsmall, sample_name):
    
    shortSampleName = sample_name.replace("HIP", "")
    trainSet = AllinsSeqsDFsmall[AllinsSeqsDFsmall.index != shortSampleName]
    testSet = AllinsSeqsDFsmall[AllinsSeqsDFsmall.index == shortSampleName]
    
    print len(trainSet.index.unique())
    print len(testSet.index.unique())
    
    return trainSet, testSet
    
#-------------------------------------------------------------------------------------------------------------

# # generate dataframe summarizing the conditional probabilities of each 
# # nt given the former nt, based on a sequence list

def gen_dinucNormDF(seqList, addPrior):
    # import pandas as pd
    nt2 = ['A', 'C', 'G', 'T']
    if addPrior == 'weak constant':
        A = [1, 1, 1, 1]
        C = [1, 1, 1, 1]
        G = [1, 1, 1, 1]
        T = [1, 1, 1, 1]
    elif addPrior == 'weak differ': 
        A = [3, 3, 2, 2]
        C = [2, 5, 2, 2]
        G = [2, 2, 4, 2]
        T = [2, 3, 2, 3]
    elif addPrior == 'strong constant':
        A = [100, 100, 100, 100]
        C = [100, 100, 100, 100]
        G = [100, 100, 100, 100]
        T = [100, 100, 100, 100]
    elif addPrior == 'strong differ':
        A = [300, 300, 200, 200]
        C = [150, 500, 200, 150]
        G = [200, 200, 400, 200]
        T = [200, 300, 200, 300]
        
        
        
    dinuc_df = pd.DataFrame({'nt2':nt2, 'A': A, 'C': C, 'G': G, 'T': T})
    dinuc_df.set_index('nt2', inplace=True)
    for seq in seqList:
        if seq is not None:
            # print seq
            for n in range(1, len(seq)):
                nt1 = seq[n - 1]
                nt2 = seq[n]
                dinuc = nt1 + nt2
                dinuc_df.loc[nt2, nt1] += 1
        # print dinuc_df
    dinuc_norm_df = dinuc_df / dinuc_df.sum()
    
    return dinuc_norm_df


# # generate a dictionary summarizing the probabilities of each nt as the first 
# # nt in the sequence, based on a sequence list
def gen_nt1FreqsDict(seqList, addPrior):
    if addPrior == 'weak constant':
        nt1freq = {'A':1, 'C':1, 'G':1, 'T':1}
    elif addPrior == 'weak differ':
        nt1freq = {'A':2, 'C':4, 'G':2, 'T':2}
    elif addPrior == 'strong constant':
        nt1freq = {'A':100, 'C':100, 'G':100, 'T':100}
    elif addPrior == 'strong differ':
        nt1freq = {'A':200, 'C':350, 'G':230, 'T':220}
        
    for seq in seqList:
        if seq is not None:
            nt1 = seq[0]
            nt1freq[nt1] += 1
    factor = 1.0 / sum(nt1freq.itervalues())
    nt1NormFreq = {k: v * factor for k, v in nt1freq.iteritems() }

    return nt1NormFreq



# # generate a series summairzing the probabilities of each ins1 length
def calc_length_prob(insSeq_df):
    lengthCount = insSeq_df.length.value_counts(normalize=True)
    return lengthCount

# # generates summaries of probs for each ins1 sequence features:  dinucleotide conditional probabilities, first nt 
# # probabilities and insertion length probability

# # inputs:
# # (1) a df containing all sequences ('sequence') and their corresponding lengths
# # ('length')
# # (2) a length value above which, all length values will be treated together 
# # (decide on this value ## using a hist of all length values in your 
# # sequence list.
# # outputs:
# # (1) a dict of dataframes, each corresponsing to different sequence length
# # (2) a dict of dictionaries, each correspond to different length
# # (3) a series summarizing the probabilities of each insertion length



def genSeqProbsForLengths(insSeq_df, threshLength, addPrior):
    dinucNormDFDict = {}
    nt1FreqsDictDict = {}
    
    for length in range(0, threshLength + 1):
        print length
        df = insSeq_df[insSeq_df['length'] == length]
        # print length
        # print len(df)
        seqListForLength = list(df['sequence'])
        dinucNormDF_forLength = gen_dinucNormDF(seqListForLength, addPrior)
        nt1freqs_forLength = gen_nt1FreqsDict(seqListForLength, addPrior)
        dinucNormDFDict[length] = dinucNormDF_forLength
        nt1FreqsDictDict[length] = nt1freqs_forLength
    print 'done with normal lengths...    '
    
    Longerdf = insSeq_df[insSeq_df['length'] > threshLength]
    seqListLonger = list(Longerdf['sequence'])
    dinucNormDF_Longer = gen_dinucNormDF(seqListLonger, addPrior)
    nt1freqs_Longer = gen_nt1FreqsDict(seqListLonger, addPrior)
    dinucNormDFDict['Longer'] = dinucNormDF_Longer
    nt1FreqsDictDict['Longer'] = nt1freqs_Longer
    lengthCount = calc_length_prob(insSeq_df)
    
    return dinucNormDFDict, nt1FreqsDictDict, lengthCount

# # This function calculate the log probability of an ins sequence, based on the probabilities of the insertion length, 
# # dinucleotide sequence and the first nucleotide:
# # inputs:
# # (1) the ins sequence
# # (2) nt1FreqsDictDict= a dict of dictionaries, each correspond to different length
# # (3) dinucNormDFDict= a  dict of dataframes, each corresponsing to different sequence length
# # (4) lengthCount = a series summarizing the probabilites of each insertion length

def calc_insSeq_prob(seq, threshLength, nt1FreqsDictDict, dinucNormDFDict, lengthCount, capping=None):
    # import numpy as np 
    if capping == 'weak':
        capping_value = -6
    elif capping == 'strong':
        capping_value = -3
        
    
    if seq is None:
        seqProb = lengthCount[0]
        logProbSeq = np.log10(seqProb)
        if capping is not None:
            if logProbSeq < capping_value:
                logProbSeq == capping_value  
    else:
        length = len(seq)
        if length in lengthCount:
            if length > threshLength:
                correct_length = 'Longer'
            else:
                correct_length = length
            
            nt1 = seq[0] 
            nt1prob = nt1FreqsDictDict[correct_length][nt1]
            lognt1prob = np.log10(nt1prob)
            lengthProb = lengthCount[length]
            logLengthProb = np.log10(lengthProb)
            if capping is not None:
                if logLengthProb < capping_value:
                    logLengthProb == capping_value 
            logProbSeq = lognt1prob + logLengthProb
            for n in range(1, len(seq)):
                nt = seq[n]
                nt_1 = seq[n - 1]
                depProbNt = dinucNormDFDict[correct_length].loc[nt, nt_1]
                logDepProbNT = np.log10(depProbNt)
                if capping is not None:
                    if logDepProbNT < capping_value:
                        logDepProbNT == capping_value
                
                logProbSeq = logProbSeq + logDepProbNT
        else:
            logProbSeq = np.nan
        
        
    return logProbSeq

def calc_observedFractLog(TestSeq, TestSet):
    TestSequenceFrac = TestSet.sequence.value_counts(normalize=True)   
    

    if TestSeq in TestSequenceFrac:
        
        TestSeqFreq = TestSequenceFrac[TestSeq]
        TestSeqFreqLog = np.log10(TestSeqFreq)
    else:
        TestSeqFreqLog = np.nan
    
    return TestSeqFreqLog

def calc_observed_value_thresh(TestSet):
    totalTestSeqNum = len(TestSet)
    threshold = np.log10(2. / totalTestSeqNum)
    return threshold



#------------------------------------------------------------------------------------------------

# # the following function summarizes all functions needed to calculate ins information:
# # count frequencies in train and test sets, calculate probs, count observed frequencies and comparing them. 
# # plotting relevant images and saving the data

# # insType = - choose 1 for ins1 and 2 for ins2

def calc_plot_save_modelParamPrediction(sample_name, basicUnit, insType, threshLength, train_fraction, axB, toShuffle=True, addPrior=True, capping=None):
    
    # # step 1: generate train and test data sets, count ins1 sequence frequencies in those sets:
    
    print 'step1: generate data sets and count parameter frequencies...'
    
    strTrainFrac = str(train_fraction).replace(".", "")
    import time
    cdate = str(time.strftime("%d%m%Y"))
    cdate
    
       
    if insType == 1:
        insSeq_df = genIns1SeqDF(sample_name, basicUnit)  # generate sequence and length df
    if insType == 2:
        insSeq_df = genIns2SeqDF(sample_name, basicUnit)
    sample_length = len(insSeq_df)
    print 'sample length is %s' % sample_length
    trainSet, TestSet = divide_data_into_TrainAndTest(insSeq_df, train_fraction, toShuffle)  # devide this df into train and test sets
    print 'trainSet start with %s' % trainSet[:5]
    print 'testSet start with %s' % TestSet[:5]


    dinucNormDFDict_train, nt1FreqsDictDict_train, lengthCount_train = genSeqProbsForLengths(trainSet, threshLength, addPrior) 
    # calculate probs in the train set (length, first nt, dinucelotides)
    
    dinucNormDFDict_test, nt1FreqsDictDict_test, lengthCount_test = genSeqProbsForLengths(TestSet, threshLength, addPrior) 
    # calculate probs in the test set (length, first nt, dinucelotides)
    
    
    # # step 2: for each sequence, calculate the predicted and observed frequencies and summarize in a df:
    
    # #loop over each sequence in the test set:
    # # (1) calculate its predicted dreuquency based on the train set probs.
    # # (2) calculate its predicted dreuquency based on the test set probs.
    # # (3) calculate its observed  freuquency in the test set. 
    # # generate a df summarizing all sequences, their predicted and observed frequencies. save it to pickles
    
    
    print 'step2: compare predicted and observed...'
    TestSeqList = []
    calcFreqTrainList = []
    calcFreqTestList = []
    TestSeqFreqLogList = []

    for TestSeq in list(TestSet['sequence']):
        # print TestSeq
        calcFreqTrain = calc_insSeq_prob(TestSeq, threshLength, nt1FreqsDictDict_train, dinucNormDFDict_train, lengthCount_train, capping)
        calcFreqTest = calc_insSeq_prob(TestSeq, threshLength, nt1FreqsDictDict_test, dinucNormDFDict_test, lengthCount_test, capping)
        TestSeqFreqLog = calc_observedFractLog(TestSeq, TestSet)
        TestSeqList.append(TestSeq)
        calcFreqTrainList.append(calcFreqTrain)
        calcFreqTestList.append(calcFreqTest)
        TestSeqFreqLogList.append(TestSeqFreqLog)

    insSeqProbsCompare = pd.DataFrame({'sequence': TestSeqList, 'calcFreqTrain':calcFreqTrainList, 'calcFreqTest':calcFreqTestList,
                                     'TestSeqFreqLog':TestSeqFreqLogList})
    insSeqProbsCompare = insSeqProbsCompare.drop_duplicates()
    observedOnceValue = insSeqProbsCompare['TestSeqFreqLog'].min()
    insSeqProbsCompare['PredictedToAppear'] = np.where(insSeqProbsCompare['calcFreqTrain'] >= observedOnceValue, 1, 0)
    perc_predicted_to_appear = float(insSeqProbsCompare['PredictedToAppear'].sum()) / len(insSeqProbsCompare)
    perc_predicted_to_appear = round(perc_predicted_to_appear * 100, 2)
    
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins1_dfs/expObsFreqDF_%s_%s_%s' 
              % (sample_name, strTrainFrac, cdate), 'wb') as f:
        pickle.dump(insSeqProbsCompare, f)
    f.close()
    
    
    # # step 3: extract only sequences which are predicted, or observed to repeat twice or more
    # # (note that in any case those sequences are only those who appear at least once in the test set)
    
    
    print 'step 3: plotting....'
    threshold = calc_observed_value_thresh(TestSet)  # #define the value which represtn two repeats in the test set:
    insSeqProbsCompareNoDupThresh = insSeqProbsCompare[(insSeqProbsCompare['TestSeqFreqLog'] >= threshold) | 
                                                       (insSeqProbsCompare['calcFreqTrain'] >= threshold)]
    
    # # define parameters and plot correlation images with p and r values :
    df = insSeqProbsCompareNoDupThresh
    
    # x_var='calcFreqTrain'
    # y_var='calcFreqTest'
    # filename='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/SuffStat Images/%s_%s_%s_%s_%s' %(sample_name, strTrainFrac, x_var,y_var,cdate)
    
    
    # figPredCorrel,axPredCorrel,rPredCorrel,pPredCorrel=draw_correlation_scatter(df[x_var], df[y_var], figsize = (6, 6), 
    #                        xticks=None, yticks=None,\
    #                             xlim = None, ylim = None, r = 'pearson', ms=4, logd = False,\
    #                             xlab = x_var, ylab = y_var, filename = filename, title = '%s_correlation between %s and %s\n (train set fraction=%s)' %(sample_name, x_var, y_var,train_fraction),
    #                            c = "b", alpha=0.2, grid = True, dpi = 300, xticklabels = None, 
    #                            add_identity=True, contour = False)

    # plt.show()
    
    x_var = 'calcFreqTrain'
    y_var = 'TestSeqFreqLog'
    filename = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/SuffStat Images/%s_%s_%s_%s_%s' % (sample_name, strTrainFrac, x_var, y_var, cdate)
    figExpToObsCorrel, axExpToObsCorrel, rExpToObsCorrel, pExpToObsCorrel = draw_correlation_scatter_forSubplots(df[x_var], df[y_var], figsize=(6, 6),
                             xticks=None, yticks=None, \
                                 xlim=None, ylim=None, r='pearson', ms=4, logd=False, \
                                 xlab=None, ylab=None, filename=filename,
                                 title='sample %s\nTrain fraction=%s\nPercent predicted to appear=%s' % (sample_name, train_fraction, perc_predicted_to_appear),
                                 c="b", alpha=0.2, grid=True, dpi=300, xticklabels=None,
                                add_identity=True, contour=False, axB=axB)

    # plt.show()
    #   "#a0a0a0"
    
    
    return figExpToObsCorrel, axExpToObsCorrel, rExpToObsCorrel, pExpToObsCorrel, sample_length

#-------------------------------------------------------------------------------------------------

  # # this function predicts ins probs in a sample, based on all other samples

# # input:



# (sample_name, basicUnit, insType, threshLength, train_fraction, axB,toShuffle=True, addPrior=True)

def predictInsoutSampleBasedOnAllOthers(sample_name, AllnsSeqsDFsmall, basicUnit, insType, axB, addPrior, capping, Trim):
    print 'make sure that ins seq DF matches the selected basicUnit!!!!'
    print sample_name
    
    
    shortSampleName = sample_name.replace("HIP", "")
    
    # generating train set and test set for each sample (the test set is all sequences in a sample and the train set is 
    # sequences from all other samples)
    
    print 'generating train and test sets...'
    trainSet = AllnsSeqsDFsmall[AllnsSeqsDFsmall.index != shortSampleName]
    testSet = AllnsSeqsDFsmall[AllnsSeqsDFsmall.index == shortSampleName]
    
    print len(trainSet.index.unique())
    print len(testSet.index.unique())
    
    threshLength = 20
    # calculating probs for the train set (long!!)
    
    # calculating probs for the test set
    print 'calculating probs for the test set...'
    dinucNormDFDict_test, nt1FreqsDictDict_test, lengthCount_test = genSeqProbsForLengths(testSet, threshLength, addPrior)
       
    print 'calculating probs for the train set (long!!)...'
    dinucNormDFDict_train, nt1FreqsDictDict_train, lengthCount_train = genSeqProbsForLengths(trainSet, threshLength, addPrior)
    # comparing predicted and observed frequencies for all ins1 sequence in the sample:

    print 'comparing predicted and observed frequencies for all ins1 sequence in the sample (long!)...'
   
    TestSeqList = []
    calcFreqTrainList = []
    calcFreqTestList = []
    TestSeqFreqLogList = []

    for TestSeq in list(testSet['sequence']):
        # print TestSeq
        calcFreqTrain = calc_insSeq_prob(TestSeq, threshLength, nt1FreqsDictDict_train, dinucNormDFDict_train, lengthCount_train, capping)
        calcFreqTest = calc_insSeq_prob(TestSeq, threshLength, nt1FreqsDictDict_test, dinucNormDFDict_test, lengthCount_test, capping)
        TestSeqFreqLog = calc_observedFractLog(TestSeq, testSet)
        TestSeqList.append(TestSeq)
        calcFreqTrainList.append(calcFreqTrain)
        calcFreqTestList.append(calcFreqTest)
        TestSeqFreqLogList.append(TestSeqFreqLog)

    insSeqProbsCompare = pd.DataFrame({'sequence': TestSeqList, 'calcFreqTrain':calcFreqTrainList, 'calcFreqTest':calcFreqTestList,
                                     'TestSeqFreqLog':TestSeqFreqLogList})
    insSeqProbsCompare = insSeqProbsCompare.drop_duplicates()
    observedOnceValue = insSeqProbsCompare['TestSeqFreqLog'].min()
    insSeqProbsCompare['PredictedToAppear'] = np.where(insSeqProbsCompare['calcFreqTrain'] >= observedOnceValue, 1, 0)
    perc_predicted_to_appear = float(insSeqProbsCompare['PredictedToAppear'].sum()) / len(insSeqProbsCompare)
    perc_predicted_to_appear = round(perc_predicted_to_appear * 100, 2)
    
    # save the expected to observed freqs dataframe to pickle:
    
    print 'savings dataframes and dictionaries...'
    f1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins_opt/expObsFreqDFopt_%s_%s_%s_%s_%s' % (sample_name, insType, basicUnit, addPrior, capping)
    insSeqProbsCompare.to_pickle(f1)
    
    # save the calculated probs to pickles
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins_opt/dinucNormDFDict_train_%s_%s_%s_%s_%s' % (sample_name, insType, basicUnit, addPrior, capping),'wb') as f2:
        pickle.dump(dinucNormDFDict_train,f2)
    f2.close()

    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins_opt/nt1FreqsDictDict_train_%s_%s_%s_%s_%s' % (sample_name, insType, basicUnit, addPrior, capping),'wb') as f3:
        pickle.dump(nt1FreqsDictDict_train,f3)
    f3.close()

    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins_opt/lengthCount_train_%s_%s_%s_%s_%s' % (sample_name, insType, basicUnit, addPrior, capping),'wb') as f4:
        pickle.dump(lengthCount_train,f4)
    f4.close()
    
    
    print 'calculatine exp to obs correlations:...'
    threshold = calc_observed_value_thresh(testSet)  # #define the value which represtn two repeats in the test set:
    insSeqProbsCompareNoDupThresh = insSeqProbsCompare[(insSeqProbsCompare['TestSeqFreqLog'] >= observedOnceValue) | 
                                                       (insSeqProbsCompare['calcFreqTrain'] >= observedOnceValue)]
    # # define parameters and plot correlation images with p and r values :
    if Trim=='without':
        df = insSeqProbsCompareNoDupThresh
    elif Trim=='trimming':
        df=insSeqProbsCompareNoDupThresh
        TrimmingValue=np.log10(1./len(trainSet))
        df['calcFreqTrain']=np.where(df['calcFreqTrain']<TrimmingValue,TrimmingValue,df['calcFreqTrain'])   
    else:
        df = insSeqProbsCompareNoDupThresh
        TrimmingValue=np.log10(1./len(trainSet))
        df['calcFreqTrain']=np.where(df['calcFreqTrain']<TrimmingValue,np.nan,df['calcFreqTrain'])
    
    
    

    from scipy.stats import pearsonr, spearmanr
    x_var = 'calcFreqTrain'
    y_var = 'TestSeqFreqLog'
    x = df[x_var]
    y = df[y_var]
    nx = np.isnan(x)
    ny = np.isnan(y)
    n = nx + ny
    newx = list(x[~n])
    newy = list(y[~n])
    r, p = pearsonr(newx, newy)
    
        
    return r, p, perc_predicted_to_appear 


#---------------------------------------------------------------------------------------

# # this function predicts ins probs in a sample, based on all other samples

# # input:





def predictInsoutSampleBasedOnAllOthers_whenExist(sample_name, AllnsSeqsDFsmall, basicUnit, insType, axB, addPrior, capping,Trim):
    print 'for sample %s, only df loading is needed' %sample_name
        
    print 'generating test sets...'
    shortSampleName = sample_name.replace("HIP", "")
    testSet = AllnsSeqsDFsmall[AllnsSeqsDFsmall.index == shortSampleName]   
    trainSet = AllnsSeqsDFsmall[AllnsSeqsDFsmall.index != shortSampleName]
    
    print 'loading expObsFreqDFopt table...'
    file1='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SuffStat/exp_and_obs_freqs_ins_opt/expObsFreqDFopt_%s_%s_%s_%s_%s' % (sample_name, insType, basicUnit, addPrior, capping)
    
    insSeqProbsCompare=pd.read_pickle(file1)
    
    print 'calculatine exp to obs correlations:...'
    perc_predicted_to_appear = float(insSeqProbsCompare['PredictedToAppear'].sum()) / len(insSeqProbsCompare)
    perc_predicted_to_appear = round(perc_predicted_to_appear * 100, 2)
    observedOnceValue = insSeqProbsCompare['TestSeqFreqLog'].min()
    threshold = calc_observed_value_thresh(testSet)  # #define the value which represtn two repeats in the test set:
    insSeqProbsCompareNoDupThresh = insSeqProbsCompare[(insSeqProbsCompare['TestSeqFreqLog'] >= observedOnceValue) | 
                                                       (insSeqProbsCompare['calcFreqTrain'] >= observedOnceValue)]
    # # define parameters and plot correlation images with p and r values :
    # df=insSeqProbsCompare
    if Trim=='without':
        df = insSeqProbsCompareNoDupThresh
    elif Trim=='trimming':
        df=insSeqProbsCompareNoDupThresh
        TrimmingValue=np.log10(1./len(trainSet))
        df['calcFreqTrain']=np.where(df['calcFreqTrain']<TrimmingValue,TrimmingValue,df['calcFreqTrain'])   
    else:
        df = insSeqProbsCompareNoDupThresh
        TrimmingValue=np.log10(1./len(trainSet))
        df['calcFreqTrain']=np.where(df['calcFreqTrain']<TrimmingValue,np.nan,df['calcFreqTrain'])
    
    from scipy.stats import pearsonr, spearmanr
    x_var = 'calcFreqTrain'
    y_var = 'TestSeqFreqLog'
    x = df[x_var]
    y = df[y_var]
    nx = np.isnan(x)
    ny = np.isnan(y)
    n = nx + ny
    newx = list(x[~n])
    newy = list(y[~n])
    r, p = pearsonr(newx, newy)
    
        
    return r, p, perc_predicted_to_appear
    
#------------------------------------------------------------------------------------------------------------------------------------------------

# # this function was taken from the group depository and adapted
# # generates scatter plot with r,p indications:


def draw_correlation_scatter_forSubplots(x, y, figsize=(3, 3), xticks=None, yticks=None, \
                             xlim=None, ylim=None, r=None, ms=4, logd=False, \
                             xlab=None, ylab=None, filename=None, title=None,
                             c="#a0a0a0", grid=True, dpi=800, xticklabels=None,
                             add_identity=None, contour=False, axB=None, **figkwargs):
    from scipy.stats import pearsonr, spearmanr
    fig = plt.figure(figsize=figsize, dpi=dpi)
    if axB is None:
       axB = fig.add_subplot(111)     
    if contour:
        print "Contour plot are experimental here"
        # axB.hist2d(x,y,bins = 40,norm=LogNorm())
    else:
        axB.plot(x, y, 'o', c=c, ms=ms, **figkwargs)
    if logd:
        axB.set_xscale('log', basex=2)
        axB.set_yscale('log', basey=2)
    if xticks is not None:
        axB.set_xticks(xticks)
        axB.set_xticklabels(xticks)
    if xticklabels is not None:
        axB.set_xticklabels(xticklabels)
    if yticks is not None:
        axB.set_yticks(yticks)
        axB.set_yticklabels(yticks)
    if xlim is not None:
        axB.set_xlim(xlim)
    if ylim is not None:
        axB.set_ylim(ylim)
    if r is not None: 
        if r == 'pearson':
            nx = np.isnan(x)
            ny = np.isnan(y)
            n = nx + ny
            newx = list(x[~n])
            newy = list(y[~n])
            r, p = pearsonr(newx, newy)
            axB.text(0.01, 0.99, "r=%.4f p=%.4f" % (r, p), transform=axB.transAxes, verticalalignment='top', ha='left', fontsize=14, color='red')
        elif r == 'spearman':
            nx = np.isnan(x)
            ny = np.isnan(y)
            n = nx + ny
            newx = list(x[~n])
            newy = list(y[~n])
            r, p = spearmanr(newx, newy)
            axB.text(0.01, 0.99, "r=%.4f p=%.4f" % (r, p), transform=axB.transAxes, verticalalignment='top', ha='left', fontsize=14, color='red')
        else:
            r = np.nan
            p = np.nan
    if xlab is not None:
        axB.set_xlabel(xlab)
    if ylab is not None:
        axB.set_ylabel(ylab)
    if title is not None:
        axB.set_title(title, fontsize=8)
    if grid:
        axB.grid()
    if add_identity:
        axB.plot(x, x, 'k-', linewidth=0.5)
    if filename is not None:
        fig.savefig(filename, bbox_inches='tight', dpi=dpi)
    axB.margins(0.1, 0.1)
    
    # axB.set_xmargin(0.2); axB.autoscale_view()
    # print r
    # print p
    
    return fig, axB, r, p
#--------------------------------------------------------------------------------

def draw_correlation_scatter_onePlotperFigure(x, y, figsize=(3, 3), xticks=None, yticks=None, \
                             xlim=None, ylim=None, r=None, ms=4, logd=False, \
                             xlab=None, ylab=None, filename=None, title=None,
                             c="#a0a0a0", grid=True, dpi=800, xticklabels=None,
                             xticklabelsSize=None, yticklabelsSize=None, xticksAlign=None,
                             contour=False, **figkwargs):
    from scipy.stats import pearsonr, spearmanr
    fig = plt.figure(figsize=figsize, dpi=dpi)
    axB = fig.add_subplot(111)
    if contour:
        print "Contour plot are experimental here"
        # axB.hist2d(x,y,bins = 40,norm=LogNorm())
    else:
        axB.plot(x, y, 'o', c=c, ms=ms, **figkwargs)
    if logd:
        axB.set_xscale('log', basex=2)
        axB.set_yscale('log', basey=2)
    if xticks is not None:
        axB.set_xticks(xticks)
        axB.set_xticklabels(xticks)
    if xticklabels is not None:
        axB.set_xticklabels(xticklabels)
    if xticklabelsSize is not None:
        axB.tick_params(axis='x', labelsize=xticklabelsSize)
        # axB.set_xticklabels(xticklabels, fontsize=xticklabelsSize)
    if yticklabelsSize is not None:
        axB.tick_params(axis='y', labelsize=yticklabelsSize)
        # axB.set_yticklabels(yticklabels, fontsize=yticklabelsSize)
    if xticksAlign is not None:
        axB.tick_params(axis='x', direction='out', labelrotation=xticksAlign)
    if yticks is not None:
        axB.set_yticks(yticks)
        axB.set_yticklabels(yticks)
    if xlim is not None:
        axB.set_xlim(xlim)
    if ylim is not None:
        axB.set_ylim(ylim)
    if r is not None: 
        if r == 'pearson':
            nx = np.isnan(x)
            ny = np.isnan(y)
            n = nx + ny
            newx = list(x[~n])
            newy = list(y[~n])
            r, p = pearsonr(newx, newy)
            axB.text(0.01, 0.99, "r=%.4f p=%.4f" % (r, p), transform=axB.transAxes, verticalalignment='top', ha='left', fontsize=14, color='red')
        elif r == 'spearman':
            nx = np.isnan(x)
            ny = np.isnan(y)
            n = nx + ny
            newx = list(x[~n])
            newy = list(y[~n])
            r, p = spearmanr(newx, newy)
            axB.text(0.01, 0.99, "r=%.4f p=%.6f" % (r, p), transform=axB.transAxes, verticalalignment='top', ha='left', fontsize=14, color='red')
    if xlab is not None:
        axB.set_xlabel(xlab)
    if ylab is not None:
        axB.set_ylabel(ylab)
    if title is not None:
        plt.title(title, fontsize=16)
    if grid:
        axB.grid()
    if filename is not None:
        fig.savefig(filename, bbox_inches='tight', dpi=dpi)
    axB.margins(0.1, 0.1)
    # axB.set_xmargin(0.2); axB.autoscale_view()
    return fig, axB
