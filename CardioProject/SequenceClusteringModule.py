######################################################################################################
# File: SequenceClusteringModule.py
# Version: 0.0
# Date: 20.12.18
# Updated: 30.12.18
# Shani Ben-Ari Fuchs, shani.ben-ari@weizmann.ac.il
# 
#
# 
# Python version: 2.7
###################################################################################################
# required improvements:
# 1. define general function for concatenations and replace all concat functions
# 2. define general functions for batching and sending to queue and replace unique batching functions
# 3. unify the 'extension' naming along all functions in the module
#
#
##############################################################################################
'''
This module takes AllUniqueWithCounts dataframe (which contains all unique aa sequences for all samples
in a cohort), finds all sequence clusters which includes sequences with a defined levenshtein distance
(default=1) from the cluster head sequence, and then finds how many sequences from each cluster appear 
in each sample in the cohort. 

running options:
1. from the begining. 
2. start with existing levenshtein distance dfs which need to be concatenate to levDistMat (specify
-onlyConcat True while running)
3. start from existing levDistMat (need to specify path while running)
4. start from existing clusterMatrix_dropped (organized clusters after dropping redundant clusters,
need to specify path while running) 
3. 

'''

import argparse
import numpy as np
import pandas as pd
# import modin.pandas as pd
import os
import random
from Levenshtein import distance as levdist
from scipy.spatial.distance import pdist, squareform
from addloglevels import sethandlers
import gc
from os import listdir, mkdir, makedirs
from os.path import isfile, join, isdir, exists
from SegalQueue.qp import qp, fakeqp
import cPickle as pickle
from time import gmtime, strftime
from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import * 
import math
import Utils
from collections import Counter


'''
generate the sequence table to use: productive/non productive/all, specific sample/random samples etc:
'''
def generate_seqTable(command_args):
    
        name_extension = ''
        
        # load sequence table:
        AllUniqueWithCounts = pd.read_pickle(command_args.path_to_AllUniqueWithCounts)
        print ('AllUniqueWithCounts.shape: ', AllUniqueWithCounts.shape)
        
        # filter for specific samples and/or sequence type:
        if command_args.seq_type == 'prod': 
            UniqueSeqs = AllUniqueWithCounts[AllUniqueWithCounts['prod_stat'] == 1]
            print 'only prod sequences are used'
            print ('UniqueSeqs.shape: ', UniqueSeqs.shape)
            name_extension = name_extension + '_prod'
        elif command_args.seq_type == 'nonProd': 
            UniqueSeqs = AllUniqueWithCounts[AllUniqueWithCounts['prod_stat'] == 0]
            print 'only nonProd sequences are used'
            print ('UniqueSeqs.shape: ', UniqueSeqs.shape)
            name_extension = name_extension + '_nonProd'
        else:
            UniqueSeqs = AllUniqueWithCounts
        
        if command_args.sample is not None:
            UniqueSeqs = UniqueSeqs[UniqueSeqs['Sample'] == command_args.sample]
            print 'only sample %s is used' % command_args.sample
            print ('UniqueSeqs.shape: ', UniqueSeqs.shape)
            name_extension = name_extension + '_sample' + command_args.sample
            
        if command_args.n_random_samples > 0:    
            samples = UniqueSeqs['Sample'].unique().tolist()
            randSamples = random.sample(samples, command_args.n_random_samples)
            randDF_list = []
            for sample in randSamples:
                    df = UniqueSeqs[UniqueSeqs['Sample'] == sample]
                    randDF_list.append(df)
            randDF = pd.concat(randDF_list)
            UniqueSeqs = randDF
            print 'only %s random samples are used' % command_args.n_random_samples
            print ('UniqueSeqs.shape: ', UniqueSeqs.shape)
            name_extension = name_extension + '_' + str(command_args.n_random_samples) + 'randSamples'
            
        return UniqueSeqs, name_extension

'''
this function takes subset of sequences,find which sequence pairs have levinshtein distance smaller or equal to max_dist
and save the resulting dfs.
this function is called by gen_levdist_table()
'''

def calc_levDist(seqs, min_seq, max_seq, max_dist, name_extension, levDistDFdir):    
    
    lev_file = 'levdist' + name_extension + '_%s_%s.dat' % (min_seq, max_seq)
    f2 = levDistDFdir + lev_file
    
    if isfile(f2):
        print 'this file already exist, breaking...'
        return
    else:
        # calculate levenstein distance:
        df = pd.DataFrame()
        count = 0
        
        for i in range(min_seq, max_seq):
    #         print ('i=',i)
            seq1 = seqs[i]
            len1 = len(seq1)
            for j in range(i + 1, len(seqs)):
                count = count + 1
                seq2 = seqs[j]                   
                len2 = len(seq2)
                if np.abs(len2 - len1) > max_dist:
                    continue
                else:
            #             print ('len1,len2= ', len(seq1), len(seq2))
                    seq_len_max = max(len1, len2)
                #             print len1,len2,seq_len_max
    
                    dist = 0
    
                    for pos in range(0, seq_len_max):
    
                        if (pos > len1 - 1) or (pos > len2 - 1): dist = dist + 1
                        elif seq1[pos] != seq2[pos]: dist = dist + 1
                #                 print pos,seq1[:pos+1],seq2[:pos+1],dist
                        if dist > max_dist: break
    
                        elif pos == seq_len_max - 1:
                            df.loc[count, 'seq1id'] = i
                            df.loc[count, 'seq2id'] = j
                            df.loc[count, 'dist'] = dist      
        
        print ('sequence analyzed: %s:%s' % (min_seq, max_seq))
        print ('number of pairs analyzed: ' , count)
#         print ('df memory usage is: ', df.memory_usage(index=True, deep=True).sum())     
        print ('number of pairs with distance smaller than %s is %s' % (max_dist, df.shape[0]))       
        
        # save df:
        
        print 'saving lev dist matrix...'
        df.to_pickle(f2)
        print 'lev distmat saved'
       
    return df
    
'''
this function:
1. generate an index file (seqs-seqIDs) and save it
2. divide the sequences into batches and send jobs with the function calc_levDist which in turn finds pairs of sequences with distance smaller 
then max_dist and save them. 
'''

def gen_levdist_table(command_args, seqs, name_extension):
    
    # get parameters:
    max_dist = command_args.levDistThresh
    nSeqsBatch = command_args.nSeqsInBatchLev
    output_dir = command_args.output_dir
    levDistDFdir = output_dir + 'levDistDFdir/'
    if not isdir(levDistDFdir):
        makedirs(levDistDFdir)
    
    
    sethandlers()
    
    if type(seqs) == list:
        print 'seqs type is a list'
    elif type(seqs) == pd.core.frame.DataFrame:
        try:
            seqs = seqs.index.tolist()
            print 'seqs type is a dataframe'
        except: 
            print 'couldnt identify seqs type'
            return
        
    seqs = list(set(seqs))  # remove dupliates from list
    seqs = [x for x in seqs if len(x) > 3]  # remove sequences shorter than 4 aas.
    print ('number of sequences: ', len(seqs))
    
    print 'generating seqIDs table'
    seqIDs = pd.Series(index=range(len(seqs)), data=seqs).rename('seq')
    seqID_file = 'seqID' + name_extension + '.dat'
    f1 = output_dir + seqID_file
    print 'saving seqID_file...'
    seqIDs.to_pickle(f1)
    print 'seqID_file saved'
       
    # divide work into batches:
    tempDir = command_args.tempDir
    if not isdir(tempDir):
        makedirs(tempDir)
    os.chdir(tempDir)
        
    with qp(jobname='levDist', q=['himem7.q'], mem_def="1G", trds_def=2, deleteCSHwithnoerr=True,
                 tryrerun=False, max_u=200) as q:
        q.startpermanentrun()
        wait_on = []
    
        min_seq = 0
        max_seq = min_seq + nSeqsBatch
        print ('expected number of batches is ', math.ceil((float(len(seqs))) / nSeqsBatch))
        
        cnt=0
        while min_seq < len(seqs):
            print 'running levdist calc- seqs %s:%s' % (min_seq, max_seq)
            if max_seq > len(seqs): max_seq = len(seqs)
        
    #         seqs_job=seqs[min_seq:max_seq]  
            wait_on.append(q.method(calc_levDist, kwargs={'seqs':seqs, 'min_seq':min_seq, 'max_seq':max_seq,
                                    'max_dist':max_dist, 'name_extension':name_extension,
                        'levDistDFdir':levDistDFdir}))
            cnt += 1 
            min_seq = min_seq + nSeqsBatch
            max_seq = max_seq + nSeqsBatch
            
            
        print ("Sent %d jobs" % cnt)    
        try:                     
            res = q.waitforresults(wait_on)
        except:
            print 'some error occured'
                
                
        print 'all sequences were analyzed!'
    
    return

'''
this function takes all levDistDFs generated by batches and concatenate them to generate 
the levDistMat
'''

def concat_distmats(command_args):
    print 'concatnating levdist dfs...'
    output_dir = command_args.output_dir
    levDistDFdir = output_dir + 'levDistDFdir/'
    
    files = [levDistDFdir + f for f in listdir(levDistDFdir) if isfile(levDistDFdir + f)]
    dfs = []
    for f in files:
        df = pd.read_pickle(f)
        dfs.append(df)
    levDistMatTotal = pd.concat(dfs, axis=0, ignore_index=True)
    print ('levDistMatTotal memory usage is: ', levDistMatTotal.memory_usage(index=True, deep=True).sum())     
    
    if command_args.saveLevDist:
        try:
            print 'saving levDistMatTotal....'
            f3 = output_dir + 'levDistMatTotal.dat'
            levDistMatTotal.to_pickle(f3)
            print 'levDistMatTotal was saved...'
        except:
            print 'a problem occured while trying to save levDistMatTotal'
        
    if command_args.deleteLevdistDFs:
        for f in files:
            os.remove(f)
        print 'removed all levdist files'
    print 'done concatenating...'
        
    return levDistMatTotal
        
'''
this function takes levDistMat, add self indications (each sequence has a pair with itself with
dist=1, and then group the table by seq1id (cluster head) and seq2id (cluster sequences).
then, it compares all cluster pairs (in which the first cluster has more than 2 sequences) to find clusters which
are contained within other clusters and thus are redundant. this comparison is done in batches,
the lists of clusters to drop are united in the end and the redundant clusters are dropped from the 
cluster matrix which is then saved.

'''    
        
def convert_levDistMat_to_ClusterMatrix(levDistMatTotal, command_args):
   
    tempDir = command_args.tempDir
    os.chdir(tempDir)
    
    print levDistMatTotal.info(memory_usage='deep')
    print 'downcasting floats...'
    levDistMatTotal = levDistMatTotal.apply(pd.to_numeric, downcast='float')
    print levDistMatTotal.info(memory_usage='deep')
    
    # add self cluster to each cluster:
    selfDF = pd.DataFrame({'seq1id': levDistMatTotal.seq1id.unique(), 'seq2id': levDistMatTotal.seq1id.unique(), 'dist':1})
    levDistMatFinal = pd.concat([levDistMatTotal, selfDF])

    # group levDistMat by cluster head sequence:
    print 'grouping table...'
    clusterMatrix_multiSeqs = levDistMatFinal.groupby(['seq1id', 'seq2id'])['dist'].count()
    print clusterMatrix_multiSeqs.head()
    
    getClustersInfo = levDistMatFinal.groupby('seq1id')['seq2id'].count().sort_values(ascending=False)
    
    return clusterMatrix_multiSeqs, getClustersInfo

def drop_clusters(command_args, getClustersInfo, clusterMatrix_multiSeqs):
    
    if command_args.fixed_filtering_perc_per_cohort is not None: 
        extension = '_cohortfiltering%s-%sperc' % (command_args.fixed_filtering_perc_per_cohort,
                                                command_args.maxSharingPerc)
        extension = extension.replace('.', '')
    elif command_args.nShared is not None: extension = '_nShared%s' % command_args.nShared
    
    clusterDroppingDir = command_args.output_dir + 'clustersToDropdir%s/' % extension
    clusterDroppingName = 'clustersToDropdir%s'
    if not isdir(clusterDroppingDir):
        makedirs(clusterDroppingDir)
    merged = pd.merge(clusterMatrix_multiSeqs.reset_index(), pd.DataFrame(getClustersInfo), how='left', left_on='seq1id', right_index=True)
    merged = merged.rename(columns={'seq2id_x':'seq2id', 'seq2id_y':'cluster_length'}).sort_values(by='cluster_length')
    merged_grouped = merged.groupby(['seq1id', 'seq2id'])['cluster_length'].mean()
    cluster_heads = merged_grouped.index.get_level_values('seq1id').unique().tolist()
    cluster_heads_reversed = cluster_heads[::-1]
    
#     clusterList=getClustersInfo.index.astype(int).tolist()
#     clusterListToUse=getClustersInfo[getClustersInfo>2].index.astype(int).tolist()
#     clusterListToUse_indices=getClustersInfo.reset_index().reset_index().set_index('seq1id').loc[clusterListToUse,:]['index'].tolist()
#     
#     clusterListFile=command_args.output_dir+'clusterList.pkl'
#     pickle.dump(clusterList, open(clusterListFile, "wb"))
#     clusterListToUseFile=command_args.output_dir+'clusterListToUse.pkl'
#     pickle.dump(clusterListToUse, open(clusterListToUseFile, "wb"))
#     clusterListToUse_indices_File=command_args.output_dir+'clusterListIndices.pkl'
#     pickle.dump(clusterListToUse_indices, open(clusterListToUse_indices_File, "wb"))
#     print ('number of clusters: ', len(clusterList))
#     print ('number of clusters to check: ',len(clusterListToUse))
#     print clusterList[:5]
#     print clusterListToUse[:5]
    
    # divide work into batches:
    batchSize = command_args.dropping_batch_size 
    sethandlers()
    with qp(jobname='dropClusters', q=['himem7.q'], mem_def="5G", trds_def=2, deleteCSHwithnoerr=True,
                 tryrerun=False, max_u=300) as q:
        q.startpermanentrun()
        wait_on = []
    
        
        min_cluster = 0
        max_cluster = min_cluster + batchSize
        nBatches = math.ceil((float(len(cluster_heads_reversed))) / batchSize)
        print ('expected number of batches is', nBatches)
        
        files = listdir(clusterDroppingDir)
        if len(files) == nBatches:
            print 'all clusterToDrop files exist, skipping clusterDropping calculation step'
        else:
        
            cnt = 0
            while min_cluster < len(cluster_heads_reversed):
                if max_cluster > len(cluster_heads_reversed): max_cluster = len(cluster_heads_reversed)
                
                doNotRun = False
                for f in files:
                    if (str(min_cluster) in f) and (str(max_cluster) in f):
                        doNotRun = True
                        break
                    
                if not doNotRun: 
                    wait_on.append(q.method(filter_identical_clusters2, kwargs={'cluster_heads_reversed':cluster_heads_reversed,
                                'merged_grouped':merged_grouped, 'min_cluster':min_cluster, 'max_cluster':max_cluster,
                                'clusterDroppingDir':clusterDroppingDir}))
                    cnt += 1
                else:
                    print 'clusters %s:%s were already analyzed' % (min_cluster, max_cluster)
                    
                min_cluster = min_cluster + batchSize
                max_cluster = max_cluster + batchSize
            
            print ("Sent %d jobs" % cnt)
            try:                     
                res = q.waitforresults(wait_on)
            except:
                print 'some error occured'      
            print 'all batches were analyzed!'
    
    print 'concatenating and saving all clusters to drop...'
    clustersToDrop = concat_clustersToDrop(command_args, clusterDroppingDir, clusterDroppingName)

    # drop rows and add indication for self sequence being in the cluster:
    # drop redundant clusters
    
    print 'dropping redundant clusters:'
    clusterMatrix_multiSeqs_dropped = clusterMatrix_multiSeqs.drop(clustersToDrop, axis=0)
    
        
    try:
        f4 = command_args.output_dir + 'clusterMatrix_multiSeqs_dropped%s.dat' % extension
        clusterMatrix_multiSeqs_dropped.to_pickle(f4)
        print 'clusterMatrix_multiSeqs_dropped was saved...'
    except:
        print 'clusterMatrix_multiSeqs_dropped was not saved...'
    
    return clusterMatrix_multiSeqs_dropped

'''
this function loads lists of clusters to drop and then unified them to include only
unique values
'''

def concat_clustersToDrop(command_args, clusterDroppingDir, clusterDroppingName):
    print 'concatnating clusters to drop...'
    
    files = [clusterDroppingDir + f for f in listdir(clusterDroppingDir) if isfile(clusterDroppingDir + f)]
    clustersToDrop = []
    for f in files:
        clusters = pd.read_pickle(f)
        clustersToDrop = clustersToDrop + clusters
    
    clustersToDrop = list(set(clustersToDrop))
    print ('final number of clusters to drop: ', len(clustersToDrop))
    
    try:
        print 'saving clustersToDrop....'
        f3 = command_args.output_dir + clusterDroppingName + 'pkl'
        pickle.dump(clustersToDrop, open(f3, "wb"))
        print 'final list of clustersToDrop was saved'
    except:
        print 'a problem occured while trying to save the final list of clustersToDrop'
       
    return clustersToDrop

'''
this function find redundant clusters:
it loops over all clusters in the clusterMatrix:
1. check whether the cluster contains more then 2 samples, otherwise continue.
2. if yes, extract the sequences in the cluster.
3. loop over all other clusters following this cluster, for each of them extract the sequences and check
whether they are contained whitin the first cluster. if so, register the cluster as clusterToDrop.


this function is called in batches by the function convert_levDistMat_to_ClusterMatrix
'''
def filter_identical_clusters2(cluster_heads_reversed, merged_grouped, min_cluster, max_cluster, clusterDroppingDir):
    clustersToDrop = []
    for n, i in enumerate(cluster_heads_reversed[min_cluster:max_cluster]):
        if n % 100 == 0: print n
        length = merged_grouped.loc[i].unique()[0]
#         print i,length
        df = merged_grouped.reset_index(level='seq2id')
        clusters_to_compare = df[(df['seq2id'] == i) & (df['cluster_length'] > length)].index.tolist()
        if len(clusters_to_compare) > 0:
            seqs1 = merged_grouped.loc[i].index.get_level_values('seq2id').tolist()
    #         print ('seqs1: ',seqs1)
            for j in clusters_to_compare:
                seqs2 = merged_grouped.loc[j].index.get_level_values('seq2id').tolist()
    #             print ('seqs2: ',seqs2)
                if set(seqs1).issubset(seqs2):
                    print ('dropping', i, seqs1, j, seqs2)
                    clustersToDrop.append(i)
#                     print 'now breaking, should not see seqs2 again before new seqs1'
                    break
                
    print ('min_caluster=%s,max_cluster=%s,number of clusters to drop=%s' % (min_cluster, max_cluster, len(clustersToDrop)))
    if not isdir(clusterDroppingDir):
        makedirs(clusterDroppingDir)
    f4 = clusterDroppingDir + 'clustersToDrop_%s_%s.dat' % (min_cluster, max_cluster)
    pickle.dump(clustersToDrop, open(f4, "wb"))
    print 'clusterMatrix_final was saved...'
    return clustersToDrop
            
    #  &(~df['seq2id'].isin(clustersToDrop)   
    #     print clusters_to_compare

# def filter_identical_clusters2(clusterMatrix_multiSeqs,clusterListFile,
#                                clusterListToUse_indices_File,min_cluster,max_cluster,output_dir):
# 
#     print 'start working on sequences %s:%s' %(min_cluster,max_cluster)
#     clustersToDrop=[]
#     clusterList=Utils.Load(clusterListFile)
#     clusterListToUse_indices=Utils.Load(clusterListToUse_indices_File)
#     for i in clusterListToUse_indices[min_cluster:max_cluster]:
#         if clusterList[i]  not in clustersToDrop:
#             seqs1=clusterMatrix_multiSeqs.loc[clusterList[i],:].index.get_level_values('seq2id').tolist()
#             for j in range(i+1,len(clusterList)):
#                 if clusterList[j] not in clustersToDrop:
#                     seqs2=clusterMatrix_multiSeqs.loc[clusterList[j],:].index.get_level_values('seq2id').tolist()
# #                     print seqs1,seqs2
#                     if set(seqs2).issubset(seqs1):
#                         print ('dropping', clusterList[i],seqs1,clusterList[j],seqs2)
#                         clustersToDrop.append(clusterList[j])
#     
#     print ('min_caluster=%s,max_cluster=%s,number of clusters to drop=%s' %(min_cluster,max_cluster,len(clustersToDrop)))
#     clusterDroppingDir=output_dir+'clustersToDropdir/'
#     if not isdir(clusterDroppingDir):
#         makedirs(clusterDroppingDir)
#     f4=clusterDroppingDir+'clustersToDrop_%s_%s.dat' %(min_cluster,max_cluster)
#     pickle.dump(clustersToDrop, open(f4, "wb"))
#     print 'clusterMatrix_final was saved...'
#                             
#     return clustersToDrop


def get_info_per_sample(samples, min_sample, max_sample, command_args, clusterMatrix_dropped_seqs,
                            tempDFdir):
        
    AllUniqueWithCounts = pd.read_pickle(command_args.path_to_AllUniqueWithCounts)
    print 'start gathering info from samples %s:%s:' % (min_sample, max_sample)
    sampleByClusterDF_temp = pd.DataFrame()
    
    
    for n, sample in enumerate(samples[min_sample:max_sample]):
#             if n%100==0: print n
        df = AllUniqueWithCounts[AllUniqueWithCounts['Sample'] == sample]
        clusterMatrix = clusterMatrix_dropped_seqs[clusterMatrix_dropped_seqs['seq2'].isin(df.index)]
        clusterMatrix_grouped = pd.DataFrame(clusterMatrix.groupby('seq1')['dist'].count())
#         print clusterMatrix_grouped.sort_values(by='dist',ascending=False).head()
        sampleByClusterDF_temp = pd.merge(sampleByClusterDF_temp, clusterMatrix_grouped, how='outer', left_index=True, right_index=True).rename(columns={'dist':sample})
    sampleByClusterDF_temp.index = sampleByClusterDF_temp.index.rename('cluster')
    sampleByClusterDF_temp = sampleByClusterDF_temp.T
    print 'downcasting...'
    sampleByClusterDF_temp = sampleByClusterDF_temp.apply(pd.to_numeric, downcast='float')
    print ('sampleByClusterDF_temp.shape: ', sampleByClusterDF_temp.shape)
    print ('sampleByClusterDF_temp memory usage: ', sampleByClusterDF_temp.memory_usage(index=True).sum())
    print 'saving dataframe to dat file...'
    f6 = tempDFdir + 'sampleByClusterDF_temp_%s_%s.dat' % (min_sample, max_sample)
    sampleByClusterDF_temp.to_pickle(f6)
    print 'saving dataframe to csv file...'
    f7 = tempDFdir + 'sampleByClusterDF_temp_%s_%s.csv' % (min_sample, max_sample)
    sampleByClusterDF_temp.to_csv(f7)
    print 'saved sampleByClusterDF_temp_%s_%s' % (min_sample, max_sample)

    return 
        
def concat_sampleDFs(command_args):
    print 'concatnating sample DFs...'
    sampleDFdir = command_args.output_dir + 'sampleDFdir/'
    partlyConcatDFdir = sampleDFdir + 'partlyConcatDFdir/'
    if not isdir(partlyConcatDFdir):
        makedirs(partlyConcatDFdir)
    
    sampleDFfiles = [sampleDFdir + f for f in listdir(sampleDFdir) if (isfile(sampleDFdir + f)) and ('.dat' in f)]
    
    # sending jobs to concat 4 dataframes each time
    sethandlers()
    with qp(jobname='concatSamples', q=['himem7.q'], mem_def="10G", trds_def=2, deleteCSHwithnoerr=True,
                 tryrerun=False, max_u=200) as q:
        q.startpermanentrun()
        wait_on = []
        
        batchSize = command_args.samples_concat_batch_size
        min_tempSampleDF = 0
        max_tempSampleDF = min_tempSampleDF + batchSize
        print ('expected number of batches for partial concatenation of sample info is ', math.ceil((float(len(sampleDFfiles))) / batchSize))
        print ('sampleDFdir: ', sampleDFdir)
        
        cnt = 0
        
        dfFiles = listdir(partlyConcatDFdir)
        while min_tempSampleDF < len(sampleDFfiles):
            if max_tempSampleDF > len(sampleDFfiles): max_tempSampleDF = len(sampleDFfiles)
#             print ('files:', files)
            doNotRun = False
            for f in dfFiles:
                if (str(min_tempSampleDF) in f) and (str(max_tempSampleDF) in f):
                    doNotRun = True
                    print f, doNotRun
                    break
            if not doNotRun:
                wait_on.append(q.method(partial_concat_samples, kwargs={'command_args':command_args,
                    'sampleDFfiles':sampleDFfiles, 'min_tempSampleDF':min_tempSampleDF,
                    'max_tempSampleDF':max_tempSampleDF, 'partlyConcatDFdir':partlyConcatDFdir}))
                cnt += 1
            else:
                print 'samples %s:%s were already analyzed' % (min_tempSampleDF, max_tempSampleDF)
        
            min_tempSampleDF = min_tempSampleDF + batchSize
            max_tempSampleDF = max_tempSampleDF + batchSize
        
        print ("Sent %d jobs" % cnt)     
        try:                     
            res = q.waitforresults(wait_on)
            print ('job result: ', res)
        except:
            print 'some error occured'   
       
    print 'all partial concatenations were done'
    
    # final concatenation:...
    print 'final concatenation...'
    sethandlers()
    with qp(jobname='finalConc', q=['himem7.q'], mem_def="10G", trds_def=4, deleteCSHwithnoerr=True,
                 tryrerun=False, max_u=200) as q:
        q.startpermanentrun()
        wait_on = []
        
        wait_on.append(q.method(final_concat, kwargs={'command_args':command_args,
                                                     'partlyConcatDFdir':partlyConcatDFdir}))
        
    try:                     
        res = q.waitforresults(wait_on)
    except:
        print 'some error occured' 
        
    return    
        
        
def final_concat(command_args, partlyConcatDFdir):
    partlyConcatDF_files = listdir(partlyConcatDFdir) 
    dfList = []
    for f in partlyConcatDF_files:
        df = pd.read_pickle(partlyConcatDFdir + f)
        dfList.append(df) 
    final_concat_df = pd.concat(dfList)
    
    f7 = command_args.output_dir + 'sampleByClusterDF.dat' 
    final_concat_df.to_pickle(f7)
    
    return
    
    
def partial_concat_samples(command_args, sampleDFfiles, min_tempSampleDF, max_tempSampleDF, partlyConcatDFdir):    
    
    dfList = []
    for n, f in enumerate(sampleDFfiles[min_tempSampleDF:max_tempSampleDF]):
        print n, f.split('/')[-1]
        sampleDF = pd.read_pickle(f)
#         print 'downcasting...'
#         sampleDF=sampleDF.apply(pd.to_numeric,downcast='float')
        print 'appending...'
        dfList.append(sampleDF)
#         sampleDF.to_csv(command_args.output_dir+'sampleByClusterDF.csv',index=False)
        del sampleDF
    
#     print ('now reading csv file')    
#     sampleByClusterDF=pd.read_csv(command_args.output_dir+'sampleByClusterDF.csv')
    print 'concatenating...'
    sampleByClusterDF_partial = pd.concat(dfList)

    print ('final shape of sampleByClusterDF: ', sampleByClusterDF_partial.shape)
    
    f7 = partlyConcatDFdir + 'sampleByClusterDF_partial_%s_%s.pickle' % (min_tempSampleDF, max_tempSampleDF)
    sampleByClusterDF_partial.to_pickle(f7)
    
    return 

'''
this function takes the final clusterMatrix, replace sequence ids with real sequences based on
the seqID table and then for each sample in the cohort, check which clusters it contains and 
how many sequences per each cluster. the result is a sample-by-cluster table in which each cell 
represent the number of sequences in the cluster in the sample.

the function also filter the table to incldue only clusters that appear in a minimal number of samples
(default=10)
'''

def get_clusters_for_samples(clusterMatrix, seqID, command_args):
    
    if 'sampleByClusterDF.dat' not in listdir(command_args.output_dir):
        
        print 'sampleByClusterDF not exist yet, start extracting cluster info for samples'
    
        print 'replacing seqIDs with real sequences...'
        clusterMatrix_2 = clusterMatrix.reset_index()
        clusterMatrix_seqs = pd.merge(clusterMatrix_2, pd.DataFrame(seqID), how='left', left_on='seq2id', right_index=True).rename(columns={'seq':'seq2'})
        clusterMatrix_seqs = pd.merge(clusterMatrix_seqs, pd.DataFrame(seqID), how='left', left_on='seq1id', right_index=True).rename(columns={'seq':'seq1'})
        clusterMatrix_seqs = clusterMatrix_seqs.drop(['seq1id', 'seq2id'], axis=1)
        
        AllUniqueWithCounts = pd.read_pickle(command_args.path_to_AllUniqueWithCounts)
        samples = AllUniqueWithCounts.Sample.unique()
        
        tempDir = command_args.tempDir
        os.chdir(tempDir)
        
        tempDFdir = command_args.output_dir + 'sampleDFdir/'
        if not isdir(tempDFdir):
            makedirs(tempDFdir)
        
        # get cluster information for each 30 samples:   
        sethandlers()
        with qp(jobname='sampleInfo', q=['himem7.q'], mem_def="4G", trds_def=2, deleteCSHwithnoerr=True,
                     tryrerun=False, max_u=200) as q:
            q.startpermanentrun()
            wait_on = []
            
            batchSize = 30
            min_sample = 0
            max_sample = min_sample + batchSize
            print ('expected number of batches is', math.ceil((float(len(samples))) / batchSize))
            print ('tempDFdir: ', tempDFdir)
            
            cnt = 0
            while min_sample < len(samples):
                if max_sample > len(samples): max_sample = len(samples)
                files = listdir(tempDFdir)
    #             print ('files:', files)
                doNotRun = False
                for f in files:
                    if (str(min_sample) in f) and (str(max_sample) in f):
                        doNotRun = True
                        print f, doNotRun
                        break
                if not doNotRun:
                    wait_on.append(q.method(get_info_per_sample, kwargs={'samples':samples,
                            'min_sample':min_sample, 'max_sample':max_sample, 'command_args':command_args,
                            'clusterMatrix_dropped_seqs':clusterMatrix_seqs, 'tempDFdir':tempDFdir}))
                    cnt += 1
                else:
                    print 'samples %s:%s were already analyzed' % (min_sample, max_sample)
            
                min_sample = min_sample + batchSize
                max_sample = max_sample + batchSize
            
            print ("Sent %d jobs" % cnt)     
            try:                     
                res = q.waitforresults(wait_on)
            except:
                print 'some error occured'      
        print 'all samples were analyzed!'
        
        sampleByClusterDF = concat_sampleDFs(command_args)
        print 'sampleByClusterDF is ready!!!'
    else:
        print 'loading sampleByClusterDF...'
        sampleByClusterDF = pd.read_pickle(command_args.output_dir + 'sampleByClusterDF.dat')
    
    if command_args.fixed_filtering_perc_per_cohort is not None:
        perc = command_args.fixed_filtering_perc_per_cohort
        print 'filtering to includes only clusters that appear in %s perc of healthy, patients or all together:' % perc
        
        # load PNP and Cardio sample lists:
        with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP530', 'rb') as fp:
            PNP530 = pickle.load(fp)
        with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/Cardio126', 'rb') as fp:
            Cardio126 = pickle.load(fp)
        
        
        # get clusters that appear in Xperc of samples in sampleGroup1:
        sampleByClusterDF_group1 = sampleByClusterDF.loc[PNP530, :]
        nShared_group1 = math.ceil(len(sampleByClusterDF_group1) * perc)
        print ('nShared_group1= ', nShared_group1)
        group1_percClusters = sampleByClusterDF_group1.loc[:, sampleByClusterDF_group1.count() >= nShared_group1].columns.tolist()
        print ('n clusters shared by %s or more samples in group1= ' % nShared_group1, len(group1_percClusters))
        
        # get clusters that appear in Xperc of samples in sampleGroup2:
        sampleByClusterDF_group2 = sampleByClusterDF.loc[Cardio126, :]
        nShared_group2 = math.ceil(len(sampleByClusterDF_group2) * perc)
        print ('nShared_group2= ', nShared_group2)
        group2_percClusters = sampleByClusterDF_group2.loc[:, sampleByClusterDF_group2.count() >= nShared_group2].columns.tolist()
        print ('n clusters shared by %s or more samples in group2= ' % nShared_group2, len(group2_percClusters))
        
        # get clusters that appear in Xperc of All samples:
        nShared_all = math.ceil(len(sampleByClusterDF) * perc)
        print ('nShared_all= ', nShared_all)
        all_percClusters = sampleByClusterDF.loc[:, sampleByClusterDF.count() >= nShared_all].columns.tolist()
        print ('n clusters shared by %s or more samples in all samples= ' % nShared_all, len(all_percClusters))
    
        # concatenate lists:
        ClusterListShared = list(set(group1_percClusters + group2_percClusters + all_percClusters))
        print ('final number of clusters to use= ', len(ClusterListShared))
        
        # filter table to have only those sequences
        sampleByClusterDF_perc = sampleByClusterDF.loc[:, ClusterListShared]
        print ('final sampleByClusterDF_perc shape is: ', sampleByClusterDF_perc.shape)
        # filter table to not include clusters that appear in more than percTooMany of samples:
        percTooMany = command_args.maxSharingPerc
        nTooMany = math.ceil(len(sampleByClusterDF) * percTooMany)
        print 'filtering out clusters that appear in more than %s (%s) samples' % (nTooMany, percTooMany)
        sampleByClusterDF_perc = sampleByClusterDF_perc.loc[:, sampleByClusterDF_perc.count() < nTooMany]
        print ('final sampleByClusterDF_perc shape is: ', sampleByClusterDF_perc.shape)

        

        try:
            f6 = command_args.output_dir + 'sampleByClusterDF_cohortFiltering%s-%sperc.dat' % (command_args.fixed_filtering_perc_per_cohort,
                                                                                        command_args.maxSharingPerc)
            f6 = f6.replace('0.', '0')
            sampleByClusterDF_perc.to_pickle(f6)
            print 'saved filtered sampleByClusterDF to pickle'
        except:
            print 'couldnt save filtered sampleByClusterDF to pickle'
        
        sampleByClusterDF = sampleByClusterDF_perc
    
    elif command_args.nShared is not None:
        print ('filtering sampleByClusterDF to include only clusters shared by %s or more samples' % command_args.nShared)
        filtered = sampleByClusterDF.loc[:, sampleByClusterDF.count() >= command_args.nShared]
        print ('filtered.shape: ', filtered.shape)
        print ('filtered memory usage: ', filtered.memory_usage(index=True).sum())
        
        # save filtered to csv and pickle:
        try:
            f7 = command_args.output_dir + 'sampleByClusterDF_nShared%s.dat' % command_args.nShared
            filtered.to_pickle(f7)
            print 'saved filtered sampleByClusterDF to pickle'
        except:
            print 'couldnt save filtered sampleByClusterDF to pickle'
            
#         try:
#             f8=command_args.output_dir+'sampleByClusterDF_nShared%s.csv' %command_args.nShared
#             filtered.to_csv(f8)
#             print 'saved filtered sampleByClusterDF to csv'
#         except:
#             print 'couldnt save filtered sampleByClusterDF to csv'
        
        sampleByClusterDF = filtered
        
    else:
        print ('not_filtered.shape: ', sampleByClusterDF.shape)
        print ('filtered memory usage: ', sampleByClusterDF.memory_usage(index=True).sum())
        try:
            f7 = command_args.output_dir + 'sampleByClusterDF.dat' 
            sampleByClusterDF.to_pickle(f7)
            print 'saved non-filtered sampleByClusterDF'
        except:
            print 'couldnt save non-filtered sampleByClusterDF'
        
    
    return sampleByClusterDF
        

def save_args_to_file(fileName, command_args):
    # save parameters to a text file:
    print ('saving parameters to file')
    ctime = strftime("%Y%m%d_%H%M%S", gmtime())
    argFile = command_args.output_dir + fileName + '_%s.txt' % ctime
    new_file = open(argFile, mode="w")
    for arg in vars(command_args): 
        new_file.write(str(arg) + ": " + str(getattr(command_args, arg)) + "\n")
    new_file.close()



def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('output_dir', help='Path to output directory', type=str, default=None)
    parser.add_argument('-path_to_AllUniqueWithCounts', help='Path to  all unique aa sequences in dataset-AllUniqueWithCounts',
                         type=str, default='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/sharingAnalysis/AllUniqueWithCounts')
    parser.add_argument('-seq_type', help=' "prod"/"nonProd"/"All" - whether to use only productive/ only non-productive or all sequences',
                         type=str, default='All')
    parser.add_argument('-sample', help='define specific sample name to analyze', type=str, default=None)
    parser.add_argument('-n_random_samples', help='number of samples to randomly sample for analysis, if <0, all samples will be used',
                         type=int, default=-1)
    parser.add_argument('-path_to_levDistMat', help='path to existing lev distance matrix, if None, will generate a new one', type=str, default=None)
    parser.add_argument('-levDistThresh', help='what is largest lev distance allowed in cluster, choose 1,2 or 3', type=int, default=1)
    parser.add_argument('-nSeqsInBatchLev', help='how many sequences per batch for levdist calculations', type=int, default=2000)
    parser.add_argument('-dropping_batch_size', help='how many clusters per batch for clustersToDrop\
calculations', type=int, default=250)
    parser.add_argument('-onlyConcat', help='if True, do not generate levdists but conctenate the ones existing in the output_dir', type=bool, default=False)
    
    parser.add_argument('-deleteLevdistDFs', help='if true, delete all levdist files after concatenating them', type=bool, default=False)
    parser.add_argument('-saveLevDist', help='if true, save full levdist table', type=bool, default=True)
    parser.add_argument('-nShared', help='number of samples to be shared by cluster', type=int, default=3)
    parser.add_argument('-path_to_clusterMatrix_dropped', help='path to clusterMatrix_dropped', type=str, default=None)
    parser.add_argument('-dropClusters', help='bolean, whether to drop clusters before assigning to sample',
                        type=bool, default=0)
    parser.add_argument('-tempDir', help='tempDir', type=str,
        default='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/temp_q_dir/')
    parser.add_argument('-samples_concat_batch_size', help='number of 30-samples dfs to concatenate\
in the first concatenating step', type=int, default=4)
    parser.add_argument('-fixed_filtering_perc_per_cohort', help='perc shared by healthy samples, patients\
or both together', type=float, default=None)
    parser.add_argument('-maxSharingPerc', help='maximal percent of samples to be shared by a cluster to be included',
                        type=float, default=0.9)
    command_args = parser.parse_args()
    
    return command_args

def load_seqID_table(command_args):
    print 'loading seqIDs table'
    for f in listdir(command_args.output_dir):
        if 'seqID' in f: seqIDs = pd.read_pickle(command_args.output_dir + f)
        
    return seqIDs
'''
this function is not called by any of the other functions in this module and should be called seperately
it takes sampleByClusterDF (or seqByClusterDF, that contains nans and not 0's!)
and finds all clusters that are over represented in sample group1 vs. sample group2 or vice-versa, using 
both fisher (forcluster presence/absence info) or MW (for number of sequences per cluster per sample info)
tests. 
'''
def compare_seq_clusters_between_populations(sampleByClusterDF, sampleGroup1, sampleGroup1name,
                                             sampleGroup2, sampleGroup2name, output_dir,
                                             percShared=0.05, percTooMany=0.95, nShared=None):

    from scipy.stats import fisher_exact, mannwhitneyu
    import math
    import datetime
    print ('nShared=', nShared)
    ####### note - sampleByClusterDF cab also be sampleBySequenceDF but make sure that the df contains nans and not 0's !!!!!#########
    
    # ##generate cluster list:
    
    # get clusters that appear in Xperc or nShared of samples in sampleGroup1:
    print ('perc Shared is ', percShared)
    sampleByClusterDF_group1 = sampleByClusterDF.loc[sampleGroup1, :]
    if nShared is None:
        nShared_group1 = math.ceil(len(sampleByClusterDF_group1) * percShared)
    else:
        nShared_group1 = nShared
    print ('nShared_group1= ', nShared_group1)
    group1_percClusters = sampleByClusterDF_group1.loc[:, sampleByClusterDF_group1.count() >= nShared_group1].columns.tolist()
    print ('n clusters shared by %s or more samples in group1= ' % nShared_group1, len(group1_percClusters))
    
    # get clusters that appear in Xperc of samples in sampleGroup2:
    sampleByClusterDF_group2 = sampleByClusterDF.loc[sampleGroup2, :]
    if nShared is None:
        nShared_group2 = math.ceil(len(sampleByClusterDF_group2) * percShared)
    else:
        nShared_group2 = nShared
    print ('nShared_group2= ', nShared_group2)
    group2_percClusters = sampleByClusterDF_group2.loc[:, sampleByClusterDF_group2.count() >= nShared_group2].columns.tolist()
    print ('n clusters shared by %s or more samples in group2= ' % nShared_group2, len(group2_percClusters))
    
    # get clusters that appear in Xperc of All samples:
    if nShared is None:
        nShared_all = math.ceil(len(sampleByClusterDF) * percShared)
    else:
        nShared_all = nShared * 2
    print ('nShared_all= ', nShared_all)
    all_percClusters = sampleByClusterDF.loc[:, sampleByClusterDF.count() >= nShared_all].columns.tolist()
    all_percClusters = []
    print ('n clusters shared by %s or more samples in all samples= ' % nShared_all, len(all_percClusters))
    all_percClusters = []
    
    
    # concatenate lists:
    ClusterListShared = list(set(group1_percClusters + group2_percClusters + all_percClusters))
    print ('total number of clusters shared by more than %s ' % percShared, len(ClusterListShared))
    
    # filter table to have only those sequences
    sampleByClusterDF_perc = sampleByClusterDF.loc[:, ClusterListShared]
    # filter table to not include clusters that appear in more than 95% of samples:
    nTooMany = math.ceil(len(sampleByClusterDF) * percTooMany)
    print 'filtering out clusters that appear in more than %s (%s) samples' % (nTooMany, percTooMany)
    sampleByClusterDF_perc = sampleByClusterDF_perc.loc[:, sampleByClusterDF_perc.count() < nTooMany]
    print ('final sampleByClusterDF_perc shape is: ', sampleByClusterDF_perc.shape)

    # ## cluster comparison between groups:
    
    comparison_df = pd.DataFrame()
    sampleByClusterDF_perc = sampleByClusterDF_perc.fillna(0)
    sampleByClusterDF_perc['cohort'] = np.where(sampleByClusterDF_perc.index.isin(sampleGroup1), 0,
                np.where(sampleByClusterDF_perc.index.isin(sampleGroup2), 1, np.nan))
    sampleByClusterDF_perc = sampleByClusterDF_perc[sampleByClusterDF_perc['cohort'].notnull()]
    print ('n samples in group1= ', len(sampleByClusterDF_perc) - sampleByClusterDF_perc.cohort.sum())
    print ('n samples in group2= ', sampleByClusterDF_perc.cohort.sum())
    sampleByClusterDF_perc_binary = (sampleByClusterDF_perc > 0).astype(int)

    print ('start comparing, time is:', datetime.datetime.now())
    for n, cluster in enumerate(sampleByClusterDF_perc_binary.columns[:-1]):
        
        if n % 500 == 0: print n
        # fisher test:
        clusterDF = sampleByClusterDF_perc_binary[[cluster, 'cohort']]

        tab = pd.crosstab(clusterDF[cluster], clusterDF['cohort'], dropna=False).fillna(0)
        try:
            OR_f, p_f = fisher_exact(tab, alternative='two-sided')
        except:
            print ('couldnt execute fisher test for seq #%s' % n)
            print (tab)
            OR_f = 999999
            p_f = 1

        comparison_df.loc[n, 'cluster_head'] = cluster
        comparison_df.loc[n, 'OR_Fisher'] = OR_f
        comparison_df.loc[n, 'p_Fisher'] = p_f
        comparison_df.loc[n, 'absent_%s' % sampleGroup1name] = tab.iloc[0, 0]
        comparison_df.loc[n, 'absent_%s' % sampleGroup2name] = tab.iloc[0, 1]
        comparison_df.loc[n, 'present_%s' % sampleGroup1name] = tab.iloc[1, 0]
        comparison_df.loc[n, 'present_%s' % sampleGroup2name] = tab.iloc[1, 1]

        # MW test:
        groups = sampleByClusterDF_perc.groupby('cohort')[cluster]
        data = {}
        for name, group in groups:
            data[name] = group.tolist()
        data0 = data[0]
        data1 = data[1]

        s_MW, p_MW = mannwhitneyu(data0, data1)

        comparison_df.loc[n, 's_MW'] = s_MW
        comparison_df.loc[n, 'p_MW'] = p_MW

    print ('finished comparing, time is:', datetime.datetime.now())
    # add corrected p-values for multiple comparisons:
    nTests = len(comparison_df)
    FDR = 0.1
    comparison_df = comparison_df.sort_values(by='p_Fisher')

    for test in ['Fisher', 'MW']:
        comparison_df = add_corrected_pValues(comparison_df, pValueColumn='p_%s' % test, nTests=nTests, FDR=FDR).\
    rename(columns={'Sig by bonferroni corrected pVal':'Sig by bonferroni corrected pVal_%s' % test, 'sig. by FDR=%s' % FDR:'sig. by FDR=%s_%s' % (FDR, test),
                   'corr p_values by FDR=%s' % FDR:'corr p_values by FDR=%s_%s' % (FDR, test)})
    
    # add sequence identity to cluster head/sequence:
    CDR3_identity = pd.read_excel('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR CDR3 sequence databases/\
combined annotation_list_clean_popped.xlsx')
    
    comparison_df = pd.merge(comparison_df, CDR3_identity, how='left', left_on='cluster_head',
                                right_index=True)
    
    if not isdir(output_dir):makedirs(output_dir)
    if nShared is None:
        comparison_df_file = output_dir + 'Fisher_MW_results_%s_%s_percShared%s_percTooMany%s.xlsx'\
% (sampleGroup1name, sampleGroup2name, percShared, percTooMany)
    else:
        comparison_df_file = output_dir + 'Fisher_MW_results_%s_%s_nShared%s_percTooMany%s.xlsx'\
% (sampleGroup1name, sampleGroup2name, nShared, percTooMany)
    comparison_df_file = comparison_df_file.replace('0.', '0')
    comparison_df.to_excel(comparison_df_file)
    print ('comparison df was saved...')
    
    return comparison_df




def compare_seq_clusters_wrapper(sampleByClusterDF, sampleGroup1, sampleGroup1name, sampleGroup2, sampleGroup2name, output_dir,
                                 percShared=0.25, percTooMany=0.95, nShared=None):
    
   # imports and definitions:
    import math
    MyPath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF'   
    f1 = '/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/TargetDFs/PNP530Cardio126_AgeGender.xlsx' 
    PNP530Cardio126_AgeGender = pd.read_excel(f1).set_index('BD')

    All = sampleGroup1 + sampleGroup2

    # calculate fisher and man-whitnney:
    sampleByClusterDF = sampleByClusterDF.loc[All, :]

    print ('starting...')
    results = compare_seq_clusters_between_populations(sampleByClusterDF, sampleGroup1,
        sampleGroup1name, sampleGroup2, sampleGroup2name, output_dir, percShared, percTooMany, nShared)

    # compare age distributions:
    Age_sampleGroup1 = PNP530Cardio126_AgeGender.loc[sampleGroup1, 'Age'].dropna().tolist()
    Age_sampleGroup2 = PNP530Cardio126_AgeGender.loc[sampleGroup2, 'Age'].dropna().tolist()
    dataList = [(sampleGroup1name, Age_sampleGroup1), (sampleGroup2name, Age_sampleGroup2)]
    
    try:
        fig, ax = plt.subplots()
        # print (dataList)
        title = '%s vs. %s - Age distribution' % (sampleGroup1name, sampleGroup2name)
        plotHistComprison(dataList, ax, title, showLegend=True, nBins=20, toAnnotate=True, alpha=None, plotType='hist')
        plt.show()
    except:
        print ('cant do plots, sorry')

    # compare gender distribution:
    Gender_sampleGroup1 = Counter(PNP530Cardio126_AgeGender.loc[sampleGroup1, 'Gender_Male'].dropna().tolist())
    Gender_sampleGroup2 = Counter(PNP530Cardio126_AgeGender.loc[sampleGroup2, 'Gender_Male'].dropna().tolist())
    
    GenderDF = pd.concat([pd.DataFrame(index=[sampleGroup1name], data=Gender_sampleGroup1),
                        pd.DataFrame(index=[sampleGroup2name], data=Gender_sampleGroup2)]).fillna(0)
    print (GenderDF)
    
    try:
        R_f, p_f = fisher_exact(GenderDF, alternative='two-sided')
        print ('Fisher r,p = ', R_f, p_f)  
    except:
        print ('couldnt execute fisher test for gender distributions')
    
    return














'''
this function is not called by any of the other functions in this module and should be called seperately
'''

def get_cluster_seqeunces_per_clusterheadID(output_dir, head, isSeq=True):
    
    
    levDistMatTotalFile = output_dir + 'levDistMatTotal.dat'
    levDistMatTotal = pd.read_pickle(levDistMatTotalFile)
    clusterMatrix = levDistMatTotal.groupby(['seq1id', 'seq2id']).count()
    
    # get seqID matrix:
    f3 = output_dir + 'seqID_prod.dat'
    seqID = pd.read_pickle(f3)
    
    print ('replacing seqIDs with real sequences...')
    clusterMatrix_2 = clusterMatrix.reset_index()
    clusterMatrix_seqs = pd.merge(clusterMatrix_2, pd.DataFrame(seqID), how='left', left_on='seq2id', right_index=True).rename(columns={'seq':'seq2'})
    clusterMatrix_seqs = pd.merge(clusterMatrix_seqs, pd.DataFrame(seqID), how='left', left_on='seq1id', right_index=True).rename(columns={'seq':'seq1'})
    
    clusterMatrix_seqs_groupByHeadID = clusterMatrix_seqs.groupby(['seq1id', 'seq1', 'seq2']).count()
    
    if isSeq: 
        level = 'seq1id'
        clusterMatrix_seqs_groupByHeadID = clusterMatrix_seqs.groupby(['seq1', 'seq1id', 'seq2']).count()
    else:
        level = 'seq1'
        clusterMatrix_seqs_groupByHeadID = clusterMatrix_seqs.groupby(['seq1id', 'seq1', 'seq2']).count()
        
    clusterhead = clusterMatrix_seqs_groupByHeadID.loc[head, :].index.get_level_values(level)[0]
    seqs = clusterMatrix_seqs_groupByHeadID.loc[head, :].index.get_level_values('seq2').tolist()
    
    return clusterhead, seqs

'''
the following function is not called by the module
it was used directly from iPython on high CPU machines to concatenate 4 dfs of 30 samples into 1 df.
the resulting dfs should be further concatenate to yield eventually the full dataframe combined from 
22 dfs of ~30 samples each.
'''

#---------------------------------------------------------
'''
the following function (and its helper function) takes a df with a column that
contains cluster heads, or takes a list of cluster heads (and then insert df=None)
and finds all the sequences in this cluster.
optionally, it summarizes the annotations of the sequences in the cluster

'''
def gen_string_from_dict(dict1):
    a=''  
    for k,v in dict1.items():
        i=':'.join([str(k),str(v)])
        a=','.join([a,i])
        
    return a

def get_annot_for_seqs_in_cluster(df,cluster_head_col,output_dir,
                        getIdentities=True,cluster_head_list=None,identity_file=None,sumIdentities=True):
   
    from scipy.stats import chi2_contingency
    from collections import Counter
    
    if df is None:
        df=pd.DataFrame({'cluster_head':cluster_head_list})
        cluster_head_col='cluster_head'
        print df
        
    if identity_file is None:
        identity_file='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR CDR3 sequence databases/combined annotation_list_clean_popped.xlsx'
    CDR3_identity=pd.read_excel(identity_file)
    df2=df.copy()
    levDistMatTotalFile=output_dir+'levDistMatTotal.dat'
    levDistMatTotal=pd.read_pickle(levDistMatTotalFile)
    
    #get seqID matrix:
    print ('getting seqID info')
    for f in listdir(output_dir):
        if 'seqID' in f: seq_ID_file=f
    f3=output_dir+seq_ID_file
    seqID=pd.read_pickle(f3)
        
    #get sequences in cluster and their identities:
    print ('getting info for each cluster...')
    phen_related_identities={};
       
    for row in df.index:
        print (row)
        cluster_head=df.loc[row,cluster_head_col]
        
       
        #get seqID for cluster head:
        try:
            head_ID=list(np.where(seqID == cluster_head)[0])[0]
        except: 
            print ('this row is repeated more than once?')
            head_ID=list(np.where(seqID == cluster_head.tolist()[0])[0])[0]
        #get IDs for all sequences in the cluster:
        cluster_seq_IDs=levDistMatTotal[levDistMatTotal['seq1id']==head_ID]['seq2id'].tolist()

        #get sequences for all sequences in the cluster:
        cluster_seq_sequences=[seqID.loc[x] for x in cluster_seq_IDs]

    #   if get identities:
        df2.loc[row,'n sequences in cluster']=len(cluster_seq_sequences)+1
        df2.loc[row,'cluster sequences']=','.join(cluster_seq_sequences)
        if getIdentities:
            cluster_seq_sequences.insert(0,cluster_head)
            annot_list=[]
            for seq in cluster_seq_sequences:
                try:
                    annot=CDR3_identity.loc[seq,'combined annotation_list_clean']
                    try:
                        print (annot.shape) #a way to detect whether a sequence has more than 1 annotation (then it will become a series here)
                        print ('sequence %s has more than 1 annotations' %seq)
                        print (annot)
                        for a1 in list(annot):
                            annot_list.append(a1)
                    except:
                        annot_list.append(annot)
                except:
                    annot_list.append(None)

            identities_count=Counter(annot_list)
            for k,v in identities_count.items():
                try:
                    phen_related_identities[k]=phen_related_identities[k]+v
                except:
                    phen_related_identities[k]=v
                       
            a=gen_string_from_dict(identities_count)
            df2.loc[row,'annotation summary']=a
            
    if sumIdentities:
        phen_related_identities={i:phen_related_identities[i] for i in phen_related_identities if i!=None}
        df2.loc['summary','annotation summary']=gen_string_from_dict(phen_related_identities)
                                 
    else:
        phen_related_identities=None
    return df2,phen_related_identities
    
#------------------------------------------------------------
'''
the followig function crop clusterDF to includes only clusters shared by more than
min_share percent of samples and less than max_share percent of samples.
'''


def get__limited_sharing(clusterDF,min_share,max_share):
    
    print ('initial shape is: ',clusterDF.shape)
    
    min_share_number=math.ceil(len(clusterDF) * min_share)
    max_share_number=math.ceil(len(clusterDF) * max_share)
    
    print ('min_share_number: ',min_share_number)
    print ('max_share_number: ',max_share_number)
    
    clusterDF_limited=clusterDF.loc[:, clusterDF.count()>=min_share_number]
    clusterDF_limited=clusterDF_limited.loc[:, clusterDF.count()<max_share_number]
    
    print ('final shape is: ',clusterDF_limited.shape)
    clusterDF_limited=clusterDF_limited.fillna(0)
    
    return clusterDF_limited
#----------------------------------------------------------
def conc4(first_sample):
    import datetime
    from os import listdir, mkdir, makedirs
    from os.path import isfile, join, isdir, exists 
    print ('start processing, time is ', datetime.datetime.now())
    now = datetime.datetime.now()

    folder = '/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/sharingAnalysis/seqClusters_allProd_maxdist1/'
    conc4_folder = folder + 'conc4_df/'
    df_folder = folder + 'sampleDFdir/'

    if not isdir(conc4_folder):
        makedirs(conc4_folder)
    sample_list = range(first_sample, first_sample + 150, 30)
#     sample_list=[656 for x in sample_list if x==660 else x]

    df_conc2_list = []
    for i in range(len(sample_list)):
        if i + 2 > (len(sample_list)) - 1:
            print ('i=%s, reached the end, quitting' % i)
            break
        if i % 2 == 0:
            print 
            f1 = df_folder + 'sampleByClusterDF_temp_%s_%s.dat' % (sample_list[i], sample_list[i + 1])
            f2 = df_folder + 'sampleByClusterDF_temp_%s_%s.dat' % (sample_list[i + 1], sample_list[i + 2])
            print ('start concatenating 2dfs, time is:', datetime.datetime.now())
            print (f1)
            print (f2)
            df1 = pd.read_pickle(f1)
            df2 = pd.read_pickle(f2)
            df_conc2 = pd.concat([df1, df2])
            df_conc2_list.append(df_conc2)
            conc2_file1 = conc4_folder + 'conc2_df_%s_%s.dat' % (sample_list[i - 2], sample_list[i])
            df_conc2.to_pickle(conc2_file1)
        if i % 4 == 2:
            print ('start concatenating 4dfs, time is :', datetime.datetime.now())
            f3 = conc4_folder + 'conc4_df_%s_%s.dat' % (sample_list[i - 2], sample_list[i + 2])
            print (f3)
           
            df_conc4 = pd.concat(df_conc2_list)
            df_conc2_list = []
            
            print ('saving conc4, time is: ', datetime.datetime.now())
            print (f3)
            try:
                df_conc4.to_pickle(f3)
            except:
                print ('couldnt save conc4')
                
    return


def main():
    print ('main')
    command_args = get_args()

    # check folder existance: 
    if not os.path.exists(command_args.path_to_AllUniqueWithCounts):
        print ("path_to_AllUniqueWithCounts doesn't exist!"); return
    #     make_dir_if_not_exists(command_args.output_dir)
    if not os.path.isdir(command_args.output_dir):
        os.makedirs(command_args.output_dir)
        print ('made new output dir')
    save_args_to_file('temp_argFile', command_args)
    tempDir = command_args.tempDir
    if not isdir(tempDir):
        makedirs(tempDir)
    
    if command_args.path_to_clusterMatrix_dropped is None:
    
        # generate/get lev distance matrix:
        if command_args.onlyConcat:
            levDistMatTotal = concat_distmats(command_args)
            seqIDs = load_seqID_table(command_args)
        elif command_args.path_to_levDistMat is not None:
            print 'loading existing lev distance matrix...'
            levDistMatTotal = pd.read_pickle(command_args.path_to_levDistMat)
            seqIDs = load_seqID_table(command_args)
        else:   
            UniqueSeqs, name_extension = generate_seqTable(command_args)
            seqIDs = gen_levdist_table(command_args, UniqueSeqs, name_extension)
            levDistMatTotal = concat_distmats(command_args)
            
        # convert lev distance matrix to clustermatrix:
        print 'convert lev distance matrix to cluster matrix...'
        clusterMatrix_multiSeqs, getClustersInfo = convert_levDistMat_to_ClusterMatrix(levDistMatTotal, command_args)
        if command_args.dropClusters:
            print 'dropping clusters...'
            clusterMatrix_multiSeqs_dropped = drop_clusters(command_args, getClustersInfo, clusterMatrix_multiSeqs)
            clusterMatrix_multiSeqs = clusterMatrix_multiSeqs_dropped
    else:
        seqIDs = load_seqID_table(command_args)
        clusterMatrix_multiSeqs, getClustersInfo = convert_levDistMat_to_ClusterMatrix(levDistMatTotal, command_args)
        print 'loading clusterMatrix_dropped...'
        clusterMatrix_multiSeqs_dropped = pd.read_pickle(command_args.path_to_clusterMatrix_dropped)
        clusterMatrix_multiSeqs = clusterMatrix_multiSeqs_dropped
        
     # use the stacked dropped table to get data per sample:
    print 'get clusters for samples, dropClusters=%s' % command_args.dropClusters
    sampleByClusterDF = get_clusters_for_samples(clusterMatrix_multiSeqs, seqIDs, command_args)
         
    if not command_args.dropClusters:
        clustersTouse = sampleByClusterDF.columns
        clusterToUseIDs = seqIDs.reset_index().reset_index().set_index('seq').loc[clustersTouse]['index'].tolist()
        print ('clusterMatrix_multiSeqs.shape[0]: ', clusterMatrix_multiSeqs.shape[0])
        print ('number of shared clusters before dropping is (should be smaller: ', len(clusterToUseIDs))
        print ('number of *unique* shared clusters before dropping is (should be equal): ', len(set(clusterToUseIDs)))
        if clusterMatrix_multiSeqs.shape[0] > len(clusterToUseIDs):
            print 'filtering clusterMatrix to include only shared clusters before dropping clusters'
            clusterMatrix_multiSeqs = clusterMatrix_multiSeqs[clusterMatrix_multiSeqs.index.isin(clusterToUseIDs, level=0)]
            getClustersInfo = getClustersInfo[getClustersInfo.index.isin(clusterToUseIDs, level=0)]
        else:
            print 'number of shared clusters is equal to number of clusters in clusterMatrix, start dropping clusters...'
        clusterMatrix_multiSeqs_dropped = drop_clusters(command_args, getClustersInfo, clusterMatrix_multiSeqs)
        final_clusters = clusterMatrix_multiSeqs_dropped.index.get_level_values('seq1id').unique().tolist()
        final_cluster_seqeunces = seqIDs.iloc[final_clusters].tolist()
        print ('number of clusters after sharing and dropping =', len(final_cluster_seqeunces))
        print ('the following list should contain sequences: ', final_cluster_seqeunces[:5])
        sampleByClusterDF_final = sampleByClusterDF[final_cluster_seqeunces]
        print ('sampleByClusterDF shape: ', sampleByClusterDF.shape)
        print ('sampleByClusterDF_final shape: ', sampleByClusterDF_final.shape)
        
        if command_args.fixed_filtering_perc_per_cohort is not None: 
            extension = '_cohortfiltering%s-%sperc' % (command_args.fixed_filtering_perc_per_cohort,
                                                    command_args.maxSharingPerc)
            extension = extension.replace('.', '')
        elif command_args.nShared is not None: extension = '_nShared%s' % command_args.nShared
        else: extension = ''
    
        f7 = command_args.output_dir + 'sampleByClusterDF%s_dropped.dat' % extension
        try:
            sampleByClusterDF_final.to_pickle(f7)
            print 'sampleByClusterDF_shared_dropped'
        except:
            print 'couldnt save sampleByClusterDF_shared_dropped'
       
    print 'thats all for now!' 
    
    return
        
               
if __name__ == "__main__":
#     sethandlers()
    main()
