'''
this code analyzes the sequences within replicates of the same samples, and is written for samples that have 
2X 1ug and 2X0.5ug replicates. 
***sould be generalized, also add option for aa sequence analysis***
'''

import pandas as pd
import numpy as np
from scipy.stats.stats import pearsonr
# from scipy.stats import mannwhitneyu
import math
# import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# import collections
# from itertools import chain
from Utils import Load, Write
from venn import venn
from scipy.stats import mannwhitneyu
import cPickle as pickle
from myplots import replicate_scatter, seq_freq_calc
# from datetime import datetime

#---------------------------------------------------------------------------------------------------

def generate_SeqLists(sample):
# # for samples in the sample list, reads the TSV data for each replicate (2X0.5ug, 2X1ug) and convert it to a seperate dataframe
    global seq_df, a05_df, b05_df, a1_df, b1_df, comb05_df, comb1_df
    a05_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sa_05ug.tsv" % sample)
    b05_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sb_05ug.tsv" % sample)
    comb05_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%s_05ug.tsv" % sample)
    a1_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sa_1ug.tsv" % sample)
    b1_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sb_1ug.tsv" % sample)
    comb1_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%s_1ug.tsv" % sample)
    a05_df = a05_df.set_index('nucleotide')
    b05_df = b05_df.set_index('nucleotide')
    comb05_df = comb05_df.set_index('nucleotide')
    a1_df = a1_df.set_index('nucleotide')
    b1_df = b1_df.set_index('nucleotide')
    comb1_df = comb1_df.set_index('nucleotide')
# # generate lists of unique nucleotide sequences for each replicate, set of all sequences in all samples, and combined list for 0.5ug and 1ug samples:
    global seq_a05, seq_b05, seq_a1, seq_b1, all_seq_set, seq_05, seq_1, seq_comb05, seq_comb1
    seq_a05 = [str(i) for i in list(a05_df.index)]
#     print len(seq_a05)
    seq_b05 = [str(i) for i in list(b05_df.index)]
#     print len(seq_b05)
    seq_comb05 = [str(i) for i in list(comb05_df.index)]
    seq_a1 = [str(i) for i in list(a1_df.index)]
#     print len(seq_a1)
    seq_b1 = [str(i) for i in list(b1_df.index)]
#     print len(seq_b1)
    seq_comb1 = [str(i) for i in list(comb1_df.index)]
    all_seq = seq_a05 + seq_b05 + seq_a1 + seq_b1
#     print len(all_seq)
    all_seq_set = list(set(all_seq))
#     print len(all_seq_set)
    seq_05 = seq_a05 + seq_b05
    seq_1 = seq_a1 + seq_b1
#------------------------------------------------------------------------------
def unique_rearrang_bar():
    quad_samples = [139, 207, 418, 415]
    dup_samples = [202, 398, 419, 148]
    total_samples = quad_samples + dup_samples
    a05_n = []
    b05_n = []
    comb05_n = []
    a05b05_n = []
    a1_n = []
    b1_n = []
    comb1_n = []
    a1b1_n = []
    
    
    for sample in total_samples:
        a05_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sa_05ug.tsv" % sample)
        b05_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sb_05ug.tsv" % sample)
        comb05_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%s_05ug.tsv" % sample)
        seq_a05 = [str(i) for i in list(a05_df['nucleotide'])]
        seq_b05 = [str(i) for i in list(b05_df['nucleotide'])]
        seq_comb05 = [str(i) for i in list(comb05_df['nucleotide'])]
#         seq_a05b05=set(seq_a05+seq_b05)
        n_seq_a05 = len(set(seq_a05))
        n_seq_b05 = len(set(seq_b05))
        n_seq_comb05 = len(set(seq_comb05))
        n_seq_a05b05 = n_seq_a05 + n_seq_b05
        a05_n.append(n_seq_a05)
        b05_n.append(n_seq_b05)
        comb05_n.append(n_seq_comb05)
        a05b05_n.append(n_seq_a05b05)
       
        
        if sample in quad_samples:
            a1_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sa_1ug.tsv" % sample)
            b1_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sb_1ug.tsv" % sample)
            comb1_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%s_1ug.tsv" % sample)
            seq_a1 = [str(i) for i in list(a1_df['nucleotide'])]
            seq_b1 = [str(i) for i in list(b1_df['nucleotide'])]
            seq_comb1 = [str(i) for i in list(comb1_df['nucleotide'])]
#             seq_a1b1=set(seq_a1+seq_b1)
            n_seq_a1 = len(set(seq_a1))
            n_seq_b1 = len(set(seq_b1))
            n_seq_comb1 = len(set(seq_comb1))
            n_seq_a1b1 = n_seq_a1 + n_seq_b1
            a1_n.append(n_seq_a1)
            b1_n.append(n_seq_b1)
            comb1_n.append(n_seq_comb1)
            a1b1_n.append(n_seq_a1b1)
        else:
            a1_n.append(0.0)
            b1_n.append(0.0)
            comb1_n.append(0.0)
            a1b1_n.append(0.0)
    print a05_n
    print b05_n
    print comb05_n
    print a05b05_n
    fig = plt.figure()
    space = 0.2
    width = (1 - space) / 8
    indeces = range(len(total_samples))
    pos_a05 = [i - 3.5 * width for i in indeces]
    pos_b05 = [i - 2.5 * width for i in indeces]
    pos_comb05 = [i - 1.5 * width for i in indeces]
    pos_a05b05 = [i - 0.5 * width for i in indeces]
    pos_a1 = [i + 0.5 * width for i in indeces]
    pos_b1 = [i + 1.5 * width for i in indeces]
    pos_comb1 = [i + 2.5 * width for i in indeces]
    pos_a1b1 = [i + 3.5 * width for i in indeces]
 
    plt.bar(pos_a05 , a05_n, align='center', label='a_0.5ug', width=width, color='blue')
    plt.bar(pos_b05 , b05_n, align='center', label='b_0.5ug', width=width, color='blue')
    plt.bar(pos_comb05 , comb05_n, align='center', label='comb_0.5ug', width=width, color='green')
    plt.bar(pos_a05b05 , a05b05_n, align='center', label='a05b05', width=width, color='black') 
    plt.bar(pos_a1 , a1_n, align='center', label='a_1ug', width=width, color='red')
    plt.bar(pos_b1 , b1_n, align='center', label='b_1ug', width=width, color='red')
    plt.bar(pos_comb1 , comb1_n, align='center', label='comb_1ug', width=width, color='yellow')
    plt.bar(pos_a1b1 , a1b1_n, align='center', label='a1b1', width=width, color='pink') 
    plt.xticks(indeces, total_samples, fontsize=9)
    # plt.yscale('log')
    # plt.ylim(100000,10000000)
    # plt.yticks(np.logspace(5,7))
    plt.title('Number of Rearrangements Per sample and DNA amount')
    plt.ylabel('Unique Rearrangements')
    plt.xlabel('Sample')
    plt.margins(0.02, 0)
    plt.legend(loc='upper right', fontsize=9)
    
    # plt.annotate('Test', xy=(1, 0), xycoords='axes fraction', fontsize=16,
    #                 horizontalalignment='right', verticalalignment='bottom')
    # plt.text(0.5, 0.5, 'Mean-0.5ug samples= %s' % mean_05, ha='center', va='center',transform=transAxes, fontsize=14)
    # plt.text(0.91, 0.91 'Mean-1ug samples=%s' % mean_1, fontsize=14)
    
    if create_pdf:
        print'hello'
        with PdfPages('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Images/Unique_Rearrange_all_samples_041216.pdf') as pdf:
            pdf.savefig(fig)
        pdf.close
    else:
        plt.show()
            
 #-------------------------------------------------------------------------------------------------------------
 
def template_number_bar():
    quad_samples = [139, 207, 418, 415]
    dup_samples = [202, 398, 419, 148]
    total_samples = quad_samples + dup_samples
    a05_n = []
    b05_n = []
    comb05_n = []
    a05b05_n = []
    a1_n = []
    b1_n = []
    comb1_n = []
    a1b1_n = []
    
    
    for sample in total_samples:
        a05_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sa_05ug.tsv" % sample)
        b05_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sb_05ug.tsv" % sample)
        comb05_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%s_05ug.tsv" % sample)
        n_temp_a05 = a05_df['count (templates/reads)'].sum()
        n_temp_b05 = b05_df['count (templates/reads)'].sum()
        n_temp_comb05 = comb05_df['count (templates/reads)'].sum()
        n_temp_a05b05 = n_temp_a05 + n_temp_b05
        a05_n.append(n_temp_a05)
        b05_n.append(n_temp_b05)
        comb05_n.append(n_temp_comb05)
        a05b05_n.append(n_temp_a05b05)
        
        
        if sample in quad_samples:
            a1_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sa_1ug.tsv" % sample)
            b1_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sb_1ug.tsv" % sample)
            comb1_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%s_1ug.tsv" % sample)
            n_temp_a1 = a1_df['count (templates/reads)'].sum()
            n_temp_b1 = b1_df['count (templates/reads)'].sum()
            n_temp_comb1 = comb1_df['count (templates/reads)'].sum()
            n_temp_a1b1 = n_temp_a1 + n_temp_b1
            a1_n.append(n_temp_a1)
            b1_n.append(n_temp_b1)
            comb1_n.append(n_temp_comb1)
            a1b1_n.append(n_temp_a1b1)
        else:
            a1_n.append(0.0)
            b1_n.append(0.0)
            comb1_n.append(0.0)
            a1b1_n.append(0.0)
    
    
    fig = plt.figure()
    space = 0.2
    width = (1 - space) / 8
    indeces = range(len(total_samples))
    pos_a05 = [i - 3.5 * width for i in indeces]
    pos_b05 = [i - 2.5 * width for i in indeces]
    pos_comb05 = [i - 1.5 * width for i in indeces]
    pos_a05b05 = [i - 0.5 * width for i in indeces]
    pos_a1 = [i + 0.5 * width for i in indeces]
    pos_b1 = [i + 1.5 * width for i in indeces]
    pos_comb1 = [i + 2.5 * width for i in indeces]
    pos_a1b1 = [i + 3.5 * width for i in indeces]
        
    plt.bar(pos_a05 , a05_n, align='center', label='a_0.5ug', width=width, color='blue')
    plt.bar(pos_b05 , b05_n, align='center', label='b_0.5ug', width=width, color='blue')
    plt.bar(pos_comb05 , comb05_n, align='center', label='comb_0.5ug', width=width, color='green')
    plt.bar(pos_a05b05 , a05b05_n, align='center', label='a05b05', width=width, color='black') 
    plt.bar(pos_a1 , a1_n, align='center', label='a_1ug', width=width, color='red')
    plt.bar(pos_b1 , b1_n, align='center', label='b_1ug', width=width, color='red')
    plt.bar(pos_comb1 , comb1_n, align='center', label='comb_1ug', width=width, color='yellow')
    plt.bar(pos_a1b1 , a1b1_n, align='center', label='a1b1', width=width, color='pink') 
    plt.xticks(indeces, total_samples, fontsize=9)
#     plt.yscale('log')
#     plt.ylim(100000,10000000)
#     plt.yticks(np.logspace(5,7))
    plt.title('Number of Templates Per Sample and DNA Amount')
    plt.ylabel('Total templates')
    plt.xlabel('Sample')
    plt.margins(0.02, 0)
    plt.legend(loc='upper right', fontsize=9)
    
#     plt.annotate('Test', xy=(1, 0), xycoords='axes fraction', fontsize=16,
#                     horizontalalignment='right', verticalalignment='bottom')
#     plt.text(0.5, 0.5, 'Mean-0.5ug samples= %s' % mean_05, ha='center', va='center',transform=transAxes, fontsize=14)
#     plt.text(0.91, 0.91 'Mean-1ug samples=%s' % mean_1, fontsize=14)
#   

    if create_pdf:
        print'hello'
        with PdfPages('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Images/Total_templates_all_samples_041216.pdf') as pdf:
            pdf.savefig(fig)
        pdf.close
    else:
        plt.show()
                       
        
        
        
#--------------------------------------------------------------------------------------------------------    

def gen_venns(sample):
# # for each sample, generate a figure with 3 subplots: 
# # (1) venn diagram of 4 ellipse, 1 for each replicate. (2) venn diagram for the 2 0.5 replicates.
# # (3) venn diagram for the 2 1 replicates.
# # currently saves each fgure in a seperate PDF
    fig = plt.figure(figsize=(5, 11))
    plt.suptitle('Sample %s' % sample, fontsize=18)
    ax = plt.subplot(3, 1, 1)
    venn(ax, [set(seq_a05), set(seq_b05), set(seq_a1), set(seq_b1)], ["a05", "b05", "a1", "b1"])
    ax = plt.subplot(3, 1, 2)
    venn(ax, [set(seq_a05), set(seq_b05)], ["a05", "b05"])
    ax = plt.subplot(3, 1, 3)
    venn(ax, [set(seq_a1), set(seq_b1)], ["a1", "b1"])
    
    
    print
    if create_pdf:
        with PdfPages('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Images/venn_diagrams_%s.pdf' % sample) as pdf:
            pdf.savefig(fig)
        pdf.close
    else:
        plt.show() 
#-------------------------------------------------------------------------------------------------------

def seq_intersections_stats(sample):
# # Calculates the percent of sequences that fill the following criteria, out of the total number of unique nucleotide sequences in the 4 replicates of the sample
# # sequences that were detected in 0.5ug samples (both together) but not in 1ug samples ,and vice versa
# # sequences that were detected only in one, but not the other, replicate of 0.5ug samples, and 1ug samples.
# # generate a dataframe summarizing this data for all samples and saves it into a csv file
    inter_a05_b05 = set(seq_a05) & set(seq_b05)
    add_a05_b05 = round(float(len(set(seq_a05) - inter_a05_b05)) / len(all_seq_set) * 100, 2)
    add_b05_a05 = round(float(len(set(seq_b05) - inter_a05_b05)) / len(all_seq_set) * 100, 2)
    inter_a1_b1 = set(seq_a1) & set(seq_b1)
    add_a1_b1 = round(float(len(set(seq_a1) - inter_a1_b1)) / len(all_seq_set) * 100, 2)
    add_b1_a1 = round(float(len(set(seq_b1) - inter_a1_b1)) / len(all_seq_set) * 100, 2)
    inter_1_05 = set(seq_05) & set(seq_1)
    add_1_05 = round(float(len(set(seq_1) - inter_1_05)) / len(all_seq_set) * 100, 2)
    add_05_1 = round(float(len(set(seq_05) - inter_1_05)) / len(all_seq_set) * 100, 2)
    
    seq_intersections_list.append({'sample' : sample, 'a05 not in b05' : add_a05_b05, 'b05 not in a05' : add_b05_a05, 'a1 not in b1' : add_a1_b1, 'b1 not in a1' : add_b1_a1, '1 samples not 05 samples': add_1_05, '05 samples not 1 samples': add_05_1})  
    print seq_intersections_list.append
    
#----------------------------------------------------------------------------------------------------------
   
def generate_SeqDF(sample):   
## ***a highly time consuming function!***
# #generates a dataframe for each samples, summarizing for each sequence in which sample(s) it exists. 
    print 'working on %s df' % sample
    
    exist_a05 = []
    exist_b05 = []
    exist_a1 = []
    exist_b1 = []
    
    [exist_a05.append(1) if seq in seq_a05 else exist_a05.append(0) for seq in all_seq_set]
    [exist_b05.append(1) if seq in seq_b05 else exist_b05.append(0) for seq in all_seq_set] 
    [exist_a1.append(1) if seq in seq_a1 else exist_a1.append(0) for seq in all_seq_set] 
    [exist_b1.append(1) if seq in seq_b1 else exist_b1.append(0) for seq in all_seq_set]   
    
        
    seq_df = pd.DataFrame({'Nucleotide': all_seq_set, 'Exist a05' : exist_a05, 'Exist b05' : exist_b05, 'Exist a1' : exist_a1, 'Exist b1' : exist_b1})
    seq_df['n samples exist'] = seq_df['Exist a05'] + seq_df['Exist b05'] + seq_df['Exist a1'] + seq_df['Exist b1']
    seq_df['exist in 0.5ug samples'] = seq_df['Exist a05'] + seq_df['Exist b05']
    seq_df['exist in 1ug samples'] = seq_df['Exist a1'] + seq_df['Exist b1']
    seq_df_file = Write('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/seq_df_%s.df' % sample, seq_df)
    
    print 'Finished working on %s df' % sample
#----------------------------------------------------------------------------------------------------------

def correlations():
    compl_sample_list = [202, 207, 398, 139, 419, 148, 418, 415]
    Lymph_status = ['High', 'Normal', 'Normal', 'High', 'High', 'Normal', 'Normal', 'High']
    Lymph_per = [44, 30, 24, 42, 47, 35, 21, 44]
    Gender = ['Female', 'Female', 'male', 'Female', 'male', 'male', 'Female', 'male']
    Age = [63, 61, 64, 29, 30, 31, 46, 39]
    
    n_seq_05 = []
    max_freq = []
    for sample in compl_sample_list:
        a05_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sa_05ug.tsv" % sample)
        b05_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/BD%sb_05ug.tsv" % sample)
        seq_a05 = [str(i) for i in list(a05_df['nucleotide'])]
        seq_b05 = [str(i) for i in list(b05_df['nucleotide'])]
        seq_05 = seq_a05 + seq_b05
        n_seq_05.append(len(set(seq_05)))
        max_freq_sample = max(list(a05_df['frequencyCount (%)']) + list(b05_df['frequencyCount (%)'])) 
        max_freq.append(max_freq_sample)   
    correl_df = pd.DataFrame({'Sample': compl_sample_list, 'Lymph Status': Lymph_status, 'Lymphocyte percentage':Lymph_per, 'Gender': Gender, 'Age': Age, 'n sequences in 05 samples': n_seq_05, 'Max frequency': max_freq})
    gender_means_nseq = correl_df.groupby(['Gender'])[['n sequences in 05 samples']].mean()
    gender_means_nseq_list = list(gender_means_nseq['n sequences in 05 samples'])
    female_nseq = correl_df.loc[correl_df['Gender'] == 'Female']['n sequences in 05 samples']
    male_nseq = correl_df.loc[correl_df['Gender'] == 'Male']['n sequences in 05 samples']
    gender_means_maxfreq = correl_df.groupby(['Gender'])[['Max frequency']].mean()
    gender_means_maxfreq_list = list(gender_means_maxfreq['Max frequency'])
    female_maxfreq = correl_df.loc[correl_df['Gender'] == 'Female']['Max frequency']
    male_maxfreq = correl_df.loc[correl_df['Gender'] == 'Male']['Max frequency']
    
    print gender_means_nseq
# #     print correl_df.groupby(['Gender']).describe()
# #     print correl_df.groupby(['Lymph Status']).describe()
#     sum_aa_seq_reads = sample_df.groupby(['aminoAcid'])[['count (reads)']].sum()
    
    fig = plt.figure(figsize=(8, 11))
    
    plt.subplot(3, 2, 1)
    plt.bar([1, 2], gender_means_nseq_list, align='center')
    plt.ylim(40000, 45000)
    plt.xticks([1, 2], ['Female', 'Male'], fontsize=10)
    plt.title('Gender effect on the number of \n unique sequences per sample', fontsize=12)
    plt.ylabel('Unique sequences')
    plt.margins(0.02, 0.02)
    MU_s_nseq, MU_p_nseq = mannwhitneyu(female_nseq, male_nseq, use_continuity=True, alternative='two-sided')
    plt.text(1, 44000, 'p-value (MW)=%s' % MU_p_nseq, color='red')
    
    plt.subplot(3, 2, 2)
    plt.bar([1, 2], gender_means_maxfreq_list, align='center')
#     plt.ylim(40000,45000)
    plt.xticks([1, 2], ['Female', 'Male'], fontsize=10)
    plt.title('Gender effect on the  \n maximal clonal frequency', fontsize=12)
    plt.ylabel('Maximal clonal frequency')
    plt.margins(0.02, 0.02)
    MU_s_maxfreq, MU_p_maxfreq = mannwhitneyu(female_maxfreq, male_maxfreq, use_continuity=True, alternative='two-sided')
    plt.text(0.8, 0.04, 'p-value (MW)=%s' % MU_p_maxfreq, color='red')
    
    plt.subplot(3, 2, 3)
    plt.scatter(Lymph_per, n_seq_05)
    plt. title('Correlation between %Lymphocytes in blood \n and number of unique sequences', fontsize=12)
    plt.xlabel('Lymphocytes %', fontsize=10)
    plt.ylabel('Unique sequences')
    plt.plot(Lymph_per, np.poly1d(np.polyfit(Lymph_per, n_seq_05, 1))(Lymph_per))
    r, p = pearsonr(Lymph_per, n_seq_05)
    r_round = round(r, 3)
    p_round = round(p, 6)
    print r_round, p
    plt.text(20, 60000, 'r2=%s' % r_round)
    plt.text(20, 55000, 'p-value=%s' % p_round, color='red')
    
    plt.subplot(3, 2, 4)
    plt.scatter(Lymph_per, max_freq)
    plt. title('Correlation between %Lymphocytes in blood \n and maximal clonal frequency', fontsize=12)
    plt.xlabel('Lymphocytes %', fontsize=10)
    plt.ylabel('Max clonal frequency')
    plt.plot(Lymph_per, np.poly1d(np.polyfit(Lymph_per, max_freq, 1))(Lymph_per))
    r, p = pearsonr(Lymph_per, max_freq)
    r_round = round(r, 3)
    p_round = round(p, 6)
    print r_round, p
    plt.text(20, 0.09, 'r2=%s' % r_round)
    plt.text(20, 0.08, 'p-value=%s' % p_round)
    
    plt.subplot(3, 2, 5)
    plt.scatter(Age, n_seq_05)
    plt.title('Correlation between age and \n number of unique sequences', fontsize=12)
    plt.xlabel('Age', fontsize=10)
    plt.ylabel('Unique sequences')
    plt.plot(Age, np.poly1d(np.polyfit(Age, n_seq_05, 1))(Age))
    r, p = pearsonr(Age, n_seq_05)
    r_round = round(r, 3)
    p_round = round(p, 6)
    print r_round, p
    plt.text(30, 25000, 'r2=%s' % r_round)
    plt.text(30, 20000, 'p-value=%s' % p_round)
    
    plt.subplot(3, 2, 6)
    plt.scatter(Age, max_freq)
    plt.title('Correlation between age and \n maximal clonal frequency', fontsize=12)
    plt.xlabel('Age', fontsize=10)
    plt.ylabel('Maximal clonal frequency')
    plt.plot(Age, np.poly1d(np.polyfit(Age, max_freq, 1))(Age))
    r, p = pearsonr(Age, max_freq)
    r_round = round(r, 3)
    p_round = round(p, 6)
    print r_round, p
    plt.text(30, 0.09, 'r2=%s' % r_round)
    plt.text(30, 0.08, 'p-value=%s' % p_round)
    
    plt.subplots_adjust(left=0.12, bottom=0.07, right=0.92, top=0.95, wspace=0.32, hspace=0.33)
      
    if create_pdf:
        print'hello'
        with PdfPages('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Images/correlations.pdf') as pdf:
            pdf.savefig(fig)
        pdf.close
    else:
        plt.show() 
#--------------------------------------------------------------------------------------------------------
''' 
(1) generates corresponding lists of template numbers for each seq in a1, b1 and comb1 samples
(adapted from the function 'template_count()'
(2) generates df from these lists, which also shows for each sequence if the sum of template number in a1 and b1 
equals the number in the combined sample or not, and then calculate statistics comparing sequences with equal
and non-equal counts. 
'''

def seq_temp_calc(sample):
    global temp_in_a1, temp_in_b1, temp_in_comb1
    print ('template number analysis for sample %s...' % sample)
    seq1_unique = list(set(seq_a1 + seq_b1 + seq_comb1))  # list of all sequences that exist in at least one of the 1ug samples
# # comment from here to avoid timewaste, when value are already stored in pickles:     
#     temp_in_a1=[] #list of template numbers for each sequence in sample a1
#     temp_in_b1=[] #list of template numbers for each sequence in sample b1
#     temp_in_comb1=[] #list of template numbers for each sequence in sample comb1
#        
#     for seq in seq1_unique:           
#         if seq in seq_a1 and seq in seq_b1 and seq in seq_comb1:
#             temp_in_a1.append(a1_df.loc[seq,'count (templates/reads)'])
#             temp_in_b1.append(b1_df.loc[seq,'count (templates/reads)'])
#             temp_in_comb1.append(comb1_df.loc[seq,'count (templates/reads)'])
#         elif seq in seq_a1 and seq in seq_b1 and seq not in seq_comb1:  
#             temp_in_a1.append(a1_df.loc[seq,'count (templates/reads)'])
#             temp_in_b1.append(b1_df.loc[seq,'count (templates/reads)'])
#             temp_in_comb1.append(0)
#         elif seq in seq_a1 and seq not in seq_b1 and seq not in seq_comb1:  
#             temp_in_a1.append(a1_df.loc[seq,'count (templates/reads)'])
#             temp_in_b1.append(0)
#             temp_in_comb1.append(0)
#         elif seq not in seq_a1 and seq in seq_b1 and seq in seq_comb1:  
#             temp_in_a1.append(0)
#             temp_in_b1.append(b1_df.loc[seq,'count (templates/reads)'])
#             temp_in_comb1.append(comb1_df.loc[seq,'count (templates/reads)'])
#         elif seq not in seq_a1 and seq not in seq_b1 and seq in seq_comb1:  
#             temp_in_a1.append(0)
#             temp_in_b1.append(0)
#             temp_in_comb1.append(comb1_df.loc[seq,'count (templates/reads)'])
#         elif seq in seq_a1 and seq not in seq_b1 and seq in seq_comb1:
#             temp_in_a1.append(a1_df.loc[seq,'count (templates/reads)'])
#             temp_in_b1.append(0)
#             temp_in_comb1.append(comb1_df.loc[seq,'count (templates/reads)'])
#         elif seq not in seq_a1 and seq in seq_b1 and seq not in seq_comb1:  
#             temp_in_a1.append(0)
#             temp_in_b1.append(b1_df.loc[seq,'count (templates/reads)'])
#             temp_in_comb1.append(0)
#             
#             
#     print ('started template count for sample %s...' %sample)
#     compare_counts=[]
#     for i in range(len(seq1_unique)):
#         if temp_in_a1[i]+temp_in_b1[i]==temp_in_comb1[i]:
#             compare_counts.append(1)
#         else:
#             compare_counts.append(0)
#       
# #since it is time consuming to generate this list, save them as pickle and next time load the pickles:
#     pickle.dump(temp_in_a1, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/temp_in_a1_%s.p' %sample,"wb" ))
#     pickle.dump(temp_in_b1, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/temp_in_b1_%s.p' %sample,"wb" ))
#     pickle.dump(temp_in_comb1, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/temp_in_comb1_%s.p' %sample,"wb" ))
#     pickle.dump(compare_counts, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/compare_counts_%s.p' %sample,"wb" ))
    temp_in_a1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/temp_in_a1_%s.p' % sample, "rb"))
    temp_in_b1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/temp_in_b1_%s.p' % sample, "rb"))
    temp_in_comb1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/temp_in_comb1_%s.p' % sample, "rb"))
    compare_counts = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/compare_counts_%s.p' % sample, "rb"))
    

    temp_counts_df = pd.DataFrame({'sequence':seq1_unique, 'a1 counts': temp_in_a1, 'b1_counts': temp_in_b1, 'comb1_counts' : temp_in_comb1, 'equal': compare_counts})
#     print temp_counts_df.groupby(['b1 counts']).describe()
    
    one_temp_only = []
    for i in range(len(seq1_unique)): 
        if temp_in_a1[i] == 1 and temp_in_b1[i] == 1 and temp_in_comb1[i]:
            one_temp_only.append(seq1_unique[i])
    
    print len(one_temp_only)
    print float(len(one_temp_only)) / len(seq1_unique) * 100
       

#---------------------------------------------------------------------------------------------------------
def generate_replicate_view(sample):
# # generate an image for each sample that contains: 
# # (1) 1ug sample comparison using scatter
# # (2) histograms of copy numbers per sequence, in each 1ug replicate
# # (3) rarefaction plots for each 1ug replicate

    f, ((ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9)) = plt.subplots(3, 3, sharey=False, figsize=(11, 8))
    plt.suptitle ('Sample %s - 1ug duplicate analysis' % sample, fontsize=16)
    


    seq_freq_calc(sample)
    replicate_scatter(sample, freq_in_a1, freq_in_b1, 'a1', 'b1', ax1)
    replicate_scatter(sample, freq_in_a1, freq_in_comb1, 'a1', 'comb1', ax2)
    replicate_scatter(sample, freq_in_b1, freq_in_comb1, 'b1', 'comb1', ax3)
    template_count(sample)
    hist_plot(sample, a1_temp_freq, 'a1', bins, xlimit, 'red', ax4)
    hist_plot(sample, b1_temp_freq, 'b1', bins, xlimit, 'green', ax5)
    hist_plot(sample, comb1_temp_freq, 'comb1', bins, xlimit, 'black', ax6)
    rarefaction_calc(sample)
    rarefaction_plot(sample, n_samp_seq_a1, 'a1', 'red', ax7)
    rarefaction_plot(sample, n_samp_seq_b1, 'b1', 'green', ax8)
    rarefaction_plot(sample, n_samp_seq_comb1, 'comb1', 'black', ax9)
    plt.subplots_adjust(left=0.1, bottom=0.06, right=0.95, top=0.88, wspace=0.43, hspace=0.48)
 
    if create_pdf:
        print'hello'
        with PdfPages('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Images/dup_analysis_with_comb_sample_%s.pdf' % sample) as pdf:
            pdf.savefig(f)
        pdf.close
    else:
        plt.show()



#----------------------------
def seq_freq_calc(sample):
# # generate scatter plots of the relative frequency of each sequence in each of the two 1ug samples
# # when the seqeunce exist only in one sample, the value in the other sample is 1e-6
    global freq_in_a1, freq_in_b1, freq_in_comb1, max_freq
    print ('replicate analysis for sample %s...' % sample)
    seq1_unique = list(set(seq_a1 + seq_b1 + seq_comb1))  # list of all sequences that exist in at least one of the 1ug samples
# # comment from here to avoid timewaste, when value are already stored in pickles:     
#     freq_in_a1=[] #list of frequencies for each sequence in sample a1
#     freq_in_b1=[] #list of frequencies for each sequence in sample b1
#     freq_in_comb1=[] #list of frequencies for each sequence in sample comb1
#       
#     for seq in seq1_unique:           
#         if seq in seq_a1 and seq in seq_b1 and seq in seq_comb1:
#             freq_in_a1.append(a1_df.loc[seq,'frequencyCount (%)'])
#             freq_in_b1.append(b1_df.loc[seq,'frequencyCount (%)'])
#             freq_in_comb1.append(comb1_df.loc[seq,'frequencyCount (%)'])
#         elif seq in seq_a1 and seq in seq_b1 and seq not in seq_comb1:  
#             freq_in_a1.append(a1_df.loc[seq,'frequencyCount (%)'])
#             freq_in_b1.append(b1_df.loc[seq,'frequencyCount (%)'])
#             freq_in_comb1.append(1e-6)
#         elif seq in seq_a1 and seq not in seq_b1 and seq not in seq_comb1:  
#             freq_in_a1.append(a1_df.loc[seq,'frequencyCount (%)'])
#             freq_in_b1.append(1e-6)
#             freq_in_comb1.append(1e-6)
#         elif seq not in seq_a1 and seq in seq_b1 and seq in seq_comb1:  
#             freq_in_a1.append(1e-6)
#             freq_in_b1.append(b1_df.loc[seq,'frequencyCount (%)'])
#             freq_in_comb1.append(comb1_df.loc[seq,'frequencyCount (%)'])
#         elif seq not in seq_a1 and seq not in seq_b1 and seq in seq_comb1:  
#             freq_in_a1.append(1e-6)
#             freq_in_b1.append(1e-6)
#             freq_in_comb1.append(comb1_df.loc[seq,'frequencyCount (%)'])
#         elif seq in seq_a1 and seq not in seq_b1 and seq in seq_comb1:
#             freq_in_a1.append(a1_df.loc[seq,'frequencyCount (%)'])
#             freq_in_b1.append(1e-6)
#             freq_in_comb1.append(comb1_df.loc[seq,'frequencyCount (%)'])
#         elif seq not in seq_a1 and seq in seq_b1 and seq not in seq_comb1:  
#             freq_in_a1.append(1e-6)
#             freq_in_b1.append(b1_df.loc[seq,'frequencyCount (%)'])
#             freq_in_comb1.append(1e-6)
      
    # since it is time consuming to generate this list, save them as pickle and next time load the pickles:
#     pickle.dump(freq_in_a1, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freq_in_a1_%s.p' %sample,"wb" ))
#     pickle.dump(freq_in_b1, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freq_in_b1_%s.p' %sample,"wb" ))
#     pickle.dump(freq_in_comb1, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freq_in_comb1_%s.p' %sample,"wb" ))
    freq_in_a1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freq_in_a1_%s.p' % sample, "rb"))
    freq_in_b1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freq_in_b1_%s.p' % sample, "rb"))
    freq_in_comb1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freq_in_comb1_%s.p' % sample, "rb"))
    print ('finished replicate analysis for sample %s...' % sample)
    
    max_freq = max(freq_in_a1 + freq_in_b1 + freq_in_comb1)
    print max_freq
#     plt.figure(figsize=(8, 8))
#     plt.suptitle('Sequences Frequency Correlation Between a1 and b1 Replicates')

    
#------------------------------------------------------------------------------------------

def seq_freq_correl(sample):
    
    freq_in_a1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freq_in_a1_%s.p' % sample, "rb"))
    freq_in_b1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freq_in_b1_%s.p' % sample, "rb"))
    freq_in_comb1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freq_in_comb1_%s.p' % sample, "rb"))
    
# #the following lines replace 1e-6 values (that represent zeros) with the other minimal value of the lists. 
# #in order to go back to normal calculation, set freq_in_a1_new=freq_in_a1 and so on.


    freq_in_a1_set = set(freq_in_a1)
    freq_in_a1_set.remove(min(freq_in_a1_set))
    min_a1 = min(freq_in_a1_set)
    freq_in_b1_set = set(freq_in_b1)
    freq_in_b1_set.remove(min(freq_in_b1_set))
    min_b1 = min(freq_in_b1_set)
    freq_in_a1_new = [min_a1 if i == 1e-06 else i for i in freq_in_a1]
    freq_in_b1_new = [min_b1 if i == 1e-06 else i for i in freq_in_b1]
    freq_in_comb1_new = freq_in_comb1
    
 # #log transformation for the values in the lists (natural log)   
    freq_in_a1_log = [np.log(i) for i in freq_in_a1_new]
    freq_in_b1_log = [np.log(i) for i in freq_in_b1_new]
    freq_in_comb1_log = [np.log(i) for i in freq_in_comb1_new]

# #generating lists of r and p values for pearson correlation. the lists are formatted into df 
# # and csv file in the main function:   
    r_a1comb1, p_a1comb1 = pearsonr(freq_in_a1_log, freq_in_comb1_log)
    r_b1comb1, p_b1comb1 = pearsonr(freq_in_b1_log, freq_in_comb1_log)
    r_a1comb1_list.append(r_a1comb1)
    p_a1comb1_list.append(float(p_a1comb1))
    r_b1comb1_list.append(r_b1comb1)
    p_b1comb1_list.append(float(p_b1comb1))
    
#----------------------------------------------------------------------------------------
def frequency_comparison(sample): 

    seq1_unique = list(set(seq_a1 + seq_b1 + seq_comb1))
    freq_in_a1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freq_in_a1_%s.p' % sample, "rb"))
    freq_in_b1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freq_in_b1_%s.p' % sample, "rb"))
    freq_in_comb1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/freq_in_comb1_%s.p' % sample, "rb"))

    seq_freq_df = pd.DataFrame({'sequence': seq1_unique, 'a1 freq': freq_in_a1, 'b1 freq': freq_in_b1, 'comb1 freq': freq_in_comb1})
    ranks = seq_freq_df.rank(axis=0, method='average', numeric_only=None, na_option='keep', ascending=False)
    seq_freq_df['a1 freq rank'] = ranks['a1 freq']
    seq_freq_df['b1 freq rank'] = ranks['b1 freq']
    seq_freq_df['comb1 freq rank'] = ranks['comb1 freq']
    seq_freq_df.sort(columns='comb1 freq rank', axis=0, ascending=True, inplace=True, kind='quicksort')
    seq_freq_df.to_csv('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/QC/Sequences frequencies in duplicates %s.csv' % sample)
    a1_freq_rank = list(seq_freq_df['a1 freq rank'])
    b1_freq_rank = list(seq_freq_df['b1 freq rank'])
    comb1_freq_rank = list(seq_freq_df['comb1 freq rank'])
    
    missed_in_b = 0
    missed_in_a = 0
    for x in range(0, len(comb1_freq_rank)):
        conditionA = a1_freq_rank[x] < 100 and b1_freq_rank[x] > 30000
        conditionB = b1_freq_rank[x] < 100 and a1_freq_rank[x] > 30000
        if conditionA:
            missed_in_b = missed_in_b + 1
        if conditionB:
            missed_in_a = missed_in_a + 1
    print missed_in_b
    print missed_in_a
    
    
def replicate_scatter(sample, data1, data2, name1, name2, ax): 
# #data1 and data2 are corresponding lists of frequencies of each sequence in the two replicates. 
# # name1 and name2 are the names of the datasets compared
# # note that the order is critical here, as these lists carry data for the *same* sequences.
   
    
    if ax is None:
        ax = plt.gca()
    scatter1 = ax.scatter(data1, data2, alpha=0.1)
    ax.set_xscale('log')
    ax.set_yscale('log')
    labels = ['', '0', '1e-5', '1e-4', '1e-3', '1e-2', '1e-1', 1] 
    ax.set_xticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1], minor=False)
    # # 1e-6 is defined as '0', as zero can't be plotted in log scale. 
    # # since I couldn't control the margins, I added e-7 so the '0' data will be spaced from the axes. 
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_yticks([1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1], minor=False)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlim(5e-7, 1)
    ax.set_ylim(5e-7, 1)
    ax.set_title('Sequences Frequency Correlation \n Between Duplicates', fontsize=12)
    ax.set_xlabel('%s (Sequence Frequency)' % name1, fontsize=9)
    ax.set_ylabel('%s (Sequence Frequency)' % name2, fontsize=9)
    ax.margins(0.2)
    return scatter1
    
#------------------------------------------------------------------------------------------------------#----------------------------------------------------------------------------------------


def template_count(sample):
   
# #generate lists of copy numbers per clone in each replicate
# #and also other parameters needed for the histogram
    global a1_temp_freq, b1_temp_freq, comb1_temp_freq, bins, xlimit
    a1_temp_freq = list(a1_df['count (templates/reads)'])
    b1_temp_freq = list(b1_df['count (templates/reads)'])
    comb1_temp_freq = list(comb1_df['count (templates/reads)'])
        
    
    
    a1_max = max(a1_temp_freq)
    b1_max = max(b1_temp_freq)
    comb1_max = max(comb1_temp_freq)
    total_max = max(a1_max, b1_max, comb1_max)
    xlimit = roundup(total_max, 100)
    print xlimit
    bins = np.linspace(0, xlimit, 50)
    print ('finished template count for sample %s...' % sample)
    
#---------------------------------------------------------------------------------------------------
def hist_plot(sample, data, name, bins, xlimit, color, ax):
# # generates histogram for the clonal copy number     
    if ax is None:
        ax = plt.gca()
    hist1 = ax.hist(a1_temp_freq, bins, label='%s' % name, color=color)
    ax.set_yscale('log')
    ax.set_xlim(0, xlimit)
    xlabels = range(0, xlimit, 100)
    ax.set_ylim(0.5, 100000)
#     ax.legend()
#     ax.set_margins(0.2,0.2)
    ax.set_xlabel('N templates per sequence', fontsize=9)
    ax.set_ylabel('N unique sequences', fontsize=9)
    ax.set_title('%s-Copy Number Distribution' % name, fontsize=12) 
    ax.set_yticks([1, 10, 100, 1000, 10000, 100000], minor=False)
    ax.set_yticklabels(['1', '10', '100', '10^3', '10^4', '10^5'], fontsize=8)
    ax.set_xticklabels(xlabels, fontsize=8) 
    return hist1

#-------------------------------------------------------------------------------------------  
def roundup(x, fold):
    return int(math.ceil(x / float(fold))) * fold   

#----------------------------------------------------------------------------------------------


def rarefaction_calc(sample):
    global samp_size, max_samp_size, n_samp_seq_a1, n_samp_seq_b1, n_samp_seq_comb1, ymax
    print ('rarefaction analysis for sample %s' % sample)

# # coment from here when values are already stored in pickles:
#     the number of repeats
    a1_n_rep = list(a1_df['count (templates/reads)'])
    b1_n_rep = list(b1_df['count (templates/reads)'])
    comb1_n_rep = list(comb1_df['count (templates/reads)'])
     
# #generate 'popped' list of sequences according to the number of repeats of each sequence:
    seq_a1_popped = []
    for i in range(0, len(seq_a1)):
        for j in range(0, a1_n_rep[i]):
            seq_a1_popped.append(seq_a1[i])
     
    seq_b1_popped = []
    for i in range(0, len(seq_b1)):
        for j in range(0, b1_n_rep[i]):
            seq_b1_popped.append(seq_b1[i])
    seq_comb1_popped = []
    for i in range(0, len(seq_comb1)):
        for j in range(0, comb1_n_rep[i]):
            seq_comb1_popped.append(seq_comb1[i])
     
     
    max_samp_size = min([len(seq_a1_popped), len(seq_b1_popped), len(seq_comb1_popped)])
    samp_size = range(0, max_samp_size, 100)
#     n_samp_seq_a1=[]
#     n_samp_seq_b1=[]
#     n_samp_seq_comb1=[]
#     for size in samp_size:

#         rand_seq_a1=np.random.choice(seq_a1_popped, size, replace=False)
#         rand_seq_b1=np.random.choice(seq_b1_popped, size, replace=False)
#         rand_seq_comb1=np.random.choice(seq_comb1_popped, size, replace=False)
#         seq_n_a1=len(set(rand_seq_a1))
#         seq_n_b1=len(set(rand_seq_b1))
#         seq_n_comb1=len(set(rand_seq_comb1))
#         n_samp_seq_a1.append(seq_n_a1)
#         n_samp_seq_b1.append(seq_n_b1)
#         n_samp_seq_comb1.append(seq_n_comb1)
# #         
# #since it is time consuming to generate theses lists, save them as pickles and next time load the pickles:
#     pickle.dump(n_samp_seq_a1, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/n_samp_seq_a1_%s.p' %sample,"wb" ))
#     pickle.dump(n_samp_seq_b1, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/n_samp_seq_b1_%s.p' %sample,"wb" ))
#     pickle.dump(n_samp_seq_comb1, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/n_samp_seq_comb1_%s.p' %sample,"wb" ))
    n_samp_seq_a1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/n_samp_seq_a1_%s.p' % sample, "rb"))
    n_samp_seq_b1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/n_samp_seq_b1_%s.p' % sample, "rb"))
    n_samp_seq_comb1 = pickle.load(open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/n_samp_seq_comb1_%s.p' % sample, "rb"))
    ymax = max(n_samp_seq_a1[-1], n_samp_seq_b1[-1], n_samp_seq_comb1[-1])
    print ymax
    print (n_samp_seq_a1[-1], n_samp_seq_b1[-1], n_samp_seq_comb1[-1])
    print ('finished rarefaction analysis for sample %s...' % sample)
    print samp_size         
    
def rarefaction_plot(sample, data1, name1, color1, ax):
    if ax is None:
        ax = plt.gca()
    curve1 = ax.plot(samp_size, data1, label=name1, color=color1)
  
#     curve2=ax.plot(samp_size,n_samp_seq_b1,label='b1',color='green')
#     curve3=ax.plot(samp_size,n_samp_seq_comb1,label='comb1',color='black')
    ax.set_title('Rarefaction Plot', fontsize=12)
    ax.set_xlabel('Sample Size', fontsize=9)
    ax.set_ylabel('Number of unique sequences', fontsize=9)
    xticks = range(0, roundup(max_samp_size, 10000) + 10000, 10000)
    yticks = range(0, roundup(ymax, 10000) + 10000, 10000)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(xticks, fontsize=8)
    ax.set_yticklabels(yticks, fontsize=8)
    ax.legend(loc='lower right')
    ax.set_xlim(0, roundup(max_samp_size, 10000))
    ax.set_ylim(0, roundup(ymax, 10000))
    return curve1
   
   
#--------------------------------------------------------------------------------------- 

    
    
    
    
    
#     seq_df=Load('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/QC/seq_df_%s.df' %sample)
#     print seq_df[0:20]
    
    

#-------------------------------------------------------------------------------

#---------------------------------------------------------

'''
main function:
(1) list the sample 
(2) define whether to generate pdf for the images
(3) define an empty list fot the sequence intersection which will be the basis for the intersection data 
dataframe and file in function 'seq_intersections_stats'
(4) calls functions to generate sequence lists, venn diagrams, statistics df and full sequence data dataframe
'''

sample_list = [139, 207, 415, 418]

# # general prepearations:
global create_pdf, seq_05_all, figlist
create_pdf = True
seq_intersections_list = []
figlist = []
r_a1comb1_list = []
p_a1comb1_list = []
r_b1comb1_list = []
p_b1comb1_list = []
    

correlations()
'''function that looks for correlation between gender, age and nLymphocytes, and data parameters of 0.5ug samples'''

# unique_rearrang_bar()
# template_number_bar()


for sample in sample_list:
    print sample
# #list of sample analysis function:
    generate_SeqLists(sample)  # required for following functions    
#     gen_venns(sample)                  
#     seq_intersections_stats(sample)
#     generate_SeqDF(sample)
    seq_temp_calc(sample)  
#     generate_replicate_view(sample) ## generate an image for each sample that contains: 
                                    # # (1) 1ug sample comparison using scatter
                                    # # (2) histograms of copy numbers per sequence, in each 1ug replicate
                                    # # (3) rarefaction plots for each 1ug replicate
#     
#     frequency_comparison(sample)
#     seq_freq_correl(sample)
# seq_intersections_df=pd.DataFrame(seq_intersections_list)
# seq_intersections_df.to_csv('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/seq_intersections.csv')

# seq_freq_correl_df=pd.DataFrame({'Sample': sample_list, 'R - a1-comb1' : r_a1comb1_list, 'P - a1-comb1' : p_a1comb1_list, 'R - b1-comb1' : r_b1comb1_list, 'P - b1-comb1' : p_b1comb1_list})
# seq_freq_correl_df.to_csv('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/QC/seq_freq_correl_removezeros.csv')
# if create_pdf:
#     with PdfPages('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Images/seq_freq_correl.pdf') as pdf:
#         for fig in figlist:
#             pdf.savefig(fig)
#     pdf.close


