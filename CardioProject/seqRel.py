'''
this is the most updated version of the sequence reliability checks!! based on thenotebook 'Sequence reliability'
'''



import os
#from os import listdir
#from os.path import isfile, join
# from Utils import Load, Write
import pandas as pd
import numpy as np
#from scipy import stats
# import math
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib import cm
# import plotly.plotly as py
from myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot
#from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
#from Bio.SeqUtils import GC
from pop_calcs import *
from pop_organize import *
from queue.qp import qp,fakeqp
from addloglevels import sethandlers
import logging 
from Utils import cacheOnDisk
from bm_preproc import BoyerMoore

#------------------------------------------------------------------------------
## functions required:
"""bm_preproc.py: Boyer-Moore preprocessing."""

__author__ = "Ben Langmead"

import unittest


def z_array(s):
    """ Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s """
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)

    # Initial comparison of s[1:] with prefix
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else:
            break

    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1

    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            # Case 2
            # Calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                # Case 2a: zkp wins
                z[k] = zkp
            else:
                # Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    """ Compile the N array (Gusfield theorem 2.2.2) from the Z array """
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    """ Compile L' array (Gusfield theorem 2.2.2) using p and N array.
        L'[i] = largest index j less than n such that N[j] = |P[i:]| """
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    """ Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]| """
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):
    """ Compile lp' array (Gusfield theorem 2.2.4) using N array. """
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    """ Return tables needed to apply good suffix rule. """
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    """ Given a mismatch at offset i, and given L/L' and l' arrays,
        return amount to shift as determined by good suffix rule. """
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    """ Given a full match of P to T, return amount to shift as
        determined by good suffix rule. """
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    """ Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table.  Table is indexed by offset
        then by character. """
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab


class BoyerMoore(object):
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    def __init__(self, p, alphabet='ACGT'):
        # Create map from alphabet characters to integers
        self.amap = {alphabet[i]: i for i in range(len(alphabet))}
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        """ Return # skips given by bad character rule at offset i """
        assert c in self.amap
        assert i < len(self.bad_char)
        ci = self.amap[c]
        return i - (self.bad_char[i][ci]-1)

    def good_suffix_rule(self, i):
        """ Given a mismatch at offset i, return amount to shift
            as determined by (weak) good suffix rule. """
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        """ Return amount to shift in case where P matches T """
        return len(self.small_l_prime) - self.small_l_prime[1]


class TestBoyerMoorePreproc(unittest.TestCase):

    def test_z_1(self):
        s = 'abb'
        #    -00
        z = z_array(s)
        self.assertEqual([3, 0, 0], z)

    def test_z_2(self):
        s = 'abababab'
        #    00604020
        z = z_array(s)
        self.assertEqual([8, 0, 6, 0, 4, 0, 2, 0], z)

    def test_z_3(self):
        s = 'abababab'
        #    00604020
        z = z_array(s)
        self.assertEqual([8, 0, 6, 0, 4, 0, 2, 0], z)

    def test_n_1(self):
        s = 'abb'
        #    01-
        n = n_array(s)
        self.assertEqual([0, 1, 3], n)

    def test_n_2(self):
        s = 'abracadabra'
        #    1004010100-
        n = n_array(s)
        self.assertEqual([1, 0, 0, 4, 0, 1, 0, 1, 0, 0, 11], n)

    def test_n_3(self):
        s = 'abababab'
        #    0204060-
        n = n_array(s)
        self.assertEqual([0, 2, 0, 4, 0, 6, 0, 8], n)

    def test_big_l_prime_1(self):
        s = 'abb'
        #    001
        big_l_prime = big_l_prime_array(s, n_array(s))
        self.assertEqual([0, 0, 2], big_l_prime)

    def test_big_l_prime_2(self):
        s = 'abracadabra'
        #    01234567890
        # L' 00000003007
        # L  00000003337
        big_l_prime = big_l_prime_array(s, n_array(s))
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8], big_l_prime)

    def test_small_l_prime_1(self):
        s = 'abracadabra'
        # N  1004010100-
        # l'           1
        # l'        4
        # l' 44444444111
        small_l_prime = small_l_prime_array(n_array(s))
        self.assertEqual([11, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1], small_l_prime)

    def test_good_suffix_match_mismatch_1(self):
        p = 'GGTAGGT'
        big_l_prime, big_l, small_l_prime = good_suffix_table(p)
        self.assertEqual([0, 0, 0, 0, 3, 0, 0], big_l_prime)
        self.assertEqual([0, 0, 0, 0, 3, 3, 3], big_l)
        self.assertEqual([7, 3, 3, 3, 3, 0, 0], small_l_prime)
        self.assertEqual(0, good_suffix_mismatch(6, big_l_prime, small_l_prime))
        self.assertEqual(0, good_suffix_mismatch(6, big_l, small_l_prime))
        #  t:      xT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(7, good_suffix_mismatch(5, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(5, big_l, small_l_prime))
        #  t:     xGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(7, good_suffix_mismatch(4, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(4, big_l, small_l_prime))
        #  t:    xGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(3, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(3, big_l, small_l_prime))
        #  t:   xAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(2, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(2, big_l, small_l_prime))
        #  t:  xTAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(1, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(1, big_l, small_l_prime))
        #  t: xGTAGGT
        #  p: GGTAGGT
        # L': -000300
        #  L: -000333
        self.assertEqual(4, good_suffix_mismatch(0, big_l_prime, small_l_prime))
        self.assertEqual(4, good_suffix_mismatch(0, big_l, small_l_prime))

    def test_good_suffix_table_1(self):
        s = 'abb'
        #    001
        big_l_prime, big_l, small_l_prime = good_suffix_table(s)
        self.assertEqual([0, 0, 2], big_l_prime)
        self.assertEqual([0, 0, 2], big_l)
        self.assertEqual([3, 0, 0], small_l_prime)

    def test_good_suffix_table_2(self):
        s = 'abracadabra'
        #    01234567890
        # L' 00000003007
        # L  00000003337
        # l' -4444444111
        big_l_prime, big_l, small_l_prime = good_suffix_table(s)
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8], big_l_prime)
        self.assertEqual([0, 0, 0, 0, 0, 0, 0, 4, 4, 4, 8], big_l)
        self.assertEqual([11, 4, 4, 4, 4, 4, 4, 4, 1, 1, 1], small_l_prime)

#if __name__ == '__main__':
#    unittest.main()



#---------
def boyer_moore_with_mm(p, p_bm, t,n):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p, n=number of mismatches allowed"""
        
    i = 0
    
    #occurrences = []
    mm_count_list=[]
    while i < len(t) - len(p) + 1:
        mm_count=0
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                mm_count+=1
                mismatched = True
                if mm_count==1:               
                    skip_bc = p_bm.bad_character_rule(j, t[i+j])
                    skip_gs = p_bm.good_suffix_rule(j)
                    shift = max(shift, skip_bc, skip_gs)
                if mm_count>n:
                    break
                            
        if not mismatched:
                #occurrences.append(i)
                skip_gs = p_bm.match_skip()
                shift = max(shift, skip_gs)
        i += shift
        mm_count_list.append(mm_count)
    min_mm=min(mm_count_list)
    return mm_count_list, min_mm

#-------------------------------------------------
#sample_list=['BD207','BD207a_05ug','BD207b_05ug','BD418','BD415a_05ug','BD415b_05ug','BD418_322939_2']
sample_list=['HIP13853','HIP02848','HIP14036']

for sample_name in sample_list:
    basePath='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/seqRel/check_seqRel_%s' %sample_name ## define the path to
    if not os.path.exists(basePath):
        print 'not exist' # generate the directory of doesn't exist
        os.makedirs(basePath)
        
    @cacheOnDisk(basePath=basePath, filename='%(sample_name)s_seqRel_%(min_seq)s_%(max_seq)s', force=True)
    def check_seqRel(sample_name,min_seq,max_seq):
        n=3
        print 'getting sample df...'
        sample_df=pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/%s.tsv" %sample_name) 
        if 'count (templates)' in sample_df.columns.values:
            column='count (templates)'
        elif 'count (templates/reads)' in sample_df.columns.values:
            column='count (templates/reads)'
        else:
            column='count (reads)'
        columns_to_keep = ['nucleotide', column]
        seq_rel_df = sample_df[columns_to_keep]
        #seq_rel_df.set_index('nucleotide', inplace=True)
        print 'finished generating sample dataframes...'
        comparison_list=[]
        if max_seq>len(seq_rel_df):
            max_seq=len(seq_rel_df)
        for seq1 in range(min_seq, max_seq):
#             if seq1%500==0:
#                 print ('count_seq1=%s' %seq1)
            p=seq_rel_df.loc[seq1,'nucleotide']
            #segment_length=round(len(p)/(n+1))
            #min_mis=1000
            #if segment_length>=2:  
            #    for i in range(n+1): # loop over partitions
            #        start=int(i*segment_length) #define start and end site for each partition
            #        end=int(min((i+1)*segment_length,len(p)))
            p_bm=BoyerMoore(p, alphabet='ACGT')
            count_seq2=0
            for seq2 in range(seq1+1,len(seq_rel_df.index)):
                #if count_seq2%50000==0:
                #    print ('count_seq2=%s' %count_seq2)
                t=seq_rel_df.loc[seq2,'nucleotide']
                mm_count_list, min_mm=boyer_moore_with_mm(p, p_bm, t,n)
                
                count_seq2=count_seq2+1
                if min_mm<=n:
                    seq1_reads=seq_rel_df.loc[seq1,column]
                    seq2_reads=seq_rel_df.loc[seq2,column]
    
                    comparison_list.append({'#mm':min_mm,'#seq1':seq1,'seq1':p,'seq1 reads':seq1_reads, '#seq2':seq2,'seq2':t,'seq2 reads':seq2_reads})
        
        seq_compar_df=pd.DataFrame(comparison_list)
        return seq_compar_df





    sample_df=pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/%s.tsv" %sample_name)
    sampleLength=len(sample_df)
    print 'the length of %s is %s' %(sample_name, sampleLength)
       
    sethandlers()
    os.chdir(basePath)
    ## the inputs for qp are: jobname, q=machine list, *** add max_r to prevent exploding the cluster!!***
    with qp('check_seqRel_job3', q=['himem7.q', 'himemint.q'], mem_def="10G", trds_def=2, deleteCSHwithnoerr=True, tryrerun=False, max_u=120) as q:
        q.startpermanentrun()
        wait_on =[]
        
    ##now define a loop that divide the job and send each part seperately:
    ## consider making a list of integer for the min/max_seq so the first jobs will be shorter than the last ones. 
        min_seq=0
        max_seq=1000
        sample_name=sample_name
        while  min_seq<sampleLength:                                     
            print min_seq
            wait_on.append(q.method(check_seqRel,kwargs={'sample_name':sample_name,'min_seq':min_seq,'max_seq':max_seq}))
                ##q.method takes the desired function with its arguments and send it to the queue.
            min_seq=min_seq+1000
            max_seq=max_seq+1000
        q.wait(wait_on)











    