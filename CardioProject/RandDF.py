'''
use this script to try new stuff
'''


import numpy as np
# import scipy as sc
import pandas as pd
# import math
# from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import random, string
import cPickle as pickle
from myplots import roundup, rounddown
import sys
import timeit
from Bio.Seq import Seq
from myplots import percentile_cut_off





# def find_decimal_fold(num):
#     s_num=str(num)
#     fraction=s_num.split('.')[1]
#     print fraction
#     fold=0
#     for i in range(len(fraction)):
#         if fraction[i]!='0':
#             fold=i+1
#             break
#     return fold
#     
# find_decimal_fold(1.00000026)
#         



# my_seq='AGGGGTTTGCAA'
# # print('reverse complement is %s' % my_seq.reverse())
# # print('reverse complement is %s' % complement(my_seq.reverse()))
# print('reverse complement is %s' % my_seq.reverse_complement())
# # print('reverse complement is %s' % reverse(my_seq.complement()))      


# my_seq = Seq("TGGGCCTCATATTTATCCTATATACCATGTTCGTATGGTGGCGCGATGTTCTACGTGAATCCACGTTCGAAGGACATCATACCAAAGTCGTACAATTAGGACCTCGATATGGTTTTATTCTGTTTATCGTATCGGAGGTTATGTTCTTTTTTGCTCTTTTTCGGGCTTCTTCTCATTCTTCTTTGGCACCTACGGTAGAG")
# my_aa = my_seq.translate()
# print my_aa

 






# def get_extension1(filename):
#     return(filename.split(".")[-1])
# 
# def get_extension2(filename):
#     import os.path
#     return(os.path.splitext(filename)[1])
# 
# def get_extension3(filename):
#     return filename[filename.rfind('.'):][1:]
# 
# # print get_extension1('myfile')
# print get_extension2('myfile')
# # print get_extension3('myfile')

# f=open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/aa','w')
# def trial1():
#     a={'a': 1, 'b': 2, 'c': 3}
#     for letter, number in a.items():
#         print(letter, number)
# 
# def trial2():
#     print 'bbbb'
# 
# trial1()
#  
# def count1(dna, base):
#     i = 0
#     for c in dna:
#         if c == base:
#             i += 1 
#     return i
# 
# def count2(dna, base):
#     i = 0 
#     for j in range(len(dna)):
#         if dna[j] == base:
#             i += 1 
#     return i 
# 
# def count3(dna, base):
#     match = [c == base for c in dna]
#     return sum(match)
# 
# def count4(dna, base):
#     return dna.count(base)
# 
# def count5(dna, base):
#     return len([i for i in range(len(dna)) if dna[i] == base])
# 
# def count6(dna,base):
#     return sum(c == base for c in dna)
# dna='aaaaattttttggggggatgagdtagatgccc'
# base='t'
# 
# print timeit.timeit('count1(dna,base)',number=1000000)
#
#  
# # data=pd.DataFrame(np.random.random(size=(100,4)),columns=list('ABCD'))
# # A=list(data['A'])
# # print np.mean(A)
# # pickle.dump(A, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/A.p',"wb" ))
# # new_a= pickle.load( open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/A.p', "rb" ) )
# # print np

# print(__version__)
# print(version)
# 
# print sys.version
# print (sys.__version__)

# 