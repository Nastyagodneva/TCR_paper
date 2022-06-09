'''
this file contains 'organizational' functions for the 'PopulationAnalysis.py' module:
getting data frim TSV file, etc.  

'''


from os import listdir
from os.path import isfile, join
# from Utils import Load, Write
import pandas as pd
import numpy as np
from scipy import stats
# import math
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib import cm
# import plotly.plotly as py
from myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot
from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
from Bio.SeqUtils import GC

#---------------------------------------------------------------------------------------

def get_sample_data(sample_name, generate_dfs): 
## this function generates dfs (general, only productive and only non productive) for each sample, and save
## as pickles. alternatively, it loads these dfs from pickles. 
    
    print 'getting sample data...'
    
    if generate_dfs:
        sample_df=pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/%s.tsv" %sample_name)  
        sample_df_prod = sample_df[sample_df['sequenceStatus'] == 'In']
        sample_df_non_prod = sample_df[sample_df['sequenceStatus'] != 'In']
        
        pickle.dump(sample_df, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_%s' %sample_name, "wb"))
        pickle.dump(sample_df_prod, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_prod_%s' %sample_name, "wb"))
        pickle.dump(sample_df_non_prod, open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_non_prod_%s' %sample_name, "wb"))
    else:
        sample_df= pickle.load( open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_%s' %sample_name,"rb" ))
        sample_df_prod= pickle.load( open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_prod_%s' %sample_name,"rb" ))
        sample_df_non_prod= pickle.load( open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/sample_df_non_prod_%s' %sample_name,"rb" ))
    
    print 'finished getting sample data'
    return sample_df, sample_df_prod, sample_df_non_prod

#____________________________________________________________________________________
## the following function generate a dict to convert the column number to the number of df and
## function used to generate the data in this column:
## the format of the items in the dictionary are {column: [df number, function number]

def generate_col_dict(df_n, function_list):
    ind=0
    col_dict={}
    for i in range(df_n):
        for j in range(len(function_list)):
            col_dict[ind]=[j, i]
            ind+=1
    return col_dict

#--------------------------------------------------------------------------------

def generate_col_title(col_dict, df_names, function_names):
## call this function with lists of df and function names (Strings) that match exactly 
## the df and function lists used to generate the result DataFrame:

    title_dict={}
    for k,v in col_dict:
        d=v[0]
        f=v[1]
        title_dict[k]=[df_names[d]+'_'+function_names[f]]
    return title_dict

#--------------------------------------------------------------------------------------


## this function checks whether a sample has a dfs, so the calculating functions will
## not get stuck trying to generate df from a "problematic" sample. 
def get_sample_with_dfs():
    print 'getting list of samples that have dfs...'
    df_file_names = [f for f in listdir("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles") if isfile(join("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles", f))]
    df_file_names = [df_file for df_file in df_file_names if df_file.startswith ('sample_df_HIP')]
    samples_with_df=[s.replace("sample_df_", "") for s in df_file_names]
    print 'finished getting list of samples that have dfs'
    return df_file_names,samples_with_df

#--------------------------------------------------------------------------------------