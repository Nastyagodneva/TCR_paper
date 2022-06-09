## main module: 
## (1) gets sample data file names from the directory - only files that start with HIP and are tsv files.
## (2) generates a csv file to summarize the stats for each sample.
## (3) generates dictionaries to summarize the stats.

from os import listdir
from os.path import isfile, join
from Utils import Load, Write
import pandas as pd
import numpy as np
import math





def generate_seq_list (data_file):
##generates nt and aa seq lists from the #1 and #2 columns of the .tsv file, respectively. cleans these list by removing the header and empty values
        global AAseqs
        lines=data_file.readlines()
        DNAseqs=[]
        AAseqs=[]
        for x in lines:
                DNAseqs.append(x.split("\t") [0])
                AAseqs.append(x.split("\t") [1])
        AAseqs=AAseqs[1::]
        AAseqs = [x for x in AAseqs if x !='']
        DNAseqs=DNAseqs[1::]
        DNAseqs = [x for x in DNAseqs if x !='']
        
        
def calculate_sample_stat(AAseqs):
##calculates the number of unique aa sequences and the maximal number of DNA seq per aa seq, for each file. 
        from collections import Counter
        global max_values, n_uniq_aa
        n_uniq_aa=len(set(AAseqs))
        aa_counter=Counter(AAseqs)
        aa_counter_values=aa_counter.values()
        max_values=max(aa_counter_values)
       
         
    
def get_sample_tags():
    global sample_tags
    ST=pd.read_csv("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SampleTags.csv")
    sample_tags=ST.set_index('Sample')
   
    

def get_age_gender(sample_tags, datafile):
    global sample_age,sample_gender, sample_age_group
    if datafile in sample_tags.index:
        sample_age=float(sample_tags.loc[datafile,'Age'])       
        sample_gender=sample_tags.loc[datafile,'Gender']     
        if sample_gender==' Female':
            sample_gender='Female'
        if sample_gender==' Male':
            sample_gender='Male'
    else:
        sample_age=np.NaN
        sample_gender=np.NaN 
    
    sample_age_group=[]
    if math.isnan(sample_age):
        sample_age_group=np.NaN
    elif sample_age <=18:
        sample_age_group='adolescent'
    elif sample_age<=35:
        sample_age_group='young adult'
    elif sample_age<=60:
        sample_age_group='adult'
    else:
        sample_age_group='elderly'
    
                                


global n_max_aa_freq_all, n_uniq_aa_all, sample_tags
onlyfiles = [f for f in listdir("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/") if isfile(join("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/", f))]
onlyfiles = [datafile for datafile in onlyfiles if datafile.startswith ('HIP') and datafile.endswith('.tsv')]

n_uniq_aa_all = {}
n_max_aa_freq_all={}

# import csv
#csvfile = open('/net/mraid08/export/genie/Microbiome/Analyses/ShaniBAF/TCR_demo_data/SampleStat', 'w')
# fieldnames = ['Sample', 'N unique aa seq','N max nt seq per aa seq']
# writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
# writer.writeheader()
countfiles=0
dp1=[] #generates the dataframe
get_sample_tags()



for datafile in onlyfiles:
        global n_uniq_aa, max_value, sample_age, sample_gender
        countfiles=countfiles+1
        data_file = open("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/%s" % datafile, 'r')
        print datafile
        print countfiles
        AAseqs=[]
        generate_seq_list (data_file)
        calculate_sample_stat(AAseqs)
        get_age_gender(sample_tags,datafile)
        
        n_uniq_aa_all[datafile]=n_uniq_aa
        n_max_aa_freq_all[datafile]=max_values
        
#         writer.writerow({'Sample': datafile, 'N unique aa seq': n_uniq_aa, 'N max nt seq per aa seq': max_values})
        dp1.append({'Sample': datafile, 'Age': sample_age, 'Age Group': sample_age_group,'Gender': sample_gender,'N unique aa seq': n_uniq_aa, 'N max nt seq per aa seq': max_values})
# csvfile.close
stat_summary=pd.DataFrame(dp1)
stat_summary=stat_summary.set_index('Sample')
stat_summary_file=Write('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Stat_summary.df',stat_summary)
stat_summary.to_csv('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/SampleStat.csv')
print stat_summary



