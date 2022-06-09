## this script is not in used and not updated!


from os import listdir
from os.path import isfile, join
from Utils import Load, Write
import pandas as pd
import numpy as np
import math
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm
import collections
# import plotly.plotly as py


  
def get_data(): 
##sample types:
    sample_type=['Rep1_0.5', 'Rep2_0.5', 'Comb_0.5', 'Rep1_1', 'Rep2_1', 'Comb_1']
##insert file names for 
    file_list=['HIP00110', 'HIP00169', 'HIP00594', 'HIP00602', 'HIP00614', 'HIP00640']
    global df_dict
    df_dict={}
    for i in range(len(sample_type)):
        sample=sample_type[i]
        print i
        print sample
        df_dict['%s_All' % sample]=pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/%s.tsv" %file_list[i])
        df_dict['%s_Prod' % sample]=df_dict['%s_All' % sample][df_dict['%s_All' % sample]['sequenceStatus']=='In']
        df_dict['%s_Non_Prod' % sample]=df_dict['%s_All' % sample][df_dict['%s_All' % sample]['sequenceStatus']!='In']
    print len(df_dict)
    print df_dict.keys()
   
   
def reads_graphs():
    df_list=df_dict.keys()
    print df_list
    reads_all={}
    reads_prod={}
    reads_non_prod={}
    dict_seq_status={}
    for i in df_list:
        if 'All' in i:
            reads_all[i]=df_dict[i]['count (reads)'].sum()
     
            sample_seq_status = df_dict[i].groupby(['sequenceStatus'])[['frequencyCount (%)']].sum()
            sample_seq_status_list=list(sample_seq_status['frequencyCount (%)'])
            sample_seq_status_list=[round(j,2) for j in sample_seq_status_list]
            dict_seq_status[i]=sample_seq_status_list
#             labels_seq_status=[str(j)+'%' for j in sample_seq_status_list]  
#             legend=sample_seq_status_list.index.tolist()
        elif 'Non_Prod' in i:
            reads_non_prod[i]=df_dict[i]['count (reads)'].sum()
        else:
            reads_prod[i]=df_dict[i]['count (reads)'].sum()    
    
    ordered_reads_all=collections.OrderedDict(sorted(reads_all.items()))
    print ordered_reads_all
    ordered_reads_prod=collections.OrderedDict(sorted(reads_prod.items()))
    print ordered_reads_prod
    ordered_reads_non_prod=collections.OrderedDict(sorted(reads_non_prod.items()))
    print ordered_reads_non_prod
    ordered_dict_seq_status=collections.OrderedDict(sorted(dict_seq_status.items()))
    print ordered_dict_seq_status
    
    plt.figure(figsize=(8, 8))
    
    plt.subplot(3,1,1)
    space=0.4
    width=(1-space)/3
    indeces = range(len(ordered_reads_prod))
    pos_all = [i-width for i in indeces]
    pos_prod = [i for i in indeces]
    pos_non_prod = [i+width for i in indeces]
    plt.bar(pos_all , ordered_reads_all.values(), align='center', label='All', width=width)
    plt.bar(pos_prod, ordered_reads_prod.values(), align='center', label='Prod', color='red', width=width)
    plt.bar(pos_non_prod, ordered_reads_non_prod.values(), align='center', label='Non-Prod', color='green', width=width)
    plt.xticks(pos_prod, ordered_reads_prod.keys(), fontsize=9)
    plt.yscale('log')
    plt.ylim(100000,10000000)
    plt.yticks(np.logspace(5,7))
    plt.title('Total reads')
    plt.ylabel('Number of reads')
    plt.margins(0.02,0.02)
    plt.legend(loc='upper right', fontsize=9)
    
    
    for i in range(0,len(ordered_dict_seq_status)):
        print i
        print ordered_dict_seq_status.values()[i]
        labels=list(str(j)+'%' for j in ordered_dict_seq_status.values()[i])
        print labels
        
        plt.subplot(3,3,i+4)
        plt.pie(ordered_dict_seq_status.values()[i], shadow=True, startangle=90,labels=labels)
        plt.title(ordered_dict_seq_status.keys()[i])
        if i==0:
            plt.legend(['In','Out','Stop'], loc='best',fontsize=8)
        plt.subplots_adjust(left=0.09,bottom=0.10, right=0.96, top=0.95, wspace=0.20,hspace=0.35)
    
#     
#     if create_pdf:
#         global fig1
#         fig1=plt.gcf()
#         figlist.append(fig1)
#     else:
    plt.show()
    
    
      
#     plt.subplot(2,3,2)
#     plt.bar(range(len(reads_prod)), reads_prod.values(), align='center')
#     plt.xticks(range(len(reads_prod)), reads_prod.keys(), rotation='vertical')
#     plt.yscale('log')
#     plt.ylim(0,10000000)
#     plt.yticks(np.logspace(0,7,7))
#     plt.title('Total reads')
#     plt.ylabel('Number of reads')
#     plt.margins(0.02,0.02)
#       
#     plt.subplot(2,3,3)
#     plt.bar(range(len(reads_non_prod)), reads_non_prod.values(), align='center')
#     plt.xticks(range(len(reads_non_prod)), reads_non_prod.keys(), rotation='vertical')
#     plt.yscale('log')
#     plt.ylim(0,10000000)
#     plt.yticks(np.logspace(0,7,7))
#     plt.title('Total reads')
#     plt.ylabel('Number of reads')
#     plt.margins(0.02,0.02)
#     plt.show()
    
       
    



    

get_data()
reads_graphs()     

