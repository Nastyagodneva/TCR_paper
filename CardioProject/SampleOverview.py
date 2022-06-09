'''
this script generates bar graphs from the SampleOverview file (exported from the analyzer)
it can generate either a graph comparing the number of unique rearrangements, or total templates
depending on the 'parameter' and 'name' variables. 
make sure to change the file name to avoid overwrite, if necessary

***THIS script generates the same graphs as the 'unique_rearrang_bar()' and the 'template_number_bar()' in
OptSeqAnalysis module generate***
'''



import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import math
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import collections




overview_df = pd.read_table("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/QC/SampleOverview_new.tsv")  
overview_df = overview_df.set_index('sample_name')
rearrangements = overview_df['total_rearrangements'].to_dict()
templates= overview_df['total_templates'].to_dict()
# fract_prod = overview_df['total_rearrangements'].to_dict()
values_list = []
# print rearrang
# print rearrang.keys()
# print rearrang.keys()[0]
sample_list = ['139', '207', '415', '418', '148', '202', '398', '419']

parameter=rearrangements  ##choose 'rearrangements/templates'
name='Rearrangements'  ##choose Rearrangements/Templates

for sample in sample_list:
    print sample
    sample_dict = {}
    for i in range(len(parameter.keys())):
        if sample in parameter.keys()[i]:
            if 'a_1ug' in parameter.keys()[i]:
                a1 = parameter.values()[i]
                sample_dict['a1'] = a1
            elif 'b_1ug' in parameter.keys()[i]:
                b1 = parameter.values()[i]
                sample_dict['b1'] = b1
            elif 'a_05ug' in parameter.keys()[i]:
                a05 = parameter.values()[i]
                sample_dict['a05'] = a05
            elif 'b_05ug' in parameter.keys()[i]:
                b05 = parameter.values()[i]
                sample_dict['b05'] = b05
            elif '05ug' in parameter.keys()[i] and 'a_05ug' not in parameter.keys()[i] and 'b_05ug' not in parameter.keys()[i]:
                comb05 = parameter.values()[i]
                sample_dict['comb05'] = comb05
            elif '1ug' in parameter.keys()[i] and 'a_1ug' not in parameter.keys()[i] and 'b_1ug' not in parameter.keys()[i]:
                comb1 = parameter.values()[i]
                sample_dict['comb1'] = comb1
                print comb1
    ordered_sample_dict = collections.OrderedDict(sorted(sample_dict.items()))
    values_list.append(ordered_sample_dict)
list_a05 = []
list_b05 = []
list_a1 = []
list_b1 = []
list_comb05=[]
list_comb1=[]
list_a05b05=[]
list_a1b1=[]
list_erron_05=[]
list_erron_1=[]
for i in range(len(values_list)):
    sampledict = values_list[i]
    a_05 = sampledict['a05']
    list_a05.append(a_05)
    b_05 = sampledict['b05']
    list_b05.append(b_05)
    comb_05 = sampledict['comb05']
    list_comb05.append(comb_05)
    list_a05b05.append(a_05+b_05)
    if 'a1' in sampledict:
        a_1 = sampledict['a1']
        list_a1.append(a_1)
    else:
        list_a1.append(0.0)
    if 'b1' in sampledict:    
        b_1 = sampledict['b1']
        list_b1.append(b_1)
    else:
        list_b1.append(0.0)
    if 'comb1' in sampledict:    
        comb_1 = sampledict['comb1']
        print comb_1
        list_comb1.append(comb_1)
    else:
        list_comb1.append(0.0)
    if 'a1' in sampledict and 'b1' in sampledict:    
        list_a1b1.append(a_1+b_1)
    else:
        list_a1b1.append(0.0)

erron_temp_df=pd.DataFrame({'sample': sample_list, 'list_a05b05': list_a05b05, 'list_comb05': list_comb05, 'list_a1b1': list_a1b1, 'list_comb1': list_comb1 })
erron_temp_df['perc_erron_05']=(erron_temp_df['list_a05b05']-erron_temp_df['list_comb05'])/erron_temp_df['list_a05b05']*100
erron_temp_df['perc_erron_1']=(erron_temp_df['list_a1b1']-erron_temp_df['list_comb1'])/erron_temp_df['list_a1b1']*100
print erron_temp_df

erron_temp_df.to_csv(('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/QC/percent_erron_%s.csv' %name))


# list05 = list_a05 + list_b05
# print list05
# list1 = list_a1 + list_b1
# print list1
# while 0.0 in list1: list1.remove(0.0)
# print list1
# 
# mean_05 = np.mean(list05)
# mean_1 = np.mean(list1)
# MU_s, MU_p = mannwhitneyu(list05, list1, use_continuity=True, alternative='two-sided')
# print mean_05, mean_1
# print MU_p


  
fig=plt.figure(figsize=(8,8))    
space = 0.2
width = (1 - space) / 8
indeces = range(len(sample_list))
pos_a05 = [i - 3.5 * width for i in indeces]
pos_b05 = [i - 2.5 * width for i in indeces]
pos_comb05 = [i - 1.5 * width for i in indeces]
pos_a05b05 = [i - 0.5 * width for i in indeces]
pos_a1= [i + 0.5 * width for i in indeces]
pos_b1 = [i + 1.5 * width for i in indeces]
pos_comb1 = [i + 2.5 * width for i in indeces]
pos_a1b1 = [i + 3.5 * width for i in indeces]

print list_comb1


  
plt.bar(pos_a05 , list_a05, align='center', label='a_0.5ug', width=width, color='blue')
plt.bar(pos_b05 , list_b05, align='center', label='b_0.5ug', width=width, color='blue')
plt.bar(pos_comb05 , list_comb05, align='center', label='comb_0.5ug', width=width, color='green')
plt.bar(pos_a05b05 , list_a05b05, align='center', label='a05b05', width=width, color='black') 
plt.bar(pos_a1 , list_a1, align='center', label='a_0.5ug', width=width, color='red')
plt.bar(pos_b1 , list_b1, align='center', label='b_0.5ug', width=width, color='red')
plt.bar(pos_comb1 , list_comb1, align='center', label='comb_1ug', width=width, color='yellow')
plt.bar(pos_a1b1 , list_a1b1, align='center', label='a1b1', width=width, color='pink') 
plt.xticks(indeces, sample_list, fontsize=9)
# plt.yscale('log')
# plt.ylim(100000,10000000)
# plt.yticks(np.logspace(5,7))
plt.title('Number of %s Per sample and DNA amount' %name)
plt.ylabel('n %s' %name)
plt.xlabel('Sample')
plt.margins(0.02, 0)
plt.legend(loc='upper right', fontsize=12)
# plt.annotate('Test', xy=(1, 0), xycoords='axes fraction', fontsize=16,
#                 horizontalalignment='right', verticalalignment='bottom')
# plt.text(0.5, 0.5, 'Mean-0.5ug samples= %s' % mean_05, ha='center', va='center',transform=transAxes, fontsize=14)
# plt.text(0.91, 0.91 'Mean-1ug samples=%s' % mean_1, fontsize=14)

plt.show()
with PdfPages('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Images/SamplesOverView_%s.pdf' %name) as pdf:
    pdf.savefig(fig)
    pdf.close


# 
# 
# fract_prod = overview_df['fraction_productive'].to_dict()
# values_list = []
# 
# for sample in sample_list:
#     print sample
#     sample_dict = {}
#     for i in range(len(fract_prod.keys())):
#         if sample in fract_prod.keys()[i]:
#             if 'a_1' in fract_prod.keys()[i]:
#                 a1 = fract_prod.values()[i]
#                 sample_dict['a1'] = a1
#             elif 'b_1' in fract_prod.keys()[i]:
#                 b1 = fract_prod.values()[i]
#                 sample_dict['b1'] = b1
#             elif 'a_05' in fract_prod.keys()[i]:
#                 a05 = fract_prod.values()[i]
#                 sample_dict['a05'] = a05
#             elif 'b_05' in fract_prod.keys()[i]:
#                 b05 = fract_prod.values()[i]
#                 sample_dict['b05'] = b05
#     ordered_sample_dict = collections.OrderedDict(sorted(sample_dict.items()))
#     values_list.append(ordered_sample_dict)
# 
# list_a05 = []
# list_b05 = []
# list_a1 = []
# list_b1 = []
# for i in range(len(values_list)):
#     sampledict = values_list[i]
#     a_05 = sampledict['a05']
#     list_a05.append(a_05)
#     b_05 = sampledict['b05']
#     list_b05.append(b_05)
#     if 'a1' in sampledict:
#         a_1 = sampledict['a1']
#         list_a1.append(a_1)
#     else:
#         list_a1.append(0.0)
#     if 'b1' in sampledict:    
#         b_1 = sampledict['b1']
#         list_b1.append(b_1)
#     else:
#         list_b1.append(0.0)
#         
# list05 = list_a05 + list_b05
# print list05
# list1 = list_a1 + list_b1
# print list1
# while 0.0 in list1: list1.remove(0.0)
# print list1
# mean_05 = np.mean(list05)
# mean_1 = np.mean(list1)
# MU_s, MU_p = mannwhitneyu(list05, list1, use_continuity=True, alternative='two-sided')
# print mean_05, mean_1
# print MU_p
# 
# 
# 
# plt.subplot(2, 1, 2)
#    
#     
# space = 0.4
# width = (1 - space) / 4
# indeces = range(len(sample_list))
# pos_a05 = [i - 1.5 * width for i in indeces]
# pos_b05 = [i - 0.5 * width for i in indeces]
# pos_a1 = [i + 0.5 * width for i in indeces]
# pos_b1 = [i + 1.5 * width for i in indeces]
#   
# plt.bar(pos_a05 , list_a05, align='center', label='a_0.5ug', width=width, color='blue')
# plt.bar(pos_b05 , list_b05, align='center', label='b_0.5ug', width=width, color='red')
# plt.bar(pos_a1 , list_a1, align='center', label='a_1ug', width=width, color='green')
# plt.bar(pos_b1 , list_b1, align='center', label='b_1ug', width=width, color='black')   
# plt.xticks(indeces, sample_list, fontsize=9)
# # plt.yscale('log')
# plt.ylim(0.7, 1)
# # plt.yticks(np.logspace(5,7))
# plt.title('Productive rearrangement fraction Per sample and DNA amount')
# plt.ylabel('Fraction Productive')
# plt.xlabel('Sample')
# plt.margins(0.02, 0)
# plt.legend(loc='upper right', fontsize=12)
# 
# plt.subplots_adjust(left=0.08, bottom=0.10, right=0.95, top=0.89, wspace=0.38, hspace=0.25)
# fig = plt.gcf()
# plt.show()
# 
# with PdfPages('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/QC/SamplesOverView.pdf') as pdf:
#     pdf.savefig(fig)
#     pdf.close
    

