import cPickle as pickle
import pandas as pd
import numpy as np
from time import gmtime, strftime
from os import makedirs
from os.path import isdir
from ShaniBA.CardioProject.Figures.GeneralFigureFunctions import remove_spines, edit_species_names


import random
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns


################################################################################
# variable to change:
FDR_t_species=0.2
top100_phen_t_FDR = 0.1
include_unknown_species=True

outcome_name = 'CV hospitalization including chest pain'
CALC_RESULT_DIR_NAME = 'CV hospitalization including chest \
pain_min_shared_species01_max_shared_species1_min_shared_cluster03_max_shared_cluster1_3/'
excel_files_dir_name = 'excel_files_new_cv_hospitalization3/'
# sample_list1_file = 'PNP530_Age4570_Gender_Male11_HbA1C4464_BMI040.pkl'
# sample_list2_file = 'Cardio126_Age3566_Gender_Male11_HbA1C4464_BMI040.pkl'

##################################################################################

#directories:
GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
SAMPLE_LIST_DIR = GENIE_BASE_DIR + 'Sample files/BD lists/'
CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5_updated_phens/'
CALC_RESULT_DIR = FIGURE_DIR + CALC_RESULT_DIR_NAME
EXCEL_FILES_DIR = FIGURE_DIR + excel_files_dir_name
CLUSTER_ANALYSIS_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/sharingAnalysis/seqClusters_allProd_maxdist1/'

if not isdir(FIGURE_DIR): makedirs(FIGURE_DIR)
if not isdir(EXCEL_FILES_DIR): makedirs(EXCEL_FILES_DIR)
PRED_RESULTS_DIR='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/'

if include_unknown_species: suf = '_inc_unknown_species'
else: suf=''

# general functions:

def convert_logp_values(df,direction_col,logp_col='-log10p_MW',suffix=None):
    df['p_to_show'] = np.where(df[direction_col] == 1, df[logp_col],df[logp_col]* -1)
    if suffix is not None:
        df = df.rename(columns={'p_to_show': 'p_to_show_%s' %suffix})
    return df


############################################################################################
#  Generation of result excel files for figure 5
############################################################################################

# get sample lists - only samples that have swab data
with open(CALC_RESULT_DIR + 'sample_list1_swab.pkl', 'rb') as fp1:
    sample_list1 = pickle.load(fp1)
with open(CALC_RESULT_DIR + 'sample_list2_swab.pkl', 'rb') as fp2:
    sample_list2 = pickle.load(fp2)



# (1) generate summary files:
summary_df = pd.read_excel(CALC_RESULT_DIR + 'summary_withAnot.xlsx')

#filterout directionOK=0
summary_df2 = summary_df[summary_df['directionOK']==1]

#get species+families:
summary_df2['species2'] = summary_df2['species'].apply(lambda x: '|'.join(x.split('|')[-6:-3]+[x.split('|')[-1]]))
summary_df2['genus'] = summary_df2['species'].apply(lambda x: x.split('|')[-5])
summary_df2['known_mb'] = np.where(summary_df2['species2'].str.contains('g__unknown|s__unknown'),0,1)

#filter out unknown species_families:
if not include_unknown_species: summary_df2 = summary_df2[summary_df2['known_mb'] == 1]

#count how many sequences
print ('number of unique clusters: ', summary_df2['cluster'].nunique())

summary_df2.to_excel(CALC_RESULT_DIR + 'summary_processed%s.xlsx' %suf)


summary_df_by_cluster =summary_df2.groupby(['cluster','species2']).mean()
summary_df_by_cluster = summary_df_by_cluster.sort_values(by=['-log10p_MW_cluster_phen','-log10p_MW_cluster_species',
                                                              '-log10p_MW_species_phen'],ascending = False)
summary_df_by_cluster.to_excel(CALC_RESULT_DIR + 'summary_df_by_cluster%s.xlsx' %suf)

#cluster by genus and check sequence similarity
summary_df_by_genus =summary_df2.groupby(['genus','cluster']).mean()
summary_df_by_genus = summary_df_by_genus.sort_values(by=['-log10p_MW_cluster_phen','-log10p_MW_cluster_species',
                                                              '-log10p_MW_species_phen'],ascending = False)
summary_df_by_genus.to_excel(CALC_RESULT_DIR + 'summary_df_by_genus%s.xlsx' %suf)

# (2) cluster-disease data:

# get top 100 list:

cluster_phen = pd.read_excel(CALC_RESULT_DIR + 'concise_summary_df_cluster_phen.xlsx')

cluster_phen = cluster_phen.sort_values(by='-log10p_MW', ascending=False)
top100 = cluster_phen['cluster'].tolist()[:100]

# get annotation info
top100_clusters = cluster_phen[cluster_phen['cluster'].isin(top100)]
direction_col = 'Associated_with_%s_positive' %outcome_name
top100_clusters = convert_logp_values(top100_clusters,direction_col,logp_col='-log10p_MW',suffix='cluster_isCardio')

annot = summary_df[['cluster','annotation summary']].drop_duplicates()
top100_clusters_with_annot = pd.merge(top100_clusters[['cluster','p_to_show_cluster_isCardio']],annot,
                                      how='left',left_on='cluster',right_on='cluster')
annot_dict = {r',None:.+':'',r' :.+':'','(':'',')':'', 'HomoSapiens':'', 'YellowFeverVirus':'YFV', 'Type':'',
              }
for k,v in annot_dict.items():
    top100_clusters_with_annot['annotation summary'] = \
        top100_clusters_with_annot['annotation summary'].str.replace(k, v)


top100_clusters_with_annot = top100_clusters_with_annot.set_index('cluster')
top100_clusters_with_annot.to_excel(EXCEL_FILES_DIR + 'top100_clusters_with_annot.xlsx')

# (3) cluster-species data (for all associations with FDR<limit)

summary_df3 = summary_df2[(summary_df2['MW_p_corrPval_FDR0.1_cluster_species'] <= FDR_t_species) &
                          (summary_df2['MW_p_corrPval_FDR0.1_species_phen'] <= FDR_t_species)]
print ('number of unique clusters in summary_df3: ', summary_df3['cluster'].nunique())
print ('number of unique species in summary_df3: ', summary_df3['species2'].nunique())
print (summary_df3['species2'].unique())

## add minus to log p-valus if association is negative
direction_col = 'species_pos_associated_with_cluster'
summary_df3 = convert_logp_values(summary_df3,direction_col,logp_col='-log10p_MW_cluster_species',suffix='cluster_species')
col_order = (summary_df3.sort_values(by='p_to_show_cluster_species', ascending=False))['species2'].unique().tolist()
summary_df3 = summary_df3.drop_duplicates(subset=['cluster','species2'],keep='first')
summary_df3_edited = summary_df3.copy()
summary_df3_edited['species2'] = edit_species_names(summary_df3_edited['species2'].tolist())
summary_df3_edited.to_excel(EXCEL_FILES_DIR + 'cluster_species_hm_fdr%s%s_longform.xlsx')

## generate pivot table:
cluster_species_hm = summary_df3.pivot(index='cluster',columns='species2',values='p_to_show_cluster_species')
cluster_species_hm = cluster_species_hm.loc[top100,col_order]
cluster_species_hm.to_excel(EXCEL_FILES_DIR + 'cluster_species_hm_fdr%s%s.xlsx' %(str(FDR_t_species).replace('.',''),
                                                                                  suf))


# (4) cluster-phens data:
# top100_phen_df = pd.DataFrame(index=top100)
# phen_list = ['isOld','Gender_Male','HbA1C_high','Obese']
#
# for phen in phen_list:
#     dir1 = FIGURE_DIR + '%s_min_shared_species01_max_shared_species1\
# _min_shared_cluster03_max_shared_cluster1/' %phen
#     df = pd.read_excel(dir1 + 'concise_summary_df_cluster_phen.xlsx')
#     df = df[df['MW_p_corrPval_FDR0.1'] <=top100_phen_t_FDR]
#     df['logp'] = np.where(df['Associated_with_%s_positive' %phen] == 1,
#                           df['-log10p_MW'],df['-log10p_MW']*-1)
#     df = df[['cluster','logp']].set_index('cluster')
#     try:
#         df = df.loc[top100,:].rename(columns={'logp':'%s_logp' %phen})
#     except:
#         print ('none of top100 clusters is associated with %s with FDR %s' %(phen,top100_phen_t_FDR))
#         df = pd.DataFrame(index=top100,data={'%s_logp' %phen: 100*[np.nan]})
#     top100_phen_df = pd.merge(top100_phen_df,df,how='left',left_index=True,right_index=True)
# top100_phen_df.to_excel(EXCEL_FILES_DIR + 'top100_phen_df_fdr%s%s.xlsx' %(str(top100_phen_t_FDR).replace('.',''),
#                                                                           suf))

# (5) merge everything and generate one excel file
df_list = [top100_clusters_with_annot,cluster_species_hm]
cluster_species_phen_in_ACS = pd.concat(df_list,axis=1)

cluster_species_phen_in_ACS.to_excel(EXCEL_FILES_DIR + 'cluster_species_phen_in_ACS_speciesFDR%s_phenFDR_%s%s.xlsx'\
                                     %(str(FDR_t_species).replace('.',''),
                                       str(top100_phen_t_FDR).replace('.',''),
                                        suf))


# (6) additional data - save each to seperate excel file


# cluster raw data - num of clusters per sample in top100 clusters
cluster_df = pd.read_pickle(CLUSTER_ANALYSIS_DIR + \
                        'sampleByClusterDF_cohortfiltering005-09perc_dropped.dat')
cluster_df_for_fig = cluster_df.loc[sample_list1+sample_list2, top100]
cluster_df_for_fig.T.fillna(0).to_excel(EXCEL_FILES_DIR + 'cluster_df_for_fig.xlsx')

# pred_proba for each sample:
Y_pred=pd.read_pickle(PRED_RESULTS_DIR+'isCardio/XGB_randomSearch_25_byRepFeatPCA10RelsVDJnocorr0999AgeGender/\
predictions_df.pkl')
Y_pred=Y_pred.loc[sample_list1+sample_list2].rename(columns={'isCardio': 'isCardio_pred_proba'})
Y_pred_real = Y_pred.copy()
Y_pred_real['isCardio']=np.where((Y_pred_real.index.str.replace('BD','').astype(int)) > 949,1,0)
Y_pred_real.sort_values(by=['isCardio','isCardio_pred_proba'],inplace=True)
Y_pred_real['isCardio_pred_proba'] = Y_pred_real['isCardio_pred_proba'].astype('float')

Y_pred_real = Y_pred_real.T.sort_index()
Y_pred_real.to_excel(EXCEL_FILES_DIR + 'Y_pred_real.xlsx')
print ('Y_pred_real.shape: ',Y_pred_real.shape)



# species X disease logp
## add minus to log p-valus if association is negative
summary_df4 = summary_df3.copy()
direction_col = 'species_pos_associated_with_phen'
summary_df4 = convert_logp_values(summary_df4,direction_col,logp_col='-log10p_MW_species_phen',suffix='species_phen')
col_order = (summary_df4.sort_values(by='p_to_show_species_phen', ascending=False))['species2'].unique().tolist()
summary_df4 = summary_df4.drop_duplicates(subset=['cluster','species2'],keep='first')

##
species_disease =summary_df4[['species2','p_to_show_species_phen']].drop_duplicates(subset='species2').\
    set_index('species2').loc[col_order]
species_disease.to_excel(EXCEL_FILES_DIR + 'species_disease_logp_fdr%s%s.xlsx' %(str(FDR_t_species).replace('.',''),
                                                                                         suf))



