import pandas as pd
import cPickle as pickle
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from os import makedirs
from os.path import isdir
from ShaniBA.CardioProject.paper_summer2019.fig5.ver2.calc_roc_auc_new_outcome import AucCalculator
from ShaniBA.CardioProject.TCRMicrobiomeBinaryPhenCalculator import get_shared_data_for_2_sample_lists


# directories:
GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
SAMPLE_LIST_DIR = GENIE_BASE_DIR + 'Sample files/BD lists/'
CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5_updated_phens/'
PRED_RESULTS_DIR='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/'
OR_DIR = FIGURE_DIR + 'or_calcs/'

if not isdir(OR_DIR): makedirs(OR_DIR)

#files and dataframe:
outcome_file = CARDIO_PHEN_DIR + 'new_outcome_df.xlsx'
outcome_df = pd.read_excel(outcome_file).set_index('BD')

pred_proba_file =PRED_RESULTS_DIR+'isCardio/XGB_randomSearch_25_byRepFeatPCA10RelsVDJnocorr0999AgeGender/\
predictions_df.pkl'

cluster_data_file = DATA_DIR + 'TCR_seqs_clusters/newX_onlySeqs_025_085_noNans.dat'
cluster_data = pd.read_pickle(cluster_data_file)

cv_df = outcome_df['CV hospitalization including chest pain'].dropna()

cluster_patients_df = cluster_data.loc[cv_df.index.tolist(),:]
cluster_patients_df_binary = (cluster_patients_df > 0).astype(int)


#sample lists:
sample_list1_name = 'Cardio126_CV hospitalization including chest pain00'
sample_list1_file = 'Cardio126_CV hospitalization including chest pain00.pkl'
sample_list2_name = 'Cardio126_CV hospitalization including chest pain11'
sample_list2_file = 'Cardio126_CV hospitalization including chest pain11.pkl'

##### CONSTANT VARIABLES: #####

sample_list_name1 = sample_list1_name
sample_list_name2 = sample_list2_name

with open(SAMPLE_LIST_DIR + sample_list1_file,'rb') as fp1: sample_list1 = pickle.load(fp1)
with open(SAMPLE_LIST_DIR + sample_list2_file,'rb') as fp2: sample_list2 = pickle.load(fp2)


#check sharing rates of interesting clusters:
cluster_list = ['CASSLQQGNTEAFF','CASSPTGSETQYF','CASSLGWGGEQYF','CASSLETGVYEQYF','CASSPTQDYGYTF']
cluster_patients_df_int = cluster_patients_df.loc[:,cluster_list]
cluster_patients_df_int_binary = cluster_patients_df_binary.loc[:,cluster_list]

cluster_count = cluster_patients_df_int_binary.sum()

shared03 = get_shared_data_for_2_sample_lists(cluster_patients_df, sample_list1, sample_list2,
                                       min_shared=0.3, max_shared=1, min_value=0)

shared05 = get_shared_data_for_2_sample_lists(cluster_patients_df, sample_list1, sample_list2,
                                       min_shared=0.5, max_shared=1, min_value=0)

for cluster in cluster_list:
    print cluster
    print cluster in shared05.columns

with open(SAMPLE_LIST_DIR + 'cv_hospt_shared05_clusters.pkl', 'wb') as fp:
    pickle.dump(shared05.columns.tolist(), fp)
