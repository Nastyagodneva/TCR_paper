# imports and definitions
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
PRED_RESULTS_DIR='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5/'
CLUSTER_ANALYSIS_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/sharingAnalysis/seqClusters_allProd_maxdist1/'


def gen_mb_tcr_matrix():
    print ('loading mb data')
    mb_df = pd.read_excel(DATA_DIR + 'Mb_data_by_BD/LargeOrNewGenusSGBs_\
5MSubsampling_00001threshold_Noneratio_sSGB_byBD_PNP530Cardio126_swab.xlsx').set_index('BD')
    print ('mb_df.shape',mb_df.shape)
    print ('loading tcr data')
    cluster_df = pd.read_pickle(DATA_DIR +'TCR_seqs_clusters/newX_onlySeqs_025_085_noNans.dat')
    print ('cluster_df.shape', cluster_df.shape)
    print ('merging')
    mb_cluster_df = pd.merge(mb_df,cluster_df,how='inner',left_index=True,right_index=True)
    mb_cluster_df.to_pickle(DATA_DIR + 'mbSwabSpecies_clusters025085_data.dat')
    print ('mb_cluster_df.shape: ',mb_cluster_df.shape)

    return

gen_mb_tcr_matrix()
# get prediction results:

## hosp1:
# hosp_pred_dir1 = PRED_RESULTS_DIR+ 'outcomes/CVhospitalizationincludingchestpain_XGB_byTCRfeatures_optByAUC/'
# hosp_Results1 = pd.read_pickle(hosp_pred_dir1+'results_df.pkl')
## hosp2:
hosp_pred_dir2 = PRED_RESULTS_DIR+ 'outcomes/CVhospitalizationincludingchestpain_XGB_byTCRfeatures_optByAUC2/'
hosp_Results2 = pd.read_pickle(hosp_pred_dir2+'results_df.pkl')

##unplannedPCI:
PCI_pred_dir = PRED_RESULTS_DIR+ 'outcomes/UnplannedPCI_XGB_byTCRfeatures_optByAUC/'
PCI_Results1 = pd.read_pickle(PCI_pred_dir+'results_df.pkl')
