# imports and definitions
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency,ttest_ind
from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import add_corrected_pValues,compare_TCR_between_groups
import os
from scipy.stats import ttest_ind, mannwhitneyu
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
PRED_RESULTS_DIR='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5/'

outcome_df = pd.read_excel(CARDIO_PHEN_DIR + 'outcomeDF.xlsx').set_index('BD')
phen_df = pd.read_excel(CARDIO_PHEN_DIR + 'PCI_LVEF_CVA_CABG_bin.xlsx').set_index('BD')

merged = pd.merge(outcome_df,phen_df, how='left',left_index=True,right_index=True)

result_df = pd.DataFrame()
n=0
for outcome in outcome_df.columns:
    for phen in phen_df.columns:
        data1 = merged[merged.index.isin(phen_df[phen_df[phen]==0].index)][outcome]
        data2 = merged[merged.index.isin(phen_df[phen_df[phen]==1].index)][outcome]
        mean_data1 = data1.mean()
        mean_data2 = data2.mean()

        result_df.loc[n, 'phen'] = phen
        result_df.loc[n, 'outcome'] = outcome
        result_df.loc[n, 'phen_0_count'] = len(data1)
        result_df.loc[n, 'phen_1_count'] = len(data2)
        result_df.loc[n, 'phen_0_sum_outcome'] = np.sum(data1)
        result_df.loc[n, 'phen_1_sum_outcome'] = np.sum(data2)


        cont_table = pd.merge(
            pd.DataFrame(data1.value_counts(dropna=False).rename(phen+'_0')),
            pd.DataFrame(data2.value_counts(dropna=False).rename(phen+'_1')),
            how='inner', left_index=True, right_index=True)
        cont_table = cont_table.fillna(0)
        try:
            chi, chi_p, dof, expected = chi2_contingency(cont_table)
        except:
            print ('could not execute chi test for phen %s' % phen)
            chi = np.nan;
            chi_p = np.nan;
            dof = np.nan;
            expected = np.nan

        result_df.loc[n, 'p_value'] = round(chi_p, 4)
        result_df.loc[n, 'test'] = 'chi_test'

        n+=1
        
result_df = add_corrected_pValues(result_df, pValueColumn='p_value',
                    nTests=len(result_df), FDR=0.1)
result_df = result_df.sort_values(by='p_value')
result_df.to_excel(os.path.join(FIGURE_DIR,'phen_outcome_associations.xlsx'))