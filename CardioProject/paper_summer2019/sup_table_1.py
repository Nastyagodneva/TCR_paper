import pandas as pd
from scipy.stats import chi2_contingency, ttest_ind

GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR='/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
PAPER_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/June2019/'

#get phen_df
BASIC_PHEN = DATA_DIR + 'phenotypes_byBD/PNP530Cardio126_Age_Gender_Male_HbA1C_SmokingStatus_Yes_BMI_HDL_Smoking_ever_nonHDL_Cholesterol.xlsx'
phen_df = pd.read_excel(BASIC_PHEN).set_index('BD')
phen_df_PNP530 = phen_df[(phen_df.index.str.replace('BD','').astype(int))<950]
phen_df_Cardio126 = phen_df[(phen_df.index.str.replace('BD','').astype(int))>=950]

#calculate mean and std:
phen_df_PNP530_stats = phen_df_PNP530.describe().loc[['mean','std']].T
phen_df_Cardio126_stats = phen_df_Cardio126.describe().loc[['mean','std']].T
phen_df_stats_merged = pd.merge(phen_df_PNP530_stats,phen_df_Cardio126_stats,
                                how='inner', left_index=True,right_index=True,
                                suffixes=['_PNP530','_Cardio126'])

#calculate ttest or chi test:

sup_table = phen_df_stats_merged
for col in phen_df.columns:

    if phen_df[col].nunique()==2:
        cont_table = pd.merge(
            pd.DataFrame(phen_df_PNP530[col].value_counts(dropna=False).rename('PNP530')),
            pd.DataFrame(phen_df_Cardio126[col].value_counts(dropna=False).rename('Cardio126')),
            how='inner',left_index=True,right_index=True)
        chi, chi_p, dof, expected = chi2_contingency(cont_table)
        sup_table.loc[col,'p_value'] = round(chi_p,4)
        sup_table.loc[col, 'test'] = 'chi_test'
    else:
        PNP530_data = phen_df_PNP530[col][phen_df_PNP530[col].notnull()].tolist()
        Cardio126_data = phen_df_Cardio126[col][phen_df_Cardio126[col].notnull()].tolist()
        s,p = ttest_ind(PNP530_data,Cardio126_data)
        sup_table.loc[col, 'p_value'] = round(p, 4)
        sup_table.loc[col, 'test'] = 't_test'
        
#save table:
sup_table.to_excel(PAPER_DIR+'sup_table1.xlsx')