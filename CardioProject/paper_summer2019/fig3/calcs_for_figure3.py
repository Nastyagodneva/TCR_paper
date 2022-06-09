import cPickle as pickle
import pandas as pd
import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns


#definitions:
GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
SAMPLE_LIST_DIR = GENIE_BASE_DIR + 'Sample files/BD lists/'
CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig3_balanced_comparison/calc_fig3/'

####################################################################
# I compared PNP and cardio sub-cohots with 49 subjects each for clusters, species and clusters-species.
# clusters shared by 0.25-1
# species shared by 0.1-1
# threshold for clusters-phen: top 100
# threshold for species-phen: FDR<0.25
# threshold for species-clusters: FDR<0.25


# get summary table:

summary_df = pd.read_excel(FIGURE_DIR + 'isCardio_min_shared_species01\
_max_shared_species1_min_shared_cluster025_max_shared_cluster1/summary_withAnot.xlsx')

#filterout directionOK=0
summary_df2 = summary_df[summary_df['directionOK']==1]

#get species+families:
summary_df2['species2'] = summary_df2['species'].apply(lambda x: '_'.join(x.split('|')[-6:-3]))
summary_df2['genus'] = summary_df2['species'].apply(lambda x: x.split('|')[-5])
summary_df2['known_mb'] = np.where(summary_df2['species2'].str.contains('g__unknown|s__unknown'),0,1)

#filter out unknown species_families:
summary_df2 = summary_df2[summary_df2['known_mb'] == 1]

#count how many sequences
print ('number of unique clusters: ', summary_df2['cluster'].nunique())

summary_df2.to_excel(FIGURE_DIR + 'isCardio_min_shared_species01\
_max_shared_species1_min_shared_cluster025_max_shared_cluster1/summary_processed.xlsx')


summary_df_by_cluster =summary_df2.groupby(['cluster','species2']).mean()
summary_df_by_cluster = summary_df_by_cluster.sort_values(by=['-log10p_MW_cluster_phen','-log10p_MW_cluster_species',
                                                              '-log10p_MW_species_phen'],ascending = False)
summary_df_by_cluster.to_excel(FIGURE_DIR + 'isCardio_min_shared_species01\
_max_shared_species1_min_shared_cluster025_max_shared_cluster1/summary_df_by_cluster.xlsx')

#cluster by genus and check sequence similarity
summary_df_by_genus =summary_df2.groupby(['genus','cluster']).mean()
summary_df_by_genus = summary_df_by_genus.sort_values(by=['-log10p_MW_cluster_phen','-log10p_MW_cluster_species',
                                                              '-log10p_MW_species_phen'],ascending = False)
summary_df_by_genus.to_excel(FIGURE_DIR + 'isCardio_min_shared_species01\
_max_shared_species1_min_shared_cluster025_max_shared_cluster1/summary_df_by_genus.xlsx')


# get top100 cluster list:
cluster_phen = pd.read_excel(FIGURE_DIR + 'isCardio_min_shared_species01\
_max_shared_species1_min_shared_cluster025_max_shared_cluster1/concise_summary_df_cluster_phen.xlsx')

cluster_phen = cluster_phen.sort_values(by='-log10p_MW', ascending=False)
top100 = cluster_phen['cluster'].tolist()[:100]

#filter by species_cluster and species_disease relations of FDR<0.1:
FDR_t_fig=0.1
summary_df3 = summary_df2[(summary_df2['MW_p_corrPval_FDR0.1_cluster_species'] <= FDR_t_fig) &
                          (summary_df2['MW_p_corrPval_FDR0.1_species_phen'] <= FDR_t_fig)]
print ('number of unique clusters in summary_df3: ', summary_df3['cluster'].nunique())
print ('number of unique species in summary_df3: ', summary_df3['species2'].nunique())
print (summary_df3['species2'].unique())


#plot preliminary heatmap:

## add minus to log p-valus if association is negative
summary_df3['cluster_species_logp'] = np.where(summary_df3['species_pos_associated_with_cluster'] == 1,
            summary_df3['-log10p_MW_cluster_species'],summary_df3['-log10p_MW_cluster_species']*-1)
summary_df3['species_phen_logp'] = np.where(summary_df3['species_pos_associated_with_phen'] == 1,
            summary_df3['-log10p_MW_species_phen'],summary_df3['-log10p_MW_species_phen']*-1)
col_order = (summary_df3.sort_values(by='species_phen_logp', ascending=False))['species2'].unique().tolist()


## generate pivot table:
cluster_species_hm = summary_df3.pivot(index='cluster',columns='species2',values='cluster_species_logp')
cluster_species_hm = cluster_species_hm.loc[top100,col_order]

## order df columns by species association with diease:


fig_width = round(0.4*len(cluster_species_hm.columns),0)

fig1,ax1 = plt.subplots(figsize=(fig_width,12))
hm1=sns.heatmap(cluster_species_hm, ax=ax1,linewidth=0.05, yticklabels=True, xticklabels=True,cmap='seismic')
# plt.tight_layout()
ax1.tick_params(labelsize='small')
fig1.savefig(FIGURE_DIR + 'cluster-species_heatmap_fdr%s.png' %str(FDR_t_fig).replace('.',''))

# plot bar plot of species-phen relations
fig2,ax2 = plt.subplots(figsize=(fig_width,2))
df2=summary_df3[['species2','species_phen_logp']].drop_duplicates(subset='species2').\
    set_index('species2').loc[col_order]
df2.plot(kind='bar',ax=ax2,legend=False,color='black')
ax2.set_xticks=[]
fig2.savefig(FIGURE_DIR + 'species_phen_logp_fdr%s.png' %str(FDR_t_fig).replace('.',''))

# generate binary phen file for cluster comparison:
phen_file = DATA_DIR + 'phenotypes_byBD/PNP530_Age_Gender_Male_HbA1C_SmokingStatus\
_Yes_BMI_HDL_Smoking_ever_nonHDL_Cholesterol.xlsx'
phen_df = pd.read_excel(phen_file).set_index('BD')

binary_phen=phen_df.copy()
binary_phen['isOld'] = np.where(binary_phen['Age'] >= 65,1,
                        np.where(binary_phen['Age'] <=30,0,np.nan))
binary_phen['HbA1C_high'] = np.where(binary_phen['HbA1C'] >6.4,1,0)
binary_phen['Obese'] = np.where(binary_phen['BMI'] >=30,1,0)
binary_phen = binary_phen[['Gender_Male','isOld','HbA1C_high','Obese']]
binary_phen.to_excel(DATA_DIR + 'phenotypes_byBD/PNP530_binary_phen.xlsx')

# heatmap for top100 clusters associated with old age, male, high HbA1C, obesity:
top100_phen_t_FDR = 0.2

## get -logp pvalues for each
top100_phen_df = pd.DataFrame(index=top100)
phen_list = ['isOld','Gender_Male','HbA1C_high','Obese']

for phen in phen_list:
    dir1 = FIGURE_DIR + '%s_min_shared_species01_max_shared_species1\
_min_shared_cluster03_max_shared_cluster1/' %phen
    df = pd.read_excel(dir1 + 'concise_summary_df_cluster_phen.xlsx')
    df = df[df['MW_p_corrPval_FDR0.1'] <=top100_phen_t_FDR]
    df['logp'] = np.where(df['Associated_with_%s_positive' %phen] == 1,
                          df['-log10p_MW'],df['-log10p_MW']*-1)
    df = df[['cluster','logp']].set_index('cluster')
    try:
        df = df.loc[top100,:].rename(columns={'logp':'%s_logp' %phen})
    except:
        print ('none of top100 clusters is associated with %s with FDR %s' %(phen,top100_phen_t_FDR))
        df = pd.DataFrame(index=top100,data={'%s_logp' %phen: 100*[np.nan]})
    top100_phen_df = pd.merge(top100_phen_df,df,how='left',left_index=True,right_index=True)
top100_phen_df.to_excel(FIGURE_DIR + 'top100_phen_df_fdr%s.xlsx' %str(top100_phen_t_FDR).replace('.',''))

fig3,ax3 = plt.subplots(figsize=(1,12))
hm3=sns.heatmap(top100_phen_df, ax=ax3,linewidth=0.05, yticklabels=False,
                xticklabels=True,cmap='PiYG')
plt.tight_layout()
ax3.tick_params(labelsize='small')
fig3.savefig(FIGURE_DIR + 'top100_phen_heatmap_fdr%s.png' %str(top100_phen_t_FDR).replace('.',''))



## merge with top 100