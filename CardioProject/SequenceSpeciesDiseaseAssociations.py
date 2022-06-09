# Goals:
# 1. calculate species/genus-disease association in balanced male cohorts
# 2. calculate species/gender-cluster associations in controls+patients
# 3. merge with cluster-disease associations in balanced male cohorts.
# 4. save to file


# imports:

import pandas as pd
import os
import time
import pickle

# Path definitions:
BASE_PATH = os.path.join('/net', 'mraid08', 'export')
GENIE_PATH = os.path.join(BASE_PATH,'genie')
JAFAR_PATH = os.path.join(BASE_PATH,'jafar','Microbiome','Analyses','ShaniBAF')
MY_PATH = os.path.join(GENIE_PATH,'Lab','Personal','ShaniBAF')

SAMPLE_LIST_DIR = os.path.join(MY_PATH,'Sample files','BD lists')
# CARDIO_PHEN_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'
CLUSTER_ANALYSIS_DIR = os.path.join(JAFAR_PATH,'sharingAnalysis',
                                    'seqClusters_allProd_maxdist1')
RESULT_DIR=os.path.join(MY_PATH,'TCR_real_data','TCR_mb_results',
                        'TCRclusters_MBsega_associations_April2019')

DATA_DIR = os.path.join(JAFAR_PATH,'Microbiome','Analyses','TCR')
# other definitions:
cdate = str(time.strftime("%d%m%Y"))

# SAMPLE LISTS:
with open(os.path.join(SAMPLE_LIST_DIR,'PNP530'),'rb') as fp:
    PNP530 = pickle.load(fp)
with open(os.path.join(SAMPLE_LIST_DIR,'Cardio126'),'rb') as fp:
    Cardio126 = pickle.load(fp)
PNP530Cardio126 = PNP530 + Cardio126



##specific parameters:

# Get data:

## get cluster-disease associations:
### get file
# Cluster_disease_associations_balanced = pd.read_excel(os.path.join(JAFAR_PATH,'sharingAnalysis',
#                             'seqClusters_allProd_maxdist1','Fisher_MW_comparisons',
# 'Fisher_MW_results_PNP530_balancedAge_males_prodClusRed_Cardio126_balancedAge_males_prodClusRed_percShared01_percTooMany085.xlsx'))
## get only with mw p corrected < 0.2
## calcualte enriched where by OR value
### get only cluster_head, p_MW, p_corrected, enriched where
### calculate -log10 p_value

print 'Done!'

## get species and genus data for relevant participants
sega_species_genotek_005095=pd.read_excel('/net/mraid08/export/jafar/Microbiome/Analyses/TCR/Mb_data_by_BD/\
LargeOrNewGenusSGBs_5MSubsampling_00001threshold_Noneratio_sSGB_byBD_PNP530_005_095_genotek.xlsx').set_index('BD')
print ('sega_species_genotek_005095.shape: ',sega_species_genotek_005095.shape)

sega_species_swab_005095=pd.read_excel('/net/mraid08/export/jafar/Microbiome/Analyses/TCR/Mb_data_by_BD/\
LargeOrNewGenusSGBs_5MSubsampling_00001threshold_Noneratio_sSGB_byBD_PNP530_005_095_swab.xlsx').set_index('BD')
print ('sega_species_swab_005095.shape: ',sega_species_swab_005095.shape)

sega_genus_genotek_005095=pd.read_excel('/net/mraid08/export/jafar/Microbiome/Analyses/TCR/Mb_data_by_BD/\
LargeOrNewGenusSGBs_5MSubsampling_00001threshold_Noneratio_gSGB_byBD_PNP530_005_095_genotek.xlsx').set_index('BD')
print ('sega_genus_genotek_005095.shape: ',sega_genus_genotek_005095.shape)

sega_genus_swab_005095=pd.read_excel('/net/mraid08/export/jafar/Microbiome/Analyses/TCR/Mb_data_by_BD/\
LargeOrNewGenusSGBs_5MSubsampling_00001threshold_Noneratio_gSGB_byBD_PNP530_005_095_swab.xlsx').set_index('BD')
print ('sega_genus_swab_005095.shape: ',sega_genus_swab_005095.shape)



## get cluster presence/absence matrix for relevant participants


# run tests:



# save to file