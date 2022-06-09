
import pandas as pd
import cPickle as pickle
import os
from scipy.stats import chi2_contingency,ttest_ind

GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR='/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
CARDIO_PHEN_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'

# get annotation frequencies
featureDF = pd.read_excel(DATA_DIR+'TCRfeatures/PNP530Cardio126_annotation_summary_\
freqSum_bd.xlsx').set_index('BD')


# get total seq numbers

#normalize annotation with total seq numbers:

#save in feature folder and change TCR comparison


#take only patients, merge with admission date


#sort by admission data and plot annotation frequencies, look for temporal pattern


#