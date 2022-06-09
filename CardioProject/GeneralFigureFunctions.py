#### imports:

# system & general:
from os import listdir, mkdir, makedirs
from os.path import isfile, join, isdir, exists
import cPickle as pickle
import os
import re
import time

# data analysis and statistics:
import pandas as pd
import numpy as np
from scipy import stats
import random
from scipy.stats import pearsonr, fisher_exact
import math
from scipy.spatial.distance import braycurtis, pdist, euclidean

# figures:
import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns 
from matplotlib.ticker import FormatStrFormatter
from matplotlib import gridspec
# from PNPChip.ForPaper.Figures.nature_guidline_utils import m2inch
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# my functions:
from ShaniBA.myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot, adjusted_roundup
from ShaniBA.MyFunctionsShani import *
from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import * 
from ShaniBA.TCR_microbiome_interactions.TCR_microbiome_interactions_functions import *
from ShaniBA.TCR_microbiome_interactions.TCR_microbiome_interactions_functions2 import *
from ShaniBA.SampleLists.SampleFileFunctions import *
from ShaniBA.PhenotypicData.PhenotypeGenerationFunctions import *
from ShaniBA.CardioProject.CardioFunctions import *
from ShaniBA.PredictionPipeline.PredictionFunctions import *
from ShaniBA.TCR_feature_generation.SubsamplingFunctions import *

#####path definitions:
MyPath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF'

with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP530', 'rb') as fp:
    PNP530 = pickle.load(fp)
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/Cardio126', 'rb') as fp:
    Cardio126 = pickle.load(fp)
PNP530Cardio126 = PNP530 + Cardio126

FIG1_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Presentations and Manuscripts/CardioTCR paper/FigureDraft_Jan19/\
Fig1_AgeGenderPrediction/'
PRED_RESULTS_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/'

##### general definitions:
pd.set_option('display.expand_frame_repr', False)
pd.set_option('display.width', 1000)
np.set_printoptions(precision=4, linewidth=200)

cdate = str(time.strftime("%d%m%Y"))
cdate

def add_subfigure_letter(ax, letter, fontsize=18, x_frac_offset=0,y_frac_offset=1.15,
                          **text_kw):
    ax.annotate(letter, xy=(x_frac_offset, y_frac_offset), xycoords='axes fraction', fontsize=fontsize,
                xytext=(0, 0), textcoords='offset points',
                ha='right', va='top', **text_kw)
    return ax

def plot_many_roc(ax,sample_list,color_map,prediction_subdir_list=None):
    Y_file='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/TargetDFs/isCardio.dat'
    Y=pd.read_pickle(Y_file).loc[sample_list,:]
    # print Y.head()

    #plot best prediction:
    d1='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/isCardio/'
    prediction_subdir_list=['XGB_randomSearch_25_byPredictedGender/', 'XGB_randomSearch_25_byPredictedAge/',
               'XGB_randomSearch_25_byPredictedAgeGender/','XGB_randomSearch_25_byRepFeatPCA10percVDJ0999PredictedAgeGender/']
    cmap=plt.get_cmap(color_map, 5)
    color_list=[cmap(1),cmap(2),cmap(3),cmap(4)]
    # print colorlist
    dataset_name_list=['Pred. Gender','Pred. Age','Pred. Age+Gender','TCR features + Pred. Age+Gender']
    for n in range(len(prediction_subdir_list)):
        d=prediction_subdir_list[n]
        color=color_list[n]
    #     print ("color=",color)

        pred=pd.read_pickle(d1+d+'predictions_df.pkl')
    #     print pred.head()
        y_true=Y.isCardio
        y_pred=pred.isCardio.loc[sample_list]
        pos_label=1
        dataset_name=dataset_name_list[n]

        ax, roc_auc=plot_roc_curve(ax,y_true,y_pred,pos_label,dataset_name=dataset_name,color=color)
        if n==2: auc_AgeGender=roc_auc
        elif n==3: auc_TCR=roc_auc
    AUC_diff=round(auc_TCR-auc_AgeGender,3)
    handles, labels = ax.get_legend_handles_labels()
    labels=[l.replace('area = ','') for l in labels]
    ax.legend(handles, labels,loc='3',fontsize='small')
    # ax.legend(loc='best',bbox_to_anchor=(1.05, 1))

    return ax,AUC_diff

def edit_feature_names(feature_list):
    # get readable feature names
    replace_tuple = (
        ('_1', '') , ('_0', '(non-prod)') , ('clonal', '') , ('normSeqNums_per2000', 'Norm seq #') , ('_NT', ' nt') , ('_aa', ' aa'),
        ('_rel_seq_count', ' seq #'), ('_mean', ' mean'), ('_std', ' std'), ('cdr3Length', 'CDR3 length'), ('Multiplesclerosis', 'MS '), ('Colorectalcancer', 'Colorectal cancer'),
        ('_rel_cum_freq(perc)', ' seq freq.'), ('Deletion', ' deletion'), ('Insertion', ' insertion'), ('top', 'Top'), ('mean_nt_per', 'Mean nt per'),
        ('totalAnnotatefreqs', 'Total annotated Seqs'), ('shannon', 'shannon '),('(perc)',''),('Influenza','Influenza-related'),
        ('cum','cumul.'),
        ('_', '-'), ('-', ' ')
        )
    feature_list_new = feature_list
    for r in replace_tuple:
        feature_list_new = [feature.replace(*r) for feature in feature_list_new] 
    feature_list_new = [feature.capitalize() for feature in feature_list_new]
    
    return feature_list_new

def remove_spines(ax=None,removeFigBorders=False,removeTicks=True,removeTicklabels=True):
    if ax is None:
        ax = plt.gca()
    if removeTicks:
        ax.tick_params(axis='both', which='both', bottom=False, top=False,
                        left=False, right=False) 
    if removeTicklabels:
        ax.tick_params(labelbottom=False, labelleft=False)
    if removeFigBorders:
        for spine in ax.spines.values():
                spine.set_edgecolor('white')
            
    return ax



