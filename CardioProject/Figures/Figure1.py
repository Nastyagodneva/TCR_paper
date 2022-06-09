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
from ShaniBA.CardioProject.Figures.GeneralFigureFunctions import *

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

#### subplot functions:
def set_fig1_definitions():
    params = {
   'axes.labelsize': 16,
   'font.size': 12,
   'legend.fontsize': 10,
   'xtick.labelsize': 12,
   'ytick.labelsize': 12,
   'text.usetex': False,
#    'figure.figsize': [m2inch(183), m2inch(247)],#[4.5, 4.5]
   'figure.dpi': 300,
   'xtick.direction':'out'}


    mpl.rcParams.update(params)
    mpl.rcParams['axes.linewidth'] = 0.5
    mpl.rcParams['axes.edgecolor'] = 'black'
    mpl.rcParams['axes.facecolor'] = 'white'
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['xtick.minor.pad'] = 4

    return

# ax1: age prediction correlation:
def plot_age_prediction_correlation(ax):
    
    #### this function can be converted to more modular function to plot phenotype prediction correlation###
    f1 = PRED_RESULTS_DIR + 'TargetDFs/PNP530_Age.xlsx'
    PNP530_Age = pd.read_excel(f1).set_index('BD')

    dir1 = PRED_RESULTS_DIR + 'PNP530_majorPhenotypes/XGB_presetHyperpar_AgebyRepFeatPCA10RelsVDJnocorr0999/'
    prediction_df_healthy_age = pd.read_pickle(dir1 + 'predictions_df.pkl')


    Age = PNP530_Age[PNP530_Age['Age'].notnull()]['Age'].rename('Real Age').astype('float')
    pred_age = prediction_df_healthy_age['Age'].loc[Age.index].rename('Predicted Age').astype('float')
    merged = pd.merge(pd.DataFrame(Age), pd.DataFrame(pred_age), how='inner', left_index=True, right_index=True)
    # print merged.head()

    # # plot data:
    x = 'Real Age'; y = 'Predicted Age'
    merged.plot(x, y, ax=ax, kind='scatter', alpha=0.5, c='darkred', s=60)
    ax.plot(np.unique(merged[x]), np.poly1d(np.polyfit(merged[x], merged[y], 1))(np.unique(merged[x])), c='black', linewidth=1)
    r, p = MyPearsonr(merged[x], merged[y])
    print r, p
    ax.annotate('r=%s\np=%.1E' % (round(r, 2), p), xy=(0.04, 0.98), xycoords='axes fraction', fontsize=mpl.rcParams['font.size'],
                    xytext=(0, 0), textcoords='offset points', fontweight='bold',
                    ha='left', va='top')
    ax.set_xlim(10, 80)
    ax.set_ylim(10, 80)
    ticks = [20, 40, 60, 80]
    ax.set_xticks(ticks); ax.set_yticks(ticks)
    
    return ax

# ## ax2, ax4: plotting shap summary 
def gen_shap_summary_to_axes(ax, fig, target_name, nTopFeatures, shap_df,
                              features_df, jitter=None, sample_num=None,
                               feature_name_list=None, scalingMethod='minmax',
                               addColorBar=True, alpha=1):
    # ##optional values for scalingMethod: 'minmax','scale','perc',None
    
    
     # calcualte feature order by sum of abs shap values and take top n values:
    if nTopFeatures is not None:
        topN_features = shap_df.abs().sum().sort_values(ascending=False)[:nTopFeatures].index.tolist()
        print topN_features
        topN_features_new = edit_feature_names(topN_features)
        print topN_features_new
    else:
        topN_features = shap_df.abs().sum().sort_values(ascending=False)
    
    # if not None, select limited number of samples to present (used to edit the function and test it quickly)
    if sample_num is None:
        sample_num = shap_df.shape[0]
       
    features_df_topN = features_df.loc[:, topN_features].iloc[:sample_num, :]
    shap_df_topN = shap_df.loc[:, topN_features].iloc[:sample_num, :]
    print ('sample_num= ', sample_num)

    # scale feature data- to enable color coding to show samples high and low in each feature:
    if scalingMethod is not None:
        if scalingMethod == 'minmax':
            # using min-max scaler 
            scaler = preprocessing.MinMaxScaler(copy=True, feature_range=(0, 1))
            features_df_topN_scaled = pd.DataFrame(index=features_df_topN.index, columns=features_df_topN.columns, data=scaler.fit_transform(features_df_topN))
        if scalingMethod == 'scale':
            # center to the mean and to get unit variance
            scaler = preprocessing.scale(copy=True)
            features_df_topN_scaled = pd.DataFrame(index=features_df_topN.index, columns=features_df_topN.columns, data=scaler.fit_transform(features_df_topN))
        elif scalingMethod == 'perc':
            features_df_topN_scaled = features_df_topN / features_df_topN.sum()
        else:
            print 'scaling method was not identified, no scaling was done'
            features_df_topN_scaled = features_df_topN
    else:
        print 'no scaling was done'
        features_df_topN_scaled = features_df_topN
    
    # transform to longform tables and merge:
    features_df_topN_scaled_long = pd.melt(features_df_topN_scaled.reset_index(),
                                           id_vars='index', value_name='feature_value')
    shap_df_topN_long = pd.melt(shap_df_topN.reset_index(), id_vars='index', value_name='shap_value')
    merged = pd.merge(features_df_topN_scaled_long, shap_df_topN_long, how='inner', left_on=['index', 'variable'], right_on=['index', 'variable']).rename(columns={'index':'Sample', 'variable':'feature'})
    
    print ('features_df_topN.shape=', features_df_topN.shape)
    print ('features_df_topN_scaled.shape=', features_df_topN_scaled.shape)
    print ('features_df_topN_scaled_long.shape=', features_df_topN_scaled_long.shape)
    print ('shap_df_topN.shape=', shap_df_topN.shape)
    print ('shap_df_topN.shape=', shap_df_topN.shape)
    print ('merged.shape=', merged.shape)
    print ('count shap values per feature:')
    print merged.head(20)
    print 'merged groupby feature - count:'
    print merged[['feature', 'shap_value']].groupby('feature').count()
    



    # ##plot:
    ax.grid(False)
    # generate unuseful scatter plot just to give colorbar to color range...
    stam = plt.figure()
    stam1 = plt.scatter(merged['shap_value'], merged['feature_value'] , c=merged['feature_value'], cmap='RdBu_r')
    plt.clf()
    
    # generate specific ax for the colorbar and plot the colorbar:
    if addColorBar:
        divider = make_axes_locatable(ax)
        ax_cb = divider.new_horizontal(size="4%", pad=0.05)
        fig.add_axes(ax_cb)
        cb1 = fig.colorbar(stam1, cax=ax_cb, orientation='vertical')
        ax_cb.yaxis.set_ticks([0, 0.5, 1])
        ax_cb.yaxis.set_ticklabels([0, 0.5, 1], fontsize=mpl.rcParams['ytick.labelsize'] - 1)
        ax_cb.set_ylabel('Norm Feature Value', fontsize=mpl.rcParams['axes.labelsize'] - 3)
        
    # plot the scatter plot itself:
    plot = sns.stripplot(x='feature', y='shap_value', hue='feature_value', data=merged,
                      palette='RdBu_r',
                     alpha=alpha, ax=ax, s=5, jitter=jitter)
        
#         .swarmplot   edgecolor='grey'  ,linewidth=2
#     else:
#         plot = sns.stripplot(x='feature', y='shap_value', hue='feature_value', data=merged,
#                       palette='RdBu_r',
#                        edgecolor='none', alpha=.4, ax=ax, s=5)
    
    plot.get_legend().set_visible(False)
    ax.set_xlabel("")
    ax.set_ylabel('Shap Value (%s)' % target_name)
    ax.set_xticklabels(topN_features_new)  # change ticklables to readable feature names
    plt.setp(ax.get_xticklabels(), rotation=90, ha="left",
         rotation_mode="anchor", fontsize=mpl.rcParams['font.size'] - 4, fontweight='bold')
    
    # color ticklables according to their content:
    for n, lab in enumerate(ax.get_xticklabels()):
        str_lab = str(lab)
#         print n,str(lab)
        if ('insertion' in str_lab) or ('deletion' in str_lab) or ('length' in str_lab):
            plt.setp(ax.get_xticklabels()[n], color='green')
        elif ('Norm' in str_lab) or ('Top' in str_lab) or ('count' in str_lab)\
or ('Shannon' in str_lab) or ('Berger' in str_lab) or ('Simpson' in str_lab) :
            plt.setp(ax.get_xticklabels()[n], color='red')
        elif ('seq #' in str_lab) or ('seq freq.' in str_lab) or ('annotate' in str_lab) :
            plt.setp(ax.get_xticklabels()[n], color='orange')
        elif ('non-prod' in str_lab) :
            plt.setp(ax.get_xticklabels()[n], fontstyle='italic')

#             plt.setp(ax.get_xticklabels()[n], r'\underline{Parameters}: ')
            
    ax.set_xticks([x + 0.25 for x in ax.get_xticks()])
    ax.tick_params(axis='x', direction="in", pad=-3, labelright=True)


    return ax, fig

# ax3: gender prediction correlation
def plot_gender_prediction_ROC_PR(ax):
    # # get data:
    dir2 = PRED_RESULTS_DIR + 'PNP530_majorPhenotypes/XGB_presetHyperpar_Gender_MalebyRepFeatPCA10RelsVDJnocorr0999/'
    gender_pred_proba = pd.read_pickle(dir2 + 'predictions_df.pkl').astype('float').rename(columns={'Gender_Male':'pred_proba'})

    # y_pred_proba_healthy_gender.head()
    gender_file = PRED_RESULTS_DIR + 'TargetDFs/PNP530_Gender.xlsx'
    gender = pd.read_excel(gender_file).set_index('BD').astype('float').rename(columns={'Gender_Male':'Gender'})

    merged = pd.merge(gender, gender_pred_proba, how='inner', left_index=True, right_index=True)


    # #plot:
    pos_label = 1
    ax, inset_axes,roc_auc, pr_auc, prevalence = plot_ROC_PR_AUC(y=pd.DataFrame(merged['Gender']), y_pred_df=pd.DataFrame(merged['pred_proba']),
                                    ax=ax, color1='darkred', color2='grey', ticklabelsize=mpl.rcParams['xtick.labelsize'],
                                      textsize=mpl.rcParams['font.size'], labelsize=mpl.rcParams['axes.labelsize'], add_texts=False)
    
    ax.annotate('ROC AUC=%s\nPR=%s\nPrevalence=%s' % (round(roc_auc, 3), round(pr_auc, 2), round(prevalence, 2)), xy=(0.04, 0.98), xycoords='axes fraction',
                 fontsize=mpl.rcParams['font.size'], xytext=(0, 0), textcoords='offset points', fontweight='bold', ha='left', va='top')
    ax.set_xlabel('False Positive Rate (Gender)')
    ax.set_ylabel('True Positive Rate (Gender)')
    
    return ax

##### plot figure:

if __name__ == "__main__":

    set_fig1_definitions()
    
    # ## adefine figure structure:
    fig = plt.figure(figsize=(18, 8))
    gs1 = gridspec.GridSpec(2, 4, wspace=0.5, hspace=0.5)
    
    # ## plot ax1:
    ax1 = plt.Subplot(fig, gs1[0, 0])
    ax1.tick_params(axis='x', direction="out")
    fig.add_subplot(ax1)
    ax1 = plot_age_prediction_correlation(ax1)
    
    # ## define parameters for ax2:
    ax2 = plt.Subplot(fig, gs1[0, 1:])
    fig.add_subplot(ax2)
    
    # load shap values:
    dir1 = PRED_RESULTS_DIR + 'PNP530_majorPhenotypes/XGB_presetHyperpar_AgebyRepFeatPCA10RelsVDJnocorr0999/'
    Age_pred_shap = pd.read_pickle(dir1 + 'shap_values.pkl')['Age']
    
    # load feature values
    Xfeatures_file = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Predictions/featureDFs/TCRfeaturesNoT_PCA10percShared_withRels_VDJ_noCorr0-999_nanFilled_noConsts.dat'
    Xfeatures = pd.read_pickle(Xfeatures_file).loc[PNP530, :]
    for col in Xfeatures.columns:
        Xfeatures[col] = Xfeatures[col].fillna(Xfeatures[col].median())
    
    nTopFeatures = 20
    shap_df = Age_pred_shap
    features_df = Xfeatures
    jitter = True
    sample_num = None
    scalingMethod = 'perc'
    target_name = 'Age'
    
    
    # ## plot ax2:
    ax2, fig = gen_shap_summary_to_axes(ax2, fig, target_name, nTopFeatures, shap_df, features_df, jitter=jitter, sample_num=sample_num, scalingMethod=scalingMethod)
    
    # ##plot ax3:
    ax3 = plt.Subplot(fig, gs1[1, 0])
    fig.add_subplot(ax3)
    ax3 = plot_gender_prediction_ROC_PR(ax3)
    
    
    # ## define parameters for ax4:
    ax4 = plt.Subplot(fig, gs1[1, 1:])
    fig.add_subplot(ax4)
    
    # load shap values:
    dir2 = PRED_RESULTS_DIR + 'PNP530_majorPhenotypes/XGB_presetHyperpar_Gender_MalebyRepFeatPCA10RelsVDJnocorr0999/'
    Gender_pred_shap = pd.read_pickle(dir2 + 'shap_values.pkl')['Gender_Male']
    
    # load feature values
    Xfeatures_file = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Predictions/featureDFs/TCRfeaturesNoT_PCA10percShared_withRels_VDJ_noCorr0-999_nanFilled_noConsts.dat'
    Xfeatures = pd.read_pickle(Xfeatures_file).loc[PNP530, :]
    for col in Xfeatures.columns:
        Xfeatures[col] = Xfeatures[col].fillna(Xfeatures[col].median())
    
    nTopFeatures = 20
    shap_df = Gender_pred_shap
    features_df = Xfeatures
    jitter = True
    sample_num = None
    scalingMethod = 'perc'
    target_name = 'Gender'
    
    # ## plot ax4:
    ax4, fig = gen_shap_summary_to_axes(ax4, fig, target_name, nTopFeatures, shap_df, features_df, jitter=jitter, sample_num=sample_num, scalingMethod=scalingMethod)
    
    # ## add letters:
    ax1 = add_subfigure_letter(ax1, 'A', fontweight='bold')
    ax2 = add_subfigure_letter(ax2, 'B', fontweight='bold')
    ax3 = add_subfigure_letter(ax3, 'C', fontweight='bold')
    ax4 = add_subfigure_letter(ax4, 'D', fontweight='bold')
    
    # ## save figure:
    # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.45, hspace=0.5)
    fig.savefig(FIG1_DIR + 'figure1_%s.png' % cdate, dpi=300)
    
    # Save axes seperately
    for subp in [('ax1', ax1), ('ax2', ax2), ('ax3', ax3), ('ax4', ax4)]:
        extent = subp[1].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
        fig.savefig(FIG1_DIR + '%s_figure_expanded_%s.png' % (subp[0], cdate), bbox_inches=extent.expanded(1.3, 1.35))
    
    print 'Finished figure 1!!'
    # plt.show()
    



