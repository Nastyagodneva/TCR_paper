### imports:
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
# mpl.use('Agg')

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
from ShaniBA.CardioProject.Figures.Figure1 import *

#ML:
from sklearn.linear_model import LogisticRegression, LinearRegression


### definitions:
#####path definitions:
MyPath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF'

with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP530', 'rb') as fp:
    PNP530 = pickle.load(fp)
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/Cardio126', 'rb') as fp:
    Cardio126 = pickle.load(fp)
PNP530Cardio126 = PNP530 + Cardio126

FIG2_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Presentations and Manuscripts/CardioTCR paper/FigureDraft_Jan19/\
Fig2_AgeGenderPrediction/'
PRED_RESULTS_DIR='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/'
FEATURES_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Predictions/featureDFs/'
SAMPLE_LIST_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/'
CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'


##### general definitions:
pd.set_option('display.expand_frame_repr', False)
pd.set_option('display.width', 1000)
np.set_printoptions(precision=4, linewidth=200)

cdate = str(time.strftime("%d%m%Y"))
cdate

#### FUNCTIONS: #######

def set_fig2_definitions():
    params = {
   'axes.labelsize': 16,
   'font.size': 12,
   'legend.fontsize': 14,
    'axes.titlesize':16,
    'axes.titleweight':'bold',
   'xtick.labelsize': 12,
   'ytick.labelsize': 12,
   'text.usetex': False,
#    'figure.figsize': [m2inch(183), m2inch(247)],#[4.5, 4.5]
#    'figure.dpi': 300,
   'xtick.direction':'out'}


    mpl.rcParams.update(params)
    mpl.rcParams['axes.linewidth'] = 0.5
    mpl.rcParams['axes.edgecolor'] = 'black'
    mpl.rcParams['axes.facecolor'] = 'white'
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['xtick.minor.pad'] = 4


    return

#wrap the previous function with definition of input variables (include shuffled)

def plot_isCardio_multi_ROC(ax):
    
    sample_list=PNP530+Cardio126
    color_map='Reds'
    
    ### plot only the interesting prediction (to get single legend for the plot):
    prediction_subdir_list1=['XGB_randomSearch_25_byRepFeatPCA10percVDJ0999PredictedAgeGender/']
    dataset_name_list1=['TCR features + Pred. Age,Gender']
    
    #plot:
    ax,AUC_diff=plot_many_roc(ax,sample_list,color_map,prediction_subdir_list=prediction_subdir_list1)
#     ,prediction_dir=None,
#                      dataset_name_list=dataset_name_list1,AUC_diff_min='byRepFeat',AUC_diff_max='byRepFeat')
    
    #generate legend:
    handles, labels = ax.get_legend_handles_labels()
    replace_list=[(dataset_name_list1[0],'AUC='),('(',''),(')',''),('area = ','')]
    for r in replace_list:
        labels=[x.replace(r[0],r[1]) for x in labels]
    #add legend manually, as its impossible to add automatically two different legends for the same plot
    first_legend=ax.legend(handles, labels,loc='lower right')  
    ax.add_artist(first_legend)
    
    ### plot ROC for all predictions,without a legend as the legend is identical to the ax2 legend:
    prediction_subdir_list2=['XGB_randomSearch_25_byPredictedGender/', 'XGB_randomSearch_25_byPredictedAge/',
                   'XGB_randomSearch_25_byPredictedAgeGender/','XGB_randomSearch_25_byOldXShuffleWithPredictedAgeGender/',
                            'XGB_randomSearch_25_byRepFeatPCA10percVDJ0999PredictedAgeGender/']
    dataset_name_list2=['Pred. Gender','Pred. Age','Pred. Age,Gender','Shuff features + Pred. Age,Gender',
                       'TCR features + Pred. Age,Gender']
    
    #plot:
    ax,AUC_diff=plot_many_roc(ax,sample_list,color_map,prediction_subdir_list=prediction_subdir_list2)
#     ,prediction_dir=None,
#                      dataset_name_list=dataset_name_list2,AUC_diff_min=None,AUC_diff_max=None)
    
    # 'remove' legend (because the legend of ax2 will be used. make sure the legend of ax1 and ax2 remains identical! otherwise, need to add seperate legend here)
    labels2=[]
    handles2=[]
    ax.legend(handles2, labels2)
    
    ax.set_title('ROC curve\ndeltaAUC=%s'%AUC_diff,fontweight='bold')
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    
    return ax

def plot_isCardio_multi_PR(ax):
    
    sample_list=PNP530+Cardio126
    color_map='Reds'
    
    ### plot only best prediction, because I want to have 2 different legends for the same plot:
    prediction_subdir_list1=['XGB_randomSearch_25_byRepFeatPCA10percVDJ0999PredictedAgeGender/']
    dataset_name_list1=['TCR features + Pred. Age,Gender']
    
    #plot:
    ax,AUC_diff=plot_many_PR(ax,sample_list,color_map,prediction_subdir_list=prediction_subdir_list1,prediction_dir=None,
                 dataset_name_list=dataset_name_list1,AUC_diff_min='byRepFeat',AUC_diff_max='byRepFeat')
    #edit legened:
    handles, labels = ax.get_legend_handles_labels()
    replace_list=[(dataset_name_list1[0],'AUC='),('(',''),(')','')]
    for r in replace_list:
        labels=[x.replace(r[0],r[1]) for x in labels]
    #add legend manually:
    first_legend=ax.legend(handles, labels,loc='lower right')  
    ax.add_artist(first_legend)
    
    ###plot ALL datasets this time:
    prediction_subdir_list2=['XGB_randomSearch_25_byPredictedGender/', 'XGB_randomSearch_25_byPredictedAge/',
                   'XGB_randomSearch_25_byPredictedAgeGender/','XGB_randomSearch_25_byOldXShuffleWithPredictedAgeGender/',
                            'XGB_randomSearch_25_byRepFeatPCA10percVDJ0999PredictedAgeGender/']
    dataset_name_list2=['Pred. Gender','Pred. Age','Pred. Age,Gender','Shuff features + Pred. Age,Gender',
                       'TCR features + Pred. Age,Gender']
    #plot:
    ax,AUC_diff=plot_many_PR(ax,sample_list,color_map,prediction_subdir_list=prediction_subdir_list2,prediction_dir=None,
                 dataset_name_list=dataset_name_list2,AUC_diff_min=None,AUC_diff_max=None)
    #edit legend:
    handles2, labels2 = ax.get_legend_handles_labels()
    labels2=[x.split('(')[0] for x in labels2]
   
    handles2=handles2[1:] # remove the first handle as it repeats the last handle
    labels2=labels2[1:] # remove the corresponsing label
    ax.legend(handles2, labels2,loc='best',bbox_to_anchor=(1.05, 1))
    
    ax.set_title('Precision-Recall curve\ndeltaAUC=%s'%AUC_diff,fontweight='bold')
    
    return ax


def shap_scatter_3d(ax,nToShow=20):
    from mpl_toolkits.mplot3d import Axes3D
    
#     ax = Axes3D(ax)
    
    ### get shap data for age,gender,isCardio and calculate sum of abs values per feature:
    shap_age=pd.read_pickle(PRED_RESULTS_DIR + 'PNP530_majorPhenotypes/XGB_presetHyperpar_AgebyRepFeatPCA10RelsVDJnocorr0999/shap_values.pkl')['Age']
    shap_age_values=pd.DataFrame(shap_age.abs().sum(),columns=['shap_age'])
    shap_gender=pd.read_pickle(PRED_RESULTS_DIR + 'PNP530_majorPhenotypes/XGB_presetHyperpar_Gender_MalebyRepFeatPCA10RelsVDJnocorr0999/shap_values.pkl')['Gender_Male']
    shap_gender_values=pd.DataFrame(shap_gender.abs().sum(),columns=['shap_gender'])
    shap_isCardio=pd.read_pickle(PRED_RESULTS_DIR + 'isCardio/XGB_randomSearch_25_byRepFeatPCA10percVDJ0999PredictedAgeGender/shap_values.pkl')['isCardio']
    shap_isCardio_values=pd.DataFrame( shap_isCardio.abs().sum(),columns=['shap_isCardio'])
        
    ###calculate multiple regression and res (isCardioShap prediction vs isCardioShap real)
    merged=pd.merge(shap_age_values,shap_gender_values,how='inner',left_index=True,right_index=True)
    merged=pd.merge(merged,shap_isCardio_values,how='inner',left_index=True,right_index=True)
    
    lm = LinearRegression()  # intercept=True, normalize=False
    lm.fit(merged[['shap_age','shap_gender']],merged['shap_isCardio'])

    shap_isCardio_pred = pd.DataFrame(index=merged.index.tolist(), data=lm.predict(merged[['shap_age','shap_gender']]),columns=['shap_isCardio_pred'])
    merged=pd.merge(merged,shap_isCardio_pred,how='inner',left_index=True,right_index=True)
     
    #add color list according to residule values:
    merged['res'] = merged['shap_isCardio'] - merged['shap_isCardio_pred']
    merged=merged.sort_values(by='res',ascending=False)
    colorList=['r' if n<nToShow else 'black' for n in range(len(merged))]
    merged['color']=colorList
    topShapFeatures=merged[merged['color']=='r'].index.tolist()
    print merged.shape
    print merged.head(20)
    
    ### plot all shap values
    ax.scatter(merged['shap_age'],merged['shap_gender'],merged['shap_isCardio'],c=merged['color'],s=40,alpha=0.4)

    ###plot prediction surface:
    xmin=-500; xmax=1500
    ymin=0; ymax=300
    x_surf = np.arange(xmin,xmax,10)               # generate a mesh
    y_surf = np.arange(ymin,ymax,10) 
    x_surf, y_surf = np.meshgrid(x_surf, y_surf)

    exog = pd.core.frame.DataFrame({'x': x_surf.ravel(), 'y': y_surf.ravel()})
    out = lm.predict(exog)
    ax.plot_surface(x_surf, y_surf,
                    out.reshape(x_surf.shape),
                    rstride=1,
                    cstride=1,
                    color='None',
                    alpha = 0.1)
    
    ###edit plot
    ax.set_xlabel('Age', labelpad=0)
    ax.set_ylabel('Gender', labelpad=0)
    ax.set_zlabel('isCardio', labelpad=0)
    
    ax.set_xlim(xmin,xmax); ax.set_xticks([])
    ax.set_ylim(ymin,ymax); ax.set_yticks([])
    ax.set_zlim(0,150); ax.set_zticks([])
    ax.set_title('Shap values',fontsize=16,fontweight='bold')
    
    return ax,topShapFeatures

def plot_balanced_age_distrbution(ax):

    #get balanced sample lists
    with open(SAMPLE_LIST_DIR+'PNP530_balancedAge_males','rb') as fp:
        PNP530_balanced=pickle.load(fp)
    with open(SAMPLE_LIST_DIR+'Cardio126_balancedAge_males','rb') as fp:
        Cardio126_balanced=pickle.load(fp)

    #get age data:
    Age=pd.read_excel('%s/TCR_real_data/PNP530Cardio126Combined/Phenotypes/\
PNP530Cardio126_Age.xlsx' %MyPath).set_index('BD')

    dataList=[('Healthy (n=%s)' %len(PNP530_balanced),Age.loc[PNP530_balanced,'Age'].tolist()),('Patients (n=%s)' %len(Cardio126_balanced),Age.loc[Cardio126_balanced,'Age'].tolist())]
    title=''

    ax,ks_p_cohort1_cohort2,t_p_cohort1_cohort2,p_Anov,filename=plotHistComprison(dataList,ax,title,showLegend=False,nBins=10,toAnnotate=False,
                          colorList=['Grey','darkred'],alpha=None,
                                               plotType='hist')
    
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels,loc='upper left',bbox_to_anchor=(-0.03,1.02))
    ax.set_title('Age Distributions\n(t-test: p=n.s.)')
    ax.set_yticks([0.1,0.2])
    ax.set_ylabel('Frequency')

    return ax

def plot_pvalues_feature_comparison(ax):
    
    #get pvalue excel data:
    feat_comp=pd.read_excel('%s/TCR_real_data/PNP530Cardio126Combined/featureSummaryDFs/\
feature_comparison_MW_PNP530Cardio126_balancedAge_males_incAnnot.xlsx' %MyPath)
    feat_comp['-log10_p']=-np.log10(feat_comp['p_MW'])

    #plot pvalues in log scale
    feat_comp['x']=1
    
    
    interesting_features=['n2Insertion_mean_1','Influenza_seq_count','top1000clonal_nt_1','Influenza_cum_freq(perc)']
    feat_comp['color']=np.where(feat_comp.index.isin(interesting_features),1,0)
    
    cpal = sns.color_palette(['black','orange'])
#     , desat=0.2
    sns.swarmplot(x='x', y='-log10_p', hue='color', data=feat_comp, edgecolor='none',  ax=ax, s=6,palette=cpal,alpha=0.7)
    
    #plot FDR line:
    ax.axhline(1.585,color='navy',linestyle='--')
    ax.text(-0.47,1.587,'FDR=0.1',color='navy',ha='left',va='bottom',fontweight='bold')
    
    #edit plot:
#     ax.set_yscale('log')
#     ax.set_ylim(0.00001,1)
    ax.legend([],[])
    ax.set_ylabel('P-value (log)')
    ax.set_xlabel('')
    ax.set_xticks([])
    ax.set_title('Feature Comparison\nResults(MW Test)')
    
    return ax

def plot_feat_comp_kde(axList):
    
    ###get data:
    feature_list=['n2Insertion_mean_1','PC1','top1000clonal_nt_1','V09_J02_0',
                  'CMV_rel_seq_count','SIV_rel_cum_freq(perc)','CalcifiedAorticStenosisdisease_rel_seq_count',
                  'Melanoma_rel_cum_freq(perc)']
    feature_name_list=edit_feature_names(feature_list)
    #get sample lists:
    with open(SAMPLE_LIST_DIR+'PNP530_balancedAge_males') as fp:
        PNP530_balancedAge_males=pickle.load(fp)
    with open(SAMPLE_LIST_DIR+'Cardio126_balancedAge_males') as fp:
        Cardio126_balancedAge_males=pickle.load(fp)
    
    
    #get feature data
    feat_annot=pd.read_pickle('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Predictions/featureDFs/X_withPredictedAgeGender.dat')
   
#     featuresFD=pd.read_excel('%s/TCR_real_data/PNP530Cardio126Combined/featureSummaryDFs/PNP530Cardio126_featureDF_topFeatures.xlsx' %MyPath)
#     annotFD=pd.read_excel('%s/TCR_real_data/PNP530Cardio126Combined/annotationFeatures/PNP530Cardio126_PNP530Cardio126_annotationFeatures.xlsx' %MyPath)
#     feat_annot=pd.merge(featuresFD,annotFD,how='outer',left_index=True,right_index=True)
    feat_annot_PNP=feat_annot.loc[PNP530_balancedAge_males,:]
    feat_annot_cardio=feat_annot.loc[Cardio126_balancedAge_males,:]
    
    pnp_data_name='Healthy'
    cardio_data_name='Patients'
    
    newAxList=[]
    for n,feature in enumerate(feature_list):
        
        ax=axList[n]
        fig.add_subplot(ax)
        
        #get data for feature:
        pnp_data=feat_annot_PNP[feature].dropna().tolist()
        cardio_data=feat_annot_cardio[feature].dropna().tolist()        
        dataList=[(pnp_data_name,pnp_data),(cardio_data_name,cardio_data)]
        
        title=''
        
#         if n==0: 
#             showLegend=True
#             handles, labels = ax.get_legend_handles_labels()
#             ax.legend(labels, bbox_to_anchor=(-0.2,0.95))
#         else: showLegend=False

        ax,ks_p_cohort1_cohort2,t_p_cohort1_cohort2,p_Anov,filename=plotHistComprison(dataList,ax,
                                    title,showLegend=False,nBins=25,toAnnotate=False,alpha=None,plotType='hist',
                                    colorList=['grey','darkred'])
#         ax.set_xticks([])
#         ax.set_yticks([])

# 		ax.tick_params(axis='both', direction='in')
        ax.legend([],[])
        ax.set_title(feature_name_list[n],fontsize='medium',y=-15)
        ax.tick_params(axis='both',direction='in')
        
#         if n==0: asterisk='**'
#         else: asterisk='*'
#         ax.text(0.01,0.98,feature_name_list[n]+'(%s)' %asterisk,transform=ax.transAxes,ha='left',va='top',
#                fontsize='small',fontweight='bold')
        ymin,ymax=ax.get_ylim()
        ax.set_ylim(ymin,ymax*1.2)
        ax.set_xticks([]); ax.set_yticks([])
        
        newAxList.append(ax)
#         if n==0: 
#             handles, labels = ax.get_legend_handles_labels()
#             ax.legend(labels, fontsize='x-large',bbox_to_anchor=(-0.1,0.95))
#     fig.subplots_adjust(top=0.99,hspace=0.25,wspace=0.25) 

    
    return newAxList


####### RUN SCRIPT TO PLOT FIGURE 2#####

set_fig2_definitions()

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(18, 18))
gs0 = gridspec.GridSpec(3, 5, wspace=0.5, hspace=0.5)
gs1 = gridspec.GridSpec(3, 4, wspace=0.5, hspace=0.5)

##################################################################
ax1 = plt.Subplot(fig, gs0[0, 0:2])
fig.add_subplot(ax1)
ax1=plot_isCardio_multi_ROC(ax1)

##################################################################
ax2 = plt.Subplot(fig, gs0[0, 2:4])
fig.add_subplot(ax2)

ax2=plot_isCardio_multi_PR(ax2)

##################################################################
# ax3 = plt.Subplot(fig, gs1[1, 0])
ax3=fig.add_subplot(gs0[1, :2],projection='3d')
nToShow=20
ax3,topShapFeatures=shap_scatter_3d(ax3,nToShow=nToShow)

print topShapFeatures

##################################################################
ax4 = plt.Subplot(fig, gs0[1, 2:])
fig.add_subplot(ax4)

shap_isCardio=pd.read_pickle(PRED_RESULTS_DIR + 'isCardio/XGB_randomSearch_25_byRepFeatPCA10percVDJ0999PredictedAgeGender/shap_values.pkl')['isCardio']
shap_isCardio=shap_isCardio.loc[:,topShapFeatures]
# load feature values
Xfeatures_file = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/Predictions/featureDFs/TCRfeaturesNoT_PCA10percShared_withRels_VDJ_noCorr0-999_nanFilled_noConsts.dat'
Xfeatures = pd.read_pickle(Xfeatures_file).loc[PNP530, :]
Xfeatures=Xfeatures.loc[:,topShapFeatures]
for col in Xfeatures.columns:
    Xfeatures[col] = Xfeatures[col].fillna(Xfeatures[col].median())

nTopFeatures = nToShow
shap_df = shap_isCardio
features_df = Xfeatures
jitter = 0.1
sample_num = None
scalingMethod = 'perc'
target_name = 'isCardio'

ax4,fig = gen_shap_summary_to_axes(ax4, fig,target_name, nTopFeatures, shap_df, features_df, jitter=jitter, sample_num=sample_num, scalingMethod=scalingMethod)
#######################################################################
ax5 = plt.Subplot(fig, gs1[2, 0])
fig.add_subplot(ax5)
# ax5=plot_balanced_age_distrbution(ax5)
ax5=plot_pvalues_feature_comparison(ax5)

########################################################################
# ax6 = plt.Subplot(fig, gs1[2, 1])
# fig.add_subplot(ax6)
# ax6=plot_pvalues_feature_comparison(ax6)

########################################################################

gs00 = gridspec.GridSpecFromSubplotSpec(2, 4, subplot_spec=gs1[2,1:],wspace=0, hspace=0)
ax11=plt.Subplot(fig, gs00[:, :])
fig.add_subplot(ax11)
ax11.set_title('Specific Feature Comparisons')
ax11.set_ylabel('Frequency')
ax11.set_xticks([]); ax11.set_yticks([])

ax7=plt.Subplot(fig, gs00[0, 0])
ax8=plt.Subplot(fig, gs00[0, 1])
ax9=plt.Subplot(fig, gs00[0, 2])
ax10=plt.Subplot(fig, gs00[0, 3])
ax12=plt.Subplot(fig, gs00[1, 0])
ax13=plt.Subplot(fig, gs00[1, 1])
ax14=plt.Subplot(fig, gs00[1, 2])
ax15=plt.Subplot(fig, gs00[1, 3])

plot_feat_comp_kde(axList=[ax7,ax8,ax9,ax10,ax12,ax13,ax14,ax15])

########################################################################


## add letters:
ax1 = add_subfigure_letter(ax1, 'A', fontweight='bold')
ax2 = add_subfigure_letter(ax2, 'B', fontweight='bold')
ax3 = add_subfigure_letter(ax3, 'C', fontweight='bold')
ax4 = add_subfigure_letter(ax4, 'D', fontweight='bold')
ax5 = add_subfigure_letter(ax5, 'E', fontweight='bold')
# ax6 = add_subfigure_letter(ax6, 'F', fontweight='bold')
ax7 = add_subfigure_letter(ax7, 'F', fontweight='bold')


# ## save figure:
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.45, hspace=0.5)
fig.savefig(FIG2_DIR + 'figure2_%s_2.png' % cdate, dpi=300)

# # Save axes seperately
# for subp in [('ax1', ax1), ('ax2', ax2), ('ax3', ax3), ('ax4', ax4),
#                 ('ax5', ax5), ('ax6', ax6)]:
#     extent = subp[1].get_window_extent().transformed(fig.dpi_scale_trans.inverted())
#     fig.savefig(FIG2_DIR + '%s_figure_expanded_%s.png' % (subp[0], cdate), bbox_inches=extent.expanded(1.3, 1.35))

print 'Finished figure 2!!'

plt.show()


