################################################################################################
# all excel files used to plot this figure are generate by:
# calc_for_figure3.py (define FDR threshold there as needed)
# the fundemental data calculations were done using TCRMicrobiomeBinaryPhenCalculator.py
# using the following sample lists:
# PNP530_Age4570_Gender_Male11_HbA1C4464_BMI040
# Cardio126_Age3566_Gender_Male11_HbA1C4464_BMI040
################################################################################################

import cPickle as pickle
import pandas as pd
import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns
from os import makedirs
from os.path import isdir
import time
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from ShaniBA.CardioProject.Figures.GeneralFigureFunctions import remove_spines, edit_species_names

if __name__ == "__main__":
    # variable to change:
    FDR_t_species=0.1
    top100_phen_t_FDR = 0.2
    include_unknown_species=True

    CALC_RESULT_DIR_NAME = 'isCardio_min_shared_species01\
_max_shared_species1_min_shared_cluster025_max_shared_cluster1/'
    pvalue_log_ticks=[-4,0,4]
    font_size_dict = {
        'very_small':5,
        'small': 6,
        'medium': 7,
        'large': 12
    }
    # sample_list1_file = 'PNP530_Age4570_Gender_Male11_HbA1C4464_BMI040.pkl'
    # sample_list2_file = 'Cardio126_Age3566_Gender_Male11_HbA1C4464_BMI040.pkl'



    #definitions:
    GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
    DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
    SAMPLE_LIST_DIR = GENIE_BASE_DIR + 'Sample files/BD lists/'
    CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'
    FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig3_balanced_comparison/calc_fig3/'
    CALC_RESULT_DIR = FIGURE_DIR + CALC_RESULT_DIR_NAME
    EXCEL_FILES_DIR = FIGURE_DIR + 'excel_files_new/'
    PLOT_DIR = FIGURE_DIR + 'plots/'
    if not isdir(PLOT_DIR): makedirs(PLOT_DIR)

    cdate = str(time.strftime("%d%m%Y"))
    if include_unknown_species: suf = '_inc_unknown_species'
    else: suf=''

    ##########################################################
    # get sample lists - only samples that have swab data

    ## healthy (n=45):
    with open(CALC_RESULT_DIR + 'sample_list1_swab.pkl', 'rb') as fp1:
        sample_list1 = pickle.load(fp1)
    ## patients (n=48):
    with open(CALC_RESULT_DIR + 'sample_list2_swab.pkl', 'rb') as fp2:
        sample_list2 = pickle.load(fp2)

    ############################################################

#functions:

def get_species(x):
        if isinstance(x, str):
            return x.split('|')[-6:-3]


# define figure parameters
def set_fig3_definitions(size_factor=1):
    """
    this function sets matplotlib parameters to be constant across the whole figure.
    :param size_factor: a factor to multiply all parameters in order to efficiently increase/decrese the figures
    while keeping the correct size relations.
    :return:
    """
    figsize=(4.75 *size_factor, 6.5*size_factor)
    global figsize
    #TODO: arrange this funciton

    changing_size_params = {
    'axes.labelsize': 5,
    'font.size': 5,
    'legend.fontsize': 5,
    'axes.titlesize':7,
    'xtick.labelsize': 4,
    'ytick.labelsize': 3.2,
    'xtick.major.size': 1.75,
    'ytick.major.size': 1.75,
    'axes.linewidth': 0.5,
    'lines.linewidth': 0.8,
    'lines.markersize': 4
    }

    other_params = {
    'axes.titleweight': 'bold',
    'text.usetex': False,
    'xtick.direction': 'out',
    'xtick.major.pad': 1,
    'ytick.major.pad': 1,
    'axes.edgecolor': 'black',
    'axes.facecolor':'white',
    'figure.dpi': 300,
    # 'axes.labelpad': 0.8,
    'legend.labelspacing': 0.4,
    # 'legend.edgecolor': '0',
    'legend.frameon': False,
        'axes.xmargin': 0.1,
        'axes.ymargin': 0.1,
    }

    if size_factor != 1:
        for k,v in changing_size_params.items():
            changing_size_params[k] = v*size_factor


    mpl.rcParams.update(changing_size_params)
    mpl.rcParams.update(other_params)

    return

def get_col_order(df,order_func = np.sum):
    healthy_cols = [col for col in df.columns if int(col.replace('BD','')) <950]
    patient_cols = [col for col in df.columns if int(col.replace('BD', '')) > 949]

    df_healthy = df.loc[:,healthy_cols]
    df_patients = df.loc[:,patient_cols]

    df_healthy_order = df_healthy.apply(lambda x: order_func(x), axis=0).sort_values(ascending=False).index.tolist()
    df_patients_order = df_patients.apply(lambda x: order_func(x), axis=0).sort_values().index.tolist()

    col_order = df_healthy_order + df_patients_order

    return col_order





def cluster_sample_heatmap(hm_ax,cb_ax):
    cluster_sample_df = pd.read_excel(EXCEL_FILES_DIR + 'cluster_df_for_fig.xlsx')

    #order cols:
    col_order = get_col_order(cluster_sample_df,order_func = np.sum)
    cluster_sample_df = cluster_sample_df.loc[:,col_order]

    #define min, max_vals:
    max_val = int(cluster_sample_df.max().max())
    min_val = 0
    dif = max_val - min_val

    #plot heatmap and colorbar
    cmap_2 = plt.get_cmap('Reds', dif + 1)
    hm2 = sns.heatmap(cluster_sample_df, cmap=cmap_2, ax=hm_ax, cbar=False, linewidth=0.05, vmin=min_val, vmax=max_val,
                      yticklabels=True, xticklabels=False)
    cbar_2 = fig.colorbar(hm_ax.collections[0], ax=cb_ax, use_gridspec=False, 
                          fraction=0.5, aspect=7, orientation = 'horizontal', anchor = (0.5, 0.2))
    # orientation = 'horizontal', anchor = (0.25, 0.5),

    #edit heatmap:
    hm_ax.text(0,1,'# present sequences per cluster',
               transform=hm_ax.transAxes, ha='left',va='bottom')
    hm_ax.set_ylabel('')
    for spine in hm_ax.spines.values():
        spine.set_visible('on')
    hm_ax.tick_params(axis='y', which='both', right='off',left='off')
    hm_ax.set_yticklabels(cluster_sample_df.index)
    # hm_ax.set_yticklabels(hm_ax.get_yticklabels())
    hm_ax.axvline(x=len(sample_list1), color='red')

    # add titles:
    x_patients = (len(sample_list1) + float(len(sample_list2)) / 2) / (len(sample_list1) + len(sample_list2))
    x_healthy = (float(len(sample_list1)) / 2) / (len(sample_list1) + len(sample_list2))
    hm_ax.text(x_patients, 1.03, 'Patients', ha='center', transform=hm_ax.transAxes, va='bottom')
    hm_ax.text(x_healthy, 1.03, 'Healthy', ha='center', transform=hm_ax.transAxes, va='bottom')

    #edit cbar:
    cbar_2.set_ticks([float(dif) / 12 + float(dif) * x / (dif + 1) for x in range(dif + 1)])
    cbar_2.set_ticklabels(range(dif + 1))
    # cbar_2.ax.set_xticklabels(cbar_2.ax.get_xticklabels())
    cb_ax = remove_spines(cb_ax,removeFigBorders=True)
    cbar_2.ax = remove_spines(cbar_2.ax,removeTicklabels=False)
    cb_ax.text(0.5, -3.3, '# present sequences per cluster', transform=cb_ax.transAxes,
               ha='center', va='top')

    return hm_ax,cb_ax

def cluster_barplot(bp_ax):
    cluster_disease_df = pd.read_excel(EXCEL_FILES_DIR + 'top100_clusters_with_annot.xlsx').\
        set_index('cluster')
    # cluster_disease_df['abs'] = cluster_disease_df['p_to_show_cluster_isCardio'].abs()
    # cluster_disease_df = cluster_disease_df.sort_values(by='abs',ascending=True)

    #plot:
    bp1 = cluster_disease_df['p_to_show_cluster_isCardio'].plot(kind='barh',color='black',
                                                                width=0.3, ax=bp_ax)

    #edit plot:
    bp_ax.invert_yaxis() #reverse cluster order to match other subplots
    bp_ax.text(0, 1, 'cluster-disease\nassociations',
               transform=bp_ax.transAxes, ha='left', va='bottom')
    bp_ax.tick_params(axis='y', which='both', left=False, right=False)
    bp_ax.tick_params(axis='x', which='both', top=False)

    bp_ax.set_yticklabels('')
    bp_ax.set_ylabel(''); bp_ax.set_xlabel('processed\nlog p-value')
    bp_ax.set_xticks(pvalue_log_ticks)
    bp_ax.set_xticklabels(pvalue_log_ticks, rotation=90)
    bp_ax.spines['top'].set_visible(False)
    bp_ax.spines['right'].set_visible(False)
    bp_ax.spines['left'].set_position(('data', 0))

    return bp_ax

def y_pred_real_heatmap(hm_ax,cb_ax):

    y_pred_real = pd.read_excel(EXCEL_FILES_DIR + 'Y_pred_real.xlsx')
    y_pred_real.index = ['Real', 'Pred']

    #plot heatmap and colorbar:
    hm1 = sns.heatmap(y_pred_real, cmap="Greys", ax=hm_ax, cbar=False, linewidth=0.05)
    cbar_1 = fig.colorbar(hm_ax.collections[0], ax=cb_ax, use_gridspec=False, location='bottom',
                              pad=0.08, fraction=0.5, aspect=7)

    #edit heatmap:
    hm_ax = remove_spines(hm_ax, removeTicklabels=False,removeFigBorders=False)
    hm_ax.set_xticklabels('')
    hm_ax.set_yticklabels(hm1.get_yticklabels(), rotation = 0)

    hm_ax.axvline(x=len(sample_list1), color='red')
    x_patients = (len(sample_list1) + float(len(sample_list2))/ 2 ) / (len(sample_list1)+len(sample_list2))
    x_healthy = (float(len(sample_list1)) / 2) / (len(sample_list1)+len(sample_list2))
    hm_ax.text(x_patients, 1.15, 'Patients', ha='center', transform=hm_ax.transAxes,va='bottom')
    hm_ax.text(x_healthy, 1.15, 'Healthy', ha='center', transform=hm_ax.transAxes,va='bottom')


    #edit colorbar:
    cbar_1.set_ticks([0, 0.5, 1])
    cb_ax = remove_spines(cb_ax,removeFigBorders=True)
    cb_ax.text(0.5, 1.75, 'Healthy (0) / Patient (1)', transform=cb_ax.transAxes, ha='center', va='top')

    return hm_ax, cb_ax




def species_heatmap(hm_ax,cb_ax):

    cluster_species_pivot = pd.read_excel(EXCEL_FILES_DIR + 'cluster_species_hm_fdr%s%s.xlsx' %(str(FDR_t_species).\
                                          replace('.',''),suf)).set_index('cluster')
    new_species = edit_species_names(cluster_species_pivot) #edit species names for plot
    cluster_species_pivot.columns = new_species

    #plot heatmap and colorbar:
    hm3 = sns.heatmap(cluster_species_pivot, cmap="seismic", ax=hm_ax, cbar=False, linewidth=0.05,vmin=-5,vmax=5,
                      yticklabels=False, xticklabels=False)
    cbar_3 = fig.colorbar(hm_ax.collections[0], ax=cb_ax, use_gridspec=False,
                               fraction=0.5, aspect=7, orientation = 'horizontal', anchor = (0.5, 0.2))
    # orientation = 'horizontal', anchor = (0.25, 0.5)
    # pad = 0.08 location='bottom',
    #edit heatmap:
    hm_ax = remove_spines(hm_ax, removeTicklabels=False)
    for spine in hm_ax.spines.values():
        spine.set_visible('on')
    hm_ax.text(0, 1, 'cluster-species associations', transform=hm_ax.transAxes,
               ha='left', va='bottom')
    hm_ax.set_xticklabels(''); hm_ax.set_yticklabels('')
    hm_ax.set_ylabel('')

    #edit colorbar:
    cbar_3.set_ticks(pvalue_log_ticks)
    # cbar_3.ax.set_xticklabels(cbar_3.ax.get_xticklabels())
    cb_ax = remove_spines(cb_ax,removeFigBorders=True)
    cbar_3.ax = remove_spines(cbar_3.ax, removeTicklabels=False)
    cb_ax.text(0.5, -3.3, 'cluster-species associations\nlog10 p-values',
               transform=cb_ax.transAxes, ha='center', va='top')

    return hm_ax, cb_ax

def species_barplot(bp_ax):
    species_disease_df = pd.read_excel(EXCEL_FILES_DIR + 'species_disease_logp_fdr%s%s.xlsx'\
                                       %(str(FDR_t_species).replace('.',''),suf)).set_index('species2')
    # species_disease_df['abs'] = species_disease_df['p_to_show_species_phen'].abs()
    # species_disease_df = species_disease_df.sort_values(by='abs', ascending=True)
    species_name_list = species_disease_df.index.tolist()
    species = edit_species_names(species_name_list)
    species = [x.replace('\n', '') for x in species]

    #plot:
    species_disease_df['p_to_show_species_phen'].plot(kind='bar', ax=bp_ax, legend=False, color='black',width=0.4)
    bp_ax.invert_xaxis()

    #edit:
    bp_ax.set_ylabel('processed\nlog p-value')
    bp_ax.spines['top'].set_visible(False)
    bp_ax.spines['right'].set_visible(False)
    bp_ax.spines['bottom'].set_position(('data', 0))
    bp_ax.xaxis.set_ticks_position('top')
    bp_ax.set_xticklabels(species, rotation=90)
    bp_ax.set_yticks(pvalue_log_ticks)
    bp_ax.tick_params(axis='y', which='major', pad=2)
    bp_ax.tick_params(top='off', bottom='off', left='off', right='off', labelleft='on')

    bp_ax.set_yticklabels(pvalue_log_ticks)
    bp_ax.set_xlabel('')
    # bp_ax = remove_spines(bp_ax, removeFigBorders=True)

    return


def phen_heatmap(hm_ax):

    cluster_phen_pivot = pd.read_excel(EXCEL_FILES_DIR + 'top100_phen_df_fdr%s%s.xlsx' %(str(top100_phen_t_FDR).\
                                       replace('.',''),suf)).drop('Obese_logp', axis=1) #drop obese because it's all nans
    cluster_phen_pivot.columns = ['Age','Gender','HbA1C']

    #plot heatmap:
    hm4 = sns.heatmap(cluster_phen_pivot, cmap="PiYG", ax=hm_ax, cbar=False, linewidth=0.05,vmin=-5,vmax=5,
                      yticklabels=False, xticklabels=True)

    #edit heatmap:
    hm_ax.text(0, 1, 'cluster-phenotype associations',
               transform=hm_ax.transAxes, ha='left', va='bottom')
    hm_ax = remove_spines(hm_ax, removeTicklabels=False)
    for spine in hm_ax.spines.values():
        spine.set_visible('on')
    hm4.xaxis.set_ticks_position('top')
    hm4.set_xticklabels(['Age','Gender','HbA1C'], rotation=90)
    # hm4.tick_params(labeltop='on')
    hm_ax.set_yticklabels('')

    return hm_ax

############################################################################

# plot:

if __name__ == "__main__":
    set_fig3_definitions(size_factor=1)


    # define figure structure
    fig = plt.figure(figsize=figsize, dpi=300)

    ##add sub-figure letters and remove spines:
    ax = plt.gca()
    plt.text(0.005, 0.95, 'A', ha='left', va='top', transform=ax.transAxes, fontsize='xx-large', fontweight='bold')
    plt.text(0.55, 0.95, 'B', ha='left', va='top', transform=ax.transAxes, fontsize='xx-large', fontweight='bold')
    # plt.text(0.895,1.03,'C',ha='left',va='bottom',transform=ax.transAxes,fontsize='xx-large',fontweight='bold')
    remove_spines(removeFigBorders=True)

    gs3 = gridspec.GridSpec(1, 1)
    gs3.update(top=0.75, bottom=0.72,left=0.61, right=0.98)
    ax3 = fig.add_subplot(gs3[0, 0])

    gs4 = gridspec.GridSpec(1, 1)
    gs4.update(top=0.7, bottom=0.09, left=0.03, right=0.06)
    ax4 = fig.add_subplot(gs4[0, 0])

    gs5 = gridspec.GridSpec(1, 1)
    gs5.update(top=0.7, bottom=0.09, left=0.21, right=0.58)
    ax5 = fig.add_subplot(gs5[0, 0])

    gs6 = gridspec.GridSpec(1, 1)
    gs6.update(top=0.7, bottom=0.09, left=0.61, right=0.98)
    ax6 = fig.add_subplot(gs6[0, 0])

    gs8 = gridspec.GridSpec(1, 1)
    gs8.update(top=0.08, bottom=0.05, left=0.21, right=0.58)
    ax8 = fig.add_subplot(gs8[0, 0])

    gs9 = gridspec.GridSpec(1, 1)
    gs9.update(top=0.08, bottom=0.05, left=0.61, right=0.98)
    ax9 = fig.add_subplot(gs9[0, 0])


    ax5,ax8 = cluster_sample_heatmap(hm_ax=ax5,cb_ax=ax8)
    # ax1, ax7 = y_pred_real_heatmap(hm_ax=ax1,cb_ax=ax7)
    ax6, ax9 = species_heatmap(hm_ax=ax6,cb_ax=ax9)
    # ax6 = phen_heatmap(hm_ax=ax6)
    ax3 = species_barplot(ax3)
    ax4 = cluster_barplot(ax4)


    # plt.subplots_adjust(left=0.01,right=0.99,bottom=0.08,top=0.9,wspace=0.7,hspace=0.7)
    # dpi=150
    plt.savefig(PLOT_DIR + 'figure3_ver4_%s.png' % cdate, dpi=300)
    print ('need to save fig from console to keep figure proportions!')
