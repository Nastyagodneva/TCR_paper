import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from os import makedirs
from os.path import isdir
import time
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from ShaniBA.CardioProject.paper_summer2019.fig5.ver2.plot_functions_for_fig5 import *
from ShaniBA.CardioProject.Figures.GeneralFigureFunctions import remove_spines
# from ShaniBA.CardioProject.Figures.GeneralFigureFunctions import remove_spines, edit_species_names

cdate = str(time.strftime("%d%m%Y"))


#directories:
GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/'
FIGURE_FILES_DIR = FIGURE_DIR + 'files_for_figure/'



def get_data():
    """
    this function extracts all dataframe required for the plot and declare them as global variables
    all the files for these dataframes were copied to a designated directory (FIGURE_FILES_DIR) so they
    need to be manually updated if updated in their original location
    """
    sup_table_4 = pd.read_excel(FIGURE_FILES_DIR+'sup table 4.xlsx')
    or_result_table = pd.read_excel(FIGURE_FILES_DIR+ 'OR_top_roc_auc_rehosp05.xlsx')
    new_outcome_df = pd.read_excel(FIGURE_FILES_DIR + 'new_outcome_df.xlsx').set_index('BD')
    # cluster_data = pd.read_pickle(DATA_DIR + 'TCR_seqs_clusters/newX_onlySeqs_025_085_noNans.dat')
    cluster_data = pd.read_excel(FIGURE_FILES_DIR + 'short_cluster_table.xlsx')
    global sup_table_4, or_result_table,new_outcome_df, cluster_data
    return


def set_fig5_definitions(size_factor=1):
    """
    this function sets matplotlib parameters to be constant across the whole figure.
    :param size_factor: a factor to multiply all parameters in order to efficiently increase/decrese the figures
    while keeping the correct size relations.
    :return:
    """
    figsize=(2.6*size_factor, 4.5*size_factor)
    global figsize
    #TODO: arrange this funciton

    changing_size_params = {
    'axes.labelsize': 7,
    'font.size': 6,
    'legend.fontsize': 7,
    'axes.titlesize':7,
    'xtick.labelsize': 5.5,
    'ytick.labelsize': 5.5,
    'xtick.major.size': 1.75,
    'ytick.major.size': 1.75,
    'axes.linewidth': 0.5,
    'lines.linewidth': 1.2,
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


def plot_volcano(ax,data_df,fdr=0.1,color_sig='orange',color_non_sig='black',
                 top5_color='red',marker_size=8, alpha=0.6):
    """
    this function generates the first part of the figure, which is a volcano plot
    of auc-roc values vs. their p-values for all clusters in data_df. colored by the significance (non-significant,
    significant with FDR<FDR parameter, and the top 5 clusters)
    :param ax: subplot axes
    :param data_df: a dataframe containing the clusters to plot with their auc roc values, p-values and FDR values, typiclly
    it is the sup_table_4
    :param fdr: fsr threshold to color clusters.
    :param color_sig: color for significant clusters
    :param color_non_sig: color for non-significcant clusters
    :param top5_color: color for the top 5 clusters
    :param marker_size:
    :param alpha:
    :return:
    """

    #divide data to top5, significant by fdr, and non-significant
    data_df['-log_p_value'] = -np.log10(data_df['p_value (200,000 permutations)'].replace(0,0.00000000001))
    top_5 = data_df.iloc[:5]
    top5_clusters = top_5['cluster'].tolist()
    global top5_clusters
    other_data = data_df.iloc[5:]
    sig_cond = (other_data['FDR corrected p-value'] < fdr) & (other_data['roc_auc'] >0.7)
    sig_data = other_data[sig_cond]
    non_sig_data = other_data[~sig_cond]

    #define color and lable for each ata group
    data_list = [
        (top_5, top5_color, 'top_5_clusters'),
        (sig_data,color_sig,'FDR<%s' %fdr),
        (non_sig_data,color_non_sig,'n.s.'),
        ]

    #plot all data with corresponding color and label:
    for item in data_list:
        ax.scatter(item[0]['roc_auc'], item[0]['-log_p_value'],
                   color=item[1], s=marker_size, label=item[2],
                   alpha=alpha, edgecolors='none')

    #edit plot:
    ax.set_xlabel('AUC (roc curve)')
    ax.set_ylabel('p-value')
    ax.legend()
    return ax

def plot_or_fig5(ax, data_df):
    """
    this funciton plots the left panel of the second part of the figure,
    which is the OR confidence interval plot for the top 5 clsuters.
    :param ax:
    :param data_df: typically or_result_table
    :return:
    """
    ax = plot_or(or_result_table=data_df, cluster_list=top5_clusters,
                 cluster_num_thresh=0, ax=ax)
    return ax

def plot_multiple_roc(ax_list):
    """
    this funciton plots the right panel of the second part of the figure,
    which is a combination of 5 distinct roc curve plots
    :param ax_list: list of 5 axes to hold the roc-auc curves
    :return:
    """

    outcome = 'CV hospitalization including chest pain'
    subplot_list = zip(ax_list, top5_clusters)
    for n,s in enumerate(subplot_list):
        ax = s[0]
        cluster = s[1]
        ax.yaxis.tick_right()
        fig_obj = Fig5SubplotBoxRocplot(ax=ax, cluster=cluster, outcome=outcome,
                                        cluster_data=cluster_data,outcome_df=new_outcome_df)
        fig_obj.gen_roc_plot()
        ax.set_ylabel('')
        ax.set_yticks([0.5,1])
        if n!=4:
            ax.set_xlabel('')
            ax.set_xticks([])
        else:
            ax.set_xticks([0.5, 1])
    return



if __name__ == '__main__':
    get_data() #extracts all dataframes required for the figure
    set_fig5_definitions(size_factor=1.4) #sets mpl.rc_params: size, color and shape definitions for differnt
                                            # aspects of the figure
    fig = plt.figure(figsize=figsize)

    # generate the general outline of the figure
    gs01 = gridspec.GridSpec(1,1) # figure upper part
    gs01.update(top=0.98, bottom=0.68,left=0.15, right=0.85)
    gs02 = gridspec.GridSpec(5,5,wspace=0, hspace=0) #figure lower part
    gs02.update(top=0.58, bottom=0.05, left=0.35,right=0.85)

    #generate all the subplots within the outline:
    ax1 = fig.add_subplot(gs01[0,0])
    ax2 = fig.add_subplot(gs02[:,:3])

    # add a virtual subplot to add a text for all 5 subplots' ylabel:
    ax3 = fig.add_subplot(gs02[:,3:])
    ax3.text(1.35, 0.5, 'TPR', ha='left', va='center', rotation=90,
             transform=ax3.transAxes, fontsize=mpl.rc_params()['axes.labelsize'])
    remove_spines(ax=ax3, removeFigBorders=True)

    ax3_0 = fig.add_subplot(gs02[0,3:])
    ax3_1 = fig.add_subplot(gs02[1, 3:])
    ax3_2 = fig.add_subplot(gs02[2, 3:])
    ax3_3 = fig.add_subplot(gs02[3, 3:])
    ax3_4 = fig.add_subplot(gs02[4, 3:])

    # plot the figure:
    ax1 = plot_volcano(ax=ax1,data_df=sup_table_4)
    ax2 = plot_or_fig5(ax=ax2, data_df=or_result_table)
    plot_multiple_roc(ax_list=[ax3_0, ax3_1, ax3_2, ax3_3, ax3_4])

    #save figure:
    fig.savefig(FIGURE_DIR + 'figure5_python_1_25_changed_size_%s_trial.png' %cdate)

