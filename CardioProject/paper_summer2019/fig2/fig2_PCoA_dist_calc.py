from __future__ import absolute_import, division, print_function

import pandas as pd
import numpy as np
import cPickle as pickle
import matplotlib as mpl
import matplotlib.pyplot as plt


from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import genDistMat,calc_PCoA
from ShaniBA.MyFunctionsShani import plotHistComprison
#############################################################################
#definitions:
#definitions:
GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
SAMPLE_LIST_DIR = GENIE_BASE_DIR + 'Sample files/BD lists/'
AllUniqueWithCounts_path='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/sharingAnalysis/\
AllUniqueWithCounts'
FIGURE2_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/June2019/calc_fig2/'

CLUSTER_FILE = GENIE_BASE_DIR + 'TCR_real_data/PNP530Cardio126Combined/sharingAnalysis/\
sharingMatrix_PNP530Cardio126_minNshared5_RA_onlyProductiveTrue'
CLUSTER_FILE_NAME = 'aa_seqs'
# CLUSTER_FILE = DATA_DIR + 'TCR_seqs_clusters/newX_onlySeqs_025_085_noNans.dat
# CLUSTER_FILE_NAME = 'clusters_025_085'
seed = 1

# get sample lists
sample_list1_name = 'PNP530_Age4570_Gender_Male11_HbA1C4464_BMI040'
sample_list1_file = 'PNP530_Age4570_Gender_Male11_HbA1C4464_BMI040.pkl'
sample_list2_name = 'Cardio126_Age3566_Gender_Male11_HbA1C4464_BMI040'
sample_list2_file =  'Cardio126_Age3566_Gender_Male11_HbA1C4464_BMI040.pkl'

with open(SAMPLE_LIST_DIR + sample_list1_file,'rb') as fp1: sample_list1 = pickle.load(fp1)
with open(SAMPLE_LIST_DIR + sample_list2_file,'rb') as fp2: sample_list2 = pickle.load(fp2)


############################################################################
# functions:

def fig2_ver3_plotPairWiseDistance(df_condensed_org, ax):
    df_condensed_org['sample1_isPatient'] = np.where(df_condensed_org['sample1'].isin(sample_list2), 1, 0)
    df_condensed_org['sample2_isPatient'] = np.where(df_condensed_org['sample2'].isin(sample_list2), 1, 0)
    df_condensed_org['samples_isPatient_sum'] = df_condensed_org[['sample1_isPatient', 'sample2_isPatient']].sum(axis=1)
    df_condensed_org['samples_isPatient_mapped'] = df_condensed_org['samples_isPatient_sum'].map({0: 'Both healthy',
                                                                                                  1: '1 patient',
                                                                                                  2: 'Both patients'})

    dataList = [
        ('Both healthy', df_condensed_org[df_condensed_org['samples_isPatient_mapped'] == 'Both healthy']['dist']),
        ('Both patients', df_condensed_org[df_condensed_org['samples_isPatient_mapped'] == 'Both patients']['dist'])]

    title = ''
    colorList = ['grey', 'darkred']
    alpha = 1
    text_kws = {'fontsize': 'large', 'fontweight': 'bold', 'color': 'red',
                'horizontalalignment': 'left', 'verticalalignment': 'top'}
    plotType = 'kde'

    ax, ks_p_cohort1_cohort2, t_p_cohort1_cohort2, p_Anov, filename,meanText = plotHistComprison(dataList, ax, title,
                                                                                        showLegend=False, nBins=20,
                                                                                        toAnnotate=False,
                                                                                        colorList=colorList,
                                                                                        alpha=alpha,
                                                                                        text_kws=text_kws,
                                                                                        plotType=plotType)
    ax.text(0.02, 0.95, 'p=%.0E (t-test)' % t_p_cohort1_cohort2, transform=ax.transAxes, ha='left', va='top',
            fontsize='medium')

    ax.set_ylabel('Density');
    ax.set_xlabel('Pair-wise sample distances')
    handles, labels = ax.get_legend_handles_labels()
    labels = [label.replace('Both ', '') for label in labels]
    ax.legend(loc='lower left', labels=labels, fontsize=14)

    return ax


def fig2_ver3_plot_PCoA_balancedMales(ax_PCoA, ax_PCoA_dist):
    print ('reading cluster matrix...')
    sampleByClusterDF = pd.read_pickle('/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/\
sharingAnalysis/seqClusters_allProd_maxdist1/sampleByClusterDF_cohortfiltering005-085perc_dropped.dat')
    print ('binarizing matrix...')
    sampleByClusterDF_binary = (sampleByClusterDF > 0).astype(int)
    # get samples lists:
    with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/sample_list1') as fp:
        sample_list1 = pickle.load(fp)
    with open(
            '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/sample_list2') as fp:
        sample_list2 = pickle.load(fp)
    sampleByClusterDF_binary_balanced = sampleByClusterDF_binary.loc[
                                        sample_list1 + sample_list2, :].fillna(0)
    print ('sampleByClusterDF_binary_balanced shape: ', sampleByClusterDF_binary_balanced.shape)

    # calculate PCoA and plot:
    df = sampleByClusterDF_binary_balanced
    metric = 'jaccard'
    sample_list_list = [('Healthy', sample_list1), ('Patients', sample_list2)]
    color_list = ['grey', 'darkred']
    pcoa_n1_toplot = 0
    pcoa_n2_toplot = 1
    toScale = False,
    toAnnotate = False

    pcoa_df, fig, ax, df_condensed_org = calc_PCoA(df, metric, sample_list_list, color_list, pcoa_n1_toplot=0,
                                                   pcoa_n2_toplot=1, toScale=False,
                                                   toAnnotate=False, ax=ax_PCoA, calculateSeperation=False)
    replacement_list = [("Text(0.5,0,u'", ""), ("Text(0,0.5,u'", ""), ("')", ""), ("PC0", "PC1")]
    Xlabel = ax_PCoA.xaxis.get_label()
    Ylabel = ax_PCoA.yaxis.get_label()
    print ('original Xlabel is ', Xlabel)
    print ('original Ylabel is ', Ylabel)
    for item in replacement_list:
        Xlabel = str(Xlabel).replace(item[0], item[1])
        Ylabel = str(Ylabel).replace(item[0], item[1])
        Ylabel = str(Ylabel).replace("PC1", "PC2")

    dist_pcoa_toplot = 0

    dataList = [('Healthy', pcoa_df.loc[sample_list1, pcoa_n1_toplot]),
                ('Patients', pcoa_df.loc[sample_list2, pcoa_n1_toplot])]
    title = ''
    colorList = ['grey', 'darkred']
    alpha = 1
    text_kws = {'fontsize': 'large', 'fontweight': 'bold', 'color': 'red',
                'horizontalalignment': 'left', 'verticalalignment': 'top'}
    plotType = 'kde'

    #     ax_PCoA_dist,ks_p_cohort1_cohort2,t_p_cohort1_cohort2,p_Anov,filename=plotHistComprison(dataList,ax_PCoA_dist,
    #                     title,showLegend=False,nBins=20,toAnnotate=False,colorList=colorList,alpha=alpha,
    #                           text_kws=text_kws,plotType=plotType)

    ax_PCoA_dist = fig2_ver3_plotPairWiseDistance(df_condensed_org, ax)

    # some more edits:
    ax_PCoA_dist.set_xticklabels('')
    ax_PCoA_dist.set_title('')
    ax_PCoA_dist.set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
    #     ax_PCoA.yaxis.set_major_locator(plt.MaxNLocator(5))
    t_p_cohort1_cohort2 = 0.01
    ax_PCoA.text(0.02, 0.95, 'p=%s (PCo%s,t-test)' % (round(t_p_cohort1_cohort2, 2), dist_pcoa_toplot + 1),
                 horizontalalignment='left',
                 verticalalignment='top', transform=ax_PCoA.transAxes, fontsize='x-large')
    #     fontsize=mpl.rcParams['axes.labelsize']-4
    ax_PCoA.set_ylabel(Ylabel, fontsize=mpl.rcParams['axes.labelsize'])
    handles, labels = ax_PCoA.get_legend_handles_labels()
    ax_PCoA.legend(loc='lower left', labels=labels, fontsize=14)

    ax_PCoA.set_xlabel(Xlabel)
    ax_PCoA_dist.set_ylabel('Density', labelpad=25, fontsize=mpl.rcParams['axes.labelsize'])
    ax_PCoA_dist.set_yticks([0, 5, 10])

    # fig.align_ylabels(axes=[ax_PCoA,ax_PCoA_dist])
    ax_PCoA.yaxis.set_label_coords(-0.1, 0.5)
    ax_PCoA_dist.yaxis.set_label_coords(-0.1, 0.5)

    return ax_PCoA, ax_PCoA_dist


def fig2_ver3_plot_PCoA_balancedMales_new(ax_PCoA_dist,sample_list1,sample_list2):
    print ('reading cluster matrix...')
    sampleByClusterDF = pd.read_pickle(CLUSTER_FILE)
    print ('binarizing matrix...')
    sampleByClusterDF_binary = (sampleByClusterDF > 0).astype(int)
    # get samples lists:

    sampleByClusterDF_binary_balanced = sampleByClusterDF_binary.loc[
                                        sample_list1 + sample_list2, :].fillna(0)
    print ('sampleByClusterDF_binary_balanced shape: ', sampleByClusterDF_binary_balanced.shape)

    # generate df_condensed_org:
    ## sample_list_list = list of tuples, each tuple is composed of a string which is the name of the sample list, and the sample list
    from skbio.stats.ordination import PCoA
    from skbio.stats.distance import DistanceMatrix

    df = sampleByClusterDF_binary_balanced
    metric = 'jaccard'
    sample_list_list = [('Healthy', sample_list1), ('Patients', sample_list2)]
    color_list = ['grey', 'darkred']
    toScale = False,
    toAnnotate = False

    ###generate distance matrix:
    print ('generating distance matrix...')
    df_condensed_org, distMat_square = genDistMat(df, metric, generateSquareform=False)

    ax_PCoA_dist = fig2_ver3_plotPairWiseDistance(df_condensed_org, ax_PCoA_dist)



    #     ax_PCoA.set_xlabel(Xlabel)
    ax_PCoA_dist.set_ylabel('Density', labelpad=25, fontsize=mpl.rcParams['axes.labelsize'])
    ax_PCoA_dist.set_yticks([0, 5, 10])

    # fig.align_ylabels(axes=[ax_PCoA,ax_PCoA_dist])
    #     ax_PCoA.yaxis.set_label_coords(-0.1, 0.5)
    ax_PCoA_dist.yaxis.set_label_coords(-0.1, 0.5)

    return ax_PCoA_dist


fig,ax = plt.subplots(figsize=(6, 4))

fig2_ver3_plot_PCoA_balancedMales_new(ax, sample_list1,sample_list2)

fig.savefig(FIGURE2_DIR + 'fig2_ver3_plot_PCoA_balancedMales_new_%s.png' %CLUSTER_FILE_NAME)

plt.show()


