from __future__ import absolute_import, division, print_function


import datetime
from os import listdir, makedirs
from os.path import isdir
import os
from time import gmtime, strftime
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import numpy as np
#### imports: #######
import pandas as pd
from GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import analyze_associations_binary_nonbinary
from TCR_feature_generation.SequenceClusteringModule import get__limited_sharing, get_annot_for_seqs_in_cluster
from CardioProject.Figures.GeneralFigureFunctions import make_custome_legend
###### paths: ########

GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
SAMPLE_LIST_DIR = GENIE_BASE_DIR + 'Sample files/BD lists/'
CLUSTER_ANALYSIS_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/sharingAnalysis/seqClusters_allProd_maxdist1/'
CARDIO_DIR = GENIE_BASE_DIR + 'TCR_real_data/CardioSamples/'
CARDIO_PHEN_DIR = CARDIO_DIR + 'phenotypicData/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5_updated_phens_repeatrun/'

start_time = datetime.datetime.now()
#### VARIABLES TO CHANGE: #######

#phenotype definition:
# phen_df should be a binary dataframe with 'BD' as index
phen_name = 'CV hospitalization including chest pain'
phen_col = 'CV hospitalization including chest pain'
phen_file = CARDIO_PHEN_DIR + 'new_outcome_df.xlsx'
phen_df = pd.DataFrame(pd.read_excel(phen_file).set_index('BD')[phen_col])
TOP_RESULT_DIR = FIGURE_DIR

#define do's and don't do's:
calc_species_phen = True
calc_cluster_phen = True
calc_cluster_species = True
merge_all_results = True
plot_pvalues = True
plot_top_pairs = True
rewrite = False

#define sharing threshold (for choosing which clusters and species to include in analysis):

min_shared_species = 0.1
max_shared_species = 1
min_shared_cluster = 0.3
max_shared_cluster = 1

#define thresholds for species and clusters association with the phen, for selecting them to test their
# association
#choose either FDR or nTOP threshold for each!!!
FDR_lim_species_phen = 0.25
nTOP_species_phen = None
FDR_lim_cluster_phen = None
nTOP_cluster_phen = 100
FDR_lim_cluster_species = 0.1
nTOP_cluster_species = None
FDR_for_plotting = 0.25

sample_list1_name = 'Cardio126_CV hospitalization including chest pain00'
sample_list1_file = 'Cardio126_CV hospitalization including chest pain00.pkl'
sample_list2_name = 'Cardio126_CV hospitalization including chest pain11'
sample_list2_file = 'Cardio126_CV hospitalization including chest pain11.pkl'

##### CONSTANT VARIABLES: #####

sample_list_name1 = sample_list1_name
sample_list_name2 = sample_list2_name

with open(SAMPLE_LIST_DIR + sample_list1_file,'rb') as fp1: sample_list1 = pickle.load(fp1)
with open(SAMPLE_LIST_DIR + sample_list2_file,'rb') as fp2: sample_list2 = pickle.load(fp2)


# sample X species datafram with relative
# abundances data, trimmed at 0.0001
original_species_file = DATA_DIR + 'Mb_data_by_BD/LargeOrNewGenusSGBs_\
5MSubsampling_00001threshold_Noneratio_sSGB_byBD_PNP530Cardio126_swab.xlsx'
original_species_df_swab_name = '5M_species_PNP530Cardio126_swab'  # string

original_cluster_file = CLUSTER_ANALYSIS_DIR + \
                        'sampleByClusterDF_cohortfiltering005-09perc_dropped.dat'


##### global functions: ######
def gen_result_dir():
    result_dir_name = phen_name + '_' + \
                      'min_shared_species%s' % min_shared_species + '_' + \
                      'max_shared_species%s' % max_shared_species + '_' + \
                      'min_shared_cluster%s' % min_shared_cluster + '_' + \
                      'max_shared_cluster%s' % max_shared_cluster + '/'
    result_dir_name = result_dir_name.replace('.','')
    RESULT_DIR = TOP_RESULT_DIR + result_dir_name
    return RESULT_DIR

def write_vars_to_file(RESULT_DIR):

    print('saving variables to file')
    ctime = strftime("%Y%m%d_%H%M%S", gmtime())
    varFile = RESULT_DIR + 'var_file_%s.txt' % ctime
    new_file = open(varFile, mode="w")
    new_file.write('phen_name' + ": " + str(phen_name) + "\n")
    new_file.write('phen_file' + ": " + str(phen_file) + "\n")
    new_file.write('calc_species_phen' + ": " + str(calc_species_phen) + "\n")
    new_file.write('calc_cluster_phen' + ": " + str(calc_cluster_phen) + "\n")
    new_file.write('calc_cluster_species' + ": " + str(calc_cluster_species) + "\n")
    new_file.write('phen_file' + ": " + str(phen_file) + "\n")
    new_file.write('original_species_file' + ": " + str(original_species_file) + "\n")
    new_file.write('original_cluster_file' + ": " + str(original_cluster_file) + "\n")
    new_file.write('min_shared_species' + ": " + str(min_shared_species) + "\n")
    new_file.write('max_shared_species' + ": " + str(max_shared_species) + "\n")
    new_file.write('min_shared_cluster' + ": " + str(min_shared_cluster) + "\n")
    new_file.write('max_shared_cluster' + ": " + str(max_shared_cluster) + "\n")
    new_file.write('FDR_lim_species_phen' + ": " + str(FDR_lim_species_phen) + "\n")
    new_file.write('nTOP_species_phen' + ": " + str(nTOP_species_phen) + "\n")
    new_file.write('FDR_lim_cluster_phen' + ": " + str(FDR_lim_cluster_phen) + "\n")
    new_file.write('nTOP_cluster_phen' + ": " + str(nTOP_cluster_phen) + "\n")
    new_file.write('sample_list1_file' + ": " + str(sample_list1_file) + "\n")
    new_file.write('sample_list2_file' + ": " + str(sample_list2_file) + "\n")

    new_file.close()

    return


def get_shared_data_for_2_sample_lists(data_df, sample_list1, sample_list2,
                                       min_shared, max_shared, min_value=0.0001):
    '''
         this function takes a dataframe and two sample_lists and for each list calculates
         all features (species / clusters etc) that are present in at least min_shared and not
         more than max_shared samples of the sample list-->shared features
         then, it takes from the original data_df all features that are shared features in
         at least one of the sample_lists
         data_df: dataframe, sample X features
         sample_list1, sample_list2 - list of BD samples
         min_shared, max_shared - floats between 0 to 1. max_shared>min_shared
         '''
    data_shared_list = [];
    samples_shared_list = []
    for n, sample_list in enumerate([sample_list1, sample_list2]):
        #    assert sample_list in data_df.index.tolist(), 'not all samples in\
        # sample list %s are included in data_df.index' % n+1
        assert max_shared > min_shared, 'max_shared is not larger than min_shared'

        print('')
        print('sample list %s' % n)

        data = data_df.loc[sample_list].dropna(how='all')
        print('data.shape= ', data.shape)
        samples = data.index.tolist()

        # if species data, replace 0.0001 with np.nan:
        if data.min().min() == min_value:
            data = data.replace(min_value, np.nan)
        # get only data shared by min_shared-max_shared fraction of samples:
        cols = get__limited_sharing(data, min_shared, max_shared).columns.tolist()

        data_shared = data.loc[:, cols]
        print('data_shared.shape: ', data_shared.shape)

        data_shared_list.append(data_shared)
        samples_shared_list.append(samples)

    data_shared_1 = data_shared_list[0]
    data_shared_2 = data_shared_list[1]
    samples_shared_1 = samples_shared_list[0]
    samples_shared_2 = samples_shared_list[1]

    # get final data df  -PNP530 AND Cardio126 together, only balanced samples, only species/seqeuences that appear in at least 5 samples in
    # at least one cohort:
    final_data = data_df.loc[samples_shared_1 + samples_shared_2, list(
        set(data_shared_1.columns.tolist() + data_shared_2.columns.tolist()))]

    print('final_data.shape: ', final_data.shape)
    print('')

    return final_data


def plot_cluster_Mb_disease_relations2(clusterDF, mbDF, cluster_name, mb_name, ax=None, use_binary_cluster=True,
                                       cluster_as_cat=True):
    '''
    new function that uses seaborn's stripplot.
    clusterDF is a sample by cluster dataframe (not binary)
    mbDF is a sample by species dataframe (not binary) probably can use also genus here.
    cluster_name  = cluster head sequence (string)
    mb_name = formal species name (string)
    use_binary_cluster - whether to binarize the cluster information (to become False/True whether there is at least
    one sequence of the cluster in this sample). use false only if cluser_as_cat is True.
    cluster_as_cat - bool, whether to use the cluster info as the x-axis variable (False means to use the isCardio info instead)

    '''
    merged = pd.merge(pd.DataFrame(clusterDF[cluster_name]), pd.DataFrame(mbDF[mb_name]), how='inner', left_index=True,
                      right_index=True)
    y_values = merged[cluster_name].unique().tolist()
    merged['isCardio'] = np.where(merged.index.str.replace('BD', '').astype(int) < 950, 'Healthy', 'Patients')

    # define cluster info column to use: original or binary version
    if use_binary_cluster:
        merged['%s_binary' % cluster_name] = merged[cluster_name] > 0
        cluster_col = '%s_binary' % cluster_name
        xlabel_cluster = 'TCRseq cluster present'
    else:
        cluster_col = cluster_name
        xlabel_cluster = '# TCRseq in cluster'

    # define which column will be the x axis column:
    if cluster_as_cat:
        x_col = cluster_col
        split_col = 'isCardio'
        xlabel = xlabel_cluster
        legend_title = ''
        color_list = ['grey', 'darkred']
    else:
        x_col = 'isCardio'
        merged['%s_binary' % cluster_name] = merged['%s_binary' % cluster_name].astype('str')
        split_col = cluster_col
        xlabel = ''
        legend_title = 'TCRseq cluster present'
        color_list = ['forestgreen', 'red']

    if ax is None:
        fig, ax = plt.subplots()

    cmap = mpl.colors.ListedColormap(color_list, name='from_list', N=2)
    mpl.cm.register_cmap('from_list', cmap)
    cpal = sns.color_palette('from_list', n_colors=2)

    ax = sns.stripplot(x=x_col, y=mb_name, hue=split_col, data=merged, split=True, ax=ax, palette=cpal, size=10,
                       alpha=0.7)
    #     ax = sns.violinplot(x=x_col, y=mb_name, hue=split_col, data=merged,ax=ax,color='red')

    ax.set_ylabel('Species Relative Abundance')
    ax.set_xlabel(xlabel)

    marker_list = ['o', 'o']
    label_list = list(merged[split_col].unique())
    bbox_to_anchor = (1.01, 0.99)
    loc = 'upper left'
    make_custome_legend(color_list, marker_list, label_list, ax, title=legend_title, loc=loc,
                        bbox_to_anchor=bbox_to_anchor)

    mb_name_species = mb_name.split('|')[-4]
    mb_name_genus = mb_name.split('|')[-5]

    ax.set_title('%s\n%s\n%s' % (cluster_name, mb_name_species, mb_name_genus), y=1.03)
    ax.set_yscale('log')
    ax.set_ylim(0.00008, 0.1)
    return ax


class Calculator():
    def __init__(self):
        self._calc_species_phen = calc_species_phen
        self._calc_cluster_phen = calc_cluster_phen
        self._calc_cluster_species = calc_cluster_species
        self._merge_all_results = merge_all_results
        self.sample_list1 = sample_list1
        self.sample_list2 = sample_list2

    def _get_suffix(self):
        ## generate suffix:
        suffix = '_'
        for var in [FDR_lim_species_phen, nTOP_species_phen, FDR_lim_cluster_phen,
                    nTOP_cluster_phen, FDR_lim_cluster_species, nTOP_cluster_species]:
            if var is not None: suffix = ''.join([suffix, str(var).replace('.', '')])
        return suffix

    def get_samples_with_swabs(self, sample_list1, sample_list2, swab_sample_list):
        '''
        this functions edits the sample lists to include only
        samples that have swab species data
        return sample_list1_swab, sample_list2_swab
        '''
        sample_list1_swab = [sample for sample in sample_list1 if sample in swab_sample_list]
        sample_list2_swab = [sample for sample in sample_list2 if sample in swab_sample_list]

        return sample_list1_swab, sample_list2_swab

    def _species_phen_calculator(self, species_df, phen_df):
        '''
        this function calls the function analyze_associations_binary_nonbinary for calculating
        species_phen associations, print the header of the results and save the results (summary df
        and concise_summary df)

        '''
        print('')
        print('Starting species_phen_calculator...')

        summary_df_species_phen, concise_species_phen = analyze_associations_binary_nonbinary(
            nonbinaryDF=species_df,
            binaryDF=phen_df, nShared_pernonbinary=0, nShared_perbinary=0, nonbinary_zero=0.0001,
            nonbinary_name='species', binary_name=phen_name, make_concise=True, DirToSave=RESULT_DIR)

        print('species_phen_calculator result head:')
        print(summary_df_species_phen.head())
        print('saving results...')
        summary_df_species_phen.to_excel(RESULT_DIR + 'summary_df_species_phen.xlsx')
        concise_species_phen.to_excel(RESULT_DIR + 'concise_summary_df_species_phen.xlsx')
        print('done!')

        return summary_df_species_phen, concise_species_phen


    def _cluster_phen_calculator(self, cluster_df, phen_df):
        '''
        this function calls the function analyze_associations_binary_nonbinary for calculating
        cluster_phen associations, print the header of the results and save the results (summary df
        and concise_summary df)

        '''
        print('')
        print('Starting cluster_phen_calculator...')

        summary_df_cluster_phen, concise_cluster_phen = analyze_associations_binary_nonbinary(
            nonbinaryDF=cluster_df,
            binaryDF=phen_df, nShared_pernonbinary=0, nShared_perbinary=0, nonbinary_zero=0,
            nonbinary_name='cluster', binary_name=phen_name, make_concise=True, DirToSave=RESULT_DIR)

        print('cluster_phen_calculator result head:')
        print(summary_df_cluster_phen.head())
        print('saving results...')
        summary_df_cluster_phen.to_excel(RESULT_DIR + 'summary_df_cluster_phen.xlsx')
        concise_cluster_phen.to_excel(RESULT_DIR + 'concise_summary_df_cluster_phen.xlsx')
        print('done!')

        return summary_df_cluster_phen, concise_cluster_phen

    def _get_cluster_and_species_data_for_assoc_test(self,
                    concise_species_phen,concise_cluster_phen,
                    FDR_lim_species_phen=None, FDR_lim_cluster_phen=None,
                    nTOP_species_phen=None, nTOP_cluster_phen=None):
        '''
        return specific clusters and species (nTOP or all with corrected p-value<FDR_lim)
        '''

        # get species list:
        if FDR_lim_species_phen is not None:
            print ('generating a list of all species associated with the phenotype with\
        FDR<%s' %FDR_lim_species_phen)
            FDR_species_col = [col for col in concise_species_phen.columns if 'MW_p_corrPval_FDR' in col][0]
            top_species_list = (concise_species_phen[concise_species_phen[FDR_species_col] <= FDR_lim_species_phen])['species'].tolist()
        elif nTOP_species_phen is not None:
            print ('generating a list of top %s species associated with the phenotype' %nTOP_species_phen)
            try:
                top_species_list = (concise_species_phen[:nTOP_species_phen])['species'].tolist()
            except:
                top_species_list = concise_species_phen['species'].tolist()
        else:
            print ('no species threshold was defined, return all species')
            top_species_list = concise_species_phen['species'].tolist()

        #get cluster list:
        if FDR_lim_cluster_phen is not None:
            print ('generating a list of all cluster associated with the phenotype with\
        FDR<%s' %FDR_lim_cluster_phen)
            FDR_cluster_col = [col for col in concise_cluster_phen.columns if 'MW_p_corrPval_FDR' in col][0]
            top_cluster_list = (concise_cluster_phen[concise_cluster_phen[FDR_cluster_col] <= FDR_lim_cluster_phen])['cluster'].tolist()
        elif nTOP_cluster_phen is not None:
            print ('generating a list of top %s cluster associated with the phenotype' %nTOP_cluster_phen)
            try:
                top_cluster_list = (concise_cluster_phen[:nTOP_cluster_phen])['cluster'].tolist()
            except:
                top_cluster_list = concise_cluster_phen['cluster'].tolist()
        else:
            print ('no cluster threshold was defined, return all cluster')
            top_cluster_list = concise_cluster_phen['cluster'].tolist()

        print('top_cluster_list length=',len(top_cluster_list))
        print('top_species_list length=', len(top_species_list))
        clusters_species_test_num = len(top_cluster_list) * len(top_species_list)
        print ('expected number of tests: ', clusters_species_test_num )
        return top_cluster_list, top_species_list, clusters_species_test_num


    def _cluster_species_calculator(self, cluster_df, species_df, top_cluster_list, top_species_list):
        '''
        1. get cluster_df and species_df data only for top_clusters and top_species
        2. calculate MW and fisher test and generate result df and concise result df:
        call the function analyze_associations_binary_nonbinary for the calculations
        ######## need to copy it and its helper function edit_summaryDF to a .py file
        '''

        cluster_df = cluster_df.loc[:, top_cluster_list] #get data only for top clusters
        cluster_df = (cluster_df > 0).astype(int) #binarize cluster_df
        species_df = species_df.loc[:, top_species_list] # get data only for top species

        print('')
        print('Starting cluster_species_calculator...')

        summary_df_cluster_species, concise_cluster_species = analyze_associations_binary_nonbinary(
            nonbinaryDF=species_df,
            binaryDF=cluster_df, nShared_pernonbinary=0, nShared_perbinary=0, nonbinary_zero=0.0001,
            nonbinary_name='species', binary_name='cluster', make_concise=True, DirToSave=RESULT_DIR)

        print('cluster_species_calculator result head:')
        print(summary_df_cluster_species.head())
        print('saving results...')
        summary_df_cluster_species.to_excel(RESULT_DIR + 'summary_df_cluster_species.xlsx')
        concise_cluster_species.to_excel(RESULT_DIR + 'concise_summary_df_cluster_species.xlsx')
        print('done!')

        return summary_df_cluster_species, concise_cluster_species

    def _merge_cluster_species_disease_data(self, concise_cluster_phen,
                                            concise_species_phen, concise_cluster_species,
                                            suffix):
        '''
        1. get only significant data from each concise_summary df (species_phen, cluster_phen, cluster_species)
        according to defined thresholds.
        2. merge, calculate directionOK (that the directions of relationships don't contradict each other)
        3. add annotations
        '''

                ### get only top cluster_species pairs:
        if FDR_lim_cluster_species is not None:
            print ('generating a list of all cluster associated with the species with\
        FDR<%s' %FDR_lim_cluster_species)
            FDR_cluster_col = [col for col in concise_cluster_species.columns if 'MW_p_corrPval_FDR' in col][0]
            concise_cluster_species = concise_cluster_species[concise_cluster_species[FDR_cluster_col] <= FDR_lim_cluster_species]
        elif nTOP_cluster_species is not None:
            print ('generating a list of top %s cluster associated with species' %nTOP_cluster_species)
            try:
                concise_cluster_species = concise_cluster_phen[:nTOP_cluster_phen]
            except:
                pass
        else:
            print ('no cluster-species threshold was defined, return all cluster-species pairs')

        ### rename columns to enable clean merge:
        concise_cluster_phen = concise_cluster_phen.rename\
            (columns = {'MW_p_corrPval_FDR0.1':'MW_p_corrPval_FDR0.1_cluster_phen',
                        '-log10p_MW':'-log10p_MW_cluster_phen',
                        'Associated_with_%s_positive' %phen_name: 'cluster_pos_associated_with_phen'})
        concise_species_phen = concise_species_phen.rename \
            (columns={'MW_p_corrPval_FDR0.1': 'MW_p_corrPval_FDR0.1_species_phen',
                      '-log10p_MW':'-log10p_MW_species_phen',
                                   'Associated_with_%s_positive' % phen_name: 'species_pos_associated_with_phen'})
        concise_cluster_species = concise_cluster_species.rename \
            (columns={'MW_p_corrPval_FDR0.1': 'MW_p_corrPval_FDR0.1_cluster_species',
                      '-log10p_MW':'-log10p_MW_cluster_species',
                                   'Associated_with_cluster_positive' : 'species_pos_associated_with_cluster'})

        cluster_species_phen_summary_df = pd.merge(concise_cluster_species,concise_cluster_phen,how='inner',
                                                   left_on='cluster',right_on='cluster')
        cluster_species_phen_summary_df = pd.merge(cluster_species_phen_summary_df, concise_species_phen, how='inner',
                                                   left_on='species', right_on='species')
        cluster_species_phen_summary_df = cluster_species_phen_summary_df.sort_values\
            (by='MW_p_corrPval_FDR0.1_cluster_species')

        print('cluster_species_phen_summary_df.shape: ',cluster_species_phen_summary_df.shape)

        # calcualte that all relation direction are ok:
        cluster_species_phen_summary_df['directionOK'] = np.where(
            cluster_species_phen_summary_df['species_pos_associated_with_cluster'] == cluster_species_phen_summary_df[
                'species_pos_associated_with_phen'] + cluster_species_phen_summary_df[
                'cluster_pos_associated_with_phen'], 0, 1)

        print(cluster_species_phen_summary_df.iloc[:4, :4])
        cluster_species_phen_summary_df.to_excel(
            RESULT_DIR + 'cluster_species_phen_summary_df_%s.xlsx' %suffix)

        # add cluster sequences and annotations:
        print ('getting cluster sequences and annotations:')
        df = cluster_species_phen_summary_df
        output_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/sharingAnalysis/seqClusters_allProd_maxdist1/'
        cluster_head_col = 'cluster'
        cluster_species_phen_summary_df_withAnot, phen_related_identities = get_annot_for_seqs_in_cluster(df, cluster_head_col,
                                            output_dir, getIdentities=True, cluster_head_list=None, identity_file=None,
                                            sumIdentities=True)
        cluster_species_phen_summary_df_withAnot.to_excel(
            RESULT_DIR + 'cluster_species_phen_summary_df_%s_withAnot.xlsx' %suffix)
        return cluster_species_phen_summary_df_withAnot

    def plot_p_value_distributions(self,FDR=None):
        '''
        this function takes the cluster_species_phen_summary_df and plot -log10 p-values for each comparison seperately
        (cluster_phen, species_phen, cluster_species) and enables comparison of results
        '''
        result_list = [
            ('TCR seq cluster-%s' %phen_name, 'concise_summary_df_cluster_phen'),
            ('Species-%s' %phen_name, 'concise_summary_df_species_phen'),
            ('Species-TCR seq cluster', 'concise_summary_df_cluster_species')
        ]
        ### generate settings:
        fig, ax = plt.subplots(figsize=(10, 6))
        color_list = ["red", "orange", "black"]
        cmap = mpl.colors.ListedColormap(color_list, name='from_list', N=3) #generate cmap from the color list
        mpl.cm.register_cmap('from_list', cmap)
        cpal = sns.color_palette('from_list', n_colors=3) #generates palette from the cmap

        total_df = pd.DataFrame()
        for n, result in enumerate(result_list):
            name = result[0]
            f_name = result[1]
            if FDR is None: FDR = 0.1
            print(n, name, f_name)
            try:
                df = pd.read_excel(RESULT_DIR + f_name + '.xlsx')
            except:
                print ('the file %s doesnt exist, cant plot p-value distributions' %f_name)
                return

            for col in df.columns:
                if 'Associated_with_' in col: col_to_calc = col  # get the name of the column that holds info 
            df['direction'] = np.where(df[col_to_calc] == 0, -1, 1)
            df['log10p_MW_direction'] = df['direction'] * df['-log10p_MW']
            df['color'] = df['MW_p_corrPval_FDR0.1'].apply(
                lambda x: 'red' if x < FDR else 'orange' if x < FDR+0.05 else 'black')
            df = df[['log10p_MW_direction', 'color']]
            df['goal'] = name
            total_df = pd.concat([total_df, df])

        # print('total_df.shape', total_df.shape)
        # print('')

        ax = sns.stripplot(x='goal', y='log10p_MW_direction',hue='color', data=total_df,
                           jitter=0.02, ax=ax, alpha=0.2, palette=cpal)

        ax.set_xlabel('')
        ax.set_ylabel('Log 10 p-value (MW test)\n(*association direction)')
        ax.set_ylim(-10, 10)

        # make legend:
        marker_list = ["o", "o", "o"]
        label_list = ['<'+str(FDR), str(FDR)+'-'+str(FDR+0.05), '>'+str(FDR+0.05)]  # need to make it used by the plotting itself
        loc = None
        bbox_to_anchor = None
        title = 'Corrected p-value'
        make_custome_legend(color_list, marker_list, label_list, ax, loc=loc,
                            bbox_to_anchor=bbox_to_anchor, title=title)

        fig.savefig(RESULT_DIR +
        'TCRseqCluster_MBspecies_disease_association_pvaluesplot_FDR%s.png' %FDR)

        return fig, ax


    def plot_5_top_triads(self,cluster_df, species_df, suffix):
        '''
        this function automatically plot the top 5 triads using the
        plot_cluster_Mb_disease_relations2 function and save the plots'
        '''

        try:
            df = pd.read_excel(RESULT_DIR + 'cluster_species_phen_summary_df_%s_withAnot.xlsx' %suffix)
        except:
            print ('merged file doesnt exist')
            return

        ### arrange merged result by descending order of sum of logs

        df['log_sum'] = df.loc[:,[col for col in df.columns if '-log' in col]].sum(axis=1)
        df = df.sort_values(by='log_sum', ascending = False)
        top5_pair_list = zip(df['cluster'],df['species'])

        ###choose top 5 pairs and run the plot function for each, including saving figs:
        for pair in top5_pair_list:
            cluster = pair[0]; species = pair[1]

            fig, ax = plt.subplots(figsize=(6, 6))
            cluster_name = cluster; mb_name = species
            mb_name_species = mb_name.split('|')[-4].replace('_', '')
            mb_name_genus = mb_name.split('|')[-5].replace('_', '')
            use_binary_cluster = True
            cluster_as_cat = True

            ax = plot_cluster_Mb_disease_relations2(cluster_df, species_df, cluster_name, mb_name, ax,
                                                    use_binary_cluster=use_binary_cluster,
                                                    cluster_as_cat=cluster_as_cat)
            fig.savefig(RESULT_DIR + 'cluster_species_pair_plots/' + 'cluster_species_disease_relations_%s_%s%s_dividedByDisease%s_TCRbinary%s.png' % (
            cluster_name, mb_name_genus, mb_name_species, str(cluster_as_cat), str(use_binary_cluster)),
            bbox_inches='tight')

        return fig, ax


##### run script: #####

if __name__ == '__main__':

    ###define calculator and its attributes:

    Calculator = Calculator()
    Calculator._suffix = Calculator._get_suffix()
    Calculator.sample_list1 = sample_list1 #this should be replaced with samples with swabs only if species data is included
    Calculator.sample_list2 = sample_list2 #this should be replaced with samples with swabs only if species data is included

    ### define result folders and generate them if necessary:
    RESULT_DIR = gen_result_dir()
    RESULT_FIG_DIR = RESULT_DIR+'figures/'

    if not isdir(RESULT_DIR):
        makedirs(RESULT_DIR)
    if not isdir(RESULT_FIG_DIR):
        makedirs(RESULT_FIG_DIR)

    ### define which calculations are necessary, according to the existance of files and global parameters:
    if 'summary_df_species_phen.xlsx' in listdir(RESULT_DIR) and not rewrite:
        Calculator._calc_species_phen = False
        print('summary_df_species_phen.xlsx already exists in RESULT_DIR and rewrite=False,\
    thus skipping species_phen_association calculator!')
    if 'summary_df_cluster_phen.xlsx' in listdir(RESULT_DIR) and not rewrite:
        print('summary_df_cluster_phen.xlsx already exists in RESULT_DIR and rewrite=False,\
    thus skipping cluster_phen_association calculator!')
        Calculator._calc_cluster_phen = False
    if 'summary_df_cluster_species.xlsx' in listdir(RESULT_DIR) and not rewrite:
        print('summary_df_cluster-species.xlsx already exists in RESULT_DIR and rewrite=False,\
    thus skipping cluster_species_association calculator!')
        Calculator._calc_cluster_species = False
    if 'cluster_species_phen_summary_df_%s_withAnot.xlsx' % Calculator._suffix \
            in listdir(RESULT_DIR) and not rewrite:
        print('cluster_species_phen_summary_df_%s_withAnot.xlsx already exists in RESULT_DIR and rewrite=False,\
    thus skipping results merging!' %Calculator._suffix)
        Calculator._merge_all_results = False

    write_vars_to_file(RESULT_DIR)


    ### get species data and list of samples with swabs if necessary:
    if Calculator._calc_species_phen or Calculator._calc_cluster_species or plot_top_pairs:
        print('species data is included in analysis,get swab samples data:')

        # load data:
        Calculator.original_species_df_swab = pd.read_excel(original_species_file).set_index('BD')

        # get lists of samples with swab data:
        Calculator.sample_list1, Calculator.sample_list2 = \
            Calculator.get_samples_with_swabs(sample_list1, sample_list2,
                                              Calculator.original_species_df_swab.index.tolist())
        with open (RESULT_DIR + 'sample_list1_swab.pkl','wb') as fp1:
            pickle.dump(Calculator.sample_list1,fp1)
        with open (RESULT_DIR + 'sample_list2_swab.pkl','wb') as fp2:
            pickle.dump(Calculator.sample_list2,fp2)

        # get shared species data (note that the df is re-filled with 0.0001)
        Calculator._species_df = get_shared_data_for_2_sample_lists(Calculator.original_species_df_swab,
                                                                    Calculator.sample_list1, Calculator.sample_list2,
                                                                    min_shared_species, max_shared_species).fillna(0.0001)
        if Calculator._calc_species_phen:
            # calculate species_phen associations:
            Calculator.summary_df_species_phen, Calculator.concise_species_phen = \
                    Calculator._species_phen_calculator(Calculator._species_df, phen_df)
        else:
            print('skipping species_phen assoication calculation')

    ###get cluster data:
    if Calculator._calc_cluster_phen or Calculator._calc_cluster_species or plot_top_pairs:
        print('cluster data is included in analysis')

        # load data:
        Calculator.original_cluster_df = pd.read_pickle(original_cluster_file)

        # get shared cluster data (note that the df is filled with 0's)
        Calculator._cluster_df = get_shared_data_for_2_sample_lists(Calculator.original_cluster_df,
                                                                    Calculator.sample_list1, Calculator.sample_list2,
                                                                    min_shared_cluster,
                                                                max_shared_cluster).fillna(0)
        if Calculator._calc_cluster_phen:
            # calculate cluster_phen associations:
            Calculator.summary_df_cluster_phen, Calculator.concise_cluster_phen = \
                Calculator._cluster_phen_calculator(Calculator._cluster_df, phen_df)
        else:
            print('skipping cluster_phen association calculation')

        if Calculator._calc_cluster_species:
            ### get data and calculate cluster-species associations only if calc_cluster_species is True
            ### and the resulting file doesn't exist or rewrite is True:

            print ('getting cluster_phen and species_phen association data:')

            #load concise summary files:
            concise_species_phen = pd.read_excel(RESULT_DIR + 'concise_summary_df_species_phen.xlsx')
            concise_cluster_phen = pd.read_excel(RESULT_DIR + 'concise_summary_df_cluster_phen.xlsx')

            Calculator.top_cluster_list, Calculator.top_species_list, clusters_species_test_num = \
                Calculator._get_cluster_and_species_data_for_assoc_test(concise_species_phen,
                        concise_cluster_phen, FDR_lim_species_phen, FDR_lim_cluster_phen,
                        nTOP_species_phen, nTOP_cluster_phen)

            #calculate cluster_species associations if numer of potential associations to test is
            # more than 0 and less than 400000:
            if (clusters_species_test_num > 0) and (clusters_species_test_num < 400000):
                Calculator.summary_df_cluster_species, Calculator.concise_cluster_species = \
                    Calculator._cluster_species_calculator(Calculator._cluster_df, Calculator._species_df,
                                                           Calculator.top_cluster_list, Calculator.top_species_list)
            else:
                print ('clusters_species_test_num = %s and therefor cluster_species associations\
    were not calculated' %clusters_species_test_num)

        else:
            print ('skipping cluster_species association calculation')

    ### merge results from all association tests:
    if Calculator._merge_all_results:
        try:
            print('merge cluster,species and phen data')
            concise_species_phen = pd.read_excel(RESULT_DIR + 'concise_summary_df_species_phen.xlsx')
            concise_cluster_phen = pd.read_excel(RESULT_DIR + 'concise_summary_df_cluster_phen.xlsx')
            concise_cluster_species = pd.read_excel(RESULT_DIR + 'concise_summary_df_cluster_species.xlsx')

            cluster_species_phen_summary_df_withAnot = Calculator._merge_cluster_species_disease_data(
                concise_cluster_phen, concise_species_phen, concise_cluster_species, Calculator._suffix)
        except:
            print ('couldnt merge all results')
    else:
        print ('skipping merging all results')


    ###plots:
    if plot_pvalues:
        print('plotting p-values')
        fig,ax = Calculator.plot_p_value_distributions(FDR_for_plotting)
        plt.show()

    if plot_top_pairs:
        print ('plotting and saving top 5 pairs')
        fig,ax = Calculator.plot_5_top_triads(Calculator._cluster_df, Calculator._species_df,
                                              Calculator._suffix)



    print('Done!')
    end_time=datetime.datetime.now()
    print ('started on %s' %start_time)
    print ('ended on %s' %end_time)
    print ('total time= %s' %(end_time - start_time))
