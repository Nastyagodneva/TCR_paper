import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from os import makedirs
from os.path import isdir
from ShaniBA.CardioProject.paper_summer2019.fig5.ver2.calc_roc_auc_new_outcome import AucCalculator

if __name__ == "__main__":
    # directories:
    GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
    DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
    SAMPLE_LIST_DIR = GENIE_BASE_DIR + 'Sample files/BD lists/'
    CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'
    FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5_updated_phens/'
    PRED_RESULTS_DIR='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/'
    OR_DIR = FIGURE_DIR + 'or_calcs/'

    if not isdir(OR_DIR): makedirs(OR_DIR)

    #files and dataframe:
    outcome_file = CARDIO_PHEN_DIR + 'new_outcome_df.xlsx'
    outcome_df = pd.read_excel(outcome_file).set_index('BD')

    pred_proba_file =PRED_RESULTS_DIR+'isCardio/XGB_randomSearch_25_byRepFeatPCA10RelsVDJnocorr0999AgeGender/\
predictions_df.pkl'

    cluster_data_file = DATA_DIR + 'TCR_seqs_clusters/newX_onlySeqs_025_085_noNans.dat'
    cluster_data = pd.read_pickle(cluster_data_file)
    # cluster_data_5 = cluster_data[['CASSPTGSETQYF','CASSLGWGGEQYF','CASSLQQGNTEAFF',
    #                  'CASSLETGVYEQYF','CASSPTQDYGYTF']]
    # cluster_data_5.to_excel(FIGURE_DIR + 'short_cluster_table.xlsx')
    # cluster_data = pd.read_excel(FIGURE_DIR + 'short_cluster_table.xlsx')

    mb_data_file = DATA_DIR + 'Mb_data_by_BD/LargeOrNewGenusSGBs_\
5MSubsampling_00001threshold_Noneratio_sSGB_byBD_PNP530Cardio126_swab.xlsx'
    mb_df = pd.read_excel(mb_data_file).set_index('BD')



class Fig5SubplotAbs(object):

    def __init(self,ax,cluster,outcome,cluster_data,outcome_df):
        self.ax = ax
        self.cluster_df = cluster_data
        self.outcome_df = outcome_df
        self.cluster = cluster
        self.outcome = outcome


    def gen_plot(self):
        pass


class Fig5SubplotBoxRocplot(Fig5SubplotAbs):

    def __init__(self, ax,cluster,outcome,cluster_data,outcome_df):
        super(self.__class__, self).__init__()
        self.ax = ax
        self.cluster = cluster
        self.outcome = outcome
        self.cluster_df = cluster_data
        self.outcome_df = outcome_df


    def _gen_data(self):
        df = pd.merge(pd.DataFrame(self.outcome_df[self.outcome]),
                      pd.DataFrame(self.cluster_df[self.cluster]),
                      how='inner', left_index=True, right_index=True)
        df = df.dropna(how='any')
        return df


    def gen_boxplot(self):
        self.data = self._gen_data()
        # self.ax = sns.boxplot(x=self.outcome, y=self.cluster, data=self.data)
        self.ax = sns.swarmplot(y=self.outcome, x=self.cluster, data=self.data, color=".25")
        self.ax.set_xlabel(self.outcome)
        self.ax.set_ylabel(self.cluster)
        return self.ax


    def gen_roc_plot(self):
        self.data = self._gen_data()
        feature_value_list = self.data[self.cluster].tolist()
        outcome_value_list = self.data[self.outcome].tolist()
        return AucCalculator(feature_value_list=feature_value_list,
                             outcome_value_list=outcome_value_list,ax=self.ax)._plot_roc_auc()

    def calc_or(self,cluster_num_thresh=0):
        self.data = self._gen_data()
        self.max_seq_in_cluster = self.data[self.cluster].max().max()
        self.data[self.cluster] = (self.data[self.cluster] > cluster_num_thresh).astype(int)
        tab = pd.crosstab(index=self.data[self.cluster],columns=self.data[self.outcome], margins=False)

        exposed_cases = tab.iloc[1, 1]
        exposed_non_cases = tab.iloc[1, 0]
        non_exposed_cases = tab.iloc[0,1]
        non_exposed_non_cases = tab.iloc[0, 0]

        self.odd_ratio = float(exposed_cases * non_exposed_non_cases) / \
                         (exposed_non_cases * non_exposed_cases)

        self.upper_95_ci = np.e ** (np.log(self.odd_ratio) +
                               1.96*(math.sqrt(1./exposed_cases + 1./exposed_non_cases +
                                               1./non_exposed_cases + 1./non_exposed_non_cases)))
        self.lower_95_ci = np.e ** (np.log(self.odd_ratio) -
                               1.96 * (math.sqrt(1. / exposed_cases + 1. / exposed_non_cases +
                                                 1. / non_exposed_cases + 1. / non_exposed_non_cases)))
        self.is_or_sig = not ((1 > self.lower_95_ci) and (1 < self.upper_95_ci))

        return self.odd_ratio, self.upper_95_ci, self.lower_95_ci, self.is_or_sig, self.max_seq_in_cluster



def summarize_or(cluster_outcome_pair_list, table_file_name='six_clusters_result_df'):
    """
    This method summarizes  odd ratios and their 95% confidence interval for having the outcome, given the presence
    level of the TCR cluster
    :param cluster_outcome_pair_list: list of tuples. each tuple contains the outcome and the cluster (str, str)
    :return: df with or results for all cluter_outcome pairs
    """

    df = pd.DataFrame()
    t_list = [0, 1, 2]
    count = 0
    for pair in cluster_outcome_pair_list:
        outcome = pair[0]
        cluster = pair[1]

        for t in t_list:

            try:
                odd_ratio, upper_95_ci, lower_95_ci, is_or_sig, max_seq_in_cluster = \
                    Fig5SubplotBoxRocplot(ax=None, cluster=cluster,
                                          outcome=outcome, cluster_data=cluster_data, outcome_df=outcome_df).calc_or(
                        cluster_num_thresh=t)

                df.loc[count, 'cluster'] = cluster
                df.loc[count, 'outcome'] = outcome
                df.loc[count, 'seq num min.'] = t + 1
                df.loc[count, 'odd_ratio'] = odd_ratio
                df.loc[count, 'upper_95_ci'] = upper_95_ci
                df.loc[count, 'lower_95_ci'] = lower_95_ci
                df.loc[count, 'is_or_sig'] = is_or_sig
                df.loc[count, 'max_seq_in_cluster'] = max_seq_in_cluster

            except:
                print ('couldnt calculate for %s %s %s' % (cluster, outcome, t))

            count += 1

    df.to_excel(OR_DIR + table_file_name + '.xlsx')
    return df


def plot_or(or_result_table,cluster_list=None, cluster_num_thresh=0,ax=None,
            markersize=None, capsize=2,markerfacecolor='orange',
            markeredgecolor='orange',ecolor='black' ):

    if cluster_list is None:
        cluster_list = or_result_table['cluster'].unique().tolist()
    if ax is None:
        fig,ax = plt.subplots()

    #prepare the table:
    or_result_table = or_result_table[(or_result_table['cluster'].isin(cluster_list)) &
                                      (or_result_table['seq num min.']==cluster_num_thresh+1)]

    #get the data:
    x = or_result_table['odd_ratio']
    y = range(len(cluster_list))
    low_e_err  = or_result_table['lower_95_ci']
    high_e_err = or_result_table['upper_95_ci'] - x

    #plot the data
    ax.errorbar(x, y, xerr=(low_e_err,high_e_err), ls='',ms=markersize, capsize=capsize, fmt = 'o',
                markerfacecolor=markerfacecolor, markeredgecolor=markeredgecolor,
                ecolor=ecolor)
    print ("min_OR: %s" %low_e_err.min())
    print ("max_OR: %s" % high_e_err.max())
    ax.set_xlabel('OR')
    ax.set_ylabel('Cluster')
    ax.set_yticks(y)
    ax.set_yticklabels(cluster_list)
    ax.set_xlim (0, round(np.max(or_result_table['upper_95_ci']))+3)
    ax.axvline(1, ls='--', c='grey')

    return ax



class Fig5SubplotScatterplot(Fig5SubplotAbs):

    def __init__(self, ax, cluster, outcome, cluster_data, outcome_df,mb_df,species):
        super(Fig5SubplotScatterplot, self).__init__(ax, cluster, outcome, cluster_data, outcome_df)
        self.mb_df = mb_df
        self.species = species


    def _gen_data(self):
        df = pd.merge(pd.DataFrame(self.outcome_df[self.outcome]),
                      pd.DataFrame(self.cluster_df[self.cluster]),
                      how='inner', left_index=True, right_index=True)
        df = pd.merge(df,
                      pd.DataFrame(self.mb_df[self.species]),
                      how='inner', left_index=True, right_index=True)
        df = df.dropna(how='any')
        return df



    def gen_plot(self):
        self.data = self._gen_data()
        self.ax = plt.scatter(x=self.data[self.species],
                              y=self.data[self.cluster],
                              c=self.data[self.outcome])
        self.ax.set_xlabel(self.species)
        self.ax.set_ylabel(self.outcome)
        return self.ax

if  __name__ == '__main__':


    updated_cluster_list = pd.read_excel(OR_DIR+'top_roc_auc_rehosp05.xlsx').iloc[:,0].tolist()
    updated_cluster_list_with_phen = zip(len(updated_cluster_list)*['CV hospitalization including chest pain'],
                                         updated_cluster_list)

    or_result_table = summarize_or(updated_cluster_list_with_phen, 'OR_top_roc_auc_rehosp05')




    cluster_outcome_pair_list = [
                ('CV hospitalization including chest pain','CASSPTGSETQYF'),
                ('CV hospitalization including chest pain','CASSLGWGGEQYF'),
                 ('CV hospitalization including chest pain','CASSLQQGNTEAFF'),
                 ('CV hospitalization including chest pain','CASSLETGVYEQYF'),
                 ('CV hospitalization including chest pain','CASSPTQDYGYTF'),
                 ]

    fig, ax = plt.subplots()
    or_result_table = summarize_or(cluster_outcome_pair_list)
    ax = plot_or(or_result_table,cluster_list=None, cluster_num_thresh=0,ax=ax)
    plt.subplots_adjust(left=0.25)
    fig.savefig(FIGURE_DIR + 'or_plot_preliminary.png', dpi=150)

    fig2,axes = plt.subplots(nrows=5,ncols=1, sharex=True, figsize = (4,20))

    for n, pair in enumerate(cluster_outcome_pair_list):
        outcome = pair[0]
        cluster = pair[1]
        ax = axes.flatten()[n]

        Fig5SubplotBoxRocplot(ax=ax,cluster=cluster,outcome=outcome,cluster_data=cluster_data,
                              outcome_df=outcome_df).gen_roc_plot()

        if n!=2:
            ax.set_ylabel('')
        if n!=4:
            ax.set_xlabel('')

    plt.subplots_adjust(left=0.3, hspace = 0.3)
    fig2.savefig(FIGURE_DIR + 'roc_curves.png', dpi=150)





