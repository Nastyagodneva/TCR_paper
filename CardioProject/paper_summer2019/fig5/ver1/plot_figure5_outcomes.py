import pandas as pd
import matplotlib.pyplot as plt
import cPickle as pickle

from paper_summer2019.fig5.ver1.calc_fig5_clustering import get_outcome
from paper_summer2019.fig5.ver1.calc_fig5 import remove_outliers_from_col
from ShaniBA.CardioProject.paper_summer2019.figure_utils import compare_numeric_data_between_2_groups, plot_boplot_with_datapoints

GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
PRED_RESULTS_DIR='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/'

# TODO: add significance mark to boxplots

with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/Cardio126', 'rb') as fp:
    Cardio126 = pickle.load(fp)
##################################

# functions:


def get_data_fig5a(outcome, outlier_t=1):
    '''
    this function plots arrange the data for part A of figure 5: comparing the delta of pred_proba values (full prediciton based on TCR info,
    and prediciton based on age and gender) between negative and positive outcome patients
    :param outcome: string, outcome name
    :param outlier_t: int/float (how many stdev below and above the mean should serve as threshold to remove outliers
    :return: Y_pred, data df
    '''

    # data required: outcome_df, pred_proba data
    any_outcome_df = get_outcome(outcome)

    pred_proba_all = Y_pred_all = pd.read_pickle(PRED_RESULTS_DIR + 'isCardio/XGB_randomSearch_25_byRepFeatPCA10RelsVDJnocorr0999AgeGender/\
predictions_df.pkl')

    pred_proba_real_age_gender = pd.read_pickle(
        PRED_RESULTS_DIR + 'isCardio/XGB_randomSearch_25_byAgeGender/predictions_df.pkl')

    # arrange data:
    Y_pred = pd.merge(pred_proba_all, pred_proba_real_age_gender, how='inner',
                      left_index=True, right_index=True, suffixes=['_by_all', '_by_real_age_gender'])
    Y_pred = pd.merge(Y_pred, any_outcome_df, how='inner',
                      left_index=True, right_index=True)

    Y_pred['delta_pred'] = Y_pred['isCardio_by_all'] - Y_pred['isCardio_by_real_age_gender']
    Y_pred = remove_outliers_from_col(Y_pred, 'delta_pred', outlier_t)
    for col in Y_pred.columns:
        Y_pred[col] = pd.to_numeric(Y_pred[col])

    return Y_pred

def plot_fig_5a(outcome,outlier_t):
    ##figure 5a:
    Y_pred = get_data_fig5a(outcome, outlier_t)

    # calculate difference:
    p_t, p_mw = compare_numeric_data_between_2_groups(Y_pred, outcome, 'delta_pred')

    # plot differences:
    fig, ax = plt.subplots(figsize=(8, 10))
    color_list = ['red', 'darkred']
    ax = plot_boplot_with_datapoints(df=Y_pred.dropna(how='any'), x_column=outcome, y_column='delta_pred', ax=ax,
                                     color_list=color_list)
    ax.set_title('prediction probability increase\nby TCR feature info',fontsize='x-large')
    ax.set_xlabel('Any outcome', fontsize='x-large')
    ax.set_ylabel('Change in prediction probability\nof being ACS patient', fontsize='x-large')
    ax.set_xticklabels(['Neg', 'Pos'])

    plt.subplots_adjust(left=0.2)
    fig.savefig(FIGURE_DIR + '5a.png')

    return ax


def plot_fig_5b(outcome):
    # get data:
    any_outcome_df = get_outcome(outcome)
    feature_df = pd.read_excel(DATA_DIR +
                               'TCRfeatures/TCRfeaturesNoT_PCA10percShared_withRels_VDJ\
_noCorr0-999_nanFilled_noConsts_Cardio126.xlsx').set_index('BD')

    CMV_patients = pd.merge(pd.DataFrame(feature_df.loc[Cardio126, 'CMV_rel_cum_freq(perc)']),
                            any_outcome_df, how='inner', left_index=True, right_index=True)

    # compare distributions:
    p_t, p_mw = compare_numeric_data_between_2_groups(CMV_patients, outcome, 'CMV_rel_cum_freq(perc)')

    # plot:
    fig, ax = plt.subplots(figsize=(8, 10))
    color_list = ['red', 'darkred']
    ax = plot_boplot_with_datapoints(df=CMV_patients.dropna(how='any'), x_column=outcome,
                                     y_column='CMV_rel_cum_freq(perc)', ax=ax,
                                     color_list=color_list)
    ax.set_title('CMV-related sequences frequencies',fontsize='x-large')
    ax.set_xlabel('Any outcome in\nfirst 21m', fontsize='x-large')
    ax.set_ylabel('cumulative relative frequency', fontsize='x-large')
    ax.set_xticklabels(['Neg','Pos'])

    plt.subplots_adjust(left=0.2)
    fig.savefig(FIGURE_DIR + '5b.png')

    return ax


##################################
# variables:
outcome_5a = 'any_outcome'
outcome_5b = 'early developers_638'
outlier_t = 1

if __name__ == '__main__':

    ## fig 5a: pred proba delta comparison:
    ax = plot_fig_5a(outcome=outcome_5a, outlier_t = outlier_t)

    ## figure 5b: CMV and HSV comparison
    ax = plot_fig_5b(outcome = outcome_5b)


##C: node and edges graph for species-cluster pairs

#data required: cluster_species_phen_summary file for any_outcome