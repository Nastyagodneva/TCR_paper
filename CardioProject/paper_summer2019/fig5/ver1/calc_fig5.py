# imports and definitions
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
PRED_RESULTS_DIR='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5/'

def new_outcomes():
    outcome_df = pd.read_excel(CARDIO_PHEN_DIR+'outcomeDF.xlsx').set_index('BD')

    outcome_df['mi_stroke_PCI'] = np.where(outcome_df['Acute MI'].isnull(),np.nan,
                                  np.where(outcome_df.iloc[:, :3].sum(axis=1) > 0, 1, 0))
    outcome_df['any_outcome'] = np.where(outcome_df['Acute MI'].isnull(),np.nan,
                                np.where(outcome_df.iloc[:,:4].sum(axis=1)>0,1,0))

    fu = pd.read_excel(CARDIO_PHEN_DIR + 'follow_up_period_length_2019-05-13.xlsx').set_index('BD')
    merged = pd.merge(outcome_df,fu,how='left',left_index=True,right_index=True)
    for fu_time in [547,638,730]:
        cond1 = merged['any_outcome']==1; cond2 = merged['follow_up_period_length']<fu_time
        merged['early developers'] = np.where(merged['Acute MI'].isnull(),np.nan,
                            (np.where(cond1 & cond2,1,0)))
        outcome_df['early developers_%s' %fu_time] = merged['early developers']

    outcome_df.to_excel(CARDIO_PHEN_DIR + 'outcomeDF.xlsx')
    for col in outcome_df.columns:
        outcome_df[col].to_excel(CARDIO_PHEN_DIR + col + '.xlsx')

    return outcome_df

def get_followup_period(follow_up_end="2019-05-13"):
    def conv_datetime(s):
        return datetime.datetime.strptime(s, '%Y-%m-%d')

    follow_up_end_str = follow_up_end.replace('_','')
    disc = pd.read_excel(CARDIO_PHEN_DIR + 'DischargeDate.xlsx')
    disc['Discharge Date'] = disc['Discharge Date'].apply(lambda x: pd.Timestamp(x))
    disc['follow_up_period_length'] = conv_datetime(follow_up_end) - disc['Discharge Date']
    disc2 = disc[['BD', 'follow_up_period_length']].set_index('BD')
    disc2.to_excel(CARDIO_PHEN_DIR + 'follow_up_period_length_%s.xlsx'%follow_up_end_str)
    return

def remove_outliers_from_col(df,col,std):
    col_mean=df[col].mean()
    col_std = df[col].std()
    lower_bound=col_mean - std*col_std
    higher_bound = col_mean + std * col_std
    df[col] = np.where((df[col]>lower_bound) &\
                       (df[col]<higher_bound),df[col],np.nan)
    return df

def gen_data():
    def calc_correct(row):
        if (row['pred_delta'] > -0.05) and (row['pred_delta'] < 0.05):
            return 'no change'
        elif (row['pred_delta'] < -0.05) and (row['isCardio'] == 0):
            return 'correct'
        elif (row['pred_delta'] > 0.05) and (row['isCardio'] == 1):
            return 'correct'
        else:
            return 'in correct'

    # get pred proba data for predicted age+gender
    Y_pred_age_gender=pd.read_pickle(PRED_RESULTS_DIR+'isCardio/XGB_randomSearch_25_byPredictedAgeGender/\
predictions_df.pkl')

    # get pred_proba for all
    Y_pred_all=pd.read_pickle(PRED_RESULTS_DIR+'isCardio/XGB_randomSearch_25_byRepFeatPCA10RelsVDJnocorr0999AgeGender/\
predictions_df.pkl')


    Y_pred_real = pd.merge(Y_pred_age_gender.rename(columns={'isCardio':'pred_proba_age_gender'}),
                           Y_pred_all.rename(columns={'isCardio':'pred_proba_all'}),
                           how='outer',left_index=True, right_index=True).dropna(how='all')
    Y_pred_real['pred_delta'] = Y_pred_real['pred_proba_all'] - Y_pred_real['pred_proba_age_gender']
    # get isCardio
    Y_pred_real['isCardio'] = np.where((Y_pred_real.index.str.replace('BD','').astype(int)) > 949,1,0)

    # calculate delta correctness:
    Y_pred_real['delta correctness'] = Y_pred_real.apply(calc_correct,axis=1)

    # get main phens:
    phen_df = pd.read_excel(DATA_DIR + 'phenotypes_byBD/PNP530Cardio126_Age_Gender\
_Male_HbA1C_SmokingStatus_Yes_BMI_HDL_Smoking_ever_nonHDL_Cholesterol.xlsx').set_index('BD')
    Y_pred_real = pd.merge(Y_pred_real,phen_df,how='left', left_index=True,right_index=True)

    #get 4 outcomes
    outcomeDF = pd.read_excel(CARDIO_PHEN_DIR+'outcomeDF.xlsx').set_index('BD')
    Y_pred_real = pd.merge(Y_pred_real,outcomeDF,
                           how='left', left_index=True,right_index=True)

    #get STEMI
    diagnosis = pd.read_excel(CARDIO_PHEN_DIR+'AdmissionDiagnosis.xlsx').set_index('BD')
    Y_pred_real = pd.merge(Y_pred_real,pd.DataFrame(diagnosis),
                           how='left', left_index=True,right_index=True)

    #get predicted_age:
    pred_age = pd.read_pickle(PRED_RESULTS_DIR + 'Cardio126_phens_basedOnHealthy/pred_and_real_Age.pkl')[
        ['average pred', 'Age']]
    Y_pred_real = pd.merge(Y_pred_real,
                           pd.DataFrame(pred_age['average pred']).rename(
                               columns={'average pred': 'pred_age'}),
                           how='left', left_index=True, right_index=True)
    #get isCardio pred proba based on real age and gender:
    pred_proba_real_age_gender = pd.read_pickle(
        PRED_RESULTS_DIR+'isCardio/XGB_randomSearch_25_byAgeGender/predictions_df.pkl')
    Y_pred_real = pd.merge(Y_pred_real,
                       pred_proba_real_age_gender.rename(
                           columns={'isCardio': 'pred_proba_real_age_gender'}),
                       how='left', left_index=True, right_index=True)
    # get isOstial phen:
    is_ostial = pd.read_excel(CARDIO_PHEN_DIR+'OstialDiseaseCorrect_binary.xlsx')\
        .set_index('BD')
    Y_pred_real = pd.merge(Y_pred_real, pd.DataFrame(is_ostial),
                           how='left', left_index=True, right_index=True)

    # get pred_proba_all_new_prediction:
    pred_proba_all_new_prediction = pd.read_pickle(PRED_RESULTS_DIR + \
                'isCardio_new_predictions/XGB_byTCRfeaturesONLY_optByAUC/predictions_df.pkl')
    Y_pred_real = pd.merge(Y_pred_real,
                   pred_proba_all_new_prediction.rename(
                                    columns={'isCardio': 'pred_proba_all_new_prediction'}),
                                   how='left', left_index=True, right_index=True)

    # get pred_proba_all_with_real_age_gender_new_prediction:
    pred_proba_all_with_real_age_gender_new_prediction = pd.read_pickle(PRED_RESULTS_DIR + \
                  'isCardio_new_predictions/XGB_byTCRfeaturesANDrealAgeGender_optByAUC/predictions_df.pkl')
    Y_pred_real = pd.merge(Y_pred_real,
                           pred_proba_all_with_real_age_gender_new_prediction.rename(
                               columns={'isCardio': 'pred_proba_all_with_real_age_gender_new_prediction'}),
                           how='left', left_index=True, right_index=True)

    Y_pred_real.to_excel(FIGURE_DIR + 'Y_pred_real_phens.xlsx')
    return Y_pred_real

def compare_descriptive_stats(phen,Y_pred_real):
    # compare delta and delta correctness for outcomes:
    for col in Y_pred_real.columns:
        try:
            Y_pred_real[col] = Y_pred_real[col].astype(float)
        except:
            print ('couldnt convert col %s to float' %col)

    outcome_list = ['Unplanned PCI','CV hospitalization including chest pain']
    cols = ['pred_proba_age_gender','pred_proba_all','pred_delta', 'delta correctness']


    for name,group in Y_pred_real.groupby(phen):
        print name
        print group[['pred_proba_age_gender','pred_proba_all','pred_delta']].describe()
    print ('')

    return


def compare_delta_cols(phen,Y_pred_real,col1='pred_proba_age_gender',
                       col2='pred_proba_all',color_by_age=False,
                       remove_outliers=None):
    '''
    this function calculate the delta between two columns, and compare it
    among two phen groups
    :param phen: string. name of binary phen.
    'Unplanned PCI'/'CV hospitalization including chest pain' etc.
    :param Y_pred_real:
    :param col1: string, name of col in Y_pred_real
    :param col2: string, name of col in Y_pred_real
    :param color_by_age: bool. if False, will color by phen (neg/pos)
    :return:
    '''
    for col in [col1, col2]:
        try:
            Y_pred_real[col] = pd.to_numeric(Y_pred_real[col])
        except:
            print ('column %s is not numeric, can perform comparison')
            return

    # compare delta_pred for Unplanned PCI negative and positive:
    df = Y_pred_real[Y_pred_real['isCardio'] == 1]

    df['delta_col'] = df[col2]-df[col1]
    if remove_outliers is not None:
            df = remove_outliers_from_col(df, 'delta_col', remove_outliers)

    #(1) calculate statistical difference between groups regarding delta_col values

    data={}
    for name,group in df.groupby(phen):
        data[name] = group['delta_col'].tolist()
    data_0=data[0]; data_1=data[1]
    t,p = ttest_ind(data_0,data_1,nan_policy='omit')
    s_mw,p_mw = mannwhitneyu(data_0,data_1)

    #(2) plot scatter plot of the two cols:

    fig,ax = plt.subplots(figsize=(10,10))
    cmap = mpl.colors.ListedColormap(['green', 'red'], name='from_list', N=2)

    if color_by_age:
        scat = ax.scatter(df[col1], df[col2],
                          c=df['Age'], s=df[phen] * 100 + 20)
    else:
        scat = ax.scatter(df[col1],df[col2],
               c=df[phen],cmap=cmap,s=df['Age'])

    plt.colorbar(scat)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.plot([0,1],[0,1],c='black',lw=3)
    ax.set_title(phen+'\n(colored by age %s' %color_by_age)
    ax.set_xlabel('pred_proba_age_gender')
    ax.set_ylabel('pred_proba_all')
    # ax.legend()
    ax.text(0.01,0.99,'p_mw=%s' %p_mw, transform=ax.transAxes, ha='left',va='top')

    fig.savefig(FIGURE_DIR+'pred_proba_comp_%s_%s_%s.png' %(phen,color_by_age,remove_outlier))

    # plot box plot:
    fig2, ax2 = plt.subplots(figsize=(10, 10))
    sns.boxplot(x=phen, y='delta_col', data=df,ax=ax2)
    ax2.set_ylabel('delta: %s-\n%s' %(col2,col1))
    ax2.text(0.01, 0.99, 'p_mw=%s' % p_mw, transform=ax2.transAxes, ha='left', va='top')
    fig2.savefig(FIGURE_DIR + '%s_%s_delta_bb_%s_%s.png' % (col1,col2,phen,remove_outlier))
    print ('p_mw=%s' % p_mw)
    print ('p_t=%s' % p)

    return fig, fig2, p_mw


def compare_one_col(phen,Y_pred_real,col,remove_outliers=None):
    '''

    :param phen:
    :param Y_pred_real:
    :param col:
    :param remove_outliers: default None. if not None, int/float that define outlier
    removal limit (how many stds from mean is the limit)
    :return:
    '''
    try:
        Y_pred_real[col] = pd.to_numeric(Y_pred_real[col])
    except:
        print ('column %s is not numeric, can perform comparison')
        return

        #(1) calculate statistical difference between groups regarding col values

    data={}
    df = Y_pred_real[Y_pred_real['isCardio'] == 1]
    if remove_outliers is not None:
        df = remove_outliers_from_col(df,col,remove_outliers)
    for name,group in df.groupby(phen):
        data[name] = group[col].tolist()
    data_0=data[0]; data_1=data[1]

    t,p = ttest_ind(data_0,data_1,nan_policy='omit')
    s_mw,p_mw = mannwhitneyu(data_0,data_1)

    #(2) plot box plot for col


    fig2, ax2 = plt.subplots(figsize=(10, 10))
    sns.boxplot(x=phen, y=col, data=df)

    ax2.text(0.01, 0.99, 'p_mw=%s' % p_mw, transform=ax2.transAxes, ha='left', va='top')
    fig2.savefig(FIGURE_DIR + '%s_bb_%s_%s.png' % (col, phen,remove_outlier))
    print ('p_mw=%s' %p_mw)
    print ('p_t=%s' %p)
    return  fig2, p_mw

def specific_calc_for_STEMI(Y_pred_real):

    Y_pred_real['AdmissionDiagnosis'].value_counts()
    data={}
    for name, group in Y_pred_real.groupby('AdmissionDiagnosis'):
        data[name] = group['pred_delta'].dropna().tolist()
    data_1 = data[1]
    data_2 = data[2]
    data_3 = data[3]
    t, p = ttest_ind(data_2, data_3,nan_policy='omit')
    s_mw, p_mw = mannwhitneyu(data_2, data_3)

    df = Y_pred_real[Y_pred_real['isCardio'] == 1]
    fig,ax = plt.subplots(figsize=(10,10))
    cmap = mpl.colors.ListedColormap(['grey','green', 'red'], name='from_list', N=3)
    ax.scatter(df['pred_proba_age_gender'],df['pred_proba_all'],
               c=df['AdmissionDiagnosis'],cmap=cmap,s=100)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.plot([0,1],[0,1],c='black',lw=3)
    ax.set_title('AdmissionDiagnosis')
    ax.set_xlabel('pred_proba_age_gender')
    ax.set_ylabel('pred_proba_all')
    ax.legend()
    ax.text(0.01,0.99,'p_mw=%s' %p_mw, transform=ax.transAxes, ha='left',va='top')
    return




###################33

if __name__ == '__main__':

    outcome_df = new_outcomes()
    Y_pred_real = gen_data()
    #
    phen_list = ['early developers_730']
    # # phen_list=['mi_stroke_PCI','any_outcome','CV hospitalization including chest pain',
    # #            'Unplanned PCI']
    # # remove_outliers_list = [None,1,3]
    remove_outliers_list = [None,1]
    #
    for phen in phen_list:
        for remove_outlier in remove_outliers_list:

            print phen,remove_outlier

            #compare pred_proba
            print ('compare pred_proba abs')
            col='pred_proba_all'
            compare_one_col(phen,Y_pred_real,col,remove_outliers=remove_outlier)

            # compare pred_proba
            print ('compare pred_proba_real_age_gender abs')
            col = 'pred_proba_real_age_gender'
            compare_one_col(phen, Y_pred_real, col, remove_outliers=remove_outlier)

            # compare pred_proba_new
            print ('compare pred_proba abs new')
            col = 'pred_proba_all_new_prediction'
            compare_one_col(phen, Y_pred_real, col, remove_outliers=remove_outlier)

            # compare delta between real age and predicted age:
            print('compare delta between predicted and real age')
            col1='Age'; col2='pred_age'
            compare_delta_cols(phen,Y_pred_real,col1=col1,col2=col2,
                               color_by_age=False,remove_outliers=remove_outlier)


            # compare delta between new prediction with features and not predicted age gender, to
            # prediction by age and gender:
            print('compare delta between new prediction and pred by age gender')
            col1 = 'pred_proba_real_age_gender';
            col2 = 'pred_proba_all_new_prediction'
            compare_delta_cols(phen, Y_pred_real, col1=col1, col2=col2,
                               color_by_age=False, remove_outliers=remove_outlier)

             #compare delta between full prediction and prediction by age and gender:
            print('compare delta between real prediction and pred by age gender')
            col1='pred_proba_real_age_gender'; col2='pred_proba_all'
            compare_delta_cols(phen,Y_pred_real,col1=col1,col2=col2,
                               color_by_age=False,remove_outliers=remove_outlier)

            print('compare delta pred_proba_all and pred_proba_age_gender')
            col1 = 'pred_proba_age_gender';
            col2 = 'pred_proba_all'
            compare_delta_cols(phen, Y_pred_real, col1=col1, col2=col2,
                               color_by_age=False, remove_outliers=remove_outlier)

