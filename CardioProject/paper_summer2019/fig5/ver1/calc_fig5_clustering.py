# imports and definitions
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, mannwhitneyu
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn import preprocessing
import os
import cPickle as pickle
from ShaniBA.GeneralFeaturePhenotypeInteractions.Feature_phenotype_functions import add_corrected_pValues,\
    compare_TCR_between_groups, calc_PCA


GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
PRED_RESULTS_DIR='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5/'

with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/Cardio126', 'rb') as fp:
    Cardio126 = pickle.load(fp)

def get_suffix(only_productive,only_annotation,only_calc_features,get_suffix):
    d = {'only_productive': only_productive,
         'only_annotation': only_annotation,
         'only_calc_features': only_calc_features,
         'filter_negatives_only': filter_negatives_only}
    suffix = ''
    for k,v in d.items():
        if v:  suffix = suffix + '_' + k
    return suffix

def get_productive_only(df):
    col_list = [col for col in df.columns if col.endswith('_1')]
    return df[col_list]

def get_annotation_only(df):
    col_list = [col for col in df.columns if 'cum_freq' in col]
    return df[col_list]

def get_calc_features_only(df):
    col_list=df.columns.tolist()
    for prefix in ['V','D','J']:
        col_list = [col for col in col_list if not col.startswith(prefix)]
        col_list = [col for col in col_list if 'rel' not in col]
    return df[col_list]

def filter_by_followup_length(df,outcome_df,filter_negatives_only=True,length_limit=547,follow_up_end="2019-05-13"):
    #length_limit = 1.5 years
    outcome = outcome_df.columns.tolist()[0]
    follow_up_end_str = follow_up_end.replace('_', '')
    fu = pd.read_excel(CARDIO_PHEN_DIR + 'follow_up_period_length_%s.xlsx' % follow_up_end_str).set_index('BD')
    df = pd.merge(df,fu,how='left',left_index=True,right_index=True)
    df = pd.merge(df,outcome_df,how='left',left_index=True,right_index=True)
    if filter_negatives_only:
        df_negatives = df[(df['follow_up_period_length'] > length_limit) & \
                          (df[outcome] ==0)]
        df_positives = df[df[outcome] == 1]
        df = pd.concat([df_negatives,df_positives])
    else:
        df = df[df['follow_up_period_length']> length_limit]
    return df


def get_patients_feature_df(patient_list,only_productive=False,only_annotation=False,only_calc_features=False):
    feature_df=pd.read_excel(DATA_DIR +
                             'TCRfeatures/TCRfeaturesNoT_PCA10percShared_withRels_VDJ\
_noCorr0-999_nanFilled_noConsts_Cardio126.xlsx').set_index('BD')

    feature_df = feature_df.loc[patient_list,:]

    if only_productive: feature_df = get_productive_only(feature_df)
    elif only_annotation: feature_df = get_annotation_only(feature_df)
    elif only_calc_features: feature_df = get_calc_features_only(feature_df)

    return feature_df

def compare_TCR_features(feature_df,outcome_df,filter_negatives_only,filter_by_fu_length=None):
    if filter_by_fu_length is not None:
        feature_df = filter_by_followup_length(feature_df, outcome_df, filter_negatives_only=filter_negatives_only,
                    length_limit=filter_by_fu_length)
    outcome = outcome_df.columns.tolist()[0]
    groupName1 = outcome + '_0'; groupName2 = outcome + '_1'
    TCRdf_binary=pd.DataFrame()

    negatives = outcome_df[outcome_df[outcome]==0].index.tolist()
    positives = outcome_df[outcome_df[outcome] == 1].index.tolist()
    sample_list1 = [BD for BD in feature_df.index if BD in negatives]
    sample_list2 = [BD for BD in feature_df.index if BD in positives]

    print ('# patients in %s is %s ' %(groupName1,len(sample_list1)))
    print ('# patients in %s is %s ' % (groupName2, len(sample_list2)))

    feature_df = feature_df.drop(outcome, axis=1)
    feature_comparison_df,TCR_comparison_df = compare_TCR_between_groups(sample_list1,groupName1,
                                sample_list2,groupName2,feature_df,TCRdf_binary,
                                compare_features=True,compare_seqs=False)
    return feature_comparison_df

def normalize_df(df):
    return pd.DataFrame(index=df.index, columns=df.columns, data=preprocessing.scale(df,copy=False))

def get_outcome(outcome):
    outcome_df=pd.read_excel(CARDIO_PHEN_DIR +'outcomeDF.xlsx').set_index('BD')
    return pd.DataFrame(outcome_df[outcome])

def get_patients_clusters(patient_list, outcome_df,do_binary = True,filter_by_fu_length=None,filter_negatives_only=False):
    print ('loading TCR cluster file...')
    cluster_df = pd.read_pickle(DATA_DIR + 'TCR_seqs_clusters/\
newX_onlySeqs_040_085_noNans.dat')

    cluster_df = cluster_df.loc[patient_list, :]

    if filter_by_fu_length is not None:
        cluster_df = filter_by_followup_length(cluster_df, outcome_df, filter_negatives_only,length_limit=filter_by_fu_length)

    if do_binary:
        print ('binarizing file...')
        cluster_df = (cluster_df>0).astype(int)
    print ('Done!')
    return cluster_df

def gen_clustermap(df,only_partial=False,filter_by_fu_length=None,
                   do_cluster=True, filter_negatives_only=False,
                   figsize_frac=None, outcome_to_color=None):

    row_cluster = True; col_cluster=True
    if not do_cluster:
        row_cluster = False
        col_cluster = True
    outcome = outcome_to_color.columns.tolist()[0]
    df = pd.merge(df, outcome_to_color, how='inner', left_index=True, right_index=True)
    color_map = {0: 'b', 1: 'r', np.nan: 'g'}
    row_colors = df.iloc[:, -1].map(color_map)

    if only_partial:
        df = df.iloc[:10,:10]
    print ('df shape before filtering by fu period:' ,df.shape)
    if filter_by_fu_length is not None:
        df = filter_by_followup_length(df, outcome_to_color, filter_negatives_only,length_limit=filter_by_fu_length)
    print ('df shape after filtering by fu period:', df.shape)
    if figsize_frac is not None:
        height = np.ceil(df.shape[0] * figsize_frac)
        width = np.ceil(df.shape[1] * figsize_frac)
    else:
        height = 10; width = 10
    # if outcome_to_color is not None:
    #     df = pd.merge(df,outcome_to_color,how='inner',left_index=True,right_index=True)
    #     color_map = {0:'b',1:'r',np.nan:'g'}
    #     row_colors = df.iloc[:,-1].map(color_map)
    #     outcome=outcome_to_color.columns.tolist()[0]
    print ('df shape after filtering by fu period:', df.shape)
    df = df.sort_values(by=outcome, ascending=False)
    if filter_by_fu_length is not None: df = df.drop('follow_up_period_length',axis=1)
    g = sns.clustermap(df.drop(outcome,axis=1),xticklabels=True,yticklabels=True,
                       figsize=(width,height),row_colors=row_colors,
                       row_cluster=row_cluster, col_cluster=col_cluster)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize='xx-small')
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), fontsize='xx-small')
    plt.subplots_adjust(bottom=0.3)
    # g.ax_heatmap.set_title('Subplot Title', y=1.25)
    return g

def gen_heatmap(df,only_partial=False,figsize_frac=None,outcome_to_color=None):
    outcome = outcome_to_color.columns.tolist()[0]
    df = pd.merge(df, outcome_to_color, how='inner', left_index=True, right_index=True)
    color_map = {0: 'b', 1: 'r', np.nan: 'g'}
    row_colors = df.iloc[:, -1].map(color_map)
    df = df.sort_values(by=outcome,ascending=False)
    if only_partial:
        df = df.iloc[:10,:10]
    if figsize_frac is not None:
        height = np.ceil(df.shape[0] * figsize_frac)
        width = np.ceil(df.shape[1] * figsize_frac)
    else:
        height = 10; width = 10


    g = sns.clustermap(df.drop(outcome,axis=1),xticklabels=True,yticklabels=True,
                       figsize=(width,height),row_colors=row_colors,
                       row_cluster=False, col_cluster=False)
    plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), fontsize='xx-small')
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), fontsize='xx-small')
    plt.subplots_adjust(bottom=0.3)
    return g

def plot_PCA(feature_df,outcome_df,suffix,to_scale=True):
    outcome = outcome_df.columns[0]
    sample_list1 = (outcome_df[outcome_df[outcome]==0]).index.tolist()
    sample_list2 = (outcome_df[outcome_df[outcome] == 1]).index.tolist()
    print ('outcome: ', outcome)
    print ('sample list 1:', sample_list1)
    print ('sample list 2:', sample_list2)

    sample_list_list=[(outcome+'_0', sample_list1),(outcome+'_1', sample_list2)]
    color_list=['blue','red']

    PCAdf,fig,ax = calc_PCA(feature_df,isSparse=False,n_comp=10,sample_list_list=sample_list_list,
              color_list=color_list,pca_n1_toplot=0,pca_n2_toplot=1,toScale=to_scale,
             toPlot=True, toAnnotate=True,fig=None,ax=None,calculateSeperation=True)
    ax.set_title('PCA_%s%s' %(outcome,suffix))
    plt.subplots_adjust(right=0.8)

    return PCAdf

if __name__ == "__main__":


    outcome = 'CV hospitalization including chest pain'
    filter_by_fu_length= None
    filter_negatives_only=False
    only_productive=False
    only_annotation=False
    only_calc_features=True
    feature_for_pca='clusters' #'clusters' / 'features'

    suffix = get_suffix(only_productive,only_annotation,only_calc_features,filter_negatives_only)


    print ('loading feature df...')
    feature_df = get_patients_feature_df(patient_list=Cardio126,only_productive=only_productive,only_annotation=only_annotation,
                                only_calc_features=only_calc_features)
    print ('loading outcome df...')
    outcome_df = get_outcome(outcome)

    print ('loading cluster df...')
    cluster_df = get_patients_clusters(patient_list=Cardio126, outcome_df=outcome_df, do_binary=True,
                                       filter_by_fu_length=filter_by_fu_length,
                                       filter_negatives_only=filter_negatives_only)


    # print ('comparing TCR features...')
    # feature_comparison_df = compare_TCR_features(feature_df,outcome_df,filter_negatives_only,filter_by_fu_length)
    # feature_comparison_df.to_excel(os.path.join(FIGURE_DIR+'feature_comparison/',
    #                     'feature_comparison_df%s_%s_%s.xlsx' %(suffix,outcome,str(filter_by_fu_length).replace('.',''))))
    #
    # print feature_comparison_df.head()
    # print ('normalizing df...')
    # feature_df_scaled = normalize_df(feature_df)

    print ('calculating PCA...')
    if feature_for_pca == 'clusters':
        suffix = suffix + '_clusters'
        PCAdf = plot_PCA(cluster_df,outcome_df,suffix,to_scale=False)
    else:
        PCAdf = plot_PCA(feature_df, outcome_df, suffix)



    # print('generating heatmapmap...')
    # do_cluster=False
    # g = gen_clustermap(feature_df_scaled,only_partial=False,filter_by_fu_length=filter_by_fu_length,
    #                    do_cluster=do_cluster, filter_negatives_only=filter_negatives_only,
    #                    figsize_frac=0.25,outcome_to_color=outcome_df)
    #
    #
    # print('generating clustermap...')
    # do_cluster=True
    # g = gen_clustermap(feature_df_scaled,only_partial=False,filter_by_fu_length=filter_by_fu_length,
    #                    do_cluster=do_cluster,
    #                    figsize_frac=0.25,outcome_to_color=outcome_df)

    # df = pd.DataFrame(np.random.randint(0,99,(4,10)))
    # g2 = sns.clustermap(df)
    # plt.setp(g2.ax_heatmap.yaxis.get_majorticklabels(), fontsize='xx-large',rotation=0)
    # g2.ax_heatmap.set_title('Subplot Title',y=1.25)


    # iris = sns.load_dataset("iris")
    # species = iris.pop("species")
    # lut = dict(zip(species.unique(), "rbg"))
    # row_colors = species.map(lut)
    # >>> g = sns.clustermap(iris, row_colors=row_colors)
    # g = sns.clustermap(iris)
