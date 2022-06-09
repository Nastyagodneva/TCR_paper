from scipy.stats import ttest_ind, mannwhitneyu
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ShaniBA.CardioProject.Figures.GeneralFigureFunctions import *


def compare_numeric_data_between_2_groups(df,phen_col,data_col):
    data = {}
    for name, group in df.groupby(phen_col):
        data[name] = group[data_col].tolist()
    data_0 = data[0];
    data_1 = data[1]
    s_t, p_t = ttest_ind(data_0, data_1, nan_policy='omit')
    s_mw, p_mw = mannwhitneyu(data_0, data_1)
    return p_t, p_mw


def gen_colormap_from_color_list(color_list):
    cmap = mpl.colors.ListedColormap(color_list, name='from_list', N=len(color_list))
    return cmap


def gen_palette_from_color_list(color_list, n_colors=2):
    cmap=gen_colormap_from_color_list(color_list)
    mpl.cm.register_cmap('from_list', cmap)
    cpal = sns.color_palette('from_list', n_colors=n_colors)
    return cpal


def plot_boplot_with_datapoints(df,x_column, y_column,color_list, ax=None, point_color='black'):
    if ax is None:
        fig,ax = plt.subplots()
    cpal = gen_palette_from_color_list(color_list)
    ax = sns.boxplot(x=x_column, y=y_column, data=df,ax=ax,palette=cpal)
    ax = sns.swarmplot(x=x_column, y=y_column, data=df, color=point_color)
    return ax


def plot_deciles(ax, y_test, y_pred, target_name='ACS', q=10, show_fold=False,
                 color='black', **lineproprs):
    """

    :param y_test: one column dataframe with binary classification (target)
    :param y_pred: dataframe with pred_proba columns
    :param q: number of quantiles (for deciles use q=10...)
    :param show_fold: bool, whether or not to show fold change in plot
    :param target_name: string, target name to appear on y label (for example, 'ACS'
    :param lineproprs: line style arguments
    :return: ax
    """
    y_test.columns=['test']
    y_pred.columns = ['pred']

    test_deciles = pd.merge(y_test, y_pred, how='inner', right_index=True, left_index=True)

    test_deciles['dec'] = pd.qcut(test_deciles.pred, q, labels=range(1, q + 1))

    test_deciles = test_deciles.groupby('dec').test.mean()

    test_deciles = test_deciles * 100

    ax.plot(test_deciles.index.astype(int), test_deciles.values, marker='o', color=color, **lineproprs)

    if show_fold & (test_deciles[0] != 0):
        fold = (test_deciles[-1] / test_deciles[0])
        ax.annotate('Fold increase: %d' % fold, xy=[1, test_deciles[0]], xytext=[0.1, 0.7], textcoords='axes fraction',
                     arrowprops={'arrowstyle': '-'}, color='black')
        ax.annotate('', xy=[10, test_deciles[-1]], xytext=[0.4, 0.76], textcoords='axes fraction',
                     arrowprops={'arrowstyle': '-'}, color='black')

    ax.set_xlim([0.5, 10.5])
    ax.set_xlabel('Prediction decile')
    ax.set_xticks(range(1, q + 1))
    ax.set_ylabel('%s (perc)' %target_name)

    return ax, test_deciles


def gen_shap_summary_to_axes(ax, fig, target_name, nTopFeatures, shap_df,
                             features_df, jitter=None, sample_num=None,
                             feature_name_list=None, scalingMethod='minmax',
                             addColorBar=True, alpha=1, **stripplot_kwargs):
    # ##optional values for scalingMethod: 'minmax','scale','perc',None

    # calcualte feature order by sum of abs shap values and take top n values:
    if nTopFeatures is not None:
        topN_features_df = shap_df.abs().sum().sort_values(ascending=False)[:nTopFeatures]
    else:
        topN_features_df = shap_df.abs().sum().sort_values(ascending=False)
    topN_features = topN_features_df.index.tolist()
    topN_features_abs_shap = topN_features_df.values.tolist()
    print topN_features
    topN_features_new = edit_feature_names(topN_features)
    print topN_features_new

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
            features_df_topN_scaled = pd.DataFrame(index=features_df_topN.index, columns=features_df_topN.columns,
                                                   data=scaler.fit_transform(features_df_topN))
        if scalingMethod == 'scale':
            # center to the mean and to get unit variance
            scaler = preprocessing.scale(copy=True)
            features_df_topN_scaled = pd.DataFrame(index=features_df_topN.index, columns=features_df_topN.columns,
                                                   data=scaler.fit_transform(features_df_topN))
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
    merged = pd.merge(features_df_topN_scaled_long, shap_df_topN_long, how='inner', left_on=['index', 'variable'],
                      right_on=['index', 'variable']).rename(columns={'index': 'Sample', 'variable': 'feature'})

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
    stam1 = plt.scatter(merged['shap_value'], merged['feature_value'], c=merged['feature_value'], cmap='RdBu_r')
    plt.clf()

    # generate specific ax for the colorbar and plot the colorbar:
    if addColorBar:
        divider = make_axes_locatable(ax)
        ax_cb = divider.new_horizontal(size="3%", pad=0.1)
        fig.add_axes(ax_cb)
        cb1 = fig.colorbar(stam1, cax=ax_cb, orientation='vertical')
        ax_cb.yaxis.set_ticks([0, 0.5, 1])
        ax_cb.yaxis.set_ticklabels([0, 0.5, 1])
        ax_cb.set_ylabel('Norm Feature Value', labelpad=0.8)

    # plot the scatter plot itself:
    plot = sns.stripplot(x='feature', y='shap_value', hue='feature_value', data=merged,
                         palette='RdBu_r',
                         alpha=alpha, ax=ax, s=3, jitter=jitter, **stripplot_kwargs)

    plot.get_legend().set_visible(False)
    ax.set_xlabel("")
    ax.set_ylabel('Shap Value (%s)' % target_name)
    ax.set_xticklabels(topN_features_new)  # change ticklables to readable feature names
    plt.setp(ax.get_xticklabels(), rotation=90, ha="left", va='center',
             rotation_mode="anchor", fontsize=mpl.rcParams['font.size']-0.5)

    # color ticklables according to their content:
    for n, lab in enumerate(ax.get_xticklabels()):
        str_lab = str(lab)
        #         print n,str(lab)
        if ('insertion' in str_lab) or ('deletion' in str_lab) or ('length' in str_lab):
            plt.setp(ax.get_xticklabels()[n], color='green')
        elif ('Norm' in str_lab) or ('Top' in str_lab) or ('count' in str_lab) \
                or ('Shannon' in str_lab) or ('Berger' in str_lab) or ('Simpson' in str_lab):
            plt.setp(ax.get_xticklabels()[n], color='red')
        elif ('seq #' in str_lab) or ('seq freq.' in str_lab) or ('annotate' in str_lab):
            plt.setp(ax.get_xticklabels()[n], color='orange')
        elif ('non-prod' in str_lab):
            plt.setp(ax.get_xticklabels()[n], fontstyle='italic')

    # ax.set_xticks([x + 0.35 for x in ax.get_xticks()])
    ax.tick_params(axis='x', direction="out", pad=43)

    return ax, fig, topN_features_abs_shap

