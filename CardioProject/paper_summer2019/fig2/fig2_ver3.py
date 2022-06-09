
#imports:
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from collections import OrderedDict

from matplotlib import gridspec
from ShaniBA.PredictionPipeline.PredictionFunctions import plot_roc_curve, plot_PR_curve
from ShaniBA.CardioProject.paper_summer2019.figure_utils import plot_deciles, gen_shap_summary_to_axes


#directories:
GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig2_isCardio_classification/'
FIGURE_FILES_DIR = FIGURE_DIR + 'files_for_figure/'


class Figure2:

    def __init__(self,datasets_to_use=[],
                 size_factor=1, cmap='Reds'):
        self._dataset_dict = OrderedDict({
        'shuffled': FIGURE_FILES_DIR + 'XGB_randomSearch_25_byOldXShuffleWithPredictedAgeGender/',
        'pred. age gender': FIGURE_FILES_DIR + 'XGB_randomSearch_25_byPredictedAgeGender/',
        'all features': FIGURE_FILES_DIR + 'XGB_randomSearch_25_byRepFeatPCA10percVDJ0999PredictedAgeGender/',
        # 'TCR features by mb' : FIGURE_FILES_DIR +'TCRfeatures_XGB_byMbSega10M005095binary_SWAB_optByExpVar/'
        })
        self._datasets_to_use = ['pred. age gender', 'all features']
        self._size_factor = size_factor
        self._cmap = cmap
        self._plot_color_list = self._get_color_list


    @property
    def _get_color_list(self):
        cmap = plt.get_cmap(self._cmap, len(self._datasets_to_use)+1)  # get seperate colors for graphs
        color_list = [cmap(x + 1) for x in range(len(self._datasets_to_use))]
        return color_list

    def _get_data(self):
        """
        this function extracts all dataframe required for the plot and declare them as global variables
        all the files for these dataframes were copied to a designated directory (FIGURE_FILES_DIR) so they
        need to be manually updated if updated in their original location
        """
        data_dict = {}

        for dataset_name, data_dir in self._dataset_dict.items():
            print ('now extracting data for dataset %s' %dataset_name)
            data_dict[dataset_name] = {}
            data_dict[dataset_name]['result_df'] = pd.read_pickle(data_dir + 'results_df.pkl')
            data_dict[dataset_name]['prediction_df'] = pd.read_pickle(data_dir + 'predictions_df.pkl')
            data_dict[dataset_name]['shap_df'] = pd.read_pickle(data_dir + 'shap_values.pkl')
        data_dict['isCardio'] = pd.read_pickle(FIGURE_FILES_DIR + 'isCardio.dat')
        data_dict['feature_df'] = pd.read_excel(DATA_DIR +
                    'TCRfeatures/TCRfeaturesNoT_PCA10percShared_withRels_VDJ_noCorr0-999_nanFilled_noConsts.xlsx').set_index('BD')

        self._data_dict = data_dict

    def _set_fig2_definitions(self):
        """
        this function sets matplotlib parameters to be constant across the whole figure.
        :param size_factor: a factor to multiply all parameters in order to efficiently increase/decrese the figures
        while keeping the correct size relations.
        :return:
        """
        self._figsize=(5*self._size_factor, 4.75*self._size_factor)

        #TODO: arrange this funciton

        changing_size_params = {
        'axes.labelsize': 5,
        'font.size': 5,
        'legend.fontsize': 4,
        'axes.titlesize':7,
        'xtick.labelsize': 5,
        'ytick.labelsize': 5,
        'xtick.major.size': 1.75,
        'ytick.major.size': 1.75,
        'axes.linewidth': 0.1,
        'lines.linewidth': 0.5,
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
        'figure.dpi': 200,
        # 'axes.labelpad': 0.8,
        'legend.labelspacing': 0.4,
        'legend.framealpha': 0.2,
        'legend.frameon': True,
        'legend.fancybox': True,
            'axes.xmargin': 0.1,
            'axes.ymargin': 0.1,
        }

        if self._size_factor != 1:
            for k,v in changing_size_params.items():
                changing_size_params[k] = v*self._size_factor


        mpl.rcParams.update(changing_size_params)
        mpl.rcParams.update(other_params)

        return

    def _plot_roc_curves(self, ax):

        for n,dataset in enumerate(self._datasets_to_use):

            y_pred = self._data_dict[dataset]['prediction_df']
            result_df = self._data_dict[dataset]['result_df']
            y_true = self._data_dict['isCardio']
            plot_random = True if n ==0 else False

            ax, roc_auc=plot_roc_curve(ax,y_true,y_pred,pos_label=1,dataset_name=dataset,
                                       color=self._plot_color_list[n], plot_random=plot_random)

        handles, labels = ax.get_legend_handles_labels()
        labels = [l.replace(' (area = ', '. auROC=') for l in labels]
        labels = [l.replace(')', '') for l in labels]
        ax.legend(handles, labels, loc='lower right', bbox_to_anchor=(0.995, 0.005))
        ax.set_title('')
        ax.set_xlabel(ax.get_xlabel(), labelpad=0.9)
        ax.set_ylabel(ax.get_ylabel(), labelpad=0.8)
        ax.tick_params(axis='both', which='both', bottom=True, top=False,
                       left=True, right=False)
        return ax

    def _plot_pr_curves(self, ax):
        for n,dataset in enumerate(self._datasets_to_use):
            y_pred = self._data_dict[dataset]['prediction_df']
            y_true = self._data_dict['isCardio']
            plot_prevalence = True if n == 0 else False

            ax,pr_auc = plot_PR_curve(ax,y_true,y_pred,pos_label=1,dataset_name=dataset,
                                      color=self._plot_color_list[n], plot_prevalence=plot_prevalence)

        handles, labels = ax.get_legend_handles_labels()
        # labels = [l.replace(' (area = ', '. auROC=') for l in labels]
        # labels = [l.replace(')', '') for l in labels]
        # ax.legend(handles, labels, loc='top right', bbox_to_anchor=(0.995, 0.995))
        ax.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.995, 0.995))
        ax.set_xlabel(ax.get_xlabel(), labelpad=0.9)
        ax.set_ylabel(ax.get_ylabel(), labelpad=0.8)
        ax.tick_params(axis='both', which='both', bottom=True, top=False,
                       left=True, right=False)
        return ax

    def _plot_deciles(self, ax):
        for n, dataset in enumerate(self._datasets_to_use):
            if dataset == 'all features':
                # dataset = 'TCR features inc. pred. age gender'
                y_pred = self._data_dict[dataset]['prediction_df']
                y_true = self._data_dict['isCardio']
                ax, test_deciles = plot_deciles(ax, y_test=y_true, y_pred=y_pred, target_name='ACS', q=10,
                                            color=self._plot_color_list[n], show_fold=True,
                                            lw=0.5)
        ax.set_xlabel(ax.get_xlabel(), labelpad=0.9)
        ax.set_ylabel(ax.get_ylabel(), labelpad=0.8)
        ax.set_ylim(0, 80)
        ax.margins(y=0.9)
        ax.tick_params(axis='both', which='both', bottom=True, top=False,
                       left=True, right=False)
        return ax

    def _plot_shap(self, ax, fig):
        feature_df = self._data_dict['feature_df']
        feature_df.index.rename('index', inplace=True)
        nTopFeatures = 20
        jitter = 0.1
        sample_num = None
        scalingMethod = 'perc'
        target_name = 'ACS prediction'
        dataset = 'all features'
        y_pred = self._data_dict[dataset]['prediction_df']
        y_true = self._data_dict['isCardio']
        shap_df = self._data_dict[dataset]['shap_df']['isCardio'].drop('Age', axis=1)
        ax, fig, topN_features_abs_shap = gen_shap_summary_to_axes(ax, fig, target_name=target_name, nTopFeatures=nTopFeatures, shap_df=shap_df,
                             features_df=feature_df, jitter=jitter, sample_num=sample_num,
                             feature_name_list=None, scalingMethod=scalingMethod,
                             addColorBar=True, alpha=1)
        ax.set_ylabel(ax.get_ylabel(), labelpad=0.8)
        ax.tick_params(axis='both', which='both', bottom=True, top=False,
                       left=True, right=False)

        self._topN_features_abs_shap = topN_features_abs_shap

        return ax

    def _plot_shap_values_barplot(self, ax):
        ax.bar(range(len(self._topN_features_abs_shap)), self._topN_features_abs_shap, color='black')
        ax.set_xlim(0, len(self._topN_features_abs_shap))
        ax.set_ylim(0, 200)
        ax.tick_params(axis='y', which='both', left=False, right=False)
        ax.tick_params(axis='x', which='both', top=False)

        ax.set_xticks([])
        ax.set_xticklabels('')
        ax.set_xlabel('');
        ax.set_ylabel('abs mean\nshap value', labelpad=0.8)
        yticks=[0,100, 200]
        ax.set_yticks(yticks)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_position(('data', 0))

        return ax

    def _generate_fig_outline(self):
        self.fig = plt.figure(figsize=self._figsize, dpi=mpl.rcParams['figure.dpi'])

        ##add sub-figure letters and remove spines:
        self._gs1 = gridspec.GridSpec(1, 1)
        self._gs1.update(top=0.95, bottom=0.66, left=0.05, right=0.31)
        self._ax1 = self.fig.add_subplot(self._gs1[0, 0])

        self._gs2 = gridspec.GridSpec(1, 1)
        self._gs2.update(top=0.95, bottom=0.66, left=0.38, right=0.64)
        self._ax2 = self.fig.add_subplot(self._gs2[0, 0])

        self._gs3 = gridspec.GridSpec(1, 1)
        self._gs3.update(top=0.95, bottom=0.66, left=0.71, right=0.97)
        self._ax3 = self.fig.add_subplot(self._gs3[0, 0])

        self._gs5 = gridspec.GridSpec(1, 1)
        self._gs5.update(top=0.55, bottom=0.23, left=0.11, right=0.96)
        self._ax5 = self.fig.add_subplot(self._gs5[0, 0])

        # self._gs6 = gridspec.GridSpec(1, 1)
        # self._gs6.update(top=0.67, bottom=0.12, left=0.68, right=0.98)
        # self._ax6 = self.fig.add_subplot(self.gs6[0, 0])

        self._gs8 = gridspec.GridSpec(1, 1)
        self._gs8.update(top=0.09, bottom=0.02, left=0.11, right=0.927)
        self._ax8 = self.fig.add_subplot(self._gs8[0, 0])

    def plot_fig_2(self):
        self._get_data()
        self._set_fig2_definitions()
        self._generate_fig_outline()
        self._ax1 = self._plot_roc_curves(self._ax1)
        self._ax2 = self._plot_pr_curves(self._ax2)
        self._ax3 = self._plot_deciles(self._ax3)
        self._ax5 = self._plot_shap(self._ax5, self.fig)
        self._ax8 = self._plot_shap_values_barplot(self._ax8)

        for ax, let in zip([self._ax1, self._ax2, self._ax3, self._ax5, self._ax8], ['A', 'B', 'C', 'D', 'E']):
            ax.text(-0.1, 1.12, let, fontsize='xx-large', fontweight='bold', transform=ax.transAxes, ha='left', va='top')


if __name__ == "__main__":
    fig2_obj = Figure2(size_factor=1.3)
    fig2_obj.plot_fig_2()
    dpi = int(mpl.rcParams['figure.dpi'])
    fig2_obj.fig.savefig(FIGURE_DIR + 'figure2_new_dpi%s_final.png' %dpi, dpi=dpi)


