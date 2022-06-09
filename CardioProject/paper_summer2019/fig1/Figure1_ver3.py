from scipy import stats
from scipy.spatial.distance import braycurtis
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from ShaniBA.myplots import roundup
from ShaniBA.CardioProject.Figures.GeneralFigureFunctions import *
from ShaniBA.myplots import *
import datetime

time_a = datetime.datetime.now()
print ('now starting on: %s' %time_a)


####path definitions:
MyPath = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
FIG1_DIR = MyPath + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/figures/'
RED_RESULTS_DIR = MyPath + 'predictions2/'


##### sample lists:
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/PNP530', 'rb') as fp:
    PNP530 = pickle.load(fp)
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Sample files/BD lists/Cardio126', 'rb') as fp:
    Cardio126 = pickle.load(fp)
PNP530Cardio126 = PNP530 + Cardio126

##### general definitions:
pd.set_option('display.expand_frame_repr', False)
pd.set_option('display.width', 1000)
np.set_printoptions(precision=4, linewidth=200)
cdate = str(time.strftime("%d%m%Y"))

#### functions:
def set_fig1_definitions():
    params = {
   'axes.labelsize': 16,
   'font.size': 12,
   'legend.fontsize': 8,
   'xtick.labelsize': 14,
   'ytick.labelsize': 14,
   'text.usetex': False,
#    'figure.figsize': [m2inch(183), m2inch(247)],#[4.5, 4.5]
   'figure.dpi': 300,
   'xtick.direction':'out'}
    mpl.rcParams.update(params)
    mpl.rcParams['axes.linewidth'] = 0.5
    mpl.rcParams['axes.edgecolor'] = 'black'
    mpl.rcParams['axes.facecolor'] = 'white'
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['xtick.minor.pad'] = 4
    return


def get_real_and_pred_data(phen, real_data_path, pred_data_file, patients_also=True):
    ## get predicted and real data for healthy:
    real_data = pd.read_excel(real_data_path).set_index('BD')[phen].dropna()
    pred_data_PNP530 = pd.read_pickle(pred_data_file)[phen].dropna()

    # merge for healthy:
    pred_and_real_PNP530 = pd.merge(pd.DataFrame(real_data.rename(phen)),
                                    pd.DataFrame(pred_data_PNP530.rename('pred ' + phen).astype(float)), how='inner',
                                    left_index=True, right_index=True)
    print ('pred_and_real_PNP530.shape: ', pred_and_real_PNP530.shape)
    print pred_and_real_PNP530.head()

    if patients_also:
        ### get predicted and real data for patients:
        pred_and_real_Cardio126 = \
        pd.read_pickle(PRED_RESULTS_DIR + 'Cardio126_phens_basedOnHealthy/pred_and_real_%s.pkl' % phen)[
            ['average pred', phen]]
        pred_and_real_Cardio126 = pred_and_real_Cardio126.rename(columns={'average pred': 'pred ' + phen})
        print ('pred_and_real_Cardio126: ', pred_and_real_Cardio126.shape)
        pred_and_real_Cardio126.head()
    else:
        pred_and_real_Cardio126 = None

    return pred_and_real_PNP530, pred_and_real_Cardio126


def plot_phen_prediction_ROC_PR_fig1(ax, phen, pred_and_real):
    y = pd.DataFrame(pred_and_real[phen])
    y_pred_df = pd.DataFrame(pred_and_real['pred ' + phen].rename('pred_proba'))
    # #plot:
    pos_label = 1
    ax, inset_axes, roc_auc, pr_auc, prevalence = plot_ROC_PR_AUC(y, y_pred_df, ax, color1='darkred', color2='grey',
                                                                  ticklabelsize=mpl.rcParams['xtick.labelsize'],
                                                                  textsize=mpl.rcParams['font.size'],
                                                                  labelsize=mpl.rcParams['axes.labelsize'],
                                                                  add_texts=False)
    ax.annotate('ROC AUC=%s\nPR AUC=%s' % (round(roc_auc, 3), round(pr_auc, 2)), xy=(0.04, 0.98),
                xycoords='axes fraction',
                fontsize=mpl.rcParams['font.size'], xytext=(0, 0), textcoords='offset points', fontweight='bold',
                ha='left', va='top')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    inset_axes.text(0.01, prevalence - 0.01, 'Prevalence=%s' % prevalence, transform=inset_axes.transAxes, ha='left',
                    va='top',
                    fontsize=mpl.rcParams['font.size'] - 2)

    return ax


def plot_phen_prediction_correlation_fig1(ax, phen, pred_and_real, outlierSTD=None, toAnnotate=True):
    merged = pred_and_real.rename(columns={phen: 'Real', 'pred ' + phen: 'Predicted'}).dropna(how='any')
    print ('merged.head(): ', merged.head())

    # # plot data:
    x = 'Real';
    y = 'Predicted'
    if outlierSTD is not None:
        x_lim1 = merged[x].mean() - outlierSTD * merged[x].std()
        x_lim2 = merged[x].mean() + outlierSTD * merged[x].std()
        print ('value limits are: ', x_lim1, x_lim2)
        print ('merged.shape before outlier removal= ', merged.shape)
        merged = merged[(merged[x] > x_lim1) & (merged[x] < x_lim2)]
        print ('merged.shape after outlier removal= ', merged.shape)

    merged.plot(x, y, ax=ax, kind='scatter', alpha=0.5, c='darkred', s=60)
    ax.plot(np.unique(merged[x]), np.poly1d(np.polyfit(merged[x], merged[y], 1))(np.unique(merged[x])), c='black',
            linewidth=1)
    r, p = MyPearsonr(merged[x], merged[y])
    print r, p
    if toAnnotate:
        ax.annotate('r=%s\np=%.1E' % (round(r, 2), p), xy=(0.04, 0.98), xycoords='axes fraction',
                    fontsize=mpl.rcParams['font.size'],
                    xytext=(0, 0), textcoords='offset points', fontweight='bold',
                    ha='left', va='top')
    return ax, r, p


def compare_healthy_patients_predictions_grouped(phen,
        real_data_path, pred_data_file, bin_edges, bin_names, ax=None, isBinary=False):

    if isBinary:
        bin_edges = [-1, 0, 1]

    ## get predicted and real data for healthy:

    real_data = pd.read_excel(real_data_path).set_index('BD')[phen].dropna()

    pred_data_PNP530 = pd.read_pickle(pred_data_file)[phen].dropna()

    # merge for healthy:
    pred_and_real_PNP530 = pd.merge(pd.DataFrame(real_data.rename(phen)),
                                    pd.DataFrame(pred_data_PNP530.rename('pred ' + phen).astype(float)), how='inner',
                                    left_index=True, right_index=True)
    print ('pred_and_real_PNP530.shape: ', pred_and_real_PNP530.shape)
    print pred_and_real_PNP530.head()

    ### get predicted and real data for patients:
    pred_and_real_Cardio126 = \
    pd.read_pickle(PRED_RESULTS_DIR + 'Cardio126_phens_basedOnHealthy/pred_and_real_%s.pkl' % phen)[['average pred', phen]]
    pred_and_real_Cardio126 = pred_and_real_Cardio126.rename(columns={'average pred': 'pred ' + phen})
    print ('pred_and_real_Cardio126: ', pred_and_real_Cardio126.shape)
    pred_and_real_Cardio126.head()

    ##plot:
    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))

    values_dict = {}

    for item in [('Healthy', pred_and_real_PNP530, 'grey', 'lightgrey'),
                 ('Patients', pred_and_real_Cardio126, 'darkred', 'red')]:
        df = item[1]
        df_name = item[0]
        color_real = item[2]
        color_pred = item[3]

        #     print df.columns.tolist()
        if isBinary:
            fig_range = range(2)
        else:
            fig_range = range(1, len(bin_edges))

        #         print df_name
        #         print 'df[phen]:'
        #         print df[phen]

        df['Binned ' + phen] = pd.cut(df[phen], bin_edges, labels=fig_range)
        grouped = df.groupby('Binned ' + phen).agg([np.mean, np.std, np.size])
        grouped.columns = [' '.join(col).strip() for col in grouped.columns.values]

        for n in fig_range:
            values_dict[df_name + str(n)] = df[(df[phen] >= bin_edges[n - 1]) & (df[phen] < bin_edges[n])][
                'pred ' + phen].tolist()

        print df_name
        print grouped

        if not isBinary:
            ax.errorbar(fig_range, grouped[phen + ' mean'], yerr=grouped[phen + ' std'], fmt='_', color=color_real, alpha=1,
                        linewidth=3, markersize=20,
                        label=df_name + ' (real)')
        ax.errorbar(fig_range, grouped['pred ' + phen + ' mean'], yerr=grouped['pred ' + phen + ' std'], fmt='o',
                    color=color_pred, alpha=0.5, linewidth=4, markersize=8,
                    label=df_name + ' (pred.)')
        ax.set_xticks(fig_range)
        ax.set_xticklabels(bin_names, rotation=45, fontsize='medium')
    min_values = [];
    max_values = []
    for n in fig_range:
        data1 = values_dict['Healthy' + str(n)]
        data2 = values_dict['Patients' + str(n)]
        #         min_values.append(np.min(data1+data2))
        #         max_values.append(np.max(data1+data2))
        MW_s, MW_p = mannwhitneyu(data1, data2)
        print ('MW test: %s bin ' % phen, bin_names[n - 1], MW_s, MW_p)

        t_s, t_p = ttest_ind(data1, data2)
        print ('t test: %s bin ' % phen, bin_names[n - 1], t_s, t_p)
        print ''

        ax.legend(loc='upper left')
        ax.set_ylabel(phen + '', fontsize='large')
        ax.set_xlabel(phen + ' Group', fontsize='large')

    #     ymin=adjusted_rounddown(np.min(min_values))
    #     ymax=adjusted_roundup(np.min(max_values))


    xmin = fig_range[0] - 0.5;
    xmax = fig_range[-1] + 0.5

    ax.set_xlim(xmin, xmax)
    #     ax.set_ylim(ymin,ymax)

    return ax



###### plotting figure 1: ##########
fig = plt.figure(figsize=(16, 12))
gs0 = gridspec.GridSpec(3, 3, wspace=0.3, hspace=0.3)

ax = plt.gca()

x1 = (1. / 3 - 0.3 / 4) / 2
x2 = 0.5
x3 = 1 - x1

plt.text(x1, 1.01, 'Age', ha='center', va='bottom', fontsize='xx-large', fontweight='bold', transform=ax.transAxes)
plt.text(x2, 1.01, 'Gender', ha='center', va='bottom', fontsize='xx-large', fontweight='bold', transform=ax.transAxes)
plt.text(x3, 1.01, 'HbA1C', ha='center', va='bottom', fontsize='xx-large', fontweight='bold', transform=ax.transAxes)

plt.text(-0.075, x3, 'Healthy', ha='left', va='center', fontsize='xx-large', fontweight='bold', transform=ax.transAxes,
         rotation=90)
plt.text(-0.075, x2, 'Patients', ha='left', va='center', fontsize='xx-large', fontweight='bold', transform=ax.transAxes,
         rotation=90)
plt.text(-0.075, x1, 'Comparison', ha='left', va='center', fontsize='xx-large', fontweight='bold',
         transform=ax.transAxes, rotation=90)

remove_spines(removeFigBorders=True)

real_data_path = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/phenotypes_byBD/\
PNP530Cardio126_Age_Gender_Male_HbA1C_SmokingStatus_Yes_BMI_HDL_Smoking_ever.xlsx'
pred_data_file = PRED_RESULTS_DIR + 'PNP530_majorPhenotypes/\
AgeGenderBMIHbA1CSmokingCholesterol_XGB_randomSearch_25_TCRfeatures_optByExpVar/predictions_df.pkl'

pred_and_real_PNP530_age, pred_and_real_Cardio126_age = get_real_and_pred_data('Age', real_data_path, pred_data_file,
                                                                               patients_also=True)
pred_and_real_PNP530_gender, pred_and_real_Cardio126_gender = get_real_and_pred_data('Gender_Male', real_data_path,
                                                                                     pred_data_file, patients_also=True)
pred_and_real_PNP530_HbA1C, pred_and_real_Cardio126_HbA1C = get_real_and_pred_data('HbA1C', real_data_path,
                                                                                   pred_data_file, patients_also=True)

### A:
ax1 = plt.Subplot(fig, gs0[0, 0])
fig.add_subplot(ax1)

ax2 = plt.Subplot(fig, gs0[0, 1])
fig.add_subplot(ax2)

ax3 = plt.Subplot(fig, gs0[0, 2])
fig.add_subplot(ax3)

ax4 = plt.Subplot(fig, gs0[1, 0])
fig.add_subplot(ax4)

### B:
ax5 = plt.Subplot(fig, gs0[1, 1])
fig.add_subplot(ax5)

ax6 = plt.Subplot(fig, gs0[1, 2])
fig.add_subplot(ax6)

ax7 = plt.Subplot(fig, gs0[2, 0])
fig.add_subplot(ax7)

ax8 = plt.Subplot(fig, gs0[2, 1])
fig.add_subplot(ax8)

# plot correlations for healthy predictions - Age, HbA1C:
phen_list = ('Age', 'HbA1C')
ax_healthy_list = [ax1, ax3]
ax_patient_list = [ax4, ax6]
pred_and_real_healthy_list = [pred_and_real_PNP530_age, pred_and_real_PNP530_HbA1C]
pred_and_real_patient_list = [pred_and_real_Cardio126_age, pred_and_real_Cardio126_HbA1C]
min_val_Age = adjusted_rounddown(
    np.min([pred_and_real_PNP530_age.min() * 0.98, pred_and_real_Cardio126_age.min() * 0.98]))
max_val_Age = adjusted_roundup(
    np.max([pred_and_real_PNP530_age.max() * 1.02, pred_and_real_Cardio126_age.max() * 1.02]))
min_val_HbA1C = adjusted_rounddown(
    np.min([pred_and_real_PNP530_HbA1C.min() * 0.98, pred_and_real_Cardio126_HbA1C.min() * 0.98]))
max_val_HbA1C = adjusted_roundup(
    np.max([pred_and_real_PNP530_HbA1C.max() * 1.02, pred_and_real_Cardio126_HbA1C.max() * 1.02]))

for n in range(2):
    phen = phen_list[n]
    print phen
    pred_and_real_healthy = pred_and_real_healthy_list[n]
    pred_and_real_patient = pred_and_real_patient_list[n]
    ax_healthy = ax_healthy_list[n]
    ax_patient = ax_patient_list[n]

    ax_healthy, r, p = plot_phen_prediction_correlation_fig1(ax_healthy, phen, pred_and_real_healthy, outlierSTD=None,
                                                             toAnnotate=True)
    ax_patient, r, p = plot_phen_prediction_correlation_fig1(ax_patient, phen, pred_and_real_patient, outlierSTD=None,
                                                             toAnnotate=True)

    if n == 0:
        ax_healthy.set_xlim(min_val_Age, max_val_Age)
        ax_healthy.set_ylim(min_val_Age, max_val_Age)
        ax_patient.set_xlim(min_val_Age, max_val_Age)
        ax_patient.set_ylim(min_val_Age, max_val_Age)
    else:
        ax_healthy.set_xlim(4, 12)
        ax_healthy.set_ylim(5, 6)
        ax_patient.set_xlim(4, 12)
        ax_patient.set_ylim(5, 6)

# plot ROC_AUC curves for healthy and patients-Gender
ax4 = plot_phen_prediction_ROC_PR_fig1(ax2, 'Gender_Male', pred_and_real_PNP530_gender)
ax5 = plot_phen_prediction_ROC_PR_fig1(ax5, 'Gender_Male', pred_and_real_Cardio126_gender)

## plot predition comparison for Age and Gender:
phen_list2 = ['Age', 'Gender_Male']
ax_list = [ax7, ax8]
# bin_edges_list = [[36, 41, 46, 51, 56, 61, 66, 71], [-1, 0, 1]]
# bin_names_list = [['36-40', '41-45', '46-50', '51-55', '56-60', '61-65', '66-70'], ['Females', 'Males']]
bin_edges_list = [[36, 46, 56, 66, 76], [-1, 0, 1]]
bin_names_list = [['36-45', '46-55', '56-66', '66-75'], ['Females', 'Males']]
binary_list = [False, True]

for n in range(2):
    ax = ax_list[n]
    phen = phen_list2[n]
    bin_edges = bin_edges_list[n]
    bin_names = bin_names_list[n]
    isBinary = binary_list[n]

    print ('phen: ', phen)
    print ('bin_edges ', bin_edges)
    print ('bin_names ', bin_names)

    ax = compare_healthy_patients_predictions_grouped(phen, real_data_path, pred_data_file, bin_edges, bin_names, ax=ax,
                                                      isBinary=False)

# plt.tight_layout()
plt.show()
print 'now saving'
fig.savefig(FIG1_DIR + 'Fig1_ver3.png')
print 'saved'
time_b = datetime.datetime.now()
print ('script running started on %s' %time_a)
print ('script running ended on %s' %time_b)