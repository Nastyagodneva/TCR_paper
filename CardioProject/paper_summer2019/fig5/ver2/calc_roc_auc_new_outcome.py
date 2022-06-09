from sklearn import metrics
import pandas as pd
import cPickle as pickle
import numpy as np
import time
import os

from ShaniBA.CardioProject.Feature_phenotype_functions import add_corrected_pValues

cdate = time.strftime("%d%m%Y")


GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
SAMPLE_LIST_DIR = GENIE_BASE_DIR + 'Sample files/BD lists/'
CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5_updated_phens_repeatrun/'
PRED_RESULTS_DIR='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/predictions2/'

outcome_file = CARDIO_PHEN_DIR + 'new_outcome_df.xlsx'
pred_proba_file =PRED_RESULTS_DIR+'isCardio/XGB_randomSearch_25_byRepFeatPCA10RelsVDJnocorr0999AgeGender/\
predictions_df.pkl'

cluster_data_file = DATA_DIR + 'TCR_seqs_clusters/newX_onlySeqs_025_085_noNans.dat'


if __name__ == "__main__":
    cluster_data = pd.read_pickle(cluster_data_file)
    top100_fig3_file = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig3_balanced_comparison/calc_fig3/excel_files_new/top100_clusters_with_annot.xlsx'
    top100_fig3 = (pd.read_excel(top100_fig3_file))['cluster'].tolist()
    top100_fig3_cluster_data = cluster_data.loc[:,top100_fig3].fillna(0)
    top100_fig3_cluster_data.to_excel(DATA_DIR + 'TCR_seqs_clusters/top100_fig3_cluster_data.xlsx')
#     top100_fig5_file = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
# June2019/fig5_outcomes/calc_fig5_updated_phens/excel_files_new/top100_clusters_with_annot.xlsx'
#     top100_fig5 = (pd.read_excel(top100_fig5_file))['cluster'].tolist()
#     top100_fig5_cluster_data = cluster_data.loc[:,top100_fig5].fillna(0)
#     top100_fig5_cluster_data.to_excel(DATA_DIR + 'TCR_seqs_clusters/top100_fig5_cluster_data.xlsx')
    top100_fig5_file_hosp = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5_updated_phens/ver1_excel_files_new_cv_hospitalization/top100_clusters_with_annot.xlsx'
    top100_fig5_hosp = (pd.read_excel(top100_fig5_file_hosp))['cluster'].tolist()
    top100_fig5_cluster_data_hosp = cluster_data.loc[:, top100_fig5_hosp].fillna(0)
    top100_fig5_cluster_data_hosp.to_excel(DATA_DIR + 'TCR_seqs_clusters/top100_fig5_cluster_data_hosp.xlsx')
    with open(SAMPLE_LIST_DIR + 'cv_hospt_shared05_clusters.pkl', 'rb') as fp:
        cv_hospt_shared05_clusters = pickle.load(fp)
    cv_hospt_shared05 = cluster_data.loc[:, cv_hospt_shared05_clusters].fillna(0)
    cv_hospt_shared05.to_excel(DATA_DIR + 'TCR_seqs_clusters/cv_hospt_shared05_fig5.xlsx')


################################################################
#parameters to change - part 1:
feature_file_name = 'cv_hospt_shared05'

################################################################

if feature_file_name == 'default':
    feature_file = DATA_DIR + 'TCRfeatures/TCRfeaturesNoT_PCA10percShared_withRels_VDJ_noCorr0-999_nanFilled_noConsts_Cardio126.xlsx'
elif feature_file_name == 'top100_fig3':
    feature_file = DATA_DIR + 'TCR_seqs_clusters/top100_fig3_cluster_data.xlsx'
elif feature_file_name == 'top100_fig5':
    feature_file = DATA_DIR + 'TCR_seqs_clusters/top100_fig5_cluster_data.xlsx'
elif feature_file_name == 'top100_fig5_hosp':
    feature_file = DATA_DIR + 'TCR_seqs_clusters/top100_fig5_cluster_data_hosp.xlsx'
elif feature_file_name == 'cv_hospt_shared05':
    feature_file = DATA_DIR + 'TCR_seqs_clusters/cv_hospt_shared05_fig5.xlsx'

################################################################
#parameters to change - part 2:
outcome_column_list=['CV hospitalization including chest pain']
n_permut=100
roc_abs_thresh = 0.2

file_name = FIGURE_DIR + 'permut_dir/feature_outcome_roc_abs_%sperms_t%s_featurefile_%s_%s_3.xlsx' %(
    n_permut, roc_abs_thresh,feature_file_name, cdate)
file_name = file_name.replace('.','').replace('xlsx','.xlsx')
if len(outcome_column_list) == 1:
    outcome_name = '_'.join(outcome_column_list[0].split( )[:2])
    file_name = file_name.replace('.xlsx','_%s.xlsx' %outcome_name)
################################################################



class AucCalculator:

    def __init__(self, feature_value_list, outcome_value_list,ax=None,plot_color='orange'):
        self.feature_value_list = feature_value_list
        self.outcome_value_list = outcome_value_list
        self.ax = ax
        self.plot_color = plot_color
        

    def _calc_auc(self):
        self.fpr, self.tpr, thresholds = metrics.roc_curve(self.outcome_value_list, self.feature_value_list)
        self.roc_auc = metrics.auc(self.fpr, self.tpr)
                
        return self.roc_auc


    def _plot_roc_auc(self, lw=1):
        self.fpr, self.tpr, thresholds = metrics.roc_curve(self.outcome_value_list, self.feature_value_list)
        self.roc_auc = metrics.auc(self.fpr, self.tpr)
        label='auc=%0.3f' % self.roc_auc
        self.ax.plot(self.fpr, self.tpr, color=self.plot_color, lw=lw, label=label)
        self.ax.plot([0, 1], [0, 1], color='black', lw=lw, linestyle='--')
        self.ax.set_xlim([0.0, 1.0])
        self.ax.set_ylim([0.0, 1.05])
        self.ax.set_xlabel('FPR')
        self.ax.set_ylabel('TPR')
        # self.ax.legend(loc=(0.45,-0.05))
        self.ax.text(0.99, 0.01, label, ha='right', va='bottom',
             transform=self.ax.transAxes)


#class that manges everything:
# get list of features
# get list of outcomes
# loop over pairs and run the auc calculation
# if auc>threshold, uses the permutated data and calculate real p-value

class AucCalcManager:
    def __init__(self, outcome_column_list=['Unplanned PCI', 'CV hospitalization including chest pain',
                'any_outcome'],file_name = None, n_permut=1000000,roc_abs_thresh = 0.2):
        self.feature_df = self._get_feature_df
        self.outcome_df = pd.read_excel(outcome_file).drop('RegistrationCode', axis=1).set_index('BD')
        self.outcome_column_list = outcome_column_list
        self.result_df = pd.DataFrame()
        self.count = 0
        self.file_name = file_name
        self.permut_dir = FIGURE_DIR + 'permut_dir/'
        if not os.path.isdir(self.permut_dir):
            os.makedirs(self.permut_dir)
        self.n_permut = n_permut
        # self.roc_thresh = roc_thresh
        self.roc_abs_thresh = roc_abs_thresh


    @property
    def _get_feature_df(self):
        feature_df = pd.read_excel(feature_file)
        try:
            feature_df = feature_df.set_index('BD')
        except:
            pass
        is_cardio_pred_proba = pd.read_pickle(pred_proba_file)
        feature_df_new = pd.merge(feature_df,is_cardio_pred_proba,how='left',
                                  left_index=True,right_index=True)
        return feature_df_new


    def _gen_feature_and_outcome_value_list(self, feature, outcome, outcome_df=None):
        if outcome_df is None:
            outcome_df = self.outcome_df
        df = pd.merge(pd.DataFrame(self.feature_df.loc[:,feature]),
                      pd.DataFrame(outcome_df.loc[:,outcome]),
                      how='inner', left_index=True, right_index=True)
        df = df.dropna(how='any')
        feature_value_list = df.loc[:,feature].tolist()
        outcome_value_list = df.loc[:,outcome].tolist()
        return feature_value_list, outcome_value_list


    def _permut_outcomes(self, outcome):
        outcome_permut_df = pd.DataFrame()
        print ('performing permutations for outcome %s, %s permutations' %(outcome, self.n_permut))
        ind = self.outcome_df.index
        for n in range(self.n_permut):
            if n%10000 == 0: print n
            outcome_permut_df[n] = self.outcome_df[outcome].sample(frac=1).reset_index(drop=True)

        outcome_permut_df.index = ind
        outcome_permut_df.to_pickle(self.permut_dir + outcome + '_%sperms.dat' %self.n_permut)
        return outcome_permut_df


    def _get_permut_df(self, outcome):
        try:
            outcome_permut_df = pd.read_pickle(self.permut_dir + outcome + '_%sperms.dat' %self.n_permut)
            print ('loaded permut df for outcome %s' %outcome)
        except:
            print ('permut_df for outcome %s was not found, performing %s permutations now...' %(outcome, self.n_permut))
            outcome_permut_df = self._permut_outcomes(outcome)
        return outcome_permut_df

    def _permut_loop(self, curr_n_permut,feature,outcome,roc_abs):

        try:
            with open(self.permut_dir + '%s_%s_%spermut_result_list.pkl' % (feature, outcome, self.n_permut),
                      'rb') as fp:
                result_list = pickle.load(fp)
                print ('loaded result list from pickle')
        except:
            result_list = []
            for perm in range(curr_n_permut):
                if perm%10000 == 0: print perm

                #get value lists:
                feature_value_list, outcome_value_list = self._gen_feature_and_outcome_value_list(feature=feature,
                                                        outcome=outcome, outcome_df=pd.DataFrame(self.permut_df.iloc[:,perm].rename(outcome)))
                #calculate roc_auc and collect:
                roc_auc_per = AucCalculator(feature_value_list=feature_value_list,
                                        outcome_value_list=outcome_value_list)._calc_auc()
                roc_abs_per = np.abs(0.5 - roc_auc_per)
                result_list.append(roc_abs_per)

        bigger_auc = [x for x in result_list if x > roc_abs]
        actual_num_perms = len(result_list)
        p_value = float(len(bigger_auc)) / actual_num_perms

        return result_list, p_value, actual_num_perms

    def _calc_real_p(self, outcome, feature, roc_abs):

        #permutation loop:
        n_perms = len(self.permut_df.columns)
        print ('now starting to calculate real p-value by %s permutations over %s-%s pair' %(n_perms,
                                                                                             outcome,feature))
        result_list=[]
        if self.n_permut >= 100:

            result_list, p_value, actual_num_perms = self._permut_loop(curr_n_permut=100,feature=feature,outcome=outcome,
                                            roc_abs=roc_abs)
            if p_value < 0.02:
                print ('preliminary p-value = %s, continue with permutations...' %p_value)
                result_list, p_value, actual_num_perms = self._permut_loop(curr_n_permut=self.n_permut, feature=feature, outcome=outcome,
                                                         roc_abs=roc_abs)
                #save result to pickle
                if self.n_permut > 10000:
                    with open(self.permut_dir + '%s_%s_%spermut_result_list.pkl' %(feature,outcome,self.n_permut), 'wb') as fp:
                        pickle.dump(result_list,fp)
        else:
            result_list, p_value, actual_num_perms = self._permut_loop(curr_n_permut=self.n_permut, feature=feature,
                                                                       outcome=outcome,  roc_abs=roc_abs)

        return p_value, actual_num_perms


    def execute_calculation(self):
        print ('number of outcome: ', len(self.outcome_column_list))
        print ('number of features: ', len(self.feature_df.columns))

        for outcome in self.outcome_column_list:

            self.permut_df = self._get_permut_df(outcome=outcome)

            for feature in self.feature_df.columns:


                feature_value_list, outcome_value_list = self._gen_feature_and_outcome_value_list(feature=feature,
                                                                                                  outcome=outcome)
                roc_auc = AucCalculator(feature_value_list=feature_value_list,
                                        outcome_value_list=outcome_value_list)._calc_auc()

                roc_abs = np.abs(0.5-roc_auc)

                #run 100 permutations for all. (if request to run permutations only for roc_abs above the
                #threshold, comment the following line and uncomment the ones that follow.
                p_value, actual_num_perms = self._calc_real_p(outcome=outcome, feature=feature, roc_abs=roc_abs)

                # if roc_abs >= self.roc_abs_thresh:
                #     p_value, actual_num_perms = self._calc_real_p(outcome=outcome, feature=feature,roc_abs=roc_abs)
                # else:
                #     p_value = 1; actual_num_perms = 0

                print ('feature: ', feature, ' outcome: ', outcome)
                print ( 'roc_auc: ', roc_auc, ' roc_abs: ', roc_abs, ' p_value: ', p_value)
                self.result_df.loc[self.count, 'feature'] = feature
                self.result_df.loc[self.count, 'outcome'] = outcome
                self.result_df.loc[self.count, 'roc_auc'] = roc_auc
                self.result_df.loc[self.count, 'roc_abs'] = roc_abs
                self.result_df.loc[self.count, 'p_value'] = p_value
                self.result_df.loc[self.count, 'actual_num_perms'] = actual_num_perms
                self.count +=1

        self.result_df.sort_values(by='roc_abs', ascending=False, inplace=True)
        self.result_df = add_corrected_pValues(resultsDF=self.result_df,pValueColumn='p_value',nTests=len(self.result_df),
                                               FDR=0.25)

        if self.file_name is not None:
            try:
                self.result_df.to_excel(self.file_name)

            except:
                print ('couldnt save results to file')
        return self.result_df


if __name__ == "__main__":
    AucCalcManager = AucCalcManager(file_name=file_name,n_permut=n_permut,
                                    roc_abs_thresh=roc_abs_thresh,outcome_column_list=outcome_column_list)
    result_df = AucCalcManager.execute_calculation()
