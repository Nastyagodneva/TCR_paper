import pandas as pd

GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
SAMPLE_LIST_DIR = GENIE_BASE_DIR + 'Sample files/BD lists/'
CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'


old_outcome_df_file = CARDIO_PHEN_DIR + 'outcomeDF.xlsx'
bd_reg_code_converter_file = DATA_DIR + 'TCR_tables_for_loader/helpers/BD_regcode_date_conversion_table.xlsx'
new_medical_file = CARDIO_PHEN_DIR + 'Medical_database_31-12-2019.xlsx'



class GenerateNewOutComeDF:

    def __init__(self,outcome_list=None):
        self.outcome_list = outcome_list
        pass


    @staticmethod
    def get_df(file_name):
        df = pd.read_excel(file_name)
        try:
            df = df.set_index('BD')
        except:
            print (file_name, ': couldnt set BD as index')
        return df


    def generate_new_outcome_df(self):

        self.old_outcome = self.get_df(old_outcome_df_file)
        self.bd_converter = self.get_df(bd_reg_code_converter_file)
        self.new_medical = self.get_df(new_medical_file)


        if self.outcome_list is None:
            self.outcome_list = self.old_outcome.columns.tolist()

        self.old_outcome_with_reg = pd.merge(self.old_outcome,pd.DataFrame(self.bd_converter['RegistrationCode']),how='left',
                                        left_index=True,right_index=True)
        for outcome in self.outcome_list:
            if outcome in self.old_outcome_with_reg.columns:
                self.old_outcome_with_reg.rename(columns={outcome: outcome+'_old'}, inplace=True)

        self.outcome_list = [col for col in self.outcome_list if col in self.new_medical.columns]
        self.new_outcome = pd.merge(self.old_outcome_with_reg.reset_index(),self.new_medical[self.outcome_list + ['RegistrationCode']], how='left',
                               left_on = 'RegistrationCode', right_on = 'RegistrationCode')
        self.new_outcome = self.new_outcome.set_index('BD')

        self.stats = pd.merge (pd.DataFrame(self.old_outcome.sum()), pd.DataFrame(self.new_outcome.sum()),
                          how='outer', left_index=True, right_index = True, suffixes=['old','new'])
        self.stats = self.stats.drop([x for x in self.stats.index if 'old' in x], axis=0)
        return


    def process_new_outcome(self):
        self.new_outcome = self.new_outcome.drop_duplicates()
        self.new_outcome = self.new_outcome.drop([x for x in self.new_outcome.columns if 'old' in x], axis=1)
        self.new_outcome['any_outcome'] = self.new_outcome.apply(lambda x: 1 if x.drop(['RegistrationCode']).sum()>0 else 0, axis=1)
        self.new_outcome.to_excel(CARDIO_PHEN_DIR + 'new_outcome_df.xlsx')


GenerateNewOutComeDF = GenerateNewOutComeDF()
GenerateNewOutComeDF.generate_new_outcome_df()
GenerateNewOutComeDF.process_new_outcome()

new_outcome = GenerateNewOutComeDF.new_outcome
final_stats = pd.merge(pd.DataFrame(GenerateNewOutComeDF.new_outcome.count()),
                       pd.DataFrame(GenerateNewOutComeDF.new_outcome.sum()),
                         how='inner',left_index=True, right_index=True,suffixes=['count','sum'])
