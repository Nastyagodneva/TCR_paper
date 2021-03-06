import pandas as pd

#get network file:
GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
SAMPLE_LIST_DIR = GENIE_BASE_DIR + 'Sample files/BD lists/'
CARDIO_PHEN_DIR='/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_real_data/CardioSamples/phenotypicData/'
FIGURE_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5_updated_phens/'
fig5_FDR010101_table_file = FIGURE_DIR + '/excel_files_new_cv_hospitalization2/cluster_species_hm_fdr%s%s_longform.xlsx'
fig5_FDR010101_table = pd.read_excel(fig5_FDR010101_table_file)

#filter for sequences
fig5_cluster_list = ['CASSLQQGNTEAFF', 'CASSPTGSETQYF',
        'CASSLGWGGEQYF', 'CASSLETGVYEQYF', 'CASSPTQDYGYTF']

predictor_clusters_table = fig5_FDR010101_table[fig5_FDR010101_table['cluster'].isin(fig5_cluster_list)]
predictor_clusters_table_nodups = predictor_clusters_table.drop_duplicates(subset = [col for col in predictor_clusters_table.\

                                                                                     columns if 'species' not in col])
predictor_clusters_table_nodups.to_excel(FIGURE_DIR + 'network_files/predictor_clusters_table_nodups.xlsx')
related_species = predictor_clusters_table_nodups['species'].unique().tolist()
predictor_cluster_extended_table = fig5_FDR010101_table[fig5_FDR010101_table['species'].isin(related_species)]
predictor_cluster_extended_table.to_excel(FIGURE_DIR + 'network_files/predictor_cluster_extended_table.xlsx')


#get the associated species


#get the associated sequences

