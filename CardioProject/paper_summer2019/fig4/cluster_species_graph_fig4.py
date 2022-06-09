from paper_summer2019.fig5.ver1.cluster_species_graphs import ClusterSpeciesGraph
import matplotlib.pyplot as plt
import pandas as pd

#### definitions:
n_top_clusters = 10
genie_base_dir = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
fig3_dir = genie_base_dir + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig3_balanced_comparison/calc_fig3/'
figure_dir = genie_base_dir + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig4_graph/'
fig3_cluster_species_df_file = 'isCardio_min_shared_species01_max_shared_species1_min_shared_cluster025_max_shared_cluster1/\
cluster_species_phen_summary_df__025100025.xlsx'
fig4_cluster_species_df_file = 'cluster_species_phen_summary_df_top%sclusters_01species.xlsx' %n_top_clusters

#### generate data file:

def gen_data_file():
    top = (pd.read_excel(fig3_dir + 'excel_files/top100_clusters_with_annot.xlsx')\
        .set_index('cluster')).index.tolist()[:n_top_clusters]

    df = pd.read_excel(fig3_dir + fig3_cluster_species_df_file).set_index('cluster')
    df_top = df.loc[top,:]
    df_top_species01 = df_top[df_top['MW_p_corrPval_FDR0.1_species_phen']<=0.1]

    df_top_species01 = df_top_species01.sort_values(by = ['-log10p_MW_cluster_phen',
             '-log10p_MW_species_phen','-log10p_MW_cluster_species'], ascending=False)

    df_top_species01.to_excel(figure_dir + fig4_cluster_species_df_file)

    return

## draw graph:
fig,ax = plt.subplots(figsize = (12,25))

produce_data_file = True

if produce_data_file:
    gen_data_file()


fig4_graph = ClusterSpeciesGraph(fig=fig, ax=ax,
                                 cluster_species_df_file=fig4_cluster_species_df_file,
                                 figure_dir=figure_dir,
                                 y_position_factor_species = 1.5, y_position_factor_cluster = 6,
                                 node_size_factor = 1200,
                                 species_x_position=70,
                                 label_font_size = 6, jitter=True,edge_width_factor = 0.3)

fig4_graph.gen_graph()
