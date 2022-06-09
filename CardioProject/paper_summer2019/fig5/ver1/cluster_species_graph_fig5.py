from paper_summer2019.fig5.ver1.cluster_species_graphs import ClusterSpeciesGraph
import matplotlib.pyplot as plt


fig,ax = plt.subplots(figsize = (12,30))

genie_base_dir = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
figure_dir = genie_base_dir + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5_2/'
cluster_species_df_file = 'any_outcome_min_shared_species01_max_shared_species1\
_min_shared_cluster03_max_shared_cluster1/cluster_species_phen_summary_df__fdr005top50rels.xlsx'


fig5_graph = ClusterSpeciesGraph(fig=fig, ax=ax,
                                 cluster_species_df_file=cluster_species_df_file,
                                 figure_dir=figure_dir,
                                 y_position_factor_species = 10000, y_position_factor_cluster = 10000,
                                 node_size_factor = 750,
                                 species_x_position=50,
                                 label_font_size = 5, jitter=True,edge_width_factor = 0.3, paper_figure_num=5)

fig5_graph.gen_graph()
