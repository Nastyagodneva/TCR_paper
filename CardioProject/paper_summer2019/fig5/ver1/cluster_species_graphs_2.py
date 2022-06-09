import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

genie_base_dir = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
data_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
figure_dir = genie_base_dir + 'Presentations and Manuscripts/CardioTCR paper/\
June2019/fig5_outcomes/calc_fig5_2/'
excel_files_dir = figure_dir + 'excel_files_new/'
cluster_species_df_file = 'any_outcome_min_shared_species01_max_shared_species1\
_min_shared_cluster03_max_shared_cluster1/cluster_species_phen_summary_df__fdr005top100rels.xlsx'




def _convert_species_names(species):
    '''
    this function convert species names to include only the species name and the species number.
    if species name is unknown, genus name is included and also 'species_unknown'
    new lines are introduced to species names to enable nice labeling of the nodes
    '''

    new_name = ''
    if 'unknown' in species.split('|')[-4]:
        part1 = species.split('|')[-5].split('__')[-1] + '|'
        part2 = '\nspecies_unknown|'
    else:
        part1 = species.split('|')[-4].split('__')[-1]
        part2 = ''
    if len(part1.split('_')) > 1:
        new_part1 = ''
        for n, part in enumerate(part1.split('_')):
            if n != len(part1.split('_')) - 1:
                link = '_\n'
            else:
                link = '|'
            new_name += part + link
        part1 = new_part1
    new_name += part1 + part2
    new_name += '\n' + species.split('|')[-1].split('__')[-1]

    return new_name


def edit_species_column(cluster_species_df):
    '''
    this function applies convert_species_names to the species column in cluster_species_df
    :param cluster_species_df:
    :return:
    '''
    cluster_species_df['species'] = cluster_species_df['species'].apply(lambda x: _convert_species_names(x))
    return cluster_species_df


def get_cluster_species_df():
    print ('loading and editing cluster species df')
    print ('files used: ', cluster_species_df_file)
    cluster_species_df = pd.read_excel(figure_dir + cluster_species_df_file)
    cluster_species_df = edit_species_column(cluster_species_df)
    return cluster_species_df





def add_all_nodes(G,cluster_species_df):

    for i, node in enumerate(cluster_species_df['cluster'].unique()):
        G.add_node(node)
    for i, node in enumerate(cluster_species_df['species'].unique()):
        G.add_node(node)
    return G

def add_all_edges(G,cluster_species_df):

    cluster_species_wide = pd.pivot_table(data=cluster_species_df, values='-log10p_MW_cluster_species',
                                          index='cluster', columns='species').fillna(9999)
    for i in cluster_species_wide.index:
        for j in cluster_species_wide.columns:
            # print cluster_species_wide.loc[i, j]
            if cluster_species_wide.loc[i, j] != 9999:
                G.add_edge(i, j, weight=cluster_species_wide.loc[i, j])
    return G

def gen_graph(G,cluster_species_df):

    #get all edges and nodes to G:
    print ('getting all nodes and edges to G...')
    G = add_all_nodes(G,cluster_species_df)
    G = add_all_edges(G,cluster_species_df)

    return G

######################################3

cluster_species_df = get_cluster_species_df()
G = nx.Graph()
G = gen_graph(G,cluster_species_df)

path = figure_dir + 'trial.graphml'
print (path)
nx.write_graphml(G, path)




