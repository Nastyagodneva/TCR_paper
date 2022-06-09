import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from os import makedirs
from os.path import isdir

class ClusterSpeciesGraph():

    def __init__(self, fig, ax, cluster_species_df_file,figure_dir, species_x_position=50,
                 pos_assoc_color='thistle', neg_assoc_color='lightblue',
                 y_position_factor_cluster=1, y_position_factor_species=1, node_size_factor=2000,
                 edge_width_factor = 2, label_font_size=8, jitter=True):

        self.fig = fig
        self.ax = ax
        self.data_dir = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
        self.figure_dir = figure_dir
        self.excel_files_dir = self.figure_dir + 'excel_files_new/'
        self.cluster_species_df_file = cluster_species_df_file
        self.cluster_species_df = self.get_cluster_species_df
        self.species_x_position = species_x_position
        self.pos_assoc_color = pos_assoc_color
        self.neg_assoc_color = neg_assoc_color
        self.G = self.set_graph
        self.y_position_factor_cluster = y_position_factor_cluster
        self.y_position_factor_species = y_position_factor_species
        self.node_size_factor = node_size_factor
        self.edge_width_factor = edge_width_factor
        self.label_font_size = label_font_size
        self.jitter = jitter


    @property
    def set_graph(self):
        return nx.Graph()

    @staticmethod
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

    @classmethod
    def edit_species_column(cls, cluster_species_df):
        '''
        this function applies convert_species_names to the species column in cluster_species_df
        :param cluster_species_df:
        :return:
        '''
        cluster_species_df['species'] = cluster_species_df['species'].apply(lambda x: cls._convert_species_names(x))
        return cluster_species_df

    @property
    def get_cluster_species_df(self):
        print ('loading and editing cluster species df')
        print ('files used: ', self.cluster_species_df_file)
        cluster_species_df = pd.read_excel(self.figure_dir + self.cluster_species_df_file)
        cluster_species_df = self.edit_species_column(cluster_species_df)
        return cluster_species_df


    def set_node_size(self, node, node_type):
        '''
        :param node: cluster/species variable
        :param node_type: string: 'cluster' or 'species'
        '''
        size = self.cluster_species_df.set_index(node_type).loc[node, '-log10p_MW_%s_phen' %node_type].max()
        return size

    def set_node_color(self, node, node_type):
        '''
        :param node: cluster/species variable
        :param node_type: string: 'cluster' or 'species'
        '''
        if self.cluster_species_df.set_index(node_type).loc[node,'%s_pos_associated_with_phen' %node_type].max() ==1:
            color = self.pos_assoc_color
        else:
            color = self.neg_assoc_color
        return color

    def set_node_x_position(self,node_type):
        '''
        :param node_type: string: 'cluster' or 'species'
        '''

        if not self.jitter:
            jitter_by = 0
        else:
            jitter_by = np.random.randint(-1 * np.floor(self.species_x_position / 5),
                                          np.floor(self.species_x_position / 5))
        if node_type == 'cluster':
            x_position = 0 + jitter_by
        else:
            x_position = self.species_x_position + jitter_by

        return x_position


    def add_all_nodes(self, node_type):
        '''
        this function add all nodes from a specific node type (species or cluster)
        each node gets a size, which is the logp MW test for the association between it (cluster/
        species) and the disease, a color which reflects whether the association is positive (purple) or
        negative (light blue) and a position (x_value is 0 for clusters and 50 for species)
        :param G: nx.Graph object
        :param node_type: string: 'cluster' or 'species'
        returns nx.Graph object with nodes
        '''
        if node_type == 'cluster': y_position_factor = self.y_position_factor_cluster
        else: y_position_factor = self.y_position_factor_species
        for i, node in enumerate(self.cluster_species_df[node_type].unique()):
            x_position = self.set_node_x_position(node_type)
            # print ('node: ',node)
            size = self.set_node_size(node, node_type)
            # print ('size: ', size)
            # color = None
            color = self.set_node_color(node, node_type)
            # print ('color: ', color)
            self.G.add_node(node,pos=(x_position, i * y_position_factor),size=size, color=color)

        return self.G

    def add_all_edges(self):
        '''
        this function add all edges, which are links between cluster and species.
        each edge has an attribute 'weight' which is equal to the logp MW test for the association
        between the cluster and the species, and will be used as the width of the edge line
        :param G: nx.Graph object
        :param cluster_species_df: dataframe
        returns nx.Graph object with edges
        '''
        cluster_species_wide = pd.pivot_table(data=self.cluster_species_df, values='-log10p_MW_cluster_species',
                                              index='cluster', columns='species').fillna(9999)
        for i in cluster_species_wide.index:
            for j in cluster_species_wide.columns:
                # print cluster_species_wide.loc[i, j]
                if cluster_species_wide.loc[i, j] != 9999:
                    self.G.add_edge(i, j, weight=cluster_species_wide.loc[i, j])
        return self.G

    def get_node_positions(self):
        '''
        generate a list of nodes position. list order is as the nodes order generated by G.nodes()
        '''
        return nx.get_node_attributes(self.G,'pos')

    def get_node_sizes(self):
        '''
        generate a list of nodes sizes. list order is as the nodes order generated by G.nodes()
        factor: string/float. the factor by which to multiply the logp value in order to get real size
        '''
        nodes = self.G.nodes()
        size_dict = nx.get_node_attributes(self.G, 'size')
        sizes = [size_dict[node] * self.node_size_factor for node in nodes]
        return sizes

    def get_node_colors(self):
        '''
        generate a list of nodes colors. list order is as the nodes order generated by G.nodes()
        factor: string/float. the factor by which to multiply the logp value in order to get real size
        '''
        nodes = self.G.nodes()
        color_dict = nx.get_node_attributes(self.G, 'color')
        colors = [color_dict[node] for node in nodes]
        return colors

    def get_edge_weights(self):
        '''
        generate a list of edge weights (widths). list order is as the edges order generated by G.edges()
        factor: string/float. the factor by which to multiply the logp value in order to get real size
        '''
        edges = self.G.edges()
        weights = [self.G[u][v]['weight'] * self.edge_width_factor for u,v in edges]
        return weights

    def gen_graph(self):

        #get all edges and nodes to G:
        print ('getting all nodes and edges to G...')
        self.G = self.add_all_nodes('cluster')
        self.G = self.add_all_nodes('species')
        self.G = self.add_all_edges()

        positions = self.get_node_positions()
        weights = self.get_edge_weights()
        sizes = self.get_node_sizes()
        colors = self.get_node_colors()

        print ('drawing network...')
        self.ax = nx.draw_networkx(self.G, pos = positions, with_labels=True,
                         width=weights, node_color=colors, alpha = 0.5,
                         node_size=sizes, font_size=self.label_font_size,
                         font_weight='bold', ax=self.ax)
        #pos=positions
        # self.ax.set_xticks([])
        # self.ax.set_yticks([])


        print ('saving files...')
        if not isdir(self.figure_dir):
            makedirs(self.figure_dir)
        num_nodes = len(self.G.nodes())
        self.fig.savefig(self.figure_dir + 'network_plot_%s' % num_nodes,
                         orientation='portrait')
        # path = self.figure_dir + 'trial.graphml'
        # print (path)
        # nx.write_graphml(self.G, path)


if __name__ == '__main__':
    fig,ax = plt.subplots(figsize = (12,30))
    # cluster_species_df_file = 'excel_files_new/cluster_species_for_fig.xlsx'
    cluster_species_df_file = 'any_outcome_min_shared_species01_max_shared_species1\
_min_shared_cluster03_max_shared_cluster1/cluster_species_phen_summary_df__fdr005top50rels.xlsx'
    fig5_graph = ClusterSpeciesGraph(fig=fig, ax=ax,
                                     cluster_species_df_file=cluster_species_df_file,
                                     y_position_factor_species = 10000, y_position_factor_cluster = 10000,
                                     node_size_factor = 750,
                                     species_x_position=50,
                                     label_font_size = 5, jitter=True,edge_width_factor = 0.3)

    fig5_graph.gen_graph()


