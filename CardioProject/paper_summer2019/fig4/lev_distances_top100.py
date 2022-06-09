import pandas as pd
import Levenshtein as lev


def gen_lev_longform(sequence_list):
    df = pd.DataFrame()
    count=0
    for i in range(len(sequence_list)):
        for j in range(i+1,len(sequence_list)):
            seq1 = sequence_list[i]; seq2 =  sequence_list[j]
            df.loc[count,'seq1'] = seq1
            df.loc[count, 'seq2'] = seq2
            df.loc[count,'lev_dist'] =lev.distance(seq1,seq2)
            count +=1
    df = df.sort_values(by='lev_dist')

    return df

top100_clusters_file_name = "top100_clusters_with_annot.xlsx"
dir1 = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Presentations and Manuscripts/CardioTCR paper/\
June2019/fig3_balanced_comparison/calc_fig3/excel_files_new/'
fig4_dir = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/Presentations and Manuscripts/CardioTCR paper/\
June2019/fig4_graph/'

top100_clusters = pd.read_excel(dir1+top100_clusters_file_name) ['cluster'].tolist()
print ('now calculating lev dist...')
cluster_lev = gen_lev_longform(top100_clusters)

cluster_lev.to_excel(fig4_dir+'cluster_lev_top100.xlsx')
