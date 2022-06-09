'''
This script calculate and plots different features across all samples
'''


#---------------------------------------------------------------------------------------------

from os import listdir
from os.path import isfile, join
# from Utils import Load, Write
import pandas as pd
import numpy as np
from scipy import stats
# import math
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
# from matplotlib import cm
# import plotly.plotly as py
from myplots import roundup, rounddown, find_decimal_fold, percentile_cut_off, rarefaction_calc, rarefaction_plot
from matplotlib.ticker import FormatStrFormatter
import cPickle as pickle
from Bio.SeqUtils import GC
from pop_calcs import *
from pop_organize import *





#___________________________________________________________________________________________




#########################
#    main function      #
#########################
figlist=[] ## generate figlist for pdf creation
create_pdf=True
n_samples=6

##generate list of sample file names:
onlyfiles = [f for f in listdir("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/CSVfiles") if isfile(join("/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/CSVfiles", f))]
onlyfiles = [datafile for datafile in onlyfiles if datafile.startswith ('HIP') and datafile.endswith('.csv')]
sample_names=[re.sub('.csv', '', datafile) for datafile in onlyfiles] ## generate list of 
                                                                      ## sample names out 
                                                                      ##of the sample file names
df_file_names,samples_with_df=get_sample_with_dfs() ## check which samples have df
generate_dfs=False ##False=load them from pickles, True=generate them. 

df_n=2 ## define number of dfs (prod, non-prod) *** check for update!***


## define function lists:
percProd_func_list=[perc_prod]
general_function_list=[unique_seq_n, unique_aa_n,norm_uniqe_nt_sequences, norm_uniqe_aa_sequences, max_nt_per_aa,mean_nt_per_aa, gc_content] ##*** check for update!***
clonality_func_list=[top_clonal_nt, top_1000clons_nt, Percentile1_clone_nt, median_clone_nt, Percentile999_clone_nt, mean_clonal_nt, std_clonal_nt, top_clonal_aa, top_1000clons_aa, Percentile1_clone_aa, median_clone_aa, Percentile999_clone_aa, mean_clonal_aa, std_clonal_aa]
diversity_func_list=[shannon_div_nt, simpson_div_nt, bpi_div_nt, shannon_div_aa, simpson_div_aa, bpi_div_aa]
length_func_list=[mean_cdr3, sv_cdr3, mean_vDeletion, sv_vDeletion, mean_n1Insertion, sv_n1Insertion, mean_d5Deletion, sv_d5Deletion, 
                  mean_d3Deletion, sv_d3Deletion, mean_n2Insertion, sv_n2Insertion, mean_jDeletion, sv_jDeletion]


## geneUsage_func_list=[]




##define dataframes for collecting features over all samples ** can use 'samples_with_df' 'sample_names':
## generate columns for each expected feature and fill them with 0's:
perProd_res_df=generate_res_DF(samples_with_df, percProd_func_list, 1)
general_res_df=generate_res_DF(samples_with_df, general_function_list, df_n)
clonality_res_df=generate_res_DF(samples_with_df, clonality_func_list, df_n)
diversity_res_df=generate_res_DF(samples_with_df, diversity_func_list, df_n)
length_res_df=generate_res_DF(samples_with_df, length_func_list, df_n)
####geneUsage_res_df=generate_res_DF(samples_with_df, geneUsage_func_list, df_n)
#### public_res_df=generate_res_DF(samples_with_df, public_func_list, df_n)




n=1
for d in range(len(samples_with_df[:n_samples])): ##***change here for more samples!***
    sample_name=samples_with_df[d]
    print n
    print sample_name
## extract prod and non-prod dfs: 
    sample_df, sample_df_prod, sample_df_non_prod=get_sample_data(sample_name, generate_dfs)
    percent_prod=perc_prod(sample_df)
##  perProd_res_df.loc[sample_name, 'perc_prod_df_0']=percent_prod
    general_res_df=calculate_res_df(general_res_df,general_function_list, sample_name, sample_df_prod, sample_df_non_prod) 
    clonality_res_df=calculate_res_df(clonality_res_df,clonality_func_list, sample_name, sample_df_prod, sample_df_non_prod)
    diversity_res_df=calculate_res_df(diversity_res_df,diversity_func_list, sample_name, sample_df_prod, sample_df_non_prod)
    length_res_df=calculate_res_df(length_res_df,length_func_list, sample_name, sample_df_prod, sample_df_non_prod)
    n=n+1



with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/general_res_df_%s_samples' %n_samples, "wb") as f1:
    pickle.dump(general_res_df, f1)
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/perProd_res_df_%s_samples' %n_samples, "wb") as f2:
    pickle.dump(perProd_res_df, f2)
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/clonality_res_df_%s_samples' %n_samples, "wb") as f3:
    pickle.dump(clonality_res_df, f3)
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/diversity_res_df_%s_samples' %n_samples, "wb") as f4:
    pickle.dump(diversity_res_df, f4)
with open('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Pickles/length_res_df_%s_samples' %n_samples, "wb") as f5:
    pickle.dump(length_res_df, f5)
    

## plot the result dataframes 
perc_prod_fig=plot_1_res_df(perProd_res_df, percProd_func_list, 'Percent Productive', 'linear')
general_res_fig_linear=plot_res_df(general_res_df, general_function_list, 'General (linear)', 'linear', 0.36)
general_res_fig_log=plot_res_df(general_res_df, general_function_list, 'General (log)', 'log', 0.36)
clonality_res_fig_linear=plot_res_df(clonality_res_df, clonality_func_list, 'Clonality (linear)', 'linear', 0.68)
clonality_res_fig_log=plot_res_df(clonality_res_df, clonality_func_list, 'Clonality (log)', 'log', 0.68)
diversity_res_fig_linear=plot_res_df(diversity_res_df, diversity_func_list, 'Diversity (linear)', 'linear', 0.36)
diversity_res_fig_log=plot_res_df(diversity_res_df, diversity_func_list, 'Diversity (log)', 'log', 0.36)
length_res_fig_linear=plot_res_df(length_res_df, length_func_list, 'length (linear)', 'linear', 0.66)
length_res_fig_log=plot_res_df(length_res_df, length_func_list, 'length (log)', 'log', 0.66)
diversity_res_fig_linear=plot_res_df(diversity_res_df, diversity_func_list, 'Diversity (linear)', 'linear', 0.36)
diversity_res_fig_log=plot_res_df(diversity_res_df, diversity_func_list, 'Diversity (log)', 'log', 0.36)

cdate=str(time.strftime("%d%m%Y"))
create_pdf=True
figlist=[perc_prod_fig, general_res_fig_linear, general_res_fig_log, clonality_res_fig_linear, clonality_res_fig_log, diversity_res_fig_linear,
        diversity_res_fig_log,length_res_fig_linear, length_res_fig_log]

if create_pdf:
    with PdfPages('/net/mraid08/export/genie/Lab/Personal/ShaniBAF/TCR_demo_data/Images/population_View_%s_samples_%s.pdf' %(n_samples, cdate)) as pdf:
        for fig in figlist:
            pdf.savefig(fig)
    pdf.close
else:
    plt.show()
