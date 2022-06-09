import cPickle as pickle
import pandas as pd
import numpy as np
import random
import matplotlib as mpl
import matplotlib.pyplot as plt


#definitions:
GENIE_BASE_DIR = '/net/mraid08/export/genie/Lab/Personal/ShaniBAF/'
DATA_DIR = '/net/mraid08/export/jafar/Microbiome/Analyses/TCR/'
SAMPLE_LIST_DIR = GENIE_BASE_DIR + 'Sample files/BD lists/'
AllUniqueWithCounts_path='/net/mraid08/export/jafar/Microbiome/Analyses/ShaniBAF/sharingAnalysis/\
AllUniqueWithCounts'
FIGURE2_DIR = GENIE_BASE_DIR + 'Presentations and Manuscripts/CardioTCR paper/June2019/calc_fig2/'

seed = None

# get sample lists
sample_list1_name = 'PNP530_Age4570_Gender_Male11_HbA1C4464_BMI040'
sample_list1_file = 'PNP530_Age4570_Gender_Male11_HbA1C4464_BMI040.pkl'
sample_list2_name = 'Cardio126_Age3566_Gender_Male11_HbA1C4464_BMI040'
sample_list2_file =  'Cardio126_Age3566_Gender_Male11_HbA1C4464_BMI040.pkl'

with open(SAMPLE_LIST_DIR + sample_list1_file,'rb') as fp1: sample_list1 = pickle.load(fp1)
with open(SAMPLE_LIST_DIR + sample_list2_file,'rb') as fp2: sample_list2 = pickle.load(fp2)


#get all unique:
AllUniqueWithCounts = AllUniqueWithCounts=pd.read_pickle(AllUniqueWithCounts_path)
AllUniqueWithCounts = AllUniqueWithCounts.drop(['frequencyCount (%)','prod_stat','isPublic'],axis=1)

AllUniqueWithCounts_filtered = AllUniqueWithCounts[AllUniqueWithCounts['Sample'].isin(sample_list1+sample_list2)]

unique_count = AllUniqueWithCounts_filtered[['Sample', 'nShared']].groupby('Sample').count()
unique_count['isCardio'] = np.where((unique_count.index.str.replace('BD','')).astype(int) > 949, 1, 0)


#count how many unique per samples and get stats
for name,group in unique_count.groupby('isCardio'):
    print (name)
    print (group.describe())
    print ('')

for seq_thresh in [5000,7500,9000]:
    print ( seq_thresh); print ( '')
    for name,group in unique_count.groupby('isCardio'):
        print ( name)
        print ( len(group[group['nShared']>=seq_thresh]))
    print ( '')







#get rid of samples with low number of seqs

def sample_uniques(unique_df,sample_list,n_seqs,n_samples,seed=None):

    if seed is not None: random.seed(seed)

    unique_df = unique_df[unique_df['Sample'].isin(sample_list)]
    seq_count = unique_df.groupby('Sample').count()
    enough_list = seq_count[seq_count['nShared'] >= n_seqs].index.tolist()

    random_samples = random.sample(enough_list,n_samples) #random selection of samples
    # unique_df_filtered = unique_df[unique_df['Sample'].isin(random_samples)]
    df_list=[]
    for sample in random_samples:
        df = unique_df[unique_df['Sample'] == sample]
        df_sampled = df.sample(n_seqs,random_state=seed)
        df_list.append(df_sampled)
    sampled_unique_df = pd.concat(df_list)

    return sampled_unique_df

pnp_sampled_5000_40 = sample_uniques(AllUniqueWithCounts_filtered,sample_list1,n_seqs=5000,n_samples=40)
cardio_sampled_5000_40 = sample_uniques(AllUniqueWithCounts_filtered,sample_list2,n_seqs=5000,n_samples=40)
pnp_sampled_9000_20 = sample_uniques(AllUniqueWithCounts_filtered,sample_list1,n_seqs=9000,n_samples=20)
cardio_sampled_9000_20 = sample_uniques(AllUniqueWithCounts_filtered,sample_list2,n_seqs=9000,n_samples=20)

# calc new sharing counts:

sum_df = pd.DataFrame()
df_list = [('pnp_sampled_5000_40',pnp_sampled_5000_40),('cardio_sampled_5000_40',cardio_sampled_5000_40),
           ('pnp_sampled_9000_20',pnp_sampled_9000_20),('cardio_sampled_9000_20',cardio_sampled_9000_20),
           ('all_sampled_5000_40',pd.concat([pnp_sampled_5000_40,cardio_sampled_5000_40]))]

for item in df_list:
    name=item[0]; df = item[1]

    df_grouped_by_seqs = df.reset_index()[['index','Sample']].groupby('index').count()
    df_grouped_by_seqs = df_grouped_by_seqs.rename(columns={'Sample':'n_sampled_shared'})
    n_seqs = df_grouped_by_seqs.index.nunique()
    mean_sharing_per_seq = df_grouped_by_seqs['n_sampled_shared'].mean()
    max_sharing_per_seq = df_grouped_by_seqs['n_sampled_shared'].max()

    sum_df.loc[name,'n_unique_seqs'] = n_seqs
    sum_df.loc[name, 'mean_sharing_per_seq'] = mean_sharing_per_seq
    sum_df.loc[name, 'max_sharing_per_seq'] = max_sharing_per_seq

    for n in [2, 4, 6, 8, 10]:
        shared_by_more=df.groupby('Sample').agg(lambda x: (float(len([i for i in x if i>=n]))) / len(x))
        mean_sharedbymore =  float(shared_by_more.mean())
        sum_df.loc[name, 'mean_sharedby%smore' %n] = mean_sharedbymore

    private = df.groupby('Sample').agg(lambda x: (float(len([i for i in x if i == 1]))) / len(x))
    mean_private = float(private.mean())
    sum_df.loc[name, 'mean_private'] = mean_private

    if name == 'all_sampled_5000_40':
        pnp_private = private[(private.index.str.replace('BD','').astype(int)) < 950]
        cardio_private = private[(private.index.str.replace('BD', '').astype(int)) >= 950]
        mean_pnp_private = float(pnp_private.mean())
        mean_cardio_private = float(cardio_private.mean())
        sum_df.loc[name, 'mean_pnp_private'] = mean_pnp_private
        sum_df.loc[name, 'mean_cardio_private'] = mean_cardio_private



sum_df.to_excel(FIGURE2_DIR + 'sharing_sum_df_seed%s.xlsx' %seed)

######################################
# conclusion: didnt see any significant different in sharing rates between pnp and cardio groups
#####################################

#calculate pairwise distances between samples