import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# load all of the neccesary dataframes
broja = pd.read_csv(snakemake.config['eca_decompositions']['broja'])
ccs = pd.read_csv(snakemake.config['eca_decompositions']['ccs'])
dep = pd.read_csv(snakemake.config['eca_decompositions']['dep'])
imin = pd.read_csv(snakemake.config['eca_decompositions']['imin'])
pm = pd.read_csv(snakemake.config['eca_decompositions']['pm'])

# Make new dataframes by PI term
df_list = [broja, ccs, dep, imin, pm]
labels = ['broja', 'ccs', 'dep', 'imin', 'pm']
terms = broja.drop(['Unnamed: 0', 'rule'], axis=1).columns
term_dfs = {}

# collect all five methods' answer for each term
for term in terms:
    term_dict = {}
    for m in range(len(df_list)):
        value = df_list[m][term]
        term_dict[labels[m]] = value

    new_df = pd.DataFrame(term_dict)
    term_dfs[term] = new_df


# Here's a function I can use to draw my correlation matrices
# I want to mark the labels of methods that have missing values for 
# each PI term
def correlation_matrix(df):
    # add asterisk to names of columns with nans
    labels = df.columns
    df_marked = df.copy()
    name_map = {}
    for lab in labels:
        if np.sum(pd.isnull(df[lab])) > 0:
            # can missing values be dropped or do we need to penalize
            # the correlation in some way?
            df_marked[lab] = df_marked[lab].dropna()
            name_map.update({lab: lab + '*'})

    # calculate the correlations
    df_marked = df.rename(columns=name_map)
    df_corr = df_marked.corr()

    # get the correlations with pandas
    mat = sns.heatmap(df_corr, annot=True, xticklabels=df_marked.columns,
                      yticklabels=df_marked.columns)

    return mat


# synergy plot
plt.figure()
mat = correlation_matrix(term_dfs['((0, 1, 2),)'])
plt.title('3-way Synergy')
plt.savefig(snakemake.output.syn)

# redundancy plot
plt.figure()
mat = correlation_matrix(term_dfs['((0,), (1,), (2,))'])
plt.title('3-way Redundancy')
plt.savefig(snakemake.output.red)

# unique plot
plt.figure()
mat = correlation_matrix(term_dfs['((1,),)'])
plt.title('Input 1 Unique Information')
plt.savefig(snakemake.output.uni)

# We also care about avg correlations between these across terms
corr_list = []
for df in term_dfs:
    corr = term_dfs[df].corr()
    corr_list.append(corr.values)

# I won't recompute missing values here instead I'll just hardcode the
# marked labels
new_labels = []
for lab in term_dfs['((1,),)'].columns:
    if lab == 'dep' or lab == 'broja':
        new_labels.append(lab + '*')
    else:
        new_labels.append(lab)

overall_correlations = np.mean(corr_list, axis=0)
plt.figure()
sns.heatmap(overall_correlations, annot=True, xticklabels=new_labels,
            yticklabels=new_labels)
plt.savefig(snakemake.output.avg)

# Now we will find the distributions of correlations for each pair
# extract upper (or lower) triangle of correlations
corr_arr = np.array(corr_list)
tri = np.array(np.triu_indices(corr_arr[0].shape[0], k=1))

# get tuples of labels and values
dist_labels = []
dist_vals = []
for ci in range(tri.shape[1]):
    label = '{};{}'.format(labels[tri[0, ci]], labels[tri[1, ci]])
    dist = corr_arr[:, tri[0, ci], tri[1, ci]]
    for d in dist:
        dist_labels.append(label)
        dist_vals.append(d)

# convert them to a dataframe for plotting
dist_dict = {'label': dist_labels, 'val': dist_vals}
dist_df = pd.DataFrame(dist_dict)

plt.figure(figsize=(7, 7))
sns.boxplot(x='label', y='val', data=dist_df)
plt.xticks(rotation=90)
plt.xlabel('PID pair')
plt.ylabel('Correlation by PI term')
plt.savefig(snakemake.output.dist)