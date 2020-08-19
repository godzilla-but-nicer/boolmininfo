import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# load the dataframes and put the in a list
# we also need labels
dfs = []
labels = []
for fin in snakemake.input:
    # we can pull the labels from the filenames
    label = fin.split('/')[-1].split('.')[0].split('_')[1]
    labels.append(label)

    df = pd.read_csv(fin)
    dfs.append(df)

# what we really want is a set a new dictionaries for unique, redundancy, 
# and synergy. we can put them in a dictionary. we will initialize 
unq = '((0,),)'
red = '((0,), (1,))'
syn = '((0, 1),)'
term_list = [unq, red, syn]

# this dictionary is keyed by terms and has that kind of information for each 
# method and rule
info_rules = {}
# for each term pull that column from each dataframe
for term in term_list:
    term_by_method = {}
    term_by_method['rule'] = dfs[0]['rule']
    for dfi, df in enumerate(dfs):
        col = df[term]
        term_by_method[labels[dfi]] = col
    
    # update the information dataframe
    term_df = pd.DataFrame(term_by_method).set_index('rule')
    info_rules[term] = term_df

# fancy labels for making nice plots
label_map = {'imin': r'$I_{min}$', 'wedge': r'$I_{\wedge}$', 
              'broja': r'$I_{BROJA}$', 'pm': r'$I_{\pm}$', 'ccs': r'$I_{CCS}$',
              'dep': r'$I_{dep}$'}


def display_df(df):
    # only use half of the rules (symmetry)
    dfT = df[:8].T
    tab = sns.heatmap(dfT, xticklabels=dfT.columns, annot=True, cmap='RdBu')
    # get some nice looking labels
    ylabels = []
    for ind in dfT.index:
        ylabels.append(label_map[ind])
    tab.set_yticklabels(ylabels, rotation=0)

    # label our axes too
    tab.set_xlabel('Rule')
    tab.set_ylabel('Method')
    return tab

plt.figure()
mat = display_df(info_rules[unq])
plt.title('Unique')
plt.savefig(snakemake.output.unq)

plt.figure()
mat = display_df(info_rules[red])
plt.title('Redundancy')
plt.savefig(snakemake.output.red)

plt.figure()
mat = display_df(info_rules[syn])
plt.title('Synergy')
plt.savefig(snakemake.output.syn)

