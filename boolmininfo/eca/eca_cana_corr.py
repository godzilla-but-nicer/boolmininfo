import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from ..plotting_functions.pid_plotter import pretty_labels_map

# load all of the neccesary dataframes, only need top half
wedge = pd.read_csv(snakemake.config['eca_decompositions']['wedge']).drop('Unnamed: 0', axis=1)
ccs = pd.read_csv(snakemake.config['eca_decompositions']['ccs']).drop('Unnamed: 0', axis=1)
imin = pd.read_csv(snakemake.config['eca_decompositions']['imin']).drop('Unnamed: 0', axis=1)
pm = pd.read_csv(snakemake.config['eca_decompositions']['pm']).drop('Unnamed: 0', axis=1)
bld = pd.read_csv(snakemake.config['eca_decompositions']['boolean']).drop('Unnamed: 0', axis=1)
# syndisc and gh not working yet. Also add them to lists below
# syndisc = pd.read_csv(snakemake.config['eca_decompositions']['syndisc'])
# gh = pd.read_csv(snakemake.config['eca_decompositions']['gh'])

# get the fancy labels for the PI atoms
pid_labels = imin.drop('rule', axis=1).columns
label_map = pretty_labels_map(pid_labels)

# load the canalization dataframe
canal = pd.read_csv(snakemake.input.cana)[:128]

# we need to fix the order of the columns
col_order = ['((0,), (1,), (2,))', '((0,), (1,))', '((0,), (2,))',
             '((1,), (2,))', '((0,), (1, 2))', '((1,), (0, 2))',
             '((2,), (0, 1))', '((0,),)', '((1,),)', '((2,),)', 
             '((0, 1), (0, 2), (1, 2))', '((0, 1), (0, 2))', 
             '((0, 1), (1, 2))', '((0, 2), (1, 2))', '((0, 1),)', '((0, 2),)', 
             '((1, 2),)', '((0, 1, 2),)', 'kr*', 'ks*', 'r(0)', 'r(1)',
             'r(2)', 's(0)', 's(1)', 's(2)']

# for each method we want to run correlations for each information atom against
# each canalization measure.
methods = [wedge, ccs, imin, pm, bld]
plain_labels = ['wedge', 'ccs', 'imin', 'pm', 'boolean']
fancy_labels = [r'$I_\wedge$', r'$I_{CCS}$', r'$I_{min}$', r'$I_\pm$', r'$BLD$']
cana_map = {'kr*': r'$k_r^*$', 'ks*': r'$k_s^*$',
            's(0)': r'$s^*(0)$', 's(1)': r'$s^*(1)$', 's(2)': r'$s^*(2)$',
            'r(0)': r'$r^*(0)$', 'r(1)': r'$r^*(1)$', 'r(2)': r'$r^*(2)$'}

for m, method in enumerate(methods):
    # combine PID and canalization dfs
    merged = method.merge(canal, on='rule').drop('rule', axis=1)

    # add a tiny bit of noise to fix correlations
    noisy = merged + np.random.uniform(low=-1e-6, high=1e-6, size=merged.shape)
    cor_mat = noisy[col_order].corr()
    pid_labels = [label_map[l] for l in col_order[:18]]
    cana_labels = [cana_map[l] for l in col_order[18:]]
    labels = pid_labels + cana_labels
    
    plt.figure()
    sns.heatmap(cor_mat, center=0, yticklabels=labels, 
                xticklabels=labels, cbar_kws={'label': 'Correlation'})
    plt.xlim((18, 28))
    plt.ylim((0, 18))
    plt.xlabel('Canalization Measure')
    plt.ylabel('PI atom')
    plt.title(fancy_labels[m])
    plt.tight_layout()
    plt.savefig(snakemake.output[plain_labels[m]])

