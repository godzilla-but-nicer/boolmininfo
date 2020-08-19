import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv(snakemake.input[0])
df = df.set_index('rule')

label_map = {'kr*': r'$k_r^*$', 'ks*': r'$k_s^*$', 'r(0)': r'$r(0)$', 
             'r(1)': r'$r(1)$', 's(0)': r'$s(0)$', 's(1)': r'$s(1)$'}

def display_df(df):
    # only use half of the rules (symmetry)
    dfT = df[:8].T
    tab = sns.heatmap(dfT, xticklabels=dfT.columns, annot=True, cmap='Blues')
    # get some nice looking labels
    ylabels = []
    for ind in dfT.index:
        ylabels.append(label_map[ind])
    tab.set_yticklabels(ylabels, rotation=0)

    # label our axes too
    tab.set_xlabel('Rule')
    tab.set_ylabel('Measure')
    return tab

plt.figure()
mat = display_df(df)
plt.title('Canalization')
plt.savefig(snakemake.output[0])