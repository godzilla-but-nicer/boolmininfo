from ..plotting_functions.pid_plotter import pid_plot
import pandas as pd
import matplotlib.pyplot as plt

# which rules do we want to look at?
rules = [254, 150, 146, 190, 240, 110, 30]

# load the datasets
imin = pd.read_csv(snakemake.input.imin, index_col=0)
bld = pd.read_csv(snakemake.input.bld, index_col=0)
ccs = pd.read_csv(snakemake.input.ccs, index_col=0)
wedge = pd.read_csv(snakemake.input.wedge, index_col=0)
pm = pd.read_csv(snakemake.input.pm, index_col = 0)

for rule in rules:
    # plot imin
    imin_vals = imin[imin['rule'] == rule]
    plt.figure()
    pid_plot(imin_vals)
    plt.suptitle(r'$I_{min}$:' + 'Rule {}'.format(rule))
    plt.savefig(snakemake.output.imin + 'rule_{}'.format(rule) + '.pdf')

    # plot iccs
    ccs_vals = ccs[ccs['rule'] == rule]
    plt.figure()
    pid_plot(ccs_vals)
    plt.suptitle(r'$I_{CCS}$:' + 'Rule {}'.format(rule))
    plt.savefig(snakemake.output.ccs + 'rule_{}'.format(rule) + '.pdf')

    # plot wedge
    wedge_vals = wedge[wedge['rule'] == rule]
    plt.figure()
    pid_plot(wedge_vals)
    plt.suptitle(r'$I_{\wedge}$:' + 'Rule {}'.format(rule))
    plt.savefig(snakemake.output.wedge + 'rule_{}'.format(rule) + '.pdf')

    # plot pm
    pm_vals = pm[pm['rule'] == rule]
    plt.figure()
    pid_plot(pm_vals)
    plt.suptitle(r'$I_{\pm}$:' + 'Rule {}'.format(rule))
    plt.savefig(snakemake.output.pm + 'rule_{}'.format(rule) + '.pdf')

    # plot my thing
    bld_vals = bld[bld['rule'] == rule]
    plt.figure()
    pid_plot(bld_vals)
    plt.suptitle(r'$BLD$:' + 'Rule {}'.format(rule))
    plt.savefig(snakemake.output.bld + 'rule_{}'.format(rule) + '.pdf')

