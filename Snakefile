from os.path import join as j

configfile: "workflow/config.yaml"

DATA_DIR = config["data_dir"]

PAPER_DIR = config["paper_dir"]
PAPER_SRC, SUPP_SRC = [j(PAPER_DIR, f) for f in ("main.tex", "supp.tex")]
PAPER, SUPP = [j(PAPER_DIR, f) for f in ("main.pdf", "supp.pdf")]

# rule all:
#     input:
#         PAPER, SUPP

wildcard_constraints:
   seqfile = '\w+'

rule paper:
    input:
        PAPER_SRC, SUPP_SRC
    params:
        paper_dir = PAPER_DIR
    output:
        PAPER, SUPP
    shell:
        "cd {params.paper_dir}; make"

# PID over all three-input boolean functions
rule eca_ccs:
    input:
        "scripts/eca/eca_analysis_ccs.py"
    output:
        ccs="data/eca_decompositions/ccs_df.csv"
    script:
        "scripts/eca/eca_analysis_ccs.py"

rule eca_imin:
    input:
        "scripts/eca/eca_analysis_imin.py"
    output:
        imin="data/eca_decompositions/imin_df.csv"
    script:
        "scripts/eca/eca_analysis_imin.py"

rule eca_pm:
    input:
        "scripts/eca/eca_analysis_pm.py"
    output:
        pm="data/eca_decompositions/pm_df.csv"
    script:
        "scripts/eca/eca_analysis_pm.py"

rule eca_wedge:
    input:
        "scripts/eca/eca_analysis_wedge.py"
    output:
        wedge="data/eca_decompositions/wedge_df.csv"
    script:
        "scripts/eca/eca_analysis_wedge.py"

# Compare the three input PIDs
rule eca_pid_correlations:
    input:
        expand("data/eca_decompositions/{method}_df.csv", method=config['eca_decompositions'])
    output:
        red="plots/eca_pid_corr/redundancy.pdf",
        uni="plots/eca_pid_corr/unique.pdf",
        syn="plots/eca_pid_corr/synergy.pdf",
        avg="plots/eca_pid_corr/average.pdf",
        dist="plots/eca_pid_corr/corr_dists.pdf"
    script:
        "scripts/eca/plot_eca_pid_correlations.py"

# compute the canalization measures for all of the eca
rule eca_canalization:
    input:
        "scripts/eca/eca_canalization.py"
    output:
        "data/eca_decompositions/canalization_df.csv"
    script:
        "scripts/eca/eca_canalization.py"

# compare all of the PID results to the canalization measures
rule cana_pid_correlations:
    input:
        expand("data/eca_decompositions/{method}_df.csv", method=config['eca_decompositions']),
        cana="data/eca_decompositions/canalization_df.csv"
    output:
        wedge="plots/eca_cana_compare/wedge_cana_corr.pdf",
        ccs="plots/eca_cana_compare/ccs_cana_corr.pdf",
        imin="plots/eca_cana_compare/imin_cana_corr.pdf",
        pm="plots/eca_cana_compare/pm_cana_corr.pdf",
        boolean="plots/eca_cana_compare/boolean_cana_corr.pdf"
    script:
        "scripts/eca/eca_cana_corr.py"


# Explore four-input boolean partial information. We have to sample from
# the ~65k possible functions
rule sample_k4_rules:
    input:
        "scripts/k4/sample_k4_rules.py"
    output:
        "data/k4_samples.txt"
    script:
        "scripts/k4/sample_k4_rules.py"

# this performs the decompositions, at least BROJA doesn't work so we'll need
# to explore the other methods
rule k4_decomposition:
    input:
        rules="data/k4_samples.txt"
    output:
        "data/k4_decompositions/{method}_df.csv"
    script:
        "scripts/k4/k4_decomposition_{wildcards.method}.py"

# PID over two-input boolean functions to establish the importance of wildcards
rule k2_imin:
    input:
        "scripts/two_input/two_input_imin.py"
    output:
        "data/k2_decompositions/k2_imin_df.csv"
    script:
        "scripts/two_input/two_input_imin.py"

rule k2_wedge:
    input:
        "scripts/two_input/two_input_wedge.py"
    output:
        "data/k2_decompositions/k2_wedge_df.csv"
    script:
        "scripts/two_input/two_input_wedge.py"

rule k2_pm:
    input:
        "scripts/two_input/two_input_pm.py"
    output:
        "data/k2_decompositions/k2_pm_df.csv"
    script:
        "scripts/two_input/two_input_pm.py"

rule k2_ccs:
    input:
        "scripts/two_input/two_input_ccs.py"
    output:
        "data/k2_decompositions/k2_ccs_df.csv"
    script:
        "scripts/two_input/two_input_ccs.py"

rule k2_broja:
    input:
        "scripts/two_input/two_input_broja.py"
    output:
        "data/k2_decompositions/k2_broja_df.csv"
    script:
        "scripts/two_input/two_input_broja.py"

rule k2_dep:
    input:
        "scripts/two_input/two_input_dep.py"
    output:
        "data/k2_decompositions/k2_dep_df.csv"
    script:
        "scripts/two_input/two_input_dep.py"

# look at the two-input decompositions
rule k2_rule_tables:
    input:
        expand("data/k2_decompositions/k2_{method}_df.csv", method=config['k2_decompositions'])
    output:
        red="plots/k2_rules_pid/redundancy_by_rule.pdf",
        unq="plots/k2_rules_pid/unique_by_rule.pdf",
        syn="plots/k2_rules_pid/synergy_by_rule.pdf"
    script:
        "scripts/two_input/k2_rule_tables.py"

# lets also run all of the canalization measures
rule k2_canalization:
    input:
        "scripts/two_input/two_input_canalization.py"
    output:
        "data/k2_decompositions/canalization_df.csv"
    script:
        "scripts/two_input/two_input_canalization.py"

# make them into a table as with the PID
rule k2_canalization_table:
    input:
        "data/k2_decompositions/canalization_df.csv"
    output:
        "plots/k2_canalization/canalization_by_rule.pdf"
    script:
        "scripts/two_input/k2_rules_canalization.py"

# we are going to explore how canalization enables different kinds of
# information by looking at some specific look up tables and their schemata.
# we'll pull 2 rules for each redundancy, unique, and synergy
rule k2_exploration_luts:
    input:
        "scripts/two_input/k2_lut_tables.py"
    output:
        redlut="plots/k2_schemata/redundancy/lut_1.pdf",
        redwc="plots/k2_schemata/redundancy/sc_1.pdf",
        unqlut="plots/k2_schemata/unique/lut_3.pdf",
        unqwc="plots/k2_schemata/unique/sc_3.pdf",
        synlut="plots/k2_schemata/synergy/lut_6.pdf",
        synwc="plots/k2_schemata/synergy/sc_6.pdf",
        alllut="plots/k2_schemata/all/lut_4.pdf",
        allwc="plots/k2_schemata/all/sc_4.pdf"
    script:
        "scripts/two_input/k2_lut_tables.py"

# plot a bunch of decompositions with my little pid_plotter
rule eca_pid_plots:
    input:
        imin="data/eca_decompositions/imin_df.csv",
        ccs="data/eca_decompositions/ccs_df.csv",
        wedge="data/eca_decompositions/wedge_df.csv",
        bld="data/eca_decompositions/boolean_df.csv",
        pm="data/eca_decompositions/pm_df.csv"
    output:
        imin=directory('plots/eca_pid/imin/'),
        ccs=directory('plots/eca_pid/ccs/'),
        wedge=directory('plots/eca_pid/wedge/'),
        bld=directory('plots/eca_pid/bld/'),
        pm=directory('plots/eca_pid/pm/')
    script:
        "scripts/eca/plot_many_decompositions.py"

# Get state transition graphs for ECA up to 16 cells
rule eca_stgs:
    input:
        ancient('scripts/eca/eca_stgs.py')
    output:
        stgs=directory('data/eca_stgs/')
    script:
        'scripts/eca/eca_stgs.py'

# calculate and store a bunch of stats about the dynamics of the STGs
rule stg_dynamics:
    input:
        stg_dir='data/eca_stgs/'
    output:
        'data/eca_dynamics.csv'
    script:
        'scripts/eca/eca_dynamics.py'

# extract regressed features from the dynamics table
# these are functions of the number of cells for attractors, periods, etc.
rule eca_regress_dynamics:
    input:
        dyn_csv='data/eca_dynamics.csv'
    output:
        'data/eca_dynamics_regressed.csv'
    script:
        'scripts/eca/eca_regress_dynamics.py'