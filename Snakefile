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
        "pyscripts/eca_analysis_ccs.py"
    output:
        ccs="data/eca_decompositions/ccs_df.csv"
    script:
        "pyscripts/eca_analysis_ccs.py"

rule eca_gh:
    input:
        "pyscripts/eca_analysis_gh.py"
    output:
        gh="data/eca_decompositions/gh_df.csv"
    script:
        "pyscripts/eca_analysis_gh.py"

rule eca_imin:
    input:
        "pyscripts/eca_analysis_imin.py"
    output:
        imin="data/eca_decompositions/imin_df.csv"
    script:
        "pyscripts/eca_analysis_imin.py"

rule eca_pm:
    input:
        "pyscripts/eca_analysis_pm.py"
    output:
        pm="data/eca_decompositions/pm_df.csv"
    script:
        "pyscripts/eca_analysis_pm.py"

rule eca_syndisc:
    input:
        "pyscripts/eca_analysis_syndisc.py"
    output:
        syndisc="data/eca_decompositions/syndisc_df.csv"
    script:
        "pyscripts/eca_analysis_syndisc.py"

rule eca_wedge:
    input:
        "pyscripts/eca_analysis_wedge.py"
    output:
        wedge="data/eca_decompositions/wedge_df.csv"
    script:
        "pyscripts/eca_analysis_wedge.py"
    
rule literal_distribution:
    input:
        "pyscripts/literal_distribution_thing.py"
    output:
        binary="data/eca_decompositions/binary_df.csv"
    script:
        "pyscripts/eca_literal_distribution.py"

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
        "pyscripts/plot_eca_pid_correlations.py"

# Collect statistics about compressibility of the LUT for the 3-input functions
rule eca_bool_min:
    input:
        "pyscripts/eca_bool_min.py"
    output:
        "data/eca_decompositions/bool_min.csv"
    script:
        "pyscripts/eca_bool_min.py"

# compare the three-input compressibilities and PIDs
rule plot_eca_imin_cana_correlations:
    input:
        "data/eca_decompositions/imin_df.csv"
    output:
        "plots/imin_cana_corr.pdf"
    script:
        "pyscripts/plot_imin_cana_corr.py"

# Explore four-input boolean partial information. We have to sample from
# the ~65k possible functions
rule sample_k4_rules:
    input:
        "pyscripts/sample_k4_rules.py"
    output:
        "data/k4_samples.txt"
    script:
        "pyscripts/sample_k4_rules.py"

# this performs the decompositions, at least BROJA doesn't work so we'll need
# to explore the other methods
rule k4_decomposition:
    input:
        rules="data/k4_samples.txt"
    output:
        "data/k4_decompositions/{method}_df.csv"
    script:
        "pyscripts/k4_decomposition_{wildcards.method}.py"