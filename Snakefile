from os.path import join as j

configfile: "workflow/config.yaml"

DATA_DIR = config["data_dir"]

PAPER_DIR = config["paper_dir"]
PAPER_SRC, SUPP_SRC = [j(PAPER_DIR, f) for f in ("main.tex", "supp.tex")]
PAPER, SUPP = [j(PAPER_DIR, f) for f in ("main.pdf", "supp.pdf")]

rule all:
    input:
        PAPER, SUPP

rule paper:
    input:
        PAPER_SRC, SUPP_SRC
    params:
        paper_dir = PAPER_DIR
    output:
        PAPER, SUPP
    shell:
        "cd {params.paper_dir}; make"


rule eca_decomposition:
    input:
        "scripts/eca_analysis_{method}.py"
    output:
        "data/eca_decompositions/{method}_df.csv"
    shell:
        "python scripts/eca_analysis_{wildcards.method}.py"

rule eca_pid_correlations:
    input:
        expand("data/eca_decompositions/{method}_df.csv", method=config['eca_decompositions'])
    output:
        red="plots/eca_pid_corr/redundancy.png",
        uni="plots/eca_pid_corr/unique.png",
        syn="plots/eca_pid_corr/synergy.png",
        avg="plots/eca_pid_corr/average.png",
        dist="plots/eca_pid_corr/corr_dists.png"
    script:
        "scripts/plot_eca_pid_correlations.py"
