import os
import pandas as pd

eca = pd.read_csv('data/eca_equiv_classes.csv', index_col=None)
rules = eca['rule']

out_str = """#!/bin/bash
#####  Constructed by HPC everywhere #####
#SBATCH --mail-user=patgwall@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3-00:59:00
#SBATCH --partition=general
#SBATCH --mail-type=FAIL
#SBATCH --job-name=eca_{0:d}

######  Module commands #####



######  Job commands go below this line #####
cd ~/boolmininfo/
python boolmininfo/eca/stat_dynamics.py {1:d}"""

for rule in rules:
    submit_str = out_str.format(rule, rule)
    print(submit_str)