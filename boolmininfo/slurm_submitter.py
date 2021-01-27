import os
import pandas as pd

eca = pd.read_csv('data/eca_equiv_classes.csv', index_col=None)
rules = eca[eca['wclass'] == 4]['rule']
submit_file = '/N/u/patgwall/BigRed3/boolmininfo/slurmy.script'

out_str = """#!/bin/bash
#####  Constructed by HPC everywhere #####
#SBATCH --mail-user=patgwall@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-08:00:00
#SBATCH --partition=general
#SBATCH --mail-type=FAIL
#SBATCH --job-name=eca_{0:d}
#SBATCH --output=logs/eca_{1}

######  Module commands #####



######  Job commands go below this line #####
cd ~/boolmininfo/
python boolmininfo/eca/stat_dynamics.py {2:d}
rm {3}"""


for rule in rules:
    submit_str = out_str.format(rule, rule, rule, submit_file)
    with open(submit_file, 'w') as fout:
        fout.write(submit_str)
    os.system('sbatch ' + submit_file)
