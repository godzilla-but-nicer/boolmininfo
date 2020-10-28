import dit
import pandas as pd
import numpy as np
from tqdm import tqdm
from ..binarize import make_dist_arr, make_input_strings

# load the sampled rules
rules = np.loadtxt(snakemake.input.rules)

# set up data structures for use later
n_inputs = 2**4
df_dict = []
input_ordering = make_input_strings(4)[::-1]

for rule in tqdm(rules):
    pis = {}  # gets filled with each pi term for this rule
    arr = make_dist_arr(rule, input_ordering, len(input_ordering))

    # use dit to calculate the PID
    dist = dit.Distribution(arr, [1/n_inputs]*n_inputs)
    pid = dit.pid.PID_WB(dist)

    # Update dictionary to contain row for this term
    pis['rule'] = rule
    for key in pid._pis.keys():
        pis[str(key)] = pid._pis[key]
    df_dict.append(pis)

df = pd.DataFrame(df_dict)
print(df.head())
df_fout = open(snakemake.output[0], 'w')
df.to_csv(df_fout)
df_fout.close()
