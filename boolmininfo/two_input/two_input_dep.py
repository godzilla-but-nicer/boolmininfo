import dit
import pandas as pd
from tqdm import tqdm
from ..binarize import make_dist_arr, make_input_strings

# set up variables
n_inputs = 2**2
n_rules = 2**(2**2)
df_dict = []
input_ordering = make_input_strings(2)[::-1]


for rule in tqdm(range(n_rules)):
    pis = {}  # becomes row of dataframe
    arr = make_dist_arr(rule, input_ordering, digits=4)

    # use dit to calculate decomposition
    dist = dit.Distribution(arr, [1/n_inputs]*n_inputs)
    pid = dit.pid.PID_dep(dist)

    # update the dictionary with the PI values
    pis['rule'] = rule
    for key in pid._pis.keys():
        pis[str(key)] = pid._pis[key]
    df_dict.append(pis)

# write out the dataframe
df = pd.DataFrame(df_dict)
df_fout = open(snakemake.output[0], 'w')
df.to_csv(df_fout)
df_fout.close()
