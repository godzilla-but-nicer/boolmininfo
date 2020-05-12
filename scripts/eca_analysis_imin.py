import dit
import pandas as pd
from tqdm import tqdm
from binarize import make_dist_arr

n_rules = 2**8
save_pid = {}
imins = {}
pis = {}
df_dict = []
pid_list = []
input_ordering = ['111', '110', '101', '100', '011', '010', '001', '000']

for rule in tqdm(range(n_rules)):
    pis = {}
    arr = make_dist_arr(rule, input_ordering)
    dist = dit.Distribution(arr, [1/8]*8)
    pid = dit.pid.PID_WB(dist)
    pid_list.append(pid)
    pis['rule'] = rule
    for key in pid._pis.keys():
        pis[str(key)] = pid._pis[key]
    df_dict.append(pis)

df = pd.DataFrame(df_dict)
print(df.head())
df_fout = open('imin_df.csv', 'w')
df.to_csv(df_fout)
df_fout.close()
