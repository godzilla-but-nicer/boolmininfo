import pandas as pd
from tqdm import tqdm
from ..binarize import to_binary
from cana.boolean_node import BooleanNode

# set up variables
n_inputs = 2**2
n_rules = 2**(2**2)
df_dict = []


for rule in tqdm(range(n_rules)):
    canal = {}  # becomes row of dataframe
    arr = to_binary(rule, digits=4)
    print(arr)

    # use dit to calculate decomposition
    bn = BooleanNode.from_output_list(outputs=arr, name=rule)
    ks = bn.input_symmetry()
    kr = bn.input_redundancy()
    sym0, sym1 = bn.input_symmetry(mode='input')
    red0, red1 = bn.input_redundancy(mode='input')


    # update the dictionary with the PI values
    canal['rule'] = rule
    canal['kr*'] = kr
    canal['ks*'] = ks
    canal['r(0)'] = red0
    canal['r(1)'] = red1
    canal['s(0)'] = sym0
    canal['s(1)'] = sym1

    df_dict.append(canal)

# write out the dataframe
df = pd.DataFrame(df_dict)
df_fout = open(snakemake.output[0], 'w')
df.to_csv(df_fout, index=False)
df_fout.close()