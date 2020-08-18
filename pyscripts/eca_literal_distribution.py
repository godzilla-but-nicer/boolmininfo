import pandas as pd
from tqdm import tqdm
from binarize import to_binary
from literal_distribution_thing import LiteralDistribution
from cana.boolean_node import BooleanNode

# set up variables
n_inputs = 2**3
n_rules = 2**8
df_dict = []

# figure out what isn't working
failing_rules = []


for rule in tqdm(range(n_rules)):
    pis = {}  # becomes row of dataframe
    output = to_binary(rule, 8)
    bn = BooleanNode.from_output_list(outputs=output)
    ld = LiteralDistribution(bn)
    pis = ld.run_distribute_literals()
    pis['rule'] = rule
    df_dict.append(pis)

# write out the dataframe
df = pd.DataFrame(df_dict)
df_fout = open('data/eca_decompositions/boolean.csv', 'w')
df.to_csv(df_fout)
df_fout.close()

with open('data/eca_decompositions/failing_rules.txt', 'w') as ff:
    text = '\n'.join(failing_rules)
    ff.write(text)