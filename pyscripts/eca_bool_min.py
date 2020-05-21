import pandas as pd
from tqdm import tqdm
from cana.boolean_node import BooleanNode
from binarize import to_binary
from boolean_minimization_stats import find_permutability, find_wildcards, powerset

# get all ECA rules
n_rules = 2**(2**3)
inputs = ['0', '1', '2']
n_inputs = len(inputs)
stat_list = []
cana_list = []

for rule in tqdm(range(n_rules)):
    node_dict = {}
    output_arr = to_binary(rule)
    new_node = BooleanNode.from_output_list(output_arr, name=str(rule),
                                            inputs=inputs)

    # calculate my measures using cana
    wildcards_normed = find_wildcards(new_node, inputs=3, norm=True)
    permutability_normed = find_permutability(new_node, inputs=3, norm=True)

    # add wildcard entries to row
    node_dict['rule'] = rule
    for i, input in enumerate(inputs):
        node_dict['r*(' + input + ')'] = wildcards_normed[i]

    # iterate over all possible permutabile sets by looping over
    # powerset of inputs where cardinality > 1
    # add the normed permutabilities to the row
    for ps in list(powerset(n_inputs))[n_inputs:]:
        node_dict['p*' + str(ps)] = permutability_normed[ps]

    stat_list.append(node_dict)

df = pd.DataFrame(stat_list)
df.to_csv(snakemake.output[0])
