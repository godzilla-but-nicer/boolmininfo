from cana.boolean_node import BooleanNode
from binarize import to_binary
from boolean_minimization_stats import find_permutability, find_wildcards, powerset

rule = 23
output = to_binary(rule)
inputs = 3
print(output)
node = BooleanNode.from_output_list(output, name=str(rule),
                                    inputs=['0', '1', '2'])

# extract redescriptions
ts = list(node.ts_coverage().values())

# find all potential group invarient enputs (subsets of cardinality > 1)
possible_ts = list(powerset(inputs))[inputs:]
ts_dict = {sub: 0 for sub in possible_ts}

for inp in ts:
    print(inp)

for inp in ts:
    if len(inp) > 0:  # the entries that cannot be redescribed are len zero
        # combine the two lists of permutable inputs, I dont care if its zero or 1
        ones = inp[0][1].copy()
        zeros = inp[0][2].copy()
        gi_enputs = ones + zeros
        # iterate over each and increase the group invariant enputs that exist
        if len(gi_enputs) >= 1:
            for gi in gi_enputs:
                ts_dict[tuple(gi)] += 1

for key in ts_dict.keys():
    ts_dict[key] = ts_dict[key] / 2**inputs

print(ts_dict)