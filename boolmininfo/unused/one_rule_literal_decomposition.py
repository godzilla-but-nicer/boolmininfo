from literal_distribution_thing import LiteralDistribution, exclude_subsets, PID_sets
from cana.boolean_node import BooleanNode
from binarize import to_binary

# set up variables
n_inputs = 2**3
rule = 131
output = to_binary(rule, n_inputs)
bn = BooleanNode.from_output_list(outputs=output)
# print(bn.input_symmetry())
# print(bn._two_symbols)

ld = LiteralDistribution(bn)
print('Output of ts_transitions:')
print(ld._get_ts_transitions())
# print('\n')
# print('Output of Literal Distribution')
dist = ld._literal_distribution()
# for t in dist:
    # print('Transitions to', t)
    # for key in dist[t]:
        # print(key, dist[t][key])
# print('information assignment:')
ld._distributed = dist
info = ld._assign_information()
# for t in info:
    # print(t, info[t])
# print(len(info))