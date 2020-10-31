import numpy as np
from itertools import chain, combinations


# this function definetely has bad normalization right now (maybe counting
# wrong in general)
def find_wildcards(bool_node, inputs=3, norm=True, ones=False):
    # we can rely on cana to redescribe things for us
    pis = list(bool_node.pi_coverage().values())
    counts = np.zeros(inputs) # store number of wildcards

    if ones:
        # find all of the one transitions
        mask = bool_node.look_up_table()['Out:'].values.astype(int)
        mask = mask.astype(bool)  # there must be a better way to do this
        n_one_transitions = mask.sum()
        pi_ones = np.array(pis)[mask]

        # for each transition to one
        for pi in pi_ones:
            for sc in pi:
                for i in range(inputs):
                    if sc[i] == '2':
                        counts[i] += 1

        # normalize by number of one transitions
        if norm and n_one_transitions > 0:
            counts /= n_one_transitions

    else:  # consider full LUT
        # for each input set look at each input and find all of the wildcards
        for pi in pis:
            for sc in pi:
                for i in range(inputs):
                    if sc[i] == '2':
                        counts[i] += 1

        # normalize them by the number of input sets
        if norm:
            counts /= 2**inputs

    return counts


# This function is a modified recipe from itertools docs. accessed from
# https://stackoverflow.com/questions/18035595/powersets-in-python-using-itertools
def powerset(n_inputs):
    "powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    iterable = range(n_inputs)
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))


def find_permutability(bool_node, inputs=3, norm=True, ones=False):
    # extract redescriptions
    ts = list(bool_node.ts_coverage().values())

    # find all potential group invarient enputs (subsets of cardinality > 1)
    possible_ts = list(powerset(inputs))[inputs:]
    ts_dict = {sub: 0 for sub in possible_ts}

    for inp in ts:
        if len(inp) > 0:  # check that redescription exists
            # combine the two lists of permutable inputs, I dont care if its zero or 1
            ones = inp[0][1].copy()
            zeros = inp[0][2].copy()
            gi_enputs = ones + zeros
            # iterate over each and increase the group invariant enputs that exist
            if len(gi_enputs) >= 1:
                for gi in gi_enputs:
                    ts_dict[tuple(gi)] += 1

    if norm:
        for key in ts_dict.keys():
            ts_dict[key] = ts_dict[key] / 2**inputs

    return ts_dict
