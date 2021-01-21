import numpy as np
import pandas as pd
from tqdm import tqdm
# from ..binarize import to_binary
# from .eca_stgs import step
# rng = np.random.RandomState(np.random.PCG64(np.random.SeedSequence(666)))
rng = np.random.default_rng()


def step(state, rule_bin):
    """ Takes a complete state vector of the system at a time point
    and provides the new value for Center """
    # calculate the encoding vector based on the rule binary
    # the sum of the element-wise product of this vector and the inputs
    # is the encoded state
    input_size = int(np.log2(rule_bin.shape[0]))
    expos = np.arange(input_size, 0, -1) - 1
    enc = 2**expos

    # we need to get new 3-vectors (Left, Center, Right) for each position
    # can use roll to "shift" the whole array left and right on a circle
    # and vstack to set them up as a list of 3-vectors
    vec = np.vstack((np.roll(state, 1), state, np.roll(state, -1))).astype(np.int8)
    encoded = vec.T.dot(enc).astype(np.int8)

    # lookup the index corresponding to the transition
    return rule_bin[sum(enc) - encoded]


def to_binary(n, digits=8):
    binary_digits = []
    for _ in range(digits):
        binary_digits.append(int(n % 2))
        n = int(n / 2)
    return np.array(binary_digits[::-1])


def array_in(a1, arr_list):
    return any((a1 == ai).all() for ai in arr_list)

# per run parameters
trials = 100
N = 102
max_step = 1000000
center_size = 20
rule = 15

# initialize in center only as in Langton 199x
# center = np.round(rng.uniform(size=20)).astype(np.int8)
# pad_size = int((N - center_size) / 2)
# pad = np.zeros(pad_size).astype(np.int8)
# state = np.hstack((pad, center, pad))

# initialize state and vectors for calculating stats
state = np.round(rng.uniform(size=N)).astype(np.int8)
periods = np.zeros(trials) - 1
transients = np.zeros(trials)

binary_rule = to_binary(rule, 8)
# if we want to print the thing

for ti in tqdm(range(trials)):
    state_history = []
    state = np.round(rng.uniform(size=N)).astype(np.int8)
    for si in range(max_step):
        # update the state
        state = step(state, binary_rule)
        # if we've seen this state before we're done
        if array_in(state, state_history):
            print(si)
            # we can calculate period and transient for this run
            for period_len, prev_state in enumerate(state_history[::-1]):
                if np.array_equal(prev_state, state):
                    periods[ti] = period_len + 1
                    transients[ti] = len(state_history) - period_len - 1
            # and move to next trial
            break
        # otherwise record state and continue
        else:
            state_history.append(state)

