import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import time
from scipy.stats import entropy
# from ..binarize import to_binary
# from .eca_stgs import step
# rng = np.random.RandomState(np.random.PCG64(np.random.SeedSequence(666)))
print('Initializing timer, rng, functions, df, params...')
start_time = time.time()
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


def word_entropy(state_vector, word_size):
    # stack the array on itself to get the word vectors
    word_vecs = state_vector.copy()
    for si in range(1, word_size):
        word_vecs = np.vstack((word_vecs, np.roll(state_vector, si)))
    
    # quick way to encode the words as numbers
    encoding = 2**(np.arange(word_size, 0, -1) - 1)
    words = encoding.dot(word_vecs)
    _, word_counts = np.unique(words, return_counts=True)

    return entropy(word_counts)


# load the dataframe to get the list of rules to simulate
eca_df = pd.read_csv('data/eca_equiv_classes.csv', index_col=None)
sim_rules = eca_df[eca_df['rule'] == 3]['rule']

# per run parameters
trials = 1
N = 500
max_step = 100000
word_size = 3
center_size = 20

# initialize in center only as in Langton 199x
# center = np.round(rng.uniform(size=20)).astype(np.int8)
# pad_size = int((N - center_size) / 2)
# pad = np.zeros(pad_size).astype(np.int8)
# state = np.hstack((pad, center, pad))

# allow for command-line specification of rules
if len(sys.argv) > 1:
    sim_rules = [int(sys.argv[1])]

print('Starting Calculations...')

for rule in sim_rules:
    print('Rule:', rule)
    # initialize vectors for calculating stats
    periods = np.zeros(trials, dtype=np.int) - 1
    transients = np.zeros(trials, dtype=np.int) + max_step + 1

    binary_rule = to_binary(rule, 8)
    # if we want to print the thing

    for ti in range(trials):
        print('\t' + str(ti + 1) + '/' + str(trials))
        # rule specific features to track
        # state_history = np.zeros((max_step, N))
        entropies = np.zeros(max_step)
        cycles = []
        # initialize state randomly
        state = np.round(rng.uniform(size=N)).astype(np.int8)
        for si in range(max_step):
            # update the state
            # state_history[si] = state
            entropies[si] = word_entropy(state, word_size)
            state = step(state, binary_rule)

        # we're going to find the cycles by building a list of word entropies
        # before it repeats
        round_digits = 20
        print('\t\t Calculating cycle...')
        entropies = np.round(entropies, round_digits)
        end_ent = np.round(word_entropy(state, word_size), round_digits)
        cycle_len = np.argmax(entropies[::-1] == end_ent) + 1 # returns idx of frst True
        cycle = entropies[-cycle_len:]

        # now we need to find the first occurance of a value in that list
        # this line returns the first false when looking backward
        transient = np.argmin(np.isin(entropies, cycle)[::-1])

        # add these values to my list if we found an attractor
        # I dont have a good way to check this but I'll say if we see less than
        # 1.5 cycles. his is probably a bad criterion
        if (max_step - transient) + 1.5 * cycle_len > max_step:
            cycle_len = -1
            transient = max_step + 1 # ensures that recorded value is '-1'
            cycle = []
        periods[ti] = cycle_len
        transients[ti] = max_step - transient
        cycles.append(cycle)



    np.savetxt('data/big_eca_runs/' + str(rule) + '_periods.csv', 
               periods, fmt='%d')
    np.savetxt('data/big_eca_runs/' + str(rule) + '_transients.csv', 
               transients, fmt='%d')
    with open('data/big_eca_runs/' + str(rule) + '_cycles.csv', 'a') as fout:
        for cyc in cycles:
            fout.write(str(cyc))
            fout.write('\n')

        # This is all incase i decide I do need to actually look at 
        # specific states and can't rely on entropy as a summary
        # now we need to try to find the cycle
        # end_state = state_history[-1]
        # in_cycle = [end_state]
        # for cyc_i in range(state_history.shape[0] - 1):
        #     if np.sum(end_state == state_history[-(cyc_i+1)]) == N:
        #         break
        #     else:
        #         in_cycle.append(state_history[-(cyc_i+1)])

# smoothing filter
print('Completed in {:.3f} minutes'.format((time.time() - start_time) / 60))

plt.plot(range(max_step), entropies, c='grey')
# plt.plot(range(max_step), np.isin(entropies, cycle), c='k', linestyle='--')
plt.axvline(max_step - transient)
plt.savefig('plots/test_stat.png')
