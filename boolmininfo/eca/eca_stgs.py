import numpy as np
import networkx as nx
import os
from tqdm import tqdm
from ..binarize import to_binary

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



min_cells = 3
max_cells = 16

for rule in tqdm(np.arange(2**8)):
    for cells in np.arange(min_cells, max_cells + 1):
        STG = nx.DiGraph()
        # get a name for our output file make a dir for them
        graph_file = str(rule) + '/' + str(cells) + '_cells.graphml'
        if not os.path.isdir(snakemake.output.stgs + str(rule)):
            os.mkdir(snakemake.output.stgs + str(rule))

        state_labels = np.arange(2**cells)

        # get rule in binary form
        brule = np.array(to_binary(rule))

        # we want to encode our state vectors in big numbers
        expos = np.arange(cells, 0, -1) - 1
        enc = 2**expos

        for dec in state_labels:
            state = np.array(to_binary(dec, digits=cells))
            next_state = step(state, brule)
            next_dec = np.sum(enc * next_state).astype(np.int)
            STG.add_edge(dec, next_dec)
        


        nx.write_graphml(STG, snakemake.output.stgs + graph_file)
