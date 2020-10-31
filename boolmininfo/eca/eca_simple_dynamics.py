import numpy as np
import networkx as nx
import pandas as pd
from boolmininfo.binarize import to_binary, to_decimal
from copy import copy, deepcopy

df_list = []

rules = range(2**8)
n_cells = [12, 14, 16]

for rule in rules:
    for n_cell in n_cells:
        row = {}
        row['rule'] = rule
        row['n_cells'] = n_cell

        # load the graph
        STG = nx.read_graphml(snakemake.input.stg_dir + str(rule) +
                              '/' + str(n_cell) + '_cells.graphml')

        # number and lengths of attractors
        row['m_attractors'] = nx.number_attracting_components(STG)

        cycles = nx.simple_cycles(STG)
        cycle_lens = [len(c) for c in cycles]

        row['max_period'] = np.max(cycle_lens)

        # transient lengths
        t_len = []

        cycles = nx.simple_cycles(STG)
        # need to do this for each attractor
        for cycle in cycles:
            # set to ignore
            cyc_nodes = set(cycle)
            for node in cycle:
                # initialize the queue and our distance counter
                queue = [node]
                level = {node: 0}
                # keep going til the queue is done
                while queue:
                    # remove the first element, thats what we check next
                    base = queue.pop(0)
                    for pred in [p for p in STG.predecessors(base) if p not in cyc_nodes]:
                        queue.append(pred)
                        # we know that this must be one step farther out than the checked node
                        level[pred] = level[base] + 1
                t_len.extend([l for l in level.values()])

        row['mean_transient'] = np.max(t_len)

        # derrida coefficient
        derrida = 0

        # we check this value over all posible states
        for state in range(2**n_cell):
            # get next state from graph
            next_state = next(STG.neighbors(str(state)))
            # binary represenatation of initial state
            encoded_init = to_binary(state, n_cell)

            # for every possible perturbation
            for bit in range(n_cell):
                # flip a bit to make the perturbation
                perturb = deepcopy(encoded_init)
                if perturb[bit] == 0:
                    perturb[bit] = 1
                elif perturb[bit] == 1:
                    perturb[bit] = 0

                # get the next step in the perturbed state
                next_perturbed = next(STG.neighbors(str(to_decimal(perturb, n_cell))))

                bin_next_state = np.array(to_binary(int(next_state), n_cell))
                bin_next_pert = np.array(to_binary(int(next_perturbed), n_cell))

                derrida += np.sum(~(bin_next_state == bin_next_pert))

        derrida /= 2 * n_cell * 2**n_cell

        row['derrida_coeff'] = derrida

        # ratio between configs in transient vs configs in attractor
        total_configs = 2**n_cell
        total_attractor = np.sum(cycle_lens)
        row['transient_to_attractor'] = (total_configs - total_attractor) / total_attractor

        # add the row to the list
        df_list.append(row)

df = pd.DataFrame(df_list)

df.to_csv(snakemake.output[0])
