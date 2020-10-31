import numpy as np
from cana.boolean_node import BooleanNode
from itertools import chain, combinations
from itertools import permutations
from dit import Distribution
from dit.shannon.shannon import entropy


# modified from itertools documentation
def powerset(iterable):
    "powerset([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))


def PID_sets(k):
    """ This function returns a list of the sets in the redundancy lattice """
    double_powerset = list(chain.from_iterable(
        [combinations(powerset(range(k)), r) for r in range(1, k+1)]))
    keep_sets = []
    for subset in double_powerset:
        contains_subset = False
        for i in subset:
            for j in subset:
                if i != j and set(i).issubset(j):
                    contains_subset = True
                    break
            if contains_subset:
                break
        if not contains_subset:
            keep_sets.append(subset)
    return keep_sets


def exclude_subsets(iter_of_sets):
    """ This function takes an interable of sets and returns a new list with all
        sets that are a subset of any other eliminated """
    keep_sets = []
    for si in iter_of_sets:
        si = set(si)
        any_subset = False
        for sj in iter_of_sets:
            sj = set(sj)
            if si != sj and si.issubset(sj):
                any_subset = True
                break
        if not any_subset:
            keep_sets.append(si)

    return keep_sets


class LiteralDistribution:
    """ 
    This class holds all of the objects and methods required for the literal
    distribution algorithm I have been working on to compare PID and
    boolean minimization

    Attributes
    ----------
    """

    def __init__(self, gate: BooleanNode):
        self.gate = gate
        self.info_sets = None
        self._distributed = None
        self._ts_transitions = None
        self._ts_coverage = None

    def _expand_gi_inputs(self, sets, n_inputs=None):
        """ This function ensures that every input is a member of a set in the
            permutation groups but assigning individuals to sets if they are not
            in any other set. """
        # this function can take a set of sets and a number of inputs or use
        # the instance's number of inputs
        if not n_inputs:
            n_inputs = self.gate.k
        
        # we want the inputs as a set so we can do set operations with it
        input_set = set(range(n_inputs))
        # in groups are subsets of inputs that are permutable
        in_group = set([])
        # new_groups are sets of single inputs that are not contained in any
        # permutable group
        new_groups = []
        # sets is a set of sets that is passed in and normally would be the
        # sets of inputs that share the same permutaation symbol
        for s in sets:
            # this allows us to track the inputs that appear in one or more
            # permutable sets
            for i in s:
                in_group.add(i)
        # this gives us the set of inputs that is not in any permutable group
        singletons = input_set - in_group
        for s in singletons:
            new_groups.append({s})
        
        # give us back the permutable groups plus any singletons
        return sets + new_groups

    def _get_ts_transitions(self):
        """
        Returns a list of two-symbol schemata as well as their associated
        output values.
        """
        # CANA gives us back two lists of tuples describing the entries to the
        # two-symbol schemata LUT. The first list describes the transitions to 
        # OFF and the second list, the transitions to ON. These take the form:
        # (str(input states with 2 = #), list(list(permutable inputs)), 
        #  list(list(same state inputs)))
        self.gate.input_symmetry()
        ts0s, ts1s = self.gate._two_symbols

        _ts_transitions = []
        # iterate over the schemata and the output labels so we can easuly
        # preserve that information
        for output, ts in zip(['0', '1'], [ts0s, ts1s]):
            for inputs, permutables, same_symbols in ts:
                # all groups is the union of the permutable and same symbol
                # input sets
                all_groups = [set(s) for s in permutables + same_symbols]
                # we want sets of inputs that are not subsets of any other set
                # of inputs. if x={1, 2, 3} and y={1, 2} we only want x
                big_groups = exclude_subsets(all_groups)

                # add "permutable groups" for single inputs not belonging
                # to another group
                all_groups = self._expand_gi_inputs(big_groups, self.gate.k)

                # make a new entry for the list these will be similar to the
                # CANA output but with the input sets we care about
                new_entry = (output, inputs, all_groups)
                _ts_transitions.append(new_entry)

        return _ts_transitions

    def _expand_wildcards(self, state_list):
        """ recursively replaces wildcards in a set of input states """
        # state_list is a list of input states that covers all of the
        # meaningfully different permutations allowed by the groups
        # determined in the previous function
        # we use a set for expanded because all we care about is cardinality
        # and it enforces uniqueness on its members
        expanded = set([])
        for i, states in enumerate(state_list):
            states = tuple(states)

            # is the given input state string has no wildcards we dont need to
            # do anything to it. When we've recursuvely replaced all of the 2s
            # this is how we break out of the loop
            if '2' not in states:
                expanded.add(states)

            # otherwise we have to expand all of the twos into ones and zeros
            else:
                for j, char in enumerate(states):
                    # if we find a 2 make new copies of the input string with
                    # the 2 replaced with 1 and 0 and add them to the list of 
                    # expanded input strings then repeat
                    if char == '2':
                        gp_ones = list(states)
                        gp_ones[j] = '1'

                        gp_zeros = list(states)
                        gp_zeros[j] = '0'

                        expanded.add(tuple(gp_zeros))
                        expanded.add(tuple(gp_ones))

                expanded = self._expand_wildcards(expanded)
        return expanded

    def _calculate_ts_coverage(self):
        """ Take the transitions and calculate coverage values for each """
        # find the two-symbol list if not already done
        if not self._ts_transitions:
            self._ts_transitions = self._get_ts_transitions()

        # Expand the two-symbol list to get coverage values. We will be making
        # a new list of tuples that includes the number of original LUT entries
        # that each reduced transition represents.
        ts_coverage = []
        for transition, inputs, permutables in self._ts_transitions:
            input_array = np.array([char for char in inputs])
            # Every input will be a member of at least one permutable group for
            # every transition
            for subset in permutables:
                # convert the input set to a numpy array so we can use it for
                # fancy indexing
                subset = np.array(list(subset))
                gi = input_array[subset]

                # we want only permutations that are meaningfully different
                # we'll recursively replace wildcards with literals.
                gi_perm = np.unique(np.array(list(permutations(gi))), axis=0)
                permuted = list(gi_perm)
                expansions = self._expand_wildcards(permuted)

            # all that we really care about is how much of the LUT is covered
            # by the schemata
            lut_coverage = len(expansions)
            ts_coverage.append((transition, inputs, permutables, lut_coverage))

        self._ts_coverage = ts_coverage

    def _literal_distribution(self):
        """ Function that assigns need-to-know literals to subsets of inputs
            pretty sure its really poorly written """
        # We are going to assign the transitions we've been tracking to
        # literals that an observer must know to predict the transition with
        # certainty.

        if not self._ts_coverage:
            self._calculate_ts_coverage()

        # we're going to make a list of dictionaries. The inner dictionaries
        # will be keyed by sets of inputs
        distributed = []

        # the set of all inputs will be useful for performing set
        # operations later
        input_set = set(range(self.gate.k))

        # check each transition and assign the literals
        for tran, inputs, pos_free, coverage in self._ts_coverage:

            # set up the new dictionary
            transition_dist = {}

            # if there are no wildcards, than all inputs must be known to
            # determine transition and no further accounting is required for
            # this transition
            if '2' not in inputs:
                transition_dist[str(tuple(input_set))] = coverage

            # otherwise we have to do lots of counting
            else:
                # convert input string to array for numpy indexing
                arr_inp = np.array([c for c in inputs])
                twos = set(np.where(arr_inp == '2')[0].flatten())

                # we'll populate this literal groups list to find places where
                # literals can exist
                literal_groups = []
                literal_singletons = []

                # check each permutable group for twos
                for group in pos_free:
                    group_twos = twos.intersection(group)
                    n_literals = len(group - group_twos)
                    # if there are literals in the group we have information!
                    if n_literals > 0:
                        # any element of this list with more than one index
                        # identifies where literals COULD be. They will have
                        # to be split up because these are instances of
                        # redundancy
                        group_dist = list(combinations(group, n_literals))
                        for expand in group_dist:
                            if len(expand) > 1:
                                literal_groups.append(expand)
                            else:
                                literal_singletons.append(expand)
                print('literal groups before adding singletons', inputs)
                print(literal_groups)
                
                # once all of the literal positions have been taken we simply
                # have to add all of the singeltons to all of the groups
                if len(literal_groups) > 0:
                    for group in literal_groups:
                        for singleton in literal_singletons:
                            list_group = list(group)
                            list_group.append(singleton)
                            group = tuple(list_group)
                else:
                    literal_groups = literal_singletons
                
                print('after')
                print(literal_groups)

                # filter out any literal group that is a subset of any other
                # literal_groups = exclude_subsets(literal_groups)
                # literal_groups = [tuple(s) for s in literal_groups]

                # put these into the dictionary of input groups
                for group in literal_groups:
                    transition_dist[str(group)] = coverage
            
            # add the new dictionary to the list
            distributed.append(transition_dist)
                
        print('distributed:')
        for i, _ in enumerate(distributed):
            for g in distributed[i]:
                print(i, g, distributed[i][g])

        return distributed

    def _assign_information(self):
        # get the PID atoms
        keep_sets = PID_sets(self.gate.k)
        string_sets = []
        # the way we're accessing the variables later means we need these as 
        # strings. This way we can ensure that all of the strings we use as
        # keys match in our final output
        for ks in keep_sets:
            new_set_entry = []
            for tup in ks:
                new_set_entry.append(str(tup))
            string_sets.append(set(new_set_entry))

        print(string_sets)

        # put these in a dictionary so we can have zeros for missing values
        info_sets = {str(k): 0 for k in keep_sets}

        # calculate entropy of the output distribution for normalization
        output_dist = Distribution(
            self.gate.outputs, [1/2**self.gate.k]*2**self.gate.k)
        output_entropy = entropy(output_dist)

        # in some cases we can end up with more coverage than we should have
        # so we will normalize by the total coverage
        coverage_sum = 0

        # ok now we can gather our inputs into redundancies
        # make a label for the redundancy atom by finding input groups
        # that appear in the transition
        for t, trans in enumerate(self._distributed):
            # only have to do this if there are, in fact, literals somewhere
            if len(trans.keys()) > 0:
                redundant = set(trans.keys())
                print('Transition', t)
                print(redundant)
                # find the index of the set that matches this one
                for si, key in enumerate(string_sets):
                    if key == redundant:
                        info_i = si

                # every coverage value matches in a single transition so we can
                # use any of them
                coverage = max(trans.values())
                info_sets[str(keep_sets[info_i])] = coverage
                coverage_sum += coverage

            # if there are no literals we just have to avoid a division by zero
            # and move on
            else:
                coverage_sum = 1
        
        # normalize to the gate entropy and by coverage
        for k in info_sets:
            info_sets[k] = info_sets[k] / coverage_sum * output_entropy

        return info_sets

    def run_distribute_literals(self, gate=None):
        """
        The big function that performs the algorithm
        """
        self._calculate_ts_coverage()
        self._distributed = self._literal_distribution()
        self.info_sets = self._assign_information()

        return self.info_sets
