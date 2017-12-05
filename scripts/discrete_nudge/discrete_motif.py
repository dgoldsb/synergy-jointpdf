"""
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This file builds a discrete GRN motif, and contains some of the common operations used.
"""

from __future__ import print_function

from copy import deepcopy
import itertools
import json
from operator import itemgetter
import os
import random
import sys

from jointpdf.jointpdf import JointProbabilityMatrix
import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.realpath(__file__)), '../..')
CONFIG = os.path.join(ROOT, 'config')

class DiscreteGrnMotif(JointProbabilityMatrix):
    """
    A motif which consists of:

    - a set of Genes, each assigned a number
    - a set of dependent PMFs that define probabilities that each gene
      is activated (FullNestedArrayOfProbabilities object)
    - a Boolean rule to determine states of the outputs between the genes

    The number of states a gene can be in is defined in self.numvalues.
    We can advance time with one timestep, and we can use measures such as the KL-divergence.
    When we advance time, each leaf of the original tree will sprout a new tree.
    If we advance time for all vars, we can sum all these trees together to get the object for t+1.
    It is an extension of the jointpdf package from Rick Quax.
    """

    def __init__(self, numvariables, numvalues, mode='uniform'):
        # call the super init
        super(DiscreteGrnMotif, self).__init__(numvariables, numvalues, mode)

        # set evolution style, network or transition_table
        self.evaluation_style = 'network'

        # global variables of a GRN network
        # TODO: refactor this away, the dict is annoying
        self.grn_vars = {}
        ## important to store the number of genes, as this is not the same as the numvariables
        self.grn_vars["gene_cnt"] = 0
        ## the rules are not part of the original framework
        self.grn_vars["rules"] = []
        ## correlation matrix for the dependent PDF of the initial conditions
        self.grn_vars["correlations"] = None
        ## set the default to downregulation dominates
        self.grn_vars["conflict_rule"] = 'totaleffect'
        ## also there has been a sequence of states, the latest is the current state
        self.states = []

        # global variables of a transition table
        self.transition_table = None

    def setter_grn_to_file(self, filename):
        """
        Set the settings for the GRN from a json.

        PARAMETERS
        ---
        filename: string
        """
        filepath = os.path.join(CONFIG, filename)
        with open(filepath, 'w') as outfile:
            json.dump(self.grn_vars, outfile)

    def getter_grn_from_file(self, filename):
        """
        Save the settings for the GRN to a json.

        PARAMETERS
        ---
        filename: string
        """
        filepath = os.path.join(CONFIG, filename)
        with open(filepath) as infile:
            self.grn_vars = json.load(infile)

    def construct_grn(self):
        """
        Creates a FullNestedArrayOfProbabilities object from the settings provided.
        """
        # reset the object
        self.reset(1, self.numvalues, self.joint_probabilities)

        # nasty thing with the states intially taking over those of the last created object
        self.states = []
        
        # check
        if self.transition_table is None:
            self.set_transition_table(self.grn_vars["conflict_rule"])

        # add all the genes
        if self.grn_vars["correlations"] is None:
            # use the default append_variable with random dependencies
            for _ in range(1, self.grn_vars["gene_cnt"]):
                self.append_variables(num_added_variables=1
                                      , added_joint_probabilities=None)
        else:
            # use the total correlation append_variable
            for i in range(1, self.grn_vars["gene_cnt"]):
                # find the correlation of the to-be-appended gene with previous
                row = self.grn_vars["correlations"][i]
                correlation = row[self.numvariables - 1]
                self.append_variable_correlation(correlation)
            # resample the probabilities, otherwise they get messed up somehow
            self.generate_random_joint_probabilities()

        # permute the transition table
        # there are n! possible  transition tables for the same motif
        # as we can always switch labels around
        # as such, we generate all permutations, sort them and pick the top one
        permutations = self.find_transition_table_permutations()
        permutations.sort()
        self.transition_table = permutations[0] 
        
        self.states.append(self.joint_probabilities.joint_probabilities)

    def append_variable_correlation(self, correlation):
        """
        Adds variables, enforcing set total correlations between them.

        I thought a bit more about this, and I think the correlation matrix
        is defined by the values under/above the diagonal.
        With these values given, you can derive the other correlations
        for a match/no match situation.
        For this reason, we only need these values to add a correlated variable.
        Even better, to add a new correlated variables we only need to know
        its correlation with the previously added variable.

        This should still be normalized, which is good!

        PARAMETERS
        ---
        correlations: defines the correlations with all existing variables
                      (list of floats)
        """
        # deepen the tree, use the correlation with the previous variable
        states = list(range(self.numvalues))
        leafcodes = list(itertools.product(states, repeat=self.numvariables))
        for _leafcode in leafcodes:
            leafcode_long = list(_leafcode) + [0]
            self.joint_probabilities.joint_probabilities =\
                self.deepen_leafcode_with_corr(leafcode_long,\
                self.joint_probabilities.joint_probabilities, correlation)

        # fix the object to a numpy array
        self.joint_probabilities.joint_probabilities =\
            np.array(self.joint_probabilities.joint_probabilities)

        # check normalization
        self.check_normalization()

        # increase the number of variables
        self.numvariables += 1

    def append_rule(self, inputs, outputs, rulefunction):
        """
        Append rule to the motif, keep in mind that the order of the inputs could matter.
        The function defines a truthtable, which is part of the motif operations.
        Rulefunctions are applied in a set order, so it does not matter in which order you add them.

        PARAMETERS
        ---
        inputs: list of gene indices (integers)
        outputs: list of gene indices (integers)
        rulefunction: function to call to evaluate this rule
        """

        # check if the gene actually exists
        for _input in inputs:
            if (_input < self.grn_vars["gene_cnt"]) is False:
                raise RuntimeError("you are referencing a non-existing gene")
        for _output in outputs:
            if (_output < self.grn_vars["gene_cnt"]) is False:
                raise RuntimeError("you are referencing a non-existing gene")

        # we do not validate if the number of inputs/outputs match the rule
        # add the rule
        rule = {}
        rule["inputs"] = inputs
        rule["outputs"] = outputs
        rule["rulefunction"] = rulefunction
        self.grn_vars["rules"].append(rule)

    def fetch_leafcode(self, leafcode, nested):
        """
        Fetches the value at the leafcode in the list.

        PARAMETERS
        ---
        leafcode: list of integers (corresponding to Boolean values)
        nested: a nested list, depth same as length of leafcode
        """
        if len(leafcode) > 1:
            return self.fetch_leafcode(leafcode[1:], nested[leafcode[0]])
        else:
            return nested[leafcode[0]]

    def add_leafcode(self, leafcode, nested, addition):
        """
        Alters the value at the leafcode in the list.

        PARAMETERS
        ---
        leafcode: list of integers (corresponding to Boolean values)
        nested: a nested list, depth same as length of leafcode
        addition: float
        """
        if len(leafcode) > 1:
            new_branches = []
            for i in range(0, self.numvalues):
                if i == leafcode[0]:
                    new_subbranch = self.add_leafcode(leafcode[1:], nested[leafcode[0]], addition)
                    new_branches.append(new_subbranch)
                else:
                    old_subbranch = nested[i]
                    new_branches.append(old_subbranch)
            return np.array(new_branches)
        else:
            new_branches = []
            for i in range(0, self.numvalues):
                if i == leafcode[0]:
                    new_subbranch = nested[leafcode[0]] + addition
                    new_branches.append(new_subbranch)
                else:
                    old_subbranch = nested[i]
                    new_branches.append(old_subbranch)
            return np.array(new_branches)

    def deepen_leafcode_with_corr(self, leafcode, nested, corr, popped=None):
        """
        Deepens a tree with copies of the parent.
        The leafcode we receive is 1 longer than the tree is deep.

        PARAMETERS
        ---
        leafcode: list of integers (corresponding to Boolean values)
        nested: a nested list, depth same as length of leafcode
        """

        if len(leafcode) > 1:
            new_branches = []
            for i in range(0, self.numvalues):
                if i == leafcode[0]:
                    new_subbranch = self.deepen_leafcode_with_corr(leafcode[1:],
                                                                   nested[leafcode[0]],
                                                                   corr, leafcode[0])
                    new_branches.append(new_subbranch)
                else:
                    old_subbranch = nested[i]
                    new_branches.append(old_subbranch)
            return list(new_branches)
        else:
            # we should have arrived at a single number
            leaf_value = nested
            new_branch = []

            # fairfrac is the fraction each leaf would get with r = 0
            fairfrac = 1/self.numvalues
            matchfrac = None
            if corr >= 0:
                matchfrac = fairfrac + (1 - fairfrac) * corr
            else:
                matchfrac = fairfrac + fairfrac * corr
            mismatchfrac = (1 - matchfrac) / (self.numvalues - 1)
            for i in range(0, self.numvalues):
                new_leaf = None
                if i == popped:
                    # they match, so correlation applies
                    new_leaf = leaf_value * matchfrac
                else:
                    new_leaf = leaf_value * mismatchfrac
                new_branch = new_branch + [new_leaf]
            return list(new_branch)

    def check_normalization(self):
        """
        Tiny normalization check.

        RETURNS
        ---
        normalization: Boolean
        """
        total_probabilities = np.array(np.minimum(np.maximum(self.joint_probabilities, 0.0), 1.0)
                                       , dtype=np.float128)
        return np.testing.assert_almost_equal(np.sum(total_probabilities), 1.0)

    def half_tree(self):
        """
        Halves the tree, leaving only the new state.
        """
        # find the indices sans the initial state
        indices = range(self.grn_vars["gene_cnt"], self.numvariables)

        # marginalize the distribution
        marginalized_object = self.marginalize_distribution(indices)

        # set the new properties
        self.joint_probabilities = marginalized_object.joint_probabilities
        self.numvariables = marginalized_object.numvariables

    def decide_outcome(self, tally, rule, old_value):
        """
        Decides what happens when multiple rules affect a gene.
        The most basic rules are majority rules, and downregulation dominates.

        PARAMETERS
        ---
        tally: how many times each outcome was decided by all rules (dict)
        rule: the way we decide the outcome
        old_value: the old state of the gene (int)

        RETURNS
        ---
        output_value: the new state of the gene (int)
        """
        if rule == 'totaleffect':
            # find the sum
            total_effect = 0
            for i in range(2*(-self.numvalues), 2*(self.numvalues + 1)):
                total_effect += i * tally[str(i)]

            # calculate the new value
            new_value = old_value + total_effect
            if new_value < 0:
                return 0
            elif new_value >= self.numvalues:
                return self.numvalues - 1
            else:
                return new_value
        if rule == 'down':
            # find the most negative effect
            effect = 0
            for i in range(2*(-self.numvalues), 2*(self.numvalues + 1)):
                if tally[str(i)] != 0:
                    effect = i
                    break

            # calculate the new value
            new_value = old_value + effect
            if new_value < 0:
                return 0
            elif new_value >= self.numvalues:
                return self.numvalues - 1
            else:
                return new_value
        elif rule == 'majority':
            # find the most common effect
            effect = 0
            votes_effect = 0
            options = [i for i in range(2*(-self.numvalues), 2*(self.numvalues + 1))]

            # to randomly pick in case of ties, we shuffle the order and pick
            # the first we encounter in case of ties
            random.shuffle(options)
            for i in options:
                if tally[str(i)] > votes_effect:
                    effect = i
                    votes_effect = tally[str(i)]

            # calculate the new value
            new_value = old_value + effect
            if new_value < 0:
                return 0
            elif new_value >= self.numvalues:
                return self.numvalues - 1
            else:
                return new_value
        elif rule == 'average':
            # find the average affect, then round to the nearest integer
            options = [i for i in range(2*(-self.numvalues), 2*(self.numvalues + 1))]

            total = 0
            divide_by = 0
            for i in options:
                total += i * tally[str(i)]
                divide_by += tally[str(i)]
            if divide_by > 0:
                return max(0, min(int(round(total/divide_by)), self.numvalues - 1))
            else:
                raise ValueError("divide by zero")
        elif rule == 'odds':
            raise ValueError("stochastic method is not re-implemented...")
            #diceroll = np.random.randint(0, sum(tally))
            #running_total = 0
            #for index in range(0, len(tally)):
            #    running_total += tally [index]
            #    if running_total > diceroll:
            #        return index
        else:
            raise ValueError("no valid rule defined!")


    def append_variable_grn(self, gene_index, transition_functions):
        """
        Can be removed if the transition table method works.

        Appends a GRN-ruled variable.

        PARAMETERS
        ---
        gene_index: index of the gene that should be added at the next timestep
        """
        # construct the leafcodes
        states = list(range(self.numvalues))
        leafcodes = list(itertools.product(states, repeat=self.numvariables))

        # make an empty state_transitions object
        state_transitions = []

        # fill the transitions table
        for _leafcode in leafcodes:
            # tally over all rules what state this should be
            # this is always deterministic
            tally = {}
            for i in range(2*(-self.numvalues), 2*(self.numvalues + 1)):
                tally[str(i)] = 0

            # some trickery to make sure you can also append a state
            # when you did not clear the initial one
            _leafcode_original = deepcopy(_leafcode)
            _leafcode = _leafcode[-self.grn_vars["gene_cnt"]:]

            # loop over all relevant rules
            for _func in transition_functions:
                # prepare inputs
                inputs = []
                for index_input in _func["inputs"]:
                    inputs.append(_leafcode[index_input])

                outputs = []
                for index_output in _func["outputs"]:
                    outputs.append(_leafcode[index_output])

                # figure out the output state
                output_value_func = _func["rulefunction"](inputs)

                # add to the tally
                tally[str(int(output_value_func))] += 1

            # decide what it will be this leafcode
            output_value = self.decide_outcome(tally,
                                               self.grn_vars["conflict_rule"],
                                               _leafcode[gene_index])

            # update the leafcode
            _leafcode_original = list(_leafcode_original) + [output_value]

            # add to the state transition
            state_transitions.append(_leafcode_original)

        # adjusting state transitions
        self.append_variables_using_state_transitions_table(np.array(state_transitions))

        # validate everything is still normalized
        self.check_normalization()

    def evaluate_motif(self, genes=None, cleanup_first=True):
        """
        At every leaf, this builds a new probability tree on the initial one.
        This new tree is created based on the rules, and is sparse in nature.
        The sparseness is because the rules are deterministic (for the most part).
        The only non-determinism **COULD** be introduced when a gene is both
        stimulated and regulated.
        FOR NOW, HOWEVER, I ASSUME THE CLOSEST TO INHIBITION DOMINATES.
        The original object at t+1 can be obtained by summing over all the
        trees that are appended to the old tree.

        PARAMETERS
        ---
        genes: list of integers (correspond to gene indices)
        """
        if cleanup_first:
            # drop old timestep
            if self.numvariables == 2 * self.grn_vars["gene_cnt"]:
                self.half_tree()
            elif self.numvariables != self.grn_vars["gene_cnt"]:
                raise ValueError("cannot do a timestep after partial timestep")
        else:
            if self.numvariables == 3 * self.grn_vars["gene_cnt"]:
                # find the indices sans the initial state
                indices = range(0, self.grn_vars["gene_cnt"]) +\
                          range(2*self.grn_vars["gene_cnt"], self.numvariables)

                # marginalize the distribution
                marginalized_object = self.marginalize_distribution(indices)

                # set the new properties
                self.joint_probabilities = marginalized_object.joint_probabilities
                self.numvariables = marginalized_object.numvariables
            elif self.numvariables % self.grn_vars["gene_cnt"] != 0:
                raise ValueError("cannot do a timestep after partial timestep")

        # we don't need to update all genes, but the default is we do
        if genes is None:
            genes = list(range(0, self.grn_vars["gene_cnt"]))

        if self.transition_table is None and self.evaluation_style == 'network':
            for _gene in genes:
                # filter rules with this gene as the target
                transition_functions = []
                for _rule in self.grn_vars["rules"]:
                    # check if gene is the target, if so add
                    if _gene in _rule["outputs"]:
                        transition_functions.append(_rule)
                self.append_variable_grn(_gene, transition_functions)
        elif self.transition_table is not None and (self.evaluation_style == 'transition_table' or self.evaluation_style == 'network'):
            if len(genes) < self.grn_vars["gene_cnt"]:
                raise ValueError("Transition tables only support full system evaluations")
            # adjusting state transitions
            self.append_variables_using_state_transitions_table(self.transition_table)
            self.check_normalization()
        else:
            raise ValueError("Invalid method of getting a transition table, should be transition_table or network")

        # add to the list of states
        ## find the indices sans the initial state
        indices = range(self.grn_vars["gene_cnt"], self.numvariables)
        ## marginalize the distribution
        marginalized_object = self.marginalize_distribution(indices)
        ## append the state
        self.states.append(marginalized_object.joint_probabilities.joint_probabilities)

    def reset_to_state(self, index):
        """
        Reset to the state at the index.

        PARAMETERS
        ---
        index: integer
        """
        if index > (len(self.states) - 1):
            raise ValueError("state does not exist")

        self.joint_probabilities.joint_probabilities = self.states[index]
        self.numvariables = self.grn_vars["gene_cnt"]

    def set_transition_table(self, rule='totaleffect'):
        """
        Find the state transition table, and set it. 
        Note that this always uses the total effect method if not otherwise specified!
        """
        if self.evaluation_style != 'network':
            raise ValueError("can only do this step for network motifs")

        # set up the table
        table = []

        # construct the leafcodes, each leafcode represents a state at T=0
        states = list(range(self.numvalues))
        leafcodes = list(itertools.product(states, repeat=self.grn_vars["gene_cnt"]))
        genes = list(range(0, self.grn_vars["gene_cnt"]))

        # loop over all states
        for _leafcode in leafcodes:
            # we build a new leafcode, which is the new state
            # we build it as a list for table purposes
            leafcode_old = [str(e) for e in list(_leafcode)]
            leafcode_new = []

            for _gene in genes:
                # filter rules with this gene as the target
                transition_functions = []
                for _rule in self.grn_vars["rules"]:
                    # check if gene is the target, if so add
                    if _gene in _rule["outputs"]:
                        transition_functions.append(_rule)

                # tally over all rules what state this should be
                # this is always deterministic
                tally = {}
                for i in range(2 * (-self.numvalues), 2 * (self.numvalues + 1)):
                    tally[str(i)] = 0

                # loop over all relevant rules
                for _func in transition_functions:
                    # prepare inputs
                    inputs = []
                    for index_input in _func["inputs"]:
                        inputs.append(_leafcode[index_input])

                    outputs = []
                    for index_output in _func["outputs"]:
                        outputs.append(_leafcode[index_output])

                    # figure out the output state
                    if set(_func["inputs"]) == set(_func["outputs"]):
                        # this means it is self-activating
                        output_value_func = _func["rulefunction"]([self.numvalues-1])
                    else:
                        output_value_func = _func["rulefunction"](inputs)

                    # add to the tally
                    tally[str(int(output_value_func))] += 1

                output_value = self.decide_outcome(tally, rule, _leafcode[_gene])
                leafcode_new.append(output_value)
            row = []
            row.extend(leafcode_old)
            row.extend(leafcode_new)

            # cast all to integers
            for i in range(0, len(row)):
                row[i] = int(row[i])
            table.append(row)

        # set the state
        self.transition_table = table

    def is_cyclical(self, max_cycle_length=10):
        """
        Tests if the motif contains a cycle of max. n timesteps,
        where n is given as a parameter "max_cycle_length".
        We use the transition table extensively to do this.

        :parameter
        :param max_cycle_length: the maximum length of the cycle (integer)

        :returns
        :return cycles: list of all cycles
        """
        # create the return variable
        cycles = []

        # we loop over all possible starting states
        for row in self.transition_table:
            half_length = len(row)/2
            state_start = row[:half_length]

            # we now loop, making a time step every time
            # we start at the start state, to preserve the start state we copy by values
            state_current = deepcopy(state_start)
            for i in range(0, max_cycle_length):
                state_next = self.match_row(state_current)

                # use numpy to compare
                a = np.array(state_start, dtype=np.int16)
                b = np.array(state_current, dtype=np.int16)
                c = np.array(state_next, dtype=np.int16)

                if np.array_equal(b, c):
                    # this is a static state, so we can quit looking
                    cycles.append([b])
                    break

                # check if we visited the state already in a cycle, if so this is part of the same cycle
                found_cycle = False
                for cycle in cycles:
                    for state in cycle:
                        if np.array_equal(c, state):
                            cycle.append(a)
                            found_cycle = True

                if found_cycle:
                    break

                if np.array_equal(a, c):
                    # this is a loop!
                    cycles.append([a])
                    break
                else:
                    state_current = state_next

        # the number of states that are in a cycle is len(cyclical_states)
        return cycles

    def match_row(self, state_current):
        """
        Find the next state, according to the transition table.

        :parameter
        :param state_current: the current state (list of integers)

        :returns
        :return state_next: the next state (list of integers)
        """
        state_next = None
        for row in self.transition_table:
            half_length = len(row)/2
            row_state_start = row[:half_length]
            row_state_next = row[half_length:]

            # use numpy to compare
            a = np.array(row_state_start, dtype=np.int16)
            b = np.array(state_current, dtype=np.int16)
            if np.array_equal(a, b):
                state_next = row_state_next
                break

        # if nothing is found, return None
        return state_next

    def find_in_sample(self, sample, strict=True):
        """
        See if an identical motif to self can be found in the sample.
        This is tricky, as we want to find every motif that in essence is the same...
        This means motif where 1 stimulates 2 is the same as a motif where 2 stimulates 1.

        :parameter
        :param sample: the set to compare to (list of DiscreteGrnMotif objects)
        :param strict: if True, we ignore matches with the same transition table,
            but different functions (bool)
        :returns
        :return found: the matches of the search (list of DiscreteGrnMotif objects)
        """
        # create the return variable
        found = []

        # loop over our sample
        for motif in sample:
            # first we check if the the functions are the same, if not it is a different motif
            if not strict or (self.grn_vars["rules"].sorted() == motif.grn_vars["rules"].sorted()):
                # find all transition table mutations
                permutations = self.find_transition_table_permutations()
                for permutation in permutations:
                    # to be sure, we compare numpy arrays
                    a = np.array(permutation, dtype=np.int16)
                    b = np.array(motif.transition_table, dtype=np.int16)
                    if np.array_equal(a, b):
                        found.append(motif)

        return found

    def find_transition_table_permutations(self):
        """
        Find all permutations of a transition table that are the same motif.

        :returns
        :return permutations: list of transition tables
        """
        # create return variable
        permutations = []

        # find all possible gene mappings (e.g. 0 > 1 and 1 > 0)
        genes = list(range(0, self.grn_vars["gene_cnt"]))
        genes_permutations = itertools.permutations(genes)

        # loop over all possible mappings, n! possibilities sadly
        for genes_permutation in genes_permutations:
            permutations.append(self.permute_transition_table(genes_permutation))

        return permutations            

    def permute_transition_table(self, mapping):
        """
        Permutes a transition table following a mapping.

        :parameter
        :param mapping: gene index becomes gene value of index (list of integers)

        :returns
        :return table_new: new transistion table (list of list of ints)
        """
        # create return value
        table_old = self.transition_table
        table_new = []

        for row_old in table_old:
            # to keep the order correct
            half_length = len(row_old)/2
            left_new = row_old[:half_length]

            # create and fill left_old
            left_old = [0] * len(left_new)
            for i in range(0, len(left_new)):
                left_old[mapping[i]] = left_new[i]

            # match the correct right_old
            right_old = self.match_row(left_old)

            # create and fill right_new
            right_new = [0] * len(right_old)
            for i in range(0, len(right_old)):
                right_new[i] = right_old[mapping[i]]

            # accumulate the new row
            row_new = np.concatenate((left_new, right_new), axis=0)
            table_new.append(row_new.tolist())

        # permutation done!
        return table_new
