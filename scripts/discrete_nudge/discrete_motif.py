"""
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This file builds a discrete GRN motif, and contains some of the common operations used.
"""

from __future__ import print_function

import itertools
import json
import os
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

    # global variables of a GRN
    grn_vars = {}
    ## important to store the number of genes, as this is not the same as the numvariables
    grn_vars["gene_cnt"] = 0
    ## the rules are not part of the original framework
    grn_vars["rules"] = []
    ## correlation matrix for the dependent PDF of the initial conditions
    grn_vars["correlations"] = None
    ## set the default to downregulation dominates
    grn_vars["conflict_rule"] = 'down'
    ## also there has been a sequence of states, the latest is the current state
    states = []

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

    def deepen_leafcode_old(self, leafcode, nested):
        """
        Can be removed if the transition table method works.

        Deepens a tree.
        The leafcode we receive is 1 longer than the tree is deep.
        Deepens in deterministic fashion: all probability goes to the last
        leafcode digit.

        PARAMETERS
        ---
        leafcode: list of integers (corresponding to Boolean values)
        nested: a nested list, depth same as length of leafcode
        """

        if len(leafcode) > 1:
            new_branches = []
            for i in range(0, self.numvalues):
                if i == leafcode[0]:
                    new_subbranch = self.deepen_leafcode_old(leafcode[1:],
                                                             nested[leafcode[0]])
                    new_branches.append(new_subbranch)
                else:
                    old_subbranch = nested[i]
                    new_branches.append(old_subbranch)
            return list(new_branches)
        else:
            # we should have arrived at a single number
            leaf_value = nested
            new_branch = []

            for i in range(0, self.numvalues):
                if int(i) == int(leafcode[0]):
                    new_branch = new_branch + [leaf_value]
                else:
                    new_branch = new_branch + [0]
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
        tally: how many times each outcome was decided by all rules (list of ints)
        rule: the way we decide the outcome
        old_value: the old state of the gene (int)

        RETURNS
        ---
        output_value: the new state of the gene (int)
        """
        if rule == 'down':
            output_value = 0
            if sum(tally) == 0:
                output_value = old_value
            else:
                for _ in range(0, self.numvalues):
                    if tally[output_value] != 0:
                        break
                    else:
                        output_value = output_value + 1
            return output_value
        elif rule == 'majority':
            return tally.index(max(tally))
        else:
            raise ValueError("no valid rule defined!")

    def append_variable_grn_old(self, gene_index, transition_functions):
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
        for _leafcode in leafcodes:
            # tally over all rules what state this should be
            # this is always deterministic
            tally = [0] * self.numvalues

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
                output_value_func = _func["rulefunction"](inputs, outputs[0])

                # add to the tally
                tally[output_value_func] = tally[output_value_func] + 1

            # decide what it will be this leafcode
            output_value = self.decide_outcome(tally, self.grn_vars["conflict_rule"],
                                               _leafcode[gene_index])

            # update the leafcode
            _leafcode = list(_leafcode) + [output_value]

            # extend tree
            self.joint_probabilities.joint_probabilities =\
                self.deepen_leafcode(_leafcode, self.joint_probabilities.joint_probabilities)

        # add the variable
        self.numvariables = self.numvariables + 1

        # very important: turn the thing back to a numpy array!
        # this was not possible before, as in the process of extending the tree is unbalanced
        self.joint_probabilities.joint_probabilities =\
            np.array(self.joint_probabilities.joint_probabilities)

        # validate everything is still normalized
        self.check_normalization()


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
            tally = [0] * self.numvalues

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
                output_value_func = _func["rulefunction"](inputs, outputs[0])

                # add to the tally
                tally[output_value_func] = tally[output_value_func] + 1

            # decide what it will be this leafcode
            output_value = self.decide_outcome(tally, self.grn_vars["conflict_rule"],
                                               _leafcode[gene_index])

            # update the leafcode
            _leafcode = list(_leafcode) + [output_value]

            # add to the state transition
            state_transitions.append(_leafcode)

        # add usting state transitions
        self.append_variables_using_state_transitions_table(np.array(state_transitions))

        # validate everything is still normalized
        self.check_normalization()

    def evaluate_motif(self, genes=None):
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
        # drop old timestep
        if self.numvariables == 2 * self.grn_vars["gene_cnt"]:
            self.half_tree()
        elif self.numvariables != self.grn_vars["gene_cnt"]:
            raise ValueError("cannot do a timestep: numvariables not equal to gene_cnt or twice")

        # we don't need to update all genes, but the default is we do
        if genes is None:
            genes = list(range(0, self.grn_vars["gene_cnt"]))

        for _gene in genes:
            # filter rules with this gene as the target
            transition_functions = []
            for _rule in self.grn_vars["rules"]:
                # check if gene is the target, if so add
                if _gene in _rule["outputs"]:
                    transition_functions.append(_rule)
            self.append_variable_grn(_gene, transition_functions)

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
