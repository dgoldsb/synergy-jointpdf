"""
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

    def getter_grn(self):
        """
        Get the settings for the GRN.
        """
        return self.grn_vars

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
        self.reset(self, 0, self.numvalues, self.joint_probabilities.joint_probabilities)

        # add all the genes
        if self.grn_vars["correlations"] is None:
            # use the default append_variable with random dependencies
            self.append_variables(num_added_variables=self.grn_vars["gene_cnt"]
                                  , added_joint_probabilities=None)
        else:
            # use the total correlation append_variable
            self.append_variables_correlation(num_added_variables=self.grn_vars["gene_cnt"]
                                              , correlations_matrix=self.grn_vars["correlations"])

    def append_variables_correlation(self, num_added_variables, correlations_matrix):
        """
        Adds variables, enforcing set total correlations between them.
        """
        print("Unimplemented...")
        sys.exit(0)

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
        """
        if len(leafcode) > 1:
            return self.fetch_leafcode(leafcode[1:], nested[leafcode[0]])
        else:
            return nested[leafcode[0]]

    def add_leafcode(self, leafcode, nested, addition):
        """
        Alters the value at the leafcode in the list.
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
            return list(new_branches)
        else:
            new_branches = []
            for i in range(0, self.numvalues):
                if i == leafcode[0]:
                    new_subbranch = nested[leafcode[0]] + addition
                    new_branches.append(new_subbranch)
                else:
                    old_subbranch = nested[i]
                    new_branches.append(old_subbranch)
            return list(new_branches)

    def check_normalization(self):
        """
        Tiny normalization check.
        """
        total_probabilities = np.array(np.minimum(np.maximum(self.joint_probabilities, 0.0), 1.0)
                                       , dtype=np.float128)
        return np.testing.assert_almost_equal(np.sum(total_probabilities), 1.0)

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
        """
        # we don't need to update all genes, but the default is we do
        if genes is None:
            genes = list(range(0, self.grn_vars["gene_cnt"]))

        for _gene in genes:
            # filter rules with this gene as the target
            rules_applicable = []
            for _rule in self.grn_vars["rules"]:
                # check if _gene is the target, if so add
                if _gene in _rule["outputs"]:
                    rules_applicable.append(_rule)

            # construct the leafcodes
            states = list(range(self.numvalues))
            leafcodes = list(itertools.product(states, repeat=self.numvariables))
            for _leafcode in leafcodes:
                # tally over all rules what state this should be
                # this is always deterministic
                tally = [0] * self.numvalues

                # loop over all relevant rules
                for _rule in rules_applicable:
                    # prepare inputs
                    inputs = []
                    for index_input in _rule["inputs"]:
                        inputs.append(_leafcode[index_input])

                    # figure out the output state
                    output_value_rule = _rule["rulefunction"](inputs)

                    # add to the tally
                    tally[output_value_rule] = tally[output_value_rule] + 1

                # decide what it will be this leafcode
                output_value = 0
                for _ in range(0, self.numvalues):
                    if tally[output_value] != 0:
                        break
                    else:
                        output_value = output_value + 1

                # extend tree
                print("Unimplemented...")
                sys.exit(0)

            # validate everything is still normalized
            self.check_normalization()
