"""
This file builds a discrete GRN motif, and contains some of the common operations used.
"""

from __future__ import print_function
#import os
#import sys
import itertools
import math

#import motif_functions as functions
from jointpdf.jointpdf import JointProbabilityMatrix
#from jointpdf.jointpdf import FullNestedArrayOfProbabilities
#sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/jointpdf_DJ")
#from code import nudge

class DiscreteGrnMotif(object):
    """
    A motif which consists of:

    - a set of Genes, each assigned a number
    - a set of dependent PMFs that define probabilities that each gene
      is activated (FullNestedArrayOfProbabilities object)
    - a Boolean rule to determine states of the outputs between the genes

    We can advance time with one timestep, and we can use measures such as the KL-divergence.
    It is based on the jointpdf package from Rick Quax.
    """
    def __init__(self, numvalues=2):
        """
        Creates a FullNestedArrayOfProbabilities object.

        Parameters
        ---
        numvalues: integer
        """
        self.state = JointProbabilityMatrix(numvariables=1, numvalues=numvalues
                                            , joint_probs='random')
        self.rules = []
        self.numgenes = 1
        self.numvalues = numvalues

    def append_gene(self):
        """
        Append a gene to the motif.
        """
        self.state.append_variables(num_added_variables=1, added_joint_probabilities=None)
        self.numgenes = self.numgenes + 1

    def append_rule(self, inputs, outputs, rulefunction):
        """
        Append rule to the motif, keep in mind that the order of the inputs could matter.
        The function defines a truthtable, which is part of the motif operations.
        Rulefunctions are applied in a set order, so it does not matter in which order you add them.

        Parameters
        ---
        inputs: list of gene indices (integers)
        outputs: list of gene indices (integers)
        rulefunction: function to call to evaluate this rule
        """

        # check if the gene actually exists
        for _input in inputs:
            if (_input < self.numgenes) is False:
                raise RuntimeError("you are referencing a non-existing gene")
        for _output in outputs:
            if (_output < self.numgenes) is False:
                raise RuntimeError("you are referencing a non-existing gene")

        # we do not validate if the number of inputs/outputs match the rule
        # add the rule
        rule = {}
        rule["inputs"] = inputs
        rule["outputs"] = outputs
        rule["rulefunction"] = rulefunction
        self.rules.append(rule)

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

    def evaluate_motif(self):
        """
        This rebuilds the probability tree based on the initial state and motif.
        Afterwards, synergy estimations can be done.
        Essentially, we add some determinism to the tree.
        With promotion, the target gene can either be always +, or the previous state.
        This depends on whether the promotor is active or not.
        """
        # make a copy of the state
        new_state = self.state.joint_probabilities.joint_probabilities.copy()

        # order rules: first do positive, then negative
        # result is that negative dominates
        for _rule in sorted(self.rules, key=lambda x: x["rulefunction"]):
            states = list(range(self.numvalues)) # for later, if I add 3 states
            leafcodes = list(itertools.product(states, repeat=self.numgenes))
            for _leafcode in leafcodes:
                # prepare inputs
                inputs = []
                for index_input in _rule["inputs"]:
                    inputs.append(_leafcode[index_input])

                # prepare outputs
                outputs = []
                for index_output in _rule["outputs"]:
                    outputs.append(_leafcode[index_output])

                # see if the outcome is now deterministic
                new_outputs = _rule["rulefunction"](inputs, outputs)
                if new_outputs != outputs:
                    # this output is apparently no longer possible
                    # add this value to the opposite state, so with flipped output
                    addition = self.fetch_leafcode(_leafcode, new_state)

                    # now make the actual state zero
                    new_state = self.add_leafcode(_leafcode, new_state, -addition)

                    # flip the outputs
                    _leafcode = list(_leafcode)
                    _leafcode_prev = list(_leafcode)
                    for index_output in _rule["outputs"]:
                        _leafcode[index_output] = int(math.fabs(_leafcode[index_output] - 1))
                    # add to the state with flipped output
                    print("Moving "+str(addition)+" from "
                          +str(_leafcode_prev)+" to "+str(_leafcode))
                    new_state = self.add_leafcode(_leafcode, new_state, addition)

                # calculate normalization
                total_probability = new_state
                for _ in range(0, len(_leafcode)):
                    total_probability_reduced = total_probability[0]
                    for j in range(1, self.numvalues):
                        total_probability_reduced = total_probability_reduced + total_probability[j]
                    total_probability = total_probability_reduced

                # check normalization
                if (total_probability < 0.99) or (total_probability > 1.01):
                    print("Total probability: "+str(total_probability))
                    print("THIS SEEMS TO BE NOT WORKING IF THERE IS MORE THAN ONE RULE, FIX")
                    #raise RuntimeError("the PDF is no longer normalized")
        
        self.state.joint_probabilities.joint_probabilities = new_state


