"""
This file builds a discrete GRN motif, and contains some of the common operations used.
"""

#import os
#import sys

#import motif_functions as functions
from jointpdf.jointpdf import JointProbabilityMatrix
#from jointpdf.jointpdf import FullNestedArrayOfProbabilities
#sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/jointpdf_DJ")
#from code import nudge

class DiscreteGrnMotif(object):
    """
    A motif which consists of:

    - a set of Genes, each assigned a number
    - a set of dependent PMFs that define probabilities that each gene is activated
    - a Boolean rule to determine states of the outputs between the genes

    We can advance time with one timestep, and we can use measures such as the KL-divergence.
    It is based on the jointpdf package from Rick Quax.
    """
    def __init__(self, numvalues=1):
        self.state = JointProbabilityMatrix(numvariables=1, numvalues=numvalues, joint_probs='random')
        self.rules = []
        self.numgenes = 1

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

        inputs: list of gene indices
        outputs: list of gene indices
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

    def evaluate(self):
        """
        Evaluate timesteps
        """
        # order rules: first do positive, then negative
        # result is that negative dominates
        return 0


