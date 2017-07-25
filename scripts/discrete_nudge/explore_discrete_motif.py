"""
This is a small test for the discrete motif code.
"""

from __future__ import print_function

import sys
from discrete_motif import DiscreteGrnMotif
import discrete_motif_functions as functions
import discrete_motif_operations as operations
import discrete_motif_measures as measures


def main():
    """
    This code will be mostly the main.
    """
    # create a network motif, we always start with 1 variable
    # as zero variables is unsupported
    motif = DiscreteGrnMotif(1, 2, 'random')

    # add some variables
    motif.grn_vars["gene_cnt"] = 2

    # add some rules
    motif.append_rule([0], [1], functions.plus)

    # construct
    print("The state empty: ")
    print(motif.joint_probabilities.joint_probabilities)
    motif.construct_grn()
    print("The state before: ")
    print(motif.joint_probabilities.joint_probabilities)

    # print the rules
    print("The rules: ")
    print(motif.grn_vars["rules"])

    # test the rule
    motif.evaluate_motif()
    print("The full state after: ")
    print(motif.joint_probabilities.joint_probabilities)

    # test the measures
    print("The absolute difference: ")
    print(measures.abs_diff(motif.states[0], motif.states[1]))
    print("Pretty different, but since we do time evolution that makes sense!")
    print("The mutual information: ")

    print("The synergistic information: ")

    print("The WMS information: ")

    # testing the nudges
    # reset to the first state

    # perform a nudge by going back to the state, and comparing the new outcome
    print("The nudge impact: ")

    # get the half tree and be ready for a next timestep
    motif.half_tree()
    print("The reduced state after: ")
    print(motif.joint_probabilities.joint_probabilities)

if __name__ == '__main__':
    main()
