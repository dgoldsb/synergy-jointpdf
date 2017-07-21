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
    # create a network motif, in our case a bifan
    motif = DiscreteGrnMotif()

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
    print("The state after: ")
    print(motif.joint_probabilities.joint_probabilities)

if __name__ == '__main__':
    main()
