"""
This is a small test for the discrete motif code.
"""

from __future__ import print_function
import sys
from discrete_motif import DiscreteGrnMotif
import discrete_motif_functions as functions


def main():
    """
    This code will be mostly the main.
    """
    # some presets
    motif_size = 4

    # create a network motif, in our case a bifan
    motif = DiscreteGrnMotif()

    # add some variables
    for _ in range(1, motif_size):
        motif.append_gene()
    print("Size of motif: "+str(motif.numgenes))
    print("The state: ")
    print(motif.state.)

    # add some rules
    motif.append_rule([0], [2], functions.plus)
    motif.append_rule([1], [3], functions.plus)
    motif.append_rule([0], [3], functions.minus)
    motif.append_rule([1], [2], functions.minus)
    print("The rules: ")
    print(motif.rules)

if __name__ == '__main__':
    main()
