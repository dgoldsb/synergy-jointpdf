"""
This is a small test for the discrete motif code.
"""

from __future__ import print_function

from discrete_motif import DiscreteGrnMotif
import discrete_motif_functions as functions
import discrete_motif_operations as operations
import discrete_motif_measures as measures
import discrete_motif_generator as generator

def main():
    """
    This code will be mostly the main.
    """
    # create a network motif, we always start with 1 variable
    # as zero variables is unsupported
    motif = DiscreteGrnMotif(1, 2, 'random')

    # add some variables
    motif.grn_vars["gene_cnt"] = 2

    # add some rule
    motif.append_rule([0, 1], [1], functions.xor)

    print("The state empty: ")
    print(motif.joint_probabilities.joint_probabilities)

    # construct
    motif.construct_grn()

    # OVERRIDE
    motifs = generator.generate_motifs(1, 2, 0.6)
    motif = motifs[0]

    print("The correlations: ")
    print(motif.grn_vars["correlations"])

    print("The state before: ")
    print(motif.joint_probabilities.joint_probabilities)

    # print the rules
    print("The rules: ")
    print(motif.grn_vars["rules"])

    # test the rule
    motif.evaluate_motif(genes=[0, 1])
    print("The full state after: ")
    print(motif.joint_probabilities.joint_probabilities)

    # test the measures
    print("The absolute difference: ")
    print(measures.abs_diff(motif.states[0], motif.states[1]))
    print("The Hellinger distance: ")
    print(measures.hellinger(motif.states[0], motif.states[1]))
    if measures.hellinger(motif.states[0], motif.states[1]) is not None:
        print("Pretty different, but since we do time evolution that makes sense!")
    print("The mutual information: ")
    print(measures.mutual_information(motif))
    print("The synergistic information: ")
    print(measures.synergy_quax(motif))
    print("The WMS information: ")
    print(measures.synergy_wms(motif))

    # testing the nudges
    # reset to the first state

    # perform a nudge by going back to the state, and comparing the new outcome
    #print("The nudge impact: ")

    # get the half tree and be ready for a next timestep
    motif.half_tree()
    print("The marginalized state after: ")
    print(motif.joint_probabilities.joint_probabilities)

    # OK! now reset and try simple X-OR

    # OK! now reset and try simple copy case (no rules)

    # Create a synergy profile

    # Create a scatterplot

    # Now try to create a correlation matrix


if __name__ == '__main__':
    main()
