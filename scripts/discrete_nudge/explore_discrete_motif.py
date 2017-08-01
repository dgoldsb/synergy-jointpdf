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
    motif_set = DiscreteGrnMotif(1, 2, 'random')

    # add some variables
    motif_set.grn_vars["gene_cnt"] = 2

    # add some rule
    motif_set.append_rule([0, 1], [1], functions.xor)

    print("The state empty: ")
    print(motif_set.joint_probabilities.joint_probabilities)

    # construct
    motif_set.construct_grn()

    # OVERRIDE because we can
    motifs, _ = generator.generate_motifs(2, 3, 2)
    print(motifs)
    motif = motifs[0]

    print("The correlations: ")
    print(motif.grn_vars["correlations"])

    print("The state before: ")
    print(motif.joint_probabilities.joint_probabilities)

    # print the rules
    print("The rules: ")
    print(motif.grn_vars["rules"])

    # test the rule
    print("The full state after: ")
    motif.evaluate_motif(genes=[0, 1])
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
    print("The WMS information: ")
    print(measures.synergy_wms(motif))
    #print("The synergistic information: ")
    #print(measures.synergy_quax(motif))

    # testing the nudges, let's use the original
    motif_set.evaluate_motif()

    # reset to the first state
    motif_set.reset_to_state(0)

    # perform a nudge by going back to the state
    operations.nudge_variable(motif_set, 1, 0)
    motif_set.evaluate_motif()

    # compare the outcome
    print("The nudge impact: ")
    print(measures.hellinger(motif_set.states[-2], motif_set.states[-1]))

    # get the half tree and be ready for a next timestep
    print("The marginalized state after: ")
    motif.half_tree()
    print(motif.joint_probabilities.joint_probabilities)

    # examine the memory
    print("The memory: ")
    motif = motifs[1]
    print(measures.mi_decay(motif))

    # Create a synergy profile and a scatterplot in a dedicated file
    print("Create a synergy profile and a scatterplot in a dedicated file...")

if __name__ == '__main__':
    main()
