'''
This is a small test for the discrete motif code.
'''

from __future__ import print_function

from discrete_motif import DiscreteGrnMotif
import discrete_motif_functions as functions
import discrete_motif_operations as operations
import discrete_motif_measures as measures
import discrete_motif_generator as generator

def main():
    '''
    This code will be mostly the main.
    '''
    # FIRST, TRY A STANDARD SYNERGY EXAMPLE
    print('EXAMPLE: and X-OR of 2 genes that targets the second gene...\n')

    # create a network motif, we always start with 1 variable
    # as zero variables is unsupported
    motif_xor = DiscreteGrnMotif(1, 2, 'random')

    # add some variables
    motif_xor.grn_vars['gene_cnt'] = 2

    # add some rule
    motif_xor.append_rule([0, 1], [1], functions.xor)

    # construct
    motif_xor.construct_grn()

    # apply the rule
    motif_xor.evaluate_motif(genes=[0, 1])

    # print to make sure it is happening right
    print('KLOPT NIET DE MARGINALISEER')
    print('First tree: ')
    print(motif_xor.states[0])
    print('Second tree: ')
    print(motif_xor.states[1])

    # test the measures
    print('The mutual information: ')
    print(measures.mutual_information(motif_xor))
    print('The WMS information: ')
    print(measures.synergy_wms(motif_xor))
    print('The Rick information: ')
    print(measures.synergy_quax(motif_xor))

    # test decay of information
    print('The memory: ')
    print(measures.mi_decay(motif_xor))

    # nudge the original
    motif_xor.reset_to_state(0)
    operations.nudge_variable(motif_xor, 2, 0.5)
    motif_xor.evaluate_motif()
    print('The nudge impact: ')
    print(measures.hellinger(motif_xor.states[1], motif_xor.states[-1]))

    # SECOND, TRY A STANDARD NO SYNERGY
    print('\nEXAMPLE: a positive and of 2 genes that targets the second gene...\n')
    motif_and = DiscreteGrnMotif(1, 2, 'random')
    motif_and.grn_vars['gene_cnt'] = 2
    motif_and.append_rule([0, 1], [1], functions.plus_and)
    motif_and.construct_grn()
    motif_and.evaluate_motif(genes=[0, 1])
    print('The mutual information: ')
    print(measures.mutual_information(motif_and))
    print('The WMS information: ')
    print(measures.synergy_wms(motif_and))
    print('The Rick information: ')
    print(measures.synergy_quax(motif_and))
    print('The memory: ')
    print(measures.mi_decay(motif_and))
    motif_and.reset_to_state(0)
    operations.nudge_variable(motif_and, 2, 0.5)
    motif_and.evaluate_motif()
    print('The nudge impact: ')
    print(measures.hellinger(motif_and.states[1], motif_and.states[-1]))

    # THIRD, TRY A RANDOM EXAMPLE (SOME SYNERGY PROBABLY)
    print('\nEXAMPLE: a random network...\n')
    motifs, _ = generator.generate_motifs(2, 2, 1)
    motif_rand = motifs[0]
    motif_rand.evaluate_motif(genes=[0, 1])
    print('The rules: ')
    print(motif_rand.grn_vars['rules'])
    print('The correlations: ')
    print(motif_rand.grn_vars['correlations'])
    print('The mutual information: ')
    print(measures.mutual_information(motif_rand))
    print('The WMS information: ')
    print(measures.synergy_wms(motif_rand))
    print('The Rick information: ')
    print(measures.synergy_quax(motif_rand))
    print('The memory: ')
    print(measures.mi_decay(motif_rand))
    motif_rand.reset_to_state(0)
    operations.nudge_variable(motif_rand, 2, 0.5)
    motif_rand.evaluate_motif()
    print('The nudge impact: ')
    print(measures.hellinger(motif_rand.states[1], motif_rand.states[-1]))

    # Create a synergy profile and a scatterplot in a dedicated file
    print('\nCreate a synergy profile and a scatterplot in a dedicated file...')


if __name__ == '__main__':
    main()
