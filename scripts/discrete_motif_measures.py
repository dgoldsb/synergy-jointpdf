'''
This file contains some of the common measures used on discrete motifs.
'''

from __future__ import print_function

from copy import deepcopy
import itertools
import math
import numpy as np
from scipy.linalg import norm
from scipy.optimize import minimize
from scipy.stats import entropy

import discrete_motif_operations as operations

def total_correlation(motif, indices):
    '''
    Returns the total correlation.

    @param motif: a DiscreteGrnMotif object
    @param indices: the indices over which to calculate the total correlation (list of integers)
    @returns: the total correlation (float)
    '''
    return_value = - motif.entropy(indices)
    for index in indices:
        return_value += motif.entropy(index)
    return return_value


def abs_diff(tree_1, tree_2):
    '''
    Finds the absolute difference between two same-shaped FullNestedArrayOfProbabilities objects.
    We chose this, as the Kullback-Leibler divergence cannot handle zeros.

    @param tree_1: FullNestedArrayOfProbabilities object
    @param tree_2: FullNestedArrayOfProbabilities object
    @returns: absolute difference (float)
    '''
    # flatten the trees
    t1_flat = np.array(np.copy(tree_1).flatten())
    t2_flat = np.array(np.copy(tree_2).flatten())

    if len(t1_flat) != len(t2_flat):
        return None

    returnval = 0
    for i in range(0, len(t1_flat)):
        returnval = math.fabs(t1_flat[i] - t2_flat[i])

    return returnval


def hellinger(tree_1, tree_2):
    '''
    Finds the absolute difference between two same-shaped FullNestedArrayOfProbabilities objects.
    We chose this, as the Kullback-Leibler divergence cannot handle zeros.

    @param tree_1: FullNestedArrayOfProbabilities object
    @param tree_2: FullNestedArrayOfProbabilities object
    @returns: hellinger divergence (float)
    '''
    # flatten the trees
    t1_flat = np.array(np.copy(tree_1).flatten())
    t2_flat = np.array(np.copy(tree_2).flatten())

    if len(t1_flat) != len(t2_flat):
        return None

    return norm(np.sqrt(t1_flat) - np.sqrt(t2_flat)) / np.sqrt(2)


def kl_div(tree_1, tree_2):
    '''
    Finds the KL-divergence for the motif after a time evaluation.

    @param motif: a DiscreteGrnMotif object
    @returns: KL divergence (float)
    '''
    # flatten the trees
    t1_flat = np.array(np.copy(tree_1).flatten())
    t2_flat = np.array(np.copy(tree_2).flatten())

    return entropy(t1_flat, t2_flat)


def mutual_information(motif, genes_t0 = None, genes_t1 = None):
    '''
    Finds the mutual information for the motif after a time evaluation.

    @param motif: a DiscreteGrnMotif object
    @param genes_t0 (list of variables to compare to genes_t1)
    @param genes_t1 (list of variables to compare to genes_t0)
    @returns: mutual information (float)
    '''
    if motif.grn_vars['gene_cnt'] == motif.numvariables:
        raise ValueError('do a time evaluation before attempting to check MI')

    # define what t0 and t1 are
    if genes_t0 is None:
        genes_t0 = range(0, motif.grn_vars['gene_cnt'])
    if genes_t1 is None:
        genes_t1 = range(motif.grn_vars['gene_cnt'], motif.numvariables)

    return motif.mutual_information(genes_t0, genes_t1)


def synergy_quax(motif, tolerance=0.2):
    '''
    Finds the mutual information for the motif after a time evaluation.

    @param motif: a DiscreteGrnMotif object
    @returns: upper bound of synergy according to Rick Quax (float)
    '''
    if motif.grn_vars['gene_cnt'] == motif.numvariables:
        raise ValueError('do a time evaluation before attempting to check synergy')

    # define what t0 and t1 are
    genes_t0 = range(0, motif.grn_vars['gene_cnt'])
    genes_t1 = range(motif.grn_vars['gene_cnt'], motif.numvariables)

    return motif.synergistic_information(genes_t1, genes_t0, tol_nonsyn_mi_frac=tolerance, verbose=False,
                                         minimize_method=None, num_repeats_per_srv_append=30,
                                         initial_guess_summed_modulo=False)


def synergy_uii(motif, tolerance=0.2):
    '''
    Finds the unique individual information

    @param motif: a DiscreteGrnMotif object
    @param tolerance: the allowed tolerance (float)
    @returns: uii according to Rick Quax (float)
    '''
    if motif.grn_vars['gene_cnt'] == motif.numvariables:
        raise ValueError('do a time evaluation before attempting to check synergy')

    # define what t0 and t1 are
    genes_t0 = range(0, motif.grn_vars['gene_cnt'])
    genes_t1 = range(motif.grn_vars['gene_cnt'], motif.numvariables)

    return motif.unique_individual_information(genes_t1, genes_t0, tol_nonunq=tolerance, verbose=False,
                                              num_repeats_per_append=3, assume_independence=False)


def synergy_wms(motif):
    '''
    Finds the mutual information for the motif after a time evaluation.

    @param motif: a DiscreteGrnMotif object
    @returns: WMS synergy (float)
    '''
    if motif.grn_vars['gene_cnt'] == motif.numvariables:
        raise ValueError('do a time evaluation before attempting to check synergy')

    # define what t0 and t1 are
    genes_t0 = range(0, motif.grn_vars['gene_cnt'])
    genes_t1 = range(motif.grn_vars['gene_cnt'], motif.numvariables)

    return motif.synergistic_information_naive(genes_t1, genes_t0)


def synergy_middleground(motif):
    '''
    A synergy approximation by taking the mean of the upper bound (single MI with whole system)
    and lower bound (WMS)

    @param motif: a DiscreteGrnMotif
    @returns: middle ground synergy approximation (float)
    '''
    if motif.grn_vars['gene_cnt'] == motif.numvariables:
        raise ValueError('do a time evaluation before attempting to check synergy')

    uppers = []
    genes = list(range(0, motif.grn_vars['gene_cnt']))
    for _gene in genes:
        uppers.append(mutual_information(motif, genes_t0=[_gene]))
    upper = mutual_information(motif) - max(uppers)
    lower = synergy_wms(motif)
    return (upper + lower) / 2


def mi_decay(motif, no_steps=8):
    '''
    Finds the mutual information for the motif after a time evaluation.

    @param motif: a DiscreteGrnMotif object
    @returns: list of MI state 0 vs. state i (list of floats)
    '''
    memory = [motif.entropy()]

    for i in range(1, no_steps + 1):
        motif.evaluate_motif(cleanup_first=False)
        first = range(0, motif.grn_vars['gene_cnt'])
        current = range(i*motif.grn_vars['gene_cnt'], motif.numvariables)
        memory.append(motif.mutual_information(first, current))

    # reset the motif
    motif.reset_to_state(0)
    motif.states = [deepcopy(motif.joint_probabilities.joint_probabilities)]

    return memory


def average_nudge_impact(motif, nudge_width, nudge_size, nudge_method):
    '''
    Find the average nudge effect over all possible nudges.

    @param motif: a motif object
    @param nudge_width: the number of genes to nudge
    @param nudge_size: the size of the nudge
    @param nudge_method: the nudge method to utilize
    @returns: the average nudge impact
    '''
    impacts = []

    # find all possible targets of this nudge_width
    targets_list = [list(t) for t in itertools.combinations(list(range(0, motif.grn_vars['gene_cnt'])), nudge_width)]

    for targets in targets_list:
        targets.sort()

        # make sure the motif is pristine
        motif.reset_to_state(0)
        motif.states = [deepcopy(motif.joint_probabilities.joint_probabilities)]

        # find the nudge impact
        motif.evaluate_motif()

        # save the unnudged state
        indices = range(motif.grn_vars['gene_cnt'], motif.numvariables)
        marginalized_object = motif.marginalize_distribution(indices)
        unnudged = marginalized_object.joint_probabilities.joint_probabilities

        # reset and find the nudged state
        motif.reset_to_state(0)
        operations.nudge_variable(motif, targets, nudge_size, nudge_method)
        motif.evaluate_motif()

        # we compare the two evolved states
        impact = hellinger(unnudged, motif.states[-1])
        impacts.append(impact)
    
    # reset the motif
    motif.reset_to_state(0)
    motif.states = [deepcopy(motif.joint_probabilities.joint_probabilities)]

    return sum(impacts)/len(impacts)


def normalized_synergy(motif, synergy_measure):
    '''
    Normalized synergy, leaves motif prestine.

    @param motif: a motif object
    @param synergy_measure: the measure to use (function)
    @returns: the synergy (float)
    '''
    # reset
    motif.reset_to_state(0)
    motif.states = [deepcopy(motif.joint_probabilities.joint_probabilities)]

    # find the synergy
    motif.evaluate_motif()
    synergy = synergy_measure(motif)

    # normalize by dividing by the total entropy (not mutual information)
    mutual_information_calculated = float(mutual_information(motif))
    if mutual_information_calculated != 0:
        synergy = float(synergy)/mutual_information_calculated
    else:
        # No mutual information means no synergy following this measure
        synergy = 0

    # reset the motif
    motif.reset_to_state(0)
    motif.states = [deepcopy(motif.joint_probabilities.joint_probabilities)]

    return synergy


def normalized_memory(motif):
    '''
    Normalized synergy, leaves motif prestine.

    @param motif: a motif object
    @returns: the normalized memory (float)
    '''
    # reset
    motif.reset_to_state(0)
    motif.states = [deepcopy(motif.joint_probabilities.joint_probabilities)]
    
    # calculate the memory
    decay = mi_decay(motif, 1)
    memory = decay[1] / decay[0]

    # reset the motif
    motif.reset_to_state(0)
    motif.states = [deepcopy(motif.joint_probabilities.joint_probabilities)]

    return memory
