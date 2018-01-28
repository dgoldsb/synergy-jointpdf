"""
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

Generates random discrete motif objects and saves them.
"""

from copy import deepcopy
from fractions import Fraction
import itertools
import math
import numpy as np
from operator import mul    # or mul=lambda x,y:x*y
from pyDOE import *
import random


from discrete_motif import DiscreteGrnMotif
import discrete_motif_functions as functions


def nCk(n,k): 
  return int( reduce(mul, (Fraction(n-i, i+1) for i in range(k)), 1) )


def generate_correlation_matrix(gene_cnt, betaparam = 10):
    """
    Generate a completely random gene correlation matrix.
    We cannot just fill the diagonal band and derive from there, as the we asume the system labels
    represent the flow of information.
    We use the vine method from a paper by Lewandowski, Kurowicka, and Joe.

    PARAMETERS
    ---
    gene_cnt: number of genes (int)

    RETURNS
    ---
    matrix: a gene_cnt x gene_cnt matrix (ndarray of floats)
    """
    P = np.zeros((gene_cnt, gene_cnt))
    S = np.eye(gene_cnt)

    for k in range(0, gene_cnt-1):
        for i in range(k + 1, gene_cnt):
            P[k,i] = np.random.beta(betaparam,betaparam)
            P[k,i] = (P[k,i]-0.5)*2
            p = P[k,i]
            for l in reversed(range(0, k-1)):
                p = p * math.sqrt((1-P[l,i]**2)*(1-P[l,k]**2)) + P[l,i]*P[l,k]
            S[k,i] = p
            S[i,k] = p

    # permute the matrix
    permutation = range(0, gene_cnt)
    random.shuffle(permutation)
    S_permuted = np.eye(gene_cnt)
    
    iterator_i = 0
    for i in permutation:
        iterator_j = 0
        for j in permutation:
            S_permuted[iterator_i, iterator_j] = S[i, j]
            iterator_j += 1
        iterator_i += 1
    
    return S_permuted.tolist()


def generate_rules(no_rules, no_nodes, p_1to1):
    """
    We use a Erdos-Renyi approach, without preferential attachment.
    We only use 1-to-1 and 2-to-1 functions, and randomly select origins
    and targets.
    We can consider every number of inputs a separate ER-network,
    with an accept probability tweaked to get the correct indegree.
    The 2-to-1 network effectively has n**2 nodes.

    The ratio 2-to-1 edges is way too high compared to the ratio of 1-to-1
    edges, the latter are much more common in nature.
    To compensate, we manually turn this to an approximately 25-75 ratio
    (configurable).

    PARAMETERS
    ---
    no_rules: the number of rules requested (int)
    no_nodes: the number of nodes requested (int)
    p_1to1: the chance that an edge is 1-to-1 (float)

    RETURNS
    ---
    rules: a list of rule dictionaries
    """
    # we store all rule-dictionaries in this list
    rules = []

    # we fetch all possible functions that can be attached to an edge
    functions_choices = functions.dictionary()

    # create a list of all possible target genes
    targets = range(0, no_nodes)

    # create list to sotre all actor gene combinations
    actors = []

    # create list of all combinations of genes up to 2
    actors_tuples = []
    for i in range(1, 3):
        actors_tuples += itertools.combinations(targets, i)
        actors = [list(t) for t in actors_tuples]

    # create a list of all possible in/out combinations
    edges = []
    for i in range(0, len(actors)):
        for j in range(0, len(targets)):
            combination = []
            combination.append(deepcopy(actors[i]))
            combination.append(deepcopy(targets[j]))
            edges.append(combination)

    # a list to store the accept probability for the number of input
    # first values is 1-to-1, second 2-to-1 etcetera
    accept_ps = [0.0] * 2
    # indegrees for both number of inputs, they are essentially separate networks we add
    ks = [0.0] * 2

    # determine the indegrees for both networks
    ks[0] = (no_rules/no_nodes) * p_1to1
    ks[1] = (no_rules/no_nodes) * (1 - p_1to1)

    # determine the probabiliti es for both networks
    accept_ps[0] = ks[0] / no_nodes
    accept_ps[1] = ks[1] / nCk(no_nodes, 2)

    # loop over the edges and accept or reject
    for edge in edges:
        # draw random number
        draw = np.random.random()

        # check if we accept
        no_inputs = len(edge[0])
        if draw < accept_ps[no_inputs - 1]:
            # pick a random rule
            # find the number of inputs and outputs required
            choice = random.choice(functions_choices[no_inputs])

            # set the values
            rule = {}
            rule["inputs"] = deepcopy(edge[0])
            rule["outputs"] = [deepcopy(edge[1])]
            rule["rulefunction"] = choice["f"]

            # append to rules
            rules.append(rule)
    return rules


def generate_motifs(samplesize, no_nodes, numvalues=2, indegrees=None, conflict_rule='totaleffect', chance_1to1=0.75):
    """
    Returns a list of objects.
    Improvable Barabasi-Albert network.
    This returns networks, which represent a part of the transition table search
    space that we think contains all biological GRN networks.

    PARAMETERS
    ---
    samplesize: number of objects in list
    no_nodes: number of genes
    indegrees: list of possible indegree desired (average indegree will tend
    towards the average of the list), if None random (list of integers)
    numvalues: number of possible values (int)
    conflict_rule: the rule for deciding what state we get when rules conflict (str)

    RETURNS
    ---
    motifs: list of DiscreteGrnMotif objects
    """
    # create empty list for the objects
    motifs = []

    # create list of all possible target genes
    targets = range(0, no_nodes)

    # create list of all single actor genes
    actors = []
    actors_tuples = []
    for i in range(1, 3):
        actors_tuples += itertools.combinations(targets, i)
    actors = [list(t) for t in actors_tuples]

    # for calculating the average indegree
    rules_total = 0

    # generate N samples
    for _ in range(0, samplesize):
        # create
        grn_vars = {}

        # set the size
        grn_vars["gene_cnt"] = no_nodes
        grn_vars["conflict_rule"] = conflict_rule
        # set the correlation matrix
        grn_vars["correlations"] = generate_correlation_matrix(no_nodes)

        # generate a random indegree if not given
        if indegrees is None:
            min_edges = len(targets) * 2
            max_edges = len(targets) * 4
            no_rules = random.randint(min_edges, max_edges)
        else:
            no_rules = np.random.choice(indegrees) * len(targets)
        rules_total += no_rules

        # the rules are not part of the original framework
        grn_vars["rules"] = generate_rules(no_rules, no_nodes, chance_1to1)

        motif = DiscreteGrnMotif(1, numvalues)
        motif.grn_vars = deepcopy(grn_vars)
        motif.construct_grn()

        motifs.append(motif)

    indegree_avg = (float(rules_total)/no_nodes)/samplesize
    return motifs, indegree_avg


def random_transition_table(no_nodes, num_values):
    """
    Generate a random transition table.

    PARAMETERS
    ---
    no_nodes: number of genes
    numvalues: number of possible values (int)

    RETURNS
    ---
    trans_table: a np array
    """

    # construct the leafcodes
    states = list(range(num_values))
    leafcodes = list(itertools.product(states, repeat=no_nodes))

    # make an empty state_transitions object
    state_transitions = []

    # fill the transitions table
    for _leafcode in leafcodes:
        for _ in range(0, no_nodes):
            random_state = np.random.randint(0, num_values)

            # update the leafcode
            _leafcode = list(_leafcode) + [random_state]

        # add to the state transition
        state_transitions.append(_leafcode)

    return np.array(state_transitions)


def LH_transition_table(no_nodes, num_values, seed):
    """
    Generate a transition table from the seed sampled using LHC scheme.

    PARAMETERS
    ---
    no_nodes: number of genes (int)
    numvalues: number of possible values (int)
    seed: random number between 0 and 1 (float)

    RETURNS
    ---
    trans_table: a np array
    """
    # calculate the maximum number of transition tables
    no_states = num_values ** (no_nodes)
    length_table = no_nodes * no_states
    no_possible_tables = num_values ** (length_table)

    # turn the seed into a number in the range 0 - max
    decimal_table = math.floor(seed * (no_possible_tables + 1))

    # turn this into a list representation of the digits in base "num_values"
    digits = []
    while decimal_table > 0:
        digits.insert(0, decimal_table % num_values)
        decimal_table  = decimal_table // num_values

    # pad with zeroes to the left until we have enough digits
    while len(digits) < length_table:
        digits = [0] + digits

    # construct the leafcodes
    states = list(range(num_values))
    leafcodes = list(itertools.product(states, repeat=no_nodes))

    # make an empty state_transitions object
    state_transitions = []
    # iterator for looping over the digits list
    iterator = 0

    # fill the transitions table
    for _leafcode in leafcodes:
        for _ in range(0, no_nodes):
            random_state = digits[iterator]
            iterator += 1

            # update the leafcode
            _leafcode = list(_leafcode) + [random_state]

        # add to the state transition
        state_transitions.append(_leafcode)

    return np.array(state_transitions)


def generate_random(samplesize, no_nodes, numvalues = 2):
    """
    This is a true random generator of transition table-based GRNs,
    and covers the entire search space, not just the bio-GRN part.


    PARAMETERS
    ---
    samplesize: number of objects in list
    no_nodes: number of genes
    numvalues: number of possible values (int)

    RETURNS
    ---
    motifs: list of DiscreteGrnMotif objects
    """
    # create empty list for the objects
    motifs = []

    # generate N samples
    for _ in range(0, samplesize):
        # create
        grn_vars = {}

        # set the size
        grn_vars["gene_cnt"] = no_nodes
        # set the correlation matrix
        grn_vars["correlations"] = generate_correlation_matrix(no_nodes)

        motif = DiscreteGrnMotif(1, numvalues)
        motif.grn_vars = deepcopy(grn_vars)
        motif.evaluation_style = 'transition_table'
        motif.transition_table = random_transition_table(no_nodes, numvalues)
        motif.construct_grn()

        motifs.append(motif)

    return motifs


def generate_random_LH(samplesize, no_nodes, numvalues = 2):
    """
    This is a true random generator of transition table-based GRNs,
    and covers the entire search space, not just the bio-GRN part.


    PARAMETERS
    ---
    samplesize: number of objects in list
    no_nodes: number of genes
    numvalues: number of possible values (int)

    RETURNS
    ---
    motifs: list of DiscreteGrnMotif objects
    """
    # create empty list for the objects
    motifs = []

    # generate N samples
    samples = lhs(n=1, samples=samplesize, criterion=None)
    for sample in samples:
        # create
        grn_vars = {}

        # set the size
        grn_vars["gene_cnt"] = no_nodes
        # set the correlation matrix
        grn_vars["correlations"] = generate_correlation_matrix(no_nodes)

        motif = DiscreteGrnMotif(1, numvalues)
        motif.grn_vars = deepcopy(grn_vars)
        motif.evaluation_style = 'transition_table'
        motif.transition_table = LH_transition_table(no_nodes, numvalues, sample)
        motif.construct_grn()

        motifs.append(motif)

    return motifs
