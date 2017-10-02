"""
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

Generates random discrete motif objects and saves them.
"""

from copy import deepcopy
import itertools
import random

import numpy as np

from discrete_motif import DiscreteGrnMotif
import discrete_motif_functions as functions

def generate_correlation_matrix(gene_cnt):
    """
    Generate a completely random gene correlation matrix.

    PARAMETERS
    ---
    gene_cnt: number of genes (int)

    RETURNS
    ---
    matrix: a gene_cnt x gene_cnt matrix (ndarray of floats)
    """
    # make an empty correlation matrix
    matrix = np.empty([gene_cnt, gene_cnt], dtype=float)

    # fill the bands around the diagonal
    for i in range(1, gene_cnt):
        random_corr = np.random.uniform(low=-1.0, high=1.0)
        matrix[i, i-1] = random_corr
        matrix[i-1, i] = matrix[i, i-1]

    # fill the remaining values, compute them from the band
    # these are the indirect correlations, e.g. gene 0 and 2
    for i in range(2, gene_cnt):
        for j in range(0, i-1):
            matrix[i, j] = matrix[i, j+1] * matrix[j+1, j]
            matrix[j, i] = matrix[i, j]

    return matrix

def generate_rules_naive(no_rules, no_nodes):
    """
    We should look into the way the networks are typically build up.
    However, for now a naive approach will do.
    We only use 1-to-1 and 2-to-1 functions, and randomly select origins
    and targets.

    PARAMETERS
    ---
    no_rules: the number of rules requested (int)

    RETURNS
    ---
    rules: a list of rule dictionaries
    """
    rules = []
    functions_choices = functions.dictionary()

    # create list of all possible target genes
    targets = range(0, no_nodes)

    # create list of all single actor genes
    actors = []
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

    for _ in range(0, no_rules):
        # create a dict
        rule = {}

        # pick a random edge to find a rule for
        random.shuffle(edges)
        edge = deepcopy(edges[0])
        edges.pop(0)

        # find the no_actors
        no_actors = len(edge[0])

        # pick a random rule
        # find the number of inputs and outputs required
        choice = random.choice(functions_choices[no_actors])

        # set the values
        rule["inputs"] = deepcopy(edge[0])
        rule["outputs"] = [deepcopy(edge[1])]
        rule["rulefunction"] = choice["f"]

        # append to rules
        rules.append(rule)

    return rules

def generate_motifs(samplesize, no_nodes, indegree=None, numvalues=2, conflict_rule = 'odds'):
    """
    Returns a list of objects.
    Improvable Barabasi-Albert network.

    PARAMETERS
    ---
    samplesize: number of objects in list
    no_nodes: number of genes
    indegree: average indegree desired, if None random (float)
    numvalues: number of possible values (int)

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
        if indegree is None:
            max_edges = len(targets) * len(actors)
            no_rules = random.randint(0, max_edges)
        else:
            no_rules = int(indegree * no_nodes)
        rules_total += no_rules

        # the rules are not part of the original framework
        grn_vars["rules"] = generate_rules_naive(no_rules, no_nodes)

        motif = DiscreteGrnMotif(1, numvalues, 'random')
        motif.grn_vars = deepcopy(grn_vars)
        motif.construct_grn()

        motifs.append(motif)

    indegree_avg = (float(rules_total)/no_nodes)/samplesize
    return motifs, indegree_avg
