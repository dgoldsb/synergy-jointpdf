"""
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

Generates random discrete motif objects and saves them.
"""

from copy import deepcopy
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

    for _ in range(0, no_rules):
        # create a dict
        rule = {}

        # pick a random rule
        # find the number of inputs and outputs required
        choice = random.choice(functions_choices)

        # draw these values from the genes
        nodes = range(0, no_nodes)
        outputs = []
        for _ in range(0, choice["o"]):
            random.shuffle(nodes)
            outputs.append(nodes[0])
            nodes.pop()
        nodes = range(0, no_nodes)
        inputs = []
        for _ in range(0, choice["i"]):
            random.shuffle(nodes)
            inputs.append(nodes[0])
            nodes.pop()

        # set the values
        rule["inputs"] = inputs
        rule["outputs"] = outputs
        rule["rulefunction"] = choice["f"]

        # append to rules
        rules.append(rule)

    return rules

def generate_motifs(samplesize, no_nodes, indegree, numvalues=2):
    """
    Returns a list of objects.
    Improvable Barabasi-Albert network.

    PARAMETERS
    ---
    samplesize: number of objects in list
    no_nodes: number of genes
    indegree: average indegree desired (float)
    numvalues: number of possible values (int)

    RETURNS
    ---
    motifs: list of DiscreteGrnMotif objects
    """
    # make list for the objects
    motifs = []

    # generate N samples
    for _ in range(0, samplesize):
        # create
        grn_vars = {}

        # set the size
        grn_vars["gene_cnt"] = no_nodes
        grn_vars["conflict_rule"] = 'down'
        # set the correlation matrix
        grn_vars["correlations"] = generate_correlation_matrix(no_nodes)

        # the rules are not part of the original framework
        no_rules = int(indegree * no_nodes)
        grn_vars["rules"] = generate_rules_naive(no_rules, no_nodes)

        motif = DiscreteGrnMotif(1, numvalues, 'random')
        motif.grn_vars = deepcopy(grn_vars)
        motif.construct_grn()

        motifs.append(motif)

    return motifs
