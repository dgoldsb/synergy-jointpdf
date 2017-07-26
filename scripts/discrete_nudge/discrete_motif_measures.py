"""
This file contains some of the common measures used on discrete motifs.
"""

from __future__ import print_function

import copy
import math
import sys

import numpy as np
from scipy.linalg import norm
from scipy.stats import entropy

def abs_diff(tree_1, tree_2):
    """
    Finds the absolute difference between two same-shaped FullNestedArrayOfProbabilities objects.
    We chose this, as the Kullback-Leibler divergence cannot handle zeros
    PARAMETERS
    ---
    tree_1: FullNestedArrayOfProbabilities object
    tree_2: FullNestedArrayOfProbabilities object

    RETURNS
    ---
    absolute difference: float
    """
    # flatten the trees
    t1_flat = np.array(np.copy(tree_1).flatten())
    t2_flat = np.array(np.copy(tree_2).flatten())

    returnval = 0
    for i in range(0, len(t1_flat)):
        returnval = math.fabs(t1_flat[i] - t2_flat[i])

    return returnval

def hellinger(tree_1, tree_2):
    """
    Finds the absolute difference between two same-shaped FullNestedArrayOfProbabilities objects.
    We chose this, as the Kullback-Leibler divergence cannot handle zeros
    PARAMETERS
    ---
    tree_1: FullNestedArrayOfProbabilities object
    tree_2: FullNestedArrayOfProbabilities object

    RETURNS
    ---
    absolute difference: float
    """
    # flatten the trees
    t1_flat = np.array(np.copy(tree_1).flatten())
    t2_flat = np.array(np.copy(tree_2).flatten())

    return norm(np.sqrt(t1_flat) - np.sqrt(t2_flat)) / np.sqrt(2)

def kl_div(tree_1, tree_2):
    """
    Finds the KL-divergence for the motif after a time evaluation.

    PARAMETERS
    ---
    motif: a DiscreteGrnMotif object

    RETURNS
    ---
    returnval: KL divergence (float)
    """
    # flatten the trees
    t1_flat = np.array(np.copy(tree_1).flatten())
    t2_flat = np.array(np.copy(tree_2).flatten())

    return entropy(t1_flat, t2_flat)

def mutual_information(motif):
    """
    Finds the mutual information for the motif after a time evaluation.

    PARAMETERS
    ---
    motif: a DiscreteGrnMotif object

    RETURNS
    ---
    returnval: mutual information (float)
    """
    if motif.grn_vars["gene_cnt"] == motif.numvariables:
        raise ValueError("do a time evaluation before attemting to check MI")

    # define what t0 and t1 are
    genes_t0 = range(0, motif.grn_vars["gene_cnt"])
    genes_t1 = range(motif.grn_vars["gene_cnt"], motif.numvariables)

    return motif.mutual_information(genes_t0, genes_t1)

def synergy_quax(motif):
    """
    Finds the mutual information for the motif after a time evaluation.

    PARAMETERS
    ---
    motif: a DiscreteGrnMotif object

    RETURNS
    ---
    returnval: upper bound of synergy according to Rick Quax (float)
    """
    if motif.grn_vars["gene_cnt"] == motif.numvariables:
        raise ValueError("do a time evaluation before attemting to check synergy")

    # define what t0 and t1 are
    genes_t0 = range(0, motif.grn_vars["gene_cnt"])
    genes_t1 = range(motif.grn_vars["gene_cnt"], motif.numvariables)

    return motif.synergistic_information(genes_t0, genes_t1)

def synergy_wms(motif):
    """
    Finds the mutual information for the motif after a time evaluation.

    PARAMETERS
    ---
    motif: a DiscreteGrnMotif object

    RETURNS
    ---
    returnval: WMS synergy (float)
    """
    if motif.grn_vars["gene_cnt"] == motif.numvariables:
        raise ValueError("do a time evaluation before attemting to check synergy")

    # define what t0 and t1 are
    genes_t0 = range(0, motif.grn_vars["gene_cnt"])
    genes_t1 = range(motif.grn_vars["gene_cnt"], motif.numvariables)

    return motif.synergistic_information_naive(genes_t0, genes_t1)

#TODO: something with SRVs?
