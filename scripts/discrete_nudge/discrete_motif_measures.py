"""
This file contains some of the common measures used on discrete motifs.
"""

from __future__ import print_function

import copy
import math
import sys

import numpy as np

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
    print("Unimplemented...")
    sys.exit(0)

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

def synergy_wms(motif, num_synergistic_variables):
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

    # make a copy of the object
    motif_copy = copy.deepcopy(motif)

    # append SRVs
    motif_copy.append_synergistic_variables(num_synergistic_variables,
                                            subject_variables=genes_t0,
                                            agnostic_about=genes_t1)
    genes_srv = range(motif.numvariables, motif_copy.numvariables)

    # compare outputs with SRVs
    return motif_copy.synergistic_information_naive(genes_srv, genes_t1)
