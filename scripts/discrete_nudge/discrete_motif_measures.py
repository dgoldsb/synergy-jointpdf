"""
This file contains some of the common measures used on discrete motifs.
"""

import math

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

def mutual_information():
    return 0

def synergy_quax():
    return 0

def synergy_wms():
    return 0