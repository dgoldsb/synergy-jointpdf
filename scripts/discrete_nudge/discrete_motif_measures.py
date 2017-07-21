"""
This file contains some of the common measures used on discrete motifs.
"""

import numpy as np

from scipy.stats import entropy

def KL_div(tree_1, tree_2):
    """
    Finds the KL-div between two same-shaped FullNestedArrayOfProbabilities objects.

    PARAMETERS
    ---
    tree_1: FullNestedArrayOfProbabilities object
    tree_2: FullNestedArrayOfProbabilities object

    RETURNS
    ---
    Kullback-Leibler divergence: float
    """
    # flatten the trees
    tree_1_flat = np.copy(tree_1).flatten()
    tree_2_flat = np.copy(tree_2).flatten()

    return entropy(tree_1_flat, tree_2_flat)

def mutual_information():
    return 0

def synergy_quax():
    return 0

def synergy_wms():
    return 0