# imports from the discrete motif package
from discrete_motif import DiscreteGrnMotif
import discrete_motif_functions as functions
import discrete_motif_measures as measures
import discrete_motif_generator as generator
import discrete_motif_operations as operations
import discrete_motif_plotting as visualize

# regular imports
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import random
import scipy.stats
import re
from sklearn.manifold import TSNE
import sys


def main():
    # XOR
    print("XOR")
    motif = generator.generate_random(samplesize=1, no_nodes=2, numvalues=2)[0]
    motif.transition_table = [[1, 1, 0, 0],[1, 0, 1, 1], [0, 1, 1, 1], [0, 0, 0, 0]]
    # experiment: set to equal chances
    motif.joint_probabilities.joint_probabilities = np.array([[0.25, 0.25], [0.25, 0.25]])
    motif.evaluate_motif()
    unnudged = motif.states[-1]
    print("The memory: ")
    motif.reset_to_state(0)
    print(measures.normalized_memory(motif))
    motif.reset_to_state(0)
    operations.nudge_variable(motif, [0], 0.25, 'DJ')
    motif.evaluate_motif()
    print("The nudge impact: ")
    print(measures.hellinger(unnudged, motif.states[-1]))
    print("The averaged nudge impact: ")
    print(measures.average_nudge_impact(motif, 1, 0.25, 'DJ'))
    print("Synergy")
    print(measures.normalized_synergy(motif, measures.synergy_middleground))

    # COPY
    print("\n\nCopy")
    motif = generator.generate_random(samplesize=1, no_nodes=2, numvalues=2)[0]
    motif.transition_table = [[1, 1, 1, 1],[1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 0, 0]]
    # experiment: set to equal chances
    motif.joint_probabilities.joint_probabilities = np.array([[0.25, 0.25], [0.25, 0.25]])
    motif.evaluate_motif()
    unnudged = motif.states[-1]
    print("The memory: ")
    motif.reset_to_state(0)
    print(measures.normalized_memory(motif))
    motif.reset_to_state(0)
    operations.nudge_variable(motif, [0], 0.25, 'DJ')
    motif.evaluate_motif()
    print("The nudge impact: ")
    print(measures.hellinger(unnudged, motif.states[-1]))
    print("The averaged nudge impact: ")
    print(measures.average_nudge_impact(motif, 1, 0.25, 'DJ'))
    print("Synergy")
    print(measures.normalized_synergy(motif, measures.synergy_middleground))

if __name__ == '__main__':
    main()