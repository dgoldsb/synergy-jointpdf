'''
This script allows for the creation of synergy profiles, the plotting of them,
and analysis on profiles.
'''

from __future__ import print_function

import sys
import itertools
import matplotlib.pyplot as plt
import numpy as np

__author__ = 'dgoldsb'

def scatterplot_synergy_nudgeimpact(motifs):
    """
    Adds variables, enforcing set total correlations between them.

    PARAMETERS
    ---
    motifs: a list of DiscreteGrnMotif objects
    """
    impacts = []
    synergies = []
    for motif in motifs:
        print("Unimplemented...")
    sys.exit(0)

def create_mi_profile(motif, mode):
    '''
    This function creates a MI profile from a system, which can be used later.
    '''
    if mode != 'maximum' and mode != 'average':
        print('Select the mode maximum or average.')
        return []

    # Create empty list to fill
    plot_data = [0]

    # Determine the system size
    system_size = motif.numvariables

    # Relabel to be sure
    labels = range(0, system_size)
    motif.set_labels(labels)

    # Calculate the system entropy, for normalization
    entropy_system = motif.entropy(labels)

    for i in range(1, system_size+1):
        combinations = itertools.combinations(labels, r=i)
        entropies = []
        for combination in combinations:
            combination = list(combination)
            entropy = motif.entropy(combination)
            print("found entropy of "+str(entropy)+" for combination "+str(combination))
            entropies.append(entropy)
        if mode == 'average':
            plot_data.append(np.average(entropies)/entropy_system)
        elif mode == 'maximum':
            plot_data.append(np.max(entropies)/entropy_system)

    return plot_data

def plot_mi_profile(motif, mode='maximum'):
    '''
    This function plots a profile, or a set of profiles
    '''
    plot_data = create_mi_profile(motif, mode)
    if isinstance(plot_data[0], list):
        print('You supplied a list')
    else:
        system_size = len(plot_data)-1
        x_ax = np.linspace(0, system_size, system_size+1)
        y_ax = x_ax/system_size
        plt.plot(x_ax, plot_data)
        plt.plot(x_ax, y_ax, ls="--")
        plt.show()
