'''
This script allows for the creation of synergy profiles, the plotting of them,
and analysis on profiles.
'''

from __future__ import print_function

import itertools
import matplotlib.pyplot as plt
import numpy as np

import discrete_motif_measures as measures
import discrete_motif_operations as operations

__author__ = 'dgoldsb'

def scatterplot_synergy_nudgeimpact(motifs, width, size, wms=False, filename=None):
    """
    Adds variables, enforcing set total correlations between them.

    PARAMETERS
    ---
    motifs: a list of DiscreteGrnMotif objects
    width: number of variables nudged (int)
    size: size of the nudge (float)
    wms: use wms or Quax' synergy (Boolean)
    filename: name of pdf to save to (string)
    """
    impacts = []
    synergies = []
    for motif in motifs:
        print("trying to find datapoint "+str(len(impacts)+1))
        # we only evaluate the variables that are targeted
        targets = []
        for rule in motif.grn_vars["rules"]:
            targets = list(set(targets + rule["outputs"]))
        print("the rules affect "+str(targets))
        try:
            # find the synergy
            motif.evaluate_motif(targets)
            if wms:
                synergy = measures.synergy_wms(motif)
            else:
                synergy = measures.synergy_quax(motif)

            # find the nudge impact
            motif.reset_to_state(0)
            operations.nudge_variable(motif, width, size)
            motif.evaluate_motif()
            impact = measures.abs_diff(motif.states[-2], motif.states[-1])

            synergies.append(synergy)
            impacts.append(impact)
        except:
            print("failed to find the synergy")
    plt.scatter(impacts, synergies)
    plt.xlabel("Nudge impact")
    plt.ylabel("Synergy")
    if filename is not None:
        plt.savefig('myfig.pdf', format='pdf')
    plt.show()

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
