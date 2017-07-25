'''
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This script allows for the creation of synergy profiles, the plotting of them,
and analysis on sets of profiles.
'''

from __future__ import print_function
import itertools
import matplotlib.pyplot as plt
from jointpdf import jointpdf
import numpy as np

__author__ = 'dgoldsb'

def create_system_discrete(numvariables=5, numvalues=2, joint_probs='unbiased'):
    '''
    This function creates a set of random systems
    '''
    system = jointpdf.JointProbabilityMatrix(numvariables=numvariables, numvalues=numvalues
                                             , joint_probs=joint_probs)
    return system

def create_systemset_discrete(numprofiles=10, numvariables=2, numvalues=2, joint_probs='unbiased'):
    '''
    This function creates a set of random systems
    '''
    systemset = []
    for i in range(0, numprofiles):
        system = jointpdf.JointProbabilityMatrix(numvariables=numvariables, numvalues=numvalues
                                                 , joint_probs=joint_probs)
        systemset.append({'ID': i, 'Instance':system})
    return systemset

def create_mi_profile(system, mode='maximum'):
    '''
    This function creates a MI profile from a system, which can be used later.
    '''
    if mode != 'maximum' and mode != 'average':
        print('Select the mode maximum or average.')
        return []

    # Create empty list to fill
    plot_data = [0]

    # Determine the system size
    system_size = len(system.get_labels())

    # Relabel to be sure
    labels = []
    for i in range(0, system_size):
        labels.append(i)
    system.set_labels(labels)

    # Calculate the system entropy, for normalization
    entropy_system = system.entropy(labels)

    for i in range(1, system_size+1):
        combinations = itertools.combinations(labels, r=i)
        entropies = []
        for combination in combinations:
            combination = list(combination)
            entropy = system.entropy(combination)
            print("Found entropy of "+str(entropy)+" for combination "+str(combination))
            entropies.append(entropy)
        if mode == 'average':
            plot_data.append(np.average(entropies)/entropy_system)
        elif mode == 'maximum':
            plot_data.append(np.max(entropies)/entropy_system)

    return plot_data

def plot_mi_profile(plot_data):
    '''
    This function plots a profile, or a set of profiles
    '''
    if isinstance(plot_data[0], list):
        print('You supplied a list')
    else:
        system_size = len(plot_data)-1
        x_ax = np.linspace(0, system_size, system_size+1)
        y_ax = x_ax/system_size
        plt.plot(x_ax, plot_data)
        plt.plot(x_ax, y_ax, ls="--")
        plt.show()
