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
    plot_data = []

    # Determine the system size
    system_size = len(system.get_labels())

    # Relabel to be sure
    labels = []
    for i in range(0, system_size):
        labels.append(str(i))
    system.set_labels(labels)

    for i in range(1, system_size+1):
        combinations = itertools.combinations(labels, r=i)
        minfos = []
        for combination in combinations:
            combination = list(combination)
            minfo = system.mutual_information_labels(combination[0], combination.pop(0))
            print("Found MI of "+str(minfo)+" for combination "+str(combination))
            minfos.append(minfo)
        if mode == 'average':
            plot_data.append(np.average(minfos))
        elif mode == 'maximum':
            plot_data.append(np.max(minfos))

    return plot_data

def plot_mi_profile(plot_data):
    '''
    This function plots a profile, or a set of profiles
    '''
    if isinstance(plot_data[0], list):
        print('You supplied a list')
    else:
        plt.plot(plot_data)
        plt.show()
