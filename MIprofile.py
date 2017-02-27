'''
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This script allows for the creation of synergy profiles, the plotting of them,
and analysis on sets of profiles.
'''

from __future__ import print_function
import jointpdf

__author__ = 'dgoldsb'

def create_random_systemset(numprofiles=10, numvariables=2, numvalues=2, joint_probs='unbiased'):
    '''
    This function creates a set of random systems
    '''
    systemset = []
    for i in range(0, numprofiles):
        system = jointpdf.JointProbabilityMatrix(numvariables=numvariables, numvalues=numvalues
                                                 , joint_probs=joint_probs)
        systemset.append({'ID': i, 'Instance':system})
    return systemset

def create_MI_profile(system, mode='maximum'):
    '''
    This function creates a MI profile from a system, which can be used later.
    '''
    if mode == 'maximum':
        profile = 1
    elif mode == 'average':
        profile = 0
    else:
        print('Select the mode maximum or average.')
    return profile

def plot_MI_profile(profileset):
    '''
    This function plots a profile, or a set of profiles
    '''
    return 0