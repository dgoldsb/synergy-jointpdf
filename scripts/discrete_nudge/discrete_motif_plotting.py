'''
This script allows for the creation of synergy profiles, the plotting of them,
and analysis on profiles.
'''

from __future__ import print_function

import itertools
import matplotlib.pyplot as plt
import numpy as np
import time

import discrete_motif_measures as measures
import discrete_motif_operations as operations

__author__ = 'dgoldsb'

def state_transition_table(motif, rule):
    """
    Find the state transition table, and return it as a list of list.
    This list of lists can be turned into a nice table using tabulate and display(HTML()).

    PARAMETERS:
    ---
    motif: a DiscreteGrnMotif object


    RETURNS:
    ---
    trans_table: a list of lists
    """
    # set up the table
    table = []
    header = []
    genes = list(range(0, motif.grn_vars["gene_cnt"]))
    for _gene in genes:
        header.append("Gene "+str(_gene)+" t=0")
    for _gene in genes:
        header.append("Gene "+str(_gene)+" t=1")
    table.append(header)

    # construct the leafcodes, each leafcode represents a state at T=0
    states = list(range(motif.numvalues))
    leafcodes = list(itertools.product(states, repeat=motif.grn_vars["gene_cnt"]))
    genes = list(range(0, motif.grn_vars["gene_cnt"]))

    # loop over all states
    for _leafcode in leafcodes:
        # we build a new leafcode, which is the new state
        # we build it as a list for table purposes
        leafcode_old = [str(e) for e in list(_leafcode)]
        leafcode_new = []

        for _gene in genes:
            # filter rules with this gene as the target
            transition_functions = []
            for _rule in motif.grn_vars["rules"]:
                # check if gene is the target, if so add
                if _gene in _rule["outputs"]:
                    transition_functions.append(_rule)

            # tally over all rules what state this should be
            # this is always deterministic
            tally = {}
            for i in range(2 * (-motif.numvalues), 2 * (motif.numvalues + 1)):
                tally[str(i)] = 0

            # loop over all relevant rules
            for _func in transition_functions:
                # prepare inputs
                inputs = []
                for index_input in _func["inputs"]:
                    inputs.append(_leafcode[index_input])

                outputs = []
                for index_output in _func["outputs"]:
                    outputs.append(_leafcode[index_output])

                # figure out the output state
                output_value_func = _func["rulefunction"](inputs)

                # add to the tally
                tally[str(int(output_value_func))] += 1

            # decide what it will be this leafcode
            #output_value = "/".join([str(i) for i, e in enumerate(tally) if e != 0])
            #if output_value == "":
            #    output_value = str(_leafcode[_gene])
            output_value = motif.decide_outcome(tally, rule, _leafcode[_gene])
            leafcode_new.append(output_value)
        row = []
        row.extend(leafcode_old)
        row.extend(leafcode_new)
        table.append(row)
    return table

def scatterplot_synergy_nudgeimpact(motifs, width, size, synergy_measure, filename=None, verbose=False):
    """
    Adds variables, enforcing set total correlations between them.

    PARAMETERS
    ---
    motifs: a list of DiscreteGrnMotif objects
    width: number of variables nudged (int)
    size: size of the nudge (float)
    synergy_measure: the synergy function to use (object)
    filename: name of pdf to save to (string)
    """
    impacts = []
    synergies = []
    for motif in motifs:
        if verbose:
            print("trying to find datapoint "+str(len(impacts)+1))

        # find the targets
        targets = []
        if motif.evaluation_style == 'network':
            for rule in motif.grn_vars["rules"]:
                targets = list(set(targets + rule["outputs"]))
            if verbose:
                print("the rules affect "+str(targets))

        # try to evaluate and find synergy
        try:
            # find the synergy
            # actually, let's use the full picture
            motif.evaluate_motif()
            time_start = time.time()
            mutual_information = measures.mutual_information(motif)
            synergy = synergy_measure(motif)
            time_end = time.time()
            time_diff = time_end - time_start
            if verbose:
                print("finding synergy took: "+str(time_diff))

            # find the nudge impact
            motif.reset_to_state(0)
            operations.nudge_variable(motif, width, size)
            motif.evaluate_motif()
            impact = measures.abs_diff(motif.states[-2], motif.states[-1])

            # normalize by dividing by the total mutual information
            synergies.append(float(synergy)/float(mutual_information))
            impacts.append(impact)
        except:
            print("failed to find the synergy")
    plt.scatter(impacts, synergies)
    plt.xlabel("Nudge impact")
    plt.ylabel("Synergy/MI")
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
