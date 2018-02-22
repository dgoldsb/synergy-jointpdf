'''
This script allows for the creation of synergy profiles, the plotting of them,
and analysis on profiles.
'''

from __future__ import print_function

from copy import deepcopy
import itertools
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import random
import scipy.stats
import string
import time
from tqdm import tqdm

import discrete_motif_measures as measures
import discrete_motif_operations as operations

__author__ = 'dgoldsb'


def plot_bar(samples, colors, labels, title, filename=None, axes_labels=None):
    '''
    Make a sample of

    @param samples: list of lists of integers
    @param colors: the color of each bar
    @param labels: the label of each bar
    @param title: the plot title
    @param filename: the file to save to
    @param axes_labels: labels for the axes
    '''
    # to determine the range on the x-axis
    total_sample = []

    # for each, count how common it is
    bar_datas = []
    for sample in samples:
        bar_data = [0] * (max(sample) + 2)
        for i in range(0, max(sample) + 1):
            bar_data[i] = sample.count(i)
        bar_datas.append(bar_data)
        total_sample = total_sample + sample

    # preprocess the data to the right format
    index = list(range(0, max(total_sample) + 2))

    # settings
    opacity = 0.4
    bar_width = 0.35

    # make plot
    if title is not None:
        plt.title(title)
    if axes_labels is not None:
        plt.xlabel(axes_labels[0])
        plt.ylabel(axes_labels[1])
    iterator = 0
    for bar_data in bar_datas:
        bar_data = bar_data + [0] * (len(index) - len(bar_data))
        index_data = [x+(iterator * bar_width) for x in index]
        plt.bar(index_data, bar_data, bar_width, alpha=opacity, color=colors[iterator], label=labels[iterator])
        iterator += 1
    group_cnt = iterator

    index_data = [(x + ((group_cnt - 1) * bar_width)/group_cnt) for x in index]
    plt.xticks(index_data, index)
    plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=8)
    plt.tight_layout()
    if filename is not None:
        plt.savefig(filename, format='pdf')
        plt.close()
    else:
        plt.show()


def plot_line(values, colors, labels, title, filename=None, axes_labels=None):
    '''
    Make a line plot

    @param values: numpy array, should be a list containing [x-value, [y-values]]
    @param colors: the color of each point
    @param labels: the label of each point
    @param title: the plot title
    @param filename: the file to save to
    @param axes_labels: labels for the axes
    '''
    if title is not None:
        plt.title(title)
    if axes_labels is not None:
        plt.xlabel(axes_labels[0])
        plt.ylabel(axes_labels[1])
    for current_color in set(colors):
        x_values = []
        means = []
        confidence_intervals = []
        col_plot = []
        lab_plot = []
        for i in range(0, len(colors)):
            if colors[i] == current_color:
                # TODO: do a normality test, if it fails return without plotting
                x_values.append(values[i][0])
                n, min_max, mean, var, skew, kurt = scipy.stats.describe(values[i][1])
                std = math.sqrt(var)
                confidence_interval = scipy.stats.norm.interval(0.95, loc=mean, scale=std)
                means.append(mean)
                width_interval = confidence_interval[1] - mean
                confidence_intervals.append(width_interval)
                col_plot.append(colors[i])
                lab_plot.append(labels[i])
        plt.plot(x_values, means, c=col_plot[0], label=lab_plot[0])
        plt.errorbar(x_values, means, yerr=confidence_intervals, fmt='x')
    plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=8)
    if filename is not None:
        plt.savefig(filename, format='pdf')
        plt.close()
    else:
        plt.show()


def plot_scatter(x_values, colors, labels, title, filename=None, axes_labels=None):
    '''
    example: http://alexanderfabisch.github.io/t-sne-in-scikit-learn.html

    @param x_values: numpy array
    @param colors: the color of each point
    @param labels: the label of each point
    @param title: the plot title
    @param filename: the file to save to
    @param axes_labels: labels for the axes
    '''
    if title is not None:
        plt.title(title)
    if axes_labels is not None:
        plt.xlabel(axes_labels[0])
        plt.ylabel(axes_labels[1])
    for current_color in set(colors):
        x_plot = []
        col_plot = []
        lab_plot = []
        for i in range(0, len(colors)):
            if colors[i] == current_color:
                x_plot.append(x_values[i])
                col_plot.append(colors[i])
                lab_plot.append(labels[i])
        x_plot = np.array(x_plot)
        if current_color == 'red':
            mark = 'x'
        else:
            mark = 'o'
        plt.scatter(x_plot[:, 0], x_plot[:, 1], c=col_plot[0], label=lab_plot[0], marker=mark)
    plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=8)
    if filename is not None:
        plt.savefig(filename, format='pdf')
        plt.close()
    else:
        plt.show()


def plot_scatter_3d(x_values, colors, labels, title, filename=None, axes_labels=None):
    '''
    example: http://alexanderfabisch.github.io/t-sne-in-scikit-learn.html

    @param x_values: numpy array
    @param colors: the color of each point
    @param labels: the label of each point
    @param title: the plot title
    @param filename: the file to save to
    @param axes_labels: labels for the axes
    '''
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if title is not None:
        plt.title(title)
    if axes_labels is not None:
        ax.set_xlabel(axes_labels[0])
        ax.set_ylabel(axes_labels[1])
        ax.set_zlabel(axes_labels[2])
    for current_color in set(colors):
        x_plot = []
        col_plot = []
        lab_plot = []
        for i in range(0, len(colors)):
            if colors[i] == current_color:
                x_plot.append(x_values[i])
                col_plot.append(colors[i])
                lab_plot.append(labels[i])
        x_plot = np.array(x_plot)
        if current_color == 'red':
            mark = 'x'
        else:
            mark = 'o'
        ax.scatter(x_plot[:, 0], x_plot[:, 1], x_plot[:, 2], c=col_plot[0], label=lab_plot[0], marker=mark)
    plt.legend(loc='upper right', numpoints=1, ncol=3, fontsize=8)
    if filename is not None:
        plt.savefig(filename, format='pdf')
        plt.close()
    else:
        plt.show()


def state_transition_table(motif, rule):
    '''
    Find the state transition table, and return it as a list of list.
    This list of lists can be turned into a nice table using tabulate and display(HTML()).

    PARAMETERS:
    ---
    motif: a DiscreteGrnMotif object


    RETURNS:
    ---
    trans_table: a list of lists
    '''
    # set up the table
    table = []
    header = []
    genes = list(range(0, motif.grn_vars['gene_cnt']))
    for _gene in genes:
        header.append('Gene '+str(_gene)+' t=0')
    for _gene in genes:
        header.append('Gene '+str(_gene)+' t=1')
    table.append(header)

    # construct the leafcodes, each leafcode represents a state at T=0
    states = list(range(motif.numvalues))
    leafcodes = list(itertools.product(states, repeat=motif.grn_vars['gene_cnt']))
    genes = list(range(0, motif.grn_vars['gene_cnt']))

    # loop over all states
    for _leafcode in leafcodes:
        # we build a new leafcode, which is the new state
        # we build it as a list for table purposes
        leafcode_old = [str(e) for e in list(_leafcode)]
        leafcode_new = []

        for _gene in genes:
            # filter rules with this gene as the target
            transition_functions = []
            for _rule in motif.grn_vars['rules']:
                # check if gene is the target, if so add
                if _gene in _rule['outputs']:
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
                for index_input in _func['inputs']:
                    inputs.append(_leafcode[index_input])

                outputs = []
                for index_output in _func['outputs']:
                    outputs.append(_leafcode[index_output])

                # figure out the output state
                output_value_func = _func['rulefunction'](inputs)

                # add to the tally
                tally[str(int(output_value_func))] += 1

            output_value = motif.decide_outcome(tally, rule, _leafcode[_gene])
            leafcode_new.append(output_value)
        row = []
        row.extend(leafcode_old)
        row.extend(leafcode_new)
        table.append(row)
    return table


def generate_id(size=7, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def append_id(filename):
    name, ext = os.path.splitext(filename)
    return '{name}_{uid}{ext}'.format(name=name, uid=generate_id(), ext=ext)


def scatterplot_synergy_nudgeimpact(motifs, width, size, synergy_measure, colors=None, labels=None,
                                    filename=None, verbose=False):
    '''
    Adds variables, enforcing set total correlations between them.

    @param motifs: a list of DiscreteGrnMotif objects
    @param width: number of variables nudged (int)
    @param size: size of the nudge (float)
    @param synergy_measure: the synergy function to use (object)
    @param colors: the color of a scatterpoint
    @param labels: the label of a scatterpoint
    @param filename: name of pdf to save to (string)
    @returns: impacts, list of the impact of the nudge per motif
    @returns: synergies, list of the synergy per motif
    '''
    impacts = []
    synergies = []
    for motif in motifs:
        if verbose:
            print('trying to find datapoint '+str(len(impacts)+1))

        # find the targets
        targets = []
        if motif.evaluation_style == 'network':
            for rule in motif.grn_vars['rules']:
                targets = list(set(targets + rule['outputs']))
            if verbose:
                print('the rules affect '+str(targets))

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
                print('finding synergy took: '+str(time_diff))

            # find the nudge impact
            motif.reset_to_state(0)
            operations.nudge_variable(motif, width, size, 'DJ')
            motif.evaluate_motif()
            impact = measures.abs_diff(motif.states[-2], motif.states[-1])

            # normalize by dividing by the total mutual information
            synergies.append(float(synergy)/float(mutual_information))
            impacts.append(impact)
        except:
            print('failed to find the synergy')
            # to make sure the labels are still correct
            synergies.append(None)
            impacts.append(None)

    # create the plot
    if colors is None:
        colors = ['red'] * len(motifs)
    if labels is None:
        labels = [motif.evaluation_style] * len(motifs)
    if filename is not None:
        filename = append_id(filename)
    x_values = zip(impacts, synergies)
    axes_labels = ['Nudge impact', 'Synergy/MI']
    title = 'Relation between synergy and nudge impact (%s variables, %s-valued\
            logic, %s nudge affecting %s variables)' % (motifs[0].grn_vars['gene_cnt'],
                                                        motifs[0].numvalues, size, width)
    plot_scatter(x_values, colors, labels, title, filename, axes_labels)

    return impacts, synergies


def create_mi_profile(motif, mode):
    '''
    This function creates a MI profile from a system, which can be used later.

    @param motif: a motif object
    @param mode: the type of profile to create (string)
    '''
    # reset and timestep
    motif.evaluate_motif()

    if mode != 'maximum' and mode != 'average':
        print('Select the mode maximum or average.')
        return []

    # Create empty list to fill
    plot_data = [0]

    # Determine the system size
    system_size = motif.grn_vars['gene_cnt']

    # Relabel to be sure
    labels = range(0, system_size)
    motif.set_labels(labels)

    # Calculate the system entropy, for normalization
    mi_system = measures.mutual_information(motif)

    for i in range(1, system_size+1):
        combinations = itertools.combinations(labels, r=i)
        mis = []
        for combination in tqdm(combinations):
            combination = list(combination)
            
            mi = measures.mutual_information(motif, combination)
            # print('found entropy of '+str(entropy)+' for combination '+str(combination))
            mis.append(mi)
        if mode == 'average':
            plot_data.append(np.average(mis)/mi_system)
        elif mode == 'maximum':
            plot_data.append(np.max(mis)/mi_system)

    # reset the motif
    motif.reset_to_state(0)
    motif.states = [deepcopy(motif.joint_probabilities.joint_probabilities)]

    return plot_data


def plot_mi_profile(motifs, title=None, mode='maximum', filename=None):
    '''
    This function plots a profile, or a set of profiles

    @param motifs: a list of motif objects
    @param title: a string
    @param mode: the type of profile to create (string)
    @param filename: a string
    '''
    plot_data = create_mi_profile(motifs[0], mode)
    system_size = len(plot_data)-1
    x_ax = np.linspace(0, system_size, system_size+1)
    y_ax = x_ax/system_size
    for motif in tqdm(motifs):
        plot_data = create_mi_profile(motif, mode)
        if isinstance(plot_data[0], list):
            print('You supplied a list')
        else:
            plt.plot(x_ax, plot_data)
    axes_labels = ['MI (%s)' % mode, 'System size']
    plt.xlabel(axes_labels[0])
    plt.ylabel(axes_labels[1])
    if title is not None:
        title = 'Complexity profile of %d motifs' % len(motifs)
    plt.plot(x_ax, y_ax, ls='--')
    if title is not None:
        plt.title(title)
    if filename is not None:
        plt.savefig(filename, format='pdf')
        plt.close()
    else:
        plt.show()
