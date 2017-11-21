"""
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

Generates random discrete motif objects and saves them.
"""

from copy import deepcopy
import itertools
import random

import numpy as np

from discrete_motif import DiscreteGrnMotif
import discrete_motif_functions as functions


def generate_correlation_matrix(gene_cnt):
    """
    Generate a completely random gene correlation matrix.

    PARAMETERS
    ---
    gene_cnt: number of genes (int)

    RETURNS
    ---
    matrix: a gene_cnt x gene_cnt matrix (ndarray of floats)
    """
    # make an empty correlation matrix
    matrix = np.empty([gene_cnt, gene_cnt], dtype=float)

    # fill the bands around the diagonal
    for i in range(1, gene_cnt):
        random_corr = np.random.uniform(low=-1.0, high=1.0)
        matrix[i, i-1] = random_corr
        matrix[i-1, i] = matrix[i, i-1]

    # fill the remaining values, compute them from the band
    # these are the indirect correlations, e.g. gene 0 and 2
    for i in range(2, gene_cnt):
        for j in range(0, i-1):
            matrix[i, j] = matrix[i, j+1] * matrix[j+1, j]
            matrix[j, i] = matrix[i, j]

    return matrix


def generate_rules(no_rules, no_nodes, chance_1to1):
    """
    We use a scale-free approach, with preferential attachment.
    We only use 1-to-1 and 2-to-1 functions, and randomly select origins
    and targets.

    The ratio 2-to-1 edges is way too high compared to the ratio of 1-to-1
    edges, the latter are much more common in nature.
    To compensate, we manually turn this to an approximately 25-75 ratio
    (configurable).


    PARAMETERS
    ---
    no_rules: the number of rules requested (int)
    no_nodes: the number of nodes requested (int)
    chance_1to1: the chance that an edge is 1-to-1 (float)

    RETURNS
    ---
    rules: a list of rule dictionaries
    """
    rules = []
    functions_choices = functions.dictionary()

    # create list of all possible target genes
    targets = range(0, no_nodes)

    # create list of all single actor genes
    actors = []
    actors_tuples = []
    for i in range(1, 3):
        actors_tuples += itertools.combinations(targets, i)
        actors = [list(t) for t in actors_tuples]

    # create a list of all possible in/out combinations
    edges = []
    for i in range(0, len(actors)):
        for j in range(0, len(targets)):
            combination = []
            combination.append(deepcopy(actors[i]))
            combination.append(deepcopy(targets[j]))
            edges.append(combination)

    # we use and adaptation of the Barabasi-Albert algorithm
    # we slightly adapt the algorithm to function with a set number of relations
    # and to utilize both 1-to-1 directed edges as 2-to-1 directed edges
    genes_referenced_cnt = [0] * no_nodes
    nodes_inserted = set([])
    for i in range(0, no_nodes):
        rules_left = no_rules - len(rules)
        nodes_left = no_nodes - len(nodes_inserted)

        # add the node
        nodes_inserted.add(i)
        
        # start adding rules
        for _ in range(0, (rules_left // nodes_left)):
            # add a rule
            rule = {}

            # we shuffle the edges and pick the first that contains only added
            # genes
            random.shuffle(edges)

            # we draw a number to  see the edge should be 1-to-1
            if i == 0:
                random_number_1to1 = 1
            else:
                random_number_1to1 = np.random.random()

            for k in range(0, len(edges)):
                # join the sets of inputs and outputs of the edge
                joint_set = set(edges[k][0]).union(set([edges[k][1]]))

                # see if we want a 1-to-1 edge
                proceed = False
                if random_number_1to1 < chance_1to1 and len(edges[k][0]) == 1:
                    proceed = True
                elif random_number_1to1 >= chance_1to1 and len(edges[k][0]) > 1:
                    proceed = True
                else:
                    proceed = False
                
                # only proceed if all inserted nodes are in the edge, and if the
                # currently added node is in the set
                if (len(joint_set - nodes_inserted) == 0) and (len(set([i]) - joint_set) == 0) and proceed:
                    # calculate the total number of edges
                    # divide the mean of the number of edges that are not the
                    # currently added edge by this number
                    relevant_edge_cnt = 0
                    for j in range(0, i+1):
                        if len(set([j]) - joint_set) == 0:
                            relevant_edge_cnt += genes_referenced_cnt[j]

                    relevant_edge_cnt = relevant_edge_cnt / (i+1)
                    if sum(genes_referenced_cnt) == 0:
                        accept_probability = 1
                    else:
                        accept_probability = relevant_edge_cnt / sum(genes_referenced_cnt)
                        
                    # the random part
                    random_number = np.random.random()
                    if random_number < accept_probability:
                        edge = deepcopy(edges[k])
                        edges.pop(k)

                        # find the no_actors
                        no_actors = len(edge[0])

                        # pick a random rule
                        # find the number of inputs and outputs required
                        choice = random.choice(functions_choices[no_actors])

                        # set the values
                        rule["inputs"] = deepcopy(edge[0])
                        rule["outputs"] = [deepcopy(edge[1])]
                        rule["rulefunction"] = choice["f"]

                        # append to rules
                        rules.append(rule)

                        # add one to all setmembers
                        for m in joint_set:
                            genes_referenced_cnt[m] += 1

                        # we found a rule, so break
                        break
    return rules


def generate_motifs(samplesize, no_nodes, numvalues=2, indegrees=None, conflict_rule='totaleffect', chance_1to1=0.75):
    """
    Returns a list of objects.
    Improvable Barabasi-Albert network.
    This returns networks, which represent a part of the transition table search
    space that we think contains all biological GRN networks.

    PARAMETERS
    ---
    samplesize: number of objects in list
    no_nodes: number of genes
    indegrees: list of possible indegree desired (average indegree will tend
    towards the average of the list), if None random (list of integers)
    numvalues: number of possible values (int)
    conflict_rule: the rule for deciding what state we get when rules conflict (str)

    RETURNS
    ---
    motifs: list of DiscreteGrnMotif objects
    """
    # create empty list for the objects
    motifs = []

    # create list of all possible target genes
    targets = range(0, no_nodes)

    # create list of all single actor genes
    actors = []
    actors_tuples = []
    for i in range(1, 3):
        actors_tuples += itertools.combinations(targets, i)
        actors = [list(t) for t in actors_tuples]

    # for calculating the average indegree
    rules_total = 0

    # generate N samples
    for _ in range(0, samplesize):
        # create
        grn_vars = {}

        # set the size
        grn_vars["gene_cnt"] = no_nodes
        grn_vars["conflict_rule"] = conflict_rule
        # set the correlation matrix
        grn_vars["correlations"] = generate_correlation_matrix(no_nodes)

        # generate a random indegree if not given
        if indegrees is None:
            max_edges = len(targets) * len(actors)
            no_rules = random.randint(0, max_edges)
        else:
            no_rules = np.random.choice(indegrees)
        rules_total += no_rules

        # the rules are not part of the original framework
        grn_vars["rules"] = generate_rules(no_rules, no_nodes, chance_1to1)

        motif = DiscreteGrnMotif(1, numvalues)
        motif.grn_vars = deepcopy(grn_vars)
        motif.construct_grn()

        motifs.append(motif)

    indegree_avg = (float(rules_total)/no_nodes)/samplesize
    return motifs, indegree_avg


def random_transition_table(no_nodes, num_values):
    """
    Generate a random transition table.

    PARAMETERS
    ---
    no_nodes: number of genes
    numvalues: number of possible values (int)

    RETURNS
    ---
    trans_table: a numpy array
    """

    # construct the leafcodes
    states = list(range(num_values))
    leafcodes = list(itertools.product(states, repeat=no_nodes))

    # make an empty state_transitions object
    state_transitions = []

    # fill the transitions table
    for _leafcode in leafcodes:
        for _ in range(0, no_nodes):
            random_state = np.random.randint(0, num_values)

            # update the leafcode
            _leafcode = list(_leafcode) + [random_state]

        # add to the state transition
        state_transitions.append(_leafcode)

    return np.array(state_transitions)


def generate_random(samplesize, no_nodes, numvalues = 2):
    """
    This is a true random generator of transition table-based GRNs,
    and covers the entire search space, not just the bio-GRN part.


    PARAMETERS
    ---
    samplesize: number of objects in list
    no_nodes: number of genes
    numvalues: number of possible values (int)

    RETURNS
    ---
    motifs: list of DiscreteGrnMotif objects
    """
    # create empty list for the objects
    motifs = []

    # generate N samples
    for _ in range(0, samplesize):
        # create
        grn_vars = {}

        # set the size
        grn_vars["gene_cnt"] = no_nodes
        # set the correlation matrix
        grn_vars["correlations"] = generate_correlation_matrix(no_nodes)

        motif = DiscreteGrnMotif(1, numvalues)
        motif.grn_vars = deepcopy(grn_vars)
        motif.evaluation_style = 'transition_table'
        motif.transition_table = random_transition_table(no_nodes, numvalues)
        motif.construct_grn()

        motifs.append(motif)

    return motifs
