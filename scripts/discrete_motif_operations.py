'''
This file contains some of the common operations used on discrete motifs.
'''

from random import shuffle
import numpy as np

def select_subset(arr, threshold):
    '''
    Import from Derkjan that broke down, so manually implementing.

    Select a subset of the arr (randomly) which is at least bigger than
    threshold. From these parts of the array we subtract probability as part
    of our nudge.

    @param arr: a 1-d numpy array, all entries should be greater or equal to zero
    @param threshold: floating point number
    @returns: minus_states, a 1-d numpy filled with zeros and ones
    '''
    minus_states = np.random.randint(0, 2, arr.shape[0])
    attempt_count = 0
    while (np.sum(arr[minus_states != 0]) < threshold or np.all(minus_states == 1)) and attempt_count < 20:
        minus_states = np.random.randint(0, 2, arr.shape[0])
        attempt_count += 1

    if attempt_count > 19:
        print('could not find satisfying minus state randomly')
        minus_states = np.ones(arr.shape[0])
        indices = list(range(arr.shape[0]))
        np.random.shuffle(indices)
        for index in indices:
            if np.sum(arr[minus_states != 0]) - arr[index] > threshold:
                # was selected, which does not exist, I hope this is correct
                minus_states[index] = 0

    return minus_states


def mutate_array_bigger_zero(arr, nudge_size, option):
    '''
    Import from Derkjan that broke down, so manually implementing.

    Mutate the arr under the constraint that every entry must be
    bigger or equal to zero. The total plus and minus change should
    both be equal to nudge size.

    @param  arr: a 1-d np array
    @param nudge_size: A (small) number
    @param option: string in set: {proportional, random}
                   how to create the vector to perform the nudge. With proportional
                   the plus and minus states are chosen randomly (within the
                   constraint that the minus states must have weight bigger than
                   nudge size). Then the nudge size in the plus and minus states are
                   totally proportional to the old probability masses in those states.
                   For random the weights for the minus states are selected using the
                   Dirichlet distribution, if this pushes the state under zero than,
                   (that part of the nudge is redistributed among the other states
                   again using Dirichlet)
    @returns: nudged_array, a 1-d numpy array

    '''
    nudged_array = np.copy(arr)
    minus_states = select_subset(nudged_array, nudge_size)
    minus_mask = minus_states == 1
    plus_mask = minus_states == 0

    if np.sum(nudged_array[minus_mask]) < nudge_size:
        raise ValueError('chosen minus states wrongly')

    if option == 'proportional':
        minus_part = nudged_array[minus_mask]
        minus_nudge = (minus_part/np.sum(minus_part)) * nudge_size
        plus_part = nudged_array[plus_mask]
        plus_nudge = (plus_part/np.sum(plus_part)) * nudge_size
    elif option == 'random':
        minus_nudge = nudge_size * np.random.dirichlet(
            [1]*np.count_nonzero(minus_states)
        )
        minus_nudge = np.minimum(nudged_array[minus_mask], minus_nudge)
        difference = abs(np.sum(minus_nudge)-nudge_size)
        count = 0
        while difference > 10**(-10) and count < 10:
            count += 1
            number_of_free_states = minus_nudge[nudged_array[minus_mask] != minus_nudge].shape[0]
            redistribute = difference * np.random.dirichlet([1]*number_of_free_states)
            minus_nudge[nudged_array[minus_mask] != minus_nudge] += redistribute
            minus_nudge = np.minimum(nudged_array[minus_mask], minus_nudge)
            difference = abs(np.sum(minus_nudge)-nudge_size)
        if count == 10:
            print('could not find nudge totally randomly, now proportional')
            free_space = nudged_array[minus_mask] - minus_nudge
            minus_nudge += (free_space/np.sum(free_space))*difference

        if abs(np.sum(minus_nudge)-nudge_size) > 10**(-8):
            raise ValueError('minus nudge not big enough')
        minus_nudge = minus_nudge
        plus_nudge = nudge_size * np.random.dirichlet(
            [1]*(minus_states.shape[0]-np.count_nonzero(minus_states))
        )
    else:
        raise ValueError('wrong option parameter')

    nudged_array[minus_mask] -= minus_nudge
    nudged_array[plus_mask] += plus_nudge
    return nudged_array


def nudge_distribution_non_local_non_causal(joint, nudge_labels, nudge_size, nudge_option='random'):
    '''
    Import from Derkjan that broke down, so manually implementing.

    Nudge the the variable with nudge label while keeping the
    marginal of the other variables constant. Thus assuming that the variable
    on which the nudge is performed does not causally impact the other variables.
    The nudge moves weight around in a random manner. Meaning that the weight
    does not go from one state to another state, but rather a random
    perturbation vector is placed on the states, its sum being equal
    to 0 and its absolute sum equal to 2*nudge_size.

    @param joint: a numpy array representing a discrete probability distribution
    @param nudge_labels: a list of integers
    @param nudge_size: a (small) floating point number
    @param nudge_option: string in set {random, proportional}
                         see mutate_array_bigger_zero option docs
    @returns: nudged_joint, a nudged version of the joint
    '''
    # copy the joint distribution
    nudged_joint = np.copy(joint)

    # get a list of indices for the nudged variables
    nudged_indices = tuple(range(len(joint.shape) - len(nudge_labels), len(joint.shape), 1))
    # reorder the variables so that the last variables in the distribution are the nduged ones
    nudged_joint = np.moveaxis(nudged_joint, nudge_labels, nudged_indices)
    # calculate the nudge size applied per state, the nudge size is the fraction of the present probability moved
    # this has the dimension of the probability matrix of all unnudged variables
    nudged_size_per_state = np.sum(nudged_joint, axis=nudged_indices) * nudge_size

    # apply nudges
    it = np.nditer(nudged_size_per_state, flags=['multi_index'])
    while not it.finished:
        if it.value == 0:
            it.iternext()
            continue

        # each nudge, we pass a joint of the nudged variables
        # we nudge by the previously computed nudge size, which is based on a fraction of the probability present in this particular joint
        flattened_dist = nudged_joint[it.multi_index].flatten()
        nudged_state = np.reshape(
            mutate_array_bigger_zero(flattened_dist, it.value, nudge_option),
            nudged_joint[it.multi_index].shape
        )
        nudged_joint[it.multi_index] = nudged_state
        it.iternext()

    # put the variables back in the original order
    nudged_joint = np.moveaxis(nudged_joint, nudged_indices, nudge_labels)
    return nudged_joint


def nudge_variable(motif, targets, nudge_size, source):
    '''
    Nudge a variable or a number of variables, using DJ his method or the jointpdf method.

    @param motif: a motif object
    @param targets: the target variables to nudge
    @param nudge_size: the total size of the nudge
    @param source: the method to use, either jointpdf (jointpdf package) or DJ (modified implementation of
                   Derkjan his method, included here)
    '''

    if source == 'joint_pdf':
        motif.nudge(targets, [], nudge_size)
    elif source == 'DJ':
        joint = motif.joint_probabilities.joint_probabilities
        new_joint = nudge_distribution_non_local_non_causal(joint, targets, nudge_size)
        motif.joint_probabilities.joint_probabilities = new_joint
    else:
        raise ValueError('Source not correctly entered')
