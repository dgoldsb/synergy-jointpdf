__author__ = 'rquax'


'''
Owner of this project:

    Rick Quax
    https://staff.fnwi.uva.nl/r.quax/index.html
    University of Amsterdam

You are free to use this package for your own non-profit, academic work. All I ask is to be suitably credited.
'''


import csv
import numpy as np
import itertools
import copy
from scipy.optimize import minimize
import warnings
from compiler.ast import flatten  # note: deprecated in Python 3, in which case: find another flatten
from collections import Sequence
import time


def maximum_depth(seq):
    """
    Helper function, e.g. maximum_depth([1,2,[2,4,[[4]]]]) == 4.
    :param seq: sequence, like a list of lists
    :rtype: int
    """
    seq = iter(seq)
    try:
        for level in itertools.count():
            seq = itertools.chain([next(seq)], seq)
            seq = itertools.chain.from_iterable(s for s in seq if isinstance(s, Sequence))
    except StopIteration:
        return level


# helper function,
# from http://stackoverflow.com/questions/2267362/convert-integer-to-a-string-in-a-given-numeric-base-in-python
def int2base(x, b, alphabet='0123456789abcdefghijklmnopqrstuvwxyz'):
    """

    :param x: int
    :type x: int
    :param b: int
    :param b: int
    :param alphabet:
    :rtype : str
    """

    # convert an integer to its string representation in a given base
    if b<2 or b>len(alphabet):
        if b==64: # assume base64 rather than raise error
            alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
        else:
            raise AssertionError("int2base base out of range")

    if isinstance(x,complex): # return a tuple
        return ( int2base(x.real,b,alphabet) , int2base(x.imag,b,alphabet) )

    if x<=0:
        if x==0:
            return alphabet[0]
        else:
            return '-' + int2base(-x,b,alphabet)

    # else x is non-negative real
    rets=''

    while x>0:
        x,idx = divmod(x,b)
        rets = alphabet[idx] + rets

    return str(rets)


def apply_permutation(lst, permutation):
    """
    Return a new list where the element at position ix in <lst> will be at a new position permutation[ix].
    :param lst: list
    :type lst: array_like
    :param permutation:
    :return:
    """
    new_list = [-1]*len(lst)

    assert len(permutation) == len(lst)

    for ix in xrange(len(permutation)):
        new_list[permutation[ix]] = lst[ix]

    if __debug__:
        if not -1 in lst:
            assert not -1 in new_list

    return new_list


# each variable in a JointProbabilityMatrix has a label, and if not provided then this label is used
_default_variable_label = 'variable'


class JointProbabilityMatrix():
    # nested list of probabilities. For instance for three binary variables it could be:
    # [[[0.15999394, 0.06049343], [0.1013956, 0.15473886]], [[ 0.1945649, 0.15122334], [0.11951818, 0.05807175]]].
    # The depth (np.ndim) should be numvariables, and the length at each dimension should be numvalues.
    joint_probabilities = []

    numvariables = 0
    numvalues = 0

    # note: at the moment I only store the labels here as courtesy, but the rest of the functions do not operate
    # on labels at all, only variable indices 0..N-1.
    labels = []

    def __init__(self, numvariables, numvalues, joint_probs=None, labels=None):
        self.numvariables = numvariables
        self.numvalues = numvalues

        if labels is None:
            self.labels = [_default_variable_label]*numvariables
        else:
            self.labels = labels

        if joint_probs is None:
            self.generate_random_joint_probabilities()
        elif isinstance(joint_probs, basestring):
            if joint_probs == 'random' or joint_probs == 'independent' or joint_probs == 'iid':
                self.generate_independent_joint_probabilities()
            else:
                raise ValueError('don\'t know what to do with joint_probs=' + str(joint_probs))
        else:
            self.joint_probabilities = joint_probs

        self.clip_all_probabilities()

        if self.numvariables > 0:
            assert np.ndim(self.joint_probabilities) == self.numvariables
            if self.numvariables > 1:
                assert len(self.joint_probabilities[0]) == numvalues, 'self.joint_probabilities[0] = ' \
                                                                      + str(self.joint_probabilities[0])
                assert len(self.joint_probabilities[-1]) == numvalues
            else:
                assert len(self.joint_probabilities) == numvalues
                assert np.isscalar(self.joint_probabilities[(0,)]), 'this is how joint_probabilities may be called'

            np.testing.assert_almost_equal(np.sum(self.joint_probabilities), 1.0)
        else:
            raise ValueError('numvariables == 0 not supported (yet?)')


    def copy(self):
        """
        Deep copy.
        :rtype : JointProbabilityMatrix
        """
        return copy.deepcopy(self)


    def generate_sample(self):
        rand_real = np.random.random()

        for input_states in self.statespace(self.numvariables):
            rand_real -= self.joint_probability(input_states)

            if rand_real < 0.0:
                return input_states

        assert False, 'should not happen? or maybe foating point roundoff errors built up.. or so'

        return self.statespace(len(self.numvariables))[-1]


    def generate_samples(self, n):
        return [self.generate_sample() for _ in xrange(n)]


    def transition_table_deterministic(self, variables_input, variables_output):
        """
        Create a transition table for the given variables, where for each possible sequence of values for the input
        variables I find the values for the output variables which have maximum conditional probability. I marginalize
        first over all other variables not named in either input or output.
        :type variables_input: list of int
        :type variables_output: list of int
        :return: a list of lists, where each sublist is [input_states, max_output_states, max_output_prob], where the
        first two are lists and the latter is a float in [0,1]
        """
        variables = list(variables_input) + list(variables_output)

        pdf_temp = self.marginalize_distribution(variables)

        assert len(pdf_temp) == len(variables)

        trans_table = []

        # for each possible sequence of input states, find the output states which have maximum probability
        for input_states in self.statespace(len(variables_input)):
            max_output_states = []
            max_output_prob = -1.0

            for output_states in self.statespace(len(variables_output)):
                prob = pdf_temp(list(input_states) + list(output_states))

                if prob > max_output_prob:
                    max_output_prob = prob
                    max_output_states = output_states

            # assert max_output_prob >= 1.0 / np.power(self.numvalues, len(variables_output)), \
            #     'how can this be max prob?: ' + str(max_output_prob) + '. I expected at least: ' \
            #     + str(1.0 / np.power(self.numvalues, len(variables_output)))

            trans_table.append([input_states, max_output_states, max_output_prob])

        return trans_table


    def generate_independent_joint_probabilities(self):
        self.joint_probabilities = np.zeros([self.numvalues]*self.numvariables)
        self.joint_probabilities = self.joint_probabilities + 1.0 / np.power(self.numvalues, self.numvariables)
        np.testing.assert_almost_equal(np.sum(self.joint_probabilities), 1.0)


    def generate_random_joint_probabilities(self):
        self.joint_probabilities = np.random.random([self.numvalues]*self.numvariables)
        self.joint_probabilities /= np.sum(self.joint_probabilities)


    def random_samples(self, n=1):
        # sample_indices = np.random.multinomial(n, [i for i in self.joint_probabilities.flat])

        # todo: I don't think this generalizes well for different numvalues per variable. Just generate 1000 value lists
        # and pick one according to their relative probabilities? Much more general and above all much more scalable!

        assert np.isscalar(self.numvalues), 'refactor this function, see todo above'

        flat_joint_probs = [i for i in self.joint_probabilities.flat]

        sample_indices = np.random.choice(len(flat_joint_probs), p=flat_joint_probs, size=n)

        values_str_per_sample = [int2base(smpl, self.numvalues).zfill(self.numvariables) for smpl in sample_indices]

        assert len(values_str_per_sample) == n
        assert len(values_str_per_sample[0]) == self.numvariables, 'values_str_per_sample[0] = ' \
                                                                   + str(values_str_per_sample[0])
        assert len(values_str_per_sample[-1]) == self.numvariables, 'values_str_per_sample[-1] = ' \
                                                                   + str(values_str_per_sample[-1])

        values_list_per_sample = [[int(val) for val in valstr] for valstr in values_str_per_sample]

        assert len(values_list_per_sample) == n
        assert len(values_list_per_sample[0]) == self.numvariables
        assert len(values_list_per_sample[-1]) == self.numvariables

        return values_list_per_sample


    def __call__(self, values):
        """
        Joint probability of a list or tuple of values, one value for each variable in order.
        :param values: list of values, each value is an integer in [0, numvalues)
        :type values: list
        :return: joint probability of the given values, in order, for all variables; in [0, 1]
        :rtype: float
        """
        return self.joint_probability(values=values)


    def set_labels(self, labels):
        self.labels = labels


    def get_labels(self):
        return self.labels


    def get_label(self, variable_index):
        return self.labels[variable_index]


    def set_label(self, variable_index, label):
        self.labels[variable_index] = label


    def has_label(self, label):
        return label in self.labels


    def retain_only_labels(self, retained_labels, dict_of_unwanted_label_values):
        unwanted_labels = list(set(self.labels).difference(retained_labels))

        # unwanted_label_values = [dict_of_unwanted_label_values[lbl] for lbl in unwanted_labels]

        unwanted_vars = self.get_variables_with_label_in(unwanted_labels)
        unwanted_var_labels = map(self.get_label, unwanted_vars)
        unwanted_var_values = [dict_of_unwanted_label_values[lbl] for lbl in unwanted_var_labels]

        self.duplicate(self.conditional_probability_distribution(unwanted_vars, unwanted_var_values))



    def clip_all_probabilities(self):
        """
        Make sure all probabilities in the joint probability matrix are in the range [0.0, 1.0], which could be
        violated sometimes due to floating point operation roundoff errors.
        """
        self.joint_probabilities = np.minimum(np.maximum(self.joint_probabilities, 0.0), 1.0)

        if np.random.random() < 0.1:
            try:
                np.testing.assert_almost_equal(np.sum(self.joint_probabilities), 1.0)
            except AssertionError as e:
                print 'error message: ' + str(e)

                print 'error: len(self.joint_probabilities) =', len(self.joint_probabilities)
                print 'error: shape(self.joint_probabilities) =', np.shape(self.joint_probabilities)
                if len(self.joint_probabilities) < 30:
                    print 'error: self.joint_probabilities =', self.joint_probabilities

                raise AssertionError(e)


    def estimate_from_data_from_file(self, filename, discrete=True, method='empirical'):
        """
        Same as estimate_from_data but read the data first from the file, then pass data to estimate_from_data.
        :param filename: string, file format will be inferred from the presence of e.g. '.csv'
        :param discrete:
        :param method:
        :raise NotImplementedError: if file format is not recognized
        """
        if '.csv' in filename:
            fin = open(filename, 'rb')
            csvr = csv.reader(fin)

            repeated_measurements = [tuple(row) for row in csvr]

            fin.close()
        else:
            raise NotImplementedError('unknown file format: ' + str(filename))

        self.estimate_from_data(repeated_measurements=repeated_measurements, discrete=discrete, method=method)


    def estimate_from_data(self, repeated_measurements, discrete=True, method='empirical'):
        """
        From a list of co-occurring values (one per variable) create a joint probability distribution by simply
        counting. For instance, the list [[0,0], [0,1], [1,0], [1,1]] would result in a uniform distribution of two
        binary variables.
        :param repeated_measurements: list of lists, where each sublist is of equal length (= number of variables)
        :type repeated_measurements: list of list
        :param discrete: do not change
        :param method: do not change
        """
        assert discrete and method == 'empirical', 'method/discreteness combination not supported'

        assert len(repeated_measurements) > 0, 'no data given'

        numvars = len(repeated_measurements[0])

        all_unique_values = list(set(np.array(repeated_measurements).flat))

        # todo: store the unique values as a member field so that later algorithms know what original values
        # corresponded to the values 0, 1, 2, ...

        numvals = len(all_unique_values)

        dict_val_to_index = {all_unique_values[valix]: valix for valix in xrange(numvals)}

        new_joint_probs = np.zeros([numvals]*numvars)

        # todo: when you generalize self.numvalues to an array then here also set the array instead of int

        for values in repeated_measurements:
            value_indices = tuple((dict_val_to_index[val] for val in values))

            try:
                new_joint_probs[value_indices] += 1.0 / len(repeated_measurements)
            except IndexError as e:
                print 'error: value_indices =', value_indices
                print 'error: type(value_indices) =', type(value_indices)

                raise IndexError(e)

        self.reset(numvars, numvals, joint_prob_matrix=new_joint_probs)


    def joint_probability(self, values):
        assert len(values) == self.numvariables, 'should specify one value per variable'
        assert values[0] < self.numvalues, 'variable can only take values 0, 1, ..., <numvalues - 1>: ' + str(values[0])
        assert values[-1] < self.numvalues, 'variable can only take values 0, 1, ..., <numvalues - 1>: ' \
                                            + str(values[-1])

        joint_prob = self.joint_probabilities[tuple(values)]

        assert 0.0 <= joint_prob <= 1.0, 'not a probability? ' + str(joint_prob)

        return joint_prob


    def get_variables_with_label_in(self, labels):
        return [vix for vix in xrange(self.numvariables) if self.labels[vix] in labels]


    def marginalize_distribution_retaining_only_labels(self, retained_labels):
        """
        Return a pdf of variables which have one of the labels in the retained_labels set; all other variables will
        be summed out.
        :param retained_labels: list of labels to retain; the rest will be summed over.
        :type retained_labels: sequence
        :rtype: JointProbabilityMatrix
        """

        variable_indices = [vix for vix in xrange(self.numvariables) if self.labels[vix] in retained_labels]

        return self.marginalize_distribution(variable_indices)


    def marginalize_distribution(self, retained_variables):
        """
        Return a pdf of only the given variable indices, summing out all others
        :param retained_variables: variables to retain, all others will be summed out and will not be a variable anymore
        :type: array_like
        :rtype : JointProbabilityMatrix
        """
        lists_of_possible_states_per_variable = [range(self.numvalues) for variable in xrange(self.numvariables)]

        assert hasattr(retained_variables, '__len__'), 'should be list or so, not int or other scalar'

        marginalized_joint_probs = np.zeros([self.numvalues]*len(retained_variables))

        # if len(variables):
        #     marginalized_joint_probs = np.array([marginalized_joint_probs])

        assert len(retained_variables) > 0, 'makes no sense to marginalize 0 variables'
        assert np.all(map(np.isscalar, retained_variables)), 'each variable identifier should be int in [0, numvalues)'
        assert len(retained_variables) <= self.numvariables, 'cannot marginalize more variables than I have'
        assert len(set(retained_variables)) <= self.numvariables, 'cannot marginalize more variables than I have'

        # not sure yet about doing this:
        # variables = sorted(list(set(variables)))  # ensure uniqueness?

        for values in itertools.product(*lists_of_possible_states_per_variable):
            marginal_values = [values[varid] for varid in retained_variables]

            marginalized_joint_probs[tuple(marginal_values)] += self.joint_probability(values)

        np.testing.assert_almost_equal(np.sum(marginalized_joint_probs), 1.0)

        marginal_joint_pdf = JointProbabilityMatrix(len(retained_variables), self.numvalues,
                                                    joint_probs=marginalized_joint_probs)

        return marginal_joint_pdf


    # helper function
    def appended_joint_prob_matrix(self, num_added_variables, values_so_far=[], added_joint_probabilities=None):
        if len(values_so_far) == self.numvariables:
            joint_prob_values = self.joint_probability(values_so_far)

            # submatrix must sum up to joint probability
            if added_joint_probabilities is None:
                added_joint_probabilities = np.random.random([self.numvalues]*num_added_variables)
                added_joint_probabilities /= np.sum(added_joint_probabilities)
                added_joint_probabilities *= joint_prob_values

                assert joint_prob_values <= 1.0
            else:
                np.testing.assert_almost_equal(np.sum(added_joint_probabilities), joint_prob_values)

                assert np.ndim(added_joint_probabilities) == num_added_variables
                assert len(added_joint_probabilities[0]) == self.numvalues
                assert len(added_joint_probabilities[-1]) == self.numvalues

            return list(added_joint_probabilities)
        elif len(values_so_far) < self.numvariables:
            if len(values_so_far) > 0:
                return [self.appended_joint_prob_matrix(num_added_variables,
                                                        values_so_far=list(values_so_far) + [val],
                                                        added_joint_probabilities=added_joint_probabilities)
                        for val in xrange(self.numvalues)]
            else:
                # same as other case but np.array converted, since the joint pdf matrix is always expected to be that
                return np.array([self.appended_joint_prob_matrix(num_added_variables,
                                                               values_so_far=list(values_so_far) + [val],
                                                               added_joint_probabilities=added_joint_probabilities)
                                 for val in xrange(self.numvalues)])
        else:
            raise RuntimeError('should not happen?')


    def append_variables(self, num_added_variables, added_joint_probabilities=None):
        assert num_added_variables > 0

        if isinstance(added_joint_probabilities, JointProbabilityMatrix):
            added_joint_probabilities = added_joint_probabilities.joint_probabilities

        new_joint_pdf = self.appended_joint_prob_matrix(num_added_variables,
                                                        added_joint_probabilities=added_joint_probabilities)

        assert np.ndim(new_joint_pdf) == self.numvariables + num_added_variables
        if self.numvariables + num_added_variables >= 1:
            assert len(new_joint_pdf[0]) == self.numvalues
            assert len(new_joint_pdf[-1]) == self.numvalues
        np.testing.assert_almost_equal(np.sum(new_joint_pdf), 1.0)

        # self.joint_probabilities = new_joint_pdf
        # self.numvariables = self.numvariables + num_added_variables
        self.reset(self.numvariables + num_added_variables, self.numvalues, new_joint_pdf)


    def append_variables_using_state_transitions_table(self, state_transitions):
        """
        Append one or more stochastic variables to this joint pdf, whose conditional pdf is defined by the provided
        'state transitions table'. In the rows of this table the first <self.numvariables> values are the values for
        the pre-existing stochastic variables; the added values are taken to be the deterministically determined
        added variable values, i.e., Pr(appended_vars = X, current_vars) = Pr(current_vars) so that
        Pr(appended_vars = X | current_vars) = 1 for whatever X you appended, where X is a list of values.
        :param state_transitions: list of lists, where each sublist is at least of length self.numvariables + 1 and
        is a list of values, each value in [0, numvalues).

        Can also provide a function f(values, num_values) which returns a list of values for the to-be-appended
        stochastic variables, where the argument <values> is a list of values for the existing variables (length
        self.numvariables).
        :type state_transitions: list or function
        """

        lists_of_possible_given_values = [range(self.numvalues) for _ in xrange(self.numvariables)]

        if hasattr(state_transitions, '__call__'):
            state_transitions = [list(existing_vars_values) + list(state_transitions(existing_vars_values,
                                                                                     self.numvalues))
                                 for existing_vars_values in itertools.product(*lists_of_possible_given_values)]

        extended_joint_probs = np.zeros([self.numvalues]*len(state_transitions[0]))

        # todo this for debugging? cycle through all possible values for self.numvariables and see if it is present
        # in the state_transitions
        # lists_of_possible_states_per_variable = [range(self.numvalues) for variable in xrange(self.numvariables)]

        # one row should be included for every possible set of values for the pre-existing stochastic variables
        assert len(state_transitions) == np.power(self.numvalues, self.numvariables)

        for states_row in state_transitions:
            assert len(states_row) > self.numvariables, 'if appending then more than self.numvariables values ' \
                                                        'should be specified'
            assert len(states_row) == len(state_transitions[0]), 'not all state rows of equal length; ' \
                                                                 'appending how many variables? Make up your mind. '

            # probability that the <self.numvariables> of myself have the values <state_transitions[:self.numvariables]>
            curvars_prob = self(states_row[:self.numvariables])

            assert 0.0 <= curvars_prob <= 1.0, 'probability not in 0-1'

            # set Pr(appended_vars = X | current_vars) = 1 for one set of values for the appended variables (X) and 0
            # otherwise (which is already by default), so I setting here
            # Pr(appended_vars = X, current_vars) = Pr(current_vars)
            extended_joint_probs[tuple(states_row)] = curvars_prob

        assert np.ndim(extended_joint_probs) == len(state_transitions[0])
        if len(state_transitions[0]) > 1:
            assert len(extended_joint_probs[0]) == self.numvalues
            assert len(extended_joint_probs[-1]) == self.numvalues
        np.testing.assert_almost_equal(np.sum(extended_joint_probs), 1.0)

        # self.joint_probabilities = extended_joint_probs
        # self.numvariables = len(state_transitions[0])
        self.reset(len(state_transitions[0]), self.numvalues, extended_joint_probs)


    def reorder_variables(self, variables):
        """
        Reorder the variables, for instance if self.numvariables == 3 then call with variables=[2,1,0] to reverse the
        order of the variables. The new joint probability matrix will be determined completely by the old matrix.
        It is also possible to duplicate variables, e.g. variables=[0,1,2,2] to duplicate the last variable.
        :param variables: sequence of int
        """
        assert len(variables) >= self.numvariables, 'I cannot reorder if you do\'nt give me the new ordering completely'
        # note: the code is mostly written to support also duplicating a variable

        num_variables_new = len(variables)

        joint_probs_new = np.zeros([self.numvalues]*num_variables_new) - 1

        lists_of_possible_states_per_variable = [range(self.numvalues) for _ in xrange(num_variables_new)]

        for values_new in itertools.product(*lists_of_possible_states_per_variable):
            values_old_order = [-1]*self.numvariables

            for new_varix in xrange(len(variables)):
                assert variables[new_varix] < self.numvariables, 'you specified the order of a variable index' \
                                                                 ' >= N (non-existent)'

                # can happen if a variable index is mentioned twice or more in 'variables' but in the current 'values_new'
                # they have different values. This is of course not possible, the two new variables should be equivalent
                # and completely redundant, so then I will set this joint prob. to zero and continue
                if values_old_order[variables[new_varix]] >= 0 \
                        and values_old_order[variables[new_varix]] != values_new[new_varix]:
                    assert len(variables) > self.numvariables, 'at least one variable index should be mentioned twice'

                    joint_probs_new[tuple(values_new)] = 0.0

                    break
                else:
                    # normal case
                    values_old_order[variables[new_varix]] = values_new[new_varix]

            if joint_probs_new[tuple(values_new)] != 0.0:
                assert not -1 in values_old_order, 'missing a variable index in variables=' + str(variables) \
                                                   + ', how should I reorder if you don\'t specify the new order of all' \
                                                   +  ' variables'

                assert joint_probs_new[tuple(values_new)] == -1, 'should still be unset element of joint prob matrix'

                joint_probs_new[tuple(values_new)] = self(values_old_order)
            else:
                pass  # joint prob already set to 0.0 in the above inner loop

        assert not -1 in joint_probs_new, 'not all new joint probs were computed'

        # change myself to the new joint pdf; will also check for being normalized etc.
        self.reset(num_variables_new, self.numvalues, joint_probs_new)


    def __eq__(self, other):  # approximate to 7 decimals
        if self.numvariables != other.numvariables or self.numvalues != other.numvalues:
            return False
        else:
            try:
                np.testing.assert_array_almost_equal(self.joint_probabilities, other.joint_probabilities)

                return True
            except AssertionError as e:
                assert 'not almost equal' in str(e), 'don\'t know what other assertion could have failed'

                return False


    # todo: implement __sub__, but for this it seems necessary to give each variable a unique ID at creation (__init__)
    # and keep track of them as you do operations such as marginalizing. Otherwise, subtraction is ambiguous,
    # i.e., subtracting a 2-variable pdf from a 5-variable pdf should result in a 3-variable conditional pdf, but
    # it is not clear which 3 variables should be chosen. Can choose this to always be the first 3 and let the
    # user be responsible to reorder them before subtraction, but seems prone to error if user does not know this or
    # forgets? Well can test equality by marginalizing the 5-variable pdf...
    def __sub__(self, other):
        assert len(self) >= len(other), 'cannot compute a conditional pdf consisting of a negative number of variables'

        if len(self) == len(other):
            assert self[range(len(other))] == other, 'my first ' + str(len(other)) + ' variables are not the same as ' \
                                                                                     'that of \'other\''

            return JointProbabilityMatrix(0, self.numvalues)  # return empty pdf...
        elif len(self) > len(other):
            assert self[range(len(other))] == other, 'my first ' + str(len(other)) + ' variables are not the same as ' \
                                                                                     'that of \'other\''

            return self.conditional_probability_distributions(range(len(other)))
        else:
            raise ValueError('len(self) < len(other), '
                             'cannot compute a conditional pdf consisting of a negative number of variables')


    # todo: implement __setitem__ for either pdf or cond. pdfs


    def __len__(self):
        return self.numvariables


    def __getitem__(self, item):
        if not hasattr(item, '__iter__'):
            return self.marginalize_distribution([item])
        else:
            return self.marginalize_distribution(item)


    def __add__(self, other):
        """

        :param other: can be JointProbabilityMatrix or a conditional distribution (dict of JointProbabilityMatrix)
        :type other: JointProbabilityMatrix or dict
        """

        self.append_variables_using_conditional_distributions(other)


    def matrix2vector(self):
        return self.joint_probabilities.flatten()


    def vector2matrix(self, list_probs):
        np.testing.assert_almost_equal(np.sum(list_probs), 1.0)

        assert np.ndim(list_probs) == 1

        self.joint_probabilities = np.reshape(list_probs, [self.numvalues]*self.numvariables)

        self.clip_all_probabilities()


    def params2matrix(self, parameters):
        assert len(parameters) == np.power(self.numvalues, self.numvariables) - 1

        vector_probs = [-1.0]*(np.power(self.numvalues, self.numvariables))

        remaining_prob_mass = 1.0

        for pix in xrange(len(parameters)):
            # note: small rounding errors will be fixed below by clipping
            assert -0.000001 <= parameters[pix] <= 1.000001, 'parameters should be in [0, 1]: ' + str(parameters[pix])

            # clip the parameter to the allowed range. If a rounding error is fixed by this in the parameters then
            # possibly a rounding error will appear in the probabilities?... Not sure though
            parameters[pix] = min(max(parameters[pix], 0.0), 1.0)

            vector_probs[pix] = remaining_prob_mass * parameters[pix]

            remaining_prob_mass = remaining_prob_mass * (1.0 - parameters[pix])

        assert vector_probs[-1] < 0.0, 'should still be unset by the above loop'

        # last parameter is irrelevant, must always be 1.0 is also a way to look at it
        vector_probs[-1] = remaining_prob_mass

        np.testing.assert_almost_equal(np.sum(vector_probs), 1.0)

        self.vector2matrix(vector_probs)


    def from_params(self, parameters):  # alternative constructor
        self.params2matrix(parameters)

        return self


    def matrix2params(self):
        vector_probs = self.matrix2vector()

        remaining_prob_mass = 1.0

        parameters = [-1.0]*(len(vector_probs) - 1)

        np.testing.assert_almost_equal(np.sum(vector_probs), 1.0)

        for pix in xrange(len(parameters)):
            if remaining_prob_mass > 0:
                assert remaining_prob_mass <= 1.0, 'remaining prob mass: ' + str(remaining_prob_mass)
                assert vector_probs[pix] <= remaining_prob_mass + 0.00001, \
                    'vector_probs[pix]=' + str(vector_probs[pix]) + ', remaining_prob_mass=' + str(remaining_prob_mass)

                parameters[pix] = vector_probs[pix] / remaining_prob_mass

                assert -0.000001 <= parameters[pix] <= 1.0000001, \
                    'parameters should be in [0, 1]: ' + str(parameters[pix])
            elif remaining_prob_mass == 0:
                parameters[pix] = 0
            else:
                if not remaining_prob_mass > -0.000001:
                    print 'debug: remaining_prob_mass =', remaining_prob_mass
                    print 'debug: pix =', pix, 'out of', len(parameters)

                    # todo: if just due to floating operation error, so small, then clip to zero and go on?
                    raise ValueError('remaining_prob_mass = ' + str(remaining_prob_mass)
                                     + ' < 0, which should not happen?')
                else:
                    # seems that it was intended to reach zero but due to floating operation roundoff it got just
                    # slightly under. Clip to 0 will do the trick.
                    remaining_prob_mass = 0.0  # clip to zero, so it will stay that way

                parameters[pix] = 0  # does not matter

            assert -0.000001 <= parameters[pix] <= 1.0000001, 'parameters should be in [0, 1]: ' + str(parameters[pix])

            parameters[pix] = min(max(parameters[pix], 0.0), 1.0)

            remaining_prob_mass -= remaining_prob_mass * parameters[pix]

        return parameters


    def __add__(self, other):
        """
        Append the variables defined by the (conditional) distributions in other.
        :type other: dict of JointProbabilityMatrix | JointProbabilityMatrix
        :rtype: JointProbabilityMatrix
        """

        pdf = self.copy()
        pdf.append_variables_using_conditional_distributions(other)

        return pdf

    def matrix2params_incremental(self, return_flattened=True, verbose=False):
        if self.numvariables > 1:
            # get the marginal pdf for the first variable
            pdf1 = self.marginalize_distribution([0])

            # first sequence of parameters, rest is added below here
            parameters = pdf1.matrix2params()

            pdf_conds = self.conditional_probability_distributions([0])

            # todo: remove this test, I am trying to find a bug now
            if np.random.randint(5) == 0:
                assert pdf1 + pdf_conds == self, 'should still be same distribution taken together'

            assert len(pdf_conds) == self.numvalues, 'should be one pdf for each value of first variable'

            for val in xrange(self.numvalues):
                pdf_cond = pdf_conds[tuple([val])]

                added_params = pdf_cond.matrix2params_incremental(return_flattened=False, verbose=verbose)

                if verbose:
                    print 'debug: matrix2params_incremental: recursed: for val=' + str(val) + ' I got added_params=' \
                          + str(added_params) + '.'
                    print 'debug: matrix2params_incremental: old parameters =', parameters

                # instead of returning a flat list of parameters I make it nested, so that the structure (e.g. number of
                # variables and number of values) can be inferred, and also hopefully it can be inferred to which
                # variable which parameters belong.
                # CHANGE123
                parameters.append(added_params)

                if verbose:
                    print 'debug: matrix2params_incremental: new parameters =', parameters

            if return_flattened:
                # flatten the tree structure to a list of scalars, which is sorted on the variable id
                parameters = self.scalars_up_to_level(parameters)

            return parameters
        elif self.numvariables == 1:
            return self.matrix2params()
        else:
            raise ValueError('no parameters for 0 variables')

    _debug_params2matrix = False  # internal variable, used to debug a debug statement, can be removed in a while


    def params2matrix_incremental(self, parameters):

        if __debug__:
            # store the original provided list of scalars
            original_parameters = list(parameters)

        # I suspect that both a tree-like input and a list of scalars should work... (add to unit test?)
        if np.all(map(np.isscalar, parameters)):
            assert min(parameters) > -0.0001, 'parameter(s) significantly out of allowed bounds [0,1]: ' \
                                              + str(parameters)
            assert min(parameters) < 1.0001, 'parameter(s) significantly out of allowed bounds [0,1]: ' \
                                              + str(parameters)

            # clip each parameter value to the allowed range. above I check already whether the error is not too large
            parameters = [min(max(pi, 0.0), 1.0) for pi in parameters]

            parameters = self.imbalanced_tree_from_scalars(parameters, self.numvalues)

            # verify that the procedure to make the tree out of the list of scalars is reversible and correct
            # (looking for bug)
            if __debug__ and self._debug_params2matrix:
                original_parameters2 = self.scalars_up_to_level(parameters)

                np.testing.assert_array_almost_equal(original_parameters, original_parameters2)

        if self.numvariables > 1:
            # first (numvalues - 1) values in the parameters tree structure should be scalars, as they will be used
            # to make the first variable's marginal distribution
            assert np.all(map(np.isscalar, parameters[:(self.numvalues - 1)]))

            ### start already by setting the pdf of the first variable...

            pdf_1 = JointProbabilityMatrix(1, self.numvalues)
            pdf_1.params2matrix(parameters[:(len(pdf_1.joint_probabilities.flat) - 1)])

            assert (len(pdf_1.joint_probabilities.flat) - 1) == (self.numvalues - 1), 'assumption directly above'

            assert len(pdf_1.joint_probabilities.flat) == self.numvalues

            assert len(flatten(parameters)) == len(self.joint_probabilities.flat) - 1, \
                'more or fewer parameters than needed: ' \
                              'need ' + str(len(self.joint_probabilities.flat) - 1) + ', got ' + str(len(parameters)) \
                              + '; #vars, #vals = ' + str(self.numvariables) + ', ' + str(self.numvalues)

            if __debug__ and self._debug_params2matrix:
                # remove this (expensive) check after it seems to work a few times?
                # note: for the conditions of no 1.0 or 0.0 prior probs, see the note in params2matrix_incremental
                if not 0.0 in pdf_1.matrix2params() and not 1.0 in pdf_1.matrix2params():
                    np.testing.assert_array_almost_equal(pdf_1.matrix2params(),
                                                         self.scalars_up_to_level(parameters[:(self.numvalues - 1)]))
                    np.testing.assert_array_almost_equal(pdf_1.matrix2params_incremental(),
                                                         self.scalars_up_to_level(parameters[:(self.numvalues - 1)]))

            # remove the used parameters from the list
            parameters = parameters[(len(pdf_1.joint_probabilities.flat) - 1):]
            assert len(parameters) == self.numvalues  # one subtree per conditional pdf

            pdf_conds = dict()

            ### now add other variables...

            for val in xrange(self.numvalues):
                # set this conditional pdf recursively as defined by the next sequence of parameters
                pdf_cond = JointProbabilityMatrix(self.numvariables - 1, self.numvalues)
                # CHANGE123

                # nte: parameters[0] is a sublist
                assert not np.isscalar(parameters[0])

                assert not np.isscalar(parameters[0])

                # todo: changing the parameters list is not necessary, maybe faster if not?

                # pdf_cond.params2matrix_incremental(parameters[:(len(pdf_cond.joint_probabilities.flat) - 1)])
                pdf_cond.params2matrix_incremental(parameters[0])

                # conditional pdf should have the same set of parameters as the ones I used to create it
                # (todo: remove this expensive check if it seems to work for  while)
                if self._debug_params2matrix:  # seemed to work for a long time...
                    try:
                        # todo: test can be probabilistic, e.g. performed with 10% chance?
                        np.testing.assert_array_almost_equal(pdf_cond.matrix2params_incremental(),
                                                             self.scalars_up_to_level(parameters[0]))
                    except AssertionError as e:
                        # print 'debug: parameters[0] =', parameters[0]
                        # print 'debug: len(pdf_cond) =', len(pdf_cond)
                        # print 'debug: pdf_cond.joint_probabilities =', pdf_cond.joint_probabilities

                        pdf_1_duplicate1 = pdf_cond.copy()
                        pdf_1_duplicate2 = pdf_cond.copy()

                        pdf_1_duplicate1._debug_params2matrix = False  # prevent endless recursion
                        pdf_1_duplicate2._debug_params2matrix = False  # prevent endless recursion

                        pdf_1_duplicate1.params2matrix_incremental(self.scalars_up_to_level(parameters[0]))
                        pdf_1_duplicate2.params2matrix_incremental(pdf_cond.matrix2params_incremental())

                        pdf_1_duplicate1._debug_params2matrix = True
                        pdf_1_duplicate2._debug_params2matrix = True

                        assert pdf_1_duplicate1 == pdf_cond
                        assert pdf_1_duplicate2 == pdf_cond

                        del pdf_1_duplicate1, pdf_1_duplicate2

                        # note: the cause seems to be as follows. If you provide the parameters e.g.
                        # [0.0, [0.37028884415935004], [0.98942830522914993]] then the middle parameter is superfluous,
                        # because it defines a conditional probability p(b|a) for which its prior p(a)=0. So no matter
                        # what parameter is here, the joint prob p(a,b) will be zero anyway. In other words, many
                        # parameter lists map to the same p(a,b). This makes the inverse problem ambiguous:
                        # going from a p(a,b) to a parameter list. So after building a pdf from the above example of
                        # parameter values I may very well get a different parameter list from that pdf, even though
                        # the pdf built is the one intended. I don't see a way around this because even if this class
                        # makes it uniform, e.g. always making parameter values 0.0 in case their prior is zero,
                        # but then still a user or optimization procedure can provide any list of parameters, so
                        # then also the uniformized parameter list will differ from the user-supplied.

                        # raise AssertionError(e)

                        # later add check. If this check fails then for sure there is something wrong. See also the
                        # original check below.
                        assert 0.0 in self.scalars_up_to_level(parameters) or \
                               1.0 in self.scalars_up_to_level(parameters), 'see story above. ' \
                                                                               'self.scalars_up_to_level(parameters) = ' \
                                                                               + str(self.scalars_up_to_level(parameters))

                        # original check. This check failed once, but the idea was to see if there are 0s or 1s in the
                        # prior probability distribution, which precedes the conditional probability distribution for which
                        # apparently the identifying parameter values have changed. But maybe I am wrong in that
                        # parameters[0] is not the prior only, and some prior prob. information is in all of parameters,
                        # I am not sure anymore so I added the above check to see whether that one is hit instead
                        # of this one (so above check is of course more stringent than this one....)
                        # assert 0.0 in self.scalars_up_to_level(parameters[0]) or \
                        #        1.0 in self.scalars_up_to_level(parameters[0]), 'see story above. ' \
                        #                                                        'self.scalars_up_to_level(parameters[0]) = ' \
                        #                                                        + str(self.scalars_up_to_level(parameters[0]))

                # right?
                np.testing.assert_almost_equal(np.sum(pdf_cond.joint_probabilities), 1.0)

                # remove the used parameters from the list
                # CHANGE123
                # parameters = parameters[(len(pdf_cond.joint_probabilities.flat) - 1):]
                parameters = parameters[1:]

                # add the conditional pdf
                pdf_conds[(val,)] = pdf_cond.copy()

            assert len(parameters) == 0, 'all parameters should be used to construct joint pdf'

            pdf_1.append_variables_using_conditional_distributions(pdf_conds)

            if __debug__ and self._debug_params2matrix:
                # remove this (expensive) check after it seems to work a few times?
                try:
                    np.testing.assert_array_almost_equal(pdf_1.matrix2params_incremental(),
                                                         self.scalars_up_to_level(original_parameters))
                except AssertionError as e:
                    ### I have the hunch that the above assertion is hit but that it is only if a parameter is 1 or 0,
                    ### so that the parameter may be different but that it does not matter. still don't understand
                    ### why it happens though...

                    pdf_1_duplicate = pdf_1.copy()

                    pdf_1_duplicate._debug_params2matrix = False  # prevent endless recursion

                    pdf_1_duplicate.params2matrix_incremental(self.scalars_up_to_level(original_parameters))

                    pdf_1_duplicate._debug_params2matrix = True

                    if not pdf_1_duplicate == pdf_1:
                        print 'error: the pdfs from the two different parameter lists are also not equivalent'

                        del pdf_1_duplicate

                        raise AssertionError(e)
                    else:
                        # warnings.warn('I found that two PDF objects can have the same joint prob. matrix but a'
                        #               ' different list of identifying parameters. This seems to be due to a variable'
                        #               ' having 0.0 probability on a certain value, making the associated conditional'
                        #               ' PDF of other variables 0 and therefore those associated parameters irrelevant.'
                        #               ' Find a way to make these parameters still uniform? Seems to happen in'
                        #               ' "pdf_1.append_variables_using_conditional_distributions(pdf_conds)"...')

                        # note: (duplicated) the cause seems to be as follows. If you provide the parameters e.g.
                        # [0.0, [0.37028884415935004], [0.98942830522914993]] then the middle parameter is superfluous,
                        # because it defines a conditional probability p(b|a) for which its prior p(a)=0. So no matter
                        # what parameter is here, the joint prob p(a,b) will be zero anyway. In other words, many
                        # parameter lists map to the same p(a,b). This makes the inverse problem ambiguous:
                        # going from a p(a,b) to a parameter list. So after building a pdf from the above example of
                        # parameter values I may very well get a different parameter list from that pdf, even though
                        # the pdf built is the one intended. I don't see a way around this because even if this class
                        # makes it uniform, e.g. always making parameter values 0.0 in case their prior is zero,
                        # but then still a user or optimization procedure can provide any list of parameters, so
                        # then also the uniformized parameter list will differ from the user-supplied.

                        del pdf_1_duplicate

            assert pdf_1.numvariables == self.numvariables
            assert pdf_1.numvalues == self.numvalues

            self.duplicate(pdf_1)  # make this object (self) be the same as pdf_1
        elif self.numvariables == 1:
            self.params2matrix(parameters)
        else:
            assert len(parameters) == 0, 'at the least 0 parameters should be given for 0 variables...'

            raise ValueError('no parameters for 0 variables')

    def duplicate(self, other_joint_pdf):
        """

        :type other_joint_pdf: JointProbabilityMatrix
        """
        self.reset(other_joint_pdf.numvariables, other_joint_pdf.numvalues, other_joint_pdf.joint_probabilities)

    def reset(self, numvariables, numvalues, joint_prob_matrix, labels=None):
        """
        This function is intended to completely reset the object, so if you add variables which determine the
        behavior of the object then add them also here and everywhere where called.
        :type numvariables: int
        :type numvalues: int
        :type joint_prob_matrix: JointProbabilityMatrix
        """

        self.numvariables = numvariables
        self.numvalues = numvalues
        self.joint_probabilities = joint_prob_matrix

        if labels is None:
            self.labels = [_default_variable_label]*numvariables
        else:
            self.labels = labels

        assert np.ndim(joint_prob_matrix) == self.numvariables, 'ndim = ' + str(np.ndim(joint_prob_matrix)) + ', ' \
                                                                'self.numvariables = ' + str(self.numvariables) + ', ' \
                                                                'joint matrix = ' + str(joint_prob_matrix)
        assert list(set(np.shape(self.joint_probabilities)))[0] == self.numvalues

        # # maybe this check should be removed, it is also checked in clip_all_* below, but after clipping, which
        # # may be needed to get this condition valid again?
        # if np.random.random() < 0.01:  # make less frequent
        #     np.testing.assert_array_almost_equal(np.sum(joint_prob_matrix), 1.0)

        self.clip_all_probabilities()


    def statespace(self, numvars='all'):
        if numvars == 'all':
            lists_of_possible_joint_values = [range(self.numvalues) for _ in xrange(self.numvariables)]
        elif type(numvars) in (int,):
            lists_of_possible_joint_values = [range(self.numvalues) for _ in xrange(numvars)]
        else:
            raise NotImplementedError('numvars=' + str(numvars))

        return itertools.product(*lists_of_possible_joint_values)


    def append_redundant_variables(self, num_redundant_variables):
        self.append_variables_using_state_transitions_table(lambda my_values_tuple, numvals:
                                                            [sum(my_values_tuple) % numvals]*num_redundant_variables)


    def append_independent_variables(self, joint_pdf):
        """

        :type joint_pdf: JointProbabilityMatrix
        """
        assert not type(joint_pdf) in (int, float, str), 'should pass a JointProbabilityMatrix object'

        self.append_variables_using_conditional_distributions({(): joint_pdf})


    def append_variables_using_conditional_distributions(self, pdf_conds, given_variables=None):
        """


        :param pdf_conds: dictionary of JointProbabilityMatrix objects, one for each possible set of values for the
        existing <self.numvariables> variables. If you provide a dict with one empty tuple as key then it is
        equivalent to just providing the pdf, or calling append_independent_variables with that joint pdf. If you
        provide a dict with tuples as keys which have fewer values than self.numvariables then the dict will be
        duplicated for all possible values for the remaining variables, and this function will call itself recursively
        with complete tuples.
        :type pdf_conds: dict of tuple | JointProbabilityMatrix
        :param given_variables: If you provide a dict with tuples as keys which have fewer values than
        self.numvariables then it is assumed that the tuples specify the values for the consecutive highest indexed
        variables (so always until variable self.numvariables-1). Unless you specify this list, whose length should
        be equal to the length of the tuple keys in pdf_conds
        :type given_variables: list of int
        """

        if isinstance(pdf_conds, JointProbabilityMatrix):
            pdf_conds = {tuple(my_values): pdf_conds for my_values in self.statespace()}
        else:
            assert type(pdf_conds) == dict

            num_conditioned_vars = len(pdf_conds.keys()[0])

            assert num_conditioned_vars <= len(self), 'makes no sense to condition on more variables than I have?'

            if num_conditioned_vars < len(self):
                # the provided conditional pdf conditions on fewer variables (m) than I currently exist of, so appending
                # the variables is strictly speaking ill-defined. What I will do is assume that the conditional pdf
                # is conditioned on the last m variables, and simply copy the same conditional pdf for the first
                # len(self) - m variables, which implies independence.

                pdf_conds_complete = dict()

                num_independent_vars = len(self) - num_conditioned_vars

                if __debug__:
                    # finding bug, check if all values of the dictionary are pdfs
                    for debug_pdf in pdf_conds.itervalues():
                        assert isinstance(debug_pdf, JointProbabilityMatrix), 'debug_pdf = ' + str(debug_pdf)

                statespace_per_independent_variable = [range(self.numvalues)
                                                       for _ in xrange(num_independent_vars)]

                if given_variables is None:
                    for indep_vals in itertools.product(*statespace_per_independent_variable):
                        pdf_conds_complete.update({(tuple(indep_vals) + tuple(key)): value
                                                   for key, value in pdf_conds.items()})
                else:
                    assert len(given_variables) == len(pdf_conds.keys()[0]), \
                        'if conditioned on ' + str(len(given_variables)) + 'then I also expect a conditional pdf ' \
                        + 'which conditions on ' + str(len(given_variables)) + ' variables.'

                    not_given_variables = np.setdiff1d(range(self.numvariables), given_variables)

                    assert len(not_given_variables) + len(given_variables) == self.numvariables

                    for indep_vals in itertools.product(*statespace_per_independent_variable):
                        pdf_conds_complete.update({tuple(apply_permutation(indep_vals + tuple(key),
                                                                           list(not_given_variables)
                                                                           + list(given_variables))): value
                                                   for key, value in pdf_conds.items()})

                assert len(pdf_conds_complete.keys()[0]) == len(self)

                # recurse once with a complete conditional pdf, so this piece of code should not be executed again:
                self.append_variables_using_conditional_distributions(pdf_conds=pdf_conds_complete)

                return

        assert type(pdf_conds) == dict

        num_added_variables = pdf_conds[(0,)*self.numvariables].numvariables

        assert num_added_variables > 0, 'makes no sense to append 0 variables?'
        assert self.numvalues == pdf_conds[(0,)*self.numvariables].numvalues, 'added variables should have same #values'

        # see if at the end, the new joint pdf has the expected list of identifying parameter values
        if __debug__ and len(self) == 1:
            _debug_parameters_before_append = self.matrix2params()

            # note: in the loop below the sublists of parameters will be added

        statespace_per_variable = [range(self.numvalues)
                                          for _ in xrange(self.numvariables + num_added_variables)]

        extended_joint_probs = np.zeros([self.numvalues]*(self.numvariables + num_added_variables))

        for values in itertools.product(*statespace_per_variable):
            existing_variables_values = values[:self.numvariables]
            added_variables_values = values[self.numvariables:]

            assert len(added_variables_values) == pdf_conds[tuple(existing_variables_values)].numvariables, 'error: ' \
                    'len(added_variables_values) = ' + str(len(added_variables_values)) + ', cond. numvariables = ' \
                    '' + str(pdf_conds[tuple(existing_variables_values)].numvariables) + ', len(values) = ' \
                    + str(len(values)) + ', existing # variables = ' + str(self.numvariables) + ', ' \
                    'num_added_variables = ' + str(num_added_variables)

            prob_existing = self(existing_variables_values)
            prob_added_cond_existing = pdf_conds[tuple(existing_variables_values)](added_variables_values)

            assert 0.0 <= prob_existing <= 1.0, 'prob not normalized'
            assert 0.0 <= prob_added_cond_existing <= 1.0, 'prob not normalized'

            extended_joint_probs[tuple(values)] = prob_existing * prob_added_cond_existing

            if __debug__ and len(self) == 1:
                _debug_parameters_before_append.append(pdf_conds[tuple(existing_variables_values)].matrix2params_incremental)

        np.testing.assert_almost_equal(np.sum(extended_joint_probs), 1.0)

        self.reset(self.numvariables + num_added_variables, self.numvalues, extended_joint_probs)

        if __debug__ and len(self) == 1:
            # of course this test depends on the implementation of matrix2params_incremental, currently it should
            # work
            np.testing.assert_array_almost_equal(self.scalars_up_to_level(_debug_parameters_before_append),
                                                 self.matrix2params_incremental(return_flattened=True))


    def conditional_probability_distribution(self, given_variables, given_values):
        """

        :param given_variables: list of integers
        :param given_values: list of integers
        :rtype: JointProbabilityMatrix
        """
        assert len(given_values) == len(given_variables)
        assert len(given_variables) < self.numvariables, 'no variables left after conditioning'

        lists_of_possible_states_per_variable = [range(self.numvalues) for variable in xrange(self.numvariables)]

        # overwrite the 'state spaces' for the specified variables, to the specified state spaces
        for gix in xrange(len(given_variables)):
            assert np.isscalar(given_values[gix]), 'assuming specific value, not list of possibilities'

            lists_of_possible_states_per_variable[given_variables[gix]] = \
                [given_values[gix]] if np.isscalar(given_values[gix]) else given_values[gix]

        conditioned_variables = [varix for varix in xrange(self.numvariables) if not varix in given_variables]

        conditional_probs = np.zeros([self.numvalues]*len(conditioned_variables))

        assert len(conditional_probs) > 0, 'you want me to make a conditional pdf of 0 variables?'

        assert len(given_variables) + len(conditioned_variables) == self.numvariables

        for values in itertools.product(*lists_of_possible_states_per_variable):
            values_conditioned_vars = [values[varid] for varid in conditioned_variables]

            assert conditional_probs[tuple(values_conditioned_vars)] == 0.0, 'right?'

            # note: here self.joint_probability(values) == Pr(conditioned_values | given_values) because the
            # given_variables == given_values constraint is imposed on the set
            # itertools.product(*lists_of_possible_states_per_variable); it is only not yet normalized
            conditional_probs[tuple(values_conditioned_vars)] += self.joint_probability(values)

        summed_prob_mass = np.sum(conditional_probs)

        # testing here if the summed prob. mass equals the marginal prob of the given variable values
        if __debug__:
            # todo: can make this test be run probabilistically, like 10% chance or so, pretty expensive?
            if np.all(map(np.isscalar, given_values)):
                pdf_marginal = self.marginalize_distribution(given_variables)

                prob_given_values = pdf_marginal(given_values)

                np.testing.assert_almost_equal(prob_given_values, summed_prob_mass)

        assert np.isscalar(summed_prob_mass), 'sum of probability mass should be a scalar, not a list or so: ' \
                                              + str(summed_prob_mass)
        assert np.isfinite(summed_prob_mass)

        assert summed_prob_mass >= 0.0, 'probability mass cannot be negative'

        if summed_prob_mass > 0.0:
            conditional_probs /= summed_prob_mass
        else:
            # note: apparently the probability of the given condition is zero, so it does not really matter
            # what I substitute for the probability mass here. I will add some probability mass so that I can
            # normalize it.

            # # I think this warning can be removed at some point...
            # warnings.warn('conditional_probability_distribution: summed_prob_mass == 0.0 (can be ignored)')

            # are all values zero?
            try:
                np.testing.assert_almost_equal(np.min(conditional_probs), 0.0)
            except ValueError as e:
                print 'debug: conditional_probs =', conditional_probs
                print 'debug: min(conditional_probs) =', min(conditional_probs)

                raise ValueError(e)

            np.testing.assert_almost_equal(np.max(conditional_probs), 0.0)

            conditional_probs *= 0
            conditional_probs += 1.0  # create some fake mass, making it a uniform distribution

            conditional_probs /= np.sum(conditional_probs)

        conditional_joint_pdf = JointProbabilityMatrix(len(conditioned_variables), self.numvalues,
                                                    joint_probs=conditional_probs)

        return conditional_joint_pdf


    def conditional_probability_distributions(self, given_variables):
        """

        :param given_variables:
        :return: dict of JointProbabilityMatrix, keys are all possible values for given_variables
        :rtype: dict of JointProbabilityMatrix
        """
        assert len(given_variables) < self.numvariables, 'no variables left after conditioning'

        lists_of_possible_given_values = [range(self.numvalues) for variable in xrange(len(given_variables))]

        return {tuple(given_values): self.conditional_probability_distribution(given_variables=given_variables,
                                                                               given_values=given_values)
                for given_values in itertools.product(*lists_of_possible_given_values)}


    # todo: try to make this entropy() function VERY efficient, or compiled (weave?), or C++ binding or something,
    # it is a central
    # function in all information-theoretic quantities and especially it is a crucial bottleneck for optimization
    # procedures such as in append_orthogonal* which call information-theoretic quantities a LOT.
    def entropy(self, variables=None):
        # def log2term(p):  # p log_2 p
        #     if 0.0 < p < 1.0:
        #         return -p * np.log2(p)
        #     else:
        #         return 0.0  # 0 log 0 == 0 is assumed

        if variables is None:
            if np.random.random() < 0.1:
                np.testing.assert_almost_equal(np.sum(self.joint_probabilities), 1.0)

            probs = self.joint_probabilities.flatten()
            probs = np.select([probs != 0], [probs], default=1)

            # joint_entropy = np.sum(map(log2term, self.joint_probabilities.flat))
            log_terms = probs * np.log2(probs)
            # log_terms = np.select([np.isfinite(log_terms)], [log_terms])
            joint_entropy = -np.sum(log_terms)

            assert joint_entropy >= 0.0
            # note: the extra 0.001 is for roundoff errors
            assert joint_entropy <= np.log2(self.numvalues) * self.numvariables + 0.001, \
                'joint_entropy=' + str(joint_entropy) + ', np.log2(self.numvalues) * self.numvariables = ' \
                + str(np.log2(self.numvalues) * self.numvariables)

            return joint_entropy
        else:
            assert hasattr(variables, '__iter__')

            if len(variables) == 0:  # hard-coded this because otherwise I have to support empty pdfs (len() = 0)
                return 0.0

            marginal_pdf = self.marginalize_distribution(retained_variables=variables)

            return marginal_pdf.entropy()


    def conditional_entropy(self, variables, given_variables=None):
        assert hasattr(variables, '__iter__'), 'variables1 = ' + str(variables)
        assert hasattr(given_variables, '__iter__') or given_variables is None, 'variables2 = ' + str(given_variables)

        assert max(variables) < self.numvariables, 'variables are 0-based indices, so <= N - 1: variables=' \
                                                   + str(variables) + ' (N=' + str(self.numvariables) + ')'

        if given_variables is None:
            # automatically set the given_variables to be the complement of 'variables', so all remaining variables

            given_variables = [varix for varix in xrange(self.numvariables) if not varix in variables]

            assert len(set(variables)) + len(given_variables) == self.numvariables, 'variables=' + str(variables) \
                                                                                    + ', given_variables=' \
                                                                                    + str(given_variables)

        # H(Y) + H(X|Y) == H(X,Y)
        condent = self.entropy(list(set(list(variables) + list(given_variables)))) - self.entropy(given_variables)

        assert np.isscalar(condent)
        assert np.isfinite(condent)

        assert condent >= 0.0, 'conditional entropy should be non-negative'

        return condent


    def mutual_information_labels(self, labels1, labels2):
        variables1 = self.get_variables_with_label_in(labels1)
        variables2 = self.get_variables_with_label_in(labels2)

        return self.mutual_information(variables1, variables2)


    def mutual_information(self, variables1, variables2):
        assert hasattr(variables1, '__iter__'), 'variables1 = ' + str(variables1)
        assert hasattr(variables2, '__iter__'), 'variables2 = ' + str(variables2)

        if len(variables1) == 0 or len(variables2) == 0:
            mi = 0  # trivial case, no computation needed
        elif len(variables1) == len(variables2):
            if np.equal(sorted(variables1), sorted(variables2)).all():
                mi = self.entropy(variables1)  # more efficient, shortcut
            else:
                # this one should be equal to the outer else-clause (below), i.e., the generic case
                mi = self.entropy(variables1) + self.entropy(variables2) \
                     - self.entropy(list(set(list(variables1) + list(variables2))))
        else:
            mi = self.entropy(variables1) + self.entropy(variables2) \
                 - self.entropy(list(set(list(variables1) + list(variables2))))

        assert np.isscalar(mi)
        assert np.isfinite(mi)

        # due to floating-point operations it might VERY sometimes happen that I get something like -4.4408920985e-16
        # here, so to prevent this case from firing an assert I clip this to 0:
        if -0.000001 < mi < 0.0:
            mi = 0.0

        assert mi >= 0, 'mutual information should be non-negative: ' + str(mi)

        return mi


    def synergistic_information_naive(self, variables_SRV, variables_X):
        return self.mutual_information(list(variables_SRV), list(variables_X)) \
               - sum([self.mutual_information(list(variables_SRV), list([var_xi])) for var_xi in variables_X])


    # todo: return not just a number but an object with more result information, which maybe if evaluated as
    # float then it will return the current return value
    def synergistic_information(self, variables_Y, variables_X, tolerance_nonsyn_mi=0.05, verbose=True,
                                minimize_method=None, num_repeats_per_srv_append=3):

        pdf_with_srvs = self.copy()

        syn_entropy = self.synergistic_entropy_upper_bound(variables_X)
        max_ent_per_var = np.log2(self.numvalues)
        max_num_srv_add_trials = int(round(syn_entropy / max_ent_per_var * 2))

        # note: currently I constrain to add SRVs which consist each of only 1 variable each. I guess in most cases
        # this will work fine, but I do not know for sure that this will always be sufficient (i.e., does there
        # ever an MSRV exist with bigger entropy, and if so, is this MSRV not decomposable into parts such that
        # this approach here still works?)

        for i in xrange(max_num_srv_add_trials):
            try:
                pdf_with_srvs.append_synergistic_variables(1, initial_guess_summed_modulo=False,
                                                           subject_variables=variables_X,
                                                           num_repeats=num_repeats_per_srv_append,
                                                           agnostic_about=range(len(self),
                                                                                len(pdf_with_srvs)),
                                                           minimize_method=minimize_method)
            except UserWarning as e:
                assert 'minimize() failed' in str(e), 'only known reason for this error'

                warnings.warn(str(e) + '. Will now skip this sample in synergistic_information.')

                continue

            # todo: can save one MI calculation here
            new_syn_info = pdf_with_srvs.synergistic_information_naive([-1], variables_X)
            total_mi = pdf_with_srvs.mutual_information([-1], variables_X)
            indiv_mi_list =  [pdf_with_srvs.mutual_information([-1], [xid]) for xid in variables_X]
            new_syn_info = total_mi - sum(indiv_mi_list)

            if new_syn_info < syn_entropy * 0.01:
                if (total_mi - new_syn_info) / total_mi > tolerance_nonsyn_mi:  # too much non-synergistic information?
                    pdf_with_srvs.marginalize_distribution(range(len(pdf_with_srvs) - 1))  # remove last srv

                    if verbose:
                        print 'debug: synergistic_information: a last SRV was found but with too much non-syn. info.'
                else:
                    if verbose:
                        print 'debug: synergistic_information: a last (small) ' \
                              'SRV was found, a good one, and now I will stop.'

                break  # assume that no more better SRVs can be found from now on, so stop (save time)
            else:
                if i == max_num_srv_add_trials - 1:
                    warnings.warn('synergistic_information: never converged to adding SRV with zero synergy')

                if (total_mi - new_syn_info) / total_mi > tolerance_nonsyn_mi:  # too much non-synergistic information?
                    pdf_with_srvs.marginalize_distribution(range(len(pdf_with_srvs) - 1))  # remove last srv

                    if verbose:
                        print('debug: synergistic_information: an SRV with new_syn_info/total_mi = '
                                  + str(new_syn_info) + ' / ' + str(total_mi) + ' = '
                                  + str((new_syn_info) / total_mi) + ' was found, which will be removed again because'
                                                                     ' it does not meet the tolerance of '
                                  + str(tolerance_nonsyn_mi))
                else:
                    if verbose:
                        print 'debug: synergistic_information: an SRV with new_syn_info/total_mi = ' \
                              + str(new_syn_info) + ' / ' + str(total_mi) + ' = ' \
                                  + str((new_syn_info) / total_mi) + ' was found, which satisfies the tolerance of ' \
                                  + str(tolerance_nonsyn_mi)

        if verbose:
            print 'debug: synergistic_information: number of SRVs:', len(xrange(len(self), len(pdf_with_srvs)))
            print 'debug: synergistic_information: entropy of SRVs:', \
                pdf_with_srvs.entropy(range(len(self), len(pdf_with_srvs)))
            print 'debug: synergistic_information: H(SRV_i) =', \
                [pdf_with_srvs.entropy([srv_id]) for srv_id in xrange(len(self), len(pdf_with_srvs))]
            print 'debug: synergistic_information: I(Y, SRV_i) =', \
                [pdf_with_srvs.mutual_information(variables_Y, [srv_id])
                        for srv_id in xrange(len(self), len(pdf_with_srvs))]

        syn_info = sum([pdf_with_srvs.mutual_information(variables_Y, [srv_id])
                        for srv_id in xrange(len(self), len(pdf_with_srvs))])

        return syn_info


    def synergistic_entropy(self, variables_X, tolerance_nonsyn_mi=0.05, verbose=True):

        assert False, 'todo: change from information to entropy here, just return different value and do not ask for Y'

        pdf_with_srvs = self.copy()

        syn_entropy = self.synergistic_entropy_upper_bound(variables_X)
        max_ent_per_var = np.log2(self.numvalues)
        max_num_srv_add_trials = int(round(syn_entropy / max_ent_per_var * 2))

        # note: currently I constrain to add SRVs which consist each of only 1 variable each. I guess in most cases
        # this will work fine, but I do not know for sure that this will always be sufficient (i.e., does there
        # ever an MSRV exist with bigger entropy, and if so, is this MSRV not decomposable into parts such that
        # this approach here still works?)

        for i in xrange(max_num_srv_add_trials):
            pdf_with_srvs.append_synergistic_variables(1, initial_guess_summed_modulo=False,
                                                       subject_variables=variables_X, num_repeats=3,
                                                       agnostic_about=range(len(variables_Y) + len(variables_X),
                                                                            len(pdf_with_srvs)))

            # todo: can save one MI calculation here
            new_syn_info = self.synergistic_information_naive([-1], variables_X)
            total_mi = self.mutual_information([-1], variables_X)

            if new_syn_info < syn_entropy * 0.01:
                if (total_mi - new_syn_info) / total_mi > tolerance_nonsyn_mi:  # too much non-synergistic information?
                    pdf_with_srvs.marginalize_distribution(range(len(pdf_with_srvs) - 1))  # remove last srv

                    if verbose:
                        print 'debug: synergistic_information: a last SRV was found but with too much non-syn. info.'
                else:
                    if verbose:
                        print 'debug: synergistic_information: a last (small) ' \
                              'SRV was found, a good one, and now I will stop.'

                break  # assume that no more better SRVs can be found from now on, so stop (save time)
            else:
                if i == max_num_srv_add_trials - 1:
                    warnings.warn('synergistic_information: never converged to adding SRV with zero synergy')

                if (total_mi - new_syn_info) / total_mi > tolerance_nonsyn_mi:  # too much non-synergistic information?
                    pdf_with_srvs.marginalize_distribution(range(len(pdf_with_srvs) - 1))  # remove last srv

                    if verbose:
                        print('debug: synergistic_information: an SRV with new_syn_info/total_mi = '
                                  + str(new_syn_info) + ' / ' + str(total_mi) + ' = '
                                  + str((new_syn_info) / total_mi) + ' was found, which will be removed again because'
                                                                     ' it does not meet the tolerance of '
                                  + str(tolerance_nonsyn_mi))
                else:
                    if verbose:
                        print 'debug: synergistic_information: an SRV with new_syn_info/total_mi = ' \
                              + str(new_syn_info) + ' / ' + str(total_mi) + ' = ' \
                                  + str((new_syn_info) / total_mi) + ' was found, which satisfies the tolerance of ' \
                                  + str(tolerance_nonsyn_mi)

        # if verbose:
        #     print 'debug: synergistic_information: number of SRVs:', len(xrange(len(self), len(pdf_with_srvs)))
        #     print 'debug: synergistic_information: entropy of SRVs:', \
        #         pdf_with_srvs.entropy(range(len(self), len(pdf_with_srvs)))
        #     print 'debug: synergistic_information: H(SRV_i) =', \
        #         [pdf_with_srvs.entropy([srv_id]) for srv_id in xrange(len(self), len(pdf_with_srvs))]
        #     print 'debug: synergistic_information: I(Y, SRV_i) =', \
        #         [pdf_with_srvs.mutual_information(variables_Y, [srv_id])
        #                 for srv_id in xrange(len(self), len(pdf_with_srvs))]

        syn_info = sum([pdf_with_srvs.mutual_information(variables_Y, [srv_id])
                        for srv_id in xrange(len(self), len(pdf_with_srvs))])

        return syn_info


    def synergistic_entropy_upper_bound(self, variables=None):
        """
        For a derivation, see Quax et al. (2015). No other stochastic variable can have more than this amount of
        synergistic information about the stochastic variables defined by this pdf object, or the given subset of
        variables in <variables>.
        :type variables: list of int
        :rtype: float
        """
        if variables is None:
            variables = range(len(self))

        indiv_entropies = [self.entropy([var]) for var in variables]

        return float(self.entropy() - max(indiv_entropies))


    # X marks the spot


    # todo: add optional numvalues? so that the synergistic variables can have more possible values than the
    # current variables (then set all probabilities where the original variables exceed their original max to 0)
    # todo: first you should then probably implement a .increase_num_values(num) or so (or only),
    # or let .numvalues be a list instead of a single value
    def append_synergistic_variables(self, num_synergistic_variables, initial_guess_summed_modulo=False, verbose=False,
                                     subject_variables=None, agnostic_about=None, num_repeats=1, minimize_method=None):
        """
        Append <num_synergistic_variables> variables in such a way that they are agnostic about any individual
        existing variable (one of self.numvariables thus) but have maximum MI about the set of self.numvariables
        variables taken together.
        :param agnostic_about: a list of variable indices to which the new synergistic variable set should be
        agnostic (zero mutual information). This can be used to find a 'different', second SRV after having found
        already a first one. If already multiple SRVs have been found then you can choose between either being agnostic
        about these previous SRVs jointly (so also possibly any synergy among them), or being 'pairwise' agnostic
        to each individual one, in which case you can pass a list of lists, then I will compute the added cost for
        each sublist and sum it up.
        :param num_repeats:
        :param subject_variables: the list of variables the <num_synergistic_variables> should be synergistic about;
        then I think the remainder variables the <num_synergistic_variables> should be agnostic about. This way I can
        support adding many UMSRVs (maybe make a new function for that) which then already are orthogonal among themselves,
        meaning I do not have to do a separate orthogonalization for the MSRVs as in the paper's theoretical part.
        :param num_synergistic_variables:
        :param initial_guess_summed_modulo:
        :param verbose:
        :return:
        """

        if not agnostic_about is None:
            if len(agnostic_about) == 0:
                agnostic_about = None  # treat as if not supplied

        if __debug__:  # looking for bug
            assert max(self.matrix2params_incremental()) < 1.0000001, 'param is out of bound: ' \
                                                                      + str(max(self.matrix2params_incremental()))
            assert min(self.matrix2params_incremental()) > -0.0000001, 'param is out of bound: ' \
                                                                      + str(min(self.matrix2params_incremental()))

        parameter_values_before = list(self.matrix2params_incremental())

        # if __debug__:
        #     debug_params_before = copy.deepcopy(parameter_values_before)

        # a pdf with XORs as appended variables (often already MSRV for binary variables), good initial guess?
        pdf_with_srvs = self.copy()
        pdf_with_srvs.append_variables_using_state_transitions_table(
            state_transitions=lambda vals, mv: [int(np.mod(np.sum(vals), mv))]*num_synergistic_variables)

        assert pdf_with_srvs.numvariables == self.numvariables + num_synergistic_variables

        parameter_values_after = pdf_with_srvs.matrix2params_incremental()

        assert len(parameter_values_after) > len(parameter_values_before), 'should be additional free parameters'
        if np.random.random < 0.1:  # reduce slowdown from this
            # see if the first part of the new parameters is exactly the old parameters, so that I know that I only
            # have to optimize the latter part of parameter_values_after
            np.testing.assert_array_almost_equal(parameter_values_before,
                                                 parameter_values_after[:len(parameter_values_before)])

        # this many parameters (each in [0,1]) must be optimized
        num_free_parameters = len(parameter_values_after) - len(parameter_values_before)

        assert num_synergistic_variables == 0 or num_free_parameters > 0

        if initial_guess_summed_modulo:
            # note: this is xor for binary variables
            initial_guess = parameter_values_after[len(parameter_values_before):]
        else:
            initial_guess = np.random.random(num_free_parameters)

        if verbose:
            debug_pdf_with_srvs = pdf_with_srvs.copy()
            debug_pdf_with_srvs.params2matrix_incremental(list(parameter_values_before) + list(initial_guess))

            # store the synergistic information before the optimization procedure (after procedure should be higher...)
            debug_before_syninfo = debug_pdf_with_srvs.synergistic_information_naive(
                variables_SRV=range(self.numvariables, pdf_with_srvs.numvariables),
                variables_X=range(self.numvariables))

        assert len(initial_guess) == num_free_parameters

        pdf_with_srvs_for_optimization = pdf_with_srvs.copy()

        if not subject_variables is None:
            pdf_subjects_syns_only = pdf_with_srvs_for_optimization.marginalize_distribution(
                list(subject_variables) + range(len(pdf_with_srvs) - num_synergistic_variables, len(pdf_with_srvs)))

            pdf_subjects_only = pdf_subjects_syns_only.marginalize_distribution(range(len(subject_variables)))

            if __debug__ and np.random.random() < 0.01:
                debug_pdf_subjects_only = pdf_with_srvs.marginalize_distribution(subject_variables)

                assert debug_pdf_subjects_only == pdf_subjects_only

            num_free_parameters_synonly = len(pdf_subjects_syns_only.matrix2params_incremental()) \
                                          - len(pdf_subjects_only.matrix2params_incremental())

            parameter_values_static = pdf_subjects_only.matrix2params_incremental()

            initial_guess = np.random.random(num_free_parameters_synonly)

            # pdf_subjects_syns_only should be the new object that fitness_func operates on, instead of
            # pdf_with_srvs_for_optimization
        else:
            # already like this, so simple renaming
            pdf_subjects_syns_only = pdf_with_srvs_for_optimization

            parameter_values_static = parameter_values_before

            num_free_parameters_synonly = num_free_parameters

            # subject_variables = range(len(self))

        upper_bound_synergistic_information = self.synergistic_entropy_upper_bound(subject_variables)
        if not agnostic_about is None:
            # upper_bound_agnostic_information is only used to normalize the cost term for non-zero MI with
            # the agnostic_variables (evidently a SRV is sought that has zero MI with these)
            if np.ndim(agnostic_about) == 1:
                upper_bound_agnostic_information = self.entropy(agnostic_about)
            elif np.ndim(agnostic_about) == 2:
                # in this case the cost term is the sum of MIs of the SRV with the sublists, so max cost term is this..
                upper_bound_agnostic_information = sum([self.entropy(ai) for ai in agnostic_about])
        else:
            upper_bound_agnostic_information = 0  # should not even be used...

        # todo: should lower the upper bounds if the max possible entropy of the SRVs is less...

        assert upper_bound_synergistic_information != 0.0, 'can never find any SRV!'

        # in the new self, these indices will identify the synergistic variables that will be added
        synergistic_variables = range(len(self), len(self) + num_synergistic_variables)

        def cost_func_subjects_only(free_params, parameter_values_before, extra_cost_rel_error=True):
            """
            This cost function searches for a Pr(S,Y,A,X) such that X is synergistic about S (subject_variables) only.
            This fitness function also taxes any correlation between X and A (agnostic_variables), but does not care
            about the relation between X and Y.
            :param free_params:
            :param parameter_values_before:
            :return:
            """
            assert len(free_params) == num_free_parameters_synonly

            if min(free_params) < -0.00001 or max(free_params) > 1.00001:
                warnings.warn('scipy\'s minimize() is violating the parameter bounds 0...1 I give it: '
                              + str(free_params))

                # high cost for invalid parameter values
                # note: maximum cost normally from this function is about 2.0
                return 10.0 + 100.0 * np.sum([p - 1.0 for p in free_params if p > 1.0]
                                             + [np.abs(p) for p in free_params if p < 0.0])

            # assert max(free_params) <= 1.00001, \
            #     'scipy\'s minimize() is violating the parameter bounds 0...1 I give it: ' + str(free_params)

            free_params = [min(max(fp, 0.0), 1.0) for fp in free_params]  # clip small roundoff errors

            pdf_subjects_syns_only.params2matrix_incremental(list(parameter_values_before) + list(free_params))

            # this if-block will add cost for the estimated amount of synergy induced by the proposed parameters,
            # and possible also a cost term for the ratio of synergistic versus non-synergistic info as extra
            if not subject_variables is None:
                assert pdf_subjects_syns_only.numvariables == len(subject_variables) + num_synergistic_variables

                # this can be considered to be in range [0,1] although particularly bad solutions can go >1
                if not extra_cost_rel_error:
                    cost = (upper_bound_synergistic_information - pdf_subjects_syns_only.synergistic_information_naive(
                        variables_SRV=range(len(subject_variables), len(pdf_subjects_syns_only)),
                        variables_X=range(len(subject_variables)))) / upper_bound_synergistic_information
                else:
                    tot_mi = pdf_subjects_syns_only.mutual_information(
                        range(len(subject_variables), len(pdf_subjects_syns_only)),
                        range(len(subject_variables)))

                    indiv_mis = [pdf_subjects_syns_only.mutual_information([var],
                                                                           range(len(subject_variables),
                                                                                 len(pdf_subjects_syns_only)))
                                 for var in range(len(subject_variables))]

                    syninfo_naive = tot_mi - sum(indiv_mis)

                    # this can be considered to be in range [0,1] although particularly bad solutions can go >1
                    cost = (upper_bound_synergistic_information - syninfo_naive) \
                           / upper_bound_synergistic_information

                    # add an extra cost term for the fraction of 'individual' information versus the total information
                    # this can be considered to be in range [0,1] although particularly bad solutions can go >1
                    if tot_mi != 0:
                        cost += sum(indiv_mis) / tot_mi
                    else:
                        cost += sum(indiv_mis)
            else:
                assert pdf_subjects_syns_only.numvariables == len(self) + num_synergistic_variables

                if not extra_cost_rel_error:
                    # this can be considered to be in range [0,1] although particularly bad solutions can go >1
                    cost = (upper_bound_synergistic_information - pdf_subjects_syns_only.synergistic_information_naive(
                        variables_SRV=range(len(self), len(pdf_subjects_syns_only)),
                        variables_X=range(len(self)))) / upper_bound_synergistic_information
                else:
                    tot_mi = pdf_subjects_syns_only.mutual_information(
                        range(len(self), len(pdf_subjects_syns_only)),
                        range(len(self)))

                    indiv_mis = [pdf_subjects_syns_only.mutual_information([var],
                                                                           range(len(self),
                                                                                 len(pdf_subjects_syns_only)))
                                 for var in range(len(self))]

                    syninfo_naive = tot_mi - sum(indiv_mis)

                    # this can be considered to be in range [0,1] although particularly bad solutions can go >1
                    cost = (upper_bound_synergistic_information - syninfo_naive) \
                           / upper_bound_synergistic_information

                    # add an extra cost term for the fraction of 'individual' information versus the total information
                    # this can be considered to be in range [0,1] although particularly bad solutions can go >1
                    if tot_mi != 0:
                        cost += sum(indiv_mis) / tot_mi
                    else:
                        cost += sum(indiv_mis)

            # this if-block will add a cost term for not being agnostic to given variables, usually (a) previous SRV(s)
            if not agnostic_about is None:
                assert not subject_variables is None, 'how can all variables be subject_variables and you still want' \
                                                      ' to be agnostic about certain (other) variables? (if you did' \
                                                      ' not specify subject_variables, do so.)'

                # make a conditional distribution of the synergistic variables conditioned on the subject variables
                # so that I can easily make a new joint pdf object with them and quantify this extra cost for the
                # agnostic constraint
                cond_pdf_syns_on_subjects = pdf_subjects_syns_only.conditional_probability_distributions(
                    range(len(subject_variables))
                )

                assert type(cond_pdf_syns_on_subjects) == dict
                assert len(cond_pdf_syns_on_subjects) > 0, 'conditioned on 0 values?'

                # if this hits then something is unintuitive with conditioning on variables...
                assert len(cond_pdf_syns_on_subjects.keys()[0]) == len(subject_variables)

                pdf_with_srvs_for_agnostic = self.copy()
                pdf_with_srvs_for_agnostic.append_variables_using_conditional_distributions(cond_pdf_syns_on_subjects,
                                                                                            subject_variables)

                if np.ndim(agnostic_about) == 1:
                    # note: cost term for agnostic is in [0,1]
                    cost += (pdf_with_srvs_for_agnostic.mutual_information(synergistic_variables, agnostic_about)) \
                            / upper_bound_agnostic_information
                else:
                    assert np.ndim(agnostic_about) == 2, 'expected list of lists, not more... made a mistake?'

                    for agn_i in agnostic_about:
                        # note: total cost term for agnostic is in [0,1]
                        cost += (1.0 / float(len(agnostic_about))) * \
                                pdf_with_srvs_for_agnostic.mutual_information(synergistic_variables, agn_i) \
                                / upper_bound_agnostic_information

            assert np.isscalar(cost)
            assert np.isfinite(cost)

            return float(cost)

        param_vectors_trace = []

        # these options are altered mainly to try to lower the computation time, which is considerable.
        minimize_options = {'ftol': 1e-6}

        if num_repeats == 1:
            optres = minimize(cost_func_subjects_only, initial_guess, bounds=[(0.0, 1.0)]*num_free_parameters_synonly,
                              callback=(lambda xv: param_vectors_trace.append(list(xv))) if verbose else None,
                              args=(parameter_values_static,), method=minimize_method, options=minimize_options)
        else:
            # optres_list = [minimize(cost_func_subjects_only, np.random.random(num_free_parameters_synonly),
            #                        bounds=[(0.0, 1.0)]*num_free_parameters_synonly,
            #                        callback=(lambda xv: param_vectors_trace.append(list(xv))) if verbose else None,
            #                        args=(parameter_values_static,)) for _ in xrange(num_repeats)]

            optres_list = []

            for _ in xrange(num_repeats):
                optres_ix = minimize(cost_func_subjects_only, np.random.random(num_free_parameters_synonly),
                                     bounds=[(0.0, 1.0)]*num_free_parameters_synonly,
                                     callback=(lambda xv: param_vectors_trace.append(list(xv))) if verbose else None,
                                     args=(parameter_values_static,), method=minimize_method, options=minimize_options)

                if verbose:
                    print 'note: finished a repeat. success=' + str(optres_ix.success) + ', cost=' \
                          + str(optres_ix.fun)

                optres_list.append(optres_ix)

            if verbose and __debug__:
                print 'debug: num_repeats=' + str(num_repeats) + ', all cost values were: ' \
                      + str([resi.fun for resi in optres_list])
                print 'debug: successes =', [resi.success for resi in optres_list]

            optres_list = [resi for resi in optres_list if resi.success]  # filter out the unsuccessful optimizations

            if len(optres_list) == 0:
                raise UserWarning('all ' + str(num_repeats) + ' optimizations using minimize() failed...?!')

            costvals = [res.fun for res in optres_list]
            min_cost = min(costvals)
            optres_ix = costvals.index(min_cost)

            assert optres_ix >= 0 and optres_ix < len(optres_list)

            optres = optres_list[optres_ix]

        if subject_variables is None:
            assert len(optres.x) == num_free_parameters
        else:
            assert len(optres.x) == num_free_parameters_synonly

        assert max(optres.x) <= 1.0000001, 'parameter bound significantly violated, ' + str(max(optres.x))
        assert min(optres.x) >= -0.0000001, 'parameter bound significantly violated, ' + str(min(optres.x))

        # todo: reuse the .append_optimized_variables (or so) instead, passing the cost function only? would also
        # validate that function.

        # clip to valid range
        optres.x = [min(max(xi, 0.0), 1.0) for xi in optres.x]

        # optimal_parameters_joint_pdf = list(parameter_values_before) + list(optres.x)
        # pdf_with_srvs.params2matrix_incremental(optimal_parameters_joint_pdf)

        assert len(pdf_subjects_syns_only.matrix2params_incremental()) == len(parameter_values_static) + len(optres.x)

        pdf_subjects_syns_only.params2matrix_incremental(list(parameter_values_static) + list(optres.x))

        if not subject_variables is None:
            cond_pdf_syns_on_subjects = pdf_subjects_syns_only.conditional_probability_distributions(
                range(len(subject_variables))
            )

            assert type(cond_pdf_syns_on_subjects) == dict
            assert len(cond_pdf_syns_on_subjects) > 0, 'conditioned on 0 values?'

            # if this hits then something is unintuitive with conditioning on variables...
            assert len(cond_pdf_syns_on_subjects.keys()[0]) == len(subject_variables)

            pdf_with_srvs = self.copy()
            pdf_with_srvs.append_variables_using_conditional_distributions(cond_pdf_syns_on_subjects,
                                                                           given_variables=subject_variables)
        else:
            pdf_with_srvs = pdf_subjects_syns_only  # all variables are subject

        assert pdf_with_srvs.numvariables == self.numvariables + num_synergistic_variables

        if verbose:
            parameter_values_after2 = pdf_with_srvs.matrix2params_incremental()

            assert len(parameter_values_after2) > len(parameter_values_before), 'should be additional free parameters'

            if not 1.0 in parameter_values_after2 and not 0.0 in parameter_values_after2:
                # see if the first part of the new parameters is exactly the old parameters, so that I know that I only
                # have to optimize the latter part of parameter_values_after
                np.testing.assert_array_almost_equal(parameter_values_before,
                                                     parameter_values_after2[:len(parameter_values_before)])
                np.testing.assert_array_almost_equal(parameter_values_after2[len(parameter_values_before):],
                                                     optres.x)
            else:
                # it can happen that some parameters are 'useless' in the sense that they defined conditional
                # probabilities in the case where the prior (that is conditioned upon) has zero probability. The
                # resulting pdf is then always the same, no matter this parameter value. This can only happen if
                # there is a 0 or 1 in the parameter list, (sufficient but not necessary) so skip the test then...
                pass

            # store the synergistic information before the optimization procedure (after procedure should be higher...)
            debug_after_syninfo = pdf_with_srvs.synergistic_information_naive(variables_SRV=range(self.numvariables,
                                                                                         pdf_with_srvs.numvariables),
                                                                              variables_X=range(self.numvariables))

            if verbose:
                print 'debug: append_synergistic_variables: I started from synergistic information =', \
                    debug_before_syninfo, 'at initial guess. After optimization it became', debug_after_syninfo, \
                    '(should be higher). Optimal params:', \
                    parameter_values_after2[len(parameter_values_before):]

            # print 'debug: optimization result:', optres

            # print 'debug: fitness of optimum (optres.x):', fitness_func(optres.x)
            # print 'debug: fitness of optimum (optres.x):', fitness_func(optres.x,
            #                                                             parameter_values_before=parameter_values_before)

            # print 'debug: parameter vector trace has', len(param_vectors_trace), 'vectors. First 20 and fitness:'
            # for xv in param_vectors_trace[:min(len(param_vectors_trace), 20)]:
            #     fitness = fitness_func(xv)
            #
            #     print 'debug: fitness =', fitness, '; vector =', xv

            # # fitness function should be the negative of the synergistic info right?
            # np.testing.assert_almost_equal(debug_after_syninfo, -fitness_func(optres.x))

        # if __debug__:
        #     # unchanged, not accidentally changed by passing it as reference? looking for bug
        #     np.testing.assert_array_almost_equal(debug_params_before, parameter_values_before)

        self.duplicate(pdf_with_srvs)


    def susceptibility_local_single(self, var_id, num_output_variables, perturbation_size=0.1, ntrials=25):

        if var_id >= len(self) - num_output_variables:
            raise ValueError('variable ' + str(var_id) + ' to perturb lies in output range, not input (first '
                             + str(len(self) - num_output_variables) + ' variables)')

        original_mi = self.mutual_information(range(len(self) - num_output_variables),
                                               range(len(self) - num_output_variables, len(self)))

        cond_pdf_rest = self.conditional_probability_distributions([var_id])

        pdf_var_only_orig = self.marginalize_distribution([var_id])
        params_orig = pdf_var_only_orig.matrix2params_incremental()

        def clip_to_unit_line(num):  # helper function, make sure all probabilities remain valid
            return max(min(num, 1), 0)

        # test if the procedure in the loop below is correct, namely that the new pdf gluued back together from pieces
        # is indeed the original <self>
        if __debug__ and np.random.random() < 0.01:  # make very infrequent, has worked for a while now
            # this test can be removed after it works for a while

            pdf_var_only = pdf_var_only_orig.copy()

            perturbation = np.array([0]*len(params_orig))  # zero perturbation

            new_params = map(clip_to_unit_line, params_orig + perturbation)

            np.testing.assert_array_almost_equal(new_params, params_orig)

            pdf_var_only.params2matrix_incremental(new_params)

            # note: if var_id > 0 then the ordering of the variables will have changed, putting the var_id to index 0
            pdf_var_only.append_variables_using_conditional_distributions(cond_pdf_rest)

            np.testing.assert_almost_equal(pdf_var_only.entropy(), self.entropy())

            new_mi = pdf_var_only.mutual_information(range(len(pdf_var_only) - num_output_variables),
                                                     range(len(pdf_var_only) - num_output_variables, len(pdf_var_only)))

            np.testing.assert_almost_equal(new_mi, original_mi)

        susceptibilities = []

        for i in xrange(ntrials):
            pdf_var_only = pdf_var_only_orig.copy()

            perturbation = np.random.random(len(params_orig))
            perturbation = perturbation / np.linalg.norm(perturbation) * perturbation_size  # normalize vector norm

            new_params = map(clip_to_unit_line, params_orig + perturbation)

            pdf_var_only.params2matrix_incremental(new_params)

            pdf_var_only.append_variables_using_conditional_distributions(cond_pdf_rest)

            new_mi = pdf_var_only.mutual_information(range(len(pdf_var_only) - num_output_variables),
                                                     range(len(pdf_var_only) - num_output_variables, len(pdf_var_only)))

            susceptibilities.append(abs(new_mi - original_mi))

        return np.mean(susceptibilities) / original_mi


    def susceptibilities_local(self, num_output_variables, perturbation_size=0.1, ntrials=25):
        return [self.susceptibility_local_single(varid, num_output_variables, perturbation_size=perturbation_size,
                                                 ntrials=ntrials) for varid in xrange(len(self) - num_output_variables)]


    def susceptibility_local(self, num_output_variables, perturbation_size=0.1, ntrials=25):
        return np.mean([self.susceptibility_local_single(varid, num_output_variables,
                                                         perturbation_size=perturbation_size,
                                                         ntrials=ntrials)
                        for varid in xrange(len(self) - num_output_variables)])


    def susceptibility_global(self, num_output_variables, perturbation_size=0.1, ntrials=25):
        """
        Perturb the current pdf Pr(X,Y) by changing Pr(X) slightly to Pr(X'), but keeping Pr(Y|X) fixed. Then
        measure the relative change in mutual information I(X:Y). Do this by Monte Carlo using <ntrials> repeats.
        :param num_output_variables: the number of variables at the end of the list of variables which are considered
        to be Y. The first variables are taken to be X. If your Y is mixed with the X, first do reordering.
        :param perturbation_size:
        :param ntrials:
        :return: expected relative change of mutual information I(X:Y) after perturbation
        :rtype: float
        """
        num_input_variables = len(self) - num_output_variables

        assert num_input_variables > 0, 'makes no sense to compute resilience with an empty set'

        original_mi = self.mutual_information(range(num_input_variables),
                                              range(num_input_variables, num_input_variables + num_output_variables))

        pdf_input_only = self.marginalize_distribution(range(num_input_variables))

        affected_params = pdf_input_only.matrix2params_incremental()  # perturb Pr(X) to Pr(X')
        static_params = self.matrix2params_incremental()[len(affected_params):]  # keep Pr(Y|X) fixed

        pdf_perturbed = self.copy()  # will be used to store the perturbed Pr(X)Pr(Y'|X)

        def clip_to_unit_line(num):  # helper function, make sure all probabilities remain valid
            return max(min(num, 1), 0)

        susceptibilities = []

        for i in xrange(ntrials):
            perturbation = np.random.random(len(affected_params))
            perturbation = perturbation / np.linalg.norm(perturbation) * perturbation_size  # normalize vector norm

            pdf_perturbed.params2matrix_incremental(map(clip_to_unit_line, np.add(affected_params, perturbation))
                                                    + static_params)

            susceptibility = pdf_perturbed.mutual_information(range(num_input_variables),
                                                              range(num_input_variables,
                                                                    num_input_variables + num_output_variables)) \
                             - original_mi
            susceptibilities.append(abs(susceptibility))

        return np.mean(susceptibilities) / original_mi


    def susceptibility_non_local(self, num_output_variables, perturbation_size=0.1, ntrials=25):
        """
        Perturb the current pdf Pr(X,Y) by changing Pr(Y|X) slightly, but keeping Pr(Y) and Pr(X) fixed. Then
        measure the relative change in mutual information I(X:Y). Do this by Monte Carlo using <ntrials> repeats.
        :param num_output_variables: the number of variables at the end of the list of variables which are considered
        to be Y. The first variables are taken to be X. If your Y is mixed with the X, first do reordering.
        :param perturbation_size:
        :param ntrials:
        :return: expected relative change of mutual information I(X:Y) after perturbation
        :rtype: float
        """
        num_input_variables = len(self) - num_output_variables

        assert num_input_variables > 0, 'makes no sense to compute resilience with an empty set'

        original_mi = self.mutual_information(range(num_input_variables),
                                              range(num_input_variables, num_input_variables + num_output_variables))

        susceptibilities = []

        for i in xrange(ntrials):
            resp = self.perturb_non_local(num_output_variables, perturbation_size)
            pdf_perturbed = resp.pdf

            susceptibility = pdf_perturbed.mutual_information(range(num_input_variables),
                                                              range(num_input_variables,
                                                                    num_input_variables + num_output_variables)) \
                             - original_mi
            susceptibilities.append(abs(susceptibility))

        return np.mean(susceptibilities) / original_mi


    class PerturbNonLocalResponse():
        pdf = None  # JointProbabilityMatrix object
        cost_same_output_marginal = None  # float, should be as close to zero as possible.
        cost_different_relation = None  # float, should be as close to zero as possible


    def perturb_non_local(self, num_output_variables, perturbation_size=0.1):
        """
        Perturb the pdf P(X,Y) by changing P(Y|X) without actually changing P(X) or P(Y), by numerical optimization.
        Y is assumed to be formed by stochastic variables collected at the end of the pdf.
        :param num_output_variables: |Y|.
        :param perturbation_size: the higher, the more different P(Y|X) will be from <self>
        :rtype: PerturbNonLocalResponse
        """
        num_input_variables = len(self) - num_output_variables

        assert num_input_variables > 0, 'makes no sense to perturb a relation with an empty set'

        original_params = self.matrix2params_incremental()

        static_params = list(self[range(num_input_variables)].matrix2params_incremental())

        num_free_params = len(original_params) - len(static_params)

        marginal_output_pdf = self[range(num_input_variables, len(self))]
        assert len(marginal_output_pdf) == num_output_variables, 'programming error'
        marginal_output_pdf_params = marginal_output_pdf.matrix2params_incremental()

        pdf_new = self.copy()  # just to create an object which I can replace everytime in cost function

        def clip_to_unit_line(num):  # helper function, make sure all probabilities remain valid
            return max(min(num, 1), 0)

        def cost_perturb_non_local(free_params, return_cost_list=False):
            new_params = static_params + list(map(clip_to_unit_line, free_params))

            pdf_new.params2matrix_incremental(new_params)

            marginal_output_pdf_new = pdf_new[range(num_input_variables, len(self))]
            marginal_output_pdf_new_params = marginal_output_pdf_new.matrix2params_incremental()

            cost_same_output_marginal = np.linalg.norm(np.subtract(marginal_output_pdf_new_params,
                                                                   marginal_output_pdf_params))
            cost_different_relation = np.linalg.norm(np.subtract(free_params, original_params[len(static_params):]))

            if not return_cost_list:
                cost = np.power(cost_same_output_marginal - 0.0, 2) \
                       + np.power(cost_different_relation - perturbation_size, 2)
            else:
                cost = [np.power(cost_same_output_marginal - 0.0, 2),
                        np.power(cost_different_relation - perturbation_size, 2)]

            return cost

        initial_guess_perturb_vec = np.random.random(num_free_params)
        initial_guess_perturb_vec /= np.linalg.norm(initial_guess_perturb_vec)
        initial_guess_perturb_vec *= perturbation_size

        initial_guess = np.add(original_params[len(static_params):],
                               initial_guess_perturb_vec)
        initial_guess = map(clip_to_unit_line, initial_guess)  # make sure stays in hypercube's unit volume

        optres = minimize(cost_perturb_non_local, initial_guess, bounds=[(0.0, 1.0)]*num_free_params)

        assert optres.success, 'scipy\'s minimize() failed'

        assert optres.fun >= -0.0001, 'cost function as constructed cannot be negative, what happened?'

        pdf_new.params2matrix_incremental(list(static_params) + list(optres.x))

        # get list of individual cost terms, for better diagnostics by caller
        cost = cost_perturb_non_local(optres.x, return_cost_list=True)

        resp = self.PerturbNonLocalResponse()
        resp.pdf = pdf_new  # final resulting pdf P'(X,Y), which is slightly different, perturbed version of <self>
        resp.cost_same_output_marginal = float(cost[0])  # cost of how different marginal P(Y) became (bad)
        resp.cost_different_relation = float(cost[1])  # cost of how different P(Y|X) is, compared to desired difference

        return resp


    def append_globally_resilient_variables(self, num_appended_variables, target_mi):

        # input_variables = [d for d in xrange(self.numvariables) if not d in output_variables]

        parameter_values_before = list(self.matrix2params_incremental())

        pdf_new = self.copy()
        pdf_new.append_variables(num_appended_variables)

        assert pdf_new.numvariables == self.numvariables + num_appended_variables

        parameter_values_after = pdf_new.matrix2params_incremental()

        # this many parameters (each in [0,1]) must be optimized
        num_free_parameters = len(parameter_values_after) - len(parameter_values_before)

        def cost_func_resilience_and_mi(free_params, parameter_values_before):
            pdf_new.params2matrix_incremental(list(parameter_values_before) + list(free_params))

            mi = pdf_new.mutual_information(range(len(self)), range(len(self), len(pdf_new)))

            susceptibility = pdf_new.susceptibility_global(num_appended_variables)

            return np.power(abs(target_mi - mi) / target_mi + susceptibility, 2)

        self.append_optimized_variables(num_appended_variables, cost_func=cost_func_resilience_and_mi,
                                        initial_guess=np.random.random(num_free_parameters))

        return


    def append_variables_with_target_mi(self, num_appended_variables, target_mi, verbose=False):

        # input_variables = [d for d in xrange(self.numvariables) if not d in output_variables]

        parameter_values_before = list(self.matrix2params_incremental())

        pdf_new = self.copy()
        pdf_new.append_variables(num_appended_variables)

        assert pdf_new.numvariables == self.numvariables + num_appended_variables

        parameter_values_after = pdf_new.matrix2params_incremental()

        # this many parameters (each in [0,1]) must be optimized
        num_free_parameters = len(parameter_values_after) - len(parameter_values_before)

        def cost_func_target_mi(free_params, parameter_values_before):
            pdf_new.params2matrix_incremental(list(parameter_values_before) + list(free_params))

            mi = pdf_new.mutual_information(range(len(self)), range(len(self), len(pdf_new)))

            return np.power(abs(target_mi - mi) / target_mi, 2)

        self.append_optimized_variables(num_appended_variables, cost_func=cost_func_target_mi,
                                        initial_guess=np.random.random(num_free_parameters),
                                        verbose=verbose)

        return


    def append_optimized_variables(self, num_appended_variables, cost_func, initial_guess=None, verbose=True):
        """
        Append variables in such a way that their conditional pdf with the existing variables is optimized in some
        sense, for instance they can be synergistic (append_synergistic_variables) or orthogonalized
        (append_orthogonalized_variables). Use the cost_func to determine the relation between the new appended
        variables and the pre-existing ones.
        :param num_appended_variables:
        :param cost_func: a function cost_func(free_params, parameter_values_before) which returns a float.
        The parameter set list(parameter_values_before) + list(free_params) defines a joint pdf of the appended
        variables together with the pre-existing ones, and free_params by itself defines completely the conditional
        pdf of the new variables given the previous. Use params2matrix_incremental to construct a joint pdf from the
        parameters and evaluate whatever you need, and then return a float. The higher the return value of cost_func
        the more desirable the joint pdf induced by the parameter set list(parameter_values_before) + list(free_params).
        :param initial_guess: initial guess for 'free_params' where you think cost_func will return a relatively
        low value. It can also be None, in which case a random point in parameter space will be chosen. It can also
        be an integer value like 10, in which case 10 optimizations will be run each starting from a random point
        in parameter space, and the best solution is selected.
        :param verbose:
        :rtype: scipy.optimize.OptimizeResult
        """

        # these parameters should be unchanged and the first set of parameters of the resulting pdf_new
        parameter_values_before = list(self.matrix2params_incremental())

        if __debug__:
            debug_params_before = copy.deepcopy(parameter_values_before)

        # a pdf with XORs as appended variables (often already MSRV for binary variables), good initial guess?
        # note: does not really matter how I set the pdf of this new pddf, as long as it has the correct number of
        # paarameters for optimization below
        pdf_new = self.copy()
        pdf_new.append_variables_using_state_transitions_table(
            state_transitions=lambda vals, mv: [int(np.mod(np.sum(vals), mv))]*num_appended_variables)

        assert pdf_new.numvariables == self.numvariables + num_appended_variables

        parameter_values_after = pdf_new.matrix2params_incremental()

        assert num_appended_variables > 0, 'makes no sense to add 0 variables'
        assert len(parameter_values_after) > len(parameter_values_before), 'should be >0 free parameters to optimize?'
        # see if the first part of the new parameters is exactly the old parameters, so that I know that I only
        # have to optimize the latter part of parameter_values_after
        np.testing.assert_array_almost_equal(parameter_values_before,
                                             parameter_values_after[:len(parameter_values_before)])

        # this many parameters (each in [0,1]) must be optimized
        num_free_parameters = len(parameter_values_after) - len(parameter_values_before)

        assert num_appended_variables == 0 or num_free_parameters > 0

        # if initial_guess is None:
        #     initial_guess = np.random.random(num_free_parameters)  # start from random point in parameter space

        param_vectors_trace = []  # storing the parameter vectors visited by the minimize() function

        if type(initial_guess) == int:
            num_repeats = int(initial_guess)
            initial_guess = None

            assert num_repeats > 0, 'makes no sense to optimize zero times?'
        else:
            num_repeats = 1

        optres = None

        for rep in xrange(num_repeats):
            if initial_guess is None:
                initial_guess_i = np.random.random(num_free_parameters)  # start from random point in parameter space
            else:
                initial_guess_i = initial_guess  # always start from supplied point in parameter space

            assert len(initial_guess_i) == num_free_parameters

            if verbose:
                print 'debug: starting minimize() #' + str(rep) \
                      + ' at params=' + str(initial_guess_i) + ' at cost_func=' \
                      + str(cost_func(initial_guess_i, parameter_values_before))

            optres_i = minimize(cost_func,
                              initial_guess_i, bounds=[(0.0, 1.0)]*num_free_parameters,
                              # callback=(lambda xv: param_vectors_trace.append(list(xv))) if verbose else None,
                              args=(parameter_values_before,))

            if optres_i.success:
                if verbose:
                    print 'debug: successfully ended minimize() #' + str(rep) \
                          + ' at params=' + str(optres_i.x) + ' at cost_func=' \
                          + str(optres_i.fun)

                if optres is None:
                    optres = optres_i
                elif optres.fun > optres_i.fun:
                    optres = optres_i
                else:
                    pass  # this solution is worse than before, so do not change optres

        if optres is None:
            # could never find a good solution, in all <num_repeats> attempts
            raise UserWarning('always failed to successfully optimize: increase num_repeats')

        assert len(optres.x) == num_free_parameters
        assert max(optres.x) <= 1.0001, 'parameter bound significantly violated: ' + str(optres.x)
        assert min(optres.x) >= -0.0001, 'parameter bound significantly violated: ' + str(optres.x)

        # clip the parameters within the allowed bounds
        optres.x = [min(max(xi, 0.0), 1.0) for xi in optres.x]

        optimal_parameters_joint_pdf = list(parameter_values_before) + list(optres.x)

        pdf_new.params2matrix_incremental(optimal_parameters_joint_pdf)

        assert len(pdf_new) == len(self) + num_appended_variables

        if __debug__:
            parameter_values_after2 = pdf_new.matrix2params_incremental()

            assert len(parameter_values_after2) > len(parameter_values_before), 'should be additional free parameters'
            # see if the first part of the new parameters is exactly the old parameters, so that I know that I only
            # have to optimize the latter part of parameter_values_after
            np.testing.assert_array_almost_equal(parameter_values_before,
                                                 parameter_values_after2[:len(parameter_values_before)])
            # note: for the if see the story in params2matrix_incremental()
            if not (0.0 in self.scalars_up_to_level(parameter_values_after2) or \
                            1.0 in self.scalars_up_to_level(parameter_values_after2)):
                try:
                    np.testing.assert_array_almost_equal(parameter_values_after2[len(parameter_values_before):],
                                                         optres.x)
                except AssertionError as e:
                    # are they supposed to be equal, but in different order?
                    print 'debug: sum params after 1 =', np.sum(parameter_values_after2[len(parameter_values_before):])
                    print 'debug: sum params after 2 =', optres.x
                    print 'debug: parameter_values_before (which IS equal and correct) =', parameter_values_before
                    # does this one below have a 1 or 0 in it? because then the error could be caused by the story in
                    # params2matrix_incremental()
                    print 'debug: parameter_values_after2 =', parameter_values_after2

                    raise AssertionError(e)
        if __debug__:
            # unchanged, not accidentally changed by passing it as reference? looking for bug
            np.testing.assert_array_almost_equal(debug_params_before, parameter_values_before)

        self.duplicate(pdf_new)

        return optres


    def append_orthogonalized_variables(self, variables_to_orthogonalize, num_added_variables_orthogonal,
                                        num_added_variables_parallel, verbose=True,
                                        num_repeats=1, randomization_per_repeat=0.01):

        """
        Let X=<variables_to_orthogonalize> and Y=complement[<variables_to_orthogonalize>]. Add two sets of variables
        X1 and X2 such that I(X1:Y)=0, I(X1:X)=H(X|Y); and I(X2:Y)=I(X2:X)=I(X:Y), and I(X1:X2)=0. In words, X is being
        decomposed into two parts: X1 (orthogonal to Y, MI=0) and X2 (parallel to Y, MI=max).

        This object itself will be expanded by <num_added_variables_orthogonal> + <num_added_variables_parallel>
        variables.

        Warning: if this pdf object also contains other variables Z which the orthogonalization can ignore, then
        first take them out of the pdf (marginalize or condition) because it blows up the number of free parameters
        that this function must optimize.

        (Do the optimization jointly, so both parallel and orthogonal together, because the function below seems not so
        effective for some reason.)
        :param variables_to_orthogonalize: list of variables which should be decomposed.
        :type variables_to_orthogonalize: list of int
        :param num_added_variables_orthogonal: |X1|
        :param num_added_variables_parallel: |X2|
        :param verbose:
        :param num_repeats: number of times to perform minimize() starting from random parameters
        and take the best result.
        :param randomization_per_repeat:
        :return: :raise ValueError:
        """
        assert num_added_variables_orthogonal > 0, 'cannot make parallel variables if I don\'t first make orthogonal'
        assert num_added_variables_parallel > 0, 'this may not be necessary, assert can be removed at some point?'

        # remove potential duplicates; and sort (for no reason)
        variables_to_orthogonalize = sorted(list(set(variables_to_orthogonalize)))

        pdf_ortho_para = self.copy()

        # these parameters should remain unchanged during the optimization, these are the variables (Y, X) where
        # X is the variable to be orthogonalized into the added (X1, X2) and Y is to be orthogonalized against.
        original_parameters = range(len(pdf_ortho_para.matrix2params_incremental(return_flattened=True)))

        original_variables = range(len(self))  # this includes variables_to_orthogonalize
        subject_variables = list(np.setdiff1d(original_variables, variables_to_orthogonalize))

        # the rest of the function assumes that the variables_to_orthogonalize are all at the end of the
        # original_variables, so that appending a partial conditional pdf (conditioned on
        # len(variables_to_orthogonalize) variables) will be conditioned on the variables_to_orthogonalize only. So
        # reorder the variables and recurse to this function once.
        if not max(subject_variables) < min(variables_to_orthogonalize):
            pdf_reordered = self.copy()

            pdf_reordered.reorder_variables(subject_variables + variables_to_orthogonalize)

            assert len(pdf_reordered) == len(self)
            assert len(range(len(subject_variables), pdf_reordered.numvariables)) == len(variables_to_orthogonalize)

            # perform orthogonalization on the ordered pdf, which it is implemented for
            pdf_reordered.append_orthogonalized_variables(range(len(subject_variables), len(self)))

            # find the mapping from ordered to disordered, so that I can reverse the reordering above
            new_order = subject_variables + variables_to_orthogonalize
            original_order = [-1] * len(new_order)
            for varix in xrange(len(new_order)):
                original_order[new_order[varix]] = varix

            assert not -1 in original_order
            assert max(original_order) < len(original_order), 'there were only so many original variables...'

            pdf_reordered.reorder_variables(original_order + range(len(original_order), len(pdf_reordered)))

            self.duplicate(pdf_reordered)

            # due to the scipy's minimize() procedure the list of parameters can be temporarily slightly out of bounds, like
            # 1.000000001, but soon after this should be clipped to the allowed range and e.g. at this point the parameters
            # are expected to be all valid
            assert min(self.matrix2params_incremental()) >= 0.0, \
                'parameter(s) out of bound: ' + str(self.matrix2params_incremental())
            assert max(self.matrix2params_incremental()) <= 1.0, \
                'parameter(s) out of bound: ' + str(self.matrix2params_incremental())

            return

        assert len(np.intersect1d(original_variables, variables_to_orthogonalize)) == len(variables_to_orthogonalize), \
            'original_variables should include the variables_to_orthogonalize'

        pdf_ortho_para.append_variables(num_added_variables_orthogonal)
        orthogonal_variables = range(len(self), len(pdf_ortho_para))

        orthogonal_parameters = range(len(original_parameters),
                                      len(pdf_ortho_para.matrix2params_incremental(return_flattened=True)))

        pdf_ortho_para.append_variables(num_added_variables_parallel)
        parallel_variables = range(len(pdf_ortho_para) - num_added_variables_parallel, len(pdf_ortho_para))

        parallel_parameters = range(len(orthogonal_parameters) + len(original_parameters),
                                    len(pdf_ortho_para.matrix2params_incremental(return_flattened=True)))

        assert len(np.intersect1d(orthogonal_parameters, parallel_parameters)) == 0, \
            'error: orthogonal_parameters = ' + str(orthogonal_parameters) + ', parallel_parameters = ' \
            + str(parallel_parameters)
        assert len(np.intersect1d(orthogonal_variables, parallel_variables)) == 0

        free_parameters = list(orthogonal_parameters) + list(parallel_parameters)

        initial_params_list = np.array(list(pdf_ortho_para.matrix2params_incremental(return_flattened=True)))

        # todo: make this more efficient by having only free parameters X1,X2 as p(X1,X2|X), not also conditioned on Y?
        # this would make it a two-step thing maybe again naively, but maybe there is a way to still do it simultaneous
        # (of course the cost function DOES depend on Y, but X1,X2=f(X) only and the point is to find an optimal f(),)

        # todo: START of optimization

        # let the X1,X2 parameterization only depend on X (=variables_to_orthogonalize) to reduce the parameter space
        # greatly
        pdf_X = self.marginalize_distribution(variables_to_orthogonalize)

        pdf_X_X1_X2 = pdf_X.copy()
        pdf_X_X1_X2.append_variables(num_added_variables_orthogonal)
        pdf_X_X1_X2.append_variables(num_added_variables_parallel)

        free_params_X1_X2_given_X = range(len(pdf_X.matrix2params_incremental()),
                                          len(pdf_X_X1_X2.matrix2params_incremental()))
        # parameter values at these parameter indices should not change:
        static_params_X = range(len(pdf_X.matrix2params_incremental()))
        static_params_X_values = list(np.array(pdf_X.matrix2params_incremental())[static_params_X])

        # check if building a complete joint pdf using the subset-conditional pdf works as expected
        if __debug__:
            cond_pdf_X1_X2_given_X = pdf_X_X1_X2.conditional_probability_distributions(range(len(variables_to_orthogonalize)))

            assert len(cond_pdf_X1_X2_given_X.values()[0]) == num_added_variables_orthogonal \
                                                              + num_added_variables_parallel

            pdf_test_Y_X_X1_X2 = self.copy()
            pdf_test_Y_X_X1_X2.append_variables_using_conditional_distributions(cond_pdf_X1_X2_given_X)

            # test if first so-many params are the same as pdf_ortho_para's
            np.testing.assert_array_almost_equal(initial_params_list[:len(original_parameters)],
                                                 pdf_test_Y_X_X1_X2.matrix2params_incremental()[:len(original_parameters)])
            # still, it would be immensely unlikely that the two pdfs are the same, since the conditional pdf(X1,X2|X)
            # is independently and randomly generated for both pdfs
            assert pdf_test_Y_X_X1_X2 != pdf_ortho_para

        # used repeatedly in the cost function below, prevent recomputing it every time
        ideal_H_X1 = self.conditional_entropy(variables_to_orthogonalize, subject_variables)
        ideal_mi_X2_Y = self.mutual_information(variables_to_orthogonalize, subject_variables)

        default_weights = (1, 1, 1, 1, 1, 1)  # used for cost_func_minimal_params and to know the number of weights


        # cost function used in scipy's minimize() procedure
        def cost_func_minimal_params(proposed_params, rel_weights=default_weights):
            assert len(proposed_params) == len(free_params_X1_X2_given_X)

            # relative weight coefficients for the different terms that contribute to the cost function below.
            wIso, wIvo, wIvp, wIsp, wHop, wIop = map(abs, rel_weights)

            pdf_X_X1_X2.params2matrix_incremental(static_params_X_values + list(proposed_params))

            cond_pdf_X1_X2_given_X = pdf_X_X1_X2.conditional_probability_distributions(range(len(variables_to_orthogonalize)))
            pdf_proposed_Y_X_X1_X2 = self.copy()
            pdf_proposed_Y_X_X1_X2.append_variables_using_conditional_distributions(cond_pdf_X1_X2_given_X)

            # test if first so-many params are the same as pdf_ortho_para's
            np.testing.assert_array_almost_equal(initial_params_list[:len(original_parameters)],
                                                 pdf_proposed_Y_X_X1_X2.matrix2params_incremental()[:len(original_parameters)])

            # note: unwanted effects should be positive terms ('cost'), and desired MIs should be negative terms
            # note: if you change the cost function here then you should also change the optimal value below
            # cost = wIso * pdf_proposed_Y_X_X1_X2.mutual_information(subject_variables, orthogonal_variables) \
            #        - wIsp * pdf_proposed_Y_X_X1_X2.mutual_information(subject_variables, parallel_variables) \
            #        - wHop * pdf_proposed_Y_X_X1_X2.entropy(orthogonal_variables) \
            #        + wIop * pdf_proposed_Y_X_X1_X2.mutual_information(orthogonal_variables, parallel_variables)

            # - wHop * pdf_proposed_Y_X_X1_X2.entropy(orthogonal_variables + parallel_variables) \

            # note: each term is intended to be normalized to [0, 1], where 0 is best and 1 is worst. Violation of this
            # is possible though, but they are really bad solutions.

            cost_terms = dict()
            cost_terms['Iso'] = \
                wIso * abs(pdf_proposed_Y_X_X1_X2.mutual_information(subject_variables, orthogonal_variables) - 0.0) \
                / ideal_mi_X2_Y
            cost_terms['Ivo'] = \
                wIvo * abs(pdf_proposed_Y_X_X1_X2.mutual_information(variables_to_orthogonalize, orthogonal_variables)
                           - ideal_H_X1) / ideal_H_X1
            # question: is the following one necessary?
            # cost_terms['Ivp'] = \
            #     wIvp * abs(pdf_proposed_Y_X_X1_X2.mutual_information(variables_to_orthogonalize, parallel_variables)
            #                - ideal_mi_X2_Y) / ideal_mi_X2_Y
            cost_terms['Ivp'] = 0.0
            cost_terms['Isp'] = \
                wIsp * abs(pdf_proposed_Y_X_X1_X2.mutual_information(subject_variables, parallel_variables)
                           - ideal_mi_X2_Y) / ideal_mi_X2_Y
            # cost_terms['Hop'] = wHop * abs(pdf_proposed_Y_X_X1_X2.entropy(orthogonal_variables) - ideal_H_X1)
            cost_terms['Hop'] = \
                wHop * 0.0  # let's see what happens now (trying to improve finding global optimum)
            cost_terms['Iop'] = \
                wIop * abs(pdf_proposed_Y_X_X1_X2.mutual_information(orthogonal_variables, parallel_variables) - 0.0) \
                / ideal_H_X1

            # sum of squared errors, or norm of vector in error space, to make a faster convergence hopefully
            cost = float(np.sum(np.power(cost_terms.values(), 2)))

            assert np.isfinite(cost)
            assert np.isscalar(cost)

            return float(cost)


        if __debug__:
            # for each term in the cost function above I determine what would be the optimal value,
            # note: if you change the cost function above then you should also change this optimal value
            debug_optimal_cost_value = 0

            # hack, I should actually perform some min()'s in this but this is just a rough guide:
            # note: assuming orthogonal_variables = variables_to_orthogonalize and H(parallel_variables) = 0
            debug_worst_cost_value = ideal_mi_X2_Y + abs(self.entropy(variables_to_orthogonalize) - ideal_H_X1) \
                                     + ideal_mi_X2_Y + ideal_mi_X2_Y \
                                     + abs(self.entropy(variables_to_orthogonalize) - ideal_H_X1)
                                                   # + ideal_mi_X2_Y \
                               # + 2 * abs(self.entropy(variables_to_orthogonalize) - ideal_H_X1) + ideal_mi_X2_Y

            debug_random_time_before = time.time()

            debug_num_random_costs = 20
            debug_random_cost_values = [cost_func_minimal_params(np.random.random(len(free_params_X1_X2_given_X)))
                                  for _ in xrange(debug_num_random_costs)]

            debug_random_time_after = time.time()

            # for trying to eye-ball whether the minimize() is doing significantly better than just random sampling,
            # getting a feel for the solution space
            debug_avg_random_cost_val = np.mean(debug_random_cost_values)
            debug_std_random_cost_val = np.std(debug_random_cost_values)
            debug_min_random_cost_val = np.min(debug_random_cost_values)
            debug_max_random_cost_val = np.max(debug_random_cost_values)

            print 'debug: cost values of random parameter vectors:', debug_random_cost_values, '-- took', \
                debug_random_time_after - debug_random_time_before, 'seconds for', debug_num_random_costs, 'vectors.'

        initial_guess_minimal = list(np.array(pdf_X_X1_X2.matrix2params_incremental())[free_params_X1_X2_given_X])

        # todo: bring optimal_cost_value to release code (not __debug__) and use it as criterion to decide if the
        # solution is acceptable? Maybe pass a tolerance value in which case I would raise an exception if exceeded?

        if verbose and __debug__:
            print 'debug: append_orthogonalized_variables: BEFORE optimization, cost func =', \
                cost_func_minimal_params(initial_guess_minimal), '(optimal=' + str(debug_optimal_cost_value) \
                + ', a worst=' + str(debug_worst_cost_value) + ', avg~=' + str(debug_avg_random_cost_val) \
                                                                 + '+-' + str(debug_std_random_cost_val) + ')'

            # mutual informations
            print 'debug: pdf_ortho_para.mutual_information(subject_variables, orthogonal_variables) =', \
                pdf_ortho_para.mutual_information(subject_variables, orthogonal_variables)
            print 'debug: pdf_ortho_para.mutual_information(subject_variables, parallel_variables) =', \
                pdf_ortho_para.mutual_information(subject_variables, parallel_variables)
            print 'debug: pdf_ortho_para.mutual_information(orthogonal_variables, parallel_variables) =', \
                pdf_ortho_para.mutual_information(orthogonal_variables, parallel_variables)
            print 'debug: pdf_ortho_para.mutual_information(variables_to_orthogonalize, orthogonal_variables) =', \
                pdf_ortho_para.mutual_information(variables_to_orthogonalize, orthogonal_variables)
            print 'debug: pdf_ortho_para.mutual_information(variables_to_orthogonalize, parallel_variables) =', \
                pdf_ortho_para.mutual_information(variables_to_orthogonalize, parallel_variables)
            print 'debug: pdf_ortho_para.mutual_information(variables_to_orthogonalize, ' \
                  'orthogonal_variables + parallel_variables) =', \
                pdf_ortho_para.mutual_information(variables_to_orthogonalize, orthogonal_variables + parallel_variables)
            print 'debug: pdf_ortho_para.mutual_information(variables_to_orthogonalize, ' \
                  'subject_variables) =', \
                pdf_ortho_para.mutual_information(variables_to_orthogonalize, subject_variables)

            # entropies
            print 'debug: pdf_ortho_para.entropy(subject_variables) =', \
                pdf_ortho_para.entropy(subject_variables)
            print 'debug: pdf_ortho_para.entropy(variables_to_orthogonalize) =', \
                pdf_ortho_para.entropy(variables_to_orthogonalize)
            print 'debug: pdf_ortho_para.entropy(orthogonal_variables) =', \
                pdf_ortho_para.entropy(orthogonal_variables)
            print 'debug: pdf_ortho_para.entropy(parallel_variables) =', \
                pdf_ortho_para.entropy(parallel_variables)

            print 'debug: (num free parameters for optimization:', len(free_params_X1_X2_given_X), ')'

            time_before = time.time()

        if num_repeats == 1:
            optres = minimize(cost_func_minimal_params, initial_guess_minimal,
                              bounds=[(0.0, 1.0)]*len(free_params_X1_X2_given_X),
                              args=((1,)*len(default_weights),))
        elif num_repeats > 1:
            # perform the minimization <num_repeats> times starting from random points in parameter space and select
            # the best one (lowest cost function value)

            # note: the args= argument passes the relative weights of the different (now 4) terms in the cost function
            # defined above. At every next iteration it allows more and more randomization around the value (1,1,1,1)
            # which means that every term would be equally important.
            optres_list = [minimize(cost_func_minimal_params, np.random.random(len(initial_guess_minimal)),
                                    bounds=[(0.0, 1.0)]*len(free_params_X1_X2_given_X),
                                    args=(tuple(1.0 + np.random.randn(len(default_weights))
                                                * randomization_per_repeat * repi),))
                           for repi in xrange(num_repeats)]

            if verbose and __debug__:
                print 'debug: num_repeats=' + str(num_repeats) + ', all cost values were: ' \
                      + str([resi.fun for resi in optres_list])
                print 'debug: successes =', [resi.success for resi in optres_list]

            optres_list = [resi for resi in optres_list if resi.success]  # filter out the unsuccessful optimizations

            assert len(optres_list) > 0, 'all ' + str(num_repeats) + ' optimizations using minimize() failed...?!'

            costvals = [res.fun for res in optres_list]
            min_cost = min(costvals)
            optres_ix = costvals.index(min_cost)

            assert optres_ix >= 0 and optres_ix < len(optres_list)

            optres = optres_list[optres_ix]
        else:
            raise ValueError('cannot repeat negative times')

        assert optres.success, 'scipy\'s minimize() failed'

        assert optres.fun >= -0.0001, 'cost function as constructed cannot be negative, what happened?'

        # build the most optimal PDF then finally:
        pdf_X_X1_X2.params2matrix_incremental(static_params_X_values + list(optres.x))
        cond_pdf_X1_X2_given_X = pdf_X_X1_X2.conditional_probability_distributions(range(len(variables_to_orthogonalize)))
        pdf_ortho_para = self.copy()
        pdf_ortho_para.append_variables_using_conditional_distributions(cond_pdf_X1_X2_given_X)

        # test if first so-many params are the same as pdf_ortho_para's
        np.testing.assert_array_almost_equal(initial_params_list[:len(original_parameters)],
                                             pdf_test_Y_X_X1_X2.matrix2params_incremental()[:len(original_parameters)])

        if verbose and __debug__:
            print 'debug: append_orthogonalized_variables: AFTER optimization, cost func =', optres.fun

            # mutual informations
            print 'debug: pdf_ortho_para.mutual_information(subject_variables, orthogonal_variables) =', \
                pdf_ortho_para.mutual_information(subject_variables, orthogonal_variables), \
                '(optimal=' + str(debug_optimal_cost_value) + ')'
            print 'debug: pdf_ortho_para.mutual_information(subject_variables, parallel_variables) =', \
                pdf_ortho_para.mutual_information(subject_variables, parallel_variables)
            print 'debug: pdf_ortho_para.mutual_information(orthogonal_variables, parallel_variables) =', \
                pdf_ortho_para.mutual_information(orthogonal_variables, parallel_variables)
            print 'debug: pdf_ortho_para.mutual_information(variables_to_orthogonalize, orthogonal_variables) =', \
                pdf_ortho_para.mutual_information(variables_to_orthogonalize, orthogonal_variables)
            print 'debug: pdf_ortho_para.mutual_information(variables_to_orthogonalize, parallel_variables) =', \
                pdf_ortho_para.mutual_information(variables_to_orthogonalize, parallel_variables)
            print 'debug: pdf_ortho_para.mutual_information(variables_to_orthogonalize, ' \
                  'orthogonal_variables + parallel_variables) =', \
                pdf_ortho_para.mutual_information(variables_to_orthogonalize, orthogonal_variables + parallel_variables)
            print 'debug: pdf_ortho_para.mutual_information(variables_to_orthogonalize, ' \
                  'subject_variables) =', \
                pdf_ortho_para.mutual_information(variables_to_orthogonalize, subject_variables)

            # entropies
            print 'debug: pdf_ortho_para.entropy(subject_variables) =', \
                pdf_ortho_para.entropy(subject_variables)
            print 'debug: pdf_ortho_para.entropy(variables_to_orthogonalize) =', \
                pdf_ortho_para.entropy(variables_to_orthogonalize)
            print 'debug: pdf_ortho_para.entropy(orthogonal_variables) =', \
                pdf_ortho_para.entropy(orthogonal_variables)
            print 'debug: pdf_ortho_para.entropy(parallel_variables) =', \
                pdf_ortho_para.entropy(parallel_variables)

            time_after = time.time()

            print 'debug: the optimization took', time_after - time_before, 'seconds in total.'

        self.duplicate(pdf_ortho_para)

        return

        # todo: END of optimization

        # note: the code below is 'old', in the sense that it should be equivalent above but using an optimization
        # step in a (much) larger parameter space

        def cost_func2(proposed_params):
            assert len(free_parameters) == len(proposed_params)

            params_list = list(initial_params_list[original_parameters]) + list(proposed_params)

            pdf_ortho_para.params2matrix_incremental(params_list)  # in-place

            # note: unwanted effects should be positive terms ('cost'), and desired MIs should be negative terms
            cost = pdf_ortho_para.mutual_information(subject_variables, orthogonal_variables) \
                   - pdf_ortho_para.mutual_information(subject_variables, parallel_variables) \
                   - pdf_ortho_para.mutual_information(variables_to_orthogonalize,
                                                       orthogonal_variables + parallel_variables)

            assert np.isfinite(cost)
            assert np.isscalar(cost)

            return float(cost)


        initial_guess = list(initial_params_list[free_parameters])

        if verbose and __debug__:
            print 'debug: append_orthogonalized_variables: BEFORE optimization, cost func =', cost_func2(initial_guess)

            # mutual informations
            print 'debug: pdf_ortho_para.mutual_information(subject_variables, orthogonal_variables) =', \
                pdf_ortho_para.mutual_information(subject_variables, orthogonal_variables)
            print 'debug: pdf_ortho_para.mutual_information(subject_variables, parallel_variables) =', \
                pdf_ortho_para.mutual_information(subject_variables, parallel_variables)
            print 'debug: pdf_ortho_para.mutual_information(orthogonal_variables, parallel_variables) =', \
                pdf_ortho_para.mutual_information(orthogonal_variables, parallel_variables)
            print 'debug: pdf_ortho_para.mutual_information(variables_to_orthogonalize, ' \
                  'orthogonal_variables + parallel_variables) =', \
                pdf_ortho_para.mutual_information(variables_to_orthogonalize, orthogonal_variables + parallel_variables)

            # entropies
            print 'debug: pdf_ortho_para.entropy(subject_variables) =', \
                pdf_ortho_para.entropy(subject_variables)
            print 'debug: pdf_ortho_para.entropy(variables_to_orthogonalize) =', \
                pdf_ortho_para.entropy(variables_to_orthogonalize)
            print 'debug: pdf_ortho_para.entropy(orthogonal_variables) =', \
                pdf_ortho_para.entropy(orthogonal_variables)
            print 'debug: pdf_ortho_para.entropy(parallel_variables) =', \
                pdf_ortho_para.entropy(parallel_variables)

        optres = minimize(cost_func2, initial_guess, bounds=[(0.0, 1.0)]*len(free_parameters))

        assert optres.success, 'scipy\'s minimize() failed'

        if verbose and __debug__:
            print 'debug: append_orthogonalized_variables: AFTER optimization, cost func =', optres.fun

            # mutual informations
            print 'debug: pdf_ortho_para.mutual_information(subject_variables, orthogonal_variables) =', \
                pdf_ortho_para.mutual_information(subject_variables, orthogonal_variables)
            print 'debug: pdf_ortho_para.mutual_information(subject_variables, parallel_variables) =', \
                pdf_ortho_para.mutual_information(subject_variables, parallel_variables)
            print 'debug: pdf_ortho_para.mutual_information(orthogonal_variables, parallel_variables) =', \
                pdf_ortho_para.mutual_information(orthogonal_variables, parallel_variables)
            print 'debug: pdf_ortho_para.mutual_information(variables_to_orthogonalize, ' \
                  'orthogonal_variables + parallel_variables) =', \
                pdf_ortho_para.mutual_information(variables_to_orthogonalize, orthogonal_variables + parallel_variables)

            # entropies
            print 'debug: pdf_ortho_para.entropy(subject_variables) =', \
                pdf_ortho_para.entropy(subject_variables)
            print 'debug: pdf_ortho_para.entropy(variables_to_orthogonalize) =', \
                pdf_ortho_para.entropy(variables_to_orthogonalize)
            print 'debug: pdf_ortho_para.entropy(orthogonal_variables) =', \
                pdf_ortho_para.entropy(orthogonal_variables)
            print 'debug: pdf_ortho_para.entropy(parallel_variables) =', \
                pdf_ortho_para.entropy(parallel_variables)

        self.duplicate(pdf_ortho_para)

        # del pdf_ortho_para


    # todo: set default entropy_cost_factor=0.1 or 0.05?
    # def append_orthogonalized_variables(self, variables_to_orthogonalize, num_added_variables_orthogonal=None,
    #                                     num_added_variables_parallel=None, entropy_cost_factor=0.1, verbose=True,
    #                                     num_repeats=3):
    #     """
    #     Orthogonalize the given set of variables_to_orthogonalize=X relative to the rest. I.e.,
    #     decompose X into two parts {X1,X2}
    #     :param num_added_variables_parallel:
    #     :param verbose:
    #     :param num_repeats: number of times that the optimization procedure of the cost function (for both
    #     the orthogonal variables and the parallel variables) is repeated, of which the best solution is then chosen.
    #     where I(X1:rest)=0 but I(X1:X)=H(X1) and I(X2:rest)=H(X2). The orthogonal set is added first and the parallel
    #     set last.
    #
    #     In theory the entropy of self as a whole should not increase, though it is not the end of the world if it does.
    #     But you can set e.g. entropy_cost_factor=0.1 to try and make it increase as little as possible.
    #     :param variables_to_orthogonalize: set of variables to orthogonalize (X)
    #     :param num_added_variables_orthogonal: number of variables in the orthogonal variable set to add. The more the
    #     better, at the combinatorial cost of memory and computation of course, though at some value the benefit should
    #     diminish to zero (or perhaps negative if the optimization procedure sucks at higher dimensions). If you also
    #     add parallel variables (num_added_variables_parallel > 0) then the quality of the parallel variable set depends
    #     on the quality of the orthogonal variable set, so then a too low number for num_added_variables_orthogonal
    #     might hurt.
    #     :type num_added_variables_orthogonal: int
    #     :type num_added_variables_parallel: int
    #     :param entropy_cost_factor: keep entropy_cost_factor < 1 or (close to) 0.0 (not negative)
    #     :type variables: list of int
    #     """
    #
    #     # remove potential duplicates; and sort (for no reason)
    #     variables_to_orthogonalize = sorted(list(set(variables_to_orthogonalize)))
    #
    #     # due to the scipy's minimize() procedure the list of parameters can be temporarily slightly out of bounds, like
    #     # 1.000000001, but soon after this should be clipped to the allowed range and e.g. at this point the parameters
    #     # are expected to be all valid
    #     assert min(self.matrix2params_incremental()) >= 0.0, \
    #         'parameter(s) out of bound: ' + str(self.matrix2params_incremental())
    #     assert max(self.matrix2params_incremental()) <= 1.0, \
    #         'parameter(s) out of bound: ' + str(self.matrix2params_incremental())
    #
    #     # will add two sets of variables, one orthogonal to the pre-existing 'rest' (i.e. MI zero) and one parallel
    #     # (i.e. MI max). How many variables should each set contain?
    #     # note: setting to len(variables_to_orthogonalize) means that each set has at least surely enough entropy to do the job,
    #     # but given the discrete nature I am not sure if it is also always optimal. Can increase this to be more sure?
    #     # At the cost of more memory usage and computational requirements of course.
    #     if num_added_variables_orthogonal is None:
    #         num_added_variables_orthogonal = len(variables_to_orthogonalize)
    #     if num_added_variables_parallel is None:
    #         num_added_variables_parallel = len(variables_to_orthogonalize)
    #
    #     # variables to orthogonalize against
    #     remaining_variables = sorted([varix for varix in xrange(self.numvariables)
    #                                   if not varix in variables_to_orthogonalize])
    #
    #     if not max(remaining_variables) < min(variables_to_orthogonalize):
    #         pdf_reordered = self.copy()
    #
    #         pdf_reordered.reorder_variables(remaining_variables + variables_to_orthogonalize)
    #
    #         assert len(pdf_reordered) == len(self)
    #         assert len(range(len(remaining_variables), pdf_reordered.numvariables)) == len(variables_to_orthogonalize)
    #
    #         # perform orthogonalization on the ordered pdf, which it is implemented for
    #         pdf_reordered.append_orthogonalized_variables(range(len(remaining_variables),
    #                                                             pdf_reordered.numvariables))
    #
    #         # find the mapping from ordered to disordered, so that I can reverse the reordering above
    #         new_order = remaining_variables + variables_to_orthogonalize
    #         original_order = [-1] * len(new_order)
    #         for varix in xrange(len(new_order)):
    #             original_order[new_order[varix]] = varix
    #
    #         assert not -1 in original_order
    #         assert max(original_order) < len(original_order), 'there were only so many original variables...'
    #
    #         pdf_reordered.reorder_variables(original_order + range(len(original_order), len(pdf_reordered)))
    #
    #         self.duplicate(pdf_reordered)
    #
    #         # due to the scipy's minimize() procedure the list of parameters can be temporarily slightly out of bounds, like
    #         # 1.000000001, but soon after this should be clipped to the allowed range and e.g. at this point the parameters
    #         # are expected to be all valid
    #         assert min(self.matrix2params_incremental()) >= 0.0, \
    #             'parameter(s) out of bound: ' + str(self.matrix2params_incremental())
    #         assert max(self.matrix2params_incremental()) <= 1.0, \
    #             'parameter(s) out of bound: ' + str(self.matrix2params_incremental())
    #     else:
    #         pdf_result = self.copy()  # used to store the result and eventually I will copy to self
    #
    #         ### first add the ORTHOGONAL part
    #
    #         # did not yet do minimize() or something like that so why would a parameter already be out of bound
    #         # even if by a small amount)? should be detected already at some earlier stage.
    #         assert min(pdf_result.matrix2params_incremental()) >= 0.0, 'parameter(s) out of bound, weird: ' \
    #                                                                    + str(pdf_result.matrix2params_incremental())
    #         assert max(pdf_result.matrix2params_incremental()) <= 1.0, 'parameter(s) out of bound, weird' \
    #                                                                    + str(pdf_result.matrix2params_incremental())
    #
    #         if num_added_variables_orthogonal > 0:
    #             pdf_for_ortho_only_optimization = self.copy()
    #             pdf_for_ortho_only_optimization.append_variables(num_added_variables_orthogonal)
    #
    #             num_free_parameters_ortho_only = len(pdf_for_ortho_only_optimization.matrix2params_incremental()) \
    #                                                     - len(self.matrix2params_incremental())
    #
    #             assert num_free_parameters_ortho_only > 0 or num_added_variables_orthogonal == 0
    #
    #             orthogonal_variables = range(len(self), len(self) + num_added_variables_orthogonal)
    #
    #             assert len(orthogonal_variables) == num_added_variables_orthogonal
    #
    #             # did not yet do minimize() or something like that so why would a parameter already be out of bound
    #             # even if by a small amount)? should be detected already at some earlier stage.
    #             assert min(pdf_for_ortho_only_optimization.matrix2params_incremental()) >= 0.0, \
    #                 'parameter(s) out of bound, weird: ' \
    #                 + str(pdf_for_ortho_only_optimization.matrix2params_incremental())
    #             assert max(pdf_for_ortho_only_optimization.matrix2params_incremental()) <= 1.0, \
    #                 'parameter(s) out of bound, weird' \
    #                 + str(pdf_for_ortho_only_optimization.matrix2params_incremental())
    #
    #             def cost_function_orthogonal_part(free_params, static_params, entropy_cost_factor=entropy_cost_factor):
    #                 # note: keep entropy_cost_factor < 1 or (close to) 0.0
    #
    #                 assert len(free_params) == num_free_parameters_ortho_only
    #
    #                 pdf_for_ortho_only_optimization.params2matrix_incremental(list(static_params) + list(free_params))
    #
    #                 # also try to minimize the total entropy of the orthogonal variable set, i.e., try to make the
    #                 # orthogonal part 'efficient' in the sense that it uses only as much entropy as it needs to do its
    #                 # job but no more
    #                 if entropy_cost_factor != 0.0:
    #                     entropy_cost = entropy_cost_factor * pdf_for_ortho_only_optimization.entropy(orthogonal_variables)
    #                 else:
    #                     entropy_cost = 0.0  # do not compute entropy if not used anyway
    #
    #                 # MI with 'remaining_variables' is unwanted, but MI with 'variables_to_orthogonalize' is wanted
    #                 cost_ortho = pdf_for_ortho_only_optimization.mutual_information(remaining_variables,
    #                                                                                 orthogonal_variables) \
    #                              - pdf_for_ortho_only_optimization.mutual_information(variables_to_orthogonalize,
    #                                                                                   orthogonal_variables) \
    #                              + entropy_cost
    #
    #                 return float(cost_ortho)
    #
    #             # if verbose and __debug__:
    #             #     static_param_values = self.matrix2params_incremental()
    #             #     free_param_values = pdf_for_ortho_only_optimization.matrix2params_incremental()[len(self.matrix2params_incremental()):]
    #             #
    #             #     print 'debug: append_orthogonalized_variables: orthogonal cost value before optimization =', \
    #             #         cost_function_orthogonal_part(free_param_values, static_param_values), \
    #             #         '(minimum=' + str(-self.entropy(variables_to_orthogonalize)) + ')'
    #
    #             pdf_result.append_optimized_variables(num_added_variables_orthogonal, cost_function_orthogonal_part,
    #                                                   initial_guess=num_repeats)
    #
    #             # why would a parameter already be out of bound
    #             # even if by a small amount)? should be detected already at some earlier stage.
    #             # append_optimized_variables should already fix this itself.
    #             assert min(pdf_result.matrix2params_incremental()) >= 0.0, \
    #                 'parameter(s) out of bound, weird: ' \
    #                 + str(pdf_result.matrix2params_incremental())
    #             assert max(pdf_result.matrix2params_incremental()) <= 1.0, \
    #                 'parameter(s) out of bound, weird' \
    #                 + str(pdf_result.matrix2params_incremental())
    #
    #             # if verbose and __debug__:
    #             #     static_param_values = self.matrix2params_incremental()
    #             #     free_param_values = pdf_result.matrix2params_incremental()[len(self.matrix2params_incremental()):]
    #             #
    #             #     # test whether the 'static' parameters were indeed kept static during optimization
    #             #     np.testing.assert_array_almost_equal(self.matrix2params_incremental(),
    #             #                                          pdf_result.matrix2params_incremental()[:len(self.matrix2params_incremental())])
    #             #
    #             #     print 'debug: append_orthogonalized_variables: orthogonal cost value after optimization =', \
    #             #         cost_function_orthogonal_part(free_param_values, static_param_values)
    #
    #         ### now add the PARALLEL part
    #
    #         if num_added_variables_parallel > 0:
    #             if num_added_variables_orthogonal == 0:
    #                 raise UserWarning('it is ill-defined to add \'parallel\' variables if I do not have any '
    #                                   '\'orthogonal\' variables to minimize MI against. Just also ask for '
    #                                   'orthogonal variables and then remove them (marginalize all other variables)?')
    #
    #             pdf_for_para_only_optimization = pdf_for_ortho_only_optimization.copy()
    #             # todo: it should be possible to let the parallel variables depend only on the orthogonal_variables
    #             # and the variables_to_orthogonalize, not also the remaining_variables, which would greatly
    #             # reduce the number of free parameters. But then you need to add this artial conditional pdf
    #             # to the complete pdf afterward, repeating it in some way.
    #             pdf_for_para_only_optimization.append_variables(num_added_variables_parallel)
    #
    #             num_free_parameters_para_only = len(pdf_for_para_only_optimization.matrix2params_incremental()) \
    #                                             - len(pdf_for_ortho_only_optimization.matrix2params_incremental())
    #
    #             assert num_free_parameters_para_only > 0 or num_added_variables_parallel == 0
    #
    #             parallel_variables = range(len(pdf_for_ortho_only_optimization),
    #                                        len(pdf_for_para_only_optimization))
    #
    #             assert len(np.intersect1d(parallel_variables, orthogonal_variables)) == 0
    #             assert len(np.intersect1d(parallel_variables, remaining_variables)) == 0
    #
    #             assert len(parallel_variables) == num_added_variables_parallel
    #
    #             # due to the scipy's minimize() procedure the list of parameters can be temporarily slightly
    #             # out of bounds, like
    #             # 1.000000001, but soon after this should be clipped to the allowed range and e.g. at this
    #             # point the parameters
    #             # are expected to be all valid
    #             assert min(pdf_for_para_only_optimization.matrix2params_incremental()) >= 0.0, \
    #                 'parameter(s) out of bound: ' + str(pdf_for_para_only_optimization.matrix2params_incremental())
    #             assert max(pdf_for_para_only_optimization.matrix2params_incremental()) <= 1.0, \
    #                 'parameter(s) out of bound: ' + str(pdf_for_para_only_optimization.matrix2params_incremental())
    #
    #             def cost_function_parallel_part(free_params, static_params, entropy_cost_factor=entropy_cost_factor):
    #                 # note: keep entropy_cost_factor < 1 or (close to) 0.0
    #
    #                 assert len(free_params) == num_free_parameters_para_only
    #                 assert len(free_params) > 0, 'makes no sense to optimize 0 parameters'
    #
    #                 pdf_for_para_only_optimization.params2matrix_incremental(list(static_params) + list(free_params))
    #
    #                 # also try to minimize the total entropy of the parallel variable set, i.e., try to make the
    #                 # parallel part 'efficient' in the sense that it uses only as much entropy as it needs to do its
    #                 # job but no more
    #                 if entropy_cost_factor != 0.0:
    #                     entropy_cost = entropy_cost_factor * pdf_for_para_only_optimization.entropy(parallel_variables)
    #                 else:
    #                     entropy_cost = 0.0  # do not compute entropy if not used anyway
    #
    #                 # MI with 'variables_to_orthogonalize' is wanted, but MI with 'orthogonal_variables' is unwanted
    #                 cost_para = - pdf_for_para_only_optimization.mutual_information(variables_to_orthogonalize,
    #                                                                                 parallel_variables) \
    #                             + pdf_for_para_only_optimization.mutual_information(parallel_variables,
    #                                                                                  orthogonal_variables) \
    #                             + entropy_cost
    #
    #                 return float(cost_para)
    #
    #             if verbose and __debug__:
    #                 static_param_values = pdf_for_ortho_only_optimization.matrix2params_incremental()
    #                 free_param_values = pdf_for_para_only_optimization.matrix2params_incremental()[len(pdf_for_ortho_only_optimization.matrix2params_incremental()):]
    #
    #                 print 'debug: append_orthogonalized_variables: parallel cost value before optimization =', \
    #                     cost_function_parallel_part(free_param_values, static_param_values), \
    #                     '(minimum=' + str(-self.entropy(variables_to_orthogonalize)) + ')'
    #
    #                 # note: this is a probabilistic check: with high probability the 0.0 value is suspected to
    #                 # lead to the asserted condition, but it is also possible that 0.0 just so happens to be due to
    #                 # randomness, however this should be with very small probability
    #                 if cost_function_parallel_part(free_param_values, static_param_values) == 0.0:
    #                     # is this because all entropy of the <variables_to_orthogonalize> is already completely
    #                     # in <orthogonal_variables>, so that for <parallel_variables> there is no more entropy left?
    #                     try:
    #                         np.testing.assert_almost_equal(
    #                             pdf_for_para_only_optimization.mutual_information(remaining_variables, orthogonal_variables)
    #                             ,
    #                             pdf_for_para_only_optimization.entropy(variables_to_orthogonalize)
    #                         )
    #                     except AssertionError as e:
    #                         print 'error: pdf_for_para_only_optimization.' \
    #                               'mutual_information(remaining_variables, orthogonal_variables) =', \
    #                             pdf_for_para_only_optimization.mutual_information(remaining_variables,
    #                                                                               orthogonal_variables)
    #                         print 'error: pdf_for_para_only_optimization.entropy(variables_to_orthogonalize) =', \
    #                             pdf_for_para_only_optimization.entropy(variables_to_orthogonalize)
    #                         print 'error: pdf_for_para_only_optimization.' \
    #                               'mutual_information(parallel_variables, orthogonal_variables) =', \
    #                             pdf_for_para_only_optimization.mutual_information(parallel_variables,
    #                                                                                orthogonal_variables)
    #                         print 'error: pdf_for_para_only_optimization.' \
    #                               'mutual_information(variables_to_orthogonalize, parallel_variables) =', \
    #                             pdf_for_para_only_optimization.mutual_information(variables_to_orthogonalize,
    #                                                                               parallel_variables)
    #                         print 'error: pdf_for_para_only_optimization.entropy(remaining_variables)', \
    #                             pdf_for_para_only_optimization.entropy(remaining_variables)
    #                         print 'error: pdf_for_para_only_optimization.entropy(orthogonal_variables) =', \
    #                             pdf_for_para_only_optimization.entropy(orthogonal_variables)
    #                         print 'error: pdf_for_para_only_optimization.entropy(parallel_variables) =', \
    #                             pdf_for_para_only_optimization.entropy(parallel_variables)
    #
    #                         raise AssertionError(e)
    #
    #             pdf_result.append_optimized_variables(num_added_variables_parallel, cost_function_parallel_part)
    #
    #             # due to the scipy's minimize() procedure the list of parameters can be temporarily slightly
    #             # out of bounds, like
    #             # 1.000000001, but soon after this should be clipped to the allowed range and e.g. at this
    #             # point the parameters
    #             # are expected to be all valid
    #             assert min(pdf_result.matrix2params_incremental()) >= 0.0, \
    #                 'parameter(s) out of bound: ' + str(pdf_result.matrix2params_incremental())
    #             assert max(pdf_result.matrix2params_incremental()) <= 1.0, \
    #                 'parameter(s) out of bound: ' + str(pdf_result.matrix2params_incremental())
    #
    #             if verbose and __debug__:
    #                 static_param_values = pdf_for_ortho_only_optimization.matrix2params_incremental()
    #                 free_param_values = pdf_result.matrix2params_incremental()[len(pdf_for_ortho_only_optimization.matrix2params_incremental()):]
    #
    #                 print 'debug: append_orthogonalized_variables: parallel cost value after optimization =', \
    #                     cost_function_parallel_part(free_param_values, static_param_values), \
    #                     '(minimum=' + str(-self.entropy(variables_to_orthogonalize)) + ')'
    #
    #                 # note: this is a probabilistic check: with high probability the 0.0 value is suspected to
    #                 # lead to the asserted condition, but it is also possible that 0.0 just so happens to be due to
    #                 # randomness, however this should be with very small probability
    #                 if cost_function_parallel_part(free_param_values, static_param_values) == 0.0:
    #                     # is this because all entropy of the <variables_to_orthogonalize> is already completely
    #                     # in <orthogonal_variables>, so that for <parallel_variables> there is no more entropy left?
    #                     np.testing.assert_almost_equal(
    #                         pdf_for_para_only_optimization.mutual_information(remaining_variables, orthogonal_variables)
    #                         ,
    #                         pdf_for_para_only_optimization.entropy(variables_to_orthogonalize)
    #                     )
    #
    #         self.duplicate(pdf_result)
    #
    #
    #
    #     # todo: add some debug tolerance measures to detect if the error is too large in orthogonalization, like >10%
    #     # of entropy of orthogonal_variables is MI with parallel_variables?


    def scalars_up_to_level(self, list_of_lists, max_level=None):
        """
        Helper function. E.g. scalars_up_to_level([1,[2,3],[[4]]]) == [1], and
        scalars_up_to_level([1,[2,3],[[4]]], max_level=2) == [1,2,3]. Will be sorted on level, with highest level
        scalars first.

        Note: this function is not very efficiently implemented I think, but my concern now is that it works at all.

        :type list_of_lists: list
        :type max_level: int
        :rtype: list
        """
        # scalars = [v for v in list_of_lists if np.isscalar(v)]
        #
        # if max_level > 1 or (max_level is None and len(list_of_lists) > 0):
        #     for sublist in [v for v in list_of_lists if not np.isscalar(v)]:
        #         scalars.extend(self.scalars_up_to_level(sublist,
        #                                                 max_level=max_level-1 if not max_level is None else None))

        scalars = []

        if __debug__:
            debug_max_depth_set = (max_level is None)

        if max_level is None:
            max_level = maximum_depth(list_of_lists)

        for at_level in xrange(1, max_level + 1):
            scalars_at_level = self.scalars_at_level(list_of_lists, at_level=at_level)

            scalars.extend(scalars_at_level)

        if __debug__:
            if debug_max_depth_set:
                assert len(scalars) == len(flatten(list_of_lists)), 'all scalars should be present, and not duplicate' \
                                                                    '. len(scalars) = ' + str(len(scalars)) \
                                                                    + ', len(flatten(list_of_lists)) = ' \
                                                                    + str(len(flatten(list_of_lists)))

        return scalars

    def scalars_at_level(self, list_of_lists, at_level=1):
        """
        Helper function. E.g. scalars_up_to_level([1,[2,3],[[4]]]) == [1], and
        scalars_up_to_level([1,[2,3],[[4]]], max_level=2) == [1,2,3]. Will be sorted on level, with highest level
        scalars first.
        :type list_of_lists: list
        :type max_level: int
        :rtype: list
        """

        if at_level == 1:
            scalars = [v for v in list_of_lists if np.isscalar(v)]

            return scalars
        elif at_level == 0:
            warnings.warn('level 0 does not exist, I start counting from at_level=1, will return [].')

            return []
        else:
            scalars = []

            for sublist in [v for v in list_of_lists if not np.isscalar(v)]:
                scalars.extend(self.scalars_at_level(sublist, at_level=at_level-1))

            assert np.ndim(scalars) == 1

            return scalars

    def imbalanced_tree_from_scalars(self, list_of_scalars, numvalues):
        """
        Helper function.
        Consider e.g. tree =
                        [0.36227870614214747,
                         0.48474422004766832,
                         [0.34019329926554265,
                          0.40787146599658614,
                          [0.11638879037422999, 0.64823088842780996],
                          [0.33155311703042312, 0.11398958845340294],
                          [0.13824154613818085, 0.42816388506114755]],
                         [0.15806602176772611,
                          0.32551465875945773,
                          [0.25748947995256499, 0.35415524846620511],
                          [0.64896559115417218, 0.65575802084978507],
                          [0.36051945555508391, 0.40134903827671109]],
                         [0.40568439663760192,
                          0.67602830725264651,
                          [0.35103999983495449, 0.59577145940649334],
                          [0.38917741342947187, 0.44327101890582132],
                          [0.075034425516081762, 0.59660319391007388]]]

        If you first call scalars_up_to_level on this you get a list [0.36227870614214747, 0.48474422004766832,
        0.34019329926554265, 0.40787146599658614, 0.15806602176772611, ...]. If you pass this flattened list through
        this function then you should get the above imbalanced tree structure back again.

        At each level in the resulting tree there will be <numvalues-1> scalars and <numvalues> subtrees (lists).
        :type list_of_scalars: list
        :type numvalues: int
        :rtype: list
        """

        num_levels = int(np.round(np.log2(len(list_of_scalars) + 1) / np.log2(numvalues)))

        all_scalars_at_level = dict()

        list_of_scalars_remaining = list(list_of_scalars)

        for level in xrange(num_levels):
            num_scalars_at_level = np.power(numvalues, level) * (numvalues - 1)

            scalars_at_level = list_of_scalars_remaining[:num_scalars_at_level]

            all_scalars_at_level[level] = scalars_at_level

            list_of_scalars_remaining = list_of_scalars_remaining[num_scalars_at_level:]

        # todo: further order the scalars in the expected way
        def tree_from_levels(all_scalars_at_level):
            if len(all_scalars_at_level) == 0:
                return []
            else:
                assert len(all_scalars_at_level[0]) == numvalues - 1

                if len(all_scalars_at_level) > 1:
                    assert len(all_scalars_at_level[1]) == numvalues * (numvalues - 1)
                if len(all_scalars_at_level) > 2:
                    assert len(all_scalars_at_level[2]) == (numvalues*numvalues) * (numvalues - 1), \
                        'len(all_scalars_at_level[2]) = ' + str(len(all_scalars_at_level[2])) + ', ' \
                        '(numvalues*numvalues) * (numvalues - 1) = ' + str((numvalues*numvalues) * (numvalues - 1))
                if len(all_scalars_at_level) > 3:
                    assert len(all_scalars_at_level[3]) == (numvalues*numvalues*numvalues) * (numvalues - 1)
                # etc.

                tree = list(all_scalars_at_level[0][:(numvalues - 1)])

                if len(all_scalars_at_level) > 1:
                    # add <numvalues> subtrees to this level
                    for subtree_id in xrange(numvalues):
                        all_scalars_for_subtree = dict()

                        for level in xrange(len(all_scalars_at_level) - 1):
                            num_scalars_at_level = len(all_scalars_at_level[level + 1])

                            assert np.mod(num_scalars_at_level, numvalues) == 0, 'should be divisible nu <numvalues>'

                            num_scalars_for_subtree = int(num_scalars_at_level / numvalues)

                            all_scalars_for_subtree[level] = \
                                all_scalars_at_level[level + 1][subtree_id * num_scalars_for_subtree
                                                                :(subtree_id + 1) * num_scalars_for_subtree]

                        subtree_i = tree_from_levels(all_scalars_for_subtree)

                        if len(all_scalars_for_subtree) > 1:
                            # numvalues - 1 scalars and numvalues subtrees
                            assert len(subtree_i) == (numvalues - 1) + numvalues, 'len(subtree_i) = ' \
                                                                                  + str(len(subtree_i)) \
                                                                                  + ', expected = ' \
                                                                                  + str((numvalues - 1) + numvalues)
                        elif len(all_scalars_for_subtree) == 1:
                            assert len(subtree_i) == numvalues - 1

                        tree.append(subtree_i)

                return tree

        tree = tree_from_levels(all_scalars_at_level)

        assert maximum_depth(tree) == len(all_scalars_at_level)  # should be numvariables if the scalars are parameters
        assert len(flatten(tree)) == len(list_of_scalars), 'all scalars should end up in the tree, and not duplicate'

        return tree



### UNIT TESTING:



def run_append_and_marginalize_test():
    pdf = JointProbabilityMatrix(3, 2)

    pdf_copy = pdf.copy()

    pdf.append_variables(4)

    pdf_old = pdf.marginalize_distribution(range(3))

    assert pdf_copy == pdf_old, 'adding and then removing variables should result in the same joint pdf'

    # old_params = pdf_copy.matrix2params_incremental()
    #
    # np.testing.assert_array_almost_equal(pdf.matrix2params_incremental()[:len(old_params)], old_params)


def run_reorder_test():
    pdf = JointProbabilityMatrix(4, 3)

    pdf_original = pdf.copy()

    pdf.reorder_variables([3,2,1,0])
    np.testing.assert_almost_equal(pdf.entropy(), pdf_original.entropy())
    assert pdf != pdf_original

    pdf.reorder_variables([3,2,1,0,0,1,2,3])
    assert len(pdf) == 2 * len(pdf_original)
    np.testing.assert_almost_equal(pdf.entropy(), pdf_original.entropy())
    # the first 4 and second 4 variables should be identical so MI should be maximum
    np.testing.assert_almost_equal(pdf.mutual_information([0,1,2,3], [4,5,6,7]), pdf_original.entropy())

    assert pdf.marginalize_distribution([0,1,2,3]) == pdf_original


def run_params2matrix_test():
    pdf = JointProbabilityMatrix(3, 2)

    pdf_copy = pdf.copy()

    pdf_copy.params2matrix(pdf.matrix2params())

    assert pdf_copy == pdf, 'computing parameter values from joint pdf and using those to construct a 2nd joint pdf ' \
                            'should result in two equal pdfs.'


def run_vector2matrix_test():
    pdf = JointProbabilityMatrix(3, 2)

    pdf_copy = pdf.copy()

    pdf_copy.vector2matrix(pdf.matrix2vector())

    assert pdf_copy == pdf, 'computing vector from joint pdf and using that to construct a 2nd joint pdf ' \
                            'should result in two equal pdfs.'


def run_conditional_pdf_test():
    pdf = JointProbabilityMatrix(3, 2)

    pdf_marginal_1 = pdf.marginalize_distribution([1])

    assert pdf_marginal_1.numvariables == 1

    pdf_cond_23_given_0 = pdf.conditional_probability_distribution([1], [0])
    pdf_cond_23_given_1 = pdf.conditional_probability_distribution([1], [1])

    assert pdf_cond_23_given_0.numvariables == 2
    assert pdf_cond_23_given_1.numvariables == 2

    prob_000_joint = pdf([0,0,0])
    prob_000_cond = pdf_marginal_1([0]) * pdf_cond_23_given_0([0,0])

    np.testing.assert_almost_equal(prob_000_cond, prob_000_joint)

    pdf_conds_23 = pdf.conditional_probability_distributions([1])

    assert pdf_conds_23[(0,)] == pdf_cond_23_given_0
    assert pdf_conds_23[(1,)] == pdf_cond_23_given_1


def run_append_using_transitions_table_and_marginalize_test():
    pdf = JointProbabilityMatrix(3, 2)

    pdf_copy = pdf.copy()

    lists_of_possible_given_values = [range(pdf.numvalues) for _ in xrange(pdf.numvariables)]

    state_transitions = [list(existing_vars_values) + list([int(np.mod(np.sum(existing_vars_values), pdf.numvalues))])
                         for existing_vars_values in itertools.product(*lists_of_possible_given_values)]

    pdf.append_variables_using_state_transitions_table(state_transitions)

    assert not hasattr(state_transitions, '__call__'), 'append_variables_using_state_transitions_table should not ' \
                                                       'replace the caller\'s variables'

    assert pdf.numvariables == 4, 'one variable should be added'

    pdf_old = pdf.marginalize_distribution(range(pdf_copy.numvariables))

    assert pdf_copy == pdf_old, 'adding and then removing variables should result in the same joint pdf'
    
    pdf_copy.append_variables_using_state_transitions_table(
        state_transitions=lambda vals, mv: [int(np.mod(np.sum(vals), mv))])

    assert pdf_copy == pdf, 'should be two equivalent ways of appending a deterministic variable'


def run_synergistic_variables_test_with_subjects_and_agnostics(num_subject_vars=2):
    num_agnostic_vars = 1

    pdf = JointProbabilityMatrix(num_subject_vars + num_agnostic_vars, 2)

    pdf_orig = pdf.copy()

    assert num_subject_vars > 0

    subject_variables = np.random.choice(range(len(pdf)), num_subject_vars, replace=False)
    agnostic_variables = np.setdiff1d(range(len(pdf)), subject_variables)

    pdf.append_synergistic_variables(1, subject_variables=subject_variables)

    assert pdf.marginalize_distribution(range(len(pdf_orig))) == pdf_orig, \
        'appending synergistic variables changed the pdf'
    
    synergistic_variables = range(len(pdf_orig), len(pdf))
    
    tol_rel_err = 0.2

    # ideally, the former is max and the latter is zero
    assert pdf.mutual_information(subject_variables, synergistic_variables) \
           > sum([pdf.mutual_information([sv], synergistic_variables) for sv in subject_variables]) \
           or pdf.entropy(synergistic_variables) < 0.01
    
    assert pdf.mutual_information(synergistic_variables, agnostic_variables) / pdf.entropy(synergistic_variables) \
           < tol_rel_err


def run_synergistic_variables_test_with_subjects(num_subject_vars=2):
    num_other_vars = 1

    pdf = JointProbabilityMatrix(num_subject_vars + num_other_vars, 2)

    pdf_orig = pdf.copy()

    assert num_subject_vars > 0

    subject_variables = np.random.choice(range(len(pdf)), num_subject_vars, replace=False)

    pdf.append_synergistic_variables(1, subject_variables=subject_variables)

    assert pdf.marginalize_distribution(range(len(pdf_orig))) == pdf_orig, \
        'appending synergistic variables changed the pdf'

    # ideally, the former is max and the latter is zero
    assert pdf.mutual_information(subject_variables, range(len(pdf_orig), len(pdf))) \
           > sum([pdf.mutual_information([sv], range(len(pdf_orig), len(pdf))) for sv in subject_variables])


# todo: make another such function but now use the subject_variables option of append_synerg*
def run_synergistic_variables_test(numvars=2):
    pdf = JointProbabilityMatrix(numvars, 2)

    pdf_syn = pdf.copy()

    assert pdf_syn == pdf

    # initial_guess_summed_modulo = np.random.choice([True, False])
    initial_guess_summed_modulo = False

    pdf_syn.append_synergistic_variables(1, initial_guess_summed_modulo=initial_guess_summed_modulo, verbose=False)

    assert pdf_syn.numvariables == pdf.numvariables + 1

    pdf_old = pdf_syn.marginalize_distribution(range(pdf.numvariables))

    # trying to figure out why I hit the assertion "pdf == pdf_old"
    np.testing.assert_almost_equal(np.sum(pdf_old.joint_probabilities), 1.0), 'all probabilities should sum to 1.0'
    np.testing.assert_almost_equal(np.sum(pdf.joint_probabilities), 1.0), 'all probabilities should sum to 1.0'

    np.testing.assert_array_almost_equal(pdf.joint_probabilities, pdf_old.joint_probabilities)
    assert pdf == pdf_old, 'adding and then removing variables should result in the same joint pdf'

    parameters_before = pdf.matrix2params_incremental()

    pdf_add_random = pdf.copy()
    pdf_add_random.append_variables(1)

    np.testing.assert_array_almost_equal(pdf_add_random.matrix2params_incremental()[:len(parameters_before)],
                                                parameters_before)
    np.testing.assert_array_almost_equal(pdf_syn.matrix2params_incremental()[:len(parameters_before)],
                                                parameters_before)

    # note: this assert is in principle probabilistic, because who knows what optimization procedure is used and
    # how much it potentially sucks. So see if you hit this more than once, if you hit it at all.
    assert pdf_add_random.synergistic_information_naive(range(pdf.numvariables, pdf_add_random.numvariables),
                                                                range(pdf.numvariables)) <= \
           pdf_syn.synergistic_information_naive(range(pdf.numvariables, pdf_add_random.numvariables),
                                                         range(pdf.numvariables)), 'surely the optimization procedure' \
                                                                                   ' in append_synergistic_variables ' \
                                                                                   'yields a better syn. info. than ' \
                                                                                   'an appended variable with simply ' \
                                                                                   'random interaction parameters?!'

    np.testing.assert_array_almost_equal(pdf_add_random.matrix2params_incremental()[:len(parameters_before)],
                                                parameters_before)
    np.testing.assert_array_almost_equal(pdf_syn.matrix2params_incremental()[:len(parameters_before)],
                                                parameters_before)

    syninfo = pdf_syn.synergistic_information_naive(range(pdf.numvariables, pdf_add_random.numvariables),
                                                    range(pdf.numvariables))

    condents = [pdf.conditional_entropy([varix]) for varix in xrange(len(pdf))]

    assert syninfo <= min(condents), 'this is a derived maximum in Quax 2015, synergy paper, right?'


def run_orthogonalization_test_null_hypothesis(num_subject_vars=2, num_ortho_vars=1, num_para_vars=1, numvals=2,
                                               verbose=True, num_repeats=5, tol_rel_error=0.05):
    """
    This is similar in spirit to run_orthogonalization_test, except now the variables to orthogonalize (X) is already
    known to be completely decomposable into two parts (X1, X2) which are orthogonal and parallel, resp. (See also
    the description of append_orthogonalized_variables().) So in this case I am sure that
    append_orthogonalized_variables should find very good solutions.
    :param tol_rel_error: A float in [0, 1) and preferably close to 0.0.
    :param num_subject_vars:
    :param num_ortho_vars:
    :param num_para_vars:
    :param numvals:
    :param verbose:
    :param num_repeats:
    """
    pdf = JointProbabilityMatrix(num_subject_vars, numvals)

    # note: I add the 'null hypothesis' parallel and orthogonal parts in reversed order compared to what
    # append_orthogonalized_variables returns, otherwise I have to insert a reordering or implement
    # append_redundant_variables for subsets of variables (not too hard but anyway)

    pdf.append_redundant_variables(num_para_vars)

    # the just added variables should be completely redundant
    np.testing.assert_almost_equal(pdf.mutual_information(range(num_subject_vars),
                                                          range(num_subject_vars, num_subject_vars + num_para_vars)),
                                   pdf.entropy(range(num_subject_vars, num_subject_vars + num_para_vars)))

    pdf.append_independent_variables(JointProbabilityMatrix(num_ortho_vars, pdf.numvalues))

    # the just added variables should be completely independent of all previous
    np.testing.assert_almost_equal(pdf.mutual_information(range(num_subject_vars, num_subject_vars + num_para_vars),
                                                          range(num_subject_vars + num_para_vars,
                                                                num_subject_vars + num_para_vars + num_ortho_vars)),
                                   0.0)

    vars_to_orthogonalize = range(num_subject_vars, num_subject_vars + num_para_vars + num_ortho_vars)

    pdf.append_orthogonalized_variables(vars_to_orthogonalize, num_ortho_vars, num_para_vars,
                                        num_repeats=num_repeats)

    assert len(pdf) == num_subject_vars + num_ortho_vars + num_para_vars + num_ortho_vars + num_para_vars

    result_ortho_vars = range(num_subject_vars + num_ortho_vars + num_para_vars,
                              num_subject_vars + num_ortho_vars + num_para_vars + num_ortho_vars)
    result_para_vars = range(num_subject_vars + num_ortho_vars + num_para_vars + num_ortho_vars,
                             num_subject_vars + num_ortho_vars + num_para_vars + num_ortho_vars + num_para_vars)
    subject_vars = range(num_subject_vars)

    if pdf.entropy(result_ortho_vars) != 0.0:
        try:
            assert pdf.mutual_information(subject_vars, result_ortho_vars) / pdf.entropy(result_ortho_vars) \
                   <= tol_rel_error
        except AssertionError as e:

            print 'debug: pdf.mutual_information(subject_vars, result_ortho_vars) =', \
                pdf.mutual_information(subject_vars, result_ortho_vars)
            print 'debug: pdf.entropy(result_ortho_vars) =', pdf.entropy(result_ortho_vars)
            print 'debug: pdf.mutual_information(subject_vars, result_ortho_vars) / pdf.entropy(result_ortho_vars) =', \
                pdf.mutual_information(subject_vars, result_ortho_vars) / pdf.entropy(result_ortho_vars)
            print 'debug: (ideal of previous quantity: 0.0)'
            print 'debug: tol_rel_error =', tol_rel_error

            raise AssertionError(e)

    if pdf.entropy(result_para_vars) != 0.0:
        assert pdf.mutual_information(subject_vars, result_para_vars) \
               / pdf.entropy(result_para_vars) >= 1.0 - tol_rel_error

    if pdf.entropy(result_ortho_vars) != 0.0:
        assert pdf.mutual_information(result_para_vars, result_ortho_vars) \
               / pdf.entropy(result_ortho_vars) <= tol_rel_error

    if pdf.entropy(result_para_vars) != 0.0:
        assert pdf.mutual_information(result_para_vars, result_ortho_vars) \
               / pdf.entropy(result_para_vars) <= tol_rel_error

    if pdf.entropy(vars_to_orthogonalize) != 0.0:
        try:
            assert pdf.mutual_information(vars_to_orthogonalize, list(result_para_vars) + list(result_ortho_vars)) \
                   / pdf.entropy(vars_to_orthogonalize) >= 1.0 - tol_rel_error, \
                'not all entropy of X is accounted for in {X1, X2}, which is of course the purpose of decomposition.'
        except AssertionError as e:
            print 'debug: pdf.mutual_information(vars_to_orthogonalize, ' \
                  'list(result_para_vars) + list(result_ortho_vars)) =', \
                pdf.mutual_information(vars_to_orthogonalize, list(result_para_vars) + list(result_ortho_vars))
            print 'debug: pdf.entropy(vars_to_orthogonalize) =', pdf.entropy(vars_to_orthogonalize)
            print 'debug: pdf.mutual_information(vars_to_orthogonalize, ' \
                  'list(result_para_vars) + list(result_ortho_vars))' \
                   '/ pdf.entropy(vars_to_orthogonalize) =', \
                pdf.mutual_information(vars_to_orthogonalize, list(result_para_vars) + list(result_ortho_vars)) \
                   / pdf.entropy(vars_to_orthogonalize)
            print 'debug: tol_rel_error =', tol_rel_error


def run_orthogonalization_test(num_subject_vars=2, num_orthogonalized_vars=1, numvals=3, verbose=True, num_repeats=1):
    """

    :param num_subject_vars:
    :param num_orthogonalized_vars:
    :param numvals:
    :param verbose:
    :param num_repeats:
    :raise AssertionError:
    """

    # note: in the code I assume that the number of added orthogonal variables equals num_orthogonalized_vars
    # note: in the code I assume that the number of added parallel variables equals num_orthogonalized_vars

    # todo: I am checking if it should be generally expected at all that any given X for any (Y, X) can be decomposed
    # into X1, X2 (orthogonal and parallel, resp.). Maybe in some cases it is not possible, due to the discrete nature
    # of the variables.

    pdf = JointProbabilityMatrix(num_subject_vars + num_orthogonalized_vars, numvals)

    pdf_original = pdf.copy()

    subject_vars = range(num_subject_vars)
    # these are the variables that will be 'orthogonalized', i.e., naming this X then the below appending procedure
    # will find two variable sets {X1,X2}=X where X1 is 'orthogonal' to original_vars (MI=0) and X2 is 'parallel'
    vars_to_orthogonalize = range(num_subject_vars, num_subject_vars + num_orthogonalized_vars)

    pdf.append_orthogonalized_variables(vars_to_orthogonalize, num_orthogonalized_vars, num_orthogonalized_vars,
                                        num_repeats=num_repeats)

    if verbose:
        print 'debug: computed first orthogonalization'

    assert len(pdf) == len(pdf_original) + num_orthogonalized_vars * 2

    # note: implicitly, the number of variables added for the orthogonal part is <num_orthogonalized_vars> here
    ortho_vars = range(num_subject_vars + num_orthogonalized_vars,
                       num_subject_vars + num_orthogonalized_vars + num_orthogonalized_vars)
    # note: implicitly, the number of variables added for the parallel part is <num_orthogonalized_vars> here
    para_vars = range(num_subject_vars + num_orthogonalized_vars + num_orthogonalized_vars,
                      num_subject_vars + num_orthogonalized_vars + num_orthogonalized_vars + num_orthogonalized_vars)

    assert para_vars[-1] == len(pdf) - 1, 'not all variables accounted for?'

    '''
    these should be high:
    pdf.mutual_information(original_vars, para_vars)
    pdf.mutual_information(vars_to_orthogonalize, ortho_vars)

    these should be low:
    pdf.mutual_information(original_vars, ortho_vars)
    pdf.mutual_information(para_vars, ortho_vars)
    '''

    tol_rel_err = 0.2  # used for all kinds of tests

    # some test of ordering.
    # Here: more % of the parallel variables' entropy should be correlated with subject_variables
    # than the same % for the orthogonal variables. This should be a pretty weak requirement.
    assert pdf.mutual_information(subject_vars, para_vars) / pdf.entropy(para_vars) \
           > pdf.mutual_information(subject_vars, ortho_vars) / pdf.entropy(ortho_vars), \
        '1: pdf.mutual_information(vars_to_orthogonalize, para_vars) = ' + str(pdf.mutual_information(vars_to_orthogonalize, para_vars)) \
        + ', pdf.mutual_information(para_vars, ortho_vars) = ' + str(pdf.mutual_information(para_vars, ortho_vars)) \
        + ', vars_to_orthogonalize=' + str(vars_to_orthogonalize) + ', para_vars=' + str(para_vars) \
        + ', ortho_vars=' + str(ortho_vars) + ', subject_vars=' + str(subject_vars) \
        + ', pdf.mutual_information(subject_vars, para_vars) = ' + str(pdf.mutual_information(subject_vars, para_vars)) \
        + ', pdf.mutual_information(subject_vars, ortho_vars) = ' + str(pdf.mutual_information(subject_vars, ortho_vars)) \
        + ', pdf.entropy(subject_vars) = ' + str(pdf.entropy(subject_vars)) \
        + ', pdf.entropy(vars_to_orthogonalize) = ' + str(pdf.entropy(vars_to_orthogonalize)) \
        + ', pdf.entropy(ortho_vars) = ' + str(pdf.entropy(ortho_vars)) \
        + ', pdf.entropy(para_vars) = ' + str(pdf.entropy(para_vars)) \
        + ', pdf.mutual_information(vars_to_orthogonalize, ortho_vars) = ' \
        + str(pdf.mutual_information(vars_to_orthogonalize, ortho_vars))

    assert pdf.mutual_information(vars_to_orthogonalize, ortho_vars) > pdf.mutual_information(para_vars, ortho_vars), \
        '2: pdf.mutual_information(vars_to_orthogonalize, para_vars) = ' + str(pdf.mutual_information(vars_to_orthogonalize, para_vars)) \
        + ', pdf.mutual_information(para_vars, ortho_vars) = ' + str(pdf.mutual_information(para_vars, ortho_vars)) \
        + ', vars_to_orthogonalize=' + str(vars_to_orthogonalize) + ', para_vars=' + str(para_vars) \
        + ', ortho_vars=' + str(ortho_vars) + ', subject_vars=' + str(subject_vars) \
        + ', pdf.mutual_information(subject_vars, para_vars) = ' + str(pdf.mutual_information(subject_vars, para_vars)) \
        + ', pdf.mutual_information(subject_vars, ortho_vars) = ' + str(pdf.mutual_information(subject_vars, ortho_vars)) \
        + ', pdf.entropy(subject_vars) = ' + str(pdf.entropy(subject_vars)) \
        + ', pdf.entropy(vars_to_orthogonalize) = ' + str(pdf.entropy(vars_to_orthogonalize)) \
        + ', pdf.entropy(ortho_vars) = ' + str(pdf.entropy(ortho_vars)) \
        + ', pdf.entropy(para_vars) = ' + str(pdf.entropy(para_vars)) \
        + ', pdf.mutual_information(vars_to_orthogonalize, ortho_vars) = ' \
        + str(pdf.mutual_information(vars_to_orthogonalize, ortho_vars))

    assert pdf.mutual_information(vars_to_orthogonalize, para_vars) > pdf.mutual_information(para_vars, ortho_vars), \
        '3: pdf.mutual_information(vars_to_orthogonalize, para_vars) = ' + str(pdf.mutual_information(vars_to_orthogonalize, para_vars)) \
        + ', pdf.mutual_information(para_vars, ortho_vars) = ' + str(pdf.mutual_information(para_vars, ortho_vars)) \
        + ', vars_to_orthogonalize=' + str(vars_to_orthogonalize) + ', para_vars=' + str(para_vars) \
        + ', ortho_vars=' + str(ortho_vars) + ', subject_vars=' + str(subject_vars) \
        + ', pdf.mutual_information(subject_vars, para_vars) = ' + str(pdf.mutual_information(subject_vars, para_vars)) \
        + ', pdf.mutual_information(subject_vars, ortho_vars) = ' + str(pdf.mutual_information(subject_vars, ortho_vars)) \
        + ', pdf.entropy(subject_vars) = ' + str(pdf.entropy(subject_vars)) \
        + ', pdf.entropy(vars_to_orthogonalize) = ' + str(pdf.entropy(vars_to_orthogonalize)) \
        + ', pdf.entropy(ortho_vars) = ' + str(pdf.entropy(ortho_vars)) \
        + ', pdf.entropy(para_vars) = ' + str(pdf.entropy(para_vars)) \
        + ', pdf.mutual_information(vars_to_orthogonalize, ortho_vars) = ' \
        + str(pdf.mutual_information(vars_to_orthogonalize, ortho_vars))

    # at least more than <tol*>% of the entropy of para_vars is not 'wrong' (correlated with ortho_vars)
    assert pdf.mutual_information(para_vars, ortho_vars) < tol_rel_err * pdf.entropy(para_vars)
    # at least more than <tol*>% of the entropy of ortho_vars is not 'wrong' (correlated with para_vars)
    assert pdf.mutual_information(para_vars, ortho_vars) < tol_rel_err * pdf.entropy(ortho_vars)

    # todo: print some numbers to get a feeling of what (verbose) kind of accuracy I get, or let
    # append_orthogonalized_variables return those (e.g. in a dict?)?

    if verbose:
        print 'debug: computed first bunch of MI checks'

    ### test if the total entropy of ortho_vars + para_vars is close enough to the entropy of vars_to_orthogonalize

    # note: cannot directly use pdf.entropy(ortho_vars + para_vars) because I did not yet include a cost term
    # for the total entropy ('efficiency') of the ortho and para vars (so far entropy_cost_factor=0.0)
    entropy_ortho_and_para = pdf.mutual_information(vars_to_orthogonalize, ortho_vars + para_vars)
    entropy_vars_to_ortho = pdf.entropy(vars_to_orthogonalize)

    if verbose:
        print 'debug: computed few more entropy things'

    if entropy_vars_to_ortho != 0.0:
        assert abs(entropy_vars_to_ortho - entropy_ortho_and_para) / entropy_vars_to_ortho <= tol_rel_err, \
            'the total entropy of the ortho and para vars is too high/low: ' + str(entropy_ortho_and_para) + ' versus' \
            + ' entropy_vars_to_ortho=' + str(entropy_vars_to_ortho) + ' (rel. err. tol. = ' + str(tol_rel_err) + ')'

    ideal_para_entropy = pdf_original.mutual_information(vars_to_orthogonalize, subject_vars)

    assert pdf.entropy(para_vars) >= (1.0 - tol_rel_err) * ideal_para_entropy
    try:
        assert pdf.mutual_information(vars_to_orthogonalize, para_vars) >= (1.0 - tol_rel_err) * ideal_para_entropy
    except AssertionError as e:
        print 'debug: pdf.mutual_information(vars_to_orthogonalize, para_vars) =', \
            pdf.mutual_information(vars_to_orthogonalize, para_vars)
        print 'debug: pdf.mutual_information(subject_vars, para_vars) =', \
            pdf.mutual_information(subject_vars, para_vars)
        print 'debug: ideal_para_entropy =', ideal_para_entropy
        print 'debug: tol_rel_err =', tol_rel_err

        raise AssertionError(e)

    if verbose:
        print 'debug: ...and more...'

    ideal_ortho_entropy = pdf_original.conditional_entropy(vars_to_orthogonalize)

    # H(X1) should be close to H(X|Y), i.e., the entropy in X which does not share information with Y
    if not pdf.entropy(ortho_vars) == 0.0:
        assert abs(pdf.entropy(ortho_vars) - ideal_ortho_entropy) / pdf.entropy(ortho_vars) <= tol_rel_err

    # I(X1:X) should be (almost) equal to H(X1), i.e., all entropy of the orthogonal X1 should be from X, nowhere
    # else
    if not pdf.mutual_information(vars_to_orthogonalize, ortho_vars) == 0.0:
        assert abs(pdf.mutual_information(vars_to_orthogonalize, ortho_vars) - ideal_ortho_entropy) \
               / pdf.mutual_information(vars_to_orthogonalize, ortho_vars) <= tol_rel_err, \
        'pdf.mutual_information(vars_to_orthogonalize, para_vars) = ' + str(pdf.mutual_information(vars_to_orthogonalize, para_vars)) \
        + ', pdf.mutual_information(para_vars, ortho_vars) = ' + str(pdf.mutual_information(para_vars, ortho_vars)) \
        + ', vars_to_orthogonalize=' + str(vars_to_orthogonalize) + ', para_vars=' + str(para_vars) \
        + ', ortho_vars=' + str(ortho_vars) + ', subject_vars=' + str(subject_vars) \
        + ', pdf.mutual_information(subject_vars, para_vars) = ' + str(pdf.mutual_information(subject_vars, para_vars)) \
        + ', pdf.mutual_information(subject_vars, ortho_vars) = ' + str(pdf.mutual_information(subject_vars, ortho_vars)) \
        + ', pdf.entropy(subject_vars) = ' + str(pdf.entropy(subject_vars)) \
        + ', pdf.entropy(vars_to_orthogonalize) = ' + str(pdf.entropy(vars_to_orthogonalize)) \
        + ', pdf.entropy(ortho_vars) = ' + str(pdf.entropy(ortho_vars)) \
        + ', pdf.entropy(para_vars) = ' + str(pdf.entropy(para_vars)) \
        + ', pdf.mutual_information(vars_to_orthogonalize, ortho_vars) = ' \
        + str(pdf.mutual_information(vars_to_orthogonalize, ortho_vars))

    if verbose:
        print 'debug: done!'



def run_append_conditional_pdf_test():
    pdf_joint = JointProbabilityMatrix(4, 3)

    pdf_12 = pdf_joint.marginalize_distribution([0, 1])
    pdf_34_cond_12 = pdf_joint.conditional_probability_distributions([0, 1])

    pdf_merged = pdf_12.copy()
    pdf_merged.append_variables_using_conditional_distributions(pdf_34_cond_12)

    assert pdf_merged == pdf_joint

    ### add a single variable conditioned on only the second existing variable

    rand_partial_cond_pdf = {(val,): JointProbabilityMatrix(1, pdf_joint.numvalues)
                             for val in xrange(pdf_joint.numvalues)}

    pdf_orig = pdf_joint.copy()

    pdf_joint.append_variables_using_conditional_distributions(rand_partial_cond_pdf, [1])

    assert pdf_orig == pdf_joint.marginalize_distribution(range(len(pdf_orig)))

    assert pdf_joint.mutual_information([0], [len(pdf_joint) - 1]) < 0.01, 'should be close to 0, not conditioned on'
    assert pdf_joint.mutual_information([1], [len(pdf_joint) - 1]) > 0.0, 'should not be 0, it is conditioned on'
    assert pdf_joint.mutual_information([2], [len(pdf_joint) - 1]) < 0.01, 'should be close to 0, not conditioned on'
    assert pdf_joint.mutual_information([3], [len(pdf_joint) - 1]) < 0.01, 'should be close to 0, not conditioned on'


def run_append_conditional_entropy_test():
    pdf_joint = JointProbabilityMatrix(4, 3)

    assert pdf_joint.conditional_entropy([1,2]) >= pdf_joint.conditional_entropy([1])
    assert pdf_joint.conditional_entropy([1,0]) >= pdf_joint.conditional_entropy([1])

    assert pdf_joint.conditional_entropy([0,2]) <= pdf_joint.entropy([0,2])

    assert pdf_joint.entropy([]) == 0, 'H(<empty-set>)=0 ... right? Yes think so'

    np.testing.assert_almost_equal(pdf_joint.conditional_entropy([1,2], [1,2]), 0.0)
    np.testing.assert_almost_equal(pdf_joint.entropy([0]) + pdf_joint.conditional_entropy([1,2,3], [0]),
                                   pdf_joint.entropy([0,1,2,3]))
    np.testing.assert_almost_equal(pdf_joint.conditional_entropy([0,1,2,3]), pdf_joint.entropy())


def run_params2matrix_incremental_test(numvars=3):
    pdf1 = JointProbabilityMatrix(numvars, 3)
    pdf2 = JointProbabilityMatrix(numvars, 3)

    params1 = pdf1.matrix2params_incremental(return_flattened=True)
    tree1 = pdf1.matrix2params_incremental(return_flattened=False)
    tree11 = pdf1.imbalanced_tree_from_scalars(params1, pdf1.numvalues)
    params1_from_tree1 = pdf1.scalars_up_to_level(tree1)
    params1_from_tree11 = pdf1.scalars_up_to_level(tree11)

    np.testing.assert_array_almost_equal(params1, params1_from_tree11)  # more a test of tree conversion itself
    np.testing.assert_array_almost_equal(params1, params1_from_tree1)

    pdf2.params2matrix_incremental(params1)

    params2 = pdf2.matrix2params_incremental()

    assert pdf1 == pdf2, 'computing parameter values from joint pdf and using those to construct a 2nd joint pdf ' \
                         'should result in two equal pdfs.\nparams1 = ' + str(params1) + '\nparms2 = ' + str(params2)

    pdf2.params2matrix_incremental(pdf2.matrix2params_incremental())

    assert pdf1 == pdf2, 'computing parameter values from joint pdf and using those to reconstruct the joint pdf ' \
                         'should result in two equal pdfs.'

    ### TEST the incrementality of the parameters

    pdf_marginal = pdf1.marginalize_distribution([0])
    params_marginal = pdf_marginal.matrix2params_incremental()
    np.testing.assert_array_almost_equal(params_marginal, pdf1.matrix2params_incremental()[:len(params_marginal)])

    pdf_marginal = pdf1.marginalize_distribution([0, 1])
    params_marginal = pdf_marginal.matrix2params_incremental()
    try:
        np.testing.assert_array_almost_equal(flatten(params_marginal),
                                             flatten(pdf1.matrix2params_incremental()[:len(params_marginal)]))
    except AssertionError as e:
        print '---------------------'
        print 'debug: params_marginal =                 ', np.round(params_marginal, decimals=4)
        print 'debug: pdf1.matrix2params_incremental() =', np.round(pdf1.matrix2params_incremental(), 4)
        print '---------------------'
        print 'debug: params_marginal =                 ', \
            pdf_marginal.matrix2params_incremental(return_flattened=False)
        print 'debug: pdf1.matrix2params_incremental() =', \
            pdf1.matrix2params_incremental(return_flattened=False)
        print '---------------------'

        raise AssertionError(e)

    if numvars >= 3:
        pdf_marginal = pdf1.marginalize_distribution([0, 1, 2])
        params_marginal = pdf_marginal.matrix2params_incremental()
        np.testing.assert_array_almost_equal(flatten(params_marginal),
                                             flatten(pdf1.matrix2params_incremental()[:len(params_marginal)]))


def run_scalars_to_tree_test():
    pdf = JointProbabilityMatrix(4, 3)

    list_of_scalars = pdf.matrix2params()  # does not matter what sequence of numbers, as long as the length is correct

    tree = pdf.imbalanced_tree_from_scalars(list_of_scalars, pdf.numvalues)
    np.testing.assert_array_almost_equal(pdf.scalars_up_to_level(tree), list_of_scalars)
    # np.testing.assert_array_almost_equal(pdf.scalars_up_to_level(tree), pdf.matrix2params_incremental())

    # another tree
    tree = pdf.matrix2params_incremental(return_flattened=False)

    # assert not np.isscalar(tree[-1]), 'you changed the matrix2params_incremental to return flat lists, which is good ' \
    #                                   'but then you should change this test, set an argument like flatten=False?'

    list_of_scalars2 = pdf.scalars_up_to_level(tree)

    tree2 = pdf.imbalanced_tree_from_scalars(list_of_scalars2, pdf.numvalues)

    np.testing.assert_array_almost_equal(flatten(tree), flatten(tree2))
    np.testing.assert_array_almost_equal(pdf.scalars_up_to_level(tree), pdf.scalars_up_to_level(tree2))


def run_all_tests(verbose=True, all_inclusive=False):
    run_append_and_marginalize_test()
    if verbose:
        print 'note: test run_append_and_marginalize_test successful.'

    run_params2matrix_test()
    if verbose:
        print 'note: test run_params2matrix_test successful.'

    run_vector2matrix_test()
    if verbose:
        print 'note: test run_vector2matrix_test successful.'

    run_conditional_pdf_test()
    if verbose:
        print 'note: test run_conditional_pdf_test successful.'

    run_append_using_transitions_table_and_marginalize_test()
    if verbose:
        print 'note: test run_append_using_transitions_table_and_marginalize_test successful.'

    run_append_conditional_pdf_test()
    if verbose:
        print 'note: test run_append_conditional_pdf_test successful.'

    run_scalars_to_tree_test()
    if verbose:
        print 'note: test run_scalars_to_tree_test successful.'

    run_params2matrix_incremental_test()
    if verbose:
        print 'note: test run_params2matrix_incremental_test successful.'

    run_synergistic_variables_test()
    if verbose:
        print 'note: test run_synergistic_variables_test successful.'

    run_synergistic_variables_test_with_subjects()
    if verbose:
        print 'note: test run_synergistic_variables_test_with_subjects successful.'

    run_synergistic_variables_test_with_subjects_and_agnostics()
    if verbose:
        print 'note: test run_synergistic_variables_test_with_subjects_and_agnostics successful.'

    run_append_conditional_entropy_test()
    if verbose:
        print 'note: test run_append_conditional_entropy_test successful.'

    run_reorder_test()
    if verbose:
        print 'note: test run_reorder_test successful.'

    if all_inclusive:
        run_orthogonalization_test()
        if verbose:
            print 'note: test run_orthogonalization_test successful.'

    run_orthogonalization_test_null_hypothesis()
    if verbose:
        print 'note: test run_orthogonalization_test_null_hypothesis successful.'

    if verbose:
        print 'note: finished. all tests successful.'


# def sum_modulo(values, modulo):  # deprecated, can be removed?
#     """
#     An example function which can be passed as the state_transitions argument to the
#     append_variables_using_state_transitions_table function of JointProbabilityMatrix.
#     :rtype : int
#     :param values: list of values which should be in integer range [0, modulo)
#     :param modulo: value self.numvalues usually, e.g. 2 for binary values
#     :return: for binary variables it is the XOR, for others it is summed modulo of integers.
#     """
#     return int(np.mod(np.sum(values), modulo))


### START of testing functions which try to determine if the synergy variables implementation behaves as expected,
### and for which I would like to see the results e.g. to plot in the paper



# returned by test_upper_bound_single_srv_entropy()
class TestUpperBoundSingleSRVEntropyResult():

    def __init__(self):
        self.num_subject_variables = None
        self.num_synergistic_variables = None
        self.num_values = None
        self.theoretical_upper_bounds = []
        self.entropies_srv = []  # actually lower bounds (same as entropies_lowerbound_srv), namely I(X:SRV) - sum_i I(X_i:SRV)
        self.pdfs_with_srv = []
        self.rel_errors_srv = []  # sum_indiv_mis / total_mi
        self.entropies_lowerbound_srv = []  # I(X:SRV) - sum_i I(X_i:SRV)
        self.entropies_upperbound_srv = []  # I(X:SRV) - max_i I(X_i:SRV)


    # parameters used to produce the results
    num_subject_variables = None
    num_synergistic_variables = None
    num_values = None

    # list of floats, each i'th number should ideally be >= to the i'th element of entropies_srv. These are the
    # theoretically predicted upper bounds of any SRV from the section "Consequential properties" in the synergy
    # paper
    theoretical_upper_bounds = []

    # list of floats, each i'th number is the estimated maximum entropy of a single SRV
    entropies_srv = []  # actually lower bounds (same as entropies_lowerbound_srv), namely I(X:SRV) - sum_i I(X_i:SRV)
    pdfs_with_srv = []
    rel_errors_srv = []  # sum_indiv_mis / total_mi

    entropies_lowerbound_srv = []  # I(X:SRV) - sum_i I(X_i:SRV)
    entropies_upperbound_srv = []  # I(X:SRV) - max_i I(X_i:SRV)

    def __str__(self):
        if len(self.entropies_srv) > 0:
            return '[H(S(x))=' + str(np.mean(self.entropies_srv)) \
                   + '=' + str(np.mean(self.entropies_srv)/np.mean(self.theoretical_upper_bounds)) + 'H_syn(X), +/- ' \
                   + str(np.mean(self.rel_errors_srv)*100.0) + '%]'
        else:
            return '[H(S(x))=nan' \
                   + '=(nan)H_syn(X), +/- ' \
                   + 'nan' + '%]'


def test_upper_bound_single_srv_entropy(num_subject_variables=2, num_synergistic_variables=2, num_values=2,
                                        num_samples=10, tol_rel_err=0.05, verbose=True, num_repeats_per_sample=5):
    """
    Measure the entropy a single SRV and compare it to the theoretical upper bound derived in the synergy paper, along
    with the relative error of the synergy estimation.

    Note: instead of the entropy of the SRV I use the MI of the SRV with the subject variables, because
    the optimization procedure in append_synergistic_variables does not (yet) try to simultaneously minimize
    the extraneous entropy in an SRV, which is entropy in an SRV which does not correlate with any other
    (subject) variable.

    :param num_subject_variables: number of X_i variables to compute an SRV for
    :param num_synergistic_variables:
    :param num_values:
    :param tol_rel_err:
    :param num_samples: number of entropy values will be returned, one for each randomly generated pdf.
    :param verbose:
    :param num_repeats_per_sample: number of optimizations to perform for each sample, trying to get the best result.
    :return: object with a list of estimated SRV entropies and a list of theoretical upper bounds of this same entropy,
    each one for a different, randomly generated PDF.
    :rtype: TestUpperBoundSingleSRVEntropyResult
    """

    result = TestUpperBoundSingleSRVEntropyResult()  # object to hold results

    # parameters used to produce the results
    result.num_subject_variables = num_subject_variables
    result.num_synergistic_variables = num_synergistic_variables
    result.num_values = num_values

    # shorthands
    synergistic_variables = range(num_subject_variables, num_subject_variables + num_synergistic_variables)
    subject_variables = range(num_subject_variables)

    time_before = time.time()

    # generate samples
    for trial in xrange(num_samples):
        pdf = JointProbabilityMatrix(num_subject_variables, numvalues=num_values)

        theoretical_upper_bound = pdf.entropy() - max([pdf.entropy([si]) for si in subject_variables])

        pdf.append_synergistic_variables(num_synergistic_variables, num_repeats=num_repeats_per_sample,
                                         verbose=verbose)

        # prevent double computations:
        indiv_mis = [pdf.mutual_information(synergistic_variables, [si]) for si in subject_variables]
        sum_indiv_mis = sum(indiv_mis)
        total_mi = pdf.mutual_information(synergistic_variables, subject_variables)

        # an error measure of the optimization procedure of finding synergistic variables, potentially it can be quite
        # bad, this one is normalized in [0,1]
        rel_error_srv = (sum_indiv_mis - max(indiv_mis)) / total_mi

        if rel_error_srv <= tol_rel_err:
            # note: instead of the entropy o the SRV I use the MI of the SRV with the subject variables, because
            # the optimization procedure in append_synergistic_variables does not (yet) try to simultaneously minimize
            # the extraneous entropy in an SRV, which is entropy in an SRV which does not correlate with any other
            # (subject) variable.
            # note: subtracting the sum_indiv_mi is done because this part of the entropy of the SRV is not
            # synergistic, so it would overestimate the synergistic entropy. But it is also an estimate, only valid
            # in case the non-synergistic MIs I_indiv(syns : subs) indeed factorizes into the sum, occurring e.g.
            # if the individual subject variables are independent, or if the MI with each subject variable is
            # non-redundant with the MIs with other subject variables.
            # entropy_srv = pdf.entropy(synergistic_variables)

            result.entropies_lowerbound_srv.append(total_mi - sum_indiv_mis)
            result.entropies_upperbound_srv.append(total_mi - max(indiv_mis))
            # take the best estimate for synergistic entropy to be the middle-way between lb and ub
            entropy_srv = (result.entropies_upperbound_srv[-1] + result.entropies_lowerbound_srv[-1]) / 2.0
            result.theoretical_upper_bounds.append(theoretical_upper_bound)
            result.entropies_srv.append(entropy_srv)
            result.pdfs_with_srv.append(pdf)
            result.rel_errors_srv.append(rel_error_srv)

            if verbose:
                print 'note: added a new sample, entropy_srv =', entropy_srv, ' and theoretical_upper_bound =', \
                    theoretical_upper_bound, ', and rel_error_srv =', rel_error_srv
        else:
            if verbose:
                print 'note: trial #' + str(trial) + ' will be discarded because the relative error of the SRV found' \
                      + ' is ' + str(rel_error_srv) + ' which exceeds ' + str(tol_rel_err) \
                      + ' (' + str(time.time() - time_before) + 's elapsed)'

    return result


def test_upper_bound_single_srv_entropy_many(num_subject_variables_list=list([2,3]),
                                             num_synergistic_variables_list=list([1,2]),
                                             num_values_list=list([2,3]),
                                             num_samples=10, tol_rel_err=0.05, verbose=True, num_repeats_per_sample=5):

    results_dict = dict()

    time_before = time.time()

    total_num_its = len(num_subject_variables_list) * len(num_synergistic_variables_list) * len(num_values_list)
    num_finished_its = 0

    for num_subject_variables in num_subject_variables_list:
        results_dict[num_subject_variables] = dict()

        for num_synergistic_variables in num_synergistic_variables_list:
            results_dict[num_subject_variables][num_synergistic_variables] = dict()

            for num_values in num_values_list:
                resobj = test_upper_bound_single_srv_entropy(num_subject_variables=num_subject_variables,
                                                             num_synergistic_variables=num_synergistic_variables,
                                                             num_values=num_values,
                                                             num_samples=num_samples,
                                                             tol_rel_err=tol_rel_err,
                                                             verbose=int(verbose)-1,
                                                             num_repeats_per_sample=num_repeats_per_sample)

                results_dict[num_subject_variables][num_synergistic_variables][num_values] = copy.deepcopy(resobj)

                num_finished_its += 1

                if verbose:
                    print 'note: finished', num_finished_its, '/', total_num_its, 'iterations after', \
                        time.time() - time_before, 'seconds. Result for (nx,ns,nv)=' \
                        + str((num_subject_variables, num_synergistic_variables, num_values)) + ': ', resobj

                del resobj  # looking for some inadvertent use of pointer

    return results_dict


class TestSynergyInRandomPdfs():
    syn_info_list = []
    total_mi_list = []
    indiv_mi_list_list = []  # list of list
    susceptibility_global_list = []
    susceptibilities_local_list = []  # list of list
    pdf_XY_list = []


def test_synergy_in_random_pdfs(num_variables_X, num_variables_Y, num_values,
                                num_samples=10, tolerance_nonsyn_mi=0.05, verbose=True, minimize_method=None,
                                perturbation_size=0.1, num_repeats_per_srv_append=3):

    """


    :param minimize_method: the default is chosen if None, which is good, but not necessarily fast. One other, faster
    option I found was "SLSQP", but I have not quantified rigorously how much better/worse in terms of accuracy it is.
    :param num_variables_X:
    :param num_variables_Y:
    :param num_values:
    :param num_samples:
    :param tolerance_nonsyn_mi:
    :param verbose:
    :rtype: TestSynergyInRandomPdfs
    """

    result = TestSynergyInRandomPdfs()

    time_before = time.time()

    try:
        for i in xrange(num_samples):
            pdf = JointProbabilityMatrix(num_variables_X + num_variables_Y, num_values)

            result.pdf_XY_list.append(pdf)

            result.syn_info_list.append(pdf.synergistic_information(range(num_variables_X, num_variables_X + num_variables_Y),
                                                          range(num_variables_X),
                                                          tolerance_nonsyn_mi=tolerance_nonsyn_mi,
                                                          verbose=bool(int(verbose)-1),
                                                          minimize_method=minimize_method,
                                                          num_repeats_per_srv_append=num_repeats_per_srv_append))

            result.total_mi_list.append(pdf.mutual_information(range(num_variables_X, num_variables_X + num_variables_Y),
                                                          range(num_variables_X)))

            indiv_mi_list = [pdf.mutual_information(range(num_variables_X, num_variables_X + num_variables_Y),
                                                          [xid]) for xid in xrange(num_variables_X)]

            result.indiv_mi_list_list.append(indiv_mi_list)

            result.susceptibility_global_list.append(pdf.susceptibility_global(num_variables_Y, ntrials=50,
                                                                               perturbation_size=perturbation_size))

            result.susceptibilities_local_list.append(pdf.susceptibilities_local(num_variables_Y, ntrials=50,
                                                                                 perturbation_size=perturbation_size))

            if verbose:
                print 'note: finished sample', i, 'of', num_samples, ', syn. info. =', result.syn_info_list[-1], 'after', \
                    time.time() - time_before, 'seconds'
    except KeyboardInterrupt as e:
        min_len = len(result.susceptibilities_local_list)  # last thing to append to in above loop, so min. length

        result.syn_info_list = result.syn_info_list[:min_len]
        result.total_mi_list = result.total_mi_list[:min_len]
        result.indiv_mi_list_list = result.indiv_mi_list_list[:min_len]
        result.susceptibility_global_list = result.susceptibility_global_list[:min_len]

        if verbose:
            print 'note: keyboard interrupt. Will stop the loop and return the', min_len, 'result I have so far.', \
                'Took me', time.time() - time_before, 'seconds.'

    return result


# returned by test_accuracy_orthogonalization()
class TestAccuracyOrthogonalizationResult():
    # parameters used to produce the results
    num_subject_variables = None
    num_orthogonalized_variables = None
    num_values = None

    # list of floats, each i'th number is the estimated entropy of the <num_subject_variables> variables
    actual_entropies_subject_variables = []
    actual_entropies_orthogonalized_variables = []
    actual_entropies_parallel_variables = []
    actual_entropies_orthogonal_variables = []

    # list of JointProbabilityMatrix objects
    joint_pdfs = []
    # list of floats, each i'th number is the relative error of the orthogonalization, in range 0..1
    rel_errors = []


def test_accuracy_orthogonalization(num_subject_variables=2, num_orthogonalized_variables=2, num_values=2,
                                        num_samples=10, verbose=True, num_repeats=5):
    """


    Note: instead of the entropy o the SRV I use the MI of the SRV with the subject variables, because
    the optimization procedure in append_synergistic_variables does not (yet) try to simultaneously minimize
    the extraneous entropy in an SRV, which is entropy in an SRV which does not correlate with any other
    (subject) variable.

    :param num_subject_variables:
    :param num_synergistic_variables:
    :param num_values:
    :param tol_rel_err:
    :return: object with a list of estimated SRV entropies and a list of theoretical upper bounds of this same entropy,
    each one for a different, randomly generated PDF.
    :rtype: TestUpperBoundSingleSRVEntropyResult
    """

    # todo: also include the option of known null hypothesis, to see what is the error if it is already known that 0.0
    # error is in fact possible

    # take these to be equal to num_orthogonalized_variables just to be sure
    num_orthogonal_variables = num_orthogonalized_variables
    num_parallel_variables = num_orthogonalized_variables

    result = TestAccuracyOrthogonalizationResult()  # object to hold results

    # parameters used to produce the results
    result.num_subject_variables = num_subject_variables
    result.num_synergistic_variables = num_orthogonalized_variables
    result.num_values = num_values

    # shorthands for variable index ranges for the different roles of variables
    orthogonalized_variables = range(num_subject_variables, num_subject_variables + num_orthogonalized_variables)
    subject_variables = range(num_subject_variables)
    orthogonal_variables = range(num_subject_variables + num_orthogonalized_variables,
                                 num_subject_variables + num_orthogonalized_variables + num_orthogonal_variables)
    parallel_variables = range(num_subject_variables + num_orthogonalized_variables + num_orthogonal_variables,
                               num_subject_variables + num_orthogonalized_variables + num_orthogonal_variables
                               + num_parallel_variables)

    time_before = time.time()

    # generate samples
    for trial in xrange(num_samples):
        pdf = JointProbabilityMatrix(num_subject_variables + num_orthogonalized_variables, numvalues=num_values)

        pdf.append_orthogonalized_variables(orthogonalized_variables, num_orthogonal_variables, num_parallel_variables,
                                            verbose=verbose, num_repeats=num_repeats)

        assert len(pdf) == num_subject_variables + num_orthogonalized_variables + num_orthogonal_variables \
                           + num_parallel_variables

        result.joint_pdfs.append(pdf)

        result.actual_entropies_subject_variables.append(pdf.entropy(subject_variables))
        result.actual_entropies_orthogonalized_variables.append(pdf.entropy(orthogonalized_variables))
        result.actual_entropies_orthogonal_variables.append(pdf.entropy(orthogonal_variables))
        result.actual_entropies_parallel_variables.append(pdf.entropy(parallel_variables))

        rel_error_1 = (pdf.mutual_information(orthogonalized_variables, subject_variables) - 0) \
                      / result.actual_entropies_subject_variables[-1]
        rel_error_2 = (pdf.mutual_information(orthogonalized_variables, subject_variables)
                      - pdf.mutual_information(parallel_variables, subject_variables)) \
                      / pdf.mutual_information(orthogonalized_variables, subject_variables)
        rel_error_3 = (result.actual_entropies_orthogonalized_variables[-1]
                      - pdf.mutual_information(orthogonal_variables + parallel_variables, orthogonalized_variables)) \
                      / result.actual_entropies_orthogonalized_variables[-1]

        # note: each relative error term is intended to be in range [0, 1]

        rel_error = rel_error_1 + rel_error_2 + rel_error_3

        # max_rel_error = result.actual_entropies_subject_variables[-1] \
        #                 + pdf.mutual_information(orthogonalized_variables, subject_variables) \
        #                 + result.actual_entropies_orthogonalized_variables[-1]
        max_rel_error = 3.0

        rel_error /= max_rel_error

        result.rel_errors.append(rel_error)

        print 'note: finished trial #' + str(trial) + ' after', time.time() - time_before, 'seconds. rel_error=', \
            rel_error, '(' + str((rel_error_1, rel_error_2, rel_error_3)) + ')'

    return result


# todo: make a test function which iteratively appends a set of SRVs until no more can be found (some cut-off MI),
# then I am interested in the statistics of the number of SRVs that are found, as well as the statistics of whether
# the SRVs are statistically independent or not, and whether they have synergy about each other (if agnostic_about
# passed to append_syn* is a list of lists) and whether indeed every two SRVs necessarily will have info about an
# individual Y_i. In the end, if the assumption of perfect orthogonal decomposition of any RVs too strict, or is
# is often not needed at all? Because it would only be needed in case there are at least two SRVs, the combination of
# which (combined RV) has more synergy about Y than the two individually summed, but also the combined RV has
# info about individual Y_i so it is not an SRV by itself, so the synergy contains both synergistic and individual info.

def test_num_srvs(num_subject_variables=2, num_synergistic_variables=2, num_values=2, num_samples=10,
                  verbose=True, num_repeats=5):
    assert False