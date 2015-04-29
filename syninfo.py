__author__ = 'rquax'


import numpy as np
import itertools
import copy
from scipy.optimize import minimize
import warnings
from compiler.ast import flatten  # note: deprecated in Python 3, in which case: find another flatten
from collections import Sequence


def maximum_depth(seq):
    """
    Helper function, e.g. depth([1,2,[2,4,[[4]]]]) == 4.
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


class JointProbabilityMatrix():
    # nested list of probabilities. For instance for three binary variables it could be:
    # [[[0.15999394, 0.06049343], [0.1013956, 0.15473886]], [[ 0.1945649, 0.15122334], [0.11951818, 0.05807175]]].
    # The depth (np.ndim) should be numvariables, and the length at each dimension should be numvalues.
    joint_probabilities = []

    numvariables = 0
    numvalues = 0

    def __init__(self, numvariables, numvalues, joint_probs=None):
        self.numvariables = numvariables
        self.numvalues = numvalues

        if joint_probs is None:
            self.generate_random_joint_probabilities()
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

    def copy(self):  # deep copy
        """

        :rtype : JointProbabilityMatrix
        """
        return copy.deepcopy(self)

    def generate_random_joint_probabilities(self):
        self.joint_probabilities = np.random.random([self.numvalues]*self.numvariables)
        self.joint_probabilities /= np.sum(self.joint_probabilities)

    def random_samples(self, n=1):
        # sample_indices = np.random.multinomial(n, [i for i in self.joint_probabilities.flat])

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
        Joint probability of a list of values.
        :param values: list of values, each value is an integer in [0, numvalues)
        :type values: list
        :return: joint probability of the given values, in order, for all variables; in [0, 1]
        :rtype: float
        """
        return self.joint_probability(values=values)

    def clip_all_probabilities(self):
        """
        Make sure all probabilities in the joint probability matrix are in the range [0.0, 1.0], which could be
        violated sometimes due to floating point operation roundoff errors.
        """
        self.joint_probabilities = np.minimum(np.maximum(self.joint_probabilities, 0.0), 1.0)

        try:
            np.testing.assert_almost_equal(np.sum(self.joint_probabilities), 1.0)
        except AssertionError as e:
            print 'error message: ' + str(e)

            print 'error: len(self.joint_probabilities) =', len(self.joint_probabilities)
            print 'error: shape(self.joint_probabilities) =', np.shape(self.joint_probabilities)
            if len(self.joint_probabilities) < 30:
                print 'error: self.joint_probabilities =', self.joint_probabilities

            raise AssertionError(e)

    def joint_probability(self, values):
        assert len(values) == self.numvariables, 'should specify one value per variable'
        assert values[0] < self.numvalues, 'variable can only take values 0, 1, ..., <numvalues - 1>: ' + str(values[0])
        assert values[-1] < self.numvalues, 'variable can only take values 0, 1, ..., <numvalues - 1>: ' \
                                            + str(values[-1])

        joint_prob = self.joint_probabilities[tuple(values)]

        assert 0.0 <= joint_prob <= 1.0, 'not a probability? ' + str(joint_prob)

        return joint_prob

    def marginalize_distribution(self, variables):
        """

        :rtype : JointProbabilityMatrix
        """
        lists_of_possible_states_per_variable = [range(self.numvalues) for variable in xrange(self.numvariables)]

        marginalized_joint_probs = np.zeros([self.numvalues]*len(variables))

        # if len(variables):
        #     marginalized_joint_probs = np.array([marginalized_joint_probs])

        assert len(variables) > 0, 'makes no sense to marginalize 0 variables'
        assert np.all(map(np.isscalar, variables)), 'each variable identifier should be int in [0, numvalues)'
        assert len(variables) <= self.numvariables, 'cannot marginalize more variables than I have'
        assert len(set(variables)) <= self.numvariables, 'cannot marginalize more variables than I have'

        # not sure yet about doing this:
        # variables = sorted(list(set(variables)))  # ensure uniqueness?

        for values in itertools.product(*lists_of_possible_states_per_variable):
            marginal_values = [values[varid] for varid in variables]

            marginalized_joint_probs[tuple(marginal_values)] += self.joint_probability(values)

        np.testing.assert_almost_equal(np.sum(marginalized_joint_probs), 1.0)

        marginal_joint_pdf = JointProbabilityMatrix(len(variables), self.numvalues,
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

        Can also provide a function f(values, max_value) which returns a list of values for the to-be-appended
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


    def __len__(self):
        return self.numvariables


    def matrix2vector(self):
        return [i for i in self.joint_probabilities.flat]


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
            assert 0.0 <= parameters[pix] <= 1.000001, 'parameters should be in [0, 1]: ' + str(parameters[pix])

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


    def matrix2params(self):  # todo: sort the parameters such that when appending variables the first some is unchanged
        vector_probs = self.matrix2vector()

        remaining_prob_mass = 1.0

        parameters = [-1.0]*(len(vector_probs) - 1)

        for pix in xrange(len(parameters)):
            parameters[pix] = vector_probs[pix] / remaining_prob_mass

            remaining_prob_mass = remaining_prob_mass * (1.0 - parameters[pix])

            assert 0.0 <= parameters[pix] <= 1.000001, 'parameters should be in [0, 1]: ' + str(parameters[pix])

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

    def params2matrix_incremental(self, parameters):

        if __debug__:
            # store the original provided list of scalars
            original_parameters = list(parameters)

        # I suspect that both a tree-like input and a list of scalars should work... (add to unit test?)
        if np.all(map(np.isscalar, parameters)):
            parameters = self.imbalanced_tree_from_scalars(parameters, self.numvalues)

            # verify that the procedure to make the tree out of the list of scalars is reversible and correct
            # (looking for bug)
            if __debug__:
                original_parameters2 = self.scalars_up_to_level(parameters)

                np.testing.assert_array_almost_equal(original_parameters, original_parameters2)

        if self.numvariables > 1:
            # first (numvalues - 1) values in the parameters tree structure should be scalars, as they will be used
            # to make the first variable's marginal distribution
            assert np.all(map(np.isscalar, parameters[:(self.numvalues - 1)]))

            pdf_1 = JointProbabilityMatrix(1, self.numvalues)
            pdf_1.params2matrix(parameters[:(len(pdf_1.joint_probabilities.flat) - 1)])

            assert len(pdf_1.joint_probabilities.flat) == self.numvalues

            assert len(flatten(parameters)) == len(self.joint_probabilities.flat) - 1, \
                'more or fewer parameters than needed: ' \
                              'need ' + str(len(self.joint_probabilities.flat) - 1) + ', got ' + str(len(parameters)) \
                              + '; #vars, #vals = ' + str(self.numvariables) + ', ' + str(self.numvalues)

            # remove the used parameters from the list
            parameters = parameters[(len(pdf_1.joint_probabilities.flat) - 1):]
            assert len(parameters) == self.numvalues  # one subtree per conditional pdf

            pdf_conds = dict()

            for val in xrange(self.numvalues):
                # set this conditional pdf recursively as defined by the next sequence of parameters
                pdf_cond = JointProbabilityMatrix(self.numvariables - 1, self.numvalues)
                # CHANGE123

                assert not np.isscalar(parameters[0])

                # todo: changing the parameters list is not necessary, maybe faster if not?

                # pdf_cond.params2matrix_incremental(parameters[:(len(pdf_cond.joint_probabilities.flat) - 1)])
                pdf_cond.params2matrix_incremental(parameters[0])

                # conditional pdf should have the same set of parameters as the ones I used to create it
                # (todo: remove this expensive check if it seems to work for  while)
                np.testing.assert_array_almost_equal(pdf_cond.matrix2params_incremental(),
                                                     self.scalars_up_to_level(parameters[0]))

                # remove the used parameters from the list
                # CHANGE123
                # parameters = parameters[(len(pdf_cond.joint_probabilities.flat) - 1):]
                parameters = parameters[1:]

                # add the conditional pdf
                pdf_conds[(val,)] = pdf_cond.copy()

            assert len(parameters) == 0, 'all parameters should be used to construct joint pdf'

            pdf_1.append_variables_using_conditional_distributions(pdf_conds)

            if __debug__:
                # remove this (expensive) check after it seems to work a few times?
                np.testing.assert_array_almost_equal(pdf_1.matrix2params_incremental(),
                                                     self.scalars_up_to_level(original_parameters))

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

    def reset(self, numvariables, numvalues, joint_prob_matrix):
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

        assert np.ndim(joint_prob_matrix) == self.numvariables, 'ndim = ' + str(np.ndim(joint_prob_matrix)) + ', ' \
                                                                'self.numvariables = ' + str(self.numvariables) + ', ' \
                                                                'joint matrix = ' + str(joint_prob_matrix)
        assert list(set(np.shape(self.joint_probabilities)))[0] == self.numvalues

        # maybe this check should be removed, it is also checked in clip_all_* below, but after clipping, which
        # may be needed to get this condition valid again?
        np.testing.assert_array_almost_equal(np.sum(joint_prob_matrix), 1.0)

        self.clip_all_probabilities()

    def statespace(self):
        lists_of_possible_joint_values = [range(self.numvalues) for _ in xrange(self.numvariables)]

        return itertools.product(*lists_of_possible_joint_values)

    def append_variables_using_conditional_distributions(self, pdf_conds):
        """

        :param pdf_conds: dictionary of JointProbabilityMatrix objects, one for each possible set of values for the
        existing <self.numvariables> variables.
        :type pdf_conds: dict of JointProbabilityMatrix | JointProbabilityMatrix
        """

        if isinstance(pdf_conds, JointProbabilityMatrix):
            pdf_conds = {tuple(my_values): pdf_conds for my_values in self.statespace()}
        else:
            assert type(pdf_conds) == dict

        num_added_variables = pdf_conds[(0,)*self.numvariables].numvariables

        assert num_added_variables > 0, 'makes no sense to append 0 variables?'
        assert self.numvalues == pdf_conds[(0,)*self.numvariables].numvalues, 'added variables should have same #values'

        lists_of_possible_joint_values = [range(self.numvalues)
                                          for _ in xrange(self.numvariables + num_added_variables)]

        extended_joint_probs = np.zeros([self.numvalues]*(self.numvariables + num_added_variables))

        for values in itertools.product(*lists_of_possible_joint_values):
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

        np.testing.assert_almost_equal(np.sum(extended_joint_probs), 1.0)

        self.reset(self.numvariables + num_added_variables, self.numvalues, extended_joint_probs)

    def conditional_probability_distribution(self, given_variables, given_values):
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

        assert summed_prob_mass >= 0.0, 'probability mass cannot be negative'

        if summed_prob_mass > 0.0:
            conditional_probs /= summed_prob_mass
        else:
            warnings.warn('conditional_probability_distribution: summed_prob_mass == 0.0')

            conditional_probs += 1.0

            conditional_probs /= summed_prob_mass

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


    def entropy(self, variables=None):
        def log2term(p):  # p log_2 p
            if 0.0 < p < 1.0:
                return -p * np.log2(p)
            else:
                return 0.0  # 0 log 0 == 0 is assumed

        if variables is None:
            np.testing.assert_almost_equal(np.sum(self.joint_probabilities), 1.0)

            joint_entropy = np.sum(map(log2term, self.joint_probabilities.flat))

            assert joint_entropy >= 0.0
            assert joint_entropy <= np.log2(self.numvalues) * self.numvariables

            return joint_entropy
        else:
            assert hasattr(variables, '__iter__')

            if len(variables) == 0:  # hard-coded this because otherwise I have to support empty pdfs (len() = 0)
                return 0.0

            marginal_pdf = self.marginalize_distribution(variables=variables)

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


    def mutual_information(self, variables1, variables2):
        assert hasattr(variables1, '__iter__'), 'variables1 = ' + str(variables1)
        assert hasattr(variables2, '__iter__'), 'variables2 = ' + str(variables2)

        if len(variables1) == 0 or len(variables2) == 0:
            mi = 0
        elif np.equal(sorted(variables1), sorted(variables2)).all():
            mi = self.entropy(variables1)
        else:
            mi = self.entropy(variables1) + self.entropy(variables2) \
                 - self.entropy(list(set(list(variables1) + list(variables2))))

        assert np.isscalar(mi)
        assert np.isfinite(mi)
        assert mi >= 0, 'mutual information should be non-negative'

        return mi


    def synergistic_information_naive(self, variables_Y, variables_X):
        return self.mutual_information(list(variables_Y), list(variables_X)) \
               - sum([self.mutual_information(list(variables_Y), list([var_xi])) for var_xi in variables_X])


    # todo: add optional numvalues, so that the synergistic variables can have more possible values than the
    # current variables (then set all probabilities where the original variables exceed their original max to 0)
    # todo: first you should then probably implement a .increase_num_values(num) or so
    def append_synergistic_variables(self, num_synergistic_variables, initial_guess_summed_modulo=True, verbose=False):
        parameter_values_before = list(self.matrix2params_incremental())

        if __debug__:
            debug_params_before = copy.deepcopy(parameter_values_before)

        # a pdf with XORs as appended variables (often already MSRV for binary variables), good initial guess?
        pdf_with_srvs = self.copy()
        pdf_with_srvs.append_variables_using_state_transitions_table(
            state_transitions=lambda vals, mv: [int(np.mod(np.sum(vals), mv))]*num_synergistic_variables)

        assert pdf_with_srvs.numvariables == self.numvariables + num_synergistic_variables

        parameter_values_after = pdf_with_srvs.matrix2params_incremental()

        assert len(parameter_values_after) > len(parameter_values_before), 'should be additional free parameters'
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
            debug_before_syninfo = debug_pdf_with_srvs.synergistic_information_naive(variables_Y=range(self.numvariables,
                                                                                         pdf_with_srvs.numvariables),
                                                                               variables_X=range(self.numvariables))

        assert len(initial_guess) == num_free_parameters

        pdf_with_srvs_for_optimization = pdf_with_srvs.copy()
        # NOTE: this fitness function does not try to minimize the (extraneous) entropy of the SRVs
        def fitness_func(free_params, parameter_values_before=parameter_values_before):
            assert len(free_params) == num_free_parameters

            pdf_with_srvs_for_optimization.params2matrix_incremental(list(parameter_values_before) + list(free_params))

            assert pdf_with_srvs_for_optimization.numvariables == self.numvariables + num_synergistic_variables

            fitness = -pdf_with_srvs_for_optimization.synergistic_information_naive(variables_Y=range(self.numvariables,
                                  pdf_with_srvs_for_optimization.numvariables), variables_X=range(self.numvariables))

            assert np.isscalar(fitness)
            assert np.isfinite(fitness)

            return float(fitness)

        param_vectors_trace = []

        optres = minimize(fitness_func, initial_guess, bounds=[(0.0, 1.0)]*num_free_parameters,
                          callback=(lambda xv: param_vectors_trace.append(list(xv))) if verbose else None)

        assert len(optres.x) == num_free_parameters
        assert max(optres.x) <= 1.0, 'parameter bound violated'
        assert min(optres.x) >= 0.0, 'parameter bound violated'

        optimal_parameters_joint_pdf = list(parameter_values_before) + list(optres.x)

        pdf_with_srvs.params2matrix_incremental(optimal_parameters_joint_pdf)

        assert pdf_with_srvs.numvariables == self.numvariables + num_synergistic_variables

        if verbose:
            parameter_values_after2 = pdf_with_srvs.matrix2params_incremental()

            assert len(parameter_values_after2) > len(parameter_values_before), 'should be additional free parameters'
            # see if the first part of the new parameters is exactly the old parameters, so that I know that I only
            # have to optimize the latter part of parameter_values_after
            np.testing.assert_array_almost_equal(parameter_values_before,
                                                 parameter_values_after2[:len(parameter_values_before)])
            np.testing.assert_array_almost_equal(parameter_values_after2[len(parameter_values_before):],
                                                 optres.x)

            # store the synergistic information before the optimization procedure (after procedure should be higher...)
            debug_after_syninfo = pdf_with_srvs.synergistic_information_naive(variables_Y=range(self.numvariables,
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

            # fitness function should be the negative of the synergistic info right?
            np.testing.assert_almost_equal(debug_after_syninfo, -fitness_func(optres.x))

        if __debug__:
            # unchanged, not accidentally changed by passing it as reference? looking for bug
            np.testing.assert_array_almost_equal(debug_params_before, parameter_values_before)

        self.duplicate(pdf_with_srvs)


    def orthogonalize(self, variables):
        """
        Orthogonalize the given set of variables=X relative to the rest. I.e., decompose X into two parts {X1,X2}
        where I(X1:rest)=0 but I(X1:X)=H(X1) and I(X2:rest)=H(X2).
        :param variables: sequence of int
        :rtype: JointProbabilityMatrix
        """

        assert False  # todo

        # todo: verify that H(X1,X2)  == H(X) (approx.)

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


def run_find_synergy_test(numvars=2):
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


def run_append_conditional_pdf_test():
    pdf_joint = JointProbabilityMatrix(4, 3)

    pdf_12 = pdf_joint.marginalize_distribution([0, 1])
    pdf_34_cond_12 = pdf_joint.conditional_probability_distributions([0, 1])

    pdf_merged = pdf_12.copy()
    pdf_merged.append_variables_using_conditional_distributions(pdf_34_cond_12)

    assert pdf_merged == pdf_joint


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


def run_all_tests(verbose=True):
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

    run_find_synergy_test()
    if verbose:
        print 'note: test run_find_synergy_test successful.'

    run_append_conditional_entropy_test()
    if verbose:
        print 'note: test run_append_conditional_entropy_test successful.'

    if verbose:
        print 'note: finished. all tests successful'


def sum_modulo(values, modulo):  # deprecated, can be removed?
    """
    An example function which can be passed as the state_transitions argument to the
    append_variables_using_state_transitions_table function of JointProbabilityMatrix.
    :rtype : int
    :param values: list of values which should be in integer range [0, modulo)
    :param modulo: value self.numvalues usually, e.g. 2 for binary values
    :return: for binary variables it is the XOR, for others it is summed modulo of integers.
    """
    return int(np.mod(np.sum(values), modulo))