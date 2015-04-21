__author__ = 'rquax'


import numpy as np
import itertools
import copy

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

    def joint_probability(self, values):
        assert len(values) == self.numvariables, 'should specify one value per variable'
        assert values[0] < self.numvalues, 'variable can only take values 0, 1, ..., <numvalues - 1>'

        joint_prob = self.joint_probabilities[tuple(values)]

        assert 0.0 <= joint_prob <= 1.0, 'not a probability?'

        return joint_prob

    def marginalize_distribution(self, variables):
        """

        :rtype : JointProbabilityMatrix
        """
        lists_of_possible_states_per_variable = [range(self.numvalues) for variable in xrange(self.numvariables)]

        marginalized_joint_probs = np.zeros([self.numvalues]*len(variables))

        # if len(variables):
        #     marginalized_joint_probs = np.array([marginalized_joint_probs])

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
        assert len(new_joint_pdf[0]) == self.numvalues
        assert len(new_joint_pdf[-1]) == self.numvalues
        np.testing.assert_almost_equal(np.sum(new_joint_pdf), 1.0)

        self.joint_probabilities = new_joint_pdf
        self.numvariables = self.numvariables + num_added_variables

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

    def matrix2vector(self):
        return [i for i in self.joint_probabilities.flat]

    def vector2matrix(self, list_probs):
        np.testing.assert_almost_equal(np.sum(list_probs), 1.0)

        assert np.ndim(list_probs) == 1

        self.joint_probabilities = np.reshape(list_probs, [self.numvalues]*self.numvariables)

    def params2matrix(self, parameters):
        assert len(parameters) == np.power(self.numvalues, self.numvariables) - 1

        vector_probs = [-1.0]*(np.power(self.numvalues, self.numvariables))

        remaining_prob_mass = 1.0

        for pix in xrange(len(parameters)):
            assert 0.0 <= parameters[pix] <= 1.0, 'parameters should be in [0, 1]'

            vector_probs[pix] = remaining_prob_mass * parameters[pix]

            remaining_prob_mass = remaining_prob_mass * (1.0 - parameters[pix])

        assert vector_probs[-1] < 0.0, 'should still be unset by the above loop'

        # last parameter is irrelevant, must always be 1.0 is also a way to look at it
        vector_probs[-1] = remaining_prob_mass

        np.testing.assert_almost_equal(np.sum(vector_probs), 1.0)

        self.vector2matrix(vector_probs)

    def matrix2params(self):
        vector_probs = self.matrix2vector()

        remaining_prob_mass = 1.0

        parameters = [-1.0]*(len(vector_probs) - 1)

        for pix in xrange(len(parameters)):
            parameters[pix] = vector_probs[pix] / remaining_prob_mass

            remaining_prob_mass = remaining_prob_mass * (1.0 - parameters[pix])

            assert 0.0 <= parameters[pix] <= 1.0, 'parameters should be in [0, 1]'

        return parameters

    def append_variables_using_state_transitions_table(self, state_transitions):
        """

        :param state_transitions: list of lists, where each sublist is at least of length self.numvariables + 1 and
        is a list of values, each value in [0, numvalues).
        :type state_transitions: list
        """

        # test

    def entropy(self):
        def log2term(p):  # p log_2 p
            if 0.0 < p < 1.0:
                return p * np.log2(p)
            else:
                return 0.0  # 0 log 0 == 0 is assumed

        joint_entropy = np.sum(map(log2term, self.joint_probabilities.flat))

        assert joint_entropy >= 0.0
        assert joint_entropy <= np.log2(self.numvalues) * self.numvariables

        return joint_entropy

    def conditional_probability_distribution(self, given_variables, given_values):
        assert len(given_values) == len(given_variables)
        assert len(given_variables) < self.numvariables, 'no variables left after conditioning'

        lists_of_possible_states_per_variable = [range(self.numvalues) for variable in xrange(self.numvariables)]

        # overwrite the 'state spaces' for the specified variables, to the specified state spaces
        for gix in xrange(len(given_variables)):
            assert np.isscalar(given_values[gix]), 'assuming specific value, not list of possibilities'

            lists_of_possible_states_per_variable[given_variables[gix]] = \
                [given_values[gix]] if np.isscalar(given_values[gix]) else given_values[gix]

        variables = [varix for varix in xrange(self.numvariables) if not varix in given_variables]

        conditional_joint_probs = np.zeros([self.numvalues]*len(variables))

        for values in itertools.product(*lists_of_possible_states_per_variable):
            marginal_values = [values[varid] for varid in variables]

            conditional_joint_probs[tuple(marginal_values)] += self.joint_probability(values)

        if __debug__:
            pass  # todo: test here if the summed prob. mass equals the marginal prob of the given variable values

        conditional_joint_probs /= np.sum(conditional_joint_probs)

        conditional_joint_pdf = JointProbabilityMatrix(len(variables), self.numvalues,
                                                    joint_probs=conditional_joint_probs)

        return conditional_joint_pdf

    def conditional_probability_distributions(self, given_variables):
        assert len(given_variables) < self.numvariables, 'no variables left after conditioning'

        lists_of_possible_given_values = [range(self.numvalues) for variable in xrange(len(given_variables))]

        return {tuple(given_values): self.conditional_probability_distribution(given_variables=given_variables,
                                                                               given_values=given_values)
                for given_values in itertools.product(*lists_of_possible_given_values)}


### UNIT TESTING:


def run_append_and_marginalize_test():
    pdf = JointProbabilityMatrix(3, 2)

    pdf_copy = pdf.copy()

    pdf.append_variables(4)

    pdf_old = pdf.marginalize_distribution(range(3))

    assert pdf_copy == pdf_old, 'adding and then removing variables should result in the same joint pdf'


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
