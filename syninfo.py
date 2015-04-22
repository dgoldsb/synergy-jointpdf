__author__ = 'rquax'


import numpy as np
import itertools
import copy
from scipy.optimize import minimize

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

    def joint_probability(self, values):
        assert len(values) == self.numvariables, 'should specify one value per variable'
        assert values[0] < self.numvalues, 'variable can only take values 0, 1, ..., <numvalues - 1>: ' + str(values[0])
        assert values[-1] < self.numvalues, 'variable can only take values 0, 1, ..., <numvalues - 1>: ' \
                                            + str(values[-1])

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
        if self.numvariables + num_added_variables >= 1:
            assert len(new_joint_pdf[0]) == self.numvalues
            assert len(new_joint_pdf[-1]) == self.numvalues
        np.testing.assert_almost_equal(np.sum(new_joint_pdf), 1.0)

        self.joint_probabilities = new_joint_pdf
        self.numvariables = self.numvariables + num_added_variables


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

        self.joint_probabilities = extended_joint_probs
        self.numvariables = len(state_transitions[0])


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

            assert 0.0 <= parameters[pix] <= 1.0, 'parameters should be in [0, 1]'

        return parameters


    def matrix2params_incremental(self):
        if self.numvariables > 1:
            # get the marginal pdf for the first variable
            pdf1 = self.marginalize_distribution([0])

            parameters = pdf1.matrix2params()  # first sequence of parameters, rest is added below here

            pdf_conds = self.conditional_probability_distributions([0])

            assert len(pdf_conds) == self.numvalues, 'should be one pdf for each value of first variable'

            for val in xrange(self.numvalues):
                pdf_cond = pdf_conds[tuple([val])]

                added_params = pdf_cond.matrix2params_incremental()

                parameters.extend(added_params)

            if __debug__:
                direct_params = self.matrix2params()

                assert len(direct_params) == len(parameters), 'two methods of parameterizing should be equal in length'

            return parameters
        elif self.numvariables == 1:
            return self.matrix2params()
        else:
            raise ValueError('no parameters for 0 variables')


    def params2matrix_incremental(self):
        assert False


    def append_variables_using_conditional_distributions(self, pdf_conds):

        """

        :param pdf_conds: dict of JointProbabilityMatrix
        :type pdf_conds: dict of JointProbabilityMatrix
        """

        num_added_variables = pdf_conds[(0,)*self.numvariables].numvariables

        assert num_added_variables > 0, 'makes no sense to append 0 variables?'
        assert self.numvalues == pdf_conds[(0,)*self.numvariables].numvalues

        lists_of_possible_joint_values = [range(self.numvalues)
                                          for _ in xrange(self.numvariables + num_added_variables)]

        extended_joint_probs = np.zeros([self.numvalues]*len(state_transitions[0]))


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
                return p * np.log2(p)
            else:
                return 0.0  # 0 log 0 == 0 is assumed

        if variables is None:
            joint_entropy = np.sum(map(log2term, self.joint_probabilities.flat))

            assert joint_entropy >= 0.0
            assert joint_entropy <= np.log2(self.numvalues) * self.numvariables

            return joint_entropy
        else:
            assert hasattr(variables, '__iter__')

            marginal_pdf = self.marginalize_distribution(variables=variables)

            return marginal_pdf.entropy()


    def mutual_information(self, variables1, variables2):
        assert hasattr(variables1, '__iter__')
        assert hasattr(variables2, '__iter__')

        if len(variables1) == 0 or len(variables2) == 0:
            mi = 0
        elif np.equal(variables1, variables2).all():
            mi = self.entropy(variables1)
        else:
            mi = self.entropy(variables1) + self.entropy(variables2) \
                 - self.entropy(list(set(list(variables1) + list(variables2))))

        assert np.isscalar(mi)
        assert np.isfinite(mi)
        assert mi >= 0, 'mutual information should be non-negative'

        return mi


    def synergistic_information_by_difference(self, variables_Y, variables_X):
        return self.mutual_information(variables_Y, variables_X) \
               - sum([self.mutual_information(variables_Y, var_xi) for var_xi in variables_X])


    # todo: add optional numvalues, so that the synergistic variables can have more possible values than the
    # current variables (then set all probabilities where the original variables exceed their original max to 0)
    # todo: first you should then probably implement a .increase_num_values(num) or so
    def append_synergistic_variables(self, num_synergistic_variables, initial_guess_summed_modulo=True):
        parameter_values_before = self.matrix2params()

        # a pdf with XORs as appended variables (often already MSRV for binary variables), good initial guess?
        pdf_with_srvs = self.copy()
        pdf_with_srvs.append_variables_using_state_transitions_table(
            state_transitions=lambda vals, mv: [int(np.mod(np.sum(vals), mv))]*num_synergistic_variables)

        assert pdf_with_srvs.numvariables == self.numvariables + num_synergistic_variables

        parameter_values_after = pdf_with_srvs.matrix2params()

        assert len(parameter_values_after) > len(parameter_values_before), 'should be additional free parameters'

        # this many parameters (each in [0,1]) must be optimized
        num_free_parameters = len(parameter_values_after) - len(parameter_values_before)

        if initial_guess_summed_modulo:
            # note: this is xor for binary variables
            initial_guess = parameter_values_after[len(parameter_values_before):]
        else:
            initial_guess = np.random.random(num_free_parameters)

        assert len(initial_guess) == num_free_parameters

        # NOTE: this fitness function does not try to minimize the (extraneous) entropy of the SRVs
        def fitness_func(free_params):
            assert len(free_params) == num_free_parameters

            pdf_with_srvs.params2matrix(parameter_values_before + list(free_params))

            fitness = -pdf_with_srvs.synergistic_information_by_difference(variables_Y=range(self.numvariables,
                                                                                          pdf_with_srvs.numvariables),
                                                                           variables_X=range(self.numvariables))

            assert np.isscalar(fitness)
            assert np.isfinite(fitness)

            return float(fitness)

        optres = minimize(fitness_func, initial_guess, bounds=[(0.0, 1.0)]*num_free_parameters,
                          options=None)

        assert len(optres.x) == num_free_parameters
        assert max(optres.x) <= 1.0, 'parameter bound violated'
        assert min(optres.x) >= 0.0, 'parameter bound violated'

        optimal_parameters_joint_pdf = list(parameter_values_before) + list(optres.x)

        self.numvariables = pdf_with_srvs.numvariables
        self.numvalues = pdf_with_srvs.numvalues
        self.params2matrix(optimal_parameters_joint_pdf)
        # todo: make a function which resets a joint pdf so that the above lines go in there, and if variables are
        # added then it is detected



### UNIT TESTING:



def run_append_and_marginalize_test():
    pdf = JointProbabilityMatrix(3, 2)

    pdf_copy = pdf.copy()

    pdf.append_variables(4)

    pdf_old = pdf.marginalize_distribution(range(3))

    assert pdf_copy == pdf_old, 'adding and then removing variables should result in the same joint pdf'

    old_params = pdf_copy.matrix2params_incremental()

    assert np.testing.assert_array_almost_equal(pdf.matrix2params_incremental()[:len(old_params)], old_params)


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


def run_find_synergy_test():
    pdf = JointProbabilityMatrix(2, 2)

    pdf_syn = pdf.copy()

    initial_guess_summed_modulo = np.random.choice([True, False])

    pdf_syn.append_synergistic_variables(1, initial_guess_summed_modulo=initial_guess_summed_modulo)

    assert pdf_syn.numvariables == pdf.numvariables + 1

    pdf_old = pdf_syn.marginalize_distribution(range(pdf.numvariables))

    assert pdf == pdf_old, 'adding and then removing variables should result in the same joint pdf'

    pdf_add_random = pdf.copy()
    pdf_add_random.append_variables(1)

    # note: this assert is in principle probabilistic, because who knows what optimization procedure is used and
    # how much it potentially sucks. So see if you hit this more than once, if you hit it at all.
    assert pdf_add_random.synergistic_information_by_difference(range(pdf.numvariables, pdf_add_random.numvariables),
                                                                range(pdf.numvariables)) <= \
           pdf_syn.synergistic_information_by_difference(range(pdf.numvariables, pdf_add_random.numvariables),
                                                         range(pdf.numvariables)), 'surely the optimization procedure' \
                                                                                   ' in append_synergistic_variables ' \
                                                                                   'yields a better syn. info. than ' \
                                                                                   'an appended variable with simply ' \
                                                                                   'random interaction parameters?!'

    parameters_before = pdf.matrix2params()

    assert np.testing.assert_array_almost_equal(pdf_add_random.matrix2params()[:len(parameters_before)],
                                                parameters_before)
    assert np.testing.assert_array_almost_equal(pdf_syn.matrix2params()[:len(parameters_before)], parameters_before)


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

    print 'note: finished. all tests successful'


def sum_modulo(values, modulo):
    """
    An example function which can be passed as the state_transitions argument to the
    append_variables_using_state_transitions_table function of JointProbabilityMatrix.
    :rtype : int
    :param values: list of values which should be in integer range [0, modulo)
    :param modulo: value self.numvalues usually, e.g. 2 for binary values
    :return: for binary variables it is the XOR, for others it is summed modulo of integers.
    """
    return int(np.mod(np.sum(values), modulo))