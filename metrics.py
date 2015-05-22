__author__ = 'rquax'

'''
In this file I will create functions that implement 'metrics' which operate on a pdf object and return
a number of (nested) list of numbers, depending on the metric.
'''

from jointpdf import JointProbabilityMatrix


# standard labels for standard concepts
class standard_labels:
    # for a time variable (guess should usually only be one of these)
    time = 'time'

    # can be a spin's state, or a voter's opinion in a voter model, or a bird's position and direction in a
    # bird flocking model, or the concentration of a chemical, etc.
    variable = 'variable'

    # (usually in combination with a 'time' variable.) The time-dependent PDFs of the variables (with label .variable)
    # may depend on a further initial condition, such as a snapshot system state that is used to generate system
    # trajectories from, or otherwise the initial concentrations of chemicals for instance. An initial condition
    # is distinct from a 'parameter' in the sense that it is assumed to influence the time-dependent PDFs of variables
    # only at time=0, and subsequently it only exterts an indirect influence by dissipation of information or e.g.
    # in case if hysteresis.
    initial_condition = 'initial_condition'


'''
Metrics.
'''


def shannon_entropy_simple(pdf, dict_of_args={}):  # instructive example, intended to be as small and simple as possible
    pdf = pdf.marginalize_distribution_retaining_only_labels([standard_labels.time])

    return pdf.entropy()


def shannon_equilibrium_entropy_of_system_state(pdf, dict_of_args={}):
    """
    More advanced example.
    :type pdf: JointProbabilityMatrix
    :type dict_of_args: dict
    :rtype: float
    """

    unsupported_labels = set(pdf.get_labels()).difference({standard_labels.variable, standard_labels.time})

    if len(unsupported_labels) > 0:
        print 'debug: unsupported_labels =', unsupported_labels

    for uns_lbl in unsupported_labels:
        # these are variable indices with labels that I don't understand
        variables_with_label = pdf.get_variables_with_label_in([uns_lbl])

        # I assume that these variables with unsupported labels are parameters, like 'temperature', for which
        # the user provides me with a specific value for which I need to compute my metric
        # todo: recurse over all values for all unsupported labels and return a dict of metric values?
        if not dict_of_args.has_key(uns_lbl):
            raise ValueError('PDF contains variable(s) with label=' + str(uns_lbl) + ' which makes no sense for'
                             + ' computing Shannon entropy. If it is a parameter like \'temperature\' then you '
                               'should provide a specific value for me to compute entropy for')

        # note: I suppose it can also be possible for the user of CE to first have the option to select specific
        # values for parameters (such as a temperature value) before he/she decides to pass the joint pdf
        # to a metric... but I can also imagine that a well-implemented metric tries its best to handle
        # unexpected additional dimensions, which I try here. A default behavior of iterating over all possible
        # values of an unexpected (and unspecified in dict_of_args) dimension could make a lot of sense, then it
        # should not be possible for any metric to crash in case extra dimensions are passed, but perhaps it could
        # be that the memory runs out if a dimension has a gazillion possible values... Perhaps this standard behavior
        # can be standardized. Like a function call
        # "I_expect_only_dimensions(['variable', 'time'], metric_i, dict_of_args)" which does what I do above:
        # try to find specific values for unexpected dimensions, and if not specified then iterate over its
        # possible values. Or so.

        # keep only a subpart of the pdf, namely where the unsupported label has the provided value
        pdf = pdf.conditional_probability_distribution(variables_with_label, dict_of_args[uns_lbl])

        assert not pdf is None

    assert len(set(pdf.get_labels())) in (1, 2), 'I expect either "variable" or "variable" and "time" as labels, ' \
                                                 'no more or less'

    # sum out over time, only keep pdf of variables (states)
    pdf = pdf.marginalize_distribution_retaining_only_labels([standard_labels.variable])

    return pdf.entropy()


def idt(pdf, dict_of_args={}):
    """

    :type pdf: JointProbabilityMatrix
    :param dict_of_args: dict of arbitrary arguments you want to pass, depending on the metric.
    Like a specific value for unsupported dimensions such as 'temperature', if any.
    """

    # todo: make a separate class ConditionalProbabilityistribution, wrapping the dict?

    pdf.retain_only_labels([standard_labels.time, standard_labels.variable, standard_labels.initial_condition])

    # these three are needed for idt
    assert standard_labels.time in pdf.get_labels()
    assert standard_labels.variable in pdf.get_labels()
    assert standard_labels.initial_condition in pdf.get_labels()

    # todo: implement a mutual_information which takes labels as arguments


'''
Unit testing.
'''


def run_equilibrium_entropy_test():
    data = [[1,2],[0,3],[0,1],[0,1]]

    pdf = JointProbabilityMatrix(2, 2)  # todo: enable this in one call
    pdf.estimate_from_data(data)

    assert pdf.entropy() == 1.5

    pdf.append_variables(1)
    pdf.set_label(2, 'time')  # this is a parameter that will be summed out to get the equilibrium entropy

    assert shannon_equilibrium_entropy_of_system_state(pdf) == 1.5
    assert shannon_entropy_simple(pdf) == 1.5


def run_all_tests():
    run_equilibrium_entropy_test()
    print 'note: run_equilibrium_entropy_test successful.'