__author__ = 'rquax'

'''
In this file I will create a list of functions that implement 'metrics' which operate on a pdf object and return
a number of (nested) list of numbers, depending on the metric.
'''

from jointpdf import JointProbabilityMatrix


def shannon_entropy_simple(pdf, dict_of_args={}):  # instructive example
    pdf = pdf.marginalize_distribution_retaining_only_labels(['variable'])

    return pdf.entropy()


def shannon_equilibrium_entropy_of_system_state(pdf, dict_of_args={}):
    """
    More advanced example.
    :type pdf: JointProbabilityMatrix
    :type dict_of_args: dict
    :rtype: float
    """

    unsupported_labels = set(pdf.get_labels()).difference({'variable', 'time'})

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

        # keep only a subpart of the pdf, namely where the unsupported label has the provided value
        pdf = pdf.conditional_probability_distribution(variables_with_label, dict_of_args[uns_lbl])

        assert not pdf is None

    assert len(pdf.get_labels()) in (1, 2), 'I expect either variable or variable and time as labels, no more or less'

    # sum out over time, only keep pdf of variables (states)
    pdf = pdf.marginalize_distribution_retaining_only_labels(['variable'])

    return pdf.entropy()