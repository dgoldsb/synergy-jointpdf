'''
This file contains some of the common operations used on discrete motifs.
'''

import math

def dictionary():
    '''
    Give a dictionary with all available functions,
    plus the number of outputs.

    @returns function_dict: dictionary of all functions with I/O requirements
    '''
    function_dict = [[], [], []]

    # plus
    plus_d = {}
    plus_d['f'] = plus
    plus_d['i'] = 1
    plus_d['o'] = 1
    function_dict[1].append(plus_d)

    # min
    minus_d = {}
    minus_d['f'] = minus
    minus_d['i'] = 1
    minus_d['o'] = 1
    function_dict[1].append(minus_d)

    '''
    # superplus
    superplus_d = {}
    superplus_d['f'] = superplus
    superplus_d['i'] = 1
    superplus_d['o'] = 1
    function_dict[1].append(superplus_d)

    # supermin
    superminus_d = {}
    superminus_d['f'] = superminus
    superminus_d['i'] = 1
    superminus_d['o'] = 1
    function_dict[1].append(superminus_d)
    '''

    # and
    plus_and_d = {}
    plus_and_d['f'] = plus_and
    plus_and_d['i'] = 2
    plus_and_d['o'] = 1
    function_dict[2].append(plus_and_d)

    # Nand
    min_and_d = {}
    min_and_d['f'] = min_and
    min_and_d['i'] = 2
    min_and_d['o'] = 1
    function_dict[2].append(min_and_d)

    '''
    # xor
    xor_d = {}
    xor_d['f'] = xor
    xor_d['i'] = 2
    xor_d['o'] = 1
    function_dict[2].append(xor_d)

    # negative xor
    neg_xor_d = {}
    neg_xor_d['f'] = neg_xor
    neg_xor_d['i'] = 2
    neg_xor_d['o'] = 1
    function_dict[2].append(neg_xor_d)

    # incl or
    incl_or_d = {}
    incl_or_d['f'] = incl_or
    incl_or_d['i'] = 2
    incl_or_d['o'] = 1
    function_dict[2].append(incl_or_d)

    # negative incl or
    neg_incl_or_d = {}
    neg_incl_or_d['f'] = neg_incl_or
    neg_incl_or_d['i'] = 2
    neg_incl_or_d['o'] = 1
    function_dict[2].append(neg_incl_or_d)

    # excl or
    excl_or_d = {}
    excl_or_d['f'] = excl_or
    excl_or_d['i'] = 2
    excl_or_d['o'] = 1
    function_dict[2].append(excl_or_d)

    # negative excl or
    neg_excl_or_d = {}
    neg_excl_or_d['f'] = neg_excl_or
    neg_excl_or_d['i'] = 2
    neg_excl_or_d['o'] = 1
    function_dict[2].append(neg_excl_or_d)
    '''

    return function_dict


# 1-to-1 functions
def superplus(input):
    '''
    Strong stimulation.

    @param input: list of inputs
    @returns: the increase/decrease of the target gene (integer)
    '''
    increase = input[0] * 2
    return increase


def plus(input):
    '''
    Simple stimulation.

    @param input: list of inputs
    @returns: the increase/decrease of the target gene (integer)
    '''
    increase = input[0] * 1
    return increase


def minus(input):
    '''
    Simple surpression.

    @param input: list of inputs
    @returns: the increase/decrease of the target gene (integer)
    '''
    increase = input[0] * -1
    return increase


def superminus(input):
    '''
    Strong surpression.

    @param input: list of inputs
    @returns: the increase/decrease of the target gene (integer)
    '''
    increase = input[0] * -2
    return increase


# 2-to-1 functions
def plus_and(inputs):
    '''
    Simple stimulation iff all inputs are stimulated.

    @param input: list of inputs
    @returns: the increase/decrease of the target gene (integer)
    '''
    increase = min(inputs) * 1
    return increase


def min_and(inputs):
    '''
    Simple stimulation iff all inputs are stimulated.

    @param input: list of inputs
    @returns: the increase/decrease of the target gene (integer)
    '''
    increase = min(inputs) * -1
    return increase


def xor(inputs):
    '''
    Simple XOR function.

    @param input: list of inputs
    @returns: the increase/decrease of the target gene (integer)
    '''
    increase = math.fabs(inputs[0] - inputs[1]) * 1
    return increase


def neg_xor(inputs):
    '''
    Simple NXOR function.

    @param input: list of inputs
    @returns: the increase/decrease of the target gene (integer)
    '''
    return -1 * xor(inputs)


def incl_or(inputs):
    '''
    Simple OR function.

    @param input: list of inputs
    @returns: the increase/decrease of the target gene (integer)
    '''
    increase = min(inputs) * 1
    return increase


def neg_incl_or(inputs):
    '''
    Simple negative OR function.

    @param input: list of inputs
    @returns: the increase/decrease of the target gene (integer)
    '''
    return incl_or(inputs) * -1


def excl_or(inputs):
    '''
    Simple exclusive OR function.

    @param input: list of inputs
    @returns: the increase/decrease of the target gene (integer)
    '''
    increase = 0
    if min(inputs) == 0:
        increase = max(inputs) + 1
    return increase


def neg_excl_or(inputs):
    '''
    Simple negative exclusive OR function.

    @param input: list of inputs
    @returns: the increase/decrease of the target gene (integer)
    '''
    return excl_or(inputs) * -1
