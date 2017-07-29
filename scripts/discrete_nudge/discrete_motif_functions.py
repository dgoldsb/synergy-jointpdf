"""
This file contains some of the common operations used on discrete motifs.
"""

def dictionary():
    """
    Give a dictionary with all available functions,
    plus the number of outputs.

    RETURNS
    ---
    function_dict: dictionary of all functions with I/O requirements
    """
    function_dict = []

    # plus
    plus_d = {}
    plus_d["f"] = plus
    plus_d["i"] = 1
    plus_d["o"] = 1
    function_dict.append(plus_d)

    # min
    minus_d = {}
    minus_d["f"] = minus
    minus_d["i"] = 1
    minus_d["o"] = 1
    function_dict.append(minus_d)

    # plus_and
    plus_and_d = {}
    plus_and_d["f"] = plus_and
    plus_and_d["i"] = 2
    plus_and_d["o"] = 1
    function_dict.append(plus_and_d)

    # min_and
    min_and_d = {}
    min_and_d["f"] = min_and
    min_and_d["i"] = 2
    min_and_d["o"] = 1
    function_dict.append(min_and_d)

    # xor
    xor_d = {}
    xor_d["f"] = xor
    xor_d["i"] = 2
    xor_d["o"] = 1
    function_dict.append(xor_d)

    return function_dict

# 1-to-1 functions
def plus(inputs, output):
    """
    Simple stimulation.
    """
    if output == 1:
        return output
    else:
        if inputs == [1]:
            return 1
        else:
            return output

def minus(inputs, output):
    """
    Simple surpression.
    """
    if output == 0:
        return output
    else:
        if inputs == [1]:
            return 0
        else:
            return output

# 2-to-1 functions
def plus_and(inputs, output):
    """
    Simple stimulation iff all inputs are stimulated.
    """
    if [0] in inputs:
        return output
    else:
        return 1

def min_and(inputs, output):
    """
    Simple stimulation iff all inputs are stimulated.
    """
    if [0] in inputs:
        return output
    else:
        return 0

def xor(inputs, output):
    """
    Simple XOR function
    """
    if inputs[0] == 1 and inputs[1] == 0:
        return 1
    elif inputs[1] == 1 and inputs[0] == 0:
        return 1
    else:
        return 0
 