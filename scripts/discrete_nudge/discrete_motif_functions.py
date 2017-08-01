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
    function_dict = [[], [], []]

    # plus
    plus_d = {}
    plus_d["f"] = plus
    plus_d["i"] = 1
    plus_d["o"] = 1
    function_dict[1].append(plus_d)

    # min
    minus_d = {}
    minus_d["f"] = minus
    minus_d["i"] = 1
    minus_d["o"] = 1
    function_dict[1].append(minus_d)

    # and
    plus_and_d = {}
    plus_and_d["f"] = plus_and
    plus_and_d["i"] = 2
    plus_and_d["o"] = 1
    function_dict[2].append(plus_and_d)

    # Nand
    min_and_d = {}
    min_and_d["f"] = min_and
    min_and_d["i"] = 2
    min_and_d["o"] = 1
    function_dict[2].append(min_and_d)

    # xor
    xor_d = {}
    xor_d["f"] = xor
    xor_d["i"] = 2
    xor_d["o"] = 1
    function_dict[2].append(xor_d)

    # Nxor
    nxor_d = {}
    nxor_d["f"] = nxor
    nxor_d["i"] = 2
    nxor_d["o"] = 1
    function_dict[2].append(nxor_d)

    # or
    or_d = {}
    or_d["f"] = por
    or_d["i"] = 2
    or_d["o"] = 1
    function_dict[2].append(or_d)

    # Nor
    nor_d = {}
    nor_d["f"] = nor
    nor_d["i"] = 2
    nor_d["o"] = 1
    function_dict[2].append(nor_d)

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
 
def nxor(inputs, output):
    """
    Simple NXOR function
    """
    if inputs[0] == 1 and inputs[1] == 0:
        return 0
    elif inputs[1] == 1 and inputs[0] == 0:
        return 0
    else:
        return 1

def por(inputs, output):
    """
    Simple OR function
    """
    if inputs[0] == 0 and inputs[1] == 0:
        return 0
    else:
        return 1

def nor(inputs, output):
    """
    Simple NOR function
    """
    if inputs[0] == 0 and inputs[1] == 0:
        return 1
    else:
        return 0
