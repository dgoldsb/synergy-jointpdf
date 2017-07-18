"""
This file contains some of the common operations used on discrete motifs.
"""

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
            return 0

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
            return 1

# many-to-1 functions
def plus_and(inputs, output):
    """
    Simple stimulation iff all inputs are stimulated.
    """
    if [0] in inputs:
        return output
    else:
        return 1
 