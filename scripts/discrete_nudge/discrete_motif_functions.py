"""
This file contains some of the common operations used on discrete motifs.
"""

# 1-to-1 functions
def plus(input, output):
    if output == 1:
        return 1
    else:
        if input == 1:
            return 1
        else:
            return 0

def minus(input, output):
    if output == 0:
        return 0
    else:
        if input == 1:
            return 0
        else:
            return 1

# many-to-1 functions
def plus_AND(inputs, output):
    """
    AND-function for gene regulation.
    """
    if 0 in inputs:
        return 0
    else:
        return 1
 