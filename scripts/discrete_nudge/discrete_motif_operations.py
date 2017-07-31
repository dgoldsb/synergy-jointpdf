"""
This file contains some of the common operations used on discrete motifs.
"""

from random import shuffle

from jointpdf_DJ_code.nudge import nudge_distribution_non_local_non_causal as nudge

def nudge_variable(motif, no_impacted, size=0.1):
    """
    Nudge a variable or a number of variables, using DJ his method.

    PARAMETERS
    ---
    motif: a DiscreteGrnMotif object
    no_impacted: number of variables impacted
    size: size of the nudge
    """
    joint = motif.joint_probabilities.joint_probabilities
    labels = []
    possible_labels = range(0, motif.numvariables)
    for _ in range(0, no_impacted):
        possible_labels = shuffle(possible_labels)
        labels.append(possible_labels.pop())
    new_joint = nudge(joint, labels, size)
    motif.joint_probabilities.joint_probabilities = new_joint
