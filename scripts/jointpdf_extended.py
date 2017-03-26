'''
Author: Dylan Goldsborough
Email:  dgoldsb@live.nl

This script extends the jointpdf package by Rick Quax.
'''

from __future__ import print_function
from jointpdf import jointpdf

__author__ = 'dgoldsb'

class JointProbabilityMatrix(jointpdf.JointProbabilityMatrix):
    """
    This class extends the JointProbilityMatrix written by Rick Quax.
    """

    # TODO: add any methods for discrete PDFs, especially things that have a continuous parallel

class JointProbabilityKNN():
    """
    Continuous extension for the code written by Rick Quax.
    Should use the same methods, to have full compatibility with MIprofile.py
    """