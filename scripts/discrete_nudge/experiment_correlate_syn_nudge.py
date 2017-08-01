"""
The aim is to correlate the impact of a nudge (on a single variable or on multiple) with the WMS synergy and the synergy as defined by Rick Quax.
We work in a discrete setting, and work with the jointpdf framework, as well as with extensions made by Derk-Jan Riesthuis.
"""

from __future__ import print_function

import sys

from discrete_motif_generator import generate_motifs
from discrete_motif_plotting import scatterplot_synergy_nudgeimpact


def main():
    """
    This code will be mostly the main.
    """
    # create a network motif
    motifs, indegree_avg = generate_motifs(50, 3)
    print("average indegree is "+str(indegree_avg))

    # make a plot
    scatterplot_synergy_nudgeimpact(motifs, 3, 0.7, False)

if __name__ == '__main__':
    main()
