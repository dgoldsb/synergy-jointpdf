# imports from the discrete motif package
from discrete_motif import DiscreteGrnMotif
import discrete_motif_functions as functions
import discrete_motif_measures as measures
import discrete_motif_generator as generator
import discrete_motif_operations as operations
import discrete_motif_plotting as visualize

# regular imports
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import random
import scipy.stats
import re
from sklearn.manifold import TSNE
import sys

# set folders
data_location = "../../data_pandas"
result_location = "../../result_pandas"
log_location = "../../log"


def main():
    samples = generator.generate_motifs(1, 3, 2, [4])[0]
    print(samples[0].grn_vars["rules"])
    print(len(samples[0].grn_vars["rules"]))

    # problem demonstrated below, the number of rules is waaaay too low
    sys.exit(0)
    with open(os.path.join(data_location, "experiment_k=5_l=2_e=0.100000_df.pkl"), 'rb') as input:
        df = pickle.load(input)
    print(df)

    # the following seems fine...
    with open(os.path.join(data_location, "samples_grn_k=4_l=2.pkl"), 'rb') as input:
        samples = pickle.load(input)

    synergies = []
    memories = []
    for sample in samples:
        if len(sample.grn_vars["rules"]) > 1:
            print(len(sample.grn_vars["rules"]))
    """
            print(sample.transition_table)
            print(sample.grn_vars["rules"])
            memories.append(measures.normalized_memory(sample))
            synergies.append(measures.normalized_synergy(sample, measures.synergy_middleground))
            print(memories)
            print(synergies)
    print("%.2f %.2f" % (max(synergies), min(synergies)))
    print("%.2f %.2f" % (max(memories), min(memories)))
    """

if __name__ == '__main__':
    main()