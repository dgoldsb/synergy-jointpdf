# Version to run without jupyter, improved by using Pandas for an easy to work with dataset

# imports from the discrete motif package
from discrete_motif import DiscreteGrnMotif
import discrete_motif_functions as functions
import discrete_motif_measures as measures
import discrete_motif_generator as generator
import discrete_motif_operations as operations
import discrete_motif_plotting as visualize

# regular imports
import logging
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle
import random
import scipy.stats
from sklearn.manifold import TSNE
import sys
import tabulate
import time
from time import gmtime, strftime

# settings for the experiment
synergy_measure = measures.synergy_middleground
nudge_method = "DJ"
sample_size = 200
network_sizes = [2, 3, 4, 5]
logic_sizes = [2, 3, 4]
nudge_sizes = [0.1, 0.25, 0.5]


def visualize_TSNE(dataframe):
    """
    We use t-SNE to visualize our sample. What we see makes perfect sense: the biological motifs seem to be a subset of the sample space. In higher valued logic systems, the space becomes larger, making the sample size insufficient. As t-SNE also does clustering, it then starts to appear that we have two seperate sample, but this is simply the separation from the very similar biological motifs from the rest.
    """

    # x_train should be an array of N vectors of length M
    # where N is the sample size and M is the length of the combined new state of the transition table
    # we combine both the random and bio samples, but give them different colors

    for i in range(0, len(samples)):
        x_train = []
        y_color = []
        labels = []
        # the random tables
        for sample in samples[i]["random_tables"]:
            table = sample[0]
            sample_vector = []
            for row in table.transition_table:
                sample_vector.extend(row[(len(row)/2):])
            y_color.append("red")
            x_train.append(sample_vector)
            labels.append("Random transition table")

        # the bio table
        for sample in samples[i]["bio_tables"]:
            motif = sample[0]
            sample_vector = []
            for row in motif.transition_table:
                sample_vector.extend(row[(len(row)/2):])
            y_color.append("blue")
            x_train.append(sample_vector)
            labels.append("Biological transition table")
    
        x_train = np.array(x_train)
        title = "Two-dimensional embedding of transition tables with motif size %s and %s-valued logic" % (samples[i]["network_size"], samples[i]["logic_size"])
        x_train_embedded = TSNE(n_components=2, perplexity=10, early_exaggeration=1, verbose=0).fit_transform(x_train)
        
        results_location = "../../data/n=%d_l=%d_e=%f_%s_%s_nosamples=%d" % (samples[i]["network_size"], samples[i]["logic_size"], samples[i]["nudge_size"], nudge_method, synergy_measure, len(samples[i]))
        if not os.path.isdir(results_location):
            os.mkdir(results_location)
        filename = "tsne2D.pdf" 
        visualize.plot_scatter(x_train_embedded, y_color, labels, title, os.path.join(results_location, filename))


def visualize_scatters(dataframe):
    """
    """

    # Synergy and Nudge Impact

    for i in range(0, len(samples)):
        for j in range(0, samples[i]["network_size"]):
            x_values = []
            colors = []
            labels = []
            title = "Synergy vs. Nudge impact with motif size %s and %s-valued logic" % \
                    (samples[i]["network_size"], samples[i]["logic_size"])
            axes_labels = ["Nudge impact", "Synergy"]
            # the random tables
            for sample in samples[i]["random_tables"]:
                if len(sample[2]) > 0:
                    if sample[1] is not None and sample[2][j][1] is not None:
                        x_values.append([sample[2][j][1], sample[1]])
                        colors.append("red")
                        labels.append("Random transition table")

            # the bio table
            for sample in samples[i]["bio_tables"]:
                if len(sample[2]) > 0:
                    if sample[1] is not None and sample[2][j][1] is not None:
                        x_values.append([sample[2][j][1], sample[1]])
                        colors.append("blue")
                        labels.append("Biological transition table")

            x_values = np.array(x_values)

            results_location = "../../data/n=%d_l=%d_e=%f_%s_%s_nosamples=%d" % (samples[i]["network_size"], samples[i]["logic_size"], samples[i]["nudge_size"], nudge_method, synergy_measure, len(samples[i]))
            filename = "scatter2D_synergy_resilience.pdf"
            visualize.plot_scatter(x_values, colors, labels, title, os.path.join(results_location, filename), axes_labels)

    # Synergy and Memory

    for i in range(0, len(samples)):
        x_values = []
        colors = []
        labels = []
        title = "Synergy vs. Memory with motif size %s and %s-valued logic, "\
                "%s genes targeted, %s nudge size" % \
                (samples[i]["network_size"], samples[i]["logic_size"], j+1, samples[i]["nudge_size"])
        axes_labels = ["Memory", "Synergy"]
        # the random tables
        for sample in samples[i]["random_tables"]:
            if len(sample[2]) > 0:
                if sample[1] is not None and sample[3] is not None:
                    x_values.append([sample[3], sample[1]])
                    colors.append("red")
                    labels.append("Random transition table")

        # the bio table
        for sample in samples[i]["bio_tables"]:
            if len(sample[2]) > 0:
                if sample[1] is not None and sample[3] is not None:
                    x_values.append([sample[3], sample[1]])
                    colors.append("blue")
                    labels.append("Biological transition table")

        x_values = np.array(x_values)

        results_location = "../../data/n=%d_l=%d_e=%f_%s_%s_nosamples=%d" % (samples[i]["network_size"], samples[i]["logic_size"], samples[i]["nudge_size"], nudge_method, synergy_measure, len(samples[i]))
        filename = "scatter2D_synergy_memory.pdf"
        visualize.plot_scatter(x_values, colors, labels, title, os.path.join(results_location, filename), axes_labels)

    # Memory and Nudge Impact

    for i in range(0, len(samples)):
        for j in range(0, samples[i]["network_size"]):
            x_values = []
            colors = []
            labels = []
            title = "Memory vs. Nudge impact with motif size %s and %s-valued logic, "\
                    "%s genes targeted, %s nudge size" % \
                    (samples[i]["network_size"], samples[i]["logic_size"], j+1, samples[i]["nudge_size"])
            axes_labels = ["Nudge impact", "Memory"]
            # the random tables
            for sample in samples[i]["random_tables"]:
                if len(sample[2]) > 0:
                    if sample[3] is not None and sample[2][j][1] is not None:
                        x_values.append([sample[2][j][1], sample[3]])
                        colors.append("red")
                        labels.append("Random transition table")

            # the bio table
            for sample in samples[i]["bio_tables"]:
                if len(sample[2]) > 0:
                    if sample[3] is not None and sample[2][j][1] is not None:
                        x_values.append([sample[2][j][1], sample[3]])
                        colors.append("blue")
                        labels.append("Biological transition table")
                
            x_values = np.array(x_values)
            
            results_location = "../../data/n=%d_l=%d_e=%f_%s_%s_nosamples=%d" % (samples[i]["network_size"], samples[i]["logic_size"], samples[i]["nudge_size"], nudge_method, synergy_measure, len(samples[i]))
            filename = "scatter2D_memory_resilience.pdf"
            visualize.plot_scatter(x_values, colors, labels, title, os.path.join(results_location, filename), axes_labels)

    # 3D scatter

    for i in range(0, len(samples)):       
        for j in range(0, samples[i]["network_size"]):
            x_values = []
            colors = []
            labels = []
            title = "Synergy vs. Nudge impact with motif size %s and %s-valued logic, "\
                    "%s genes targeted, %s nudge size" % \
                    (samples[i]["network_size"], samples[i]["logic_size"], j+1, samples[i]["nudge_size"])
            axes_labels = ["Nudge impact", "Synergy", "Memory"]
            # the random tables
            for sample in samples[i]["random_tables"]:
                if len(sample[2]) > 0:
                    if sample[1] is not None and sample[2][j][1] is not None:
                        x_values.append([sample[2][j][1], sample[1], sample[3]])
                        colors.append("red")
                        labels.append("Random transition table")

            # the bio table
            for sample in samples[i]["bio_tables"]:
                if len(sample[2]) > 0:
                    if sample[1] is not None and sample[2][j][1] is not None:
                        x_values.append([sample[2][j][1], sample[1], sample[3]])
                        colors.append("blue")
                        labels.append("Biological transition table")
                
            x_values = np.array(x_values)
            
            results_location = "../../data/n=%d_l=%d_e=%f_%s_%s_nosamples=%d" % (samples[i]["network_size"], samples[i]["logic_size"], samples[i]["nudge_size"], nudge_method, synergy_measure, len(samples[i]))
            filename = "scatter3D_memory_synergy_resilience.pdf"
            
            visualize.plot_scatter_3d(x_values, colors, labels, title, os.path.join(results_location, filename), axes_labels)



def visualize_profile(dataframe):
    """
    """

    
def visualize_memory(dataframe):
    """
    """

    
def test_synergy(dataframe):
    """
    """

    
def test_memory(dataframe):
    """
    """

    
def test_resilience(dataframe):
    """
    """

    
def loop_impacts(network_size, nudge_size, motif):
    # find the memory
    memory = measures.normalized_memory(motif)
    
    # try to find the synergy
    synergy = measures.normalized_synergy(motif, synergy_measure)
    
    # try to loop
    impacts = []
    for nudge_width in range(1, network_size + 1):
        if not (nudge_method == 'joint_pdf' and nudge_width == network_size):
            # we compare the two evolved states
            impact = measures.average_nudge_impact(motif, nudge_width, nudge_size, nudge_method)
            impact_tuple = (nudge_width, impact)
            impacts.append(impact_tuple)
    
    return synergy, impacts, memory

def draw_sample(sample_size, network_size, logic_size, nudge_size):
    # create dataframe
    dataframe = pd.DataFrame(columns=["system_size", "logic_size", "nudge_size", "type", "motif", "synergy", "memory", "impacts"])

    # we generate our samples
    samples_random = generator.generate_random(sample_size, network_size, logic_size)
    samples_grn = generator.generate_motifs(sample_size, network_size, logic_size)[0]

    # draw the samples
    print("Sampling with %s nodes and %s-valued logic" % (network_size, logic_size))
    loc_counter = 0
    for motif in samples_grn:
        # get the basic inputs in the dataframe row
        df_row = []
        df_row.append(network_size)
        df_row.append(logic_size)
        df_row.append(nudge_size)
        df_row.append("GRN")
        df_row.append(motif)
        motif.set_transition_table()

        # compute outputs of experiments
        synergy, impacts, memory = loop_impacts(network_size, nudge_size, motif)

        # add to dataframe
        df_row.append(synergy)
        df_row.append(memory)
        df_row.append(impacts)

        dataframe.loc[loc_counter] = df_row
        loc_counter += 1

    # enrich the random tables
    for motif in samples_random:
        # get the basic inputs in the dataframe row
        df_row = []
        df_row.append(network_size)
        df_row.append(logic_size)
        df_row.append(nudge_size)
        df_row.append("random")
        df_row.append(motif)
        motif.set_transition_table()

        # compute outputs of experiments
        synergy, impacts, memory = loop_impacts(network_size, nudge_size, motif)

        # add to dataframe
        df_row.append(synergy)
        df_row.append(memory)
        df_row.append(impacts)

        dataframe.loc[loc_counter] = df_row
        loc_counter += 1
    
    return dataframe, samples_grn, samples_random

def main():
    # set folders
    data_location = "../../data"
    log_location = "../../log"

    # logger
    mylogger = logging.getLogger('mylogger')
    handler1 = logging.FileHandler(filename=os.path.join(log_location, 'experiments_%s.log' % strftime("%Y-%m-%d %H:%M:%S", gmtime())), mode='w')
    handler1.setLevel(logging.DEBUG)
    handler1.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
    mylogger.addHandler(handler1)

    directory = "../../data"
    if not os.path.isdir(directory):
        os.mkdir(directory)

    # draw a few completely random samples, with different parameters
    for network_size in network_sizes:
        for logic_size in logic_sizes:
            for nudge_size in nudge_sizes:
                mylogger.debug("sampling %d nodes, %d logic size, %f nudge size, %s as nudge_method, %s as synergy measure" % (network_size, logic_size, nudge_size, nudge_method, synergy_measure))
                start = time.time()
                dataframe, samples_grn, samples_random = draw_sample(sample_size, network_size, logic_size, nudge_size)
                
                # save the data for future use/reruns
                name_df = "experiment_k=%d_l=%d_e=%f_%s_%s_n=%d_df.pkl" % (network_size, logic_size, nudge_size, nudge_method, synergy_measure, sample_size)
                name_grn = "samples_grn_k=%d_l=%d_n=%d.pkl" % (network_size, logic_size, sample_size)
                name_random = "samples_random_k=%d_l=%d_n=%d.pkl" % (network_size, logic_size, sample_size)
                with open(os.path.join(data_location, name_df), 'wb') as output:
                    pickle.dump(dataframe, output, pickle.HIGHEST_PROTOCOL)
                with open(os.path.join(data_location, name_grn), 'wb') as output:
                    pickle.dump(samples_grn, output, pickle.HIGHEST_PROTOCOL)
                with open(os.path.join(data_location, name_random), 'wb') as output:
                    pickle.dump(samples_random, output, pickle.HIGHEST_PROTOCOL)

                # log the finished experiment
                end = time.time()
                mylogger.debug("sampled %d motifs" % sample_size)
                mylogger.debug("sample took %d seconds" % (end - start))

    # load data back, append the dataframes
    dataframe = pd.DataFrame(columns=["system_size", "logic_size", "nudge_size", "type", "motif", "synergy", "memory", "impacts"])
    for root, dirnames, filenames in os.walk(data_location):
        for filename in filenames:
            if filename.endswith('df.pkl'):
                with open(os.path.join(root, filename), 'rb') as input:
                    dataframe.append(pickle.load(input))

    # now let's get rolling and do our tests!
    visualize_TSNE(dataframe)
    visualize_scatters(dataframe)
    visualize_profile(dataframe)
    visualize_memory(dataframe)
    test_synergy(dataframe)
    test_memory(dataframe)
    test_resilience(dataframe)
