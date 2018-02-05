#%latex imports for tables
#\usepackage{array}
#\usepackage{diagbox}
#\usepackage{multirow}
#\usepackage[table]{xcolor}

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
from sklearn.linear_model import LinearRegression
from sklearn.manifold import TSNE
import sys
from tqdm import tqdm

# set folders
data_location = "../../data_pandas"
result_location = "../../result_pandas"
log_location = "../../log"


#TODO maak virtualenv
def outliers_iqr(ys):
    """
    Credit for this function to http://colingorrie.github.io/outlier-detection.html
    """
    quartile_1, quartile_3 = np.percentile(ys, [25, 75])
    iqr = quartile_3 - quartile_1
    lower_bound = quartile_1 - (iqr * 1.5)
    upper_bound = quartile_3 + (iqr * 1.5)
    return np.where((ys > upper_bound) | (ys < lower_bound))


def histogram_cycles(dataframe):
    """
    We use t-SNE to visualize our sample. What we see makes perfect sense: the biological motifs seem to be a subset of the sample space. In higher valued logic systems, the space becomes larger, making the sample size insufficient. As t-SNE also does clustering, it then starts to appear that we have two seperate sample, but this is simply the separation from the very similar biological motifs from the rest.
    """
    # always start by finding all experiments
    experiments_dataframe = dataframe[["system_size", "logic_size"]]
    experiments = experiments_dataframe.drop_duplicates().values.tolist()
    # loop over the experiments
    for experiment in experiments:
        # fill the variables
        system_size = experiment[0]
        logic_size = experiment[1]

        # get the wanted data
        a = dataframe["logic_size"] == logic_size
        b = dataframe["system_size"] == system_size
        selected_dataframe = dataframe[a & b]
        samples = selected_dataframe[["type", "motif"]].values.tolist()

        random_cycles = []
        bio_cycles = []
        for sample in tqdm(samples):
            lengths = [len(x) for x in sample[1].is_cyclical(int(logic_size ** system_size))]
            if len(lengths) == 0:
                lengths = [0]
            if int(np.sum(lengths)) > (logic_size ** system_size):
                print(sample[1].is_cyclical(int(logic_size ** system_size)))
                print(sample[1].transition_table)
                print("impossible")
                sys.exit(1)
            max_length = int(np.max(lengths))
            if sample[0] == "random":
                random_cycles.append(max_length)
            elif sample[0] == "GRN":
                bio_cycles.append(max_length)
        
        colors = ["r", "b"]
        labels = ["random", "biological random"]
        #title = "Distribution of cycle lengths (0 if no cycle is present)"
        title = None
        axes_labels = ["Maximum length of cycle in motif", "Number of occurences"]
        foldername = "k=%d_l=%d" % (system_size, logic_size)
        folderpath = os.path.join(result_location, foldername)
        if not os.path.isdir(folderpath):
            os.mkdir(folderpath)
        filename = "cycles.pdf" 
        visualize.plot_bar([random_cycles, bio_cycles], colors, labels, title, os.path.join(folderpath, filename), axes_labels)


def visualize_TSNE(dataframe):
    """
    We use t-SNE to visualize our sample. What we see makes perfect sense: the biological motifs seem to be a subset of the sample space. In higher valued logic systems, the space becomes larger, making the sample size insufficient. As t-SNE also does clustering, it then starts to appear that we have two seperate sample, but this is simply the separation from the very similar biological motifs from the rest.
    """
    # always start by finding all experiments
    experiments_dataframe = dataframe[["system_size", "logic_size"]]
    experiments = experiments_dataframe.drop_duplicates().values.tolist()
    # loop over the experiments
    for experiment in experiments:
        # fill the variables
        system_size = experiment[0]
        logic_size = experiment[1]

        # get the wanted data
        a = dataframe["logic_size"] == logic_size
        b = dataframe["system_size"] == system_size
        selected_dataframe = dataframe[a & b]
        samples = selected_dataframe[["type", "motif"]].values.tolist()

        # x_train should be an array of N vectors of length M
        # where N is the sample size and M is the length of the combined new state of the transition table
        # we combine both the random and bio samples, but give them different colors
        x_train = []
        y_color = []
        labels = []

        for sample in tqdm(samples):
            table = sample[1]
            sample_vector = []
            for row in table.transition_table:
                sample_vector.extend(row[(len(row)/2):])
            x_train.append(sample_vector)
            if sample[0] == 'random':
                y_color.append("red")
                labels.append("Random transition table")
            elif sample[0] == 'GRN':
                y_color.append("blue")
                labels.append("Biological transition table")

        # run TSNE
        x_train = np.array(x_train)
        x_train_embedded = TSNE(n_components=2, perplexity=10, early_exaggeration=1, verbose=0).fit_transform(x_train)
        
        # save the results
        title = None
        # title = "Two-dimensional embedding of transition tables with motif size %s and %s-valued logic" % (system_size, logic_size)
        foldername = "k=%d_l=%d" % (system_size, logic_size)
        folderpath = os.path.join(result_location, foldername)
        if not os.path.isdir(folderpath):
            os.mkdir(folderpath)
        filename = "tsne2D.pdf" 
        visualize.plot_scatter(x_train_embedded, y_color, labels, title, os.path.join(folderpath, filename))


def visualize_scatters(dataframe):
    """
    Make a series of scatterplots that visualize the spreads of our experiments
    """
    # first we do experiments that include nudge size as a factor
    # always start by finding all experiments
    experiments_dataframe = dataframe[["system_size", "logic_size", "nudge_size"]]
    experiments = experiments_dataframe.drop_duplicates().values.tolist()
    # loop over the experiments
    for experiment in experiments:
        # fill the variables
        system_size = experiment[0]
        logic_size = experiment[1]
        nudge_size = experiment[2]

        # get the wanted data
        a = dataframe["logic_size"] == logic_size
        b = dataframe["system_size"] == system_size
        c = dataframe["nudge_size"] == nudge_size
        selected_dataframe = dataframe[a & b & c]
        samples = selected_dataframe[["type", "synergy", "impacts", "memory"]].values.tolist()

        # PLOT: synergy - impact
        for i in range(0, int(system_size)):
            x_values = []
            colors = []
            labels = []

            # read the samples
            for sample in samples:
                if len(sample[2]) > 0:
                    if sample[1] is not None and sample[2][i][1] is not None:
                        x_values.append([sample[1], sample[2][i][1]])
                        if sample[0] == "random":
                            colors.append("red")
                            labels.append("Random transition table")
                        elif sample[0] == "GRN":
                            colors.append("blue")
                            labels.append("Biological transition table")
            x_values = np.array(x_values)

            # plot and save
            foldername = "k=%d_l=%d_e=%f" % (system_size, logic_size, nudge_size)
            folderpath = os.path.join(result_location, foldername)
            if not os.path.isdir(folderpath):
                os.mkdir(folderpath)
            title = None
            axes_labels = ["Synergy", "Nudge impact"]
            filename = "scatter2D_synergy_resilience_width=%s.pdf" % str(i+1)
            visualize.plot_scatter(x_values, colors, labels, title, os.path.join(folderpath, filename), axes_labels)

        # PLOT: memory - impact
        for i in range(0, int(system_size)):
            x_values = []
            colors = []
            labels = []

            # read the samples
            for sample in samples:
                if len(sample[2]) > 0:
                    if sample[3] is not None and sample[2][i][1] is not None:
                        x_values.append([sample[3], sample[2][i][1]])
                        if sample[0] == "random":
                            colors.append("red")
                            labels.append("Random transition table")
                        elif sample[0] == "GRN":
                            colors.append("blue")
                            labels.append("Biological transition table")
            x_values = np.array(x_values)

            # plot and save
            title = None
            axes_labels = ["Memory", "Nudge impact"]
            filename = "scatter2D_memory_resilience_width=%s.pdf" % str(i+1)
            visualize.plot_scatter(x_values, colors, labels, title, os.path.join(folderpath, filename), axes_labels)

        # PLOT: 3D scatter
        for j in range(0, int(system_size)):
            x_values = []
            colors = []
            labels = []

            # read the samples
            for sample in samples:
                if len(sample[2]) > 0:
                    if sample[1] is not None and sample[2][i][1] is not None:
                        x_values.append([sample[1], sample[3], sample[2][j][1]])
                        if sample[0] == "random":
                            colors.append("red")
                            labels.append("Random transition table")
                        elif sample[0] == "GRN":
                            colors.append("blue")
                            labels.append("Biological transition table")
            x_values = np.array(x_values)
            
            title = None
            axes_labels = ["Synergy", "Memory", "Nudge impact"]
            filename = "scatter3D_memory_synergy_resilience.pdf"
            visualize.plot_scatter_3d(x_values, colors, labels, title, os.path.join(folderpath, filename), axes_labels)

    
    # now those that don't
    # always start by finding all experiments
    experiments_dataframe = dataframe[["system_size", "logic_size"]]
    experiments = experiments_dataframe.drop_duplicates().values.tolist()
    # loop over the experiments
    for experiment in experiments:
        # fill the variables
        system_size = experiment[0]
        logic_size = experiment[1]
        # get the wanted data
        a = dataframe["logic_size"] == logic_size
        b = dataframe["system_size"] == system_size
        selected_dataframe = dataframe[a & b]
        samples = selected_dataframe[["type", "synergy", "impacts", "memory"]].values.tolist()

        # PLOT: synergy - memory
        x_values = []
        colors = []
        labels = []
        
        # read the samples
        for sample in samples:
            if sample[1] is not None and sample[3] is not None:
                x_values.append([sample[1], sample[3]])
                if sample[0] == "random":
                    colors.append("red")
                    labels.append("Random transition table")
                elif sample[0] == "GRN":
                    colors.append("blue")
                    labels.append("Biological transition table")

        x_values = np.array(x_values)

        # plot and save
        foldername = "k=%d_l=%d" % (system_size, logic_size)
        folderpath = os.path.join(result_location, foldername)
        if not os.path.isdir(folderpath):
            os.mkdir(folderpath)
        title = None
        axes_labels = ["Synergy", "Memory"]
        filename = "scatter2D_synergy_memory.pdf"
        visualize.plot_scatter(x_values, colors, labels, title, os.path.join(folderpath, filename), axes_labels)


def visualize_profile(dataframe):
    """
    The MI-profile, as discussed in the paper.
    """
    # always start by finding all experiments
    experiments_dataframe = dataframe[["system_size", "logic_size"]]
    experiments = experiments_dataframe.drop_duplicates().values.tolist()
    # loop over the experiments
    for experiment in experiments:
        # fill the variables
        system_size = experiment[0]
        logic_size = experiment[1]

        # get the wanted data
        a = dataframe["logic_size"] == logic_size
        b = dataframe["system_size"] == system_size
        selected_dataframe = dataframe[a & b]
        samples = selected_dataframe[["type", "motif"]].values.tolist()

        # saving location
        foldername = "k=%d_l=%d" % (system_size, logic_size)
        folderpath = os.path.join(result_location, foldername)
        if not os.path.isdir(folderpath):
            os.mkdir(folderpath)

        random_motifs = []
        for sample in samples:
            if sample[0] == "random":
                random_motifs.append(sample[1])
        #title = "MI Profile (random tables) with motif size %s and %s-valued logic, nudge size of %s" % (system_size, logic_size, nudge_size)
        title = None
        filename = os.path.join(folderpath, "MIprofile_random.pdf")
        visualize.plot_mi_profile(random_motifs, title, mode='maximum', filename = os.path.join(folderpath, filename))

        grn_motifs = []
        for sample in samples:
            if sample[0] == "GRN":
                grn_motifs.append(sample[1])
        #title = "MI Profile (GRN motifs) with motif size %s and %s-valued logic, nudge size of %s" % (system_size, logic_size, nudge_size)
        title = None
        filename = os.path.join(folderpath, "MIprofile_GRN.pdf")
        visualize.plot_mi_profile(grn_motifs, title, mode='maximum', filename = os.path.join(folderpath, filename))


def visualize_impacts(dataframe):
    """
    """
    # always start by finding all experiments
    experiments_dataframe = dataframe[["system_size", "logic_size", "nudge_size"]]
    experiments = experiments_dataframe.drop_duplicates().values.tolist()
    # loop over the experiments
    for experiment in experiments:
        # fill the variables
        system_size = experiment[0]
        logic_size = experiment[1]
        nudge_size = experiment[2]

        # get the wanted data
        a = dataframe["logic_size"] == logic_size
        b = dataframe["system_size"] == system_size
        c = dataframe["nudge_size"] == nudge_size
        selected_dataframe = dataframe[a & b & c]
        samples = selected_dataframe[["type", "synergy", "impacts", "memory"]].values.tolist()

        # plot parameters
        values = []
        colors = []
        labels = []
        for j in tqdm(range(0, int(system_size))):
            x_value = j
            y_values_random = []
            y_values_grn = []
            # the random tables
            for sample in samples:
                if sample[2][j][1] is not None:
                    if sample[0] == "random":
                        y_values_random.append(float(sample[2][j][1]))
                    elif sample[0] == "GRN":
                        y_values_grn.append(float(sample[2][j][1]))

            values.append([x_value, y_values_random])
            colors.append("red")
            labels.append("Random transition table")
            values.append([x_value, y_values_grn])
            colors.append("blue")
            labels.append("Biological transition table")
    
        # save the plot
        foldername = "k=%d_l=%d_e=%f" % (system_size, logic_size, nudge_size)
        folderpath = os.path.join(result_location, foldername)
        if not os.path.isdir(folderpath):
            os.mkdir(folderpath)
        filename = os.path.join(folderpath, "impacts.pdf")
        #title = "Nudge impact vs. Nudge width with motif size %s and %s-valued logic, nudge size of %s" % (system_size, logic_size, nudge_size)
        title = None
        axes_labels = ["Nudge width","Nudge impact"]
        visualize.plot_line(values, colors, labels, title, filename=filename, axes_labels=axes_labels)


def test_synergy(dataframe):
    """
    """
    # always start by finding all experiments
    experiments_dataframe = dataframe[["system_size", "logic_size"]]
    experiments = experiments_dataframe.drop_duplicates().values.tolist()

    # loop over the experiments
    for experiment in experiments:
        # fill the variables
        system_size = experiment[0]
        logic_size = experiment[1]

        # get the wanted data
        a = dataframe["logic_size"] == logic_size
        b = dataframe["system_size"] == system_size
        selected_dataframe = dataframe[a & b]
        samples = selected_dataframe[["type", "synergy"]].values.tolist()

        random_synergies = []
        bio_synergies = []
        for sample in samples:
            if sample[0] == "random":
                random_synergies.append(float(sample[1]))
            elif sample[0] == "GRN":
                bio_synergies.append(float(sample[1]))

        t, prob = scipy.stats.wilcoxon(random_synergies, bio_synergies)
        print("test Z = %.22f and p = %.22f" % (t, prob)) 
        random_mean = np.average(random_synergies)
        bio_mean = np.average(bio_synergies)
        
        args = (system_size, logic_size, random_mean, bio_mean, t, prob)
        if prob < 0.05:
            result = "Using %s nodes and %s-valued logic, we found a significant difference between the mean synergy"\
                    " in random transition tables (%s) and biological transition"\
                    " table (%s), with Z=%s and p=%s.\n" % args
        else:
            result = "Using %s nodes and %s-valued logic, we found no significant difference between the mean synergy"\
                    " in random transition tables (%s) and biological transition"\
                    " table (%s), with Z=%s and p=%s.\n" % args
        
        # do some normality testing
        both_normal = True
        _, pnormal = scipy.stats.normaltest(random_synergies)
        if pnormal > 0.05:
            print("non-normal data encountered")
            result += "Random synergies are not normal (thankfully we use Wilcox rank test)\n"
            both_normal = False
        _, pnormal = scipy.stats.normaltest(bio_synergies)
        if pnormal > 0.05:
            print("non-normal data encountered")
            result += "GRN synergies are not normal (thankfully we use Wilcox rank test)\n"
            both_normal = False
        if both_normal:
            mu1, std1 = scipy.stats.norm.fit(random_synergies)
            mu2, std2 = scipy.stats.norm.fit(bio_synergies)
            result += "found sigmas %.2f (random) and %.2f (GRN)\n" % (std1, std2)
            print("found sigmas %.2f (random) and %.2f (GRN)" % (std1, std2))
        outs = outliers_iqr(random_synergies)
        if len(outs) > 0:
            result += "found outliers in random synergies %s\n" % str(outs)
            print("found outliers %s" % str(outs))
        outs = outliers_iqr(bio_synergies)
        if len(outs) > 0:
            result += "found outliers in GRNS synergies %s\n" % str(outs)
            print("found outliers %s" % str(outs))

            
        foldername = "k=%d_l=%d" % (system_size, logic_size)
        folderpath = os.path.join(result_location, foldername)
        if not os.path.isdir(folderpath):
            os.mkdir(folderpath)
        with open(os.path.join(folderpath, "more_synergy.txt"), "w") as text_file:
            text_file.write(result)


def test_memory(dataframe):
    """
    As a very simple memory measure, we use the mutual information between the first and second state, as a fraction of the largest of the two's entropies.
    """
    # always start by finding all experiments
    experiments_dataframe = dataframe[["system_size", "logic_size"]]
    experiments = experiments_dataframe.drop_duplicates().values.tolist()

    # loop over the experiments
    for experiment in experiments:
        # fill the variables
        system_size = experiment[0]
        logic_size = experiment[1]

        # get the wanted data
        a = dataframe["logic_size"] == logic_size
        b = dataframe["system_size"] == system_size
        selected_dataframe = dataframe[a & b]
        samples = selected_dataframe[["type", "memory"]].values.tolist()

        random_memories = []
        bio_memories = []
        for sample in samples:
            if sample[0] == "random":
                random_memories.append(float(sample[1]))
            elif sample[0] == "GRN":
                bio_memories.append(float(sample[1]))

        random_mean = np.average(random_memories)
        bio_mean = np.average(bio_memories)
        t, prob = scipy.stats.wilcoxon(random_memories, bio_memories)
        if t == 0.0:
            print("apparently we get z = 0 if there is no overlap between the two distributions")
            random_memories.sort()
            bio_memories.sort()
            print(random_memories[:10])
            print(bio_memories[-10:])
        print("test Z = %.22f and p = %.22f" % (t, prob)) 
        
        args = (system_size, logic_size, random_mean, bio_mean, t, prob)
        if prob < 0.05:
            result = "Using %s nodes and %s-valued logic, we found a significant difference between the mean memory"\
                    " in random transition tables (%s) and biological transition"\
                    " table (%s), with Z=%s and p=%s.\n" % args
        else:
            result = "Using %s nodes and %s-valued logic, we found no significant difference between the mean memory"\
                    " in random transition tables (%s) and biological transition"\
                    " table (%s), with Z=%s and p=%s.\n" % args
        
        # do some normality testing
        both_normal = True
        _, pnormal = scipy.stats.normaltest(random_memories)
        if pnormal > 0.05:
            print("non-normal data encountered")
            result += "Random memories are not normal (thankfully we use Wilcox rank test)\n"
            both_normal = False
        _, pnormal = scipy.stats.normaltest(bio_memories)
        if pnormal > 0.05:
            print("non-normal data encountered")
            result += "GRN memories are not normal (thankfully we use Wilcox rank test)\n"
            both_normal = False
        if both_normal:
            mu1, std1 = scipy.stats.norm.fit(random_memories)
            mu2, std2 = scipy.stats.norm.fit(bio_memories)
            result += "found sigmas %.2f (random) and %.2f (GRN)\n" % (std1, std2)
            print("found sigmas %.2f (random) and %.2f (GRN)" % (std1, std2))
        outs = outliers_iqr(random_memories)
        if len(outs) > 0:
            result += "found outliers in random memories %s\n" % str(outs)
            print("found outliers %s" % str(outs))
        outs = outliers_iqr(bio_memories)
        if len(outs) > 0:
            result += "found outliers in GRNS memories %s\n" % str(outs)
            print("found outliers %s" % str(outs))
            
        foldername = "k=%d_l=%d" % (system_size, logic_size)
        folderpath = os.path.join(result_location, foldername)
        if not os.path.isdir(folderpath):
            os.mkdir(folderpath)
        with open(os.path.join(folderpath, "more_memory.txt"), "w") as text_file:
            text_file.write(result)


def test_resilience(dataframe):
    """
    """
    # always start by finding all experiments
    experiments_dataframe = dataframe[["system_size", "logic_size", "nudge_size"]]
    experiments = experiments_dataframe.drop_duplicates().values.tolist()
    # loop over the experiments
    for experiment in experiments:
        # fill the variables
        system_size = experiment[0]
        logic_size = experiment[1]
        nudge_size = experiment[2]

        # get the wanted data
        a = dataframe["logic_size"] == logic_size
        b = dataframe["system_size"] == system_size
        c = dataframe["nudge_size"] == nudge_size
        selected_dataframe = dataframe[a & b & c]
        samples = selected_dataframe[["type", "impacts"]].values.tolist()
            
        # check if folder exists
        foldername = "k=%d_l=%d_e=%f" % (system_size, logic_size, nudge_size)
        folderpath = os.path.join(result_location, foldername)
        if not os.path.isdir(folderpath):
            os.mkdir(folderpath)

        # EXPERIMENT: single resilience
        random_resiliences = []
        bio_resiliences = []
        for sample in samples:
            if sample[0] == "random":
                random_resiliences.append(float(sample[1][0][1]))
            elif sample[0] == "GRN":
                bio_resiliences.append(float(sample[1][0][1]))
        t, prob = scipy.stats.wilcoxon(random_resiliences, bio_resiliences)
        print("test Z = %.22f and p = %.22f" % (t, prob)) 
        random_mean = np.average(random_resiliences)
        bio_mean = np.average(bio_resiliences)
        args = (system_size, logic_size, nudge_size, "1", random_mean, bio_mean, t, prob)
        if prob < 0.05:
            result = "Using %s nodes, %s-valued logic, and %s-epsilon %s-target nudge, we found a significant"\
                    " difference between the mean nudge impact"\
                    " in random transition tables (%s) and biological transition"\
                    " table (%s), with Z=%s and p=%s.\n" % args
        else:
            result = "Using %s nodes, %s-valued logic, and %s-epsilon %s-target nudge, we found no significant"\
                    " difference between the mean nudge impact"\
                    " in random transition tables (%s) and biological transition"\
                    " table (%s), with Z=%s and p=%s.\n" % args
        
        # do some normality testing
        both_normal = True
        _, pnormal = scipy.stats.normaltest(random_resiliences)
        if pnormal > 0.05:
            print("non-normal data encountered")
            result += "Random resiliences are not normal (thankfully we use Wilcox rank test)\n"
            both_normal = False
        _, pnormal = scipy.stats.normaltest(bio_resiliences)
        if pnormal > 0.05:
            print("non-normal data encountered")
            result += "GRN resiliences are not normal (thankfully we use Wilcox rank test)\n"
            both_normal = False
        if both_normal:
            mu1, std1 = scipy.stats.norm.fit(random_resiliences)
            mu2, std2 = scipy.stats.norm.fit(bio_resiliences)
            result += "found sigmas %.2f (random) and %.2f (GRN)\n" % (std1, std2)
            print("found sigmas %.2f (random) and %.2f (GRN)" % (std1, std2))
        outs = outliers_iqr(random_resiliences)
        if len(outs) > 0:
            result += "found outliers in random resiliences %s\n" % str(outs)
            print("found outliers %s" % str(outs))
        outs = outliers_iqr(bio_resiliences)
        if len(outs) > 0:
            result += "found outliers in GRNS resiliences %s\n" % str(outs)
            print("found outliers %s" % str(outs))

        with open(os.path.join(folderpath, "more_resilience_single.txt"), "w") as text_file:
            text_file.write(result)

        # EXPERIMENT: multiple resilience
        random_resiliences = []
        bio_resiliences = []
        for sample in samples:
            if sample[0] == "random":
                random_resiliences.append(float(sample[1][-1][1]))
            elif sample[0] == "GRN":
                bio_resiliences.append(float(sample[1][-1][1]))
        t, prob = scipy.stats.wilcoxon(random_resiliences, bio_resiliences)
        print("test Z = %.22f and p = %.22f" % (t, prob))
        random_mean = np.average(random_resiliences)
        bio_mean = np.average(bio_resiliences)
        args = (system_size, logic_size, nudge_size, str(len(samples[0][1])), random_mean, bio_mean, t, prob)
        if prob < 0.05:
            result = "Using %s nodes, %s-valued logic, and %s-epsilon %s-target nudge, we found a significant"\
                    " difference between the mean nudge impact"\
                    " in random transition tables (%s) and biological transition"\
                    " table (%s), with Z=%s and p=%s.\n" % args
        else:
            result = "Using %s nodes, %s-valued logic, and %s-epsilon %s-target nudge, we found no significant"\
                    " difference between the mean nudge impact"\
                    " in random transition tables (%s) and biological transition"\
                    " table (%s), with Z=%s and p=%s.\n" % args
        
        # do some normality testing
        both_normal = True
        _, pnormal = scipy.stats.normaltest(random_resiliences)
        if pnormal > 0.05:
            print("non-normal data encountered")
            result += "Random resiliences are not normal (thankfully we use Wilcox rank test)\n"
            both_normal = False
        _, pnormal = scipy.stats.normaltest(bio_resiliences)
        if pnormal > 0.05:
            print("non-normal data encountered")
            result += "GRN resiliences are not normal (thankfully we use Wilcox rank test)\n"
            both_normal = False
        if both_normal:
            mu1, std1 = scipy.stats.norm.fit(random_resiliences)
            mu2, std2 = scipy.stats.norm.fit(bio_resiliences)
            result += "found sigmas %.2f (random) and %.2f (GRN)\n" % (std1, std2)
            print("found sigmas %.2f (random) and %.2f (GRN)" % (std1, std2))
        outs = outliers_iqr(random_resiliences)
        if len(outs) > 0:
            result += "found outliers in random resiliences %s\n" % str(outs)
            print("found outliers %s" % str(outs))
        outs = outliers_iqr(bio_resiliences)
        if len(outs) > 0:
            result += "found outliers in GRNS resiliences %s\n" % str(outs)
            print("found outliers %s" % str(outs))


        with open(os.path.join(folderpath, "more_resilience_multiple.txt"), "w") as text_file:
            text_file.write(result)


def create_table_3v(df, experiment, caption):
    outfile = open(os.path.join(result_location, experiment + ".tex"), "w")
    outfile.write(r"\begin{table}[h]" + "\n")
    
    # find the column definitions
    column_cnt = 2 + len(list(df.nudge_size.unique()))
    column_definition = r"{|"
    for _ in range(0, column_cnt):
        column_definition += r"c|"
    column_definition += r"}"
    
    # begin tabular
    outfile.write("\n" + r"\begin{tabular}" + column_definition + "\n")
    
    # find the unique values in order for all things
    nudge_sizes = list(set(df.nudge_size))
    system_sizes = list(set(df.system_size))
    logic_sizes = list(set(df.logic_size))
    nudge_sizes.sort()
    system_sizes.sort()
    logic_sizes.sort()
    
    # header row
    outfile.write(r"\hline" + "\n")
    header = r"\# nodes & \diagbox{\# states}{$\epsilon$} "
    for nudge_size in nudge_sizes:
        header += r" & " + str(nudge_size)
    header += r"\\"
    outfile.write(header + "\n")
    outfile.write(r"\hline" + "\n")
    
    # loop over rows
    for system_size in system_sizes:
        first = True
        for logic_size in logic_sizes:
            # new row
            row = ""
            if first:
                first = False
                row += r"\multirow{" + str(len(logic_sizes)) + r"}{*}{" + str(system_size) + r"}"
            else:
                row += " "
                
            row += " & " + str(logic_size)
            for nudge_size in nudge_sizes:
                # get the p-value
                query = "logic_size == " + str(logic_size) + " | system_size == " + str(system_size) + " | nudge_size == " + str(nudge_size)
                a = df["logic_size"] == logic_size
                b = df["system_size"] == system_size
                c = df["nudge_size"] == nudge_size
                d = df["experiment"] == experiment
                print_row = len(list(df[a & b & c & d].p_value.unique())) != 0
                p_value = list(df[a & b & c & d].p_value.unique())
                z_value = list(df[a & b & c & d].statistic_value.unique())
                color = list(df[a & b & c & d].color.unique())
                if len(p_value) > 0:
                    # truncate
                    p_value = "%.3e" % p_value[0]
                    color = color[0]
                    
                    # add column
                    row += r" & %.2f" % z_value[0]
                    
                    # add stars based on significance
                    if color is not None:
                        if float(p_value) < 0.001:
                            row += r"*** \cellcolor{%s!20}" % color
                        elif float(p_value) < 0.005:
                            row += r"** \cellcolor{%s!20}" % color
                        elif float(p_value) < 0.05:
                            row += r"* \cellcolor{%s!20}" % color
                    else:
                        if float(p_value) < 0.001:
                            row += r"*** "
                        elif float(p_value) < 0.005:
                            row += r"** " 
                        elif float(p_value) < 0.05:
                            row += r"* " 
                else:
                    row += r" & "
            outfile.write(row + r"\\" + "\n")
            outfile.write(r"\cline{2-5}" + "\n")
        outfile.write(r"\hline" + "\n")
        
    outfile.write(r"\end{tabular}" + "\n")
    outfile.write(r"\centering" + "\n")
    outfile.write(caption)
    outfile.write(r"\label{%s}" % experiment)
    outfile.write("\n" + r"\end{table}"+ "\n")
    outfile.close()
    return 0


def create_table_2v(df, experiment, caption):
    outfile = open(os.path.join(result_location, experiment + ".tex"), "w")
    outfile.write(r"\begin{table}[h]" + "\n")
    
    # find the column definitions
    column_cnt = 1 + len(list(df.logic_size.unique()))
    column_definition = r"{|"
    for _ in range(0, column_cnt):
        column_definition += r"c|"
    column_definition += r"}"
    
    # begin tabular
    outfile.write("\n" + r"\begin{tabular}" + column_definition + "\n")
    
    # find the unique values in order for all things
    system_sizes = list(set(df.system_size))
    logic_sizes = list(set(df.logic_size))
    system_sizes.sort()
    logic_sizes.sort()
    
    # header row
    outfile.write(r"\hline" + "\n")
    header = r"\diagbox{\# nodes }{\# states} "
    for logic_size in logic_sizes:
        header += r" & " + str(logic_size)
    header += r"\\"
    outfile.write(header + "\n")
    outfile.write(r"\hline" + "\n")
    
    # loop over rows
    for system_size in system_sizes:
        # new row
        row = ""
            
        row += str(system_size)
        for logic_size in logic_sizes:
            # get the p-value
            query = "logic_size == " + str(logic_size) + " | system_size == " + str(system_size)
            a = df["logic_size"] == logic_size
            b = df["system_size"] == system_size
            d = df["experiment"] == experiment
            print_row = len(list(df[a & b & d].p_value.unique())) != 0
            p_value = list(df[a & b & d].p_value.unique())
            z_value = list(df[a & b & d].statistic_value.unique())
            color = list(df[a & b & d].color.unique())
            if len(p_value) > 0:
                # truncate
                p_value = "%.3e" % p_value[0]
                color = color[0]
                
                # add column
                row += r" & %.2f" % z_value[0]
                
                # add stars based on significance
                if color is not None:
                    if float(p_value) < 0.001:
                        row += r"*** \cellcolor{%s!20}" % color
                    elif float(p_value) < 0.005:
                        row += r"** \cellcolor{%s!20}" % color
                    elif float(p_value) < 0.05:
                        row += r"* \cellcolor{%s!20}" % color
                else:
                    if float(p_value) < 0.001:
                        row += r"*** "
                    elif float(p_value) < 0.005:
                        row += r"** " 
                    elif float(p_value) < 0.05:
                        row += r"* " 
            else:
                row += r" & "
        outfile.write(row + r"\\" + "\n")
        outfile.write(r"\hline" + "\n")
        
    outfile.write(r"\end{tabular}" + "\n")
    outfile.write(r"\centering" + "\n")
    outfile.write(caption)
    outfile.write(r"\label{%s}" % experiment)
    outfile.write("\n" + r"\end{table}"+ "\n")
    outfile.close()
    return 0


def main():
    if not os.path.isdir(result_location):
        os.mkdir(result_location)

    # get the parameters, unpack and store
    sample_size = -1
    with open(os.path.join(data_location, "parameters.pkl"), 'rb') as input:
        parameters = pickle.load(input)
        outfile = open(os.path.join(result_location, "parameters.txt"), "w")
        print(parameters)
        outfile.write("sample size: %d \n" % parameters["sample_size"])
        outfile.write("system sizes: %s \n" % str(parameters["network_sizes"]))
        outfile.write("logic sizes: %s \n" % str(parameters["logic_sizes"]))
        outfile.write("nudge sizes: %s \n" % str(parameters["nudge_sizes"]))
        outfile.write("nudge method: %s \n" % str(parameters["nudge_method"]))
        outfile.write("synergy measure: %s \n" % str(parameters["synergy_measure"]))
        sample_size = parameters["sample_size"]

    # load data back, append the dataframes
    dataframe = pd.DataFrame(columns=["system_size", "logic_size", "nudge_size", "type", "motif", "synergy", "memory", "impacts"])
    for root, dirnames, filenames in os.walk(data_location):
        for filename in filenames:
            if filename.endswith('df.pkl'):
                with open(os.path.join(root, filename), 'rb') as input:
                    print("opening %s" % input)
                    dataframe_new = pickle.load(input)
                    dataframe = pd.concat([dataframe, dataframe_new])

    # now let's get rolling and do our tests!
    print("doing impacts")
    visualize_impacts(dataframe)
    print("doing profile")
    visualize_profile(dataframe)
    print("doing cyclesearch")
    histogram_cycles(dataframe)
    sys.exit(1)
    print("proceeding to LaTeX")
    print("doing TSNE")
    visualize_TSNE(dataframe)
    print("doing tests")
    test_synergy(dataframe)
    test_memory(dataframe)
    test_resilience(dataframe)
    print("doing scatters")
    visualize_scatters(dataframe)

    # now also do our spearman tests
    # immediatiely put in LaTeX
    # we use the same table as later, to reuse tables
    df_2 = pd.DataFrame(columns=["experiment", "system_size", "logic_size", "p_value", "statistic_value", "color"])
    loc_df2 = 0
    df_3 = pd.DataFrame(columns=["experiment", "system_size", "logic_size", "nudge_size", "p_value", "statistic_value", "color"])
    loc_df3 = 0
    
    # do the loop
    nudge_sizes = list(set(dataframe.nudge_size))
    system_sizes = list(set(dataframe.system_size))
    logic_sizes = list(set(dataframe.logic_size))
    for system_size in system_sizes:
        for logic_size in logic_sizes:
            # correlate memory synergy
            a = dataframe["logic_size"] == logic_size
            b = dataframe["system_size"] == system_size
            c = dataframe["type"] == "random"
            selected_dataframe = dataframe[a & b & c]
            samples = selected_dataframe[["type", "synergy", "memory"]].values.tolist()
            synergies = []
            memories = []
            for sample in samples:
                if sample[0] == "random":
                    synergies.append(sample[1])
                    memories.append(sample[2])
            rho, p_value = scipy.stats.spearmanr(synergies, memories)
            df_2.loc[loc_df2] = ["random_rho_syn_mem", system_size, logic_size, p_value, rho, None]
            loc_df2 += 1

            # now do the nudge-involving experiments
            for nudge_size in nudge_sizes:
                # correlate synergy single nudge
                a = dataframe["logic_size"] == logic_size
                b = dataframe["system_size"] == system_size
                c = dataframe["type"] == "random"
                d = dataframe["nudge_size"] == nudge_size
                selected_dataframe = dataframe[a & b & d]
                samples = selected_dataframe[["type", "synergy", "impacts"]].values.tolist()
                synergies = []
                impacts = []
                for sample in samples:
                    if sample[0] == "random":
                        synergies.append(sample[1])
                        impacts.append(sample[2][0][1])
                rho, p_value = scipy.stats.spearmanr(synergies, impacts)
                df_3.loc[loc_df3] = ["random_rho_syn_singleimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1

                # correlate synergy multiple nudge
                samples = selected_dataframe[["type", "synergy", "impacts"]].values.tolist()
                synergies = []
                impacts = []
                for sample in samples:
                    if sample[0] == "random":
                        synergies.append(sample[1])
                        impacts.append(sample[2][-1][1])
                rho, p_value = scipy.stats.spearmanr(synergies, impacts)
                df_3.loc[loc_df3] = ["random_rho_syn_multimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1

                # correlate memory single nudge
                samples = selected_dataframe[["type", "memory", "impacts"]].values.tolist()
                memories = []
                impacts = []
                for sample in samples:
                    if sample[0] == "random":
                        memories.append(sample[1])
                        impacts.append(sample[2][0][1])
                rho, p_value = scipy.stats.spearmanr(memories, impacts)
                df_3.loc[loc_df3] = ["random_rho_mem_singleimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1

                # correlate memory multiple nudge
                samples = selected_dataframe[["type", "memory", "impacts"]].values.tolist()
                memories = []
                impacts = []
                for sample in samples:
                    if sample[0] == "random":
                        memories.append(sample[1])
                        impacts.append(sample[2][-1][1])
                rho, p_value = scipy.stats.spearmanr(memories, impacts)
                df_3.loc[loc_df3] = ["random_rho_mem_multimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1

                # do partial correlation synergy single nudge
                # do partial correlation synergy multiple nudge
                samples = selected_dataframe[["type", "synergy", "memory", "impacts"]].values.tolist()
                synergies = []
                memories = []
                impacts_single = []
                impacts_multiple = []
                for sample in samples:
                    if sample[0] == "random":
                        synergies.append(sample[1])
                        memories.append(sample[2])
                        impacts_single.append(sample[3][0][1])
                        impacts_multiple.append(sample[3][-1][1])
                # do linear models for residuals
                if len(synergies) == 0:
                    continue
                lr1 = LinearRegression()
                lr1.fit(np.asarray(memories).reshape(len(synergies), 1), synergies)
                lr2 = LinearRegression()
                lr2.fit(np.asarray(memories).reshape(len(impacts_single), 1), impacts_single)
                lr3 = LinearRegression()
                lr3.fit(np.asarray(memories).reshape(len(impacts_multiple), 1), impacts_multiple)
                # calculate residuals
                residuals_synergy = []
                residuals_single_impact = []
                residuals_multiple_impact = []
                for i in range(0, len(synergies)):
                    residuals_synergy.append(synergies[i] - lr1.predict(memories[i])[0])
                    residuals_single_impact.append(impacts_single[i] - lr2.predict(memories[i])[0])
                    residuals_multiple_impact.append(impacts_multiple[i] - lr3.predict(memories[i])[0])
                # do partial correlations
                rho, p_value = scipy.stats.spearmanr(residuals_synergy, residuals_single_impact)
                df_3.loc[loc_df3] = ["random_rho_partial_synergy_singleimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1
                rho, p_value = scipy.stats.spearmanr(residuals_synergy, residuals_multiple_impact)
                df_3.loc[loc_df3] = ["random_rho_partial_synergy_multimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1

                # do partial correlation memory single nudge
                # do partial correlation memory multiple nudge
                # do linear models for residuals
                if len(synergies) == 0:
                    continue
                lr1 = LinearRegression()
                lr1.fit(np.asarray(synergies).reshape(len(memories), 1), memories)
                lr2 = LinearRegression()
                lr2.fit(np.asarray(synergies).reshape(len(impacts_single), 1), impacts_single)
                lr3 = LinearRegression()
                lr3.fit(np.asarray(synergies).reshape(len(impacts_multiple), 1), impacts_multiple)
                # calculate residuals
                residuals_memory = []
                residuals_single_impact = []
                residuals_multiple_impact = []
                for i in range(0, len(synergies)):
                    residuals_memory.append(memories[i] - lr1.predict(synergies[i])[0])
                    residuals_single_impact.append(impacts_single[i] - lr2.predict(synergies[i])[0])
                    residuals_multiple_impact.append(impacts_multiple[i] - lr3.predict(synergies[i])[0])
                # do partial correlations
                rho, p_value = scipy.stats.spearmanr(residuals_memory, residuals_single_impact)
                df_3.loc[loc_df3] = ["random_rho_partial_memory_singleimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1
                rho, p_value = scipy.stats.spearmanr(residuals_memory, residuals_multiple_impact)
                df_3.loc[loc_df3] = ["random_rho_partial_memory_multimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1


                # do partial correlation synergy single nudge
                # do partial correlation synergy multiple nudge
                samples = selected_dataframe[["type", "synergy", "memory", "impacts"]].values.tolist()
                synergies = []
                memories = []
                impacts_single = []
                impacts_multiple = []
                for sample in samples:
                    if sample[0] == "GRN":
                        synergies.append(sample[1])
                        memories.append(sample[2])
                        impacts_single.append(sample[3][0][1])
                        impacts_multiple.append(sample[3][-1][1])
                # do linear models for residuals
                if len(synergies) == 0:
                    continue
                lr1 = LinearRegression()
                lr1.fit(np.asarray(memories).reshape(len(synergies), 1), synergies)
                lr2 = LinearRegression()
                lr2.fit(np.asarray(memories).reshape(len(impacts_single), 1), impacts_single)
                lr3 = LinearRegression()
                lr3.fit(np.asarray(memories).reshape(len(impacts_multiple), 1), impacts_multiple)
                # calculate residuals
                residuals_synergy = []
                residuals_single_impact = []
                residuals_multiple_impact = []
                for i in range(0, len(synergies)):
                    residuals_synergy.append(synergies[i] - lr1.predict(memories[i])[0])
                    residuals_single_impact.append(impacts_single[i] - lr2.predict(memories[i])[0])
                    residuals_multiple_impact.append(impacts_multiple[i] - lr3.predict(memories[i])[0])
                # do partial correlations
                rho, p_value = scipy.stats.spearmanr(residuals_synergy, residuals_single_impact)
                df_3.loc[loc_df3] = ["GRN_rho_partial_synergy_singleimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1
                rho, p_value = scipy.stats.spearmanr(residuals_synergy, residuals_multiple_impact)
                df_3.loc[loc_df3] = ["GRN_rho_partial_synergy_multimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1

                # do partial correlation memory single nudge
                # do partial correlation memory multiple nudge
                # do linear models for residuals
                if len(synergies) == 0:
                    continue
                lr1 = LinearRegression()
                lr1.fit(np.asarray(synergies).reshape(len(memories), 1), memories)
                lr2 = LinearRegression()
                lr2.fit(np.asarray(synergies).reshape(len(impacts_single), 1), impacts_single)
                lr3 = LinearRegression()
                lr3.fit(np.asarray(synergies).reshape(len(impacts_multiple), 1), impacts_multiple)
                # calculate residuals
                residuals_memory = []
                residuals_single_impact = []
                residuals_multiple_impact = []
                for i in range(0, len(synergies)):
                    residuals_memory.append(memories[i] - lr1.predict(synergies[i])[0])
                    residuals_single_impact.append(impacts_single[i] - lr2.predict(synergies[i])[0])
                    residuals_multiple_impact.append(impacts_multiple[i] - lr3.predict(synergies[i])[0])
                # do partial correlations
                rho, p_value = scipy.stats.spearmanr(residuals_memory, residuals_single_impact)
                df_3.loc[loc_df3] = ["GRN_rho_partial_memory_singleimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1
                rho, p_value = scipy.stats.spearmanr(residuals_memory, residuals_multiple_impact)
                df_3.loc[loc_df3] = ["GRN_rho_partial_memory_multimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1

    for system_size in system_sizes:
        for logic_size in logic_sizes:
            # correlate memory synergy
            a = dataframe["logic_size"] == logic_size
            b = dataframe["system_size"] == system_size
            c = dataframe["type"] == "GRN"
            selected_dataframe = dataframe[a & b & c]
            samples = selected_dataframe[["type", "synergy", "memory"]].values.tolist()
            synergies = []
            memories = []
            for sample in samples:
                if sample[0] == "GRN":
                    synergies.append(sample[1])
                    memories.append(sample[2])
            rho, p_value = scipy.stats.spearmanr(synergies, memories)
            df_2.loc[loc_df2] = ["GRN_rho_syn_mem", system_size, logic_size, p_value, rho, None]
            loc_df2 += 1

            # now do the nudge-involving experiments
            for nudge_size in nudge_sizes:
                # correlate synergy single nudge
                a = dataframe["logic_size"] == logic_size
                b = dataframe["system_size"] == system_size
                c = dataframe["type"] == "GRN"
                d = dataframe["nudge_size"] == nudge_size
                selected_dataframe = dataframe[a & b & c & d]
                samples = selected_dataframe[["type", "synergy", "impacts"]].values.tolist()
                synergies = []
                impacts = []
                for sample in samples:
                    if sample[0] == "GRN":
                        synergies.append(sample[1])
                        impacts.append(sample[2][0][1])
                rho, p_value = scipy.stats.spearmanr(synergies, impacts)
                df_3.loc[loc_df3] = ["GRN_rho_syn_singleimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1

                # correlate synergy multiple nudge
                samples = selected_dataframe[["type", "synergy", "impacts"]].values.tolist()
                synergies = []
                impacts = []
                for sample in samples:
                    if sample[0] == "GRN":
                        synergies.append(sample[1])
                        impacts.append(sample[2][-1][1])
                rho, p_value = scipy.stats.spearmanr(synergies, impacts)
                df_3.loc[loc_df3] = ["GRN_rho_syn_multimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1

                # correlate memory single nudge
                samples = selected_dataframe[["type", "memory", "impacts"]].values.tolist()
                memories = []
                impacts = []
                for sample in samples:
                    if sample[0] == "GRN":
                        memories.append(sample[1])
                        impacts.append(sample[2][0][1])
                rho, p_value = scipy.stats.spearmanr(memories, impacts)
                df_3.loc[loc_df3] = ["GRN_rho_mem_singleimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1

                # correlate memory multiple nudge
                samples = selected_dataframe[["type", "memory", "impacts"]].values.tolist()
                memories = []
                impacts = []
                for sample in samples:
                    if sample[0] == "GRN":
                        memories.append(sample[1])
                        impacts.append(sample[2][-1][1])
                rho, p_value = scipy.stats.spearmanr(memories, impacts)
                df_3.loc[loc_df3] = ["GRN_rho_mem_multimpact", system_size, logic_size, nudge_size, p_value, rho, None]
                loc_df3 += 1
            
    # create the tables
    compression_factor = len(list(set(df_3.nudge_size)))
    for experiment in list(df_3.experiment.unique()):
        caption = re.sub("_", " ", r"\caption{Experiment %s, $r_S$ and significance per experiment for random/GRN tables (* implies $p<0.05$, ** $p<0.005$, *** $p<0.0005$) with n=%s.}" % (experiment, sample_size))
        create_table_3v(df_3, experiment, caption)
    for experiment in list(df_2.experiment.unique()):
        sample_size_2v = compression_factor * sample_size
        caption = re.sub("_", " ", r"\caption{Experiment %s, $r_S$ and significance per experiment for random/GRN tables (* implies $p<0.05$, ** $p<0.005$, *** $p<0.0005$) with n=%s.}" % (experiment, sample_size_2v))
        create_table_2v(df_2, experiment, caption)

    # now create our latex tables of the wilcoxon rank test
    df_3 = pd.DataFrame(columns=["experiment", "system_size", "logic_size", "nudge_size", "p_value", "statistic_value", "color"])
    df_2 = pd.DataFrame(columns=["experiment", "system_size", "logic_size", "p_value", "statistic_value", "color"])

    loc_counter = 0
    
    for root, dirs, files in os.walk(result_location):
        for file in files:
            if file.endswith('.txt') and not file.endswith('parameters.txt') and (re.search('k=[0-9]+_l=[0-9]+_e=[0-9]+.[0-9]+', root) is not None):
                # this is an experiment, to be added to our tables
                cgs = list(re.findall('k=([0-9]+)_l=([0-9]+)_e=([0-9]+.[0-9]+)', root)[0])
                experiment = re.findall('more_([a-z _]+)', file)
                with open(os.path.join(root, file), 'r') as f:
                    result = f.readline()
                    z_value = re.findall('Z=([0-9]+.[0-9]+e?-?[0-9]*)', result)
                    p_value = re.findall('p=([0-9]+.[0-9]+e?-?[0-9]*)', result)
                    color = None
                    if len(re.findall("no significant", result)) > 0:
                        color = None
                    else:
                        value_bio = float(re.findall('table \(([0-9]+.[0-9]+)\)', result)[0])
                        value_ran = float(re.findall('tables \(([0-9]+.[0-9]+)\)', result)[0])
                        if value_ran > value_bio:
                            color = "yellow"
                        else:
                            color = "green"
                    df_3.loc[loc_counter] = [experiment[0], int(cgs[0]), int(cgs[1]), float(cgs[2]), float(p_value[0]), float(z_value[0]), color]
                    loc_counter += 1
            elif file.endswith('.txt') and not file.endswith('parameters.txt') and (re.search('k=[0-9]+_l=[0-9]+', root) is not None):
                # this is an experiment, to be added to our tables
                cgs = list(re.findall('k=([0-9]+)_l=([0-9]+)', root)[0])
                experiment = re.findall('more_([a-z _]+)', file)
                with open(os.path.join(root, file), 'r') as f:
                    result = f.readline()
                    z_value = re.findall('Z=([0-9]+.[0-9]+e?-?[0-9]*)', result)
                    p_value = re.findall('p=([0-9]+.[0-9]+e?-?[0-9]*)', result)
                    color = None
                    if len(re.findall("no significant", result)) > 0:
                        color = None
                    else:
                        value_bio = float(re.findall('table \(([0-9]+.[0-9]+)\)', result)[0])
                        value_ran = float(re.findall('tables \(([0-9]+.[0-9]+)\)', result)[0])
                        # TODO: change this to dashed/solid background
                        if value_ran > value_bio:
                            color = "yellow"
                        else:
                            color = "green"
                    df_2.loc[loc_counter] = [experiment[0], int(cgs[0]), int(cgs[1]), float(p_value[0]), float(z_value[0]), color]
                    loc_counter += 1
    
    # create the 4 tables
    compression_factor = len(list(set(df_3.nudge_size)))
    for experiment in list(df_3.experiment.unique()):
        caption = re.sub("_", " ", r"\caption{Experiment %s, Z-value and significance per experiment (* implies $p<0.05$, ** $p<0.005$, *** $p<0.0005$) with n=%s. Green background implies higher mean in the biological network, yellow higher mean in the random network.}" % (experiment, sample_size))
        create_table_3v(df_3, experiment, caption)
    for experiment in list(df_2.experiment.unique()):
        sample_size_2v = compression_factor * sample_size
        caption = re.sub("_", " ", r"\caption{Experiment %s, Z-value and significance per experiment (* implies $p<0.05$, ** $p<0.005$, *** $p<0.0005$) with n=%s. Green background implies higher mean in the biological network, yellow higher mean in the random network.}" % (experiment, sample_size_2v))
        create_table_2v(df_2, experiment, caption)
    
if __name__ == '__main__':
    main()
