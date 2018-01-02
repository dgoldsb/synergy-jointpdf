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
import numpy as np
import pandas as pd
import pickle
import shutil
import sys
import time
from time import gmtime, strftime

# settings for the experiment
synergy_measure = measures.synergy_middleground
nudge_method = "DJ"
sample_size = 150 # in practice this is times two, we draw a random and a GRN sample
network_sizes = [2, 3, 4, 5]
logic_sizes = [2, 3, 4]
nudge_sizes = [0.1, 0.25, 0.5]

# set folders
data_location = "../../data_pandas"
archive_location = "../../archive_pandas"
log_location = "../../log"
    
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
    # logger
    mylogger = logging.getLogger('mylogger')
    handler1 = logging.FileHandler(filename=os.path.join(log_location, 'experiments_%s.log' % strftime("%Y-%m-%d %H:%M:%S", gmtime())), mode='w')
    handler1.setLevel(logging.DEBUG)
    handler1.setFormatter(logging.Formatter('%(asctime)s [%(levelname)s] %(message)s'))
    mylogger.addHandler(handler1)

    # create the file locations and archive
    if not os.path.isdir(data_location):
        os.mkdir(data_location)
    else:
        if not os.path.isdir(archive_location):
            os.mkdir(archive_location)
        archive_folder = os.path.join(archive_location, str(int(time.time())))
        shutil.copytree(data_location, archive_folder)
        shutil.rmtree(data_location)
        os.mkdir(data_location)
        shutil.copy(os.path.join(archive_folder, ".gitignore"), os.path.join(data_location, ".gitignore"))

    # save the experiment parameters
    parameters = {}
    parameters["synergy_measure"] = synergy_measure
    parameters["nudge_method"] = nudge_method
    parameters["sample_size"] = sample_size * 2 # because we take two samples
    parameters["network_sizes"] = network_sizes
    parameters["logic_sizes"] = logic_sizes
    parameters["nudge_sizes"] = nudge_sizes
    with open(os.path.join(data_location, "parameters.pkl"), 'wb') as output:
        pickle.dump(parameters, output, pickle.HIGHEST_PROTOCOL)

    # draw a few completely random samples, with different parameters
    for network_size in network_sizes:
        for logic_size in logic_sizes:
            for nudge_size in nudge_sizes:
                mylogger.debug("sampling %d nodes, %d logic size, %f nudge size, %s as nudge_method, %s as synergy measure" % (network_size, logic_size, nudge_size, nudge_method, synergy_measure))
                start = time.time()
                dataframe, samples_grn, samples_random = draw_sample(sample_size, network_size, logic_size, nudge_size)
                
                # save the data for future use/reruns
                name_df = "experiment_k=%d_l=%d_e=%f_df.pkl" % (network_size, logic_size, nudge_size)
                name_grn = "samples_grn_k=%d_l=%d.pkl" % (network_size, logic_size)
                name_random = "samples_random_k=%d_l=%d.pkl" % (network_size, logic_size)
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

if __name__ == '__main__':
    main()