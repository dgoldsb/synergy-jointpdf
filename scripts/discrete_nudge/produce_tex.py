#\usepackage{array}
#\usepackage{diagbox}
#\usepackage{multirow}
#\usepackage[table]{xcolor}

import os
import pandas as pd
import re



# set folders
data_location = "../../data_pandas"
result_location = "../../result_pandas"
log_location = "../../log"

def visualize_TSNE(dataframe):
    """
    We use t-SNE to visualize our sample. What we see makes perfect sense: the biological motifs seem to be a subset of the sample space. In higher valued logic systems, the space becomes larger, making the sample size insufficient. As t-SNE also does clustering, it then starts to appear that we have two seperate sample, but this is simply the separation from the very similar biological motifs from the rest.
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
        samples = selected_dataframe[["type", "motif"]].values.tolist()

        # now the real specific deal
        # x_train should be an array of N vectors of length M
        # where N is the sample size and M is the length of the combined new state of the transition table
        # we combine both the random and bio samples, but give them different colors
        x_train = []
        y_color = []
        labels = []

        for sample in samples:
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
        title = "Two-dimensional embedding of transition tables with motif size %s and %s-valued logic" % (system_size, logic_size)
        x_train_embedded = TSNE(n_components=2, perplexity=10, early_exaggeration=1, verbose=0).fit_transform(x_train)
        
        # save the results
        foldername = "n=%d_l=%d_e=%f_%s_%s_nosamples=%d" % (system_size, logic_size, nudge_size, nudge_method, synergy_measure_name, len(samples))
        folderpath = os.path.join(result_location, foldername)
        if not os.path.isdir(folderpath):
            os.mkdir(folderpath)
        filename = "tsne2D.pdf" 
        visualize.plot_scatter(x_train_embedded, y_color, labels, title, os.path.join(folderpath, filename))


def visualize_scatters(dataframe):
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
        samples = selected_dataframe[["type", "synergy", "impacts"]].values.tolist()

        for i in range(0, int(system_size)):
            x_values = []
            colors = []
            labels = []
            title = "Synergy vs. Nudge impact with motif size %s, %s genes targeted, and %s-valued logic" % \
                    (system_size, i+1, logic_size)
            axes_labels = ["Nudge impact", "Synergy"]
           
            # read the samples
            for sample in samples:
                if len(sample[2]) > 0:
                    if sample[1] is not None and sample[2][i][1] is not None:
                        x_values.append([sample[2][i][1], sample[1]])
                        if sample[0] == "random":
                            colors.append("red")
                            labels.append("Random transition table")
                        elif sample[0] == "GRN":
                            colors.append("blue")
                            labels.append("Biological transition table")

            x_values = np.array(x_values)

            foldername = "n=%d_l=%d_e=%f_%s_%s_nosamples=%d" % (system_size, logic_size, nudge_size, nudge_method, synergy_measure_name, len(samples))
            folderpath = os.path.join(result_location, foldername)
            if not os.path.isdir(folderpath):
                os.mkdir(folderpath)

            filename = "scatter2D_synergy_resilience_width_%s.pdf" % str(i+1)
            visualize.plot_scatter(x_values, colors, labels, title, os.path.join(folderpath, filename), axes_labels)

def visualize_profile(dataframe):
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
        samples = selected_dataframe[["type", "motif"]].values.tolist()
def visualize_memory(dataframe):
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
        samples = selected_dataframe[["type", "motif"]].values.tolist()
def test_synergy(dataframe):
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
        samples = selected_dataframe[["type", "motif"]].values.tolist()
def test_memory(dataframe):
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
        samples = selected_dataframe[["type", "motif"]].values.tolist()
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
        samples = selected_dataframe[["type", "motif"]].values.tolist()


def create_table(df, experiment, samplesize):
    outfile = open(experiment + ".txt", "w")
    outfile.write(r"\begin{figure}[h]" + "\n")
    outfile.write(r"\label{%s}" % experiment)
    
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
                color = list(df[a & b & c & d].color.unique())
                if len(p_value) > 0:
                    # truncate
                    p_value = "%.3e" % p_value[0]
                    color = color[0]
                    
                    # add column
                    row += r" & " + p_value
                    
                    # add stars based on significance
                    if float(p_value) < 0.001:
                        row += r"*** \cellcolor{%s!60}" % color
                    elif float(p_value) < 0.005:
                        row += r"** \cellcolor{%s!40}" % color
                    elif float(p_value) < 0.05:
                        row += r"* \cellcolor{%s!20}" % color
                else:
                    row += r" & "
            outfile.write(row + r"\\" + "\n")
            outfile.write(r"\cline{2-5}" + "\n")
        outfile.write(r"\hline" + "\n")
        
    outfile.write(r"\end{tabular}" + "\n")
    outfile.write(r"\centering" + "\n")
    outfile.write(r"\caption{Experiment %s, green implies higher value in biological networks, yellow higher in random networks (n=%s)}" % (experiment, samplesize))
    outfile.write("\n" + r"\end{figure}"+ "\n")
    outfile.close()
    return 0

def main():
    if not os.path.isdir(result_location):
        os.mkdir(result_location)

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
    visualize_TSNE(dataframe)
    visualize_scatters(dataframe)
    #visualize_profile(dataframe)
    #visualize_memory(dataframe)
    #test_synergy(dataframe)
    #test_memory(dataframe)
    #test_resilience(dataframe)

    # now create our latex tables
    df = pd.DataFrame(columns=["experiment", "system_size", "logic_size", "nudge_size", "p_value", "color"])

    loc_counter = 0
    samplesize = "0"
    
    for root, dirs, files in os.walk("../../data"):
        for file in files:
            if file.endswith('.txt'):
                # this is an experiment, to be added to our tables
                cgs = list(re.findall('n=([0-9]+)_l=([0-9]+)_e=([0-9]+.[0-9]+)', root)[0])
                experiment = re.findall('more_([a-z _]+)', file)
                with open(os.path.join(root, file), 'r') as f:
                    result = f.readline()
                    p_value = re.findall('p=([0-9]+.[0-9]+e?-?[0-9]+)', result)
                    color = None
                    if len(re.findall("no significant", result)) > 0:
                        color = "white"
                    else:
                        value_bio = float(re.findall('table \(([0-9]+.[0-9]+)\)', result)[0])
                        value_ran = float(re.findall('tables \(([0-9]+.[0-9]+)\)', result)[0])
                    if value_ran > value_bio:
                        color = "yellow"
                    else:
                        color = "green"
                    df.loc[loc_counter] = [experiment[0], int(cgs[0]), int(cgs[1]), float(cgs[2]), float(p_value[0]), color]
                    loc_counter += 1
                
                samplesize = re.findall('nosamples=([0-9]+)', root)[0]
                    
    
    # create the 4 tables
    for experiment in list(df.experiment.unique()):
        caption = "(n=" + samplesize + ")"
        create_table(df, experiment, samplesize)
    
if __name__ == '__main__':
    main()