#\usepackage{array}
#\usepackage{diagbox}
#\usepackage{multirow}
#\usepackage[table]{xcolor}

import os
import pandas as pd
import re

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