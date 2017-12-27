PANDAS
4 plots
ook de image met labels enzo eruit halen, store per experiment ook een image locatie, selecteer de juiste

import os
import pandas as pd
import re

def create_table(table, experiment, caption):
	outfile = open(experiment + ".txt", "w")
	outfile.write("\begin(figure)[h]")
	outfile.write("\label(%s)" % experiment)
	
	# find the column definitions
	column_cnt = 2 + len(list(df.nudge_size.unique()))
	column_definition = "{|"
	for _ in range(0, column_cnt):
		column_definition += "c|"
	column_definition += "}"
	
	# begin tabular
	outfile.write("\begin(tabular)" + colums)
	
	# find the unique values in order for all things
	nudge_sizes = list(df.nudge_size.unique()).sorted()
	system_sizes = list(df.system_size.unique()).sorted()
	logic_sizes = list(df.nudge_size.unique()).sorted()
	
	# header row
	header = "# nodes & \diagbox{# states}{$\epsilon$} "
	for nudge_size in nudge_sizes:
		header += " & " + nudge_size
	header += "\\"
	
	# loop over rows
	for system_size in system_sizes:
		first = True
		for logic_size in logic_size:
			# new row
			row = ""
			if first:
				first = False
				row += "\multirow{" + str(len(logic_sizes)) + "}{*}{" + system_size + "}"
			else:
				row += " "
				
			row += " & " + logic_size 
			for nudge_size in nudge_size:
				# get the p-value
				query = "logic_size == " + logic_size + " | system_size == " + system_size + " | nudge_size == " + nudge_size
				p_value = list(df.query(query).p_value.unique())[0]
				
				# truncate
				p_value = p_value[0:6]
				
				# add column
				row += " & " + p_value
				
				# add stars based on significance
				if float(p_value) < 0.001:
					row += "***"
				elif float(p_value) < 0.005:
					row += "**"
				elif float(p_value) < 0.05:
					row += "*"
				
			outfile.write(row + "\\")
			outfile.write("\hline")
		
	outfile.write("\end(tabular)")
	outfile.write("\caption(%s)" % caption)
	outfile.write("\end(figure)")
	outfile.close()
	return 0

def main():
	df = pd.DataFrame(columns=["experiment", "system_size", "logic_size", "nudge_size", "p-value"])

	loc_counter = 0
	samplesize = "0"
	
	for root, dirs, files in os.walk("../../data"):
		for file in files:
			if file.endswith('.txt'):
				# this is an experiment, to be added to our tables
				cgs = re.findall('n=([0-9]+)_l=([0-9]+)_e=([0-9]+.[0-9]+)', root)
				experiment = re.findall('more_([a-z _]+)', file)
				with f as open(os.path.join(root, file), 'r'):
					result = f.readline()
					p_value = re.findall('p=([0-9]+.[0-9]+)', result)
					
					df.loc[loc_counter] = experiment + cgs + p_value
					loc_counter += 1
				
				samplesize = re.findall('nosamples=([0-9]+)', root)[0]
					
	
	# create the 4 tables
	for experiment in list(df.experiment.unique()):
		caption = "(n=" + samplesize + ")"
		create_table(df, experiment, samplesize)
	
if __name__ == '__main__':
	main()