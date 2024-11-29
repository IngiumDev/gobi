#!/usr/bin/python3
import matplotlib.pyplot as plt
import os
import pandas as pd

RESULTS = ""
'''The report should contain
at least two cumulative plots (see definition on the internet or an R tutorial) showing the
1Guideline: a bunch of plots without a clear explanation of what is shown, without description of what
can be observed, or without a statement what can be concluded is not a sufficient report and will be graded
accordingly.
3
distributions of the maximum number of skipped exons and skipped bases, respectively,
per ES-SE for every GTF file. The plots have to be saved also into your output directory
named skipped_exons.jpg and skipped_bases.jpg. '''

if __name__ == '__main__':
    # List to store the file names
    tsv_data = {}

    # Iterate over all files in the directory
    for file_name in os.listdir(RESULTS):
        # only if .tsv
        if file_name[-4:] == ".tsv":
            data = pd.read_csv(file_name, sep="\t")
            tsv_data[file_name] = data
    
    # Create the figures
    fig, ax = plt.subplots()
    
    
        
