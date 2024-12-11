#!/usr/bin/python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.interpolate import interp1d
import re


def main(annotation_name, annotation_file, lengths_file, output_path):
    lengths_df = pd.read_csv(lengths_file, sep='\t')
    length_dict = lengths_df.set_index(lengths_df.columns[0]).T.to_dict()

    genic_count = 0
    split_inconsistent_count = 0
    intergenic_count = 0
    transcriptomic_count = 0
    merged_transcriptomic_count = 0
    intronic_count = 0
    gene_dict = {}
    pcr_index_count = {}
    total_annotated_count = 0
    gdists = []

    process_annotation(annotation_file, gdists, gene_dict, genic_count, intergenic_count, intronic_count,
                       merged_transcriptomic_count, pcr_index_count, split_inconsistent_count, total_annotated_count,
                       transcriptomic_count)
    rpkm_df = calculate_all_rpkm(gene_dict, genic_count, length_dict)
    plot_rpkm_histogram(annotation_name, output_path, rpkm_df)
    plot_log_rpkm_histogram(annotation_name, output_path, rpkm_df)
    plot_rpkm_scatter(annotation_name, output_path, rpkm_df)
    plot_total_vs_unique_rpkm(annotation_name, output_path, rpkm_df)

    plot_log_total_vs_unique_rpkm(annotation_name, output_path, rpkm_df)
    plot_annotation_categories(annotation_name, genic_count, intergenic_count, intronic_count,
                               merged_transcriptomic_count, output_path, split_inconsistent_count,
                               total_annotated_count, transcriptomic_count)


def plot_annotation_categories(annotation_name, genic_count, intergenic_count, intronic_count,
                               merged_transcriptomic_count, output_path, split_inconsistent_count,
                               total_annotated_count, transcriptomic_count):
    counts = {
        'All Annotated': total_annotated_count,
        'Genic': genic_count,
        'Transcriptomic': transcriptomic_count,
        'Merged Transcriptomic': merged_transcriptomic_count,
        'Intronic': intronic_count,
        'Split Inconsistent': split_inconsistent_count,
        'Intergenic': intergenic_count
    }
    # Create a DataFrame from the counts dictionary
    counts_df = pd.DataFrame(list(counts.items()), columns=['Category', 'Count'])
    # Create the bar plot
    plt.figure(figsize=(12, 6))
    bar_plot = sns.barplot(x='Category', y='Count', data=counts_df, palette='viridis')
    # Set plot title and labels
    bar_plot.set_title("Number of Reads by Annotation Category")
    bar_plot.set_xlabel('Category')
    bar_plot.set_ylabel('Count')
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45)
    # Show the plot
    plt.tight_layout()
    plt.savefig(f'{output_path}/{annotation_name}_annotation_counts.png', dpi=900)


def plot_log_total_vs_unique_rpkm(annotation_name, output_path, rpkm_df):
    # Create a scatter plot of log(total_rpkm) vs log(unique_rpkm), colored by chromosome
    plt.figure(figsize=(12, 8))
    scatter_plot = sns.scatterplot(
        data=rpkm_df,
        x=np.log1p(rpkm_df['total_rpkm']),
        y=np.log1p(rpkm_df['unique_rpkm']),
        hue='chromosome',
        palette='tab20',  # Use a colormap with many distinct colors
        alpha=0.7
    )
    # Add y=x line
    max_value = max(np.log1p(rpkm_df['total_rpkm']).max(), np.log1p(rpkm_df['unique_rpkm']).max())
    plt.plot([0, max_value], [0, max_value], color='red', linestyle='--', linewidth=2)
    # Set plot title and labels
    scatter_plot.set_title('Log(Unique RPKM) is generally lower than Log(Total RPKM)')
    scatter_plot.set_xlabel('Log(Total RPKM)')
    scatter_plot.set_ylabel('Log(Unique RPKM)')
    # Adjust legend
    scatter_plot.legend(title='Chromosome', bbox_to_anchor=(1.05, 1), loc='upper left')
    # Show the plot
    plt.tight_layout()
    plt.show()
    plt.savefig(f'{output_path}/{annotation_name}_log_total_rpkm_vs_unique_scatter.png', dpi=900)


def plot_total_vs_unique_rpkm(annotation_name, output_path, rpkm_df):
    plt.figure(figsize=(12, 8))
    scatter_plot = sns.scatterplot(
        data=rpkm_df,
        x='total_rpkm',
        y='unique_rpkm',
        hue='chromosome',
        palette='tab20',  # Use a colormap with many distinct colors
        alpha=0.7
    )
    # Add y=x line
    max_value = max(rpkm_df['total_rpkm'].max(), rpkm_df['unique_rpkm'].max())
    plt.plot([0, max_value], [0, max_value], color='red', linestyle='--', linewidth=2)
    # Set plot title and labels
    scatter_plot.set_title('Unique RPKM is generally lower than Total RPKM')
    scatter_plot.set_xlabel('Total RPKM')
    scatter_plot.set_ylabel('Unique RPKM')
    # Adjust legend
    scatter_plot.legend(title='Chromosome', bbox_to_anchor=(1.05, 1), loc='upper left')
    # Show the plot
    plt.tight_layout()
    plt.savefig(f'{output_path}/{annotation_name}_total_rpkm_vs_unique_scatter.png', dpi=900)


def plot_rpkm_scatter(annotation_name, output_path, rpkm_df):
    rpkm_df_copy = rpkm_df.sort_values(
        by='chromosome',
        key=lambda col: col.map(custom_key)
    )
    # Create a new column for the x-axis positions
    chromosomes = rpkm_df_copy['chromosome'].unique()
    x_positions = {chrom: i for i, chrom in enumerate(chromosomes)}
    rpkm_df_copy['x_pos'] = rpkm_df_copy['chromosome'].map(x_positions)
    # Add jitter to x positions to separate points within the same chromosome
    jitter_width = 0.1
    rpkm_df_copy['x_pos_jitter_total'] = rpkm_df_copy['x_pos'] + np.random.uniform(-jitter_width, jitter_width, size=len(rpkm_df_copy))
    rpkm_df_copy['x_pos_jitter_unique'] = rpkm_df_copy['x_pos'] + np.random.uniform(-jitter_width, jitter_width,
                                                                                    size=len(rpkm_df_copy))
    # Create a figure with two subplots side by side
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))
    # Plot scatter plot for total_rpkm
    axes[0].scatter(rpkm_df_copy['x_pos_jitter_total'], rpkm_df_copy['total_rpkm'], color='blue', label='Total RPKM')
    axes[0].set_title('Scatter Plot of Total RPKM (all Reads) by Chromosome')
    axes[0].set_xlabel('Chromosome')
    axes[0].set_ylabel('RPKM')
    axes[0].set_xticks(ticks=range(len(chromosomes)))
    axes[0].set_xticklabels(chromosomes)
    axes[0].legend()
    # TODO: Mark outliers
    # Plot scatter plot for unique_rpkm
    axes[1].scatter(rpkm_df_copy['x_pos_jitter_unique'], rpkm_df_copy['unique_rpkm'], color='green', label='Unique RPKM')
    axes[1].set_title('Scatter Plot of Unique RPKM (non-PCR duplicated reads) by Chromosome')
    axes[1].set_xlabel('Chromosome')
    axes[1].set_ylabel('RPKM')
    axes[1].set_xticks(ticks=range(len(chromosomes)))
    axes[1].set_xticklabels(chromosomes)
    axes[1].legend()
    # Adjust layout and show the plot
    plt.tight_layout()
    plt.savefig(f'{output_path}/{annotation_name}_rpkm_scatter.png')


def plot_log_rpkm_histogram(annotation_name, output_path, rpkm_df):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    # Plot histogram and KDE for log(total_rpkm)
    sns.histplot(np.log1p(rpkm_df['total_rpkm']), bins='auto', kde=True, ax=axes[0], color='blue', alpha=0.7)
    axes[0].set_title('Log Distribution of Total RPKM (all Reads)')
    axes[0].set_xlabel('Log(RPKM)')
    axes[0].set_ylabel('Frequency')
    # Plot histogram and KDE for log(unique_rpkm)
    sns.histplot(np.log1p(rpkm_df['unique_rpkm']), bins='auto', kde=True, ax=axes[1], color='green', alpha=0.7)
    axes[1].set_title('Log Distribution of Unique RPKM (non-PCR duplicated reads)')
    axes[1].set_xlabel('Log(RPKM)')
    axes[1].set_ylabel('Frequency')
    # Adjust layout and show the plot
    plt.tight_layout()
    plt.savefig(f'{output_path}/{annotation_name}_log_rpkm_distribution.png', dpi=900)


def plot_rpkm_histogram(annotation_name, output_path, rpkm_df):
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
    # Plot histogram and KDE for log(total_rpkm)
    sns.histplot(rpkm_df['total_rpkm'], bins='auto', kde=True, ax=axes[0], color='blue', alpha=0.7)
    axes[0].set_title('Distribution of Total RPKM (all Reads)')
    axes[0].set_xlabel('RPKM')
    axes[0].set_ylabel('Frequency')
    # Plot histogram and KDE for log(unique_rpkm)
    sns.histplot(rpkm_df['unique_rpkm'], bins='auto', kde=True, ax=axes[1], color='green', alpha=0.7)
    axes[1].set_title('Distribution of Unique RPKM (non-PCR duplicated reads)')
    axes[1].set_xlabel('Unique RPKM')
    axes[1].set_ylabel('Frequency')
    # Adjust layout and show the plot
    plt.tight_layout()
    plt.savefig(f'{output_path}/{annotation_name}_rpkm_distribution.png', dpi=900)


def calculate_all_rpkm(gene_dict, genic_count, length_dict):
    rows = []
    # Iterate through the gene_dict
    for gene_id, vals in gene_dict.items():
        length = length_dict[gene_id]["length"]
        chr = length_dict[gene_id]["chromosome"]
        total_rpkm = calc_rpkm(vals[0], genic_count, length)
        unique_rpkm = calc_rpkm(vals[1], genic_count, length)

        # Append a dictionary to the list
        rows.append({
            'chromosome': chr,
            'gene_id': gene_id,
            'total_rpkm': total_rpkm,
            'unique_rpkm': unique_rpkm
        })
    # Create the DataFrame from the list of rows
    return pd.DataFrame(rows, columns=['chromosome', 'gene_id', 'total_rpkm', 'unique_rpkm'])

def roman_to_int(roman):
    """Convert Roman numeral to an integer."""
    roman_values = {
        'I': 1, 'V': 5, 'X': 10, 'L': 50, 'C': 100,
        'D': 500, 'M': 1000
    }
    prev_value = 0
    total = 0
    for char in reversed(roman):
        value = roman_values[char]
        if value < prev_value:
            total -= value
        else:
            total += value
        prev_value = value
    return total


def is_roman(s):
    """Check if a string is a valid Roman numeral."""
    roman_regex = re.compile(r'^(?=[MDCLXVI])M{0,4}(CM|CD|D?C{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})$')
    return bool(roman_regex.match(s))


def custom_key(chromosome):
    """Custom sorting key function."""
    if chromosome.isdigit():
        # Pure numeric chromosome
        return (0, int(chromosome))
    elif is_roman(chromosome):
        # Roman numeral chromosome
        return (1, roman_to_int(chromosome))
    else:
        # Alphabetical chromosome (e.g., X, Y)
        return (2, chromosome)


def custom_sort(chromosomes):
    """Sort a list of chromosome labels using the custom key."""
    return sorted(chromosomes, key=custom_key)


def process_annotation(annotation_file, gdists, gene_dict, genic_count, intergenic_count, intronic_count,
                       merged_transcriptomic_count, pcr_index_count, split_inconsistent_count, total_annotated_count,
                       transcriptomic_count):
    with open(annotation_file, 'r') as file:
        for line in file:
            parts = line.strip().split("\t")
            if (parts[1][0] == 'm'):
                if parts[4][7] != '0':
                    genic_count += 1
                    # Now check whether it's transcriptomic, merged or intronic
                    splitmatch = parts[5].split("|")
                    for matching_gene in splitmatch:
                        matching_gene_parts = matching_gene.split(":")
                        gene_id = matching_gene_parts[0].split(",")[0]
                        pcr_index = parts[6].split(" ")[1]
                        add_to_pcr_index(pcr_index_count, pcr_index)
                        add_to_gene_dict(gene_dict, gene_id, pcr_index)
                        if matching_gene_parts[1].startswith("MERGED"):
                            merged_transcriptomic_count += 1
                        elif matching_gene_parts[1].startswith("INTRON"):
                            intronic_count += 1
                        else:
                            transcriptomic_count += 1

                else:
                    intergenic_count += 1
                    gdists.append(parts[5].split(":")[1])
            else:
                split_inconsistent_count += 1
            total_annotated_count += 1

def calc_rpkm(num_reads, total_genes_matched, gene_length):
    return (num_reads * 1_000 * 1_000_000.0) / (total_genes_matched * gene_length)
def add_to_gene_dict(gene_dict, gene_id, pcr_index):
    if gene_id in gene_dict:
        gene_dict[gene_id] = (gene_dict[gene_id][0] + 1, gene_dict[gene_id][1] + (1 if pcr_index == '0' else 0))
    else:
        gene_dict[gene_id] = (1, 1 if pcr_index == '0' else 0)


def add_to_pcr_index(pcr_index_count, pcr_index):
    if pcr_index in pcr_index_count:
        pcr_index_count[pcr_index] += 1
    else:
        pcr_index_count[pcr_index] = 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze BAM features annotation')
    parser.add_argument('--annotation', required=True, help='Name of the annotation')
    parser.add_argument('--input', required=True, help='Input Annotation file')
    parser.add_argument('--lengths', required=True, help='Lengths file')
    parser.add_argument('--output', required=True, help='Output folder')
    args = parser.parse_args()
    main(args.annotation, args.input, args.lengths, args.output)
