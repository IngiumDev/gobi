import base64  # Import the base64 module for encoding binary data
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os
import pandas as pd
import sys
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
from scipy.stats import binom


# Function to check if all regions have a length >= 5
def check_region_length(regvec):
    # Split the regvec by '|'
    regions = regvec.split('|')

    # Check each region
    for region in regions:
        start, end = region.split('-')
        start = int(start)
        end = int(end)

        # Check if the length of the region is at least 5
        if end - start < 5:
            return False

    return True


def save_plot(plots, fig, pdf, html, plot_title, ouput):
    # Save the plot to PNG file
    plot_filename = f'{plot_title.replace(" ", "_").lower()}.png'
    fig.savefig(os.path.join(ouput, plot_filename), dpi=600)

    # Add the plot to the list of plots for PDF
    pdf.savefig(fig)  # This saves the figure to the PDF

    # Convert plot to image for HTML
    img_stream = BytesIO()
    fig.savefig(img_stream, format='png')
    img_stream.seek(0)

    # Base64 encode the image data
    img_base64 = base64.b64encode(img_stream.getvalue()).decode('utf-8')

    # Generate the HTML tag to embed the image
    img_tag = f'<img src="data:image/png;base64,{img_base64}" alt="{plot_title}" style="width:100%;"/>'
    html.write(f'<h2>{plot_title}</h2>')
    html.write(img_tag)

    plots.append(fig)  # Add the plot to the list to keep track of all generated plots


def main(readmapping_file, output_path):
    # Read the mapping data
    mapping_info_df = pd.read_csv(readmapping_file, sep='\t')

    html, html_path, pdf, pdf_path = prepare_output_files(output_path)

    # Plot Fragment Length Distribution
    fig = create_fragment_length_distribution(mapping_info_df)

    save_plot([], fig, pdf, html, "Fragment Length Distribution", output_path)
    fig = create_cumulative_mutation_count_distribution(fig, mapping_info_df)
    save_plot([], fig, pdf, html, "Cumulative Distribution of Mutation Counts", output_path)

    fig = create_read_metrics_bar_plot(fig, mapping_info_df)
    save_plot([], fig, pdf, html, "Bar Plot Read Metrics", output_path)

    fig = create_mutation_position_distribution(fig, mapping_info_df)

    save_plot([], fig, pdf, html, "Histogram of Mutation Positions", output_path)

    close_output_files(html, html_path, pdf, pdf_path)


def create_mutation_position_distribution(fig, mapping_info_df):
    # Distribution of what position the mutations are in
    # Extract mutation positions for fw and rw
    fw_mutation_positions = []
    rw_mutation_positions = []
    for fw_mut in mapping_info_df['fw_mut']:
        if isinstance(fw_mut, str) and fw_mut:
            fw_mutation_positions.extend(map(int, fw_mut.split(',')))
    for rw_mut in mapping_info_df['rw_mut']:
        if isinstance(rw_mut, str) and rw_mut:
            rw_mutation_positions.extend(map(int, rw_mut.split(',')))
    # Create histograms for fw and rw mutation positions
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    # FW mutation positions histogram
    ax1.hist(fw_mutation_positions, bins=np.unique(fw_mutation_positions) + 1, color='blue', edgecolor='black',
             alpha=0.7)
    ax1.set_title('Distribution of FW Mutation Positions')
    ax1.set_xlabel('Position')
    ax1.set_ylabel('Frequency')
    ax1.grid(True)
    ax1.xaxis.set_major_locator(MultipleLocator(5))  # Set x-axis ticks to step by 5
    # RW mutation positions histogram
    ax2.hist(rw_mutation_positions, bins=np.unique(rw_mutation_positions) + 1, color='red', edgecolor='black',
             alpha=0.7)
    ax2.set_title('Distribution of RW Mutation Positions')
    ax2.set_xlabel('Position')
    ax2.set_ylabel('Frequency')
    ax2.grid(True)
    ax2.xaxis.set_major_locator(MultipleLocator(5))  # Set x-axis ticks to step by 5
    plt.tight_layout()
    return fig


def create_read_metrics_bar_plot(fig, mapping_info_df):
    all_reads = len(mapping_info_df)
    # Number of non-split reads (fw and rw)
    non_split_read_fw = len([fw for fw in mapping_info_df['fw_regvec'] if '|' not in fw])
    non_split_read_rw = len([rw for rw in mapping_info_df['rw_regvec'] if '|' not in rw])
    # Number of non-split reads with no mismatches
    non_split_read_fw_no_mismatches = len(mapping_info_df[
                                              ~mapping_info_df['fw_regvec'].str.contains(r'\|', na=False) & pd.isna(
                                                  mapping_info_df['fw_mut'])])
    non_split_read_rw_no_mismatches = len(mapping_info_df[
                                              ~mapping_info_df['rw_regvec'].str.contains(r'\|', na=False) & pd.isna(
                                                  mapping_info_df['rw_mut'])])
    # Number of split reads
    split_read_fw = len([fw for fw in mapping_info_df['fw_regvec'] if '|' in fw])
    split_read_rw = len([rw for rw in mapping_info_df['rw_regvec'] if '|' in rw])
    # Number of split reads with no mismatches
    split_read_fw_no_mismatches = len(mapping_info_df[
                                          mapping_info_df['fw_regvec'].str.contains(r'\|', na=False) & pd.isna(
                                              mapping_info_df['fw_mut'])])
    split_read_rw_no_mismatches = len(mapping_info_df[
                                          mapping_info_df['rw_regvec'].str.contains(r'\|', na=False) & pd.isna(
                                              mapping_info_df['rw_mut'])])
    # Number of split reads with no mismatches where all regions are at least 5 basepairs long
    split_read_fw_valid_length = len(mapping_info_df[
                                         mapping_info_df['fw_regvec'].str.contains(r'\|', na=False) & mapping_info_df[
                                             'fw_regvec'].apply(
                                             check_region_length) & pd.isna(mapping_info_df['fw_mut'])])
    split_read_rw_valid_length = len(mapping_info_df[
                                         mapping_info_df['rw_regvec'].str.contains(r'\|', na=False) & mapping_info_df[
                                             'rw_regvec'].apply(
                                             check_region_length) & pd.isna(mapping_info_df['rw_mut'])])
    metrics = [('All Reads', all_reads, all_reads),  # For 'All Reads', both fw and rw should be equal
               ('Non-split Reads', non_split_read_fw, non_split_read_rw),
               ('Non-split Reads (No Mismatches)', non_split_read_fw_no_mismatches, non_split_read_rw_no_mismatches),
               ('Split Reads', split_read_fw, split_read_rw),
               ('Split Reads (No Mismatches)', split_read_fw_no_mismatches, split_read_rw_no_mismatches),
               ('Split Reads (Valid Length >5)', split_read_fw_valid_length, split_read_rw_valid_length)]
    df = pd.DataFrame(metrics, columns=['Metric', 'fw', 'rw'])
    fig, ax = plt.subplots(figsize=(10, 6))
    df.set_index('Metric').plot(kind='bar', ax=ax, color=['blue', 'red'], width=0.8)
    ax.yaxis.set_major_formatter(ticker.StrMethodFormatter('{x:,.0f}'))
    ax.set_xlabel('Metrics')
    ax.set_ylabel('Count')
    ax.set_title('Barplot of Read Metrics by Forward (fw) and Reverse (rw)')
    ax.legend(title='Direction', labels=['Forward (fw)', 'Reverse (rw)'], loc='upper right')
    plt.tight_layout()
    return fig


def close_output_files(html, html_path, pdf, pdf_path):
    # Close PDF and HTML
    pdf.close()
    html.write('</body></html>\n')
    html.close()
    print(f'Plots saved to {pdf_path} and {html_path}')


def prepare_output_files(output_path):
    # Prepare the PDF and HTML output
    pdf_path = os.path.join(output_path, 'plots.pdf')
    html_path = os.path.join(output_path, 'plots.html')
    pdf = PdfPages(pdf_path)  # Create a PDF object to save the plots
    html = open(html_path, 'w')  # Open an HTML file to save the plots
    html.write('<html><body>\n')
    html.write('<h1>Simulation Plots</h1>\n')
    return html, html_path, pdf, pdf_path


def create_fragment_length_distribution(mapping_info_df):
    fragment_lengths = [int(rw.split('-')[1]) - int(fw.split('-')[0]) + 1 for rw, fw in
                        zip(mapping_info_df['t_rw_regvec'], mapping_info_df['t_fw_regvec'])]
    mean_length = np.mean(fragment_lengths)
    std_length = np.std(fragment_lengths)
    # Explicitly create the figure
    fig = plt.figure(figsize=(10, 6))
    plt.hist(fragment_lengths, bins=np.unique(fragment_lengths), color='darkblue', edgecolor="blue", alpha=1)
    plt.axvline(mean_length, color='red', linestyle='dashed', linewidth=1, label='Mean Length')
    plt.axvline(75, color='green', linestyle='dashed', linewidth=1, label='Read Length')
    plt.text(mean_length + std_length, plt.ylim()[1] * 0.9, f'Mean: {mean_length:.2f}', color='red')
    plt.text(mean_length + std_length, plt.ylim()[1] * 0.85, f'SD: {std_length:.2f}', color='red')
    plt.xlabel('Fragment Length')
    plt.ylabel('Frequency')
    plt.title('Fragment Length Distribution')
    plt.grid(True)
    return fig


def create_cumulative_mutation_count_distribution(fig, mapping_info_df):
    fw_mutation_counts = [len(fw.split(',')) if isinstance(fw, str) and fw else 0 for fw in mapping_info_df['fw_mut']]
    rw_mutation_counts = [len(rw.split(',')) if isinstance(rw, str) and rw else 0 for rw in mapping_info_df['rw_mut']]
    mutation_counts = fw_mutation_counts + rw_mutation_counts
    # Sort the mutation counts for fw and rw
    sorted_fw_mutation_counts = np.sort(fw_mutation_counts)
    sorted_rw_mutation_counts = np.sort(rw_mutation_counts)
    # Calculate cumulative counts
    cumulative_fw_mutation_counts = np.arange(1, len(sorted_fw_mutation_counts) + 1)
    cumulative_rw_mutation_counts = np.arange(1, len(sorted_rw_mutation_counts) + 1)
    # Normalize cumulative counts to get probabilities
    cumulative_fw_probabilities = cumulative_fw_mutation_counts / len(fw_mutation_counts)
    cumulative_rw_probabilities = cumulative_rw_mutation_counts / len(rw_mutation_counts)
    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(sorted_fw_mutation_counts, cumulative_fw_probabilities, label='FW Cumulative Probability', color='blue')
    ax.plot(sorted_rw_mutation_counts, cumulative_rw_probabilities, label='RW Cumulative Probability', color='red',
            linestyle='dashed')
    # Add the CDF of the binomial distribution
    n = 75
    p = 0.01
    x = np.arange(0, max(max(sorted_fw_mutation_counts), max(sorted_rw_mutation_counts)) + 1)
    binom_cdf = binom.cdf(x, n, p)
    ax.plot(x, binom_cdf, label='Binomial CDF', color='green', linestyle='dashed')
    # Add labels, title, and legend
    ax.set_title('Cumulative Distribution of Mutation Counts with Binomial(75,0.01)')
    ax.set_xlabel('Number of Mutations')
    ax.set_ylabel('Cumulative Probability')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
    ax.grid(True)
    ax.legend()
    plt.tight_layout()
    return fig


if __name__ == '__main__':
    readmapping_file = sys.argv[1]
    output_path = sys.argv[2]
    # readmapping_file = "/Users/simon/IdeaProjects/gobi/data/readsimulator/testouput/read.mappinginfo"
    # output_path = "/Users/simon/IdeaProjects/gobi/data/readsimulator/testouput"
    main(readmapping_file, output_path)
