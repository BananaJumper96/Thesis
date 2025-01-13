import os
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import argparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans, DBSCAN as dbscan
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
import re
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description='Generate visualizations from GEMs stats.')
    parser.add_argument('--bins',nargs='+', required=True, help='Path to the directory containing bin files')
    parser.add_argument('--bin_stats', nargs='+', required=True, help='Paths to the all_bin_stats_coverage.tsv files')
    parser.add_argument('--GEM_stats', required=True, help='Path to the GEMs.stats file')
    parser.add_argument('--output', required=True, help='Path to the output directory')
    return parser.parse_args()
# Function to count proteins in .faa files
def count_unique_proteins_in_file(file_path):
    """
    Count the number of unique protein entries in a .faa file.
    Each protein entry starts with a '>' in the file, and unique entries are identified by their IDs.

    Args:
        file_path (str): Path to the .faa file.

    Returns:
        int: Number of unique protein entries in the file.
    """
    unique_ids = set()
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    protein_id = line.split()[0]  # Extract the first part of the header as the ID
                    unique_ids.add(protein_id)
        return len(unique_ids)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        return 0

def plot_gene_protein_scatter(dataframe):
    """
    Create a scatter plot showing genes vs protein count
    """
    # Create figure
    plt.figure(figsize=(10, 10))
    ax = plt.gca()
    
    # Find min and max for axes
    max_val = max(dataframe['Genes'].max(), dataframe['ProteinCount'].max())
    min_val = min(dataframe['Genes'].min(), dataframe['ProteinCount'].min())
    
    # Create scatter plot
    ax.scatter(dataframe['ProteinCount'], 
              dataframe['Genes'],
              s=100,  # Size of points
              color='#3498db',  # Nice blue color
              edgecolor='black',
              alpha=0.6)
    
    # Add y=x line (ratio = 1)
    ax.plot([min_val, max_val], [min_val, max_val], 
            'k--', alpha=0.5, label='1:1 ratio')
    
    # Customize plot
    ax.set_xlabel('Protein Count')
    ax.set_ylabel('Gene Count')
    ax.set_title('Gene Count in GEMs vs Protein Count in bins')
    
    # Set equal aspect ratio and limits
    ax.set_aspect('equal')
    ax.set_xlim(min_val*0.9, max_val*1.1)
    ax.set_ylim(min_val*0.9, max_val*1.1)
    
    # Add grid
    ax.grid(True, linestyle='--', alpha=0.3, color='grey')
    ax.set_axisbelow(True)
    
    # Add summary statistics
    stats_text = ('Summary Statistics:\n'
                 f'Average ratio: {(dataframe["Genes"] / dataframe["ProteinCount"]).mean():.3f}\n'
                 f'Min ratio: {(dataframe["Genes"] / dataframe["ProteinCount"]).min():.3f}\n'
                 f'Max ratio: {(dataframe["Genes"] / dataframe["ProteinCount"]).max():.3f}\n'
                 f'Number of genomes: {len(dataframe)}')
    
    plt.text(1.1, 0.5, stats_text,
             transform=ax.transAxes,
             bbox=dict(facecolor='white', edgecolor='black'),
             verticalalignment='center')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    plt.savefig('results/gene_protein_scatter.png',
                dpi=300,
                bbox_inches='tight')
    plt.close()

def process_and_merge_data(args):
        """
        Process FAA files, bin stats, and GEM stats and merge them into a single dataframe.
        
        Args:
            args: Parsed command line arguments containing bins, bin_stats, GEM_stats, and output paths
        
        Returns:
            pd.DataFrame: Merged statistics dataframe
        """
        sample_faa_files = {}
        # Process FAA files
        for folder in args.bins:
            for root, dirs, files in os.walk(folder):
                sample = root.split('/')[-2]
                for file in files:
                    if file.endswith(".faa"):
                        bin_id = re.search(r'bin\.(\d+)\.faa$', file)
                        if bin_id:
                            bin_id = f"bin.{bin_id.group(1)}"
                            if sample not in sample_faa_files:
                                sample_faa_files[sample] = {}
                            sample_faa_files[sample][bin_id] = os.path.join(root, file)

        # Process bin stats files
        all_bin_stats_path = {}
        for folder in args.bin_stats:
            for root, dirs, files in os.walk(folder):
                sample = root.split('/')[-1]
                for file in files:
                    if file.endswith("all_bin_stats_coverage.tsv"):
                        if sample not in all_bin_stats_path:
                            all_bin_stats_path[sample] = os.path.join(root, file)

        # Count proteins
        protein_counts = defaultdict(lambda: defaultdict(int))
        for sample, bins in sample_faa_files.items():
            for bin_id, file_path in bins.items():
                protein_counts[sample][bin_id] = count_unique_proteins_in_file(file_path)

        # Combine bin stats
        all_bin_stats = pd.DataFrame()
        for sample, file_path in all_bin_stats_path.items():
            bin_stats = pd.read_csv(file_path, sep='\t')
            bin_stats["Samples"] = sample
            bin_stats.rename(columns={"Bin": "BinID"}, inplace=True)
            all_bin_stats = pd.concat([all_bin_stats, bin_stats])

        # Create output directory
        os.makedirs(args.output, exist_ok=True)

        # Process GEM stats
        gems_stats = pd.read_csv(args.GEM_stats, sep=r'\s+')
        gems_stats.columns = ["Samples", "BinID", "Metabolites", "Reactions", "Genes"]
        gems_stats["ProteinCount"] = gems_stats.apply(
            lambda row: protein_counts[row["Samples"]].get(row["BinID"], None), axis=1
        )
        gems_stats["GenesToProteinCount"] = gems_stats.apply(
            lambda row: row["Genes"] / row["ProteinCount"] if row["ProteinCount"] else None, axis=1
        )

        # Merge stats
        merged_stats = pd.merge(
            gems_stats,
            all_bin_stats,
            on=['Samples', 'BinID'],
            how='left',
            validate='1:1'
        )

        # Save merged stats
        merged_stats.to_csv(os.path.join(args.output, "merged_stats.tsv"), sep='\t', index=False)
        
        return merged_stats

def create_quality_bubble_subplot(dataframe):
    """
    Create separate bubble plots for each bacterial family showing genome quality metrics
    with abundance information.
    """
    # Extract family level from classification
    dataframe['family'] = dataframe['classification'].apply(
        lambda x: x.split(';')[4].replace('f__', '')
    )

    # count cag-508 genomes
    cag508 = dataframe[dataframe['family'] == 'CAG-508']
    print(cag508.value_counts())

    # Get unique families and calculate number of rows/columns needed
    unique_families = sorted(dataframe['family'].unique())
    n_families = len(unique_families)
    n_cols = min(3, n_families)  # Maximum 3 columns
    n_rows = (n_families + n_cols - 1) // n_cols  # Ceiling division
    
    # Create figure and subplots with extra space on right for legend
    fig, axes = plt.subplots(n_rows, n_cols, 
                            figsize=(6*n_cols + 2, 5*n_rows),  # Added extra width for legend
                            squeeze=False)  # squeeze=False ensures axes is always 2D
    
    # Flatten axes for easier iteration
    axes_flat = axes.flatten()
    
    # Create color palette for all plots to use
    colors = plt.cm.Set3(np.linspace(0, 1, n_families))
    
    # Plot each family
    for idx, (family, color) in enumerate(zip(unique_families, colors)):
        ax = axes_flat[idx]
        family_data = dataframe[dataframe['family'] == family]
        
        # Create bubble plot
        size_scale = 1000  # Adjust this value to make bubbles larger or smaller
        ax.scatter(family_data['Completeness'],
                  family_data['Contamination'],
                  s=family_data['mg_rel_abundance'] * size_scale,
                  c=[color],  # Use the same color for all points in this family
                  alpha=0.6)
        
        # Customize each subplot
        ax.set_xlabel('Completeness (%)')
        ax.set_ylabel('Contamination (%)')
        ax.set_title(f'{family}')
        
        # Add number of genomes in this family
        n_genomes = len(family_data)
        ax.text(0.02, 0.98,
                f'n = {n_genomes}',
                transform=ax.transAxes,
                bbox=dict(facecolor='white', edgecolor='gray', alpha=0.8),
                va='top', ha='left',
                fontsize=9)
        
        # Set consistent axis limits across all subplots
        ax.set_xlim(dataframe['Completeness'].min() - 1, 
                   dataframe['Completeness'].max() + 1)
        ax.set_ylim(dataframe['Contamination'].min() - 0.5,
                   dataframe['Contamination'].max() + 0.5)
    
    # Remove empty subplots if any
    for idx in range(n_families, len(axes_flat)):
        fig.delaxes(axes_flat[idx])
    
    # Add size legend to the figure
    # Create custom legend for bubble sizes
    abundance_levels = [
        dataframe['mg_rel_abundance'].max(),
        dataframe['mg_rel_abundance'].max()/2,
        dataframe['mg_rel_abundance'].min()
    ]
    
    legend_elements = [plt.scatter([], [], 
                                 c='gray',
                                 s=abundance * size_scale,
                                 alpha=0.3,
                                 label=f'Abundance: {abundance:.2f}')
                      for abundance in abundance_levels]
    
    # Add legend with adjusted position and size
    fig.legend(handles=legend_elements,
              title='Relative Abundance',
              bbox_to_anchor=(1.1, 0.9),  # Moved further right
              loc='upper left',
              borderaxespad=0,
              handletextpad=6,
              labelspacing=8,  # Increased spacing between legend entries
              prop={'size': 10},  # Adjust font size if needed
              borderpad=5
              )  # Added padding around legend
    # Adjust layout with tighter margins
    plt.tight_layout()
    
    # Save the plot with extra right margin for legend
    plt.savefig('results/quality_bubble_subplots.svg', 
                dpi=300, 
                bbox_inches='tight',
                pad_inches=0.5)  # Added padding
    plt.close()

def plot_singlesample_plots(dataframe):
    # Load the dataset
    data = dataframe.copy()
    # Filter for sample M11-07
    data = data[data['Samples'] == 'M11-07']
    data['genus'] = data['classification'].str.extract(r'g__([\w\s]+)')
    data['phylum'] = data['classification'].str.extract(r'p__([\w\s]+)')
    data = data.dropna(subset=['genus', 'phylum'])
    genus_completeness = data.groupby('genus')['Completeness'].mean().sort_values(ascending=False)
    plt.figure(figsize=(12, 8))
    genus_completeness.plot(kind='bar', color='skyblue', edgecolor='black', alpha=0.7)
    plt.xlabel('Genus', fontsize=12)
    plt.ylabel('Average Completeness', fontsize=12)
    plt.title('Average Completeness by Genus', fontsize=14)
    plt.xticks(rotation=90, fontsize=10)
    plt.tight_layout()
    plt.grid(axis='y', alpha=0.3)
    plt.savefig('average_completeness_by_genus.png')  # Save the plot as an image.
    plt.show()
    # Plot 2: Average Completeness by Phylum
    phylum_completeness = data.groupby('phylum')['Completeness'].mean().sort_values(ascending=False)
    plt.figure(figsize=(12, 8))
    phylum_completeness.plot(kind='bar', color='lightgreen', edgecolor='black', alpha=0.7)
    plt.xlabel('Phylum', fontsize=12)
    plt.ylabel('Average Completeness', fontsize=12)
    plt.title('Average Completeness by Phylum', fontsize=14)
    plt.xticks(rotation=90, fontsize=10)
    plt.tight_layout()
    plt.grid(axis='y', alpha=0.3)
    plt.savefig('average_completeness_by_phylum.png')  # Save the plot as an image.
    plt.show()

    # Plot 1: Average mg_rel_abundance by Genus
    if 'mg_rel_abundance' in data.columns:
        genus_abundance = data.groupby('genus')['mg_rel_abundance'].mean().sort_values(ascending=False)
        plt.figure(figsize=(12, 8))
        genus_abundance.plot(kind='bar', color='skyblue', edgecolor='black', alpha=0.7)
        plt.xlabel('Genus', fontsize=12)
        plt.ylabel('Average MG Relative Abundance', fontsize=12)
        plt.title('Average MG Relative Abundance by Genus', fontsize=14)
        plt.xticks(rotation=90, fontsize=10)
        plt.tight_layout()
        plt.grid(axis='y', alpha=0.3)
        plt.savefig('average_mg_rel_abundance_by_genus.png')  # Save the plot as an image.
        plt.show()
        # Save the plot as an image
        
        

    # Plot 2: Average mg_rel_abundance by Phylum
    if 'mg_rel_abundance' in data.columns:
        phylum_abundance = data.groupby('phylum')['mg_rel_abundance'].mean().sort_values(ascending=False)
        plt.figure(figsize=(12, 8))
        phylum_abundance.plot(kind='bar', color='lightgreen', edgecolor='black', alpha=0.7)
        plt.xlabel('Phylum', fontsize=12)
        plt.ylabel('Average MG Relative Abundance', fontsize=12)
        plt.title('Average MG Relative Abundance by Phylum', fontsize=14)
        plt.xticks(rotation=90, fontsize=10)
        plt.tight_layout()
        plt.grid(axis='y', alpha=0.3)
        # Save the plot as an image
        plt.savefig('average_mg_rel_abundance_by_phylum.png')
        plt.show()
        
    # Plot 3: Average mt_rel_expression by Genus
    if 'mt_rel_expression' in data.columns:
        genus_expression = data.groupby('genus')['mt_rel_expression'].mean().sort_values(ascending=False)
        plt.figure(figsize=(12, 8))
        genus_expression.plot(kind='bar', color='skyblue', edgecolor='black', alpha=0.7)
        plt.xlabel('Genus', fontsize=12)
        plt.ylabel('Average MT Relative Expression', fontsize=12)
        plt.title('Average MT Relative Expression by Genus', fontsize=14)
        plt.xticks(rotation=90, fontsize=10)
        plt.tight_layout()
        plt.grid(axis='y', alpha=0.3)
        # Save the plot as an image
        plt.savefig('average_mt_rel_expression_by_genus.png')
        plt.show()

    # Plot 4: Average mt_rel_expression by Phylum
    if 'mt_rel_expression' in data.columns:
        phylum_expression = data.groupby('phylum')['mt_rel_expression'].mean().sort_values(ascending=False)
        plt.figure(figsize=(12, 8))
        phylum_expression.plot(kind='bar', color='lightgreen', edgecolor='black', alpha=0.7)
        plt.xlabel('Phylum', fontsize=12)
        plt.ylabel('Average MT Relative Expression', fontsize=12)
        plt.title('Average MT Relative Expression by Phylum', fontsize=14)
        plt.xticks(rotation=90, fontsize=10)
        plt.tight_layout()
        plt.grid(axis='y', alpha=0.3)
        # Save the plot as an image
        plt.savefig('average_mt_rel_expression_by_phylum.png')
        plt.show()

def plot_heatmap(dataframe):
    #Prepare data for single sample
    dataframe = dataframe[dataframe['Samples'] == 'M11-07']
    dataframe['Genus'] = dataframe['classification'].str.extract(r'g__([^;]+)')
    features = ['mg_rel_abundance', 'mt_rel_expression', 'Completeness', 'Contamination']
    genus_pivot = dataframe.groupby('Genus')[features].mean()
    # Scale the data
    scaler = StandardScaler()
    scaled_genus_data = scaler.fit_transform(genus_pivot)
    # Create the clustermap
    heatmap = sns.clustermap(
        scaled_genus_data,
        cmap="coolwarm",
        figsize=(12, 10),
        method='ward',
        metric='euclidean',
        row_cluster=True,
        col_cluster=False
    )
    # Reposition the color bar to the left
    heatmap.cax.set_position([0.05, 0.2, 0.02, 0.6])  # [x, y, width, height]
    # Center the title
    heatmap.ax_heatmap.set_title('Heatmap with Clustering by Genus', loc='center', y=1.05)
    # save the plot
    plt.savefig('heatmap_single_sample.png', dpi=300, bbox_inches='tight')
    #Data for all samples
    dataframe = dataframe
    dataframe['Genus'] = dataframe['classification'].str.extract(r'g__([^;]+)')
    features = ['mg_rel_abundance', 'mt_rel_expression', 'Completeness', 'Contamination']
    genus_pivot = dataframe.groupby('Genus')[features].mean()
    # Scale the data
    scaler = StandardScaler()
    scaled_genus_data = scaler.fit_transform(genus_pivot)
    # Create the clustermap
    heatmap = sns.clustermap(
        scaled_genus_data,
        cmap="coolwarm",
        figsize=(12, 10),
        method='ward',
        metric='euclidean',
        row_cluster=True,
        col_cluster=False
    )
    # Reposition the color bar to the left
    heatmap.cax.set_position([0.05, 0.2, 0.02, 0.6])  # [x, y, width, height]
    # Center the title
    heatmap.ax_heatmap.set_title('Heatmap with Clustering by Genus', loc='center', y=1.05)
    # save the plot
    plt.savefig('heatmap_single_sample.png', dpi=300, bbox_inches='tight')

def completeness_to_contamination(dataframe):
    """
    Plot the completeness vs contamination curve for a dataframe with reference lines
    and percentage information
    """
    # Create a color map for unique samples
    unique_samples = dataframe['Samples'].unique()
    color_dict = dict(zip(unique_samples, range(len(unique_samples))))
    numerical_colors = [color_dict[sample] for sample in dataframe['Samples']]
    
    # Calculate percentage of high-quality points
    high_quality = dataframe[
        (dataframe['Completeness'] >= 96) & 
        (dataframe['Contamination'] <= 2.5)
    ]
    percentage = (len(high_quality) / len(dataframe)) * 100

    print(f'High-quality genomes (≥96% complete, ≤2.5% contaminated): {percentage:.1f}%')
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    
    # Add reference lines with ticks (dashed style)
    plt.axvline(x=96, color='gray', linestyle='--', alpha=0.5)
    plt.axhline(y=2.5, color='gray', linestyle='--', alpha=0.5)
    
    # Create the scatter plot
    scatter = plt.scatter(dataframe['Completeness'], 
                         dataframe['Contamination'], 
                         c=numerical_colors,
                         cmap='viridis')
    
    # Customize the plot
    plt.xlabel('Completeness (%)')
    plt.ylabel('Contamination (%)')
    plt.title('Completeness vs Contamination by Sample')
    
    # Add legend with sample names
    legend_elements = [plt.scatter([], [], c=plt.cm.viridis(i/(len(unique_samples)-1)), 
                                 label=sample) 
                      for i, sample in enumerate(unique_samples)]
    plt.legend(handles=legend_elements)
    
    # Add percentage information in a box
    # plt.text(0.02, 0.98, f'High-quality genomes\n(≥96% complete, ≤2.5% contaminated):\n{percentage:.1f}%', 
    #          transform=plt.gca().transAxes,  # Use relative coordinates
    #          bbox=dict(facecolor='white', edgecolor='gray', alpha=0.8),
    #          va='top', ha='left',
    #          fontsize=9)
    
    # Save the plot
    plt.savefig('results/completeness_vs_contamination.png', dpi=300, bbox_inches='tight')
    plt.close()
def boxplots_allsamples(dataframe):
    # Extracting necessary data for plotting
    data = dataframe[['Samples', 'Completeness']]

    # Creating the boxplot
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=data, x='Samples', y='Completeness', palette='Set2')
    plt.title('Completeness Across Samples', fontsize=14)
    plt.xlabel('Samples', fontsize=12)
    plt.ylabel('Completeness (%)', fontsize=12)
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Save the plot as an image
    plt.savefig('completeness_boxplot.png')
def scatter_violin_plot(dataframe):
    # Define bar_width for overlapping bars
    data = dataframe.copy()
    # # Prepare data for overlapping bar plot by genus
    x_genus = np.arange(len(data))
    # plt.figure(figsize=(14, 6))
    # # Base ProteinCount bar (light blue)
    # plt.bar(
    #     x_genus, 
    #     grouped_data['ProteinCount'], 
    #     width=0.6, 
    #     color='lightblue', 
    #     alpha=1.0, 
    #     label='Protein Count'
    # )
    # # Calculated_ProteinCount bar stacked on top (light red, more transparent)
    # plt.bar(
    #     x_genus, 
    #     grouped_data['Calculated_ProteinCount'], 
    #     width=0.6, 
    #     color='red', 
    #     alpha=0.3, 
    #     bottom=grouped_data['ProteinCount'],
    #     label='Calculated Protein Count'
    # )
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(os.path.join(plots_dir, 'genus_protein_counts.png'))
    # plt.close() # Close the plot to free memory

    # Scatter plot for Protein Count
    plt.figure(figsize=(10, 6))
    plt.scatter(data['Completeness'], data['ProteinCount'], color='lightblue', alpha=0.7, label='Protein Count')
    plt.scatter(data['Completeness'], data['Genes'], color='red', alpha=0.7, label='GEM Gene Count')
    plt.title('Protein Count and GEM Gene Count vs Completeness for all Samples')
    plt.xlabel('Completeness (%)')
    plt.ylabel('Counts')
    plt.grid(True)
    plt.legend()
    plt.savefig(os.path.join(plots_dir, 'completeness_scatter2.png'))
    plt.close() # Close the plot to free memory

    # Set up bins and labels for the Completeness_Range column
    bins = [0, 70, 85, 90, 100]
    labels = ['0-70%', '70-85%', '85-90%', '90-100%']
    data['Completeness_Range'] = pd.cut(data['Completeness'], bins=bins, labels=labels, include_lowest=True)
    
    plt.figure(figsize=(14, 6))
    
    # Violin plot for Protein Count
    plt.subplot(1, 2, 1)
    sns.violinplot(data=data, x='Completeness_Range', y='ProteinCount', palette='Blues')
    plt.title('Protein Count by Completeness Range')
    plt.xlabel('Completeness Range')
    plt.ylabel('Protein Count')
    plt.xticks(rotation=45)
    
    # Violin plot for GEM Gene Count
    plt.subplot(1, 2, 2)
    sns.violinplot(data=data, x='Completeness_Range', y='Genes', palette='Reds')
    plt.title('GEM Gene Count by Completeness Range')
    plt.xlabel('Completeness Range')
    plt.ylabel('GEM Gene Count')
    plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'completeness_violin.png'))
    plt.close()

def plot_gene_protein_scatter_dotted(dataframe):
    """
    Create a scatter plot showing genes vs protein count
    """
    # Create figure
    plt.figure(figsize=(10, 10))
    ax = plt.gca()
    
    # Find min and max for axes
    max_val = max(dataframe['Genes'].max(), dataframe['ProteinCount'].max())
    min_val = min(dataframe['Genes'].min(), dataframe['ProteinCount'].min())
    
    # Create scatter plot
    ax.scatter(dataframe['ProteinCount'], 
              dataframe['Genes'],
              s=100,  # Size of points
              color='#3498db',  # Nice blue color
              edgecolor='black',
              alpha=0.6)
    
    # Add y=x line (ratio = 1)
    ax.plot([min_val, max_val], [min_val, max_val], 
            'k--', alpha=0.5, label='1:1 ratio')
    
    # Customize plot
    ax.set_xlabel('Protein Count')
    ax.set_ylabel('Gene Count')
    ax.set_title('Gene Count in GEMs vs Protein Count in bins')
    
    # Set equal aspect ratio and limits
    ax.set_aspect('equal')
    ax.set_xlim(min_val*0.9, max_val*1.1)
    ax.set_ylim(min_val*0.9, max_val*1.1)
    
    # Add grid
    ax.grid(True, linestyle='--', alpha=0.3, color='grey')
    ax.set_axisbelow(True)
    
    # Add summary statistics
    stats_text = ('Summary Statistics:\n'
                 f'Average ratio: {(dataframe["Genes"] / dataframe["ProteinCount"]).mean():.3f}\n'
                 f'Min ratio: {(dataframe["Genes"] / dataframe["ProteinCount"]).min():.3f}\n'
                 f'Max ratio: {(dataframe["Genes"] / dataframe["ProteinCount"]).max():.3f}\n'
                 f'Number of genomes: {len(dataframe)}')
    
    plt.text(1.1, 0.5, stats_text,
             transform=ax.transAxes,
             bbox=dict(facecolor='white', edgecolor='black'),
             verticalalignment='center')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    plt.savefig('results/gene_protein_scatter.png',
                dpi=300,
                bbox_inches='tight')
    plt.close()
def main():
    # load the args
    args = parse_args()
    # Process and merge data
    merged_stats = process_and_merge_data(args)
    # Create explicit copies to avoid SettingWithCopyWarning
    gems_stats = merged_stats.copy()
    output = args.output
    plots_dir = os.path.join(output, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    # Move into the plots directory
    os.chdir(plots_dir)
    # Remove rows with NaN values
    gems_stats = gems_stats.dropna() # Drop rows with NaN values
    # Visualizations ------------------------------------------------------------------------------------------------------------------------------------ 
    data = gems_stats.copy() # Use the merged dataframe for visualizations
    # Create output directory for plots
    plot_singlesample_plots(gems_stats)
    plot_heatmap(gems_stats)
    create_quality_bubble_subplot(gems_stats)
    boxplots_allsamples(gems_stats)
    completeness_to_contamination(gems_stats)
    scatter_violin_plot(gems_stats)
    # Extract genus information from the 'classification' column
    if 'classification' in data.columns: # Check if the column exists
        data.loc[:, 'Genus'] = data['classification'].str.split(';g__').str[1].str.split(';').str[0] # Extract the genus
    else:
        print("Warning: 'classification' column not found") 
        data.loc[:, 'Genus'] = 'Unknown'

    # Plot the 'GenesToProteinCount' ratio with updated label names
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=data, x='Samples', y='GenesToProteinCount', palette="viridis")
    plt.title('Distribution of Genes Divided by Bin Protein Ratio Across Samples')
    plt.xlabel('Samples')
    plt.ylabel('Genes Divided by Bin Protein Ratio')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
    # Use `ProteinCount` instead of `GeneToProteinRatio` where applicable
    data.loc[:, 'Calculated_ProteinCount'] = data['Genes'] / data['GenesToProteinCount'] 
    data['Calculated_ProteinCount'] = data['Genes'] / data['GenesToProteinCount']
    # print(data)
    # Group by genus and sum the counts for each genus
    grouped_data = data.groupby('Genus')[['ProteinCount', 'Calculated_ProteinCount']].sum().reset_index()
    # print(grouped_data)
    # Sort the data by total protein count (ProteinCount + Calculated_ProteinCount)
    grouped_data['TotalProteinCount'] = grouped_data['ProteinCount'] + grouped_data['Calculated_ProteinCount']
    grouped_data = grouped_data.sort_values(by='TotalProteinCount', ascending=False)
    

    # # Clustering Analysis
    # numerical_features = [
    #     'Metabolites', 'Reactions', 'Genes', 'ProteinCount',
    #     'mg_rel_abundance', 'mt_rel_expression',
    #     'Completeness', 'Contamination'
    # ]
    # scaler = StandardScaler()
    # scaled_data = scaler.fit_transform(data[numerical_features])

    
    # plt.title('Clustered Heatmap of Pairwise Correlation of All Enzymes')
    # plt.show()
    # plt.savefig(os.path.join(plots_dir, 'enzyme_correlation_heatmap.png'))
    # plt.close()
    # normalized_heatmap_data = minmax_scaler.fit_transform(heatmap_data)
    # normalized_heatmap_data_df = pd.DataFrame(normalized_heatmap_data, index=heatmap_data.index, columns=heatmap_data.columns)

    # Generate a clustered heatmap with normalized data
    sns.clustermap(
        normalized_heatmap_data_df,
        cmap='coolwarm',
        figsize=(12, 10),
        method='ward',
        cbar_kws={'label': 'Normalized Value (0-1)'}
    )
    plt.title('Clustered Heatmap of Genera (Normalized Data)')
    plt.xlabel('Attributes')
    plt.ylabel('Genus Clusters')
    plt.show()

if __name__ == "__main__":
    main()