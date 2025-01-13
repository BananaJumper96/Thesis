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
def main():
    # load the args
    args = parse_args()
    # Define the samples folder path
    sample_faa_folder = args.bins
    # Define the bin_stats folder path
    all_bin_stats_folder = args.bin_stats
    # Define the Gem stats path
    gems_stats_path = args.GEM_stats
    # Define the output directory
    output = args.output
    # print("This is the output for:sample_faa_folder", sample_faa_folder)
    # print("This is the output for:all_bin_stats_folder", all_bin_stats_folder)
    # print("This is the output for:GEM_stats", gems_stats_path)
    # Get all .faa files in the sample_faa_folder
    sample_faa_files = {}
    for folder in sample_faa_folder:
        for root, dirs, files in os.walk(folder):
            sample = root.split('/')[-2] # Get the sample name from path
            for file in files:
                if file.endswith(".faa"):
                    bin_id = re.search(r'bin\.(\d+)\.faa$', file) # Extract bin ID from the file name
                    if bin_id:
                        bin_id = f"bin.{bin_id.group(1)}" # Format the bin ID
                        if sample not in sample_faa_files: # Check if the sample is already in the dictionary
                            sample_faa_files[sample] = {} # Initialize the dictionary for the sample
                        sample_faa_files[sample][bin_id] = os.path.join(root, file)
    # print("this is the output for:sample_faa_files", sample_faa_files)
    # Get all all_bin_stats_coverage.tsv files in the all_bin_stats_folder
    all_bin_stats_path = {}
    for folder in all_bin_stats_folder:
        for root, dirs, files in os.walk(folder): # Walk through the directory
            sample = root.split('/')[-1]  # Get the sample name from path
            for file in files:
                if file.endswith("all_bin_stats_coverage.tsv"):
                    if sample not in all_bin_stats_path: # Initialize the dictionary for the sample
                        all_bin_stats_path[sample] = os.path.join(root, file) # Store the file path
    # print("This is the output for:all_bin_stats_path", all_bin_stats_path)
    # Iterate over the sample_faa_files dictionary 
    protein_counts = defaultdict(lambda: defaultdict(int))
    # Iterate over the list of file paths
    for sample, bins in sample_faa_files.items():
        for bin_id, file_path in bins.items(): # Iterate over the bin IDs and file paths
            protein_count = count_unique_proteins_in_file(file_path) # Count the proteins in the file
            protein_counts[sample][bin_id] = protein_count # Store the protein count in the dictionary
    # print("This is the output for:protein_counts", protein_counts)
    # Define the paths to the .faa files, all_bin_stats_coverage.tsv, and GEMs.stats files
    # Initialize an empty DataFrame to store all bin stats
    all_bin_stats = pd.DataFrame()
    # Read and combine all bin stats files
    for sample, file_path in all_bin_stats_path.items():
        bin_stats = pd.read_csv(file_path, sep='\t')
        print("This is the output for:bin_stats", bin_stats)
        # Add the sample name to the DataFrame
        bin_stats["Samples"] = sample
        # Rename the columns to match the GEMs.stats file
        bin_stats.rename(columns={"Bin": "BinID"}, inplace=True)
        all_bin_stats = pd.concat([all_bin_stats, bin_stats])
    # print("This is the output for:all_bin_stats", all_bin_stats)
    # Check if there are NaN values in the DataFrame
    # Configure pandas to display all rows and columns
    pd.set_option('display.max_rows', None)  # Show all rows
    pd.set_option('display.max_columns', None)  # Show all columns
    pd.set_option('display.width', None)  # No truncation in width
    # print("Rows with NaN values:")
    # print(all_bin_stats[all_bin_stats.isna().any(axis=1)])
    # Ensure output directory exists
    os.makedirs(output, exist_ok=True)
    # Read GEMs.stats and extract necessary columns
    gems_stats = pd.read_csv(gems_stats_path, sep=r'\s+')
    gems_stats.columns = ["Samples", "BinID", "Metabolites", "Reactions", "Genes"]
    # Merge protein counts with GEMs.stats
    # Add the ProteinCount column
    gems_stats["ProteinCount"] = gems_stats.apply(lambda row: protein_counts[row["Samples"]].get(row["BinID"], None), axis=1
    )
    # Add the GenesToProteinCount column
    gems_stats["GenesToProteinCount"] = gems_stats.apply(
    lambda row: row["Genes"] / row["ProteinCount"] if row["ProteinCount"] else None, axis=1
    )
    # Check columns in all_bin_stats
    print("Columns in all_bin_stats:", all_bin_stats.columns)
    # Check columns in gems_stats
    print("Columns in gems_stats:", gems_stats.columns)

    # Ensure the merge columns are present in both DataFrames
    # Merge the DataFrames using both Samples and BinID as keys
    merged_stats = pd.merge(
        gems_stats,
        all_bin_stats,
        on=['Samples', 'BinID'],
        how='left',
        validate='1:1'  # Ensures one-to-one merge
    )

    # Save the merged dataframe
    merged_stats.to_csv(os.path.join(output, "merged_stats.tsv"), sep='\t', index=False)
    
    # Print info about the merge
    # print(f"Original GEMs stats rows: {len(gems_stats)}")
    # print(f"Original bin stats rows: {len(all_bin_stats)}")
    # print(f"Merged stats rows: {len(merged_stats)}")
    # print(merged_stats)
    
    
    # Update gems_stats with merged data
    # Create explicit copies to avoid SettingWithCopyWarning
    gems_stats = merged_stats.copy()

    # Remove rows with NaN values
    gems_stats = gems_stats.dropna() # Drop rows with NaN values
    # Visualizations ------------------------------------------------------------------------------------------------------------------------------------ 
    data = gems_stats.copy() # Use the merged dataframe for visualizations
    # Create output directory for plots
    plots_dir = os.path.join(output, 'plots')
    os.makedirs(plots_dir, exist_ok=True)
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

    # print(grouped_data)
    # Define bar_width for overlapping bars
    bar_width = 0.6

    # Prepare data for overlapping bar plot by genus
    x_genus = np.arange(len(grouped_data))
    plt.figure(figsize=(14, 6))
    # Base ProteinCount bar (light blue)
    plt.bar(
        x_genus, 
        grouped_data['ProteinCount'], 
        width=0.6, 
        color='lightblue', 
        alpha=1.0, 
        label='Protein Count'
    )
    # Calculated_ProteinCount bar stacked on top (light red, more transparent)
    plt.bar(
        x_genus, 
        grouped_data['Calculated_ProteinCount'], 
        width=0.6, 
        color='red', 
        alpha=0.3, 
        bottom=grouped_data['ProteinCount'],
        label='Calculated Protein Count'
    )
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'genus_protein_counts.png'))
    plt.close() # Close the plot to free memory

    # Scatter Plot
    plt.xlabel('Genus')
    plt.ylabel('Protein Count')
    plt.xticks(x_genus + bar_width / 2, grouped_data['Genus'], rotation=90, fontsize=8) 
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(os.path.join(plots_dir, 'completeness_scatter.png'))
    plt.close() # Close the plot to free memory

    # Viol
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

    # Violin Plot
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

    # Clustering Analysis
    numerical_features = [
        'Metabolites', 'Reactions', 'Genes', 'ProteinCount',
        'mg_rel_abundance', 'mt_rel_expression',
        'Completeness', 'Contamination'
    ]
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data[numerical_features])

    # Perform hierarchical clustering on the pairwise enzyme correlations
    linkage_matrix = linkage(enzyme_correlation, method='ward')

    # Plot the clustered heatmap
    plt.figure(figsize=(14, 12))
    sns.clustermap(
        enzyme_correlation,
        cmap="viridis",
        figsize=(14, 12),
        row_linkage=linkage_matrix,
        col_linkage=linkage_matrix,
        cbar_kws={'label': 'Correlation'},
        xticklabels=False,  # Hide labels for clarity in dense plots
        yticklabels=False,
    )
    plt.title('Clustered Heatmap of Pairwise Correlation of All Enzymes')
    plt.show()
    plt.savefig(os.path.join(plots_dir, 'enzyme_correlation_heatmap.png'))
    plt.close()
    normalized_heatmap_data = minmax_scaler.fit_transform(heatmap_data)
    normalized_heatmap_data_df = pd.DataFrame(normalized_heatmap_data, index=heatmap_data.index, columns=heatmap_data.columns)

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