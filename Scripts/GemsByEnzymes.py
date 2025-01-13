import pandas as pd
import numpy as np
import seaborn as sns
import os
import matplotlib.pyplot as plt
import argparse
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

def parse_args():
    parser = argparse.ArgumentParser(description='Create a heatmap of the enzymes per gems')
    parser.add_argument('--ec_summary', required=True, help='Path to the EC summary file')
    parser.add_argument('--ec_dir', required=True, help='Path to the directory containing the EC files')
    parser.add_argument('--output', required=True, help='Path to the output directory')
    return parser.parse_args()

def main():
    args = parse_args()
    ec_summary = pd.read_csv(args.ec_summary, sep='\\s+', header=None, names=['Count', 'EC_Number'])
    # Get all unique EC numbers
    ec_numbers = ec_summary['EC_Number']
    # Get the directory containing the EC files
    ec_dir = args.ec_dir
    # Get a list of .xml.ec files in the directory
    ec_files = [f for f in os.listdir(ec_dir) if f.endswith('.xml.ec')]
    # Trim the first column of the ec_files
    # Initialize a dataframe to store the presence of each EC number in each file
    ec_presence = pd.DataFrame(0, index=ec_files, columns=ec_numbers)

    # Iterate through each .xml.ec file
    for ec_file in ec_files:
        # Read the file
        with open(os.path.join(ec_dir, ec_file), 'r') as file:
            content = file.read()
            # Check for the presence of each EC number in the file
            for ec_number in ec_numbers:
                if ec_number in content:
                    ec_presence.at[ec_file, ec_number] = 1
            ec_vector = ec_presence.loc[ec_file]
    # Perform clustering on EC presence
    clustering = AgglomerativeClustering(n_clusters=6, affinity='euclidean', linkage='ward')
    clusters = clustering.fit_predict(ec_presence)
    # Add cluster labels to the DataFrame
    ec_presence['Cluster'] = clusters
    
    # Calculate enzyme correlations
    enzyme_correlation = ec_presence.drop(columns=['Cluster']).corr()

    # Visualize the data as a heatmap
    # Visualize the data as a heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(ec_presence.drop(columns=['Cluster']), cmap='viridis', cbar=True)
    plt.title('Heatmap of EC Presence/Absence')
    plt.xlabel('Enzyme EC Numbers')
    plt.ylabel('Bins (Files)')
    plt.tight_layout()
    # Save the heatmap
    plt.savefig(os.path.join(args.output, 'ec_heatmap.pdf'))
    plt.show()
    
    # Perform hierarchical clustering on the pairwise enzyme correlations
    
    distance_matrix = squareform(1 - enzyme_correlation)
    linkage_matrix = linkage(distance_matrix, method='ward')

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
    # Save the clustered heatmap
    plt.savefig(os.path.join(args.output, 'enzyme_correlation_heatmap.pdf'))
    # Save the resulting EC presence matrix with clusters
    ec_presence.to_csv(os.path.join(args.output, 'ec_presence_with_clusters.csv'), index=True)

if __name__ == '__main__':
    main()