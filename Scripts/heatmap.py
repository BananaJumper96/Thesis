import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage


# Clustering by genus
def parse_args():
    parser = argparse.ArgumentParser(description='Take all_bin_stats_coverage.tsv file and analyze features.')
    parser.add_argument('--input_file', required=True, help='Path to the all_bin_stats_coverage.tsv file')
    parser.add_argument('--output', required=True, help='Path to the saved files')
    return parser.parse_args()
# Load dataset
def load_data(file_path):
    data = pd.read_csv(file_path, sep='\t')
    return data

# Extract genus and prepare data for clustering
def prepare_data(data):
    # Extract genus from the classification column
    data['Genus'] = data['classification'].str.extract(r'g__([^;]+)')
    # Group by genus and compute the mean for numerical attributes
    grouped_by_genus = data.groupby('Genus').mean(numeric_only=True)
    # Normalize the data for heatmap visualization
    normalized_data = (grouped_by_genus - grouped_by_genus.min()) / (grouped_by_genus.max() - grouped_by_genus.min())
    # Perform hierarchical clustering on the normalized data
    linkage_matrix = linkage(normalized_data, method="ward")
    return normalized_data, linkage_matrix



# Main function to execute the workflow
def main():
    args = parse_args()
    data = load_data(args.input_file)
    normalized_data, linkage_matrix = prepare_data(data)
    # Plot the heatmap
    sns.clustermap(
    normalized_data,
    cmap="coolwarm",
    linewidths=1,
    linecolor="white",
    row_linkage=linkage_matrix,
    col_cluster=False,  # Only cluster rows (genera)
    )
    plt.title("Genus Clusters", loc='center')
    plt.xlabel("Attributes")
    plt.ylabel("Genus")
    plt.show()
    try:
        plt.savefig(args.output, format='pdf', bbox_inches='tight')
        print(f"Plot saved to {args.output}")
    except Exception as e:
        print(f"Error saving plot: {e}")
    

if __name__ == "__main__":
    main()