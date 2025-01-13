import pandas as pd
import matplotlib.pyplot as plt
import os
import argparse

def arg_parser():
    parser = argparse.ArgumentParser(description='Calculate and plot the average total scores for each sample')
    parser.add_argument('--file_paths', nargs='+', required=True, help='List of file paths to MultiMemote CSV files')
    return parser.parse_args()

def main():
    args = arg_parser()
    file_paths = args.file_paths
    # List of file paths to MultiMemote CSV files
    data = pd.read_csv(file_paths)
    # Dictionary to store averages for each sample
    averages = {}

    # Calculate averages for each sample based on file names
    for file_path in file_paths:
        sample_name = os.path.basename(file_path).split('_')[0]  # Extract sample name (e.g., M05-04)
        df = pd.read_csv(file_path)
        average_score = df['total_score'].mean()
        averages[sample_name] = average_score

    # Convert averages to a DataFrame
    average_df = pd.DataFrame(list(averages.items()), columns=['Sample', 'Average_Score']).sort_values(by='Sample')

    # Normalize scores to 100% and cap them
    average_df['Normalized_Score'] = average_df['Average_Score'] * 100
    average_df['Capped_Score'] = average_df['Normalized_Score'].clip(upper=100)

    # Plotting the horizontal bar graph
    plt.figure(figsize=(12, 8))
    plt.barh(average_df['Sample'], average_df['Capped_Score'], color='skyblue')
    plt.xlabel('Average Total Score (%)', fontsize=12)
    plt.ylabel('Sample', fontsize=12)
    plt.title('Capped Normalized Average Total Scores for Each Sample', fontsize=14)
    plt.xlim(0, 100)  # Set the range to go up to 100%
    plt.tight_layout()
    plt.show()
