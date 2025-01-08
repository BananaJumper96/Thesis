import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Generate visualizations from GEMs stats.')
    parser.add_argument('--stats', required=True, help='Path to the GEMs.stats file')
    parser.add_argument('--output', required=True, help='File path to save the output plot (e.g., output.pdf)')
    return parser.parse_args()

def main():
    args = parse_args()

    # Read the GEMs.stats file into a pandas DataFrame
    gems = pd.read_csv(args.stats, sep=' ', header=None)
    gems.columns = ['bin', 'mets', 'rxns', 'genes']
    gems['sample'] = gems['bin'].str.replace('_.*$', '', regex=True)

    # Convert numerical columns to integers
    gems['mets'] = gems['mets'].astype(int)
    gems['rxns'] = gems['rxns'].astype(int)
    gems['genes'] = gems['genes'].astype(int)

    # Prepare data for the samples plot
    sample_counts = gems['sample'].value_counts().reset_index()
    sample_counts.columns = ['sample', 'n']
    sample_counts = sample_counts.sort_values('n', ascending=False)

    # Create the figure and specify the grid layout
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])

    # Left plot: Bar plot of the number of GEMs across samples
    ax0 = plt.subplot(gs[0])
    sns.barplot(x='n', y='sample', data=sample_counts, orient='h', ax=ax0, color="#1f77b4")
    ax0.set_title('Number of GEMs across samples')
    ax0.set_xlabel('Number of GEMs carved')
    ax0.set_ylabel('Sample ID')
    ax0.invert_yaxis()  # To match the order in the R plot

    # Right plots: Density plots arranged vertically
    gs_right = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[1])

    # Metabolites density plot
    ax1 = plt.subplot(gs_right[0])
    sns.kdeplot(gems['mets'], fill=True, color="#7fc97f", ax=ax1)
    ax1.set_title('Unique metabolites across GEMs')
    ax1.set_ylabel('')
    ax1.set_yticks([])

    # Reactions density plot
    ax2 = plt.subplot(gs_right[1])
    sns.kdeplot(gems['rxns'], fill=True, color="#beaed4", ax=ax2)
    ax2.set_title('Reactions across GEMs')
    ax2.set_ylabel('')
    ax2.set_yticks([])

    # Genes density plot
    ax3 = plt.subplot(gs_right[2])
    sns.kdeplot(gems['genes'], fill=True, color="#fdc086", ax=ax3)
    ax3.set_title('Genes across GEMs')
    ax3.set_ylabel('')
    ax3.set_yticks([])

    # Adjust layout to prevent overlap
    plt.tight_layout()

    try:
        plt.savefig(args.output, format='pdf', bbox_inches='tight')
        print(f"Plot saved to {args.output}")
    except Exception as e:
        print(f"Error saving plot: {e}")


    # Optional: Display the plot
    # plt.show()

if __name__ == '__main__':
    main()
