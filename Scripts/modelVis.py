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

    # Read stats: columns are [sample, bin, mets, rxns, genes]
    df = pd.read_csv(args.stats, sep=' ', header=None, names=['sample','bin','mets','rxns','genes'])

    # Convert numerical columns
    df['mets'] = df['mets'].astype(int)
    df['rxns'] = df['rxns'].astype(int)
    df['genes'] = df['genes'].astype(int)

    # Group by sample to get number of bins per sample
    sample_counts = df.groupby('sample')['bin'].count().reset_index(name='num_bins')
    sample_counts = sample_counts.sort_values(by='num_bins', ascending=False)

    # Set up figure
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1])

    # Left plot: number of GEMs across samples
    ax0 = fig.add_subplot(gs[0])
    sns.barplot(data=sample_counts, x='num_bins', y='sample', orient='h', ax=ax0, color="#1f77b4")
    ax0.set_title('Number of GEMs across samples')
    ax0.set_xlabel('Number of GEMs')
    ax0.set_ylabel('Sample')
    ax0.invert_yaxis()  # optional if you prefer top->bottom

    # Right side with density plots
    gs_right = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs[1])

    # Metabolites density
    ax1 = fig.add_subplot(gs_right[0])
    sns.kdeplot(df['mets'], fill=True, color="#7fc97f", ax=ax1)
    ax1.set_title('Unique metabolites across GEMs')
    ax1.set_yticks([])

    # Reactions density
    ax2 = fig.add_subplot(gs_right[1])
    sns.kdeplot(df['rxns'], fill=True, color="#beaed4", ax=ax2)
    ax2.set_title('Reactions across GEMs')
    ax2.set_yticks([])

    # Genes density
    ax3 = fig.add_subplot(gs_right[2])
    sns.kdeplot(df['genes'], fill=True, color="#fdc086", ax=ax3)
    ax3.set_title('Genes across GEMs')
    ax3.set_yticks([])

    # Final styling and save
    plt.tight_layout()
    try:
        plt.savefig(args.output, format='pdf', bbox_inches='tight')
        print(f"Plot saved to {args.output}")
    except Exception as e:
        print(f"Error saving plot: {e}")

if __name__ == '__main__':
    main()
