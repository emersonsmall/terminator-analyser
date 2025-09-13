import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

NUE_X_MIN = -35
NUE_X_MAX = -5
CE_X_MIN = -15
CE_X_MAX = 15
CS_POSITION = 1

def plot_signal_distribution(ranked_kmers: list, counts_data: dict, region: str, out_file: str) -> None:
    """
    Generates and saves a line plot using a list of ranked k-mers and their positional counts.

    Args:
        ranked_kmers (list): List of dictionaries containing k-mer information.
        counts_data (dict): Dictionary with k-mer positional counts.
        region (str): Name of the region being analyzed.
        out_file (str): Path to save the output plot image.
    """
    top_n = len(ranked_kmers)

    top_kmers = [item['kmer'] for item in ranked_kmers]

    # Prepare data for plotting
    plot_data = {kmer: counts_data.get(kmer, {}) for kmer in top_kmers}
    df = pd.DataFrame(plot_data).fillna(0).astype(int)
    if df.empty:
        print("No data to plot.")
        return
    df = df.reindex(sorted(df.index))

    # Create plot
    fig, ax = plt.subplots(figsize=(12, 7))
    for kmer in df.columns:
        ax.plot(df.index, df[kmer], marker='o', linestyle='-', markersize=4, label=kmer)
    
    # Style plot
    ax.set_xlabel('Position relative to CS', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'Distribution of Top {top_n} Signals in {region}', fontsize=14, weight='bold')
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.set_ylim(bottom=0)

    ax.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.15),
        ncol=5
    )

    if region.upper() == "NUE":
        ax.set_xlim(NUE_X_MIN, NUE_X_MAX)
    elif region.upper() == "CE":
        ax.set_xlim(CE_X_MIN, CE_X_MAX)
        ax.axvline(CS_POSITION, linestyle='--', color='red')
    
    plt.tight_layout(rect=[0, 0, 0.9, 1])

    # Save plot
    plt.savefig(out_file, dpi=300)
    print(f"Plot saved to {out_file}")
    plt.close()
