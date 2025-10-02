import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CS_POSITION = 1
LINESTYLES = ['-', '--', '-.', ':']
MARKERS = ['o', 's', '^', 'D', 'v', 'P', '*', 'X']

def plot_signal_distribution(ranked_kmers: list, counts_data: dict, region: str, x_min: int, x_max: int, out_file: str) -> None:
    """
    Generates and saves a line plot using a list of ranked k-mers and their positional counts.

    Args:
        ranked_kmers (list): List of dictionaries containing k-mer information.
        counts_data (dict): Dictionary with k-mer positional counts.
        region (str): Name of the region being analyzed.
        out_file (str): Path to save the output plot image.
    """
    tick_interval = 5
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
    fig, ax = plt.subplots(figsize=(12, 8))
    for i, kmer in enumerate(df.columns):
        linestyle = LINESTYLES[i % len(LINESTYLES)]
        marker = MARKERS[i % len(MARKERS)]
        ax.plot(
            df.index,
            df[kmer],
            label=kmer,
            linestyle=linestyle,
            marker=marker,
            markersize=5
        )
    
    # Style plot
    ax.set_xlabel('Position relative to CS', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title(f'Distribution of Top {top_n} Signals in {region}', fontsize=14, weight='bold')
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    ax.set_ylim(bottom=0)

    ax.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.15),
        ncols=5
    )

    ax.set_xlim(x_min, x_max)
    if region.upper() == "CE":
        ticks = list(range(x_min, x_max + 1, 1))
        ticks.remove(0)
        ticks.append(1)
        ticks.sort()
        ax.set_xticks(ticks)
    else:
        ax.set_xticks(range(x_min, x_max + 1, tick_interval))
    
    plt.tight_layout()

    # Save plot
    plt.savefig(out_file, dpi=300)
    print(f"Plot saved to {out_file}")
    plt.close()
