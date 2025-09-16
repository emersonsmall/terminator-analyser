import sys
import os
import argparse
import glob
from collections import defaultdict
import heapq
import statistics

from plot_utils import plot_signal_distribution

import pyfaidx

# Region coordinates for k-mer count window
# -1 is the last nt of the 3'UTR, +1 is the first nt of the downstream region
NUE_START = -50
NUE_END = -1
CE_START = -10 # Loke paper uses -15, but -10 gives better results. If using -15, AAUAAA and other NUE signals dominate
CE_END = 20

# Plot constants
NUE_X_MIN = -35
NUE_X_MAX = -5
CE_X_MIN = -10
CE_X_MAX = 15

def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyses the NUE and CE regions of the given terminator sequences."
    )
    parser.add_argument(
        "input_path",
        help="Path to the terminator sequence FASTA file/s (filepath or directory path).")
    parser.add_argument(
        "-d", 
        "--downstream-nts", 
        type=int, 
        default=50,
        help="The number of nucleotides downstream of the CS included in the terminators (default: 50)."
    )
    parser.add_argument(
        "-n",
        "--top-n",
        type=int,
        default=20,
        help="The number of k-mers to report (default: 20)."
    )
    parser.add_argument(
        "-k", 
        "--kmer-size", 
        type=int, 
        default=6,
        help="The k-mer size (default: 6)."
    )
    parser.add_argument(
        "-m",
        "--min-3utr-length",
        type=int,
        default=100,
        help="The minimum length of 3'UTRs to be included in the analysis (default: 100)."
    )
    parser.add_argument(
        "-s",
        "--step-size",
        type=int,
        default=1,
        help="The step size for k-mer counting: 1=overlapping (default), kmer_size=non-overlapping."
    )
    args = parser.parse_args()

    if not os.path.exists(args.input_path):
        parser.error(f"Input path '{args.input_path}' does not exist.")

    if not args.min_3utr_length >= abs(NUE_START):
        parser.error(f"Minimum 3'UTR length must be at least {abs(NUE_START)}.")

    if not args.downstream_nts >= CE_END:
        parser.error(f"Downstream nts must be at least {CE_END}.")

    return args


def get_kmer_counts(sequences: list[str], region_start: int, region_end: int, kmer_size: int, downstream_nts: int, step_size: int = 1) -> dict:
    """
    Counts k-mers at each position in the specified region across all sequences.
    If region_start and region_end are 0, counts k-mers across the whole sequence.
    Args:
        sequences (list[str]): List of sequences to analyze.
        region_start (int): Start position of the region (where last nt of 3'UTR = -1, e.g., -50).
        region_end (int): End position of the region.
        kmer_size (int): Size of the k-mers to count.
        step_size (int): Step size for k-mer counting (=1: overlapping k-mers. =kmer_size: non-overlapping k-mers)
    """
    positional_counts = defaultdict(lambda: defaultdict(int))

    is_global = (region_start == 0 and region_end == 0)

    for seq in sequences:
        seq_len = len(seq)
        utr_len = seq_len - downstream_nts

        for i in range(0, seq_len - kmer_size + 1, step_size):
            pos = (i + kmer_size - 1) - utr_len # anchored to rightmost nt of kmer
            if pos >= 0:
                pos += 1 # +1 is the first nt of downstream region
            
            if is_global or (region_start <= pos <= region_end):
                kmer = seq[i : i + kmer_size]
                positional_counts[kmer][pos] += 1

    return positional_counts


def rank_kmers_by_delta(kmer_counts: dict, top_n: int) -> list:
    """
    Ranks k-mers by the difference between their peak and median counts.
    Returns the top N k-mers with the highest delta.
    """
    heap = []
    for kmer, positions in kmer_counts.items():
        if not positions:
            continue

        # Find peak count for this kmer
        peak_count = 0
        peak_pos = 0
        for pos, count in positions.items():
            if count > peak_count:
                peak_count = count
                peak_pos = pos
        
        # Calculate median count for this kmer
        counts = list(positions.values())
        median_count = statistics.median(counts)

        delta = peak_count - median_count

        item = (delta, peak_count, kmer, {
            "kmer": kmer,
            "delta": delta,
            "peak_count": peak_count,
            "peak_pos": peak_pos,
            "median_count": median_count
        })

        if len(heap) < top_n:
            heapq.heappush(heap, item)
        else:
            heapq.heappushpop(heap, item)
    
    top_items = sorted([item[3] for item in heap], 
                       key=lambda x: (x["delta"], x["peak_count"]), 
                       reverse=True)

    return top_items


def print_report(region_name: str, kmers: list, kmer_size: int):
    print("\n" + "=" * 40)
    print(f"Top {len(kmers)} K-mers for {region_name}")
    print("=" * 40)
    print(f"{'Rank':<5} | {'K-mer':<{kmer_size + 2}} | {'Delta':>15} | {'Peak Count':>12} | {'Peak Pos':>15}")
    print("-" * 40)
    for i, item in enumerate(kmers):
        rank = i + 1
        print(f"{rank:<5} | {item['kmer']:<{kmer_size + 2}} | {item['delta']:>15.1f} | {item['peak_count']:>12,} | {item['peak_pos']:>15}")


def main():
    args = get_args()

    # Find all fasta files
    fasta_files = []
    if os.path.isdir(args.input_path):
        fasta_files = glob.glob(os.path.join(args.input_path, "*_terminators.fa"))
        print(f"Found {len(fasta_files)} FASTA files")
    elif os.path.isfile(args.input_path):
        fasta_files = [args.input_path]
    
    if not fasta_files:
        print(f"ERROR: No valid FASTA file/s found at path '{args.input_path}'", file=sys.stderr)
        return 1
    
    terminators = []
    skipped = 0
    num_terminators = 0

    for file in fasta_files:
        fa_records = pyfaidx.Fasta(file, as_raw=True)
        for record in fa_records:
            num_terminators += 1

            seq = str(record).upper()
            total_len = len(seq)
            utr_len = total_len - args.downstream_nts

            if utr_len < args.min_3utr_length:
                skipped += 1
                continue
            
            terminators.append(seq)
    
    nue_counts = get_kmer_counts(terminators, NUE_START, NUE_END, args.kmer_size, args.downstream_nts, args.step_size)
    ranked_nue_kmers = rank_kmers_by_delta(nue_counts, args.top_n)

    ce_counts = get_kmer_counts(terminators, CE_START, CE_END, args.kmer_size, args.downstream_nts, args.step_size)
    ranked_ce_kmers = rank_kmers_by_delta(ce_counts, args.top_n)

    print(f"\n{skipped} of {num_terminators} ({(skipped/num_terminators * 100):.2f}%) terminator sequences skipped due to insufficient 3'UTR length")
    print_report("CE", ranked_ce_kmers, args.kmer_size)
    print_report("NUE", ranked_nue_kmers, args.kmer_size)

    plot_signal_distribution(
        ranked_nue_kmers, 
        nue_counts, 
        "NUE",
        NUE_X_MIN, 
        NUE_X_MAX,
        "NUE_signals_plot.png"
    )

    plot_signal_distribution(
        ranked_ce_kmers, 
        ce_counts, 
        "CE", 
        CE_X_MIN, 
        CE_X_MAX,
        "CE_signals_plot.png"
    )

    print("\nFINISHED")


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
