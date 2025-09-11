import sys
import os
import argparse
import glob
from collections import defaultdict
import heapq
import statistics

import pyfaidx

# -1 is the last nt of the 3'UTR
NUE_START = -50
NUE_END = -5
CE_START = -15
CE_END = 20

def get_kmer_counts(sequences: list[str], region_start: int, region_end: int, kmer_size: int, downstream_nts: int, step_size: int = 1) -> dict:
    """
    Counts k-mers at each position in the specified region across all sequences.

    Args:
        sequences (list[str]): List of sequences to analyze.
        region_start (int): Start position of the region (where last nt of 3'UTR = -1, e.g., -50).
        region_end (int): End position of the region.
        kmer_size (int): Size of the k-mers to count.
        step_size (int): Step size for k-mer counting (=1: overlapping k-mers. =kmer_size: non-overlapping k-mers)
    """
    positional_counts = defaultdict(lambda: defaultdict(int))

    for seq in sequences:
        seq_len = len(seq)
        utr_len = seq_len - downstream_nts

        for i in range(0, seq_len - kmer_size + 1, step_size):
            pos = (i + kmer_size - 1) - utr_len # anchored to rightmost nt of kmer

            if region_start <= pos <= region_end:
                kmer = seq[i : i + kmer_size]
                positional_counts[kmer][pos] += 1

    return positional_counts


def rank_kmers_by_delta(kmer_counts: dict, top_n: int) -> list:
    """
    Ranks k-mers by the difference between their peak and median positional counts.
    Returns the top N k-mers with the highest delta.
    """
    heap = []
    for kmer, positions in kmer_counts.items():
        if not positions:
            continue

        counts = list(positions.values())

        peak_count = 0
        peak_pos = 0
        for pos, count in positions.items():
            if count > peak_count:
                peak_count = count
                peak_pos = pos
        
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
    parser = argparse.ArgumentParser(
        description="Analyses the NUE and CE regions of the given set of terminator sequences."
    )
    parser.add_argument("input_dir", help="Path to the directory containing terminator sequence FASTA files.")
    parser.add_argument(
        "-d", "--downstream-nts", type=int, default=50,
        help="Number of downstream nucleotides included in the terminators (default: 50)."
    )
    parser.add_argument(
        "-n", "--top-n", type=int, default=50,
        help="The number of most common k-mers to report (default: 50)."
    )
    parser.add_argument(
        "-k", "--kmer-size", type=int, default=6,
        help="The size of the k-mers to analyze (default: 6)."
    )
    parser.add_argument(
        "-m", "--min-3utr-length", type=int, default=50,
        help="Minimum length of the 3'UTR required to be included in the analysis (default: 50)."
    )
    parser.add_argument(
        "-s", "--step-size", type=int, default=1,
        help="Step size for k-mer counting: 1=overlapping (default), kmer_size for non-overlapping."
    )
    args = parser.parse_args()
    # TODO validate args
    # step_size should be <= kmer_size
    # min_3utr_length should be >= NUE_START absolute value
    # input_dir should exist
    # downstream_nts should be >= CE_END absolute value. or maybe not if want to analyse only 3'UTR

    # Find all fasta files
    fasta_files = glob.glob(os.path.join(args.input_dir, "*_terminators.fa"))
    if not fasta_files:
        print(f"ERROR: No '*_terminators.fa' files found in directory '{args.input_dir}'", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(fasta_files)} FASTA files")
    
    sequences = []
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
            
            sequences.append(seq)
            
    nue_pos_counts = get_kmer_counts(sequences, NUE_START, NUE_END, args.kmer_size, args.downstream_nts, args.step_size)
    ranked_nue_kmers = rank_kmers_by_delta(nue_pos_counts, args.top_n)

    ce_pos_counts = get_kmer_counts(sequences, CE_START, CE_END, args.kmer_size, args.downstream_nts, args.step_size)
    ranked_ce_kmers = rank_kmers_by_delta(ce_pos_counts, args.top_n)

    print(f"\n{skipped} of {num_terminators} ({(skipped/num_terminators * 100):.2f}%) terminator sequences skipped due to insufficient 3'UTR length")
    print_report("CE", ranked_ce_kmers, args.kmer_size)
    print_report("NUE", ranked_nue_kmers, args.kmer_size)

    print("\nFINISHED")


if __name__ == "__main__":
    main()
