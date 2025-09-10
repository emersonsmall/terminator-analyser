import sys
import os
import argparse
import glob
import textwrap
from collections import defaultdict
import heapq
import statistics

import pyfaidx

# -1 is the last nt of the 3'UTR
NUE_START = -35
NUE_END = -5
CE_START = -15
CE_END = 15

def get_kmer_counts(sequences: list[str], region_start: int, kmer_size: int, step_size: int = 1) -> dict:
    """
    Counts k-mers at each position in the specified region across all sequences.
    - step_size = 1: overlapping k-mers
    - step_size = kmer_size: non-overlapping k-mers
    """
    positional_counts = defaultdict(lambda: defaultdict(int))

    for seq in sequences:
        for i in range(0, len(seq) - kmer_size + 1, step_size):
            kmer = seq[i : i + kmer_size]

            pos = region_start + i + kmer_size - 1 # Counts anchored to rightmost nt of kmer

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
        description=textwrap.dedent("""
        Analyses the NUE and CE regions of the given set of terminator sequences.
        Finds the top N most common k-mers across all sequences using jellyfish.
        """)
    )
    parser.add_argument("input_dir", help="Path of the directory containing FASTA files.")
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

    # Find all fasta files
    fasta_files = glob.glob(os.path.join(args.input_dir, "*_terminators.fa"))
    if not fasta_files:
        print(f"ERROR: No '*_terminators.fa' files found in directory '{args.input_dir}'", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(fasta_files)} FASTA files")
    
    print(f"Extracting CE and NUE regions")
    ce_seqs = []
    nue_seqs = []
    skipped = 0
    num_terminators = 0

    for fasta_file in fasta_files:
        fa_records = pyfaidx.Fasta(fasta_file, as_raw=True)
        for record in fa_records:
            num_terminators += 1

            seq = str(record).upper()
            total_len = len(seq)
            utr_len = total_len - args.downstream_nts

            if utr_len < args.min_3utr_length:
                skipped += 1
                continue

            if utr_len >= abs(NUE_START) and args.downstream_nts >= CE_END:
                # Slice NUE window
                nue_start = utr_len - abs(NUE_START)
                nue_end = utr_len - abs(NUE_END)
                nue_seq = seq[nue_start - 1 : nue_end]
                nue_seqs.append(nue_seq)

                # Slice CE window
                ce_start = utr_len - abs(CE_START)
                ce_end = utr_len + CE_END
                ce_seq = seq[ce_start - 1 : ce_end]
                ce_seqs.append(ce_seq)
            else:
                skipped += 1
                continue

    nue_pos_counts = get_kmer_counts(nue_seqs, NUE_START, args.kmer_size, args.step_size)
    ranked_nue_kmers = rank_kmers_by_delta(nue_pos_counts, args.top_n)

    ce_pos_counts = get_kmer_counts(ce_seqs, CE_START, args.kmer_size, args.step_size)
    ranked_ce_kmers = rank_kmers_by_delta(ce_pos_counts, args.top_n)

    print(f"\n{skipped} of {num_terminators} ({(skipped/num_terminators * 100):.2f}%) terminator sequences skipped due to insufficient 3'UTR length")
    print_report("CE", ranked_ce_kmers, args.kmer_size)
    print_report("NUE", ranked_nue_kmers, args.kmer_size)

    print("FINISHED")


if __name__ == "__main__":
    main()
