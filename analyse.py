# Built-in libraries
import sys
import os
import argparse
import glob
from collections import defaultdict
import heapq
import statistics

# External libraries
import pyfaidx

# Local modules
from extract import TERMINATOR_FILE_SUFFIX
from plots import plot_signal_distribution


# TODO: add more visualisations, like Loke Figure 1C, 1D and Figure 3.
# TODO: add metrics to measure conservation - % occurance
# TODO: add feature to identify unique signals in a target gene. Compare gene terminator vs rest of genome vs whole genus
# TODO: Filter in one place - include all terminators, and then skip in analysis?
# TODO: validate args comprehensively and in full pipeline as well as standalone
# TODO: parallelise analysis (currently only extract is parallelised?)
# TODO: terminate gracefully
# TODO: generalise CE and NUE as 'analysis windows'

# Coordinates: -1 is the last nt of the 3'UTR, +1 is the first nt of the downstream region

# Analysis window coordinates
NUE_START = -50
NUE_END = -1
CE_START = (
    -10
)  # Loke paper uses -15, but -10 gives better results. AAUAAA and other NUE signals dominate if using -15
CE_END = 20

PLOT_NUE_X_MIN = -35
PLOT_NUE_X_MAX = -5
PLOT_CE_X_MIN = -10
PLOT_CE_X_MAX = 15
SIGNALS_PLOT_FILENAME = "_signals_plot.png"


def run_analysis(args: argparse.Namespace) -> int:
    """Runs a full analysis of terminator sequences.

    Args:
        args: Parsed command-line arguments.

    Returns:
        int: Exit code (0 for success, 1 for failure).
    """

    try:
        # Find all fasta files
        fasta_files = []
        if os.path.isdir(args.input_path):
            fasta_files = glob.glob(
                os.path.join(args.input_path, f"*{TERMINATOR_FILE_SUFFIX}")
            )
            print(f"Found {len(fasta_files)} FASTA files")
        elif os.path.isfile(args.input_path):
            fasta_files = [args.input_path]

        if not fasta_files:
            print(
                f"ERROR: No valid FASTA file/s found at path '{args.input_path}'",
                file=sys.stderr,
            )
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

        nue_counts = _count_kmers(
            terminators,
            NUE_START,
            NUE_END,
            args.kmer_size,
            args.downstream_nts,
            args.step_size,
        )
        ce_counts = _count_kmers(
            terminators,
            CE_START,
            CE_END,
            args.kmer_size,
            args.downstream_nts,
            args.step_size,
        )

        ranked_nue_kmers = _rank_kmers(nue_counts, args.top_n)
        ranked_ce_kmers = _rank_kmers(ce_counts, args.top_n)

        print(
            f"\n{skipped} of {num_terminators} ({(skipped/num_terminators * 100):.2f}%) terminator sequences skipped due to insufficient 3'UTR length"
        )
        _save_report(
            "CE",
            ranked_ce_kmers,
            args.kmer_size,
            os.path.join(args.output_dir, "CE_report.txt"),
        )
        _save_report(
            "NUE",
            ranked_nue_kmers,
            args.kmer_size,
            os.path.join(args.output_dir, "NUE_report.txt"),
        )

        os.makedirs(args.output_dir, exist_ok=True)

        nue_plot_path = os.path.join(args.output_dir, "NUE" + SIGNALS_PLOT_FILENAME)
        ce_plot_path = os.path.join(args.output_dir, "CE" + SIGNALS_PLOT_FILENAME)

        plot_signal_distribution(
            ranked_nue_kmers,
            nue_counts,
            "NUE",
            PLOT_NUE_X_MIN,
            PLOT_NUE_X_MAX,
            nue_plot_path,
        )

        plot_signal_distribution(
            ranked_ce_kmers, ce_counts, "CE", PLOT_CE_X_MIN, PLOT_CE_X_MAX, ce_plot_path
        )

        return 0
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


def add_analyse_args(parser: argparse.ArgumentParser, standalone: bool = True) -> None:
    """Adds command-line arguments for the `analyse` command to the given parser.

    Args:
        parser: The argument parser to which the arguments will be added.
        standalone: Whether to include standalone execution arguments. Defaults to True.
    """

    parser.add_argument(
        "-n",
        "--top-n",
        type=int,
        default=20,
        help="Number of k-mers to report (default: 20).",
    )
    parser.add_argument(
        "-k", "--kmer-size", type=int, default=6, help="K-mer size (default: 6)."
    )
    parser.add_argument(
        "-m",
        "--min-3utr-length",
        type=int,
        default=100,
        help="Minimum length of 3'UTRs to be included in the analysis (default: 100).",
    )
    parser.add_argument(
        "-s",
        "--step-size",
        type=int,
        default=1,
        help="Step size for k-mer counting (default: 1).",
    )

    if standalone:
        parser.add_argument(
            "input_path",
            help="Path to the terminator sequence FASTA file/s (filepath or directory path).",
        )
        parser.add_argument(
            "-o",
            "--output-dir",
            default=os.path.join("out", "plots"),
            help="Path to the output directory (default: ./out/plots).",
        )
        parser.add_argument(
            "-d",
            "--downstream-nts",
            type=int,
            default=50,
            help="Number of nucleotides downstream of the CS included in the terminators (default: 50).",
        )


def _get_args() -> argparse.Namespace:
    """Parses and validates command-line arguments for standalone execution.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """

    parser = argparse.ArgumentParser(
        description="Analyses the NUE and CE regions of the given terminator sequences."
    )
    add_analyse_args(parser)

    args = parser.parse_args()
    if not os.path.exists(args.input_path):
        parser.error(f"'{args.input_path}' does not exist.")
    if not args.min_3utr_length >= abs(NUE_START):
        parser.error(f"Minimum 3'UTR length must be at least {abs(NUE_START)}.")
    if not args.downstream_nts >= CE_END:
        parser.error(f"Downstream nts must be at least {CE_END}.")
    return args


def _count_kmers(
    sequences: list[str],
    region_start: int,
    region_end: int,
    kmer_size: int,
    downstream_nts: int,
    step_size: int = 1,
) -> dict:
    """
    Counts k-mers at each position within the specified region across all sequences.
    If region_start and region_end are 0, counts k-mers across the whole sequence.

    Args:
        sequences: List of sequences to analyse.
        region_start: Start position of the region (-1 is last nt of 3'UTR, +1 is first nt of downstream region).
        region_end: End position of the region.
        kmer_size: Size of the k-mers.
        step_size: Step size for k-mer counting (1: overlapping k-mers. kmer_size: non-overlapping k-mers)

    Returns:
        dict: Dictionary of k-mer counts by position { kmer: { pos1: count, pos2: count } }.
    """

    assert isinstance(sequences, list) and all(
        isinstance(s, str) for s in sequences
    ), "Sequences must be a list of strings."
    assert isinstance(region_start, int) and isinstance(
        region_end, int
    ), "Region start and end must be integers."
    assert (
        region_start <= region_end
    ), "Region start must be less than or equal to region end."
    assert (
        isinstance(kmer_size, int) and kmer_size > 0
    ), "K-mer size must be a positive integer."
    assert (
        isinstance(step_size, int) and step_size > 0
    ), "Step size must be a positive integer."

    positional_counts = defaultdict(lambda: defaultdict(int))

    is_global = region_start == 0 and region_end == 0

    for seq in sequences:
        seq_len = len(seq)
        utr_len = seq_len - downstream_nts

        for i in range(0, seq_len - kmer_size + 1, step_size):
            pos = (i + kmer_size - 1) - utr_len  # anchored to rightmost nt of kmer
            if pos >= 0:
                pos += 1  # +1 is the first nt of downstream region

            if is_global or (region_start <= pos <= region_end):
                kmer = seq[i : i + kmer_size]
                positional_counts[kmer][pos] += 1

    return positional_counts


def _rank_kmers(kmer_counts: dict, top_n: int) -> list:
    """
    Ranks k-mers by the difference between their peak and median counts.
    Returns the top N k-mers with the highest delta.

    Args:
        kmer_counts: Dictionary of k-mer counts by position.
        top_n: Number of top k-mers to return.
    """

    assert isinstance(kmer_counts, dict)
    assert isinstance(top_n, int) and top_n > 0, "Top N must be a positive integer."

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

        item = (
            delta,
            peak_count,
            kmer,
            {
                "kmer": kmer,
                "delta": delta,
                "peak_count": peak_count,
                "peak_pos": peak_pos,
                "median_count": median_count,
            },
        )

        if len(heap) < top_n:
            heapq.heappush(heap, item)
        else:
            heapq.heappushpop(heap, item)

    top_items = sorted(
        [item[3] for item in heap],
        key=lambda x: (x["delta"], x["peak_count"]),
        reverse=True,
    )

    return top_items


# --- HELPER FUNCTIONS ---
def _save_report(region_name: str, kmers: list, kmer_size: int, out_fpath: str) -> None:
    """Writes a formatted report to a file.

    Args:
        region_name: Name of the region (e.g., "NUE" or "CE").
        kmers: List of top k-mers with their statistics.
        kmer_size: Size of the k-mers.
    """

    lines = []
    lines.append("=" * 50 + "\n")
    lines.append(f"Top {len(kmers)} K-mers for {region_name}\n")
    lines.append("=" * 50 + "\n")
    lines.append(
        f"{'Rank':<5} | {'K-mer':<{kmer_size + 2}} | {'Delta':>8} | {'Peak Count':>10} | {'Peak Pos':>8}\n"
    )
    lines.append("-" * 50 + "\n")
    for i, item in enumerate(kmers, start=1):
        lines.append(
            f"{i:<5} | {item['kmer']:<{kmer_size + 2}} | {item['delta']:>8.1f} | {item['peak_count']:>10,} | {item['peak_pos']:>8}\n"
        )

    with open(out_fpath, "w") as f:
        f.writelines(lines)


# --- STANDALONE EXECUTION ---
def main():
    return run_analysis(_get_args())


if __name__ == "__main__":
    sys.exit(main())
