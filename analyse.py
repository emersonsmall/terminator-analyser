# Built-in libraries
import sys
import os
import argparse
import glob
import heapq
import statistics
from multiprocessing import Pool
from functools import partial
import datetime

# External libraries
import pyfaidx

# Local modules
from plots import plot_signal_distribution

# Coordinates: -1 is the last nt of the 3'UTR, +1 is the first nt of the downstream region

# Analysis window coordinates
NUE_START = -50
NUE_END = -1
CE_START = (
    -10
)  # Loke et al. uses -15, but -10 gives better results. NUE signals dominate if using -15
CE_END = 20

PLOT_NUE_X_MIN = -35
PLOT_NUE_X_MAX = -5
PLOT_CE_X_MIN = -10
PLOT_CE_X_MAX = 15
PLOT_FILE_SUFFIX = "_signals_plot.png"


def main():
    run_analysis(_get_args())


def run_analysis(args: argparse.Namespace) -> None:
    """Runs an analysis on terminator sequence fasta files.

    Args:
        args: Parsed command-line arguments.
    """

    try:
        # Find all fasta files
        fasta_files = []
        if os.path.isdir(args.input_path):
            included_accessions = getattr(args, "included_accessions", None)
            fasta_files = _get_fasta_files(args.input_path, included_accessions)
            print(f"Found {len(fasta_files)} FASTA files")
        elif os.path.isfile(args.input_path):
            fasta_files = [args.input_path]

        if not fasta_files:
            print(
                f"ERROR: No valid FASTA file/s found at path '{args.input_path}'",
                file=sys.stderr,
            )

        # process each file in parallel
        worker = partial(_process_terminator_fasta, args=args)
        total_nue_counts = {}
        total_ce_counts = {}
        total_skipped = 0
        total_terminators = 0
        accession_stats = {}

        with Pool() as pool:
            for (
                nue_counts,
                ce_counts,
                num_skipped,
                num_terminators,
                accession
            ) in pool.imap_unordered(worker, fasta_files):
                _merge_counts(total_nue_counts, nue_counts)
                _merge_counts(total_ce_counts, ce_counts)
                total_skipped += num_skipped
                total_terminators += num_terminators
                accession_stats[accession] = {
                    "total": num_terminators,
                    "included": num_terminators - num_skipped,
                    "skipped": num_skipped,
                }

        ranked_nue_kmers = _rank_kmers(total_nue_counts, args.num_kmers)
        ranked_ce_kmers = _rank_kmers(total_ce_counts, args.num_kmers)

        print(
            f"\n{total_skipped} of {total_terminators} ({(total_skipped/total_terminators * 100):.2f}%) terminators skipped"
        )

        os.makedirs(args.results_dir, exist_ok=True)

        # Prepare metadata for reports
        taxon = getattr(args, "taxon", None)
        report_metadata = {
            "taxon": taxon,
            "generated_at": datetime.datetime.now(),
            "kmer_size": args.kmer_size,
            "min_3utr_length": args.min_3utr_length,
            "num_downstream_nts": args.num_downstream_nts,
            "step_size": args.step_size,
            "total_terminators": total_terminators,
            "skipped_terminators": total_skipped,
            "included_terminators": total_terminators - total_skipped,
            "num_accessions": len(accession_stats),
            "accession_stats": accession_stats,
        }

        _save_report(
            "CE",
            ranked_ce_kmers,
            os.path.join(args.results_dir, "CE_report.txt"),
            report_metadata,
            CE_START,
            CE_END,
        )
        _save_report(
            "NUE",
            ranked_nue_kmers,
            os.path.join(args.results_dir, "NUE_report.txt"),
            report_metadata,
            NUE_START,
            NUE_END,
        )

        nue_plot_path = os.path.join(args.results_dir, "NUE" + PLOT_FILE_SUFFIX)
        ce_plot_path = os.path.join(args.results_dir, "CE" + PLOT_FILE_SUFFIX)

        plot_signal_distribution(
            ranked_nue_kmers,
            total_nue_counts,
            "NUE",
            PLOT_NUE_X_MIN,
            PLOT_NUE_X_MAX,
            nue_plot_path,
        )
        plot_signal_distribution(
            ranked_ce_kmers,
            total_ce_counts,
            "CE",
            PLOT_CE_X_MIN,
            PLOT_CE_X_MAX,
            ce_plot_path,
        )

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)


def add_analyse_args(
    parser: argparse.ArgumentParser, is_standalone: bool = True
) -> None:
    """Adds command-line arguments for the `analyse` command to the given parser.

    Args:
        parser: The argument parser to which the arguments will be added.
        is_standalone: Whether to include standalone execution arguments. Defaults to True.
    """

    parser.add_argument(
        "-n",
        "--num-kmers",
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
    parser.add_argument(
        "--results-dir",
        default=os.path.join("out", "results"),
        help="Path to the output directory (default: ./out/results).",
    )

    if is_standalone:
        parser.add_argument(
            "input_path",
            help="Path to the terminator sequence FASTA file/s (filepath or directory path).",
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
    
    if not args.num_downstream_nts >= CE_END:
        parser.error(f"Downstream nts must be at least {CE_END}.")
    
    return args


def _count_kmers(
    sequences: list[str],
    region_start: int,
    region_end: int,
    kmer_size: int,
    num_downstream_nts: int,
    step_size: int,
) -> dict:
    """
    Counts k-mers at each position within the specified region across all sequences.
    If region_start and region_end are 0, counts k-mers across the whole sequence.

    Args:
        sequences: List of sequences to analyse.
        region_start: Start position of the region.
        region_end: End position of the region.
        kmer_size: Size of the k-mers.
        step_size: Step size for k-mer counting.

    Returns:
        Dictionary of k-mer counts by position { kmer: { pos1: count, pos2: count } }.
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

    positional_counts = {}

    is_global = region_start == 0 and region_end == 0

    for seq in sequences:
        seq_len = len(seq)
        utr_len = seq_len - num_downstream_nts

        for i in range(0, seq_len - kmer_size + 1, step_size):
            pos = (i + kmer_size - 1) - utr_len  # anchored to rightmost nt of kmer
            if pos >= 0:
                pos += 1  # +1 is the first nt of downstream region

            if is_global or (region_start <= pos <= region_end):
                kmer = seq[i : i + kmer_size]
                if kmer not in positional_counts:
                    positional_counts[kmer] = {}
                if pos not in positional_counts[kmer]:
                    positional_counts[kmer][pos] = 0
                positional_counts[kmer][pos] += 1

    return positional_counts


def _process_terminator_fasta(fasta_fpath: str, args: argparse.Namespace) -> tuple:
    """
    Processes a single terminator FASTA file, counting kmers in NUE and CE regions.

    Args:
        fasta_fpath: Path to the FASTA file.
        args: Parsed command-line arguments.

    Returns:
        A tuple containing (nue_counts, ce_counts, num_skipped, num_terminators)
    """

    terminators = []
    num_skipped = 0
    num_terminators = 0

    base = os.path.basename(fasta_fpath)
    accession = os.path.splitext(base)[0]

    fa_records = pyfaidx.Fasta(fasta_fpath, as_raw=True)
    for record in fa_records:
        num_terminators += 1

        seq = str(record).upper()
        total_len = len(seq)
        utr_len = total_len - args.num_downstream_nts

        if utr_len < args.min_3utr_length:
            num_skipped += 1
            continue

        terminators.append(seq)

    nue_counts = _count_kmers(
        terminators,
        NUE_START,
        NUE_END,
        args.kmer_size,
        args.num_downstream_nts,
        args.step_size,
    )
    ce_counts = _count_kmers(
        terminators,
        CE_START,
        CE_END,
        args.kmer_size,
        args.num_downstream_nts,
        args.step_size,
    )

    return nue_counts, ce_counts, num_skipped, num_terminators, accession


def _merge_counts(target: dict, src: dict) -> None:
    """
    Merges k-mer counts from source into target.

    Args:
        target: Target dictionary to merge into.
        src: Source dictionary to merge from.
    """

    for kmer, pos_map in src.items():
        if kmer not in target:
            target[kmer] = {}
        for pos, count in pos_map.items():
            if pos not in target[kmer]:
                target[kmer][pos] = 0
            target[kmer][pos] += count


def _rank_kmers(kmer_counts: dict, n: int) -> list:
    """
    Ranks k-mers by the difference between their peak and median counts.
    Returns the top N k-mers with the highest delta.

    Args:
        kmer_counts: Dictionary of k-mer counts by position.
        n: Number of top k-mers to return.

    Returns:
        list: List of top N k-mers with their statistics.
    """

    assert isinstance(kmer_counts, dict)
    assert isinstance(n, int) and n > 0, "Top N must be a positive integer."

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

        if len(heap) < n:
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
def _get_fasta_files(input_dir: str, included_accessions: set[str] | None) -> list[str]:
    assert os.path.isdir(input_dir), f"Input path '{input_dir}' is not a directory."

    candidates = glob.glob(os.path.join(input_dir, "*.fa"))

    results = []
    for fpath in candidates:
        base = os.path.basename(fpath)
        accession = os.path.splitext(base)[0]
        # if filter set provided, include only those accessions
        if included_accessions and accession not in included_accessions:
            continue
        results.append(fpath)
    
    return results


def _save_report(
    region_name: str,
    ranked_kmers: list,
    out_fpath: str,
    metadata: dict,
    region_start: int,
    region_end: int,
) -> None:
    """Writes a formatted analysis report to a file.

    Args:
        region_name: Name of the region (e.g., "NUE" or "CE").
        ranked_kmers: List of top N k-mers with their statistics.
        out_fpath: Path to the output file.
        metadata: Dictionary containing analysis metadata.
        region_start: Start position of the analysis region.
        region_end: End position of the analysis region.
    """

    lines = []

    # Header
    lines.append("=" * 80 + "\n")
    lines.append(f"Terminator Analysis Report - {region_name} Region\n")
    lines.append("=" * 80 + "\n\n")

    # Metadata section
    lines.append("ANALYSIS METADATA\n")
    lines.append("-" * 80 + "\n")
    if metadata.get("taxon"):
        lines.append(f"Taxon:                  {metadata['taxon']}\n")
    lines.append(
        f"Generated:              {metadata['generated_at'].strftime('%Y-%m-%d %H:%M:%S')}\n"
    )
    lines.append(f"Analysis Window:        {region_start} to {region_end}\n")
    lines.append(f"K-mer Size:             {metadata['kmer_size']}\n")
    lines.append(f"Step Size:              {metadata['step_size']}\n")
    lines.append(f"Min 3'UTR Length:       {metadata['min_3utr_length']} nt\n")
    lines.append(f"Num Downstream nt:      {metadata['num_downstream_nts']} nt\n")
    skip_pct = (
        (metadata["skipped_terminators"] / metadata["total_terminators"] * 100)
        if metadata["total_terminators"] > 0
        else 0
    )
    lines.append(f"Skip Rate:              {skip_pct:.2f}%\n")
    lines.append(f"Number of Accessions:   {metadata['num_accessions']}\n")
    lines.append("\n")

    # Per-accession breakdown table
    if metadata["accession_stats"]:
        lines.append("PER-ACCESSION BREAKDOWN\n")
        lines.append("-" * 80 + "\n")

        accession_stats = metadata["accession_stats"]
        
        # Totals for TOTAL row
        total_all = metadata["total_terminators"]
        total_included = metadata["included_terminators"]
        total_skipped = metadata["skipped_terminators"]

        header = {
            "accession": "Accession",
            "total": "Total",
            "included": "Included",
            "skipped": "Skipped",
            "pct": "% of Included",
        }

        # Calculate column widths including header, data rows, AND TOTAL row
        max_widths = {
            "accession": len(header["accession"]),
            "total": len(header["total"]),
            "included": len(header["included"]),
            "skipped": len(header["skipped"]),
            "pct": len(header["pct"]),
        }

        # Include accession data rows
        for accession, stats in accession_stats.items():
            max_widths["accession"] = max(max_widths["accession"], len(accession))
            max_widths["total"] = max(max_widths["total"], len(f"{stats['total']:,}"))
            max_widths["included"] = max(
                max_widths["included"], len(f"{stats['included']:,}")
            )
            max_widths["skipped"] = max(
                max_widths["skipped"], len(f"{stats['skipped']:,}")
            )

        # Include TOTAL row in width calculations
        max_widths["accession"] = max(max_widths["accession"], len("TOTAL"))
        max_widths["total"] = max(max_widths["total"], len(f"{total_all:,}"))
        max_widths["included"] = max(max_widths["included"], len(f"{total_included:,}"))
        max_widths["skipped"] = max(max_widths["skipped"], len(f"{total_skipped:,}"))

        # Now build format string with final widths
        fmt_string = (
            f"{{accession:<{max_widths['accession']}}} | "
            f"{{total:>{max_widths['total']}}} | "
            f"{{included:>{max_widths['included']}}} | "
            f"{{skipped:>{max_widths['skipped']}}} | "
            f"{{pct:>{max_widths['pct']}}}"
        )

        separator = (
            f"{'-' * max_widths['accession']} | "
            f"{'-' * max_widths['total']} | "
            f"{'-' * max_widths['included']} | "
            f"{'-' * max_widths['skipped']} | "
            f"{'-' * max_widths['pct']}"
        )

        # Print header
        lines.append(fmt_string.format(**header) + "\n")
        lines.append(separator + "\n")

        # Print data rows sorted by total count descending
        sorted_accessions = sorted(
            accession_stats.items(), key=lambda x: x[1]["total"], reverse=True
        )
        for accession, stats in sorted_accessions:
            pct = (
                (stats["included"] / metadata["included_terminators"] * 100)
                if metadata["included_terminators"] > 0
                else 0
            )
            lines.append(
                fmt_string.format(
                    accession=accession,
                    total=f"{stats['total']:,}",
                    included=f"{stats['included']:,}",
                    skipped=f"{stats['skipped']:,}",
                    pct=f"{pct:.2f}%",
                )
                + "\n"
            )

        # Bottom-line totals
        lines.append(separator + "\n")
        lines.append(
            fmt_string.format(
                accession="TOTAL",
                total=f"{total_all:,}",
                included=f"{total_included:,}",
                skipped=f"{total_skipped:,}",
                pct="100.00%",
            )
            + "\n"
        )
        lines.append("\n")

    # Top K-mers section
    lines.append(f"TOP {len(ranked_kmers)} K-MERS\n")
    lines.append("-" * 80 + "\n")

    if not ranked_kmers:
        lines.append("No k-mers found.\n")
    else:
        kmer_size = len(ranked_kmers[0]["kmer"])

        # Calculate column widths
        header_kmers = {
            "kmer": "K-mer",
            "delta": "Delta",
            "peak": "Peak Count",
            "median": "Median Count",
            "pos": "Peak Pos",
        }

        max_widths_kmers = {
            "kmer": max(kmer_size, len(header_kmers["kmer"])),
            "delta": len(header_kmers["delta"]),
            "peak": len(header_kmers["peak"]),
            "median": len(header_kmers["median"]),
            "pos": len(header_kmers["pos"]),
        }

        for item in ranked_kmers:
            max_widths_kmers["delta"] = max(
                max_widths_kmers["delta"], len(f"{item['delta']:.1f}")
            )
            max_widths_kmers["peak"] = max(
                max_widths_kmers["peak"], len(f"{item['peak_count']:,}")
            )
            max_widths_kmers["median"] = max(
                max_widths_kmers["median"], len(f"{item['median_count']:.1f}")
            )
            max_widths_kmers["pos"] = max(
                max_widths_kmers["pos"], len(str(item["peak_pos"]))
            )

        fmt_string_kmers = (
            f"{{kmer:<{max_widths_kmers['kmer']}}} | "
            f"{{delta:>{max_widths_kmers['delta']}}} | "
            f"{{peak:>{max_widths_kmers['peak']}}} | "
            f"{{median:>{max_widths_kmers['median']}}} | "
            f"{{pos:>{max_widths_kmers['pos']}}}"
        )

        separator_kmers = (
            f"{'-' * max_widths_kmers['kmer']} | "
            f"{'-' * max_widths_kmers['delta']} | "
            f"{'-' * max_widths_kmers['peak']} | "
            f"{'-' * max_widths_kmers['median']} | "
            f"{'-' * max_widths_kmers['pos']}"
        )

        lines.append(fmt_string_kmers.format(**header_kmers) + "\n")
        lines.append(separator_kmers + "\n")

        for item in ranked_kmers:
            lines.append(
                fmt_string_kmers.format(
                    kmer=item["kmer"],
                    delta=f"{item['delta']:.1f}",
                    peak=f"{item['peak_count']:,}",
                    median=f"{item['median_count']:.1f}",
                    pos=str(item["peak_pos"]),
                )
                + "\n"
            )

    lines.append("\n")
    lines.append("=" * 80 + "\n")

    with open(out_fpath, "w") as f:
        f.writelines(lines)

    print(f"Report saved to '{out_fpath}'")


# --- STANDALONE EXECUTION ---
if __name__ == "__main__":
    sys.exit(main())
