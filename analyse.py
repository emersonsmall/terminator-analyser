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
import traceback

# External libraries
import pyfaidx

# Local modules
from plots import plot_signal_distribution

# Coordinates: -1 is the last nt of the 3'UTR, +1 is the first nt of the downstream region

# Analysis window coordinates - used for counting kmers
NUE_ANALYSIS_WINDOW_START = -50
NUE_ANALYSIS_WINDOW_END = -1
CE_ANALYSIS_WINDOW_START = (
    -10
)  # Loke et al. uses -15, but NUE signals dominate if using -15
CE_ANALYSIS_WINDOW_END = 20

# Expected region coordinates - used for % occurrence of signals
NUE_EXPECTED_START = -30
NUE_EXPECTED_END = -13
CE_EXPECTED_START = -10
CE_EXPECTED_END = 10

PLOT_NUE_X_MIN = -35
PLOT_NUE_X_MAX = -5
PLOT_CE_X_MIN = -10
PLOT_CE_X_MAX = 15
PLOT_FILE_SUFFIX = "_signals_plot.png"


def main():
    run_analysis(_get_args())


def run_analysis(args: argparse.Namespace) -> None:
    try:
        # Find all fasta files
        fasta_files = []
        if os.path.isdir(args.input_path):
            genomes = getattr(args, "genomes", {})
            fasta_files = _get_fasta_files(args.input_path, list(genomes.keys()))
            print(f"Found {len(fasta_files)} FASTA files")
        elif os.path.isfile(args.input_path):
            fasta_files = [args.input_path]

        if not fasta_files:
            print(
                f"ERROR: No valid FASTA file/s found at path '{args.input_path}'",
                file=sys.stderr,
            )
            sys.exit(1)

        # process each file in parallel
        worker = partial(_process_terminator_fasta, args=args)
        total_nue_counts = {}
        total_ce_counts = {}
        total_nue_occurrence_counts = {}
        total_ce_occurrence_counts = {}
        total_skipped = 0
        total_terminators = 0
        accession_info = {}

        # TODO: do window occurrence counts need to be separate? - probably because the % occurrence is based on the window

        with Pool() as pool:
            for result in pool.imap_unordered(worker, fasta_files):
                (
                    nue_counts,
                    ce_counts,
                    nue_occurrence_counts,
                    ce_occurrence_counts,
                    num_skipped,
                    num_terminators,
                    accession,
                ) = result

                # merge positional counts
                _merge_counts(total_nue_counts, nue_counts)
                _merge_counts(total_ce_counts, ce_counts)

                # merge occurrence counts
                for kmer, count in nue_occurrence_counts.items():
                    total_nue_occurrence_counts[kmer] = total_nue_occurrence_counts.get(kmer, 0) + count
                for kmer, count in ce_occurrence_counts.items():
                    total_ce_occurrence_counts[kmer] = total_ce_occurrence_counts.get(kmer, 0) + count

                total_skipped += num_skipped
                total_terminators += num_terminators
                accession_info[accession] = {
                    "organism_name": args.genomes.get(accession, "Unknown"),
                    "total": num_terminators,
                    "included": num_terminators - num_skipped,
                    "skipped": num_skipped,
                }

        top_nue_kmers = _rank_kmers(total_nue_counts, args.num_kmers)
        top_ce_kmers = _rank_kmers(total_ce_counts, args.num_kmers)

        # calculate % occurrence for top kmers
        total_included = total_terminators - total_skipped
        for kmer_info in top_nue_kmers:
            occurrence_count = total_nue_occurrence_counts.get(kmer_info["kmer"], 0)
            kmer_info["occurrence_count"] = occurrence_count
            kmer_info["pct_occurrence"] = (
                (occurrence_count / total_included * 100) if total_included > 0 else 0
            )

        for kmer_info in top_ce_kmers:
            occurrence_count = total_ce_occurrence_counts.get(kmer_info["kmer"], 0)
            kmer_info["occurrence_count"] = occurrence_count
            kmer_info["pct_occurrence"] = (
                (occurrence_count / total_included * 100) if total_included > 0 else 0
            )

        os.makedirs(args.results_dir, exist_ok=True)

        taxon = getattr(args, "taxon", None)
        extraction_stats = getattr(args, "extraction_stats", None)

        report_metadata = {
            "taxon": taxon,
            "generated_at": datetime.datetime.now(),
            "kmer_size": args.kmer_size,
            "min_3utr_length": args.min_3utr_length,
            "num_downstream_nt": args.num_downstream_nt,
            "step_size": args.step_size,
            "total_terminators": total_terminators,
            "skipped_terminators": total_skipped,
            "included_terminators": total_terminators - total_skipped,
            "num_accessions": len(accession_info),
            "accession_info": accession_info,
            "extraction_stats": extraction_stats,
        }

        _save_overview_report(
            os.path.join(args.results_dir, "analysis_overview.txt"),
            report_metadata,
        )
        _save_region_report(
            os.path.join(args.results_dir, "CE_report.txt"),
            "CE",
            top_ce_kmers,
            (CE_ANALYSIS_WINDOW_START, CE_ANALYSIS_WINDOW_END),
            args.kmer_size,
        )
        _save_region_report(
            os.path.join(args.results_dir, "NUE_report.txt"),
            "NUE",
            top_nue_kmers,
            (NUE_ANALYSIS_WINDOW_START, NUE_ANALYSIS_WINDOW_END),
            args.kmer_size,
        )

        nue_plot_path = os.path.join(args.results_dir, "NUE" + PLOT_FILE_SUFFIX)
        ce_plot_path = os.path.join(args.results_dir, "CE" + PLOT_FILE_SUFFIX)

        plot_signal_distribution(
            top_nue_kmers,
            total_nue_counts,
            "NUE",
            PLOT_NUE_X_MIN,
            PLOT_NUE_X_MAX,
            nue_plot_path,
        )
        plot_signal_distribution(
            top_ce_kmers,
            total_ce_counts,
            "CE",
            PLOT_CE_X_MIN,
            PLOT_CE_X_MAX,
            ce_plot_path,
        )

    except Exception:
        print(f"ERROR: {traceback.format_exc()}", file=sys.stderr)


def add_analyse_args(
    parser: argparse.ArgumentParser, is_standalone: bool = False
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
            "--num-downstream-nt",
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

    add_analyse_args(parser, is_standalone=True)

    args = parser.parse_args()

    if not os.path.exists(args.input_path):
        parser.error(f"'{args.input_path}' does not exist.")

    if not args.min_3utr_length >= abs(NUE_ANALYSIS_WINDOW_START):
        parser.error(
            f"Minimum 3'UTR length must be at least {abs(NUE_ANALYSIS_WINDOW_START)}."
        )

    if not args.num_downstream_nt >= CE_ANALYSIS_WINDOW_END:
        parser.error(f"Downstream nts must be at least {CE_ANALYSIS_WINDOW_END}.")

    return args


def _count_kmers(
    sequences: list[str],
    analysis_window: tuple[int, int],
    kmer_size: int,
    num_downstream_nt: int,
    step_size: int,
    expected_region: tuple[int, int],
) -> tuple[dict, dict]:
    """
    Counts k-mers at each position within the specified region across all sequences.
    If region_start and region_end are 0, counts k-mers across the whole sequence.
    Also counts binary presence (max 1 per sequence per kmer) of k-mers within the expected window.

    Args:
        sequences: List of sequences to analyse.
        region_start: Start position of the region.
        region_end: End position of the region.
        kmer_size: Size of the k-mers.
        step_size: Step size for k-mer counting.

    Returns:
        Tuple of (positional_counts, presence_counts):
            - positional_counts: Dictionary of k-mer counts by position { kmer: { pos1: count, pos2: count } }.
            - presence_counts: Dictionary of k-mer presence counts { kmer: count }
    """

    window_start, window_end = analysis_window
    region_start, region_end = expected_region

    assert isinstance(sequences, list) and all(
        isinstance(s, str) for s in sequences
    ), "Sequences must be a list of strings."
    assert isinstance(window_start, int) and isinstance(
        window_end, int
    ), "Region start and end must be integers."
    assert (
        window_start <= window_end
    ), "Region start must be less than or equal to region end."
    assert (
        isinstance(kmer_size, int) and kmer_size > 0
    ), "K-mer size must be a positive integer."
    assert (
        isinstance(step_size, int) and step_size > 0
    ), "Step size must be a positive integer."

    positional_counts = {}
    presence_counts = {}
    is_global = window_start == 0 and window_end == 0

    for seq in sequences:
        seq_len = len(seq)
        utr_len = seq_len - num_downstream_nt
        seen_in_seq = set()

        for i in range(0, seq_len - kmer_size + 1, step_size):
            pos = (i + kmer_size - 1) - utr_len  # anchored to rightmost nt of kmer
            if pos >= 0:
                pos += 1  # +1 is the first nt of downstream region

            if is_global or (window_start <= pos <= window_end):
                kmer = seq[i : i + kmer_size]

                # positional counts (analysis window)
                if kmer not in positional_counts:
                    positional_counts[kmer] = {}
                if pos not in positional_counts[kmer]:
                    positional_counts[kmer][pos] = 0
                positional_counts[kmer][pos] += 1

                # presence counts (expected region)
                if region_start <= pos <= region_end:
                    if kmer not in seen_in_seq:
                        seen_in_seq.add(kmer)
                        if kmer not in presence_counts:
                            presence_counts[kmer] = 0
                        presence_counts[kmer] += 1

    return positional_counts, presence_counts


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

        seq = str(record)
        total_len = len(seq)
        utr_len = total_len - args.num_downstream_nt

        if utr_len < args.min_3utr_length:
            num_skipped += 1
            continue

        terminators.append(seq)

    nue_counts, nue_occurrence_counts = _count_kmers(
        terminators,
        (NUE_ANALYSIS_WINDOW_START, NUE_ANALYSIS_WINDOW_END),
        args.kmer_size,
        args.num_downstream_nt,
        args.step_size,
        (NUE_EXPECTED_START, NUE_EXPECTED_END),
    )
    ce_counts, ce_occurrence_counts = _count_kmers(
        terminators,
        (CE_ANALYSIS_WINDOW_START, CE_ANALYSIS_WINDOW_END),
        args.kmer_size,
        args.num_downstream_nt,
        args.step_size,
        (CE_EXPECTED_START, CE_EXPECTED_END),
    )

    return (
        nue_counts,
        ce_counts,
        nue_occurrence_counts,
        ce_occurrence_counts,
        num_skipped,
        num_terminators,
        accession,
    )


def _merge_counts(target: dict, src: dict) -> None:
    """
    Merges k-mer counts from src into target.
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
def _get_fasta_files(input_dir: str, accessions: list[str]) -> list[str]:
    assert os.path.isdir(input_dir), f"Input path '{input_dir}' is not a directory."

    candidates = glob.glob(os.path.join(input_dir, "*.fa"))

    results = []
    for fpath in candidates:
        base = os.path.basename(fpath)
        accession = os.path.splitext(base)[0]
        # if filter set provided, include only those accessions
        if accessions and accession not in accessions:
            continue
        results.append(fpath)

    return results


def _save_overview_report(out_fpath: str, metadata: dict) -> None:
    separator_width = 100

    lines = []
    lines.append("=" * separator_width + "\n")
    lines.append("Terminator Analysis Overview\n")
    lines.append("=" * separator_width + "\n\n")

    lines.append("ANALYSIS METADATA\n")
    lines.append("-" * separator_width + "\n")
    if metadata.get("taxon"):
        lines.append(f"Taxon:                  {metadata['taxon']}\n")
    lines.append(
        f"Generated:              {metadata['generated_at'].strftime('%Y-%m-%d %H:%M:%S')}\n"
    )
    lines.append(f"K-mer Size:             {metadata['kmer_size']}\n")
    lines.append(f"Step Size:              {metadata['step_size']}\n")
    lines.append(f"Min 3'UTR Length:       {metadata['min_3utr_length']} nt\n")
    lines.append(f"Num Downstream nt:      {metadata['num_downstream_nt']} nt\n")
    lines.append(f"Total Terminators:      {metadata['total_terminators']:,}\n")
    lines.append(f"Included Terminators:   {metadata['included_terminators']:,}\n")
    lines.append(f"Skipped Terminators:    {metadata['skipped_terminators']:,}\n")
    skip_pct = (
        (metadata["skipped_terminators"] / metadata["total_terminators"] * 100)
        if metadata["total_terminators"] > 0
        else 0
    )
    lines.append(f"Skip Rate:              {skip_pct:.2f}%\n")
    lines.append(f"Number of Accessions:   {metadata['num_accessions']}\n")
    lines.append("\n")

    accession_info = metadata.get("accession_info") or {}
    extraction_stats = metadata.get("extraction_stats") or {}

    if accession_info:
        lines.append("TERMINATOR COUNTS BY ACCESSION\n")
        lines.append("-" * separator_width + "\n")

        header = {
            "accession": "Accession",
            "organism": "Organism name",
            "total": "Total",
            "included": "Included",
            "skipped": "Skipped",
            "pct": "% of Included",
        }
        total_all = metadata["total_terminators"]
        total_included = metadata["included_terminators"]
        total_skipped = metadata["skipped_terminators"]

        max_widths = {
            "accession": len(header["accession"]),
            "organism": len(header["organism"]),
            "total": len(header["total"]),
            "included": len(header["included"]),
            "skipped": len(header["skipped"]),
            "pct": len(header["pct"]),
        }

        for accession, info in accession_info.items():
            max_widths["accession"] = max(max_widths["accession"], len(accession))
            max_widths["organism"] = max(max_widths["organism"], len(info["organism_name"]))

        max_widths["total"] = max(max_widths["total"], len(f"{total_all:,}"))
        max_widths["included"] = max(max_widths["included"], len(f"{total_included:,}"))
        max_widths["skipped"] = max(max_widths["skipped"], len(f"{total_skipped:,}"))

        fmt = (
            f"{{accession:<{max_widths['accession']}}} | "
            f"{{organism:<{max_widths['organism']}}} | "
            f"{{total:>{max_widths['total']}}} | "
            f"{{included:>{max_widths['included']}}} | "
            f"{{skipped:>{max_widths['skipped']}}} | "
            f"{{pct:>{max_widths['pct']}}}"
        )
        sep = (
            f"{'-' * max_widths['accession']} | "
            f"{'-' * max_widths['organism']} | "
            f"{'-' * max_widths['total']} | "
            f"{'-' * max_widths['included']} | "
            f"{'-' * max_widths['skipped']} | "
            f"{'-' * max_widths['pct']}"
        )

        lines.append(fmt.format(**header) + "\n")
        lines.append(sep + "\n")

        sorted_accessions = sorted(
            accession_info.items(), key=lambda x: x[1]["total"], reverse=True
        )
        for accession, info in sorted_accessions:
            pct = (
                (info["included"] / total_included * 100) if total_included > 0 else 0
            )
            lines.append(
                fmt.format(
                    accession=accession,
                    organism=info["organism_name"],
                    total=f"{info['total']:,}",
                    included=f"{info['included']:,}",
                    skipped=f"{info['skipped']:,}",
                    pct=f"{pct:.2f}%",
                )
                + "\n"
            )

        lines.append(sep + "\n")
        lines.append(
            fmt.format(
                accession="TOTAL",
                organism=" " * max_widths["organism"],
                total=f"{total_all:,}",
                included=f"{total_included:,}",
                skipped=f"{total_skipped:,}",
                pct="100.00%",
            )
            + "\n"
        )
        lines.append("\n")

    # Extraction skip reasons per accession (if available)
    if extraction_stats:
        lines.append("TRANSCRIPT COUNTS BY ACCESSION\n")
        lines.append("-" * separator_width + "\n")

        # Collect all reason keys across accessions to build columns
        all_reasons = set()
        for info in extraction_stats.values():
            all_reasons.update(info.get("skip_reasons", {}).keys())
        all_reasons = sorted(all_reasons)

        if not all_reasons:
            lines.append("No transcripts skipped during extraction.\n\n")
        else:
            header = {
                "accession": "Accession",
                "total": "Total",
                "extracted": "Extracted",
                "skipped": "Skipped",
                **{r: r for r in all_reasons},
            }

            # Compute column widths
            maxw = {
                "accession": len(header["accession"]),
                "total": len(header["total"]),
                "extracted": len(header["extracted"]),
                "skipped": len(header["skipped"]),
            }
            for r in all_reasons:
                maxw[r] = len(r)

            # Accumulate totals
            total_extracted = 0
            total_skipped = 0
            total_transcripts = 0
            total_reasons = {r: 0 for r in all_reasons}

            for accession, info in extraction_stats.items():
                num_extracted = info.get("num_extracted", 0)
                skip_reasons = info.get("skip_reasons", {})
                acc_skipped = sum(skip_reasons.values())
                acc_transcripts = num_extracted + acc_skipped

                total_extracted += num_extracted
                total_skipped += acc_skipped
                total_transcripts += acc_transcripts
                for r in all_reasons:
                    total_reasons[r] += skip_reasons.get(r, 0)

                maxw["accession"] = max(maxw["accession"], len(accession))
                maxw["total"] = max(maxw["total"], len(f"{acc_transcripts:,}"))
                maxw["extracted"] = max(maxw["extracted"], len(f"{num_extracted:,}"))
                maxw["skipped"] = max(maxw["skipped"], len(f"{acc_skipped:,}"))
                for r in all_reasons:
                    maxw[r] = max(maxw[r], len(f"{skip_reasons.get(r,0):,}"))

            # Also consider TOTAL row widths
            maxw["accession"] = max(maxw["accession"], len("TOTAL"))
            maxw["total"] = max(maxw["total"], len(f"{total_transcripts:,}"))
            maxw["extracted"] = max(maxw["extracted"], len(f"{total_extracted:,}"))
            maxw["skipped"] = max(maxw["skipped"], len(f"{total_skipped:,}"))
            for r in all_reasons:
                maxw[r] = max(maxw[r], len(f"{total_reasons[r]:,}"))

            # Build format strings
            fmt_parts = [
                f"{{accession:<{maxw['accession']}}}",
                f"{{total:>{maxw['total']}}}",
                f"{{extracted:>{maxw['extracted']}}}",
                f"{{skipped:>{maxw['skipped']}}}",
            ]
            sep_parts = [
                "-" * maxw["accession"],
                "-" * maxw["total"],
                "-" * maxw["extracted"],
                "-" * maxw["skipped"],
            ]
            for r in all_reasons:
                fmt_parts.append(f"{{{r}:>{maxw[r]}}}")
                sep_parts.append("-" * maxw[r])

            fmt = " | ".join(fmt_parts)
            sep = " | ".join(sep_parts)

            # Write header
            lines.append(fmt.format(**header) + "\n")
            lines.append(sep + "\n")

            # Rows per accession (sorted by accession)
            for accession, info in sorted(
                extraction_stats.items(), key=lambda x: x[0]
            ):
                num_extracted = info.get("num_extracted", 0)
                skip_reasons = info.get("skip_reasons", {})
                acc_skipped = sum(skip_reasons.values())
                row = {
                    "accession": accession,
                    "total": f"{num_extracted + acc_skipped:,}",
                    "extracted": f"{num_extracted:,}",
                    "skipped": f"{acc_skipped:,}",
                }
                for r in all_reasons:
                    row[r] = f"{skip_reasons.get(r, 0):,}"
                lines.append(fmt.format(**row) + "\n")

            # Separator and TOTAL row
            lines.append(sep + "\n")
            total_row = {
                "accession": "TOTAL",
                "total": f"{total_transcripts:,}",
                "extracted": f"{total_extracted:,}",
                "skipped": f"{total_skipped:,}",
            }
            for r in all_reasons:
                total_row[r] = f"{total_reasons[r]:,}"
            lines.append(fmt.format(**total_row) + "\n\n")

    lines.append("=" * separator_width + "\n")

    with open(out_fpath, "w") as f:
        f.writelines(lines)

    print(f"Overview report saved to '{out_fpath}'")


def _save_region_report(
    out_fpath: str,
    region_name: str,
    ranked_kmers: list,
    region_window: tuple[int, int],
    kmer_size: int,
) -> None:

    separator_width = 100

    lines = []
    lines.append("=" * separator_width + "\n")
    lines.append(f"{region_name} Region Report\n")
    lines.append("=" * separator_width + "\n\n")

    lines.append(f"Analysis Window:        {region_window[0]} to {region_window[1]}\n")

    # Add expected region info
    if region_name == "NUE":
        lines.append(
            f"Expected Region:        {NUE_EXPECTED_START} to {NUE_EXPECTED_END}\n"
        )
    elif region_name == "CE":
        lines.append(
            f"Expected Region:        {CE_EXPECTED_START} to {CE_EXPECTED_END}\n"
        )

    lines.append("\n")

    lines.append(f"TOP {len(ranked_kmers)} K-MERS\n")
    lines.append("-" * separator_width + "\n")

    if not ranked_kmers:
        lines.append("No k-mers found.\n")
    else:
        header = {
            "kmer": "K-mer",
            "delta": "Delta",
            "peak": "Peak Count",
            "median": "Median Count",
            "pos": "Peak Pos",
            "occurrences": "Occurrence Count",
            "pct": "% Occurrence",
        }

        maxw = {
            "kmer": max(kmer_size, len(header["kmer"])),
            "delta": len(header["delta"]),
            "peak": len(header["peak"]),
            "median": len(header["median"]),
            "pos": len(header["pos"]),
            "occurrences": len(header["occurrences"]),
            "pct": len(header["pct"]),
        }

        for item in ranked_kmers:
            maxw["delta"] = max(maxw["delta"], len(f"{item['delta']:.1f}"))
            maxw["peak"] = max(maxw["peak"], len(f"{item['peak_count']:,}"))
            maxw["median"] = max(maxw["median"], len(f"{item['median_count']:.1f}"))
            maxw["pos"] = max(maxw["pos"], len(str(item["peak_pos"])))
            maxw["occurrences"] = max(maxw["occurrences"], len(f"{item['occurrence_count']:,}"))
            maxw["pct"] = max(maxw["pct"], len(f"{item['pct_occurrence']:.1f}%"))

        fmt = (
            f"{{kmer:<{maxw['kmer']}}} | "
            f"{{delta:>{maxw['delta']}}} | "
            f"{{peak:>{maxw['peak']}}} | "
            f"{{median:>{maxw['median']}}} | "
            f"{{pos:>{maxw['pos']}}} | "
            f"{{occurrences:>{maxw['occurrences']}}} | "
            f"{{pct:>{maxw['pct']}}}"
        )
        sep = (
            f"{'-' * maxw['kmer']} | "
            f"{'-' * maxw['delta']} | "
            f"{'-' * maxw['peak']} | "
            f"{'-' * maxw['median']} | "
            f"{'-' * maxw['pos']} | "
            f"{'-' * maxw['occurrences']} | "
            f"{'-' * maxw['pct']}"
        )

        lines.append(fmt.format(**header) + "\n")
        lines.append(sep + "\n")
        for item in ranked_kmers:
            lines.append(
                fmt.format(
                    kmer=item["kmer"],
                    delta=f"{item['delta']:.1f}",
                    peak=f"{item['peak_count']:,}",
                    median=f"{item['median_count']:.1f}",
                    pos=str(item["peak_pos"]),
                    occurrences=f"{item['occurrence_count']:,}",
                    pct=f"{item['pct_occurrence']:.1f}%",
                )
                + "\n"
            )

    lines.append("\n")
    lines.append("=" * separator_width + "\n")

    with open(out_fpath, "w") as f:
        f.writelines(lines)
    print(f"{region_name} report saved to '{out_fpath}'")


# --- STANDALONE EXECUTION ---
if __name__ == "__main__":
    sys.exit(main())
