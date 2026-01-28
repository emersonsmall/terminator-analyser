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
from reports import save_overview_report, save_region_report
from constants import CE_EXPECTED_END, CE_EXPECTED_START, NUE_EXPECTED_END, NUE_EXPECTED_START

# Coordinates: -1 is the last nt of the 3'UTR, +1 is the first nt of the downstream region

# Analysis window coordinates
NUE_ANALYSIS_WINDOW_START = -50
NUE_ANALYSIS_WINDOW_END = -1
CE_ANALYSIS_WINDOW_START = (
    -10
)  # Loke et al. uses -15, but NUE signals dominate if using -15
CE_ANALYSIS_WINDOW_END = 20

PLOT_NUE_X_MIN = -35
PLOT_NUE_X_MAX = -5
PLOT_CE_X_MIN = -10
PLOT_CE_X_MAX = 15


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

        save_overview_report(
            os.path.join(args.results_dir, "analysis_overview.txt"),
            report_metadata,
        )
        save_region_report(
            os.path.join(args.results_dir, "CE_report.txt"),
            "CE",
            top_ce_kmers,
            (CE_ANALYSIS_WINDOW_START, CE_ANALYSIS_WINDOW_END),
            args.kmer_size,
        )
        save_region_report(
            os.path.join(args.results_dir, "NUE_report.txt"),
            "NUE",
            top_nue_kmers,
            (NUE_ANALYSIS_WINDOW_START, NUE_ANALYSIS_WINDOW_END),
            args.kmer_size,
        )

        plot_signal_distribution(
            top_nue_kmers,
            total_nue_counts,
            "NUE",
            PLOT_NUE_X_MIN,
            PLOT_NUE_X_MAX,
            args.results_dir,
        )
        plot_signal_distribution(
            top_ce_kmers,
            total_ce_counts,
            "CE",
            PLOT_CE_X_MIN,
            PLOT_CE_X_MAX,
            args.results_dir,
        )

    except Exception:
        print(f"ERROR: {traceback.format_exc()}", file=sys.stderr)


def _window_coords_arg(val: str) -> tuple[tuple[int, int], ...]:
    coord_pair_strs = [v.strip() for v in val.split(";") if v.strip()]
    coord_pairs = []

    for coord_pair in coord_pair_strs:
        try:
            start_str, end_str = coord_pair.split(",")
            start = int(start_str)
            end = int(end_str)
            if start >= end and start != 0 and end != 0:
                raise ValueError
            coord_pairs.append((start, end))
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"invalid coordinate pair '{coord_pair}', must be in format '<start>,<end>' with start < end"
            )

    if not coord_pairs:
        raise argparse.ArgumentTypeError("requires at least one coordinate pair")
    
    return tuple(coord_pairs)


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
    parser.add_argument(
        "--filter-consecutive-a",
        type=int,
        default=6,
        help="Filter out terminators with this many consecutive 'A's in the first <filter-window-size> downstream nts (default: 6). Set to 0 to disable.",
    )
    parser.add_argument(
        "--filter-window-a",
        type=int,
        default=8,
        help="Filter out terminators with this many 'A's in the first <filter-window-size> downstream nts (default: 8). Set to 0 to disable.",
    )
    parser.add_argument(
        "--filter-window-size",
        type=int,
        default=10,
        help="Number of downstream nts to check for internal priming artifacts (default: 10).",
    )
    parser.add_argument(
        "--analysis-window-coords",
        type=_window_coords_arg,
        default=None,
        help="Semicolon-separated list of start,end coordinate pairs defining analysis windows (default: ((-50, -1), (-10, 20))).",
    )
    parser.add_argument(
        "--expected-window-coords",
        type=_window_coords_arg,
        default=None,
        help="Semicolon-separated list of start,end coordinate pairs defining expected regions for analysis windows (default: ((-30, -13), (-10, 10))).",
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


def validate_analyse_args(args: argparse.Namespace) -> None:
    """
    Validates the command-line arguments for the `analyse` command.
    """
    if args.min_3utr_length < abs(NUE_ANALYSIS_WINDOW_START):
        raise ValueError(
            f"Minimum 3'UTR length must be at least {abs(NUE_ANALYSIS_WINDOW_START)}."
        )

    if args.num_downstream_nt < CE_ANALYSIS_WINDOW_END:
        raise ValueError(f"Downstream nts must be at least {CE_ANALYSIS_WINDOW_END}.")

    has_analysis = args.analysis_window_coords is not None
    has_expected = args.expected_window_coords is not None

    # both or neither must be provided
    if has_analysis != has_expected:
        raise ValueError("Both --analysis-window-coords and --expected-window-coords must be provided together.")

    # apply defaults if neither provided
    args.analysis_window_coords = (
        (NUE_ANALYSIS_WINDOW_START, NUE_ANALYSIS_WINDOW_END),
        (CE_ANALYSIS_WINDOW_START, CE_ANALYSIS_WINDOW_END),
    )
    args.expected_window_coords = (
        (NUE_EXPECTED_START, NUE_EXPECTED_END),
        (CE_EXPECTED_START, CE_EXPECTED_END),
    )


def _get_args() -> argparse.Namespace:
    """Parses and validates command-line arguments for standalone execution.

    Returns:
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Analyses the specified regions of the given terminator sequences."
    )

    add_analyse_args(parser, is_standalone=True)
    args = parser.parse_args()

    if not os.path.exists(args.input_path):
        parser.error(f"'{args.input_path}' does not exist.")

    try:
        validate_analyse_args(args)
    except ValueError as e:
        parser.error(str(e))

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
    Also counts binary occurrences (max 1 per sequence per kmer) of k-mers within the expected window.

    Args:
        sequences: List of sequences to analyse.
        region_start: Start position of the region.
        region_end: End position of the region.
        kmer_size: Size of the k-mers.
        step_size: Step size for k-mer counting.

    Returns:
        Tuple of (positional_counts, occurrence_counts):
            - positional_counts: Dictionary of k-mer counts by position { kmer: { pos1: count, pos2: count } }.
            - occurrence_counts: Dictionary of k-mer occurrence counts { kmer: count }
    """

    window_start, window_end = analysis_window
    region_start, region_end = expected_region

    positional_counts = {}
    occurrence_counts = {}
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
                        if kmer not in occurrence_counts:
                            occurrence_counts[kmer] = 0
                        occurrence_counts[kmer] += 1

    return positional_counts, occurrence_counts


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

        if not _terminator_passes_filter(seq, args):
            num_skipped += 1
            continue

        terminators.append(seq)

    # TODO: for window in args.analysis_window_coords, expected_region in args.expected_window_coords etc

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


def _is_internal_priming_artifact(
    downstream_seq: str, consecutive_a: int, total_a: int, window_size: int
) -> bool:
    """
    Checks if the given downstream sequence indicates an internal priming artifact.
    Methodology based on Beaudong et al. DOI: 10.1101/gr.10.7.1001

    Args:
        downstream_sequence: The downstream sequence to check.
        consecutive_a: The number of consecutive 'A's to classify sequence as an artifact.
        total_a: The total number of 'A's in the window to classify sequence as an artifact.
        window_size: The number of nucleotides to check from the start of the sequence.

    Returns:
        True if the sequence is an internal priming artifact, False otherwise.
    """

    if consecutive_a == 0 and total_a == 0:
        return False

    region_to_check = downstream_seq[:window_size]

    if consecutive_a > 0 and "A" * consecutive_a in region_to_check:
        return True

    if total_a > 0 and region_to_check.count("A") >= total_a:
        return True

    return False


def _terminator_passes_filter(term_seq: str, args: argparse.Namespace) -> bool:
    """
    Applies filters to the terminator sequence.

    Returns:
        True, if the terminator sequence passes the filters, False otherwise.
    """

    if len(term_seq) - args.num_downstream_nt < args.min_3utr_length:
        return False

    downstream_seq = term_seq[-args.num_downstream_nt:]
    if len(downstream_seq) < args.filter_window_size:
        return False

    if _is_internal_priming_artifact(
        downstream_seq,
        args.filter_consecutive_a,
        args.filter_window_a,
        args.filter_window_size,
    ):
        return False

    return True


# --- STANDALONE EXECUTION ---
if __name__ == "__main__":
    sys.exit(main())
