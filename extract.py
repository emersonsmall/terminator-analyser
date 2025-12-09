# Built-in libraries
import os
import sys
import argparse
import textwrap
from collections import defaultdict
from multiprocessing import Pool
from functools import partial

# External libraries
import pyfaidx  # https://anaconda.org/bioconda/pyfaidx
import gffutils  # https://anaconda.org/bioconda/gffutils

# Internal modules
from get_genomes import VALID_FASTA_EXTS, VALID_ANNOTATION_EXTS

# UTR refers to 3'UTR unless otherwise specified

TERMINATOR_FILE_SUFFIX = "_terminators.fa"


def run_extraction(args: argparse.Namespace) -> int:
    try:
        file_pairs = _find_files(args.input_path)

        # process each genome in parallel
        worker = partial(_extract_all_terminators, args=args)
        with Pool() as pool:
            pool.map(worker, file_pairs)

        return 0
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


def add_extract_args(
    parser: argparse.ArgumentParser, is_standalone: bool = True
) -> None:
    """Adds command-line arguments for the `extract` command to the given parser.

    Args:
        parser: The argument parser to which the arguments will be added.
        standalone: Whether to include standalone execution arguments. Defaults to True.
    """

    parser.add_argument(
        "-r",
        "--raw-dna",
        action="store_true",
        help="Output raw DNA sequences instead of transcribed RNA.",
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
        "-d",
        "--downstream-nts",
        type=int,
        default=50,
        help="Number of nucleotides downstream of the CS to extract (default: 50).",
    )

    if is_standalone:
        parser.add_argument(
            "input_path",
            help="Path to the input folder containing FASTA and GFF files.",
        )
        parser.add_argument(
            "-o",
            "--output-dir",
            default=os.path.join("out", "terminators"),
            help="Path to the output directory (default: './out/terminators').",
        )


def _get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract terminators for the given genomes/genus."
    )
    add_extract_args(parser)

    args = parser.parse_args()
    if not os.path.isdir(args.input_path):
        parser.error(
            f"Input directory '{args.input_path}' does not exist or is not a directory."
        )
    return args


def _is_internal_priming_artifact(
    downstream_sequence: str, consecutive_a: int, total_a: int, window_size: int
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
        bool: True if the sequence is an internal priming artifact, False otherwise.
    """

    assert (
        len(downstream_sequence) >= window_size
    ), "Downstream sequence length must be greater than or equal to window size."
    assert consecutive_a >= 0, "Filter for consecutive A's must be non-negative."
    assert total_a >= 0, "Filter for total A's in window must be non-negative."
    assert window_size > 0, "Window size must be greater than 0."

    if consecutive_a == 0 and total_a == 0:
        return False

    region_to_check = downstream_sequence[:window_size].upper()

    if consecutive_a > 0 and "A" * consecutive_a in region_to_check:
        return True

    if total_a > 0 and region_to_check.count("A") >= total_a:
        return True

    return False


def _extract_terminator(
    tscript: gffutils.Feature,
    fasta: pyfaidx.Fasta,
    cds_features: list,
    exon_features: list,
    downstream_nts: int,
) -> str | None:
    """
    Extracts the terminator sequence for a given transcript feature.
    Handles strand sense and reverse complementing.
    """

    # All coordinates 1-based until modified in slice operations
    # Python slices are 0-based, end-exclusive -> subtract 1 only for start coords
    utr_parts = []

    if tscript.strand == "+":
        cds_end = cds_features[-1].end

        # Get all exons after the CDS end
        for exon in exon_features:
            if exon.end > cds_end:
                utr_start = max(
                    exon.start, cds_end + 1
                )  # One exon can contain CDS and 3'UTR content
                utr_end = exon.end
                utr_parts.append(fasta[tscript.chrom][utr_start - 1 : utr_end].seq)

        downstream_start = tscript.end + 1
        downstream_end = tscript.end + downstream_nts

    elif tscript.strand == "-":
        cds_start = cds_features[0].start

        # Get all exons before the CDS start
        for exon in exon_features:
            if exon.start < cds_start:
                utr_start = exon.start
                utr_end = min(exon.end, cds_start - 1)
                utr_parts.append(fasta[tscript.chrom][utr_start - 1 : utr_end].seq)

        downstream_start = max(1, tscript.start - downstream_nts)
        downstream_end = tscript.start - 1

    else:
        print(
            f"WARNING: unknown strand '{tscript.strand}' for feature '{tscript.id}'",
            file=sys.stderr,
        )
        return None

    if not utr_parts:
        return None

    full_utr = "".join(utr_parts)
    downstream_seq = ""
    if downstream_start <= downstream_end:
        downstream_seq = fasta[tscript.chrom][downstream_start - 1 : downstream_end].seq

    # assemble final sequence
    if tscript.strand == "+":
        term_seq = full_utr + downstream_seq  # 5' to 3'
    else:
        term_seq = downstream_seq + full_utr
        term_seq = pyfaidx.Sequence(seq=term_seq).reverse.complement.seq

    return term_seq


def _filter_sequence(term_seq: str, args: argparse.Namespace) -> bool:
    """
    Applies filters to the terminator sequence.

    Returns:
        bool: True if the sequence passes the filters, False otherwise.
    """

    if args.raw_dna:
        return True

    final_downstream_seq = term_seq[-args.downstream_nts :]
    if len(final_downstream_seq) < args.filter_window_size:
        return False

    if _is_internal_priming_artifact(
        final_downstream_seq,
        args.filter_consecutive_a,
        args.filter_window_a,
        args.filter_window_size,
    ):
        return False

    return True


def _process_transcript(
    tscript: gffutils.Feature,
    db: gffutils.FeatureDB,
    fasta: pyfaidx.Fasta,
    args: argparse.Namespace,
) -> str | None:
    """Processes a single transcript feature to extract and format its terminator sequence."""

    try:
        cds_features = list(db.children(tscript, featuretype="CDS", order_by="start"))
        exon_features = list(db.children(tscript, featuretype="exon", order_by="start"))

        # some GFF formats have CDS and exon features as children of gene, not transcript
        if not cds_features or not exon_features:
            try:
                parent_gene = list(db.parents(tscript))[0]
                if not cds_features:
                    cds_features = list(db.children(parent_gene, featuretype="CDS", order_by="start"))
                if not exon_features:
                    exon_features = list(db.children(parent_gene, featuretype="exon", order_by="start"))
            except (IndexError, StopIteration):
                pass

        if not cds_features or not exon_features:
            return None


        term_seq = _extract_terminator(
            tscript, fasta, cds_features, exon_features, args.downstream_nts
        )

        if not term_seq or not _filter_sequence(term_seq, args):
            return None

        return _format_fasta_record(tscript, term_seq, args.raw_dna)

    except Exception as e:
        print(
            f"WARNING: could not process feature '{tscript.id}': {e}", file=sys.stderr
        )
        return None


def _extract_all_terminators(file_pair: tuple, args: argparse.Namespace) -> None:
    """
    Extracts terminator sequences from the given fasta and gff files.

    Args:
        file_pair: A tuple containing (fasta_fpath, annotation_fpath).
        args: Parsed command-line arguments.
    """

    fasta_fpath, annotation_fpath = file_pair

    assert os.path.isfile(fasta_fpath), f"Fasta file '{fasta_fpath}' does not exist."
    assert os.path.isfile(
        annotation_fpath
    ), f"Annotation file '{annotation_fpath}' does not exist."

    base = os.path.basename(fasta_fpath)
    fname = os.path.splitext(base)[0]
    print(f"Processing genome '{fname}'")

    db_dir = os.path.join(
        "out", "FeatureDBs"  # TODO: make this adapt to diff output folders
    )  # Allows reuse of DBs between different taxons
    os.makedirs(db_dir, exist_ok=True)
    os.makedirs(args.output_dir, exist_ok=True)

    db_fpath = os.path.join(db_dir, f"{fname}.db")

    terminators_fpath = os.path.join(
        args.output_dir, f"{fname}{TERMINATOR_FILE_SUFFIX}"
    )

    fasta = pyfaidx.Fasta(fasta_fpath)
    db = _create_gff_db(annotation_fpath, db_fpath)

    transcript_labels = ("mRNA", "transcript")  # GFF uses "mRNA", GTF uses "transcript"
    transcripts = []
    for t_label in transcript_labels:
        try:
            if db.count_features_of_type(t_label) > 0:
                transcripts = db.features_of_type(t_label, order_by="start")
                break
        except:
            continue

    if not transcripts:
        print(
            f"WARNING: No transcripts found with labels {transcript_labels} in '{annotation_fpath}', skipping",
            file=sys.stderr,
        )
        return

    out_records = []
    num_skipped = 0
    for tscript in transcripts:
        record = _process_transcript(tscript, db, fasta, args)
        if record:
            out_records.append(record)
        else:
            num_skipped += 1

    with open(terminators_fpath, "w") as out_f:
        out_f.writelines(out_records)

    print(f"Extracted {len(out_records)} terminator sequences to '{terminators_fpath}'")
    print(f"skipped {num_skipped} transcripts")


# --- HELPER FUNCTIONS ---
def _format_fasta_record(
    tscript: gffutils.Feature, term_seq: str, raw_dna: bool
) -> str:
    """
    Formats the given terminator sequence into a FASTA record header.

    Returns:
        The formatted FASTA record.
    """

    display_id = tscript.id
    if "orig_protein_id" in tscript.attributes:
        raw_id = tscript.attributes["orig_protein_id"][0]
        display_id = raw_id.split("|")[-1]

    header = f">{display_id} | {tscript.chrom}:{tscript.start}-{tscript.end}({tscript.strand})"

    term_seq = term_seq.upper()
    if not raw_dna:
        term_seq = term_seq.replace("T", "U")
    wrapped_seq = textwrap.fill(term_seq, width=80)

    return f"{header}\n{wrapped_seq}\n"


def _find_files(dir: str) -> list[tuple[str, str]]:
    """
    Searches the given directory for matching pairs of FASTA and GFF files.

    Args:
        dir: The path to the directory to search.

    Returns:
        A list of tuples, where each tuple is a pair of FASTA and GFF file paths.
    """

    assert os.path.isdir(dir), f"Folder '{dir}' does not exist or is not a directory."

    files_by_basename = defaultdict(dict)

    for filename in os.listdir(dir):
        basename, ext = os.path.splitext(filename)
        ext = ext.lower()

        if ext in VALID_FASTA_EXTS:
            files_by_basename[basename]["fasta"] = os.path.join(dir, filename)
        elif ext in VALID_ANNOTATION_EXTS:
            files_by_basename[basename]["gff"] = os.path.join(dir, filename)

    file_pairs = []
    for basename, paths in files_by_basename.items():
        if "fasta" in paths and "gff" in paths:
            file_pairs.append((paths["fasta"], paths["gff"]))
        else:
            print(
                f"WARNING: Incomplete file pair for '{basename}'. Skipping. Found: {list(paths.keys())}",
                file=sys.stderr,
            )

    if not file_pairs:
        raise Exception(
            f"No valid FASTA and GFF file pairs found in directory '{dir}'."
        )

    return file_pairs


def _create_gff_db(annotation_fpath: str, db_fpath: str) -> gffutils.FeatureDB:
    """
    Creates a gffutils FeatureDB from the given annotation file if it does not already exist.

    Args:
        annotation_fpath: Filepath to the input annotation file.
        db_fpath: Filepath to the output database file.

    Returns:
        A gffutils.FeatureDB object.
    """

    assert os.path.isfile(annotation_fpath), f"'{annotation_fpath}' does not exist."

    if not os.path.isfile(db_fpath):
        print(f"Creating FeatureDB at '{db_fpath}'")
        gffutils.create_db(
            annotation_fpath,
            dbfn=db_fpath,
            keep_order=True,
            merge_strategy="create_unique",
            sort_attribute_values=True,
            disable_infer_transcripts=True,
            disable_infer_genes=True,
        )
        print(f"FeatureDB created at '{db_fpath}'")
    else:
        print(f"Using existing FeatureDB at '{db_fpath}'")

    return gffutils.FeatureDB(db_fpath)


# --- STANDALONE EXECUTION ---
def main() -> int:
    return run_extraction(_get_args())


if __name__ == "__main__":
    sys.exit(main())
