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
from constants import VALID_FASTA_EXTS, VALID_ANNOTATION_EXTS

# utr refers to 3'UTR unless otherwise specified


class NoUTRException(Exception):
    """Raised when a transcript has no implied 3' UTR based on CDS and exon coordinates."""

    pass


class InvalidStrandException(Exception):
    """Raised when a transcript has an invalid or unknown strand."""

    pass


class NoCDSException(Exception):
    """Raised when a transcript has no CDS features."""

    pass


class NoExonException(Exception):
    """Raised when a transcript has no exon features."""

    pass


def main() -> None:
    run_extraction(_get_args())


def run_extraction(args: argparse.Namespace) -> None:
    try:
        genomes = getattr(args, "genomes", {})
        file_pairs = _get_file_pairs(args.input_dir, list(genomes.keys()))
        extraction_stats = {}

        # process each genome in parallel
        worker = partial(_extract_all_terminators, args=args)
        with Pool() as pool:
            for result in pool.imap_unordered(worker, file_pairs):
                if result:
                    extraction_stats[result["accession"]] = {
                        "num_extracted": result["num_extracted"],
                        "skip_reasons": result["skip_reasons"],
                    }

        args.extraction_stats = extraction_stats  # store in args for downstream use

    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)


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
        "-d",
        "--num-downstream-nt",
        type=int,
        default=50,
        help="Number of nucleotides downstream of the CS to extract (default: 50).",
    )
    parser.add_argument(
        "--terminators-dir",
        default=os.path.join("out", "terminators"),
        help="Path to the output directory (default: './out/terminators').",
    )

    if is_standalone:
        parser.add_argument(
            "input_dir",
            help="Path to the input folder containing FASTA and GFF files.",
        )


def validate_extract_args(args: argparse.Namespace) -> None:
    if args.num_downstream_nt < 0:
        raise ValueError("Number of downstream nucleotides must be non-negative.")


def _get_args() -> argparse.Namespace:
    """Parses and validates command-line arguments for standalone execution.

    Returns:
        Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Extract terminators for the given set of genomes."
    )

    add_extract_args(parser)
    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        parser.error(
            f"Input directory '{args.input_dir}' does not exist or is not a directory."
        )

    try:
        validate_extract_args(args)
    except ValueError as e:
        parser.error(str(e))

    return args


def _extract_terminator(
    tscript: gffutils.Feature,
    fasta: pyfaidx.Fasta,
    cds_features: list,
    exon_features: list,
    num_downstream_nt: int,
) -> str:
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
        downstream_end = tscript.end + num_downstream_nt

    elif tscript.strand == "-":
        cds_start = cds_features[0].start

        # Get all exons before the CDS start
        for exon in exon_features:
            if exon.start < cds_start:
                utr_start = exon.start
                utr_end = min(exon.end, cds_start - 1)
                utr_parts.append(fasta[tscript.chrom][utr_start - 1 : utr_end].seq)

        downstream_start = max(1, tscript.start - num_downstream_nt)
        downstream_end = tscript.start - 1

    else:
        raise InvalidStrandException(
            f"Transcript '{tscript.id}' has invalid strand '{tscript.strand}'."
        )

    if not utr_parts:
        raise NoUTRException(f"Transcript '{tscript.id}' has no implied 3' UTR.")

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

    return term_seq.upper()


def _process_transcript(
    tscript: gffutils.Feature,
    db: gffutils.FeatureDB,
    fasta: pyfaidx.Fasta,
    args: argparse.Namespace,
) -> str | None:
    """Processes a single transcript feature to extract and format its terminator sequence."""

    cds_features = list(db.children(tscript, featuretype="CDS", order_by="start"))
    exon_features = list(db.children(tscript, featuretype="exon", order_by="start"))

    # some GFF formats have CDS and exon features as children of gene, not transcript
    if not cds_features or not exon_features:
        try:
            parent_gene = list(db.parents(tscript))[0]
            if not cds_features:
                cds_features = list(
                    db.children(parent_gene, featuretype="CDS", order_by="start")
                )
            if not exon_features:
                exon_features = list(
                    db.children(parent_gene, featuretype="exon", order_by="start")
                )
        except (IndexError, StopIteration):
            pass

    if not exon_features:
        raise NoExonException(f"Transcript '{tscript.id}' has no exon features.")
    if not cds_features:
        raise NoCDSException(f"Transcript '{tscript.id}' has no CDS features.")

    term_seq = _extract_terminator(
        tscript, fasta, cds_features, exon_features, args.num_downstream_nt
    )

    return _format_fasta_record(tscript, term_seq, args.raw_dna)


def _extract_all_terminators(file_pair: tuple, args: argparse.Namespace) -> dict | None:
    """
    Extracts terminator sequences from the given fasta and gff files.

    Args:
        file_pair: A tuple containing (fasta_fpath, annotation_fpath).
        args: Parsed command-line arguments.
    """

    fasta_fpath, annotation_fpath = file_pair

    base = os.path.basename(fasta_fpath)
    accession = os.path.splitext(base)[0]

    os.makedirs(args.terminators_dir, exist_ok=True)
    terminators_fpath = os.path.join(args.terminators_dir, f"{accession}.fa")

    # skip if terminator fasta already exists
    force = getattr(args, "force", False)
    if not force and os.path.isfile(terminators_fpath):
        print(
            f"{accession} terminators already exist at '{terminators_fpath}', skipping"
        )
        return None

    print(f"Processing genome {accession}")

    db_dir = os.path.join(
        "out", "FeatureDBs"
    )  # Allows reuse of DBs between different taxons
    os.makedirs(db_dir, exist_ok=True)

    db_fpath = os.path.join(db_dir, f"{accession}.db")

    fasta = pyfaidx.Fasta(fasta_fpath)
    db = _create_gff_db(annotation_fpath, db_fpath)

    transcript_labels = ("mRNA", "transcript")  # GFF has "mRNA", GTF has "transcript"
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
        return None

    out_records = []
    skip_reasons = {}

    for tscript in transcripts:
        try:
            record = _process_transcript(tscript, db, fasta, args)
            out_records.append(record)
        except NoUTRException:
            skip_reasons["no_utr"] = skip_reasons.get("no_utr", 0) + 1
        except InvalidStrandException:
            skip_reasons["invalid_strand"] = skip_reasons.get("invalid_strand", 0) + 1
        except NoCDSException:
            skip_reasons["no_cds"] = skip_reasons.get("no_cds", 0) + 1
        except NoExonException:
            skip_reasons["no_exon"] = skip_reasons.get("no_exon", 0) + 1
        except Exception as e:
            print(
                f"WARNING: Unexpected error processing '{tscript.id}': {e}",
                file=sys.stderr,
            )
            skip_reasons["unexpected_error"] = (
                skip_reasons.get("unexpected_error", 0) + 1
            )

    with open(terminators_fpath, "w") as out_f:
        out_f.writelines(out_records)

    print(f"Extracted {len(out_records)} terminator sequences to '{terminators_fpath}'")

    return {
        "accession": accession,
        "num_extracted": len(out_records),
        "skip_reasons": skip_reasons,
    }


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

    if not raw_dna:
        term_seq = term_seq.replace("T", "U")
    wrapped_seq = textwrap.fill(term_seq, width=80)

    return f"{header}\n{wrapped_seq}\n"


def _get_file_pairs(dir: str, accessions: list[str]) -> list[tuple[str, str]]:
    """
    Searches the given directory for matching pairs of FASTA and GFF files.

    Args:
        dir: The path to the directory to search.

    Returns:
        A list of tuples, where each tuple is a pair of FASTA and GFF file paths.
    """

    files_by_basename = defaultdict(dict)

    for filename in os.listdir(dir):
        basename, ext = os.path.splitext(filename)

        # if 'full' execution, only add files in accessions list. Otherwise, include all valid files in input directory
        if accessions and basename not in accessions:
            continue

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
if __name__ == "__main__":
    main()
