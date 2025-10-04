import os
import sys
import concurrent.futures
import argparse
import textwrap
from collections import defaultdict

# External libraries
import pyfaidx  # https://anaconda.org/bioconda/pyfaidx
import gffutils # https://anaconda.org/bioconda/gffutils

# UTR refers to 3'UTR unless otherwise specified


def main() -> int:
    return run_extraction(_get_args())

if __name__ == "__main__":
    sys.exit(main())


def run_extraction(args: argparse.Namespace) -> int:
    try:
        # find all file pairs in the input directory
        file_pairs = _find_files(args.input_path)
        tasks = [(file_pair[0], file_pair[1], args) for file_pair in file_pairs]

        # process each genome in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
            list(executor.map(_worker, tasks))

        return 0
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


def _worker(args: tuple) -> None:
    """
    Extracts the terminator sequences of a single genome.

    Args:
        args (tuple): A tuple containing (fasta_fpath, gff_fpath, cli_args).
    """
    fasta_fpath, gff_fpath, cli_args = args
    _extract_all_terminators(fasta_fpath, gff_fpath, cli_args)


def add_args_to_parser(parser: argparse.ArgumentParser, standalone: bool = True) -> None:
    parser.add_argument(
        "-r",
        "--raw-dna",
        action="store_true",
        help="Output raw DNA sequences instead of transcribed RNA."
    )
    parser.add_argument(
        "--filter-consecutive-a",
        type=int,
        default=6,
        help="Filter out terminators with this many consecutive 'A's in the first <filter-window-size> downstream nts (default: 6). Set to 0 to disable."
    )
    parser.add_argument(
        "--filter-window-a",
        type=int,
        default=8,
        help="Filter out terminators with this many 'A's in the first <filter-window-size> downstream nts (default: 8). Set to 0 to disable."
    )
    parser.add_argument(
        "--filter-window-size",
        type=int,
        default=10,
        help="Number of downstream nts to check for internal priming artifacts (default: 10)."
    )
    parser.add_argument(
        "-d",
        "--downstream-nts",
        type=int,
        default=50,
        help="Number of nucleotides downstream of the CS to extract (default: 50)."
    )

    if standalone:
        parser.add_argument(
            "input_path",
            help="Path to the input folder containing FASTA and GFF files."
        )
        parser.add_argument(
            "-o",
            "--output-dir",
            default=os.path.join("out", "terminators"),
            help="Path to the output directory (default: './out/terminators')."
        )

def _get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract terminators for the given genomes/genus."
    )
    add_args_to_parser(parser)

    args = parser.parse_args()
    if not os.path.isdir(args.input_path):
        parser.error(f"Input directory '{args.input_path}' does not exist or is not a directory.")
    return args


def _is_internal_priming_artifact(sequence: str, consecutive_a: int = 6, total_a: int = 8, window_size: int = 10) -> bool:
    """
    Checks if the given downstream sequence is likely to be an internal priming artifact.
    Methodology based on Beaudong et al. DOI: 10.1101/gr.10.7.1001

    Args:
        sequence (str): The downstream sequence to check.
        consecutive_a (int): The number of consecutive 'A's to consider the sequence an artifact (default: 6).
        window_a (int): The total number of 'A's in the window to consider the sequence an artifact (default: 8).
        window_size (int): The number of nucleotides to check from the start of the sequence (default: 10).

    Returns:
        bool: True if the sequence is an internal priming artifact, False otherwise.
    """
    assert len(sequence) >= window_size, "Downstream sequence length must be greater than or equal to window size."
    assert consecutive_a >= 0, "Filter for consecutive A's must be non-negative."
    assert total_a >= 0, "Filter for total A's in window must be non-negative."
    assert window_size > 0, "Window size must be greater than 0."

    if consecutive_a == 0 and total_a == 0:
        return False
    
    region_to_check = sequence[:window_size].upper()

    if consecutive_a > 0 and 'A' * consecutive_a in region_to_check:
        return True
    
    if total_a > 0 and region_to_check.count('A') >= total_a:
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
    Extracts the terminator sequence (3'UTR + downstream region) for a given transcript feature.
    Handles strand and reverse complementing as needed.
    """
    # All coordinates 1-based until modified in slice operations
    # Python slices are 0-based, end-exclusive -> subtract 1 only for start coords
    utr_parts = []

    if tscript.strand == "+":
        cds_end = cds_features[-1].end

        # Get all exons after the CDS end
        for exon in exon_features:
            if exon.end > cds_end:
                utr_start = max(exon.start, cds_end + 1) # One exon can contain CDS and 3'UTR content
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
        print(f"WARNING: unknown strand '{tscript.strand}' for feature '{tscript.id}'", file=sys.stderr)
        return None
    
    if not utr_parts:
        return None

    full_utr = ''.join(utr_parts)
    downstream_seq = ""
    if downstream_start <= downstream_end:
        downstream_seq = fasta[tscript.chrom][downstream_start - 1 : downstream_end].seq

    # Assemble final sequence
    if tscript.strand == "+":
        term_seq = full_utr + downstream_seq # 5' to 3'
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

    final_downstream_seq = term_seq[-args.downstream_nts:]
    if len(final_downstream_seq) < args.filter_window_size:
        return False

    if _is_internal_priming_artifact(
        final_downstream_seq, 
        args.filter_consecutive_a, 
        args.filter_window_a, 
        args.filter_window_size
    ):
        return False

    return True


def _process_transcript(tscript: gffutils.Feature, db: gffutils.FeatureDB, fasta: pyfaidx.Fasta, args: argparse.Namespace) -> str | None:
    """Processes a single transcript feature to extract and format its terminator sequence."""
    try:
        cds_features = list(db.children(tscript, featuretype="CDS", order_by="start"))
        exon_features = list(db.children(tscript, featuretype="exon", order_by="start"))
        if not cds_features or not exon_features:
            return None
        
        term_seq = _extract_terminator(tscript, fasta, cds_features, exon_features, args.downstream_nts)

        if not term_seq or not _filter_sequence(term_seq, args):
            return None
        
        return _format_fasta_record(tscript, term_seq, args.raw_dna)

    except Exception as e:
        print(f"WARNING: could not process feature '{tscript.id}': {e}", file=sys.stderr)
        return None


def _extract_all_terminators(fasta_fpath: str, gff_fpath: str, args: argparse.Namespace) -> None:
    """
    Extracts terminator sequences (3'UTR + downstream region) from the given fasta and gff files.

    Args:
        fasta_fpath (str): Filepath to the input FASTA file.
        gff_fpath (str): Filepath to the input GFF file.
        args (argparse.Namespace): Parsed command-line arguments.
    """
    assert os.path.isfile(fasta_fpath), f"Fasta file '{fasta_fpath}' does not exist."
    assert os.path.isfile(gff_fpath), f"GFF file '{gff_fpath}' does not exist."

    fname = os.path.basename(fasta_fpath).split('.')[0]
    print(f"Processing genome '{fname}'")

    db_dir = os.path.join("out", "gff_dbs") # Allows reuse of DBs between different taxons
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(db_dir, exist_ok=True)

    db_fpath = os.path.join(db_dir, f"{fname}.db")
    terminators_fpath = os.path.join(args.output_dir, f"{fname}_terminators.fa")

    fasta = pyfaidx.Fasta(fasta_fpath)
    db = _create_gff_db(gff_fpath, db_fpath)
    output_records = []
    skipped_count = 0

    for tscript in db.features_of_type("mRNA", order_by="start"):
        record = _process_transcript(tscript, db, fasta, args)
        if record:
            output_records.append(record)
        else:
            skipped_count += 1
        
    with open(terminators_fpath, 'w') as out_f:
        out_f.writelines(output_records)
    
    print(f"Extracted {len(output_records)} terminator sequences to '{terminators_fpath}'")
    print(f"skipped {skipped_count} transcripts")


# --- HELPER FUNCTIONS ---

def _format_fasta_record(tscript: gffutils.Feature, term_seq: str, raw_dna: bool) -> str:
    """
    Formats the given terminator sequence into a FASTA record header.

    Returns:
        str: The formatted FASTA record.
    """
    display_id = tscript.id
    if "orig_protein_id" in tscript.attributes:
        raw_id = tscript.attributes["orig_protein_id"][0]
        display_id = raw_id.split('|')[-1]
    
    header = f">{display_id} | {tscript.chrom}:{tscript.start}-{tscript.end}({tscript.strand})"

    term_seq = term_seq.upper()
    if not raw_dna:
        term_seq = term_seq.replace('T', 'U')
    wrapped_seq = textwrap.fill(term_seq, width=80)

    return f"{header}\n{wrapped_seq}\n"


def _find_files(dir: str) -> list[tuple[str, str]]:
    """
    Searches the given directory for matching pairs of FASTA and GFF files.

    Args:
        dir (str): The path to the directory to search.
    
    Returns:
        list[tuple[str, str]]: A list of tuples, where each tuple is a pair of FASTA and GFF file paths.
    """
    assert os.path.isdir(dir), f"Folder '{dir}' does not exist or is not a directory."

    valid_fasta_exts = ("fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn")
    valid_gff_exts = ("gff", "gff3")

    files_by_basename = defaultdict(dict)

    for filename in os.listdir(dir):
        basename, ext = os.path.splitext(filename)
        ext = ext.lstrip('.').lower()
        
        if ext in valid_fasta_exts:
            files_by_basename[basename]['fasta'] = os.path.join(dir, filename)
        elif ext in valid_gff_exts:
            files_by_basename[basename]['gff'] = os.path.join(dir, filename)
    
    file_pairs = []
    for basename, paths in files_by_basename.items():
        if 'fasta' in paths and 'gff' in paths:
            file_pairs.append((paths['fasta'], paths['gff']))
        else:
            print(f"WARNING: Incomplete file pair for '{basename}'. Skipping. Found: {list(paths.keys())}", file=sys.stderr)
    
    if not file_pairs:
        raise Exception(f"No valid FASTA and GFF file pairs found in directory '{dir}'.")

    return file_pairs


def _create_gff_db(gff_fpath: str, db_fpath: str) -> gffutils.FeatureDB:
    """
    Creates a gffutils FeatureDB from the given GFF file if it does not already exist.

    Args:
        gff_fpath (str): Filepath to the input GFF file.
        db_fpath (str): Filepath to the output database file.

    Returns:
        gffutils.FeatureDB: The created FeatureDB object.
    """
    assert os.path.isfile(gff_fpath), f"GFF file '{gff_fpath}' does not exist."

    if not os.path.isfile(db_fpath):
        print(f"Creating GFF database at '{db_fpath}'")
        gffutils.create_db(
            gff_fpath,
            dbfn=db_fpath,
            keep_order=True,
            merge_strategy="create_unique",
            sort_attribute_values=True,
        )
    else:
        print(f"Using existing GFF database at '{db_fpath}'")
    
    return gffutils.FeatureDB(db_fpath)
