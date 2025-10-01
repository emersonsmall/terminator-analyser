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
        parser.add_argument(
            "-d",
            "--downstream-nts",
            type=int,
            default=50,
            help="Number of nucleotides downstream of the CS to extract (default: 50)."
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


def is_internal_priming_artifact(sequence: str, consecutive_a: int = 6, total_a: int = 8, window_size: int = 10) -> bool:
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

    # check for consecutive A's
    if consecutive_a > 0 and 'A' * consecutive_a in region_to_check:
        return True
    
    # check for total A's in window
    if total_a > 0 and region_to_check.count('A') >= total_a:
        return True
    
    return False


def extract_terminators(fasta_fpath: str, gff_fpath: str, args: argparse.Namespace) -> None:
    """
    Extracts terminator sequences (3'UTR + downstream region) from the given fasta and gff files.

    Args:
        fasta_fpath (str): Filepath to the input FASTA file.
        gff_fpath (str): Filepath to the input GFF file.
        args (argparse.Namespace): Parsed command-line arguments.
    """
    assert os.path.isfile(fasta_fpath), f"Fasta file '{fasta_fpath}' does not exist."
    assert os.path.isfile(gff_fpath), f"GFF file '{gff_fpath}' does not exist."

    terminators_dir = args.output_dir
    db_dir = os.path.join(os.path.dirname(terminators_dir), "gff_dbs")
    os.makedirs(terminators_dir, exist_ok=True)
    os.makedirs(db_dir, exist_ok=True)

    fname = os.path.basename(fasta_fpath).split('.')[0]
    db_fpath = os.path.join(db_dir, f"{fname}.db")
    terminators_fpath = os.path.join(terminators_dir, f"{fname}_terminators.fa")

    if not os.path.isfile(db_fpath):
        try:
            print(f"Creating GFF database for '{fname}'")
            gffutils.create_db(
                gff_fpath,
                dbfn=db_fpath,
                keep_order=True,
                merge_strategy="create_unique",
                sort_attribute_values=True,
            )
        except Exception as e:
            raise Exception(f"Error creating GFF database for '{fname}': {e}")
    else:
        print(f"Using existing GFF database for '{fname}': '{db_fpath}'")
    
    fasta = pyfaidx.Fasta(fasta_fpath)
    db = gffutils.FeatureDB(db_fpath)
    output_records = []
    skipped_count = 0

    # All coordinates 1-based until modified in slice operations
    # Python slices are 0-based, end-exclusive. -1 only for start coords
    for tscript in db.features_of_type("mRNA", order_by="start"):
        try:
            cds_features = list(db.children(tscript, featuretype="CDS", order_by="start"))
            if not cds_features:
                skipped_count += 1
                continue
            
            exon_features = list(db.children(tscript, featuretype="exon", order_by="start"))
            if not exon_features:
                skipped_count += 1
                continue

            term_seq = ""
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
                downstream_end = tscript.end + args.downstream_nts

            elif tscript.strand == "-":
                cds_start = cds_features[0].start

                # Get all exons before the CDS start
                for exon in exon_features:
                    if exon.start < cds_start:
                        utr_start = exon.start
                        utr_end = min(exon.end, cds_start - 1)
                        utr_parts.append(fasta[tscript.chrom][utr_start - 1 : utr_end].seq)
                
                downstream_start = max(1, tscript.start - args.downstream_nts)
                downstream_end = tscript.start - 1

            else:
                print(f"WARNING: unknown strand '{tscript.strand}' for feature '{tscript.id}' in '{fname}'", file=sys.stderr)
                skipped_count += 1
                continue
            
            if not utr_parts:
                skipped_count += 1
                continue

            full_utr = ''.join(utr_parts)

            downstream_seq = ""
            if downstream_start <= downstream_end:
                downstream_seq = fasta[tscript.chrom][downstream_start - 1 : downstream_end].seq

            if tscript.strand == "+":
                term_seq = full_utr + downstream_seq # 5' to 3'
            else:
                term_seq = downstream_seq + full_utr
                term_seq = pyfaidx.Sequence(seq=term_seq).reverse.complement.seq

            if not args.raw_dna:
                final_downstream_seq = term_seq[-args.downstream_nts:]
                if len(final_downstream_seq) < args.filter_window_size:
                    skipped_count += 1
                    continue
            
                if is_internal_priming_artifact(final_downstream_seq, args.filter_consecutive_a, args.filter_window_a, args.filter_window_size):
                    skipped_count += 1
                    continue
            
            if term_seq:
                display_id = tscript.id
                if "orig_protein_id" in tscript.attributes:
                    raw_id = tscript.attributes["orig_protein_id"][0]
                    display_id = raw_id.split('|')[-1]
                header = f">{display_id} | {tscript.chrom}:{tscript.start}-{tscript.end}({tscript.strand})"

                term_seq = term_seq.upper()
                if not args.raw_dna:
                    term_seq = term_seq.replace('T', 'U')
                
                wrapped_seq = textwrap.fill(term_seq, width=80)
                output_records.append(f"{header}\n{wrapped_seq}\n")
            else:
                skipped_count += 1
                continue

        except Exception as e:
            print(f"WARNING: could not process feature '{tscript.id}' in '{fname}': {e}", file=sys.stderr)
            continue
        
    with open(terminators_fpath, 'w') as out_f:
        out_f.writelines(output_records)
    
    print(f"Extracted {len(output_records)} terminator sequences from '{fname}'\nskipped {skipped_count} transcripts")


def find_files(dir: str) -> list[tuple[str, str]]:
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


def worker(args: tuple) -> None:
    """
    Extracts the terminator sequences of a single genome.

    Args:
        args (tuple): A tuple containing (fasta_fpath, gff_fpath, cli_args).
    """
    fasta_fpath, gff_fpath, cli_args = args
    fname = os.path.basename(fasta_fpath).split('.')[0]

    print(f"Processing genome '{fname}'")
    extract_terminators(fasta_fpath, gff_fpath, cli_args)


def run_extraction(args: argparse.Namespace) -> int:
    file_pairs = find_files(args.input_path)
    tasks = [(file_pair[0], file_pair[1], args) for file_pair in file_pairs]

    # process each genome in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
        list(executor.map(worker, tasks))

    print("\nTerminator extraction finished")
    return 0


def main() -> int:
    """Standalone execution entry point."""
    try:
        return run_extraction(_get_args())
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())
