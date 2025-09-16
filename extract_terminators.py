import os
import sys
import concurrent.futures
import argparse
import textwrap
from collections import defaultdict

# External libraries
import pyfaidx  # https://anaconda.org/bioconda/pyfaidx
import gffutils # https://anaconda.org/bioconda/gffutils

# utr always refers to 3'UTR unless otherwise specified

OUT_DIR = "out"
TERMINATORS_DIR = os.path.join(OUT_DIR, "terminators")
DB_DIR = os.path.join(OUT_DIR, "gff_dbs")

# --- Exceptions ---
class TerminatorExtractionError(Exception):
    """Base exception for this script."""
    pass

class GFFParsingError(TerminatorExtractionError):
    """Exception for errors parsing GFF files."""
    pass

class FileProcessingError(TerminatorExtractionError):
    """Exception for errors during file processing."""
    pass


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract terminators for the given genomes/genus."
    )
    parser.add_argument("input_dir", help="Path to the input folder containing FASTA and GFF files.")
    parser.add_argument(
        "--downstream-nts",
        type=int,
        default=50,
        help="Number of nucleotides to extract downstream of the terminator (default: 50)."
    )
    parser.add_argument(
        "--equal-length",
        action="store_true",
        help="Buffer terminator sequences so that they are all equal length."
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
        "--raw-dna",
        action="store_true",
        help="Output raw DNA sequences (5' to 3') instead of transcribed RNA."
    )

    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        parser.error(f"Input directory '{args.input_dir}' does not exist or is not a directory.")
    
    return args


def is_internal_priming_artifact(downstream_seq: str, filter_consecutive_a: int = 6, filter_window_a: int = 8, window_size: int = 10) -> bool:
    """
    Checks if the given downstream sequence is likely to be an internal priming artifact.
    Methodology based on Beaudong et al. DOI: 10.1101/gr.10.7.1001 
    Args:
        downstream_seq (str): The downstream sequence to check.
        filter_consecutive_a (int): The minimum number of consecutive 'A's to consider it an artifact (default: 6).
        filter_window_a (int): The minimum number of 'A's in the window to consider it an artifact (default: 8).
        window_size (int): The number of nucleotides to check from the start of the sequence (default: 10).
    """
    assert len(downstream_seq) >= window_size, "Downstream sequence length must be greater than or equal to window size."
    assert filter_consecutive_a >= 0, "Filter for consecutive A's must be non-negative."
    assert filter_window_a >= 0, "Filter for total A's in window must be non-negative."
    assert window_size > 0, "Window size must be greater than 0."

    if filter_consecutive_a == 0 and filter_window_a == 0:
        return False
    
    region_to_check = downstream_seq[:window_size].upper()

    # check for consecutive A's
    if filter_consecutive_a > 0 and 'A' * filter_consecutive_a in region_to_check:
        return True
    
    # check for total A's in window
    if filter_window_a > 0 and region_to_check.count('A') >= filter_window_a:
        return True
    
    return False


def extract_terminators(fasta_fpath: str, gff_fpath: str, out_fpath: str, args: argparse.Namespace) -> None:
    """
    Extracts terminator sequences (3'UTR + downstream region) from the given fasta and gff files.

    Args:
        fasta_fpath (str): Path to the fasta file.
        gff_fpath (str): Path to the gff file.
        out_fpath (str): Path to the output fasta file to write terminator sequences to.
        dstream_nts (int): Number of nucleotides to extract downstream of the 3'UTR.
    """
    assert os.path.isfile(fasta_fpath), f"Fasta file '{fasta_fpath}' does not exist."
    assert os.path.isfile(gff_fpath), f"GFF file '{gff_fpath}' does not exist."

    genome_name = os.path.basename(fasta_fpath).split('.')[0]
    os.makedirs(DB_DIR, exist_ok=True)
    db_path = os.path.join(DB_DIR, f"{genome_name}.db")

    if not os.path.isfile(db_path):
        try:
            print(f"Creating GFF database for '{genome_name}'")
            gffutils.create_db(
                gff_fpath,
                dbfn=db_path,
                keep_order=True,
                merge_strategy="create_unique",
                sort_attribute_values=True,
            )
        except Exception as e:
            raise GFFParsingError(f"Error creating GFF database for '{genome_name}': {e}")
    else:
        print(f"Using existing GFF database for '{genome_name}': '{db_path}'")
    
    fasta = pyfaidx.Fasta(fasta_fpath)
    db = gffutils.FeatureDB(db_path)
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
                print(f"WARNING: unknown strand '{tscript.strand}' for feature '{tscript.id}' in '{genome_name}'", file=sys.stderr)
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
            print(f"WARNING: could not process feature '{tscript.id}' in '{genome_name}': {e}", file=sys.stderr)
            continue
        
    with open(out_fpath, 'w') as out_f:
        out_f.writelines(output_records)
    
    print(f"Extracted {len(output_records)} terminator sequences from '{genome_name}'\nskipped {skipped_count} transcripts")


def find_files(dir: str) -> list[tuple[str, str]]:
    """
    Searches the given directory for matching pairs of FASTA and GFF files.

    Args:
        dir (str): The path to the directory to search.
    
    Returns:
        list[tuple[str, str]]: A list of tuples, where each tuple is a pair of fasta and gff file paths.
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
        raise FileProcessingError(f"No valid FASTA and GFF file pairs found in directory '{dir}'.")

    return file_pairs


def worker(args):
    """
    Extracts the terminator sequences of a single genome.
    Args:
        args (tuple): A tuple containing (fasta_fpath, gff_fpath, dstream_nts).
    """
    fasta_fpath, gff_fpath, cli_args = args
    genome_name = os.path.basename(fasta_fpath).split('.')[0]

    print(f"Processing genome '{genome_name}'")
    os.makedirs(TERMINATORS_DIR, exist_ok=True)
    terminators_fpath = os.path.join(TERMINATORS_DIR, f"{genome_name}_terminators.fa")

    extract_terminators(fasta_fpath, gff_fpath, terminators_fpath, cli_args)


def main() -> int:
    try:
        args = get_args()

        file_pairs = find_files(args.input_dir)

        tasks = [(file_pair[0], file_pair[1], args) for file_pair in file_pairs]

        # process each genome in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
            list(executor.map(worker, tasks))

    except (GFFParsingError, FileProcessingError, OSError) as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted by user.", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}", file=sys.stderr)
        return 1
    
    print("\nFINISHED")
    return 0


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
