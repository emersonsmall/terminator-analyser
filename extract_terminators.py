# gffread https://ccb.jhu.edu/software/stringtie/gff.shtml
# pyfastx https://pypi.org/project/pyfastx/
# LAPTOP: /mnt/c/Users/emers/Documents/QUT/2025/EGH400 Local/project-code

# TODO: replace gffread with a python library (e.g. gffutils, pybedtools, BioPython) these may have more robust parsing and would make the code cleaner
# need to verify they have the same functionality

import os
import sys
import subprocess
import functools
import concurrent.futures
import argparse
import re
import textwrap

# External libraries
import gffutils
import pyfastx

# --- Custom Exceptions ---
class TerminatorExtractionError(Exception):
    """Base exception for this script."""
    pass

class GFFParsingError(TerminatorExtractionError):
    """Exception for errors in GFF parsing."""
    pass

class FileProcessingError(TerminatorExtractionError):
    """Exception for errors during file processing."""
    pass


def get_args() -> argparse.Namespace:
    """
    Parses command line arguments using argparse.
    """
    parser = argparse.ArgumentParser(
        description="Extract terminators from GFF files and their corresponding FASTA files."
    )
    parser.add_argument("input_dir", help="Path to the input folder containing FASTA and GFF files.")
    parser.add_argument(
        "--enable-filtering",
        action="store_true",
        help="Enable filtering of GFF files if gffread encounters errors."
    )
    parser.add_argument(
        "--max-features",
        type=int,
        default=10,
        help="Maximum number of features allowed to be removed when filtering GFFs (default: 10)."
    )
    parser.add_argument(
        "--downstream-nucleotides",
        type=int,
        default=50,
        help="Number of nucleotides to extract downstream of the terminator (default: 50)."
    )
    parser.add_argument(
        "--separate-3utrs",
        action="store_true",
        help="Extract 3' UTRs separately (default: False)."
    )
    parser.add_argument(
        "--equal-length",
        action="store_true",
        help="Buffer terminator sequences so that they are all equal length."
    )

    args = parser.parse_args()

    if not os.path.isdir(args.input_dir):
        parser.error(f"Input directory '{args.input_dir}' does not exist or is not a directory.")
    
    return args


def parse_attributes(attributes: str) -> dict[str, str]:
    """
    Parses a GFF file attribute string into a dictionary.
    Args:
        attributes (str): The attribute string from a GFF file.
    Returns:
        dict[str, str]: A dictionary mapping attribute names to their values.
    """
    assert isinstance(attributes, str), f"Invalid type for parameter 'attributes'"

    attrs = {}
    for attr in attributes.split(';'):
        if '=' in attr:
            key, val = attr.split('=', 1)
            attrs[key.strip()] = val.strip()
    return attrs


def build_feature_map(gff_lines: list[str]) -> dict[str, list[dict[str, int | str | None]]]:
    """
    Creates a dictionary mapping each feature ID to its line index and parent. Handles duplicate feature IDs.
    Args:
        gff_lines (list[str]): A list of strings where each string is a line from a gff file.
    Returns:
        dict[str, list[dict[str, int | str | None]]]: A dictionary mapping each feature ID to a list of
        dictionaries, where each dictionary represents a line with said feature ID.
    """
    assert isinstance(gff_lines, list), f"Invalid type for parameter 'gff_lines'"

    feature_map = {}
    for i, line in enumerate(gff_lines):
        if line.startswith('#') or not line.strip():
            continue # skip comments and blank lines

        try:
            # 9th column contains attribute string
            attrs = parse_attributes(line.split('\t')[8])
            f_id = attrs.get("ID")

            if f_id:
                f_data = {
                    "idx": i,
                    "parent": attrs.get("Parent"),
                }

                if f_id not in feature_map:
                    feature_map[f_id] = []
                
                feature_map[f_id].append(f_data)
        
        except IndexError:
            continue # skip lines that do not have an attribute column

    return feature_map


def find_related_features(gff_lines: list[str], feature_id: str) -> tuple[str, list[int]]:
    """
    Searches a GFF file for the given ID and returns the line indices for the feature and any related features.
    Args:
        gff_lines (list[str]): A list of strings where each string is a line from a gff file.
        start_id (str): The ID of the feature to search for.
    Returns:
        tuple[str, list[int]]: The 1st element is the root id of the feature, and the 2nd element is a list containing 
                               the line indices of the lines related to the feature.
    """
    assert isinstance(gff_lines, list), f"Invalid type for parameter 'gff_lines'"
    assert isinstance(feature_id, str), f"Invalid type for parameter 'start_id'"

    feature_map = build_feature_map(gff_lines)

    @functools.lru_cache(maxsize=None)
    def get_root(f_id: str) -> str:
        if f_id not in feature_map:
            return f_id
        
        parent = feature_map[f_id][0].get("parent")

        if parent:
            return get_root(parent)

        return f_id # if no parent found, root ID has been reached
    
    root_id = get_root(feature_id)

    all_ids = [
        f_id for f_id in feature_map
        if get_root(f_id) == root_id
    ]

    all_idxs = []
    for f_id in all_ids:
        for part in feature_map[f_id]:
            all_idxs.append(part["idx"])

    return root_id, sorted(all_idxs)


def filter_gff(feature_id: str, in_fpath: str, out_fpath: str) -> str:
    """
    Filters a GFF file to remove all lines related to a given feature ID.
    Args:
        f_id (str): The ID of the feature to remove.
        in_path (str): The path to the input GFF file.
        out_path (str): The path to the output GFF file.
    Returns:
        str: The root ID (top-level feature ID i.e., gene ID) of the removed feature.
    """
    assert os.path.isfile(in_fpath), f"Input file {in_fpath} does not exist."
    assert isinstance(feature_id, str), f"Invalid type for parameter 'f_id'"
    assert isinstance(out_fpath, str), f"Invalid type for parameter 'out_path'"

    with open(in_fpath, 'r') as f:
        lines = f.readlines()
    
    root_id, line_idxs = find_related_features(lines, feature_id)
    
    lines_to_keep = [line for i, line in enumerate(lines) if i not in line_idxs]

    with open(out_fpath, 'w') as f:
        f.writelines(lines_to_keep)

    print(f"{len(line_idxs)} line/s removed from '{in_fpath}':")
    
    return root_id


def get_transcripts(fasta_fpath: str, gff_fpath: str, out_dir: str, out_fpath: str) -> None:
    """
    Extracts transcripts for the specified fasta file and gff file.
    """
    assert os.path.isfile(fasta_fpath), f"FASTA file {fasta_fpath} does not exist."
    assert os.path.isfile(gff_fpath), f"GFF file {gff_fpath} does not exist."
    assert os.path.isdir(out_dir), f"Output directory {out_dir} does not exist or is not a directory."
    assert isinstance(out_fpath, str), f"Invalid type for parameter 'out_fname'"

    gffread_errs = ("Error parsing", "GffObj::getSpliced() error: improper genomic coordinate")
    max_iterations = 10 # TODO: use cmd line arg

    genome_name = os.path.basename(fasta_fpath).split('.')[0]

    filtered_fname = genome_name + "_filtered.gff"
    filtered_gffs_dir = os.path.join(out_dir, "filtered_gffs")
    filtered_gff_fpath = os.path.join(filtered_gffs_dir, filtered_fname)
    features_removed: list[str] = []

    for i in range(max_iterations):
        if i > 0:
            print(f"\nRe-running gffread with '{filtered_fname}'...")
        if i == 1:
            gff_fpath = filtered_gff_fpath

        try:
            subprocess.run(
                [
                    "gffread",
                    "-w", out_fpath,
                    "-g", fasta_fpath,
                    gff_fpath
                ],
                capture_output=True, 
                text=True, 
                check=True)
            
            print(f"Transcripts successfully extracted to '{out_fpath}'")

            if features_removed:
                as_str = "\n".join(features_removed)
                print(f"Features removed from '{genome_name}':\n{as_str}")
            
            break
        except subprocess.CalledProcessError as err:
            print(f"\nError running gffread for '{genome_name}':\n{err.stderr}")

            is_parsing_error = err.stderr.startswith(gffread_errs[0])
            is_improper_coord_error = err.stderr.startswith(gffread_errs[1])

            if is_parsing_error:
                gff_line = err.stderr.split('\n')[1]
                attrs = parse_attributes(gff_line.split('\t')[-1])
                id = attrs.get("ID")
            elif is_improper_coord_error:
                id = err.stderr.split(' ')[-1].strip()
            else:
                raise GFFParsingError(
                    f"Unexpected error parsing GFF file '{gff_fpath}': {err.stderr}"
                )

            os.makedirs(filtered_gffs_dir, exist_ok=True)
            feature_removed = filter_gff(id, gff_fpath, filtered_gff_fpath)
            
            features_removed.append(feature_removed)


def get_post_cds(tscript_fpath: str, out_fpath: str, line_width: int = 70) -> None:
    """
    
    """
    assert os.path.isfile(tscript_fpath), f"FASTA file '{tscript_fpath}' does not exist."
    assert isinstance(out_fpath, str), f"Invalid type for parameter 'out_fpath'"

    tscripts = pyfastx.Fasta(tscript_fpath)
    utr_count = 0

    with open(out_fpath, 'w') as out_f:
        for record in tscripts:
            cds_match = re.search(r'CDS=(\d+)-(\d+)', record.description)

            if cds_match:
                cds_end_pos = int(cds_match.group(2))
                utr_seq = record.seq[cds_end_pos:]

                if utr_seq:
                    id = record.name.split()[0]  # Get the ID from the FASTA header

                    wrapped_utr_seq = textwrap.fill(utr_seq, width=line_width)

                    out_f.write(f">{id}_3utr\n{wrapped_utr_seq}\n")
                    utr_count += 1
    
    print(f"Extracted {utr_count} 3'UTRs from a total of {len(tscripts)} transcripts from '{tscript_fpath}' to '{out_fpath}'.")


def process_genome(file_pair: tuple[str, str], n: int, num_genomes: int, out_dir: str) -> None:
    """
    Processes a single genome.
    Args:
        file_pair (tuple[str, str]): A tuple containing the fasta and gff file paths for a single genome.
        n (int): The current genome number being processed.
        num_genomes (int): The total number of genomes to process.
        out_dir (str): The output directory where results will be saved.
    """
    assert isinstance(file_pair, tuple), f"Invalid type for parameter 'file_pair'"
    assert len(file_pair) == 2, f"Invalid length for parameter 'file_pair', expected 2 elements."
    assert isinstance(n, int), f"Invalid type for parameter 'n'"
    assert isinstance(num_genomes, int), f"Invalid type for parameter 'num_genomes'"

    fasta, gff = file_pair
    genome_name = os.path.basename(fasta).split('.')[0]

    print(f"\nProcessing genome '{genome_name}' ({n} of {num_genomes})...")

    print(f"Extracting transcripts for '{genome_name}'...")
    tscripts_dir = os.path.join(out_dir, "transcripts")
    os.makedirs(tscripts_dir, exist_ok=True)
    tscript_fpath = os.path.join(tscripts_dir, genome_name + "_transcripts.fa")
    get_transcripts(fasta, gff, out_dir, tscript_fpath)

    # extract 3'UTRs from transcripts
    utrs_dir = os.path.join(out_dir, "3utrs")
    os.makedirs(utrs_dir, exist_ok=True)
    utrs_fpath = os.path.join(utrs_dir, genome_name + "_3utrs.fa")
    get_post_cds(tscript_fpath, utrs_fpath)

    
    # extract X nucleotides from the end of the sequence for full terminator
    # USE --w-add <N> option. this modifies CDS= accordingly, can find 3'UTR by subtracting from the end


def find_files(dir: str) -> list[tuple[str, str]]:
    """
    Searches the given directory for valid fasta and gff files.
    Args:
        dir (str): The path to the directory to search.
    Returns:
        list[tuple[str, str]]: A list of tuples, where each tuple is a pair of fasta and gff file paths.
    """
    assert os.path.isdir(dir), f"Folder '{dir}' does not exist or is not a directory."

    valid_fasta_exts = ("fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn")
    valid_gff_exts = ("gff", "gff3")

    # iterate through all files in the input folder
    fasta_files: list[str] = []
    gff_files: list[str] = []
    for file in sorted(os.listdir(dir)):
        file_ext = file.split('.')[-1].lower()
        if file_ext in valid_fasta_exts:
            fasta_files.append(file)
        elif file_ext in valid_gff_exts:
            gff_files.append(file)
    
    if len(fasta_files) != len(gff_files):
        raise FileProcessingError(
            f"Number of fasta files ({len(fasta_files)}) does not match number of gff files ({len(gff_files)})."
        )

    files = list(zip(fasta_files, gff_files))
    i = 0
    for fasta, gff in files:
        if fasta.split('.')[0] != gff.split('.')[0]:
            raise FileProcessingError(f"Fasta file '{fasta}' and GFF file '{gff}' do not match in name.")
        else:
            files[i] = (os.path.join(dir, fasta), os.path.join(dir, gff))
        i += 1

    return files


def get_annotated_utrs(fasta_fpath: str, gff_fpath: str) -> list[str]:
    """
    Extracts annotated 3' UTRs from a GFF file and returns them as a list of strings.
    
    Uses pyfastx to fetch sequences from the FASTA file.
    Args:
        fasta_fpath (str): The path to the FASTA file.
        gff_fpath (str): The path to the GFF file.
    Returns:
        list[str]: A list of 3'UTR sequences.
    """
    assert os.path.isfile(fasta_fpath), f"FASTA file '{fasta_fpath}' does not exist."
    assert os.path.isfile(gff_fpath), f"GFF file '{gff_fpath}' does not exist."

    db_fpath = os.path.basename(gff_fpath) + ".db"

    if not os.path.exists(db_fpath):
        print(f"Creating GFF database for '{gff_fpath}'...")
        gffutils.create_db(
            gff_fpath,
            dbfn=db_fpath,
            force=True,
            merge_strategy="merge",
        )
    
    print(f"Loading GFF database from '{db_fpath}'...")
    db = gffutils.FeatureDB(db_fpath)

    genome = pyfastx.Fasta(fasta_fpath)

    seqs = []
    for utr in db.features_of_type("three_prime_UTR"):
        try:
            seq = genome[utr.chrom][utr.start - 1 : utr.end].seq

            if utr.strand == "-":
                seq = seq.reverse.complement
            
            seqs.append(seq)
        
        except (KeyError, IndexError) as err:
            print(f"Warning: could not find feature {utr.id} in FASTA file. Error: {err}", file=sys.stderr)
            continue

    return seqs

def main() -> int:
    exit_success = 0
    exit_failure = 1

    out_dir = "out" # TODO: move into process_genome function?
    terminators_dir = os.path.join(out_dir, "terminators")

    try:
        args = get_args()
        input_dir: str = args.input_dir

        files = find_files(input_dir)

        os.makedirs(terminators_dir, exist_ok=True)

        # worker function for multiprocessing
        worker_func = functools.partial(
            process_genome,
            num_genomes=len(files),
            out_dir=out_dir
        )

        # COMPARE TO ANNOTATED THALIANA
        ground_truth_utrs = get_annotated_utrs(
            os.path.join(input_dir, "thaliana.fna"),
            os.path.join(out_dir, "filtered_gffs", "thaliana_filtered.gff")
        )

        print(f"\nExtracted {len(ground_truth_utrs)} annotated 3' UTRs from the ground truth GFF file.")
        print(f"Ground truth 3' UTRs: {ground_truth_utrs}")

        # process each genome in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
            results = executor.map(worker_func, files, range(1, len(files) + 1))
            list(results) # consume iterator - raises an exception if any worker process raised an exception

    except GFFParsingError as err:
        print(f"\nGFF Parsing Error: {err}", file=sys.stderr)
        return exit_failure
    except FileProcessingError as err:
        print(f"\nFile Processing Error: {err}", file=sys.stderr)
        return exit_failure
    except OSError as err:
        print(f"\nOS Error: {err.strerror}", file=sys.stderr)
        return exit_failure
    except KeyboardInterrupt:
        print("\nInterrupted by user.", file=sys.stderr)
        return exit_failure
    except Exception as err:
        print(f"\nAn unexpected error occurred: {err}", file=sys.stderr)
        return exit_failure
    
    print("\nFINISHED")
    return exit_success


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
