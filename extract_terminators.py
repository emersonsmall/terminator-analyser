# gffread https://ccb.jhu.edu/software/stringtie/gff.shtml
# pyfastx https://pypi.org/project/pyfastx/

# TAIR 3'UTR source https://www.arabidopsis.org/download/list?dir=Sequences%2FTAIR10_blastsets
# gene counts: https://www.ncbi.nlm.nih.gov/datasets/gene/taxon/81972/

# TODO add arg for max tolerance for 3'UTR differences for diff isoforms?

# TODO: explore gffread options for CDS=True annotation for more robust parsing?

# TODO: use ncbi API to retrieve given genomes/genus. still provide option to
#       specify local files. check if files exist, if not, download them. Can check against filenames that api provides

# TODO: pre-process/filter GFF files to remove problematic features before running gffread

# TODO: runs out of memory for thaliana when including 50 nts downstream

import os
import sys
import subprocess
import functools
import concurrent.futures
import argparse
import re
import textwrap

# External libraries
import pyfastx

OUT_DIR = os.path.join(os.getcwd(), "out")
FILTERED_GFFS_DIR = os.path.join(OUT_DIR, "filtered_gffs")
TSCRIPTS_DIR = os.path.join(OUT_DIR, "transcripts")
TERMINATORS_DIR = os.path.join(OUT_DIR, "terminators")
UTRS_DIR = os.path.join(OUT_DIR, "3utrs")

# --- Custom Exceptions ---
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
    """
    Parses command line arguments using argparse.
    """
    parser = argparse.ArgumentParser(
        description="Extract terminators for the given genomes/genus."
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
        help="Extract 3'UTRs separately (default: False)."
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


def parse_gff_attrs(attrs: str) -> dict[str, str]:
    """
    Parses a GFF feature attribute string into a dictionary.
    Args:
        attrs (str): The attribute string from a GFF file.
    Returns:
        dict[str, str]: A dictionary mapping attribute names to their values.
    """
    assert isinstance(attrs, str), f"Invalid type for parameter 'attrs'"

    attrs_dict = {}
    for attr in attrs.split(';'):
        if '=' in attr:
            key, val = attr.split('=', 1)
            attrs_dict[key.strip()] = val.strip()
    return attrs_dict


def build_feature_map(gff_lines: list[str]) -> dict[str, list[dict[str, int | str | None]]]:
    """
    Creates a dictionary mapping each feature ID to its line index and parent feature ID. Handles duplicate feature IDs.
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
            attrs = parse_gff_attrs(line.split('\t')[8])
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


def find_related_features(gff_lines: list[str], feature_id: str, feature_map: dict[str, list[dict[str, int | str | None]]]) -> tuple[str, list[int]]:
    """
    Searches the given GFF lines for the given feature ID and returns the root ID (gene ID) and line indices for the feature and any related features.
    Args:
        gff_lines (list[str]): A list of strings where each string is a line from a gff file.
        feature_id (str): The ID of the feature to search for.
        feature_map (dict[str, list[dict[str, int | str | None]]]): A dictionary mapping each feature ID to a list of lines.
    Returns:
        tuple[str, list[int]]: The root id of the feature, and a sorted list of all line indices related to the feature.
    """
    assert isinstance(gff_lines, list), f"Invalid type for parameter 'gff_lines'"
    assert isinstance(feature_id, str), f"Invalid type for parameter 'start_id'"
    assert isinstance(feature_map, dict), f"Invalid type for parameter 'feature_map'"

    @functools.lru_cache(maxsize=None)
    def get_root(f_id: str) -> str:
        if f_id not in feature_map:
            return f_id
        
        parent = feature_map[f_id][0].get("parent")

        if parent:
            return get_root(parent)

        return f_id # if no parent found, root ID reached
    
    root_id = get_root(feature_id)

    f_ids_with_common_root = [
        f_id for f_id in feature_map
        if get_root(f_id) == root_id
    ]

    line_idxs = []
    for f_id in f_ids_with_common_root:
        for feature in feature_map[f_id]:
            line_idxs.append(feature["idx"])

    return root_id, sorted(line_idxs)


def filter_gff(gff_lines: list[str], feature_id: str, in_fpath: str, out_fpath: str, feature_map: dict) -> str:
    """
    Filters a GFF file to remove all lines related to a given feature ID.
    Args:
        f_id (str): The ID of the feature to remove.
        in_path (str): The path to the input GFF file.
        out_path (str): The path to the output GFF file.
    Returns:
        str: The root ID (top-level feature ID i.e., gene ID) of the removed feature.
    """
    assert isinstance(feature_id, str), f"Invalid type for parameter 'f_id'"
    assert isinstance(out_fpath, str), f"Invalid type for parameter 'out_path'"
    
    root_id, line_idxs = find_related_features(gff_lines, feature_id, feature_map)
    
    lines_to_keep = [line for i, line in enumerate(gff_lines) if i not in line_idxs]

    with open(out_fpath, 'w') as f:
        f.writelines(lines_to_keep)
    
    return root_id


def get_transcripts(fasta_fpath: str, gff_fpath: str, out_fpath: str, max_iterations: int, dstream_nts: int) -> None:
    """
    Extracts transcripts for the specified fasta file and gff file.
    """
    assert os.path.isfile(fasta_fpath), f"FASTA file {fasta_fpath} does not exist."
    assert os.path.isfile(gff_fpath), f"GFF file {gff_fpath} does not exist."
    assert isinstance(out_fpath, str), f"Invalid type for parameter 'out_fname'"

    gffread_errs = ("Error parsing", "GffObj::getSpliced() error: improper genomic coordinate")

    genome_name = os.path.basename(fasta_fpath).split('.')[0]

    filtered_fname = genome_name + "_filtered.gff"
    filtered_gff_fpath = os.path.join(FILTERED_GFFS_DIR, filtered_fname)
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
                    #"--w-add", str(dstream_nts),
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

            is_parsing_err = err.stderr.startswith(gffread_errs[0])
            is_improper_coord_err = err.stderr.startswith(gffread_errs[1])

            if is_parsing_err:
                gff_line = err.stderr.split('\n')[1]
                attrs = parse_gff_attrs(gff_line.split('\t')[-1])
                id = attrs.get("ID")
            elif is_improper_coord_err:
                id = err.stderr.split(' ')[-1].strip()
            else:
                raise GFFParsingError(
                    f"Unexpected error parsing GFF file '{gff_fpath}': {err.stderr}"
                )

            os.makedirs(FILTERED_GFFS_DIR, exist_ok=True)
            with open(gff_fpath, 'r') as f:
                gff_lines = f.readlines()
            
            f_map = build_feature_map(gff_lines)
            feature_removed = filter_gff(gff_lines, id, gff_fpath, filtered_gff_fpath, f_map)
            
            features_removed.append(feature_removed)


def get_post_cds(tscripts_fpath: str, out_fpath: str) -> None:
    """
    
    """
    assert os.path.isfile(tscripts_fpath), f"FASTA file '{tscripts_fpath}' does not exist."
    assert isinstance(out_fpath, str), f"Invalid type for parameter 'out_fpath'"

    tscripts = pyfastx.Fasta(tscripts_fpath)
    count = 0

    output_records = []
    for record in tscripts:
        cds_match = re.search(r'CDS=(\d+)-(\d+)', record.description)

        if cds_match:
            cds_end_pos = int(cds_match.group(2))
            terminator_seq = record.seq[cds_end_pos:]

            if terminator_seq:
                id = record.name
                wrapped_seq = textwrap.fill(terminator_seq, width=80)
                output_records.append(f">{id}\n{wrapped_seq}\n")
                count += 1
    
    print(f"Extracted {count} terminators from {len(tscripts)} transcripts to '{out_fpath}'.")

    with open(out_fpath, 'w') as out_f:
        out_f.writelines(output_records)


def process_genome(file_pair: tuple[str, str], n: int, num_genomes: int, max_iterations: int, dstream_nts: int) -> None:
    """
    Processes a single genome.
    Args:
        file_pair (tuple[str, str]): A tuple containing the fasta and gff file paths for a single genome.
        n (int): The current genome number being processed.
        num_genomes (int): The total number of genomes to process.
    """
    assert isinstance(file_pair, tuple), f"Invalid type for parameter 'file_pair'"
    assert len(file_pair) == 2, f"Invalid length for parameter 'file_pair', expected 2 elements."
    assert isinstance(n, int), f"Invalid type for parameter 'n'"
    assert isinstance(num_genomes, int), f"Invalid type for parameter 'num_genomes'"

    fasta, gff = file_pair
    genome_name = os.path.basename(fasta).split('.')[0]

    print(f"\nProcessing genome '{genome_name}' ({n} of {num_genomes})...")

    print(f"Extracting transcripts for '{genome_name}'...")
    os.makedirs(TSCRIPTS_DIR, exist_ok=True)
    tscripts_fpath = os.path.join(TSCRIPTS_DIR, genome_name + "_transcripts.fa")
    get_transcripts(fasta, gff, tscripts_fpath, max_iterations, dstream_nts)

    # extract terminators from transcripts
    print(f"Extracting terminator sequences for '{genome_name}'...")
    os.makedirs(TERMINATORS_DIR, exist_ok=True)
    terminators_fpath = os.path.join(TERMINATORS_DIR, genome_name + "_terminators.fa")
    get_post_cds(tscripts_fpath, terminators_fpath)


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


def main() -> int:
    exit_success = 0
    exit_failure = 1

    try:
        args = get_args()
        input_dir = args.input_dir
        max_iterations = args.max_features
        dstream_nts = args.downstream_nucleotides

        files = find_files(input_dir)

        # worker function for multiprocessing
        worker_func = functools.partial(
            process_genome,
            num_genomes=len(files),
            max_iterations=max_iterations,
            dstream_nts=dstream_nts
        )

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
