# gffread https://ccb.jhu.edu/software/stringtie/gff.shtml
# DESKTOP: /mnt/c/Users/Emerson/Documents/QUT/2025/EGH400 Local/project-code
# LAPTOP: /mnt/c/Users/emers/Documents/QUT/2025/EGH400 Local/project-code

# TODO: replace gffread with a python library (e.g. gffutils, pybedtools, BioPython) these may have more robust parsing and would make the code cleaner
# need to verify they have the same functionality

import os
import sys
import subprocess
import functools
import concurrent.futures
import argparse

# GLOBALS
EXIT_SUCCESS: int = 0
EXIT_FAILURE: int = 1

def parse_cmd_line_args() -> argparse.Namespace:
    """
    Parses command line arguments using argparse.
    """
    parser = argparse.ArgumentParser(
        description="Extract terminators from GFF files and their corresponding FASTA files."
    )
    parser.add_argument("input_dir", help="Path to the input folder containing FASTA and GFF files.")
    parser.add_argument(
        "--filter",
        action="store_true",
        help="Enable filtering of GFF files if gffread encounters errors (default: False)."
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
        help="Buffer terminator sequences so that they are all the same length (default: False)."
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


def find_related_features(gff_lines: list[str], start_id: str) -> tuple[str, list[int]]:
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
    assert isinstance(start_id, str), f"Invalid type for parameter 'start_id'"

    feature_map = build_feature_map(gff_lines)

    if start_id not in feature_map:
        print(f"Error: ID {start_id} not found.")
        return tuple()

    @functools.lru_cache(maxsize=None)
    def get_root(f_id: str) -> str:
        if f_id not in feature_map:
            return f_id
        
        parent = feature_map[f_id][0].get("parent")

        if parent:
            return get_root(parent)

        return f_id
    
    root_id = get_root(start_id)

    all_ids = [
        f_id for f_id in feature_map
        if get_root(f_id) == root_id
    ]

    all_idxs = []
    for f_id in all_ids:
        for part in feature_map[f_id]:
            all_idxs.append(part["idx"])

    if not all_idxs:
        return tuple()

    return root_id, sorted(all_idxs)


def filter_gff(f_id: str, in_path: str, out_path: str) -> str:
    """
    Filters a GFF file to remove all lines related to a given feature ID.
    Args:
        f_id (str): The ID of the feature to remove.
        in_path (str): The path to the input GFF file.
        out_path (str): The path to the output GFF file.
    Returns:
        str: The root ID (top-level feature ID i.e., gene ID) of the removed feature.
    """
    assert os.path.isfile(in_path), f"Input file {in_path} does not exist."
    assert isinstance(f_id, str), f"Invalid type for parameter 'f_id'"
    assert isinstance(out_path, str), f"Invalid type for parameter 'out_path'"

    with open(in_path, 'r') as f:
        lines = f.readlines()
    
    root_id, line_idxs = find_related_features(lines, f_id)
    
    lines_to_keep = [line for i, line in enumerate(lines) if i not in line_idxs]

    try:
        with open(out_path, 'w') as f:
            f.writelines(lines_to_keep)
    except OSError as e:
        print(f"Error: Could not create/open output file \"{out_path}\": {e}")
        raise OSError

    print(f"{len(line_idxs)} line/s removed from {in_path}:")
    for i in line_idxs:
        print(lines[i], end="")
    
    return root_id


def create_folder(path: str) -> None:
    """
    Creates a folder at the specified path if it does not already exist.
    Args:
        path (str): The path of the folder to create.
    """
    assert isinstance(path, str), f"Invalid type for parameter 'path'"

    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Created directory: {path}")


def extract_transcripts(fasta_path: str, gff_path: str, out_path: str) -> None:
    """
    Extracts transcripts for the specified fasta file and gff file.
    """
    assert os.path.isfile(fasta_path), f"FASTA file {fasta_path} does not exist."
    assert os.path.isfile(gff_path), f"GFF file {gff_path} does not exist."
    assert isinstance(out_path, str), f"Invalid type for parameter 'out_path'"

    GFFREAD_ERRORS = ("Error parsing", "GffObj::getSpliced() error: improper genomic coordinate")
    MAX_ITERATIONS = 10

    name = os.path.basename(fasta_path).split('.')[0]
    filtered_name = name + "_filtered.gff"

    filtered_gff_path = os.path.join(FILTERED_GFFS_DIR, filtered_name)
    features_removed: list[str] = []

    for i in range(MAX_ITERATIONS):
        if i > 0:
            print(f"\nRe-running gffread with \"{filtered_name}\"...")
        if i == 1:
            gff_path = filtered_gff_path

        try:
            subprocess.run(
                [
                    "gffread",
                    "-w", out_path,
                    "-g", fasta_path,
                    gff_path
                ],
                capture_output=True, 
                text=True, 
                check=True)
            
            print(f"Transcripts successfully extracted to \"{out_path}\"")

            if features_removed:
                print(f"Features removed from \"{name}\":\n{'\n'.join(features_removed)}")
            
            break
        except subprocess.CalledProcessError as err:
            print(f"\nError running gffread for {name}:")
            print(err.stderr)

            is_parsing_error = err.stderr.startswith(GFFREAD_ERRORS[0])
            is_improper_coord_error = err.stderr.startswith(GFFREAD_ERRORS[1])

            if is_parsing_error:
                gff_line = err.stderr.split('\n')[1]
                attrs = parse_attributes(gff_line.split('\t')[-1])
                id = attrs.get("ID")
            elif is_improper_coord_error:
                id = err.stderr.split(' ')[-1].strip()
            else:
                print("\nError: unable to handle gffread error.")
                raise Exception

            try:
                create_folder(FILTERED_GFFS_DIR)
                feature_removed = filter_gff(id, gff_path, filtered_gff_path)
            except OSError as e:
                print(f"Error: {e}")
                raise Exception
            
            features_removed.append(feature_removed)


def process_genome(file_pair: tuple[str, str], n: int, num_genomes: int, transcripts_dir: str) -> None:
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
    name = os.path.basename(fasta).split('.')[0]

    print(f"\nProcessing genome '{name}' ({n} of {num_genomes})...")

    print(f"Extracting transcripts for '{name}'...")
    create_folder(transcripts_dir)
    tscript_path = os.path.join(transcripts_dir, name + "_transcripts.fa")
    extract_transcripts(fasta, gff, tscript_path)

    # extract 3'UTRs from transcripts
    # COMPARE TO ANNOTATED THALIANA only a limited number available
    
    # extract X nucleotides from the end of the sequence
    # USE --w-add <N> option. this modifies CDS= accordingly, can find 3'UTR by subtracting from the end

    # add to 3'UTR to form full terminator


def find_files(dir: str) -> list[tuple[str, str]]:
    """
    Searches the given directory for valid fasta and gff files.
    Args:
        dir (str): The path to the directory to search.
    Returns:
        list[tuple[str, str]]: A list of tuples, where each tuple is a pair of fasta and gff file paths.
    """
    assert os.path.isdir(dir), f"Folder '{dir}' does not exist or is not a directory."

    VALID_FASTA_EXTS = ("fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn")
    VALID_GFF_EXTS = ("gff", "gff3")

    # iterate through all files in the input folder
    fasta_files: list[str] = []
    gff_files: list[str] = []
    for file in os.listdir(dir):
        file_ext = file.split('.')[-1].lower()
        if file_ext in VALID_FASTA_EXTS:
            fasta_files.append(file)
        elif file_ext in VALID_GFF_EXTS:
            gff_files.append(file)
    
    if len(fasta_files) != len(gff_files):
        print("Error: The number of fasta files and gff files must match.")
        return []

    files = list(zip(fasta_files, gff_files))
    i = 0
    for fasta, gff in files:
        if fasta.split('.')[0] != gff.split('.')[0]:
            print(f"Error: Mismatched filenames: {fasta} and {gff}")
            return []
        else:
            files[i] = (os.path.join(dir, fasta), os.path.join(dir, gff))
        i += 1

    return files


def extract_3utrs(fasta_path: str, out_path: str) -> None:
    """
    
    """
    assert os.path.isfile(fasta_path), f"FASTA file {fasta_path} does not exist."
    assert isinstance(out_path, str), f"Invalid type for parameter 'out_path'"
    


def main() -> int:
    OUT_DIR: str = "out"
    UTRS_DIR: str = os.path.join(OUT_DIR, "3utrs")
    TRANSCRIPTS_DIR: str = os.path.join(OUT_DIR, "transcripts")
    FILTERED_GFFS_DIR: str = os.path.join(OUT_DIR, "filtered_gffs")
    TERMINATORS_DIR: str = os.path.join(OUT_DIR, "terminators")
    DIRS = (TERMINATORS_DIR, UTRS_DIR, TRANSCRIPTS_DIR)

    try:
        args = parse_cmd_line_args()
        input_dir: str = args.input_dir

        FILES = find_files(input_dir)
        if not FILES:
            print(f"Error: No valid fasta or gff files found in {input_dir}.")
            return EXIT_FAILURE


        for dir in DIRS:
            try:
                create_folder(dir)
            except OSError as err:
                print(f"Error creating folder {dir}: {err}")
                return EXIT_FAILURE

        # set up worker function for multiprocessing
        worker_func = functools.partial(
            process_genome,
            num_genomes=len(FILES),
            transcripts_dir=TRANSCRIPTS_DIR,
        )

        # process each genome in parallel
        with concurrent.futures.ProcessPoolExecutor(max_workers=None) as executor:
            results = executor.map(worker_func, FILES, range(1, len(FILES) + 1))
            list(results) # consume iterator object - raises an exception if any worker process failed
    
    except KeyboardInterrupt:
        print("\nProcess interrupted by user.", file=sys.stderr)
        return EXIT_FAILURE
    except Exception as err:
        print(f"\nAn unexpected error occurred: {err}", file=sys.stderr)
        return EXIT_FAILURE

    
    print("\nFINISHED")
    return EXIT_SUCCESS


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
