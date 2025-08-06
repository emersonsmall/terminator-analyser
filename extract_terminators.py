# gffread https://ccb.jhu.edu/software/stringtie/gff.shtml
# /mnt/c/Users/Emerson/Documents/QUT/2025/EGH400 Local/project-code
# /mnt/c/Users/emers/Documents/QUT/2025/EGH400 Local/project-code

import os
import sys
import subprocess
import functools

SUCCESS_EXIT_CODE = 0
FAILURE_EXIT_CODE = 1

GFFREAD_ERRORS = ("Error parsing", "GffObj::getSpliced() error: improper genomic coordinate")
MAX_GFFREAD_ITERATIONS = 10

def parse_cmd_line_args(args: list[str]) -> list[str]:
    """
    Parses command line arguments.
    """
    MIN_CMD_LINE_ARGS = 1 # path to input folder

    if len(sys.argv) != MIN_CMD_LINE_ARGS + 1:
        print(f"Usage: {sys.argv[0]} <path to input folder>")
        return []
    
    # TODO add cmd line arg for num nucleotides to extract to the right of CS
    # add arg for whether to filter GFFs or not, or set max features allowed to be removed

    input_dir = sys.argv[1]
    if not os.path.isdir(input_dir):
        print(f"Error: {input_dir} is not a valid directory.")
        return []
    
    return [input_dir]


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
    feature_map = {}
    for i, line in enumerate(gff_lines):
        if line.startswith('#') or not line.strip():
            continue # skip comments and blank lines

        try:
            # 9th column contains attributes string
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
        in_path (str): The path to the input GFF file.
        f_id (str): The ID of the feature to filter out.
        out_path (str): The path to the output GFF file.
    Returns:
        str: The root ID (top-level feature ID - i.e., gene ID) of the removed feature.
    """
    assert os.path.isfile(in_path), f"Input file {in_path} does not exist."
    assert isinstance(f_id, str), f"Feature ID must be a string"

    try:
        with open(in_path, 'r') as f:
            lines = f.readlines()
    except OSError as e:
        print(f"Error: Could not open input file '{in_path}': {e}")
        raise OSError
    
    root_id, line_idxs = find_related_features(lines, f_id)
    
    lines_to_keep = [line for i, line in enumerate(lines) if i not in line_idxs]

    try:
        with open(out_path, 'w') as f:
            f.writelines(lines_to_keep)
    except OSError as e:
        print(f"Error: Could not create/open output file '{out_path}': {e}")
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
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Created directory: {path}")


def extract_transcripts(in_dir: str, fasta_fname: str, gff_fname: str) -> int:
    """
    Extracts transcripts using gffread.
    Returns:
        int: SUCCESS_EXIT_CODE if successful, FAILURE_EXIT_CODE otherwise.
    """
    assert os.path.isdir(in_dir), f"Input folder {in_dir} does not exist or is not a directory."
    assert os.path.isfile(fasta_fname), f"FASTA file {fasta_fname} does not exist."
    assert os.path.isfile(gff_fname), f"GFF file {gff_fname} does not exist."

    return SUCCESS_EXIT_CODE


def extract_3utrs() -> int:
    """
    Extracts 3' UTRs from transcripts.
    Returns:
        int: SUCCESS_EXIT_CODE if successful, FAILURE_EXIT_CODE otherwise.
    """
    # TODO implement this function
    return SUCCESS_EXIT_CODE


def find_files(dir: str) -> list[tuple[str, str]]:
    """
    Searches the given directory for fasta and gff files and enforces gffread requirements.
    Args:
        dir (str): The path to the directory to search.
    Returns:
        list[tuple[str, str]]: A list of tuples, where each tuple is a pair of a fasta and gff filenames
    """
    assert os.path.isdir(dir), f"Folder {dir} does not exist or is not a directory."

    FASTA_EXTENSIONS = ("fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn")
    GFF_EXTENSIONS = ("gff", "gff3")

    # iterate through all files in the input folder
    fasta_files = []
    gff_files = []
    for file in os.listdir(dir):
        file_ext = file.split('.')[-1]
        if file_ext in FASTA_EXTENSIONS:
            fasta_files.append(file)
        elif file_ext in GFF_EXTENSIONS:
            gff_files.append(file)
    
    if len(fasta_files) != len(gff_files):
        print("Error: The number of fasta files and gff files must match.")
        return []
    
    files = []
    for i in range(len(fasta_files)):
        if fasta_files[i].split('.')[0] != gff_files[i].split('.')[0]:
            print(f"Error: Mismatched fasta and gff filenames: {fasta_files[i]} and {gff_files[i]}")
            files = []
            break
        files.append((fasta_files[i], gff_files[i]))

    return files


def main() -> int:
    OUT_DIR = "out"
    FILTERED_GFFS_DIR = os.path.join(OUT_DIR, "filtered_gffs")
    TRANSCRIPTS_DIR = os.path.join(OUT_DIR, "gffread_transcripts")
    TERMINATORS_DIR = os.path.join(OUT_DIR, "terminators")
    UTRS_DIR = os.path.join(OUT_DIR, "3utrs")
    ALL_DIRS = (FILTERED_GFFS_DIR, TRANSCRIPTS_DIR, TERMINATORS_DIR, UTRS_DIR)

    args = parse_cmd_line_args(sys.argv)
    if not args:
        return FAILURE_EXIT_CODE
    input_dir = args[0]

    FILES = find_files(input_dir)
    if not FILES:
        return FAILURE_EXIT_CODE

    for dir in ALL_DIRS:
        try:
            create_folder(dir)
        except OSError as e:
            print(f"Error creating folder {dir}")
            return FAILURE_EXIT_CODE

    #TODO run in parallel
    n = 1
    for fasta, gff in FILES:
        genome = fasta.split('.')[0]
        print(f"\nProcessing genome \"{genome}\" ({n} of {len(FILES)})...")

        print(f"Extracting transcripts for \"{genome}\" with gffread...")

        # extract transcripts for every gene with gffread
        tscript_path = os.path.join(TRANSCRIPTS_DIR, genome + "_transcripts.fa")
        fasta_path = os.path.join(input_dir, fasta)
        gff_path = os.path.join(input_dir, gff)

        first_run = True
        filtered_gff_path = os.path.join(FILTERED_GFFS_DIR, genome + "_filtered.gff")
        features_removed: list[str] = []

        for i in range(MAX_GFFREAD_ITERATIONS):
            if first_run:
                first_run = False
            else:
                print(f"\nRe-running gffread with {filtered_gff_path}...")
                gff_path = filtered_gff_path

            try:
                subprocess.run(
                    [
                        "gffread",
                        "-w", tscript_path, 
                        "-g", fasta_path, 
                        gff_path
                    ],
                    capture_output=True, 
                    text=True, 
                    check=True)
                
                print(f"Transcripts successfully extracted to {tscript_path}")

                if features_removed:
                    print(f"Features removed from \"{genome}\":\n{'\n'.join(features_removed)}")
                
                break
            except subprocess.CalledProcessError as err:
                print(err.stderr)

                is_parsing_error = err.stderr.startswith(GFFREAD_ERRORS[0])
                is_improper_coord_error = err.stderr.startswith(GFFREAD_ERRORS[1])
                if is_parsing_error:
                    # remove all lines related to problem feature
                    gff_line = err.stderr.split('\n')[1]
                    attrs = parse_attributes(gff_line.split('\t')[-1])
                    id = attrs.get("ID")
                elif is_improper_coord_error:
                    id = err.stderr.split(' ')[-1].strip()
                else:
                    print("\nError: unable to handle gffread error.")
                    return FAILURE_EXIT_CODE

                try:
                    feature_removed = filter_gff(id, gff_path, filtered_gff_path)
                except OSError as e:
                    return FAILURE_EXIT_CODE
                features_removed.append(feature_removed)
                
        # extract 3'UTRs from transcripts
        # COMPARE TO ANNOTATED THALIANA??? only a limited number available
        
        # extract X nucleotides from the end of the sequence
        # USE --w-add <N> gffread option

        # add to 3'UTR to form full terminator


        n += 1
    
    print("Finished")
    return SUCCESS_EXIT_CODE


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
