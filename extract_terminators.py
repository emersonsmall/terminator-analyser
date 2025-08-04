# gffread https://ccb.jhu.edu/software/stringtie/gff.shtml
# /mnt/c/Users/Emerson/Documents/QUT/2025/EGH400 Local/project-code

import os
import sys
import subprocess

SUCCESS_EXIT_CODE = 0
FAILURE_EXIT_CODE = 1

GFFREAD_ERRORS = ("Error parsing", "GffObj::getSpliced() error: improper genomic coordinate")

def parse_cmd_line_args(args: list[str]) -> list[str]:
    """
    Parses command line arguments.
    """
    MIN_CMD_LINE_ARGS = 1 # path to input folder

    if len(sys.argv) != MIN_CMD_LINE_ARGS + 1:
        print(f"Usage: {sys.argv[0]} <path to input folder>")
        return []
    
    # TODO add cmd line arg for num nucleotides to extract to the right of CS

    input_folder = sys.argv[1]
    if not os.path.isdir(input_folder):
        print(f"Error: {input_folder} is not a valid directory.")
        return []
    
    return [input_folder]


def parse_attributes(attr_string: str) -> dict[str, str]:
    """
    Parses a GFF attribute string into a dictionary.
    Args:
        attr_string (str): The attribute string from a GFF file.
    Returns:
        dict[str, str]: A dictionary mapping attribute names to their values.
    """
    attrs = {}
    for attr in attr_string.split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attrs[key.strip()] = value.strip()
    return attrs


def find_related_features(gff_lines: list[str], id: str) -> tuple[str, list[int]]:
    """
    Searches a GFF file for the given ID and finds the range of lines that contain the feature and its related features.
    ASSUMPTIONS: related features form a block in the GFF file, i.e. they are consecutive lines.
    Args:
        lines (list[str]): A list of strings where each string is a line from a gff file.
        id (str): The ID of the feature to search for.
    Returns:
        list[int]: A list containing the start and end indices (inclusive) of the lines related to the feature
    """
    start_idx = -1
    for i, line in enumerate(gff_lines):
        if f"ID={id};" in line:
            start_idx = i
            break
    
    if start_idx == -1:
        print(f"Error: ID {id} not found in GFF file contents.")
        return tuple()

    root_id = id
    min_idx = max_idx = start_idx
    top_found = False
    bottom_found = False

    while not (top_found and bottom_found):
        if not top_found:
            # check line above current block
            above_idx = min_idx - 1
            if above_idx >= 0:
                curr_line = gff_lines[min_idx]
                curr_attrs = parse_attributes(curr_line.split('\t')[-1])
                curr_id = curr_attrs.get("ID")

                line_above = gff_lines[above_idx]
                attrs_above = parse_attributes(line_above.split('\t')[-1])
                id_above = attrs_above.get("ID")
                parent_above = attrs_above.get("Parent")

                if parent_above is None and curr_id != id_above:
                    root_id = curr_id
                    top_found = True
                    continue
                
                min_idx -= 1
        
        if not bottom_found:
            # check line below current block
            below_idx = max_idx + 1
            if below_idx < len(gff_lines):
                curr_line = gff_lines[max_idx]
                curr_attrs = parse_attributes(curr_line.split('\t')[-1])
                curr_id = curr_attrs.get("ID")

                line_below = gff_lines[below_idx]
                attrs_below = parse_attributes(line_below.split('\t')[-1])
                id_below = attrs_below.get("ID")
                parent_below = attrs_below.get("Parent")

                if parent_below is None and curr_id != id_below:
                    bottom_found = True
                    continue

                max_idx += 1
    
    return (root_id, [min_idx, max_idx])


def create_folder(name: str) -> bool:
    """
    Creates a folder in the current directory if it does not already exist.
    Args:
        name (str): The name of the folder to create.
    """
    path = os.path.join(os.getcwd(), name)
    if not os.path.exists(path):
        try:
            os.makedirs(path, exist_ok=True)
            print(f"Created directory: {path}")
        except OSError as e:
            print(f"Error creating directory {path}: {e}")
            return False
    
    return True


def extract_transcripts(input_folder: str, fasta_fname: str, gff_fname: str) -> int:
    """
    Extracts transcripts using gffread.
    Returns:
        int: SUCCESS_EXIT_CODE if successful, FAILURE_EXIT_CODE otherwise.
    """
    assert os.path.isdir(input_folder), f"Input folder {input_folder} does not exist or is not a directory."
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


def ingest_files(input_folder: str) -> list[tuple[str, str]]:
    assert os.path.isdir(input_folder), f"Folder {input_folder} does not exist or is not a directory."

    FASTA_EXTENSIONS = ("fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn")
    GFF_EXTENSIONS = ("gff", "gff3")

    # iterate through all files in the input folder
    fasta_files = []
    gff_files = []
    for file in os.listdir(input_folder):
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
    OUT_FOLDER = "out"
    TRANSCRIPTS_FOLDER = os.path.join(OUT_FOLDER, "gffread_transcripts")
    TERMINATORS_FOLDER = os.path.join(OUT_FOLDER, "terminators")
    UTRS_FOLDER = os.path.join(OUT_FOLDER, "3utrs")
    FILTERED_GFFS_FOLDER = os.path.join(OUT_FOLDER, "filtered_gffs")

    args = parse_cmd_line_args(sys.argv)
    if not args:
        return FAILURE_EXIT_CODE
    input_folder = args[0]

    files = ingest_files(input_folder)
    if not files:
        return FAILURE_EXIT_CODE

    if not create_folder(TRANSCRIPTS_FOLDER) or not create_folder(UTRS_FOLDER) or not create_folder(TERMINATORS_FOLDER):
        return FAILURE_EXIT_CODE

    #TODO run in parallel
    n = 1
    for fasta, gff in files:
        genome = fasta.split('.')[0]
        print(f"Processing genome \"{genome}\" ({n} of {len(files)})...")

        print(f"Extracting transcripts for \"{genome}\" with gffread...")

        # extract transcripts for every gene with gffread
        tscript_path = os.path.join(TRANSCRIPTS_FOLDER, genome + "_transcripts.fa")
        fasta_path = os.path.join(input_folder, fasta)
        gff_path = os.path.join(input_folder, gff)

        first_run = True
        filtered_fname = genome + "_filtered.gff"
        features_removed: list[str] = []
        create_folder(FILTERED_GFFS_FOLDER)

        # max iterations?
        while True:
            if first_run:
                first_run = False
            else:
                print("\nRe-running gffread with filtered GFF file...")
                gff_path = os.path.join(FILTERED_GFFS_FOLDER, filtered_fname)

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
                
                print(f"\nTranscripts successfully extracted to {tscript_path}")

                if features_removed:
                    print(f"Features removed from GFF file: {'\n'.join(features_removed)}")
                
                break
            except subprocess.CalledProcessError as err:
                print(err.stderr)

                if err.stderr.startswith(GFFREAD_ERRORS[0]): # parsing error

                    # remove all lines related to problem feature
                    gff_line = err.stderr.split('\n')[1]
                    attrs = parse_attributes(gff_line.split('\t')[-1])
                    id = attrs.get("ID")


                    with open(gff_path, 'r') as f:
                        lines = f.readlines()
                    
                    feature = find_related_features(lines, id)
                    lines_to_remove = feature[1]

                    with open(os.path.join(FILTERED_GFFS_FOLDER, filtered_fname), 'w') as f:
                        for i, line in enumerate(lines):
                            if i not in range(lines_to_remove[0], lines_to_remove[1] + 1):
                                f.write(line)
                    
                    print(f"Lines {lines_to_remove[0]} to {lines_to_remove[1]} removed from {gff_path}:")
                    for line_idx in range(lines_to_remove[0], lines_to_remove[1] + 1):
                        print(f"{lines[line_idx].strip()}")
                    
                    features_removed.append(feature[0])
                elif err.stderr.startswith(GFFREAD_ERRORS[1]): # improper genomic coordinate
                    print("IMPROPER GENOMIC COORDINATE ERROR")


                    return FAILURE_EXIT_CODE
                else:
                    print("\nError: unable to handle gffread error.")
                    return FAILURE_EXIT_CODE

        # extract 3'UTRs from transcripts
        
        # extract X nucleotides from the end of the sequence
        # USE --w-add <N> gffread option

        # add to 3'UTR to form full terminator


        n += 1
    
    print("Finished")
    return SUCCESS_EXIT_CODE


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
