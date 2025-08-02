# gffread https://ccb.jhu.edu/software/stringtie/gff.shtml
# /mnt/c/Users/Emerson/Documents/QUT/2025/EGH400 Local/project-code

import os
import sys
import subprocess

MIN_CMD_LINE_ARGS = 1 # path to input folder

SUCCESS_EXIT_CODE = 0
FAILURE_EXIT_CODE = 1

OUT_FOLDER = "output"
TRANSCRIPTS_FOLDER = os.path.join(OUT_FOLDER, "gffread_transcripts")
TERMINATORS_FOLDER = os.path.join(OUT_FOLDER, "terminators")
UTRS_FOLDER = os.path.join(OUT_FOLDER, "3utrs")
FILTERED_GFFS_FOLDER = os.path.join(OUT_FOLDER, "filtered_gffs")

FASTA_EXTENSIONS = {"fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn"}
GFF_EXTENSIONS = {"gff", "gff3"}


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

#TODO: fix. keep going down until you find a feature with a diff id and it's not a child feature of known features
# keep going up until you find a gene feature with no parent
# use index range to extract lines
def find_related_features(gff_file_path: str, id: str) -> list[str]:
    with open(gff_file_path, 'r') as f:
        lines = f.readlines()
    
    start_idx = -1
    for i, line in enumerate(lines):
        if f"ID={id};" in line:
            start_idx = i
            break
    
    if start_idx == -1:
        return []

    known_ids = {id}
    related_indices = {start_idx}

    print("start idx, id", start_idx, id)

    while True:
        print("iterating through related features...")
        added_feature = False

        min_idx = min(related_indices)
        max_idx = max(related_indices)

        # check line above current block
        above_idx = min_idx - 1
        if above_idx >= 0:
            line_above = lines[above_idx]
            attrs_above = parse_attributes(line_above.split('\t')[-1])
            id_above = attrs_above.get("ID")
            if id_above in known_ids:
                related_indices.add(above_idx)
            line_to_check = lines[min_idx]
            parent_id = parse_attributes(line_to_check.split('\t')[-1]).get("Parent")
            if id_above == parent_id:
                related_indices.add(above_idx)
                known_ids.add(id_above)
                added_feature = True
        
        # check line below current block
        below_idx = max_idx + 1
        if below_idx < len(lines):
            line_below = lines[below_idx]
            attrs_below = parse_attributes(line_below.split('\t')[-1])
            id_below = attrs_below.get("ID")
            if id_below in known_ids:
                related_indices.add(below_idx)
            parent_below = attrs_below.get("Parent")
            if parent_below in known_ids:
                related_indices.add(below_idx)
                id_below = attrs_below.get("ID")
                if id_below:
                    known_ids.add(id_below)
                added_feature = True
        
        # stop if no new relatives found
        if not added_feature:
            break
    
    # collect all lines
    # RETURN INDEX RANGE INSTEAD
    return [lines[i] for i in sorted(list(related_indices))]


def create_folder(name: str) -> bool:
    """
    Creates a folder in the current directory if it does not already exist.
    Args:
        name (str): The name of the folder to create.
    """

    path = os.path.join(os.getcwd(), name)
    if not os.path.exists(path):
        try:
            os.mkdir(path)
            print(f"Created directory: {path}")
        except OSError as e:
            print(f"Error creating directory {path}: {e}")
            return False
    
    return True


def extract_transcripts(input_folder: str, fasta_filename: str, gff_filename: str) -> int:
    """
    Extracts transcripts using gffread.
    Returns:
        int: SUCCESS_EXIT_CODE if successful, FAILURE_EXIT_CODE otherwise.
    """
    
    assert os.path.isdir(input_folder), f"Input folder {input_folder} does not exist or is not a directory."
    assert os.path.isfile(fasta_filename), f"FASTA file {fasta_filename} does not exist."
    assert os.path.isfile(gff_filename), f"GFF file {gff_filename} does not exist."

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
    assert os.path.isdir(input_folder), f"Input folder {input_folder} does not exist or is not a directory."

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
        print("Error: The number of genome files and annotation files must match.")
        return []
    
    files = []
    for i in range(len(fasta_files)):
        if fasta_files[i].split('.')[0] != gff_files[i].split('.')[0]:
            print(f"Error: Mismatched genome and annotation filenames: {fasta_files[i]} and {gff_files[i]}")
        files.append((fasta_files[i], gff_files[i]))
    
    if not files:
        print(f"Error: No valid genome and annotation files found in folder {input_folder}.")

    return files


def main() -> int:
    # read input folder name from command line argument
    if len(sys.argv) != MIN_CMD_LINE_ARGS + 1:
        print(f"Usage: {sys.argv[0]} <path to input folder>")
        return FAILURE_EXIT_CODE
    
    # TODO add cmd line arg for num nucleotides to extract to the right of CS

    input_folder = sys.argv[1]
    if not os.path.isdir(input_folder):
        print(f"Error: {input_folder} is not a valid directory.")
        return FAILURE_EXIT_CODE
    
    files = ingest_files(input_folder)
    if not files:
        return FAILURE_EXIT_CODE
    NUM_GENOMES = len(files)

    if not create_folder(TRANSCRIPTS_FOLDER) or not create_folder(UTRS_FOLDER) or not create_folder(TERMINATORS_FOLDER):
        return FAILURE_EXIT_CODE

    #TODO run this in parallel
    n = 1
    for fasta_file, gff_file in files:
        genome_name = fasta_file.split('.')[0]
        print(f"Processing genome \"{genome_name}\" ({n} of {NUM_GENOMES})...")

        # extract transcripts for every gene with gffread
        transcript_path = os.path.join(os.getcwd(), TRANSCRIPTS_FOLDER, genome_name + "_transcripts.fa")
        fasta_path = os.path.join(input_folder, fasta_file)
        gff_path = os.path.join(input_folder, gff_file)
        print(f"Extracting transcripts for \"{genome_name}\" with gffread...")
        try:
            subprocess.run(["gffread", "-w", transcript_path, "-g", fasta_path, gff_path], capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            print(e.stderr)
            if (e.stderr.startswith("Error parsing")):
                # TODO: re-run gffread with filtered GFF file in while loop until no parsing errors (max iterations?)
                # remove all lines related to problem gene from gff file
                gff_line = e.stderr.split('\n')[1]
                attrs = parse_attributes(gff_line.split('\t')[-1])
                id = attrs.get("ID")

                lines_to_remove = find_related_features(gff_path, id)
                print(lines_to_remove)

                with open(gff_path, 'r') as f:
                    lines = f.readlines()

                create_folder(FILTERED_GFFS_FOLDER)
                with open(os.path.join(FILTERED_GFFS_FOLDER, genome_name + "_filtered.gff"), 'w') as filtered_gff:
                    for line in lines:
                        if line not in lines_to_remove: # SIMPLY USE INDEX RANGE INSTEAD
                            filtered_gff.write(line)

                print(f"{len(lines_to_remove)} lines removed from {gff_file}.")

            else:
                return FAILURE_EXIT_CODE


        # extract 3'UTRs from transcripts


        # search annotation file for gene and find global coords
        

        # extract X nucleotides from the end of the sequence


        # add to 3'UTR to form full terminator


        n += 1
    
    print("Finished")
    return SUCCESS_EXIT_CODE


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
