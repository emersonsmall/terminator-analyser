# gffread https://ccb.jhu.edu/software/stringtie/gff.shtml
# /mnt/c/Users/Emerson/Documents/QUT/2025/EGH400 Local/project-code

import os
import sys
import subprocess

MIN_CMD_LINE_ARGS = 1 # path to input folder

SUCCESS_EXIT_CODE = 0
FAILURE_EXIT_CODE = 1

OUTPUT_FOLDER = "output"
TRANSCRIPTS_FOLDER = os.path.join(OUTPUT_FOLDER, "gffread_transcripts")
TERMINATORS_FOLDER = os.path.join(OUTPUT_FOLDER, "terminators")
UTRS_FOLDER = os.path.join(OUTPUT_FOLDER, "3utrs")
FILTERED_GFF_FOLDER = os.path.join(OUTPUT_FOLDER, "filtered_gffs")

FASTA_EXTENSIONS = {"fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn"}
GFF_EXTENSIONS = {"gff", "gff3"}


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
    i = 1
    for fasta_file, gff_file in files:
        genome_name = fasta_file.split('.')[0]
        print(f"Processing genome \"{genome_name}\" ({i} of {NUM_GENOMES})...")

        # use gffread to extract transcripts for every gene
        transcript_path = os.path.join(os.getcwd(), TRANSCRIPTS_FOLDER, genome_name + "_transcripts.fa")
        fasta_path = os.path.join(input_folder, fasta_file)
        gff_path = os.path.join(input_folder, gff_file)
        print(f"Extracting transcripts for \"{genome_name}\" with gffread...")
        try:
            subprocess.run(["gffread", "-w", transcript_path, "-g", fasta_path, gff_path], capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running gffread: {e}")
            print(e.stderr)
            if (e.stderr.startswith("Error parsing")):
                # TODO remove gene from gff file (create filtered gff file with only valid genes)
                gff_line = e.stderr.split('\n')[1]
                print(gff_line)
                columns = gff_line.split('\t')
                attributes = columns[8]
                id = attributes.split(';')[0].split('=')[1]
                child_id = id

                # filter out lines with the above id and any children of that id
                create_folder(FILTERED_GFF_FOLDER)
                with open(gff_path, 'r') as f:
                    lines = f.readlines()
                    lines_written = 0

                    for line in lines:
                        if line.find(f"Parent={id}"):
                            columns = line.split('\t')
                            attributes = columns[8]
                            child_id = attributes.split(';')[0].split('=')[1]
                        if not line.startswith('#') and id not in line and child_id not in line:
                            # write line to filtered gff file
                            path = os.path.join(FILTERED_GFF_FOLDER, genome_name + "_filtered.gff")
                            with open(path, 'w') as filtered:
                                filtered.write(line)
                            lines_written += 1

                print(f"{len(lines) - lines_written} lines removed from {gff_file}.")
                # try running gffread again with filtered gff file



            else:
                return FAILURE_EXIT_CODE


        # extract 3'UTRs from transcripts


        # search annotation file for gene and find global coords
        

        # extract X nucleotides from the end of the sequence


        # add to 3'UTR to form full terminator


        i += 1
    
    print("Finished")
    return SUCCESS_EXIT_CODE


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
