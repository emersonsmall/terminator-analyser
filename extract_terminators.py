# gffread https://ccb.jhu.edu/software/stringtie/gff.shtml

import os
import sys
import subprocess

NUM_CMD_LINE_ARGS = 1

SUCCESS_EXIT_CODE = 0
FAILURE_EXIT_CODE = 1

OUTPUT_FOLDER = "output"
TRANSCRIPTS_FOLDER = os.path.join(OUTPUT_FOLDER, "gffread_transcripts")
TERMINATORS_FOLDER = os.path.join(OUTPUT_FOLDER, "terminators")
UTRS_FOLDER = os.path.join(OUTPUT_FOLDER, "3utrs")

ALLOWED_FASTA_EXTENSIONS = {"fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn"}
ALLOWED_GFF_EXTENSIONS = {"gff", "gff3"}

def create_folder(path: str) -> None:
    """
    Creates a folder at the specified path if it does not already exist.
    Args:
        path (str): The path to the folder to create.
    """

    full_path = os.path.join(os.getcwd(), path)
    if not os.path.exists(full_path):
        try:
            os.mkdir(full_path)
            print(f"Created directory: {full_path}")
        except OSError as e:
            print(f"Error creating directory {full_path}: {e}")
            sys.exit(FAILURE_EXIT_CODE)


def main():
    # read input folder name from command line argument
    if len(sys.argv) != NUM_CMD_LINE_ARGS + 1:
        print(f"Usage: {sys.argv[0]} <path to input folder>")
        return FAILURE_EXIT_CODE
    
    input_folder = sys.argv[1]
    if not os.path.isdir(input_folder):
        print(f"Error: {input_folder} is not a valid directory.")
        return FAILURE_EXIT_CODE
    
    # iterate through files in the input folder
    fasta_files = []
    gff_files = []
    for file in os.listdir(input_folder):
        file_ext = file.split('.')[1]
        if file_ext in ALLOWED_FASTA_EXTENSIONS:
            fasta_files.append(file)
        elif file_ext in ALLOWED_GFF_EXTENSIONS:
            gff_files.append(file)

    if len(fasta_files) != len(gff_files):
        print("Error: The number of genome files and annotation files must match.")
        return FAILURE_EXIT_CODE
    elif len(fasta_files) == 0:
        print("Error: A genome sequence file (FASTA format) and an annotation file (GFF format) must be provided.")
        return FAILURE_EXIT_CODE
    
    files = []
    for i in range(len(fasta_files)):
        if fasta_files[i].split('.')[0] != gff_files[i].split('.')[0]:
            print(f"Error: Mismatched genome and annotation filenames: {fasta_files[i]} and {gff_files[i]}")
            return FAILURE_EXIT_CODE
        files.append((fasta_files[i], gff_files[i]))
    
    NUM_GENOMES = len(fasta_files)

    create_folder(TRANSCRIPTS_FOLDER)
    create_folder(UTRS_FOLDER)
    create_folder(TERMINATORS_FOLDER)

    print("Extracting transcripts with gffread...")
    i = 1
    for genome, annotation in files:
        # use gffread to extract transcripts for every gene
        transcript_path = os.path.join(os.getcwd(), TRANSCRIPTS_FOLDER, genome.split('.')[0] + "_transcripts" + ".fa")
        genome_path = os.path.join(input_folder, genome)
        annotation_path = os.path.join(input_folder, annotation)
        print(f"Processing genome {i} of {NUM_GENOMES}...")
        subprocess.run(["gffread", "-w", transcript_path, "-g", genome_path, annotation_path])

        print("Extracting transcripts finished...")

        # extract 3'UTRs from transcripts


        # search annotation file for gene and find global coords
        

        # extract X nucleotides from the end of the sequence


        # add to 3'UTR to form full terminator
        i += 1
    
    print("Script finished")
    return SUCCESS_EXIT_CODE


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)