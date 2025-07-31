import os
import sys
import subprocess

NUM_CMD_LINE_ARGS = 1
TRANSCRIPTS_FOLDER_NAME = "gffread_transcripts"
OUTPUT_FOLDER_NAME = "terminator_sequences"

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
            sys.exit(1)

def main():
    # read input folder name from command line argument
    if len(sys.argv) != NUM_CMD_LINE_ARGS + 1:
        print(f"Usage: {sys.argv[0]} <path to input folder>")
        return
    
    input_folder = sys.argv[1]
    if not os.path.isdir(input_folder):
        print(f"Error: {input_folder} is not a valid directory.")
        return
    
    # iterate through files in the input folder and add filenames to list
    # TODO: make work with all types of fasta file and gff file extensions
    fasta_files = []
    annotation_files = []
    for filename in os.listdir(input_folder):
        if filename.endswith('.fa') or filename.endswith('.fasta') or filename.endswith('.fna'):
            fasta_files.append(filename)
        elif filename.endswith('.gff'):
            annotation_files.append(filename)

    if len(fasta_files) != len(annotation_files):
        print("Error: The number of genome files and annotation files must match.")
        return
    elif len(fasta_files) == 0:
        print("Error: No genome or annotation files found in the specified folder.")
        return
    
    filenames = []
    for i in range(len(fasta_files)):
        if fasta_files[i].split('.')[0] != annotation_files[i].split('.')[0]:
            print(f"Error: Mismatched genome and annotation filenames: {fasta_files[i]} and {annotation_files[i]}")
            return
        filenames.append((fasta_files[i], annotation_files[i]))
    
    create_folder(TRANSCRIPTS_FOLDER_NAME)
    create_folder(OUTPUT_FOLDER_NAME)

    for genome, annotation in filenames:
        # use gffread to extract transcripts for every gene (all exons combined)
        # gffread usage https://ccb.jhu.edu/software/stringtie/gff.shtml
        name = genome.split('.')[0]
        transcript_path = os.path.join(os.getcwd(), TRANSCRIPTS_FOLDER_NAME, name + ".fa")
        genome_path = os.path.join(input_folder, genome)
        annotation_path = os.path.join(input_folder, annotation)
        subprocess.run(["gffread", "-w", transcript_path, "-g", genome_path, annotation_path])

        # extract 3'UTRs from the transcripts


        # search annotation file for gene and extract global coords


        # extract X nucleotides from the end of the sequence

        # add to 3'UTR to form full terminator

if __name__ == "__main__":
    main()