import os

def main():
    # read input folder name from command line argument
    path_to_input_folder = ""
    # if full path not given, assume folder is in current directory
    # if len(args) < X:
        #print(Usage...)
        #return

    # iterate through files in the input folder and add filenames to list
    annotation_files = []

    # create a folder for gffread transcript files if it doesn't exist
    # TODO: make OS-agnostic
    gffread_transcripts_folder = r".\gffread_transcripts"
    if not os.path.exists(gffread_transcripts_folder):
        try:
            os.makedirs(gffread_transcripts_folder)
            print(f"Created gffread transcripts directory: {gffread_transcripts_folder}")
        except OSError as err:
            print(f"Error creating directory {gffread_transcripts_folder}: {err}")
            return

    # create an output folder if it doesn't exist
    output_folder = r".\output"
    if not os.path.exists(output_folder):
        try:
            os.makedirs(output_folder)
            print(f"Created output directory: {output_folder}")
        except OSError as err:
            print(f"Error creating directory {output_folder}: {err}")
            return
    
    for genome_file in annotation_files:
        full_path = os.path.join(path_to_input_folder, genome_file)

        # use gffread to extract transcripts for every gene (all exons combined)

        # extract 3'UTRs from the transcripts

        pass

    pass

if __name__ == "__main__":
    main()