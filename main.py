import re
import os

def parse_fasta(filepath):
    """
    Parses a FASTA file and yields (header, sequence) tuples.
    Args:
        filepath (str): The path to the FASTA file.
    Yields:
        tuple: (header_string, sequence_string)
    """
    header = None
    sequence_parts = []
    try:
        with open(filepath, 'r') as file_obj:
            for line in file_obj:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue
                if line.startswith(">"):
                    if header:  # If a previous sequence was being processed
                        yield header, "".join(sequence_parts)
                    header = line  # Store the new header
                    sequence_parts = []  # Reset sequence parts
                else:
                    if header:  # Only append if we are currently inside a FASTA record
                        sequence_parts.append(line)
            
            if header:  # Yield the last record in the file
                yield header, "".join(sequence_parts)
    except FileNotFoundError:
        print(f"Error: File not found at {filepath}")
        return
    except Exception as e:
        print(f"An error occurred while reading {filepath}: {e}")
        return

def extract_3utr(full_header, sequence):
    """
    Extracts the 3' UTR from a sequence based on CDS information in the header.
    Args:
        full_header (str): The full FASTA header line (e.g., ">rna-XYZ CDS=10-100").
        sequence (str): The nucleotide sequence.
    Returns:
        tuple: (new_header_string, utr_sequence_string) or None if CDS info is not found,
               is invalid, or if the 3' UTR is empty.
    """
    if "CDS=" not in full_header:
        return None

    # Use regex to find CDS coordinates robustly, e.g., CDS=123-456
    match = re.search(r"CDS=(\d+)-(\d+)", full_header)
    if not match:
        # print(f"Warning: Could not parse CDS coordinates from header: {full_header}")
        return None

    try:
        # The coordinates from gffread are 1-based.
        cds_end_coord = int(match.group(2))
    except ValueError:
        # print(f"Warning: Invalid CDS coordinates in header: {full_header}")
        return None

    # The sequence string is 0-indexed.
    # If the CDS ends at 1-based position 'cds_end_coord',
    # the 3' UTR starts at 0-based index 'cds_end_coord'.
    # For example, if CDS is 1-10, it occupies sequence[0] through sequence[9].
    # The 3' UTR starts at sequence[10].
    utr_sequence = sequence[cds_end_coord:]

    if not utr_sequence: # If the UTR sequence is empty
        return None

    # Create a new header for the 3' UTR
    # Remove '>' from the beginning of the original header for modification
    base_header = full_header[1:] 
    new_header = f">{base_header}_3UTR"
    
    return new_header, utr_sequence

def format_sequence_for_fasta(sequence, line_width=70):
    """
    Formats a sequence string into lines of a specified width.
    Args:
        sequence (str): The nucleotide sequence.
        line_width (int): The maximum number of characters per line.
    Returns:
        str: The formatted sequence with newlines.
    """
    return "\n".join(sequence[i:i+line_width] for i in range(0, len(sequence), line_width))

def main():
    # read filenames/folder from command line arguments
    # assumes folder is in current directory if only folder name is given
    path_to_genome_folder = ""

    # iterate through files in the specified folder and add full paths to list
    annotation_files = []

    # create an 'output_3utrs' folder in current directory if it doesn't exist
    output_path_prefix = r".\output_3utrs"

    if not os.path.exists(output_path_prefix):
        try:
            os.makedirs(output_path_prefix)
            print(f"Created output directory: {output_path_prefix}")
        except OSError as e:
            print(f"Error creating output directory {output_path_prefix}: {e}")
            return

    line_width = 70 # Default line width for FASTA output

    for genome in annotation_files:
        # use gffread to extract transcripts for every gene (all exons combined)


        # Construct output filename, e.g., "arenosa_3utrs.fa"
        base, ext = os.path.splitext(filename)
        if base.endswith("_transcripts"):
            output_base = base[:-len("_transcripts")] + "_3utrs"
        output_filename = output_base + ext
        output_filepath = os.path.join(output_path_prefix, output_filename)

        print(f"\nProcessing file: {filename}...")
        
        processed_count = 0
        utr_count = 0

        try:
            with open(output_filepath, 'w') as outfile:
                # parse_fasta yields, so we iterate through it.
                # It will handle FileNotFoundError for input_filepath internally.
                records_generator = parse_fasta(input_filepath)
                if records_generator is None: # Happens if parse_fasta had an early exit due to file error
                    print(f"Skipping {input_filepath} due to read error.")
                    continue

                for header, seq in records_generator:
                    processed_count += 1
                    result = extract_3utr(header, seq)
                    if result:
                        utr_header, utr_seq = result
                        utr_count += 1
                        outfile.write(utr_header + "\n")
                        outfile.write(format_sequence_for_fasta(utr_seq, line_width) + "\n")
            
            if processed_count > 0 or utr_count > 0 : # only print summary if file was attempted
                 print(f"Finished processing {input_filepath}.")
                 print(f"  Records processed in this file: {processed_count}")
                 print(f"  3' UTRs extracted to {output_filepath}: {utr_count}")
            elif os.path.exists(input_filepath): # If file exists but no records or UTRs, still note it.
                print(f"  No records found or no 3' UTRs extracted from {input_filepath}.")


        except IOError as e:
            print(f"Error writing to output file {output_filepath}: {e}")
        except Exception as e:
            print(f"An unexpected error occurred while processing {filename}: {e}")

if __name__ == "__main__":
    main()
