import sys
import os
import argparse
import subprocess
import glob
import textwrap
import pyfaidx

# TODO: add arg for minimum 3'UTR length??

NUE_START = -35 # -1 is the last nt of the 3'UTR
NUE_END = -5
CE_WIDTH = 30

def run_command(command: list[str]):
    print("Running command:", " ".join(command))
    try:
        process = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True
        )
        return process
    except FileNotFoundError:
        print(f"ERROR: Command not found: {command[0]}", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Command '{' '.join(command)}' failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)


def get_top_kmers(input_fasta: str, k_size: int, top_n: int) -> list:
    """Runs the jellyfish k-mer counting tool and returns the top N kmers."""
    jf_db_path = f"{input_fasta}.jf"
    dump_path = f"{input_fasta}.tsv"

    cpu_cores = os.cpu_count() or 1

    count_cmd = [
        "jellyfish", "count",
        "-m", str(k_size),
        "-s", "100M",
        "-t", str(cpu_cores),
        "-o", jf_db_path,
        "--canonical",
        input_fasta
    ]
    run_command(count_cmd)

    dump_cmd = [
        "jellyfish", "dump",
        "-c", jf_db_path,
        "-o", dump_path
    ]
    run_command(dump_cmd)

    # Parse and sort
    kmer_counts = []
    with open(dump_path, "r") as f:
        for line in f:
            kmer, count = line.strip().split()
            kmer_counts.append((kmer, int(count)))
    kmer_counts.sort(key=lambda x: x[1], reverse=True)

    # Cleanup
    os.remove(jf_db_path)
    os.remove(dump_path)

    return kmer_counts[:top_n]


def print_report(region_name: str, kmer_counts: list, k_size: int):
    print("\n" + "=" * 40)
    print(f"Top {len(kmer_counts)} K-mers for {region_name}")
    print("=" * 40)
    print(f"{'Rank':<5} | {'K-mer':<{k_size + 2}} | {'Count':>10}")
    print("-" * 40)
    for i, (kmer, count) in enumerate(kmer_counts):
        rank = i + 1
        print(f"{rank:<5} | {kmer:<{k_size + 2}} | {count:>10,}")


def generate_tiling_kmers(sequences: list[str], kmer_size: int) -> list[str]:
    """
    Generates k-mers from a list of sequences using a non-overlapping tiling approach.
    Returns a list of formatted FASTA records.
    """
    kmer_records = []
    kmer_id_counter = 0
    for seq in sequences:
        for i in range(0, len(seq) - kmer_size + 1, kmer_size):
            kmer = seq[i : i + kmer_size]
            kmer_records.append(f">kmer_{kmer_id_counter}\n{kmer}\n")
            kmer_id_counter += 1
    return kmer_records


def main():
    parser = argparse.ArgumentParser(
        description=textwrap.dedent("""
        Analyses the NUE and CE regions of the given set of terminator sequences.
        Finds the top N most common k-mers across all sequences using jellyfish.
        """)
    )
    parser.add_argument("input_dir", help="Path of the directory containing FASTA files.")
    parser.add_argument(
        "--downstream-nts", type=int, default=50,
        help="Number of downstream nucleotides included in the terminators (default: 50)."
    )
    parser.add_argument(
        "-N", "--top-n", type=int, default=50,
        help="The number of most common k-mers to report (default: 50)."
    )
    parser.add_argument(
        "-k", "--kmer-size", type=int, default=6,
        help="The size of the k-mers to analyze (default: 6)."
    )
    args = parser.parse_args()

    # Find all fasta files
    fasta_files = glob.glob(os.path.join(args.input_dir, "*_terminators.fa"))
    if not fasta_files:
        print(f"ERROR: No '*_terminators.fa' files found in directory '{args.input_dir}'", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(fasta_files)} FASTA files")
    
    # Extract all CE and NUE regions into separate temp files
    ce_fasta_path = "all_CEs.fa"
    nue_fasta_path = "all_NUEs.fa"
    print(f"Extracting CE and NUE regions")

    ce_seqs = []
    nue_seqs = []
    skipped = 0
    num_terminators = 0

    for fasta_file in fasta_files:
        fa_records = pyfaidx.Fasta(fasta_file, as_raw=True)
        for record in fa_records:
            num_terminators += 1

            seq = str(record).upper()
            total_len = len(seq)
            utr_len = total_len - args.downstream_nts

            # Slice NUE window
            if utr_len >= abs(NUE_START):
                nue_start = utr_len - abs(NUE_START)
                nue_end = utr_len - abs(NUE_END)
                nue_seq = seq[nue_start - 1 : nue_end]
                nue_seqs.append(nue_seq)
            else:
                skipped += 1
                continue

            # Slice CE window
            ce_start = utr_len - CE_WIDTH // 2
            ce_end = utr_len + CE_WIDTH // 2
            ce_seq = seq[ce_start - 1 : ce_end]
            ce_seqs.append(ce_seq)

    ce_records = generate_tiling_kmers(ce_seqs, args.kmer_size)
    nue_records = generate_tiling_kmers(nue_seqs, args.kmer_size)

    with open(ce_fasta_path, "w") as f:
        f.writelines(ce_records)
    with open(nue_fasta_path, "w") as f:
        f.writelines(nue_records)

    top_ce_kmers = get_top_kmers(ce_fasta_path, args.kmer_size, args.top_n)
    top_nue_kmers = get_top_kmers(nue_fasta_path, args.kmer_size, args.top_n)

    print(f"\n{skipped} of {num_terminators} ({(skipped/num_terminators * 100):.2f}%) terminator sequences skipped due to insufficient 3'UTR length")
    print_report("CE", top_ce_kmers, args.kmer_size)
    print_report("NUE", top_nue_kmers, args.kmer_size)

    os.remove(ce_fasta_path)
    os.remove(nue_fasta_path)

    print("FINISHED")

if __name__ == "__main__":
    main()
