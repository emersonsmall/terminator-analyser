# TAIR 3'UTR source: https://www.arabidopsis.org/download/list?dir=Sequences%2FAraport11_blastsets

import sys
import argparse
import pyfaidx
import os
from collections import defaultdict
from textwrap import fill
import re


def _get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compares two FASTA files record by record based on their IDs.",
    )
    parser.add_argument("reference_fasta", help="Path to the reference FASTA file.")
    parser.add_argument("query_fasta", help="Path to the query FASTA file.")

    args = parser.parse_args()
    if not os.path.isfile(args.reference_fasta):
        parser.error(f"Reference FASTA file '{args.reference_fasta}' does not exist.")
    if not os.path.isfile(args.query_fasta):
        parser.error(f"Query FASTA file '{args.query_fasta}' does not exist.")
    
    return args


def _parse_fasta(fasta_fpath, is_query) -> dict[str, str]:
    """
    Parses a FASTA file and returns a dictionary mapping IDs to sequences.
    - Key: Record ID (e.g., 'AT5G56260.2')
    - Value: Uppercase sequence
    """
    data = {}
    try:
        fasta_records = pyfaidx.Fasta(fasta_fpath, as_raw=True, read_long_names=True)
        for record in fasta_records:
            record_id = record.name.split()[0]
            seq = str(record).upper()
            header = record.name
            strand = "n/a"

            if is_query:
                match = re.search(r'\(([+-?])\)$', header)
                if match:
                    strand = match.group(1)
                
            data[record_id] = (seq, strand)

    except pyfaidx.FastaIndexingError as e:
        print(f"ERROR: could not parse FASTA file '{fasta_fpath}': {e}", file=sys.stderr)
        sys.exit(1)
    
    return data


def _run_comparison(args: argparse.Namespace) -> int:
    print(f"Parsing reference FASTA")
    reference_records = _parse_fasta(args.reference_fasta, False)
    print(f"Found {len(reference_records)} records in reference FASTA")

    print(f"Parsing query FASTA")
    query_records = _parse_fasta(args.query_fasta, True)
    print(f"Found {len(query_records)} records in query FASTA")

    reference_ids = set(reference_records.keys())
    query_ids = set(query_records.keys())

    common_ids = reference_ids.intersection(query_ids)
    missing_from_query = reference_ids - query_ids
    extra_in_query = query_ids - reference_ids

    matches = 0
    mismatches = defaultdict(tuple)

    for record_id in common_ids:
        ref, _ = reference_records[record_id]
        query, query_strand = query_records[record_id]
        if ref == query:
            matches += 1
        else:
            mismatches[record_id] = (fill(ref, width=80), fill(query, width=80), query_strand)
    
    num_mismatches = len(mismatches)
    num_missing_from_q = len(missing_from_query)
    num_common = len(common_ids)

    print("\nSUMMARY")
    print(f"Records found in both files: {num_common}")
    print(f" - Matches: {matches} ({matches / num_common * 100:.2f}%)")
    print(f" - Mismatches: {num_mismatches} ({num_mismatches / num_common * 100:.2f}%)")
    print(f"Records only in reference FASTA: {num_missing_from_q}")
    print(f"Records only in query FASTA: {len(extra_in_query)}")
    
    if mismatches:
        print("\nMismatched Records:")
        for record_id in sorted(mismatches)[:20]:
            print("-" * 30)
            print(f" RECORD {record_id}")
            print(f"   Reference: \n{mismatches[record_id][0]}")
            print(f"   Query    : \n{mismatches[record_id][1]}")
            print(f"Query Strand: \n{mismatches[record_id][2]}")
            print("-" * 30)
        if num_mismatches > 20:
            print(f" ... and {num_mismatches - 20} more mismatches.")
        
    if missing_from_query:
        print("\nRecords only in reference FASTA:")
        for record_id in sorted(list(missing_from_query))[:20]:
            print(f" - {record_id}")
        if num_missing_from_q > 20:
            print(f" ... and {num_missing_from_q - 20} more records.")
    
    return 0


# --- STANDALONE EXECUTION ---
def main():
    return _run_comparison(_get_args())


if __name__ == "__main__":
    sys.exit(main())
