import sys
import argparse
import pyfaidx
import os
from collections import defaultdict
from textwrap import fill

def parse_fasta(fasta_fpath) -> dict[str, str]:
    """
    Parses a FASTA file and returns a dictionary mapping IDs to sequences.
    - Key: Record ID (e.g., 'AT5G56260.2')
    - Value: Uppercase sequence
    """
    data = {}
    try:
        fasta_records = pyfaidx.Fasta(fasta_fpath, as_raw=True)
        for record in fasta_records:
            record_id = record.name.split()[0]
            seq = str(record).upper()
            data[record_id] = seq
    except pyfaidx.FastaIndexingError as e:
        print(f"ERROR: could not parse FASTA file '{fasta_fpath}': {e}", file=sys.stderr)
        sys.exit(1)
    
    return data


def main():
    parser = argparse.ArgumentParser(
        description="Compares two FASTA files record by record based on their IDs."
    )
    parser.add_argument("reference_fasta", help="Path to the reference FASTA file.")
    parser.add_argument("query_fasta", help="Path to the query FASTA file.")
    args = parser.parse_args()

    assert os.path.isfile(args.reference_fasta), f"Reference FASTA file '{args.reference_fasta}' does not exist."
    assert os.path.isfile(args.query_fasta), f"Query FASTA file '{args.query_fasta}' does not exist."

    print(f"Parsing reference FASTA")
    reference_records = parse_fasta(args.reference_fasta)
    print(f"Found {len(reference_records)} records in reference FASTA")

    print(f"Parsing query FASTA")
    query_records = parse_fasta(args.query_fasta)
    print(f"Found {len(query_records)} records in query FASTA")

    reference_ids = set(reference_records.keys())
    query_ids = set(query_records.keys())

    common_ids = reference_ids.intersection(query_ids)
    missing_from_query = reference_ids - query_ids
    extra_in_query = query_ids - reference_ids

    matches = 0
    mismatches = defaultdict(tuple)

    for record_id in common_ids:
        ref = reference_records[record_id]
        query = query_records[record_id]
        if ref == query:
            matches += 1
        else:
            mismatches[record_id] = (fill(ref, width=80), fill(query, width=80))
    
    # Summary
    print("\nSUMMARY")
    print(f"Records found in both files: {len(common_ids)}")
    print(f" - Matches: {matches}")
    print(f" - Mismatches: {len(mismatches)}")
    print(f"Records only in reference FASTA: {len(missing_from_query)}")
    print(f"Records only in query FASTA: {len(extra_in_query)}")
    
    if mismatches:
        print("\nMismatched Records:")
        for record_id in sorted(mismatches)[:20]:
            print("-" * 30)
            print(f" RECORD {record_id}")
            print(f"   Reference: \n{mismatches[record_id][0]}")
            print(f"   Query    : \n{mismatches[record_id][1]}")
            print("-" * 30)
        if len(mismatches) > 20:
            print(f" ... and {len(mismatches) - 20} more mismatches.")
        
    if missing_from_query:
        print("\nRecords only in reference FASTA:")
        for record_id in sorted(list(missing_from_query))[:20]:
            print(f" - {record_id}")
        if len(missing_from_query) > 20:
            print(f" ... and {len(missing_from_query) - 20} more records.")

if __name__ == "__main__":
    main()
