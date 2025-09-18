import argparse
import sys

from extract_terminators import run_extraction, get_args as get_extract_args
from analyse_terminators import run_analysis, get_args as get_analyse_args
from get_genomes import run_get_genomes, get_args as get_genomes_args

def main() -> int:
    """
    Main function for the terminator analysis pipeline.
    Uses subparsers to handle different stages: extraction, analysis, or full run.
    """
    parser = argparse.ArgumentParser(
        description="Terminator Analysis Pipeline"
    )
    subparsers = parser.add_subparsers(dest="command")
    subparsers.add_parser("get", help="Get all reference genomes for a given taxon.")
    subparsers.add_parser("extract", help="Extract terminator sequences from FASTA and GFF files.")
    subparsers.add_parser("analyse", help="Analyse NUEs and CEs of terminator sequences.")
    subparsers.add_parser("full", help="Runs the full pipeline: get, extract, and analyse")

    known, remaining = parser.parse_known_args()
    if not known.command:
        parser.print_help()
        return 1

    try:
        if known.command == "get":
            genomes_parser = get_genomes_args(return_parser=True)
            args = genomes_parser.parse_args(remaining)
            return run_get_genomes(args)

        elif known.command == "extract":
            extract_parser = get_extract_args(return_parser=True)
            args = extract_parser.parse_args(remaining)
            return run_extraction(args)
        
        elif known.command == "analyse":
            analyse_parser = get_analyse_args(return_parser=True)
            args = analyse_parser.parse_args(remaining)
            return run_analysis(args)

        elif known.command == "full":
            genomes_parser = get_genomes_args(return_parser=True)
            genomes_args, _ = genomes_parser.parse_known_args(remaining)
            exit_code = run_get_genomes(genomes_args)
            if exit_code != 0:
                print("\nERROR: Genome retrieval step failed. Aborting", file=sys.stderr)
                return exit_code

            extract_parser = get_extract_args(return_parser=True)
            extract_args, _ = extract_parser.parse_known_args(remaining)

            exit_code = run_extraction(extract_args)
            if exit_code != 0:
                print("\nERROR: Extraction step failed. Aborting", file=sys.stderr)
                return exit_code

            analyse_parser = get_analyse_args(return_parser=True)
            analyse_args, _ = analyse_parser.parse_known_args(remaining)

            analyse_args = argparse.Namespace(**vars(analyse_args)) # copy before modifying
            analyse_args.input_path = getattr(extract_args, 'output_dir', None)
            if analyse_args.input_path is None:
                print("\nERROR: Could not determine input path for analysis step. Aborting", file=sys.stderr)
                return 1

            return run_analysis(analyse_args)

    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        return 1

    print("\nFINISHED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
