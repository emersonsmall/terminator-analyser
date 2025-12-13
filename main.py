# Built-in libraries
import argparse
import sys
import os

# Local modules
from get_genomes import run_get_genomes, add_get_args
from extract import run_extraction, add_extract_args
from analyse import run_analysis, add_analyse_args


def main():
    parser = argparse.ArgumentParser(description="Terminator Analysis Pipeline")
    subparsers = parser.add_subparsers(dest="command", required=True)

    get_parser = subparsers.add_parser(
        "get",
        help="Get all reference genomes for a given taxon (genus or organism name).",
    )
    add_get_args(get_parser, is_standalone=True)

    extract_parser = subparsers.add_parser(
        "extract", help="Extract terminator sequences."
    )
    add_extract_args(extract_parser, is_standalone=True)

    analyse_parser = subparsers.add_parser(
        "analyse", help="Analyse NUEs and CEs of terminator sequences."
    )
    add_analyse_args(analyse_parser, is_standalone=True)

    full_parser = subparsers.add_parser(
        "full", help="Run the full pipeline: get, extract, and analyse."
    )
    full_parser.add_argument(
        "-o",
        "--output-dir",
        default="out",
        help="Path to the output directory (default: ./out).",
    )

    # Add unique args from all steps
    add_get_args(full_parser, is_standalone=False)
    add_extract_args(full_parser, is_standalone=False)
    add_analyse_args(full_parser, is_standalone=False)

    args = parser.parse_args()

    try:
        if args.command == "get":
            run_get_genomes(args)

        elif args.command == "extract":
            run_extraction(args)

        elif args.command == "analyse":
            run_analysis(args)

        elif args.command == "full":
            print("\nGETTING GENOMES")
            included_accessions = run_get_genomes(args)
            if not included_accessions:
                print("\nERROR: No accessions to process.", file=sys.stderr)

            args.input_dir = args.genomes_dir
            args.included_accessions = included_accessions

            print("\nEXTRACTING TERMINATORS")
            run_extraction(args)

            args.input_path = args.terminators_dir
            args.results_dir = os.path.join(args.results_dir, args.taxon.replace(" ", "_"))

            print("\nANALYSING TERMINATORS")
            run_analysis(args)

    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
