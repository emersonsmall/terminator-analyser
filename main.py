import argparse
import sys
import os

from get_genomes import run_get_genomes, add_get_args
from extract import run_extraction, add_extract_args
from analyse import run_analysis, add_analyse_args


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Terminator Analysis Pipeline"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    get_parser = subparsers.add_parser("get", help="Get all reference genomes (FASTA and GFF files) for a given taxon.")
    add_get_args(get_parser, standalone=True)

    extract_parser = subparsers.add_parser("extract", help="Extract terminator sequences from FASTA and GFF files.")
    add_extract_args(extract_parser, standalone=True)

    analyse_parser = subparsers.add_parser("analyse", help="Analyse NUEs and CEs of terminator sequences.")
    add_analyse_args(analyse_parser, standalone=True)

    full_parser = subparsers.add_parser("full", help="Run the full pipeline: get, extract, and analyse.")
    full_parser.add_argument(
        "-o", "--output-dir", default="out",
        help="Path to the output directory (default: ./out)."
    )
    
    # Add unique args from all steps
    add_get_args(full_parser, standalone=False)
    add_extract_args(full_parser, standalone=False)
    add_analyse_args(full_parser, standalone=False)

    args = parser.parse_args()

    try:
        if args.command == "get":
            return run_get_genomes(args)

        elif args.command == "extract":
            return run_extraction(args)
        
        elif args.command == "analyse":
            return run_analysis(args)

        elif args.command == "full":
            print("\nGETTING GENOMES")
            exit_code = run_get_genomes(args)
            if exit_code != 0:
                print("\nERROR: Genome retrieval failed.", file=sys.stderr)
                return exit_code

            genomes_out_dir = getattr(args, "input_path", None) # Set during run_get_genomes()
            if not genomes_out_dir or not os.path.isdir(genomes_out_dir):
                print("\nERROR: Genome retrieval did not produce a valid genomes directory.", file=sys.stderr)
                return 1

            print("\nEXTRACTING TERMINATORS")
            extract_args = argparse.Namespace(**vars(args)) # copy
            extract_args.input_path = genomes_out_dir
            taxon_dir = os.path.dirname(genomes_out_dir)
            extract_args.output_dir = os.path.join(taxon_dir, "terminators")

            exit_code = run_extraction(extract_args)
            if exit_code != 0:
                print("\nERROR: Terminator extraction failed.", file=sys.stderr)
                return exit_code

            print("\nANALYSING TERMINATORS")
            analyse_args = argparse.Namespace(**vars(args))
            analyse_args.input_path = extract_args.output_dir
            if not os.path.isdir(analyse_args.input_path):
                print(f"\nERROR: '{analyse_args.input_path}' not found.", file=sys.stderr)
                return 1
            analyse_args.output_dir = os.path.join(taxon_dir, "plots")

            return run_analysis(analyse_args)
        else:
            print("Invalid command. Use -h for help.", file=sys.stderr)
            return 1

    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
