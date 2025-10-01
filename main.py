import argparse
import sys
import os

from extract_terminators import run_extraction, _get_args as get_extract_args
from analyse_terminators import run_analysis, _get_args as get_analyse_args
from get_genomes import run_get_genomes, _get_args as get_genomes_args

def main() -> int:
    """
    Main function for the terminator analysis pipeline.
    Uses subparsers to handle different stages: getting genomes, extracting terminators, 
    analysing terminators, or full pipeline.
    """
    parser = argparse.ArgumentParser(
        description="Terminator Analysis Pipeline"
    )
    subparsers = parser.add_subparsers(dest="command")
    subparsers.add_parser("get", help="Get all reference genomes for a given taxon.")
    subparsers.add_parser("extract", help="Extract terminator sequences from FASTA and GFF files.")
    subparsers.add_parser("analyse", help="Analyse NUEs and CEs of terminator sequences.")
    subparsers.add_parser("full", help="Runs the full pipeline: get, extract, and analyse.")

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
            # Create unified parser that accepts all args from all steps
            full_parser = argparse.ArgumentParser(
                description="Runs the full pipeline: get, extract, and analyse"
            )

            full_parser.add_argument(
                "-o", "--output-dir", default="out",
                help="Path to the output directory (default: ./out)."
            )

            added_dests = {"help", "output_dir", "input_path"}

            all_parsers = [
                get_genomes_args(return_parser=True),
                get_extract_args(return_parser=True),
                get_analyse_args(return_parser=True)
            ]

            for p in all_parsers:
                for action in p._actions:
                    if action.dest not in added_dests:
                        full_parser.add_argument(*action.option_strings, dest=action.dest, default=action.default, help=action.help, type=action.type)
                        added_dests.add(action.dest)
            
            args = full_parser.parse_args(remaining) # Parse all args

            # Step 1: Get Genomes
            # Uses 'taxon' from args and then adds 'input_path' for next step
            print("\nSTEP 1: GET GENOMES")
            genome_args = argparse.Namespace(**vars(args))
            exit_code = run_get_genomes(genome_args)
            if exit_code != 0:
                print("\nERROR: Genome retrieval step failed. Aborting", file=sys.stderr)
                return exit_code

            # Step 2: Extract Terminators
            # args now has 'input_path' from previous step pointing to genomes dir
            print("\nSTEP 2: EXTRACT TERMINATORS")
            extract_args = argparse.Namespace(**vars(args))
            extract_args.input_path = getattr(genome_args, 'input_path', None)
            if not extract_args.input_path:
                print("\nERROR: No input path from genome retrieval step. Aborting", file=sys.stderr)
                return 1
            
            taxon_dir = os.path.dirname(extract_args.input_path)
            extract_args.output_dir = os.path.join(taxon_dir, "terminators")

            exit_code = run_extraction(extract_args)
            if exit_code != 0:
                print("\nERROR: Terminator extraction step failed. Aborting", file=sys.stderr)
                return exit_code

            # Step 3: Analyse Terminators
            # input_path for this step is output from previous step
            print("\nSTEP 3: ANALYSE TERMINATORS")
            analyse_args = argparse.Namespace(**vars(args))
            analyse_args.input_path = extract_args.output_dir
            if not os.path.isdir(analyse_args.input_path):
                print(f"\nERROR: Output dir '{analyse_args.input_path}' from extraction not found. Aborting", file=sys.stderr)
                return 1
            
            taxon_dir = os.path.dirname(extract_args.output_dir)
            analyse_args.output_dir = os.path.join(taxon_dir, "plots")

            return run_analysis(analyse_args)

    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        return 1

    print("\nFINISHED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
