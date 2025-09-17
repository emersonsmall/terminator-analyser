import os
import sys
import argparse
import requests
import zipfile
import io

API_BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
GENOMES_DIR = os.path.join("out", "genomes")

class GenomeRetrievalError(Exception):
    """Base exception for this script."""
    pass

def get_args(return_parser: bool = False) -> argparse.Namespace | argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Retrieve all reference genomes from NCBI Datasets API for the given taxon."
    )
    parser.add_argument(
        "--taxon",
        help="Taxon name (e.g., 'Homo sapiens' or 'Arabidopsis')."
    )
    parser.add_argument(
        "--api-key",
        default=os.environ.get("NCBI_API_KEY"),
        help="NCBI API key. Can also be set via NCBI_API_KEY environment variable."
    )
    parser.add_argument(
        "--genomes-dir",
        default=GENOMES_DIR,
        help="Parent directory to save downloaded genomes to (default: ./out/genomes)"
    )

    if return_parser:
        return parser
    
    return parser.parse_args()


def get_genomes_by_taxon(taxon: str, api_key: str, genomes_dir: str) -> str:
    """
    Finds, downloads, and extracts all reference genomes for a given taxon.
    """
    print(f"Searching for reference genomes for taxon '{taxon}'")
    report_url = f"{API_BASE_URL}/genome/taxon/{requests.utils.quote(taxon)}/dataset_report?filters.reference_only=true"
    report_res = _api_request(report_url, api_key)

    reports = report_res.get("reports")
    if not reports:
        raise GenomeRetrievalError(f"No reference genomes found for taxon '{taxon}'")

    num_genomes = len(reports)
    print(f"Found {num_genomes} reference genomes")

    taxon_dir_name = taxon.replace(" ", "_")
    taxon_path = os.path.join(genomes_dir, taxon_dir_name)
    os.makedirs(taxon_path, exist_ok=True)

    for i, report in enumerate(reports):
        accession = report.get("accession")
        if not accession:
            print(f"Skipping genome with missing accession: {report}")
            continue

        organism_name = report.get("organism").get("organism_name", "N/A")
        print(f"\nProcessing genome {i+1}/{num_genomes}: {accession} ({organism_name})")

        existing_files = os.listdir(taxon_path)
        has_fasta = any(f.startswith(accession) and f.endswith((".fna", ".fa", ".fasta")) for f in existing_files)
        has_gff = any(f.startswith(accession) and f.endswith(".gff") for f in existing_files)

        if has_fasta and has_gff:
            print(f"Genome {accession} already downloaded. Skipping.")
            continue

        download_url = f"{API_BASE_URL}/genome/accession/{accession}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF"
        _download_and_unzip(download_url, api_key, taxon_path, accession)
    
    return taxon_path


def run_get_genomes(args: argparse.Namespace) -> int:
    if args.taxon:
        get_genomes_by_taxon(args.taxon, args.api_key, args.genomes_dir)
        return 0
    
    elif hasattr(args, 'input_path') and args.input_path:
        if not os.path.isdir(args.input_path):
            raise FileNotFoundError(f"Input path '{args.input_path}' does not exist or is not a directory.")
        print(f"Using local genome files from: '{args.input_path}'")
        return 0
    else:
        raise ValueError("An input source is required.")


# Helper functions
def _api_request(url: str, api_key: str) -> dict:
    headers = {"Accept": "application/json", "api-key": api_key} if api_key else {"Accept": "application/json"}
    try:
        res = requests.get(url, headers=headers)
        res.raise_for_status()
        return res.json()
    except requests.RequestException as e:
        raise GenomeRetrievalError(f"API request failed for URL {url}: {e}")


def _download_and_unzip(url: str, api_key: str, out_dir: str, accession: str) -> None:
    headers = {"api-key": api_key} if api_key else {}
    print(f"Downloading files to '{out_dir}'")
    try:
        with requests.get(url, headers=headers, stream=True) as res:
            res.raise_for_status()
            with zipfile.ZipFile(io.BytesIO(res.content)) as z:
                for file_info in z.infolist():
                    filename = os.path.basename(file_info.filename)

                    if filename.endswith(('.fna', '.fa', '.fasta')):
                        with z.open(file_info) as src, open(os.path.join(out_dir, f"{accession}.fna"), 'wb') as dest:
                            dest.write(src.read())
                            
                    elif filename.endswith('.gff'):
                        with z.open(file_info) as src, open(os.path.join(out_dir, f"{accession}.gff"), 'wb') as dest:
                            dest.write(src.read())
                    
    except requests.RequestException as e:
        raise GenomeRetrievalError(f"Download failed for URL {url}: {e}")
    except zipfile.BadZipFile as e:
        raise GenomeRetrievalError(f"Failed to unzip downloaded content from {url}: {e}")


def main() -> int:
    """Standalone execution entry point."""
    args = get_args()
    try:
        return run_get_genomes(args)
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
