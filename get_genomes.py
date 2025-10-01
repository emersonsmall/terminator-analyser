import os
import sys
import argparse
import requests
import zipfile
from typing import Optional
import tempfile
import shutil
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

GENBANK_API_BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
RETRIES = 3
TIMEOUT_IN_SECONDS = 30
CHUNK_SIZE = 32 * 1024 # 32KB
BACKOFF_FACTOR = 0.3

class GenomeRetrievalError(Exception):
    """Base exception for this script."""
    pass


def add_args_to_parser(parser: argparse.ArgumentParser, standalone: bool = True) -> None:
    """Adds genome retrieval arguments to the given parser."""
    parser.add_argument(
        "taxon",
        help="Taxon name (e.g., 'Arabidopsis' or 'Arabidopsis thaliana')."
    )
    parser.add_argument(
        "--api-key",
        default=os.environ.get("NCBI_API_KEY"),
        help="NCBI API key. Can be set via NCBI_API_KEY environment variable."
    )
    parser.add_argument(
        "--max-genomes",
        type=int,
        default=None,
        help="Maximum number of genomes to download (default: None)."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing files with the same filenames."
    )

    if standalone:
        parser.add_argument(
            "-o",
            "--output-dir",
            default="out",
            help="Path to the output directory (default: /out)."
        )


def _get_args() -> argparse.Namespace:
    """Gets arguments for standalone script execution."""
    parser = argparse.ArgumentParser(
        description="Gets FASTA and GFF files for the given taxon using the Genbank API"
    )
    add_args_to_parser(parser)
    return parser.parse_args()


def get_genomes_by_taxon(
        taxon: str,
        api_key: str,
        output_dir: str,
        session: requests.Session,
        max_genomes: Optional[int] = None,
        force: bool = False,
    ) -> str:
    """
    Finds, downloads, and extracts all reference genomes for a given taxon.
    """
    taxon = taxon.strip().lower()
    print(f"Searching for reference genomes for taxon '{taxon}'")
    report_url = f"{GENBANK_API_BASE_URL}/genome/taxon/{requests.utils.quote(taxon)}/dataset_report?filters.reference_only=true"
    report_res = _api_request(session, report_url, api_key)
    reports = report_res.get("reports")
    if not reports:
        raise GenomeRetrievalError(f"No reference genomes found for taxon '{taxon}'")

    if max_genomes is not None:
        reports = reports[:max_genomes]

    num_genomes = len(reports)
    print(f"\nFound {num_genomes} reference genomes:")
    print("Organism name | Common name | Accession")
    for report in reports:
        accession = report.get("accession", "N/A")
        organism_name = report.get("organism", {}).get("organism_name", "N/A")
        common_name = report.get("organism", {}).get("common_name", "N/A")
        print(f"{organism_name} | {common_name} | {accession}")

    taxon_dir_name = taxon.replace(" ", "_")
    genomes_dir_path = os.path.join(output_dir, "taxons", taxon_dir_name, "genomes")
    os.makedirs(genomes_dir_path, exist_ok=True)

    for i, report in enumerate(reports, start=1):
        accession = report.get("accession")
        if not accession:
            print(f"Skipping genome with missing accession: {report}")
            continue

        organism_name = report.get("organism").get("organism_name", "N/A")

        existing_files = os.listdir(genomes_dir_path)
        has_fasta = any(f.startswith(accession) and f.endswith((".fna", ".fa", ".fasta")) for f in existing_files)
        has_gff = any(f.startswith(accession) and f.endswith((".gff", ".gff3")) for f in existing_files)

        if has_fasta and has_gff:
            print(f"\n{accession} genome already downloaded. Skipping.")
            continue

        print(f"\nDownloading genome {i}/{num_genomes}: {accession} ({organism_name})")
        download_url = f"{GENBANK_API_BASE_URL}/genome/accession/{accession}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF"
        _download_and_unzip(session, download_url, api_key, genomes_dir_path, accession, force=force)
    
    return genomes_dir_path


def run_get_genomes(args: argparse.Namespace) -> int:
    if not args.taxon:
        print("ERROR: Taxon name is required to download genomes.", file=sys.stderr)
        return 1

    session = _build_session()
    try:
        try:
            taxon_path = get_genomes_by_taxon(
                taxon=args.taxon,
                api_key=args.api_key,
                output_dir=args.output_dir,
                session=session,
                max_genomes=args.max_genomes,
                force=args.force,
            )

            args.input_path = taxon_path
            print(f"\nDownloaded genomes saved to: {taxon_path}")
            return 0
        except GenomeRetrievalError as e:
            print(f"ERROR: {e}", file=sys.stderr)
            return 1
    finally:
        try:
            session.close()
        except Exception:
            pass


# Helper functions
def _build_session() -> requests.Session:
    """
    Return a requests.Session configured with retry logic.
    Retries restricted to GET
    Caller responsible for calling session.close()
    """
    s = requests.Session()

    s.headers.update({
        "Accept": "application/zip, application/octet-stream, */*"
    })
    
    # only retry for GET
    allowed = frozenset(['GET'])
    try:
        retry = Retry(
            total=RETRIES,
            backoff_factor=BACKOFF_FACTOR,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=allowed,
            raise_on_status=False
        )
    except TypeError:
        # older urllib3 versions use method_whitelist
        retry = Retry(
            total=RETRIES,
            backoff_factor=BACKOFF_FACTOR,
            status_forcelist=[429, 500, 502, 503, 504],
            method_whitelist=allowed,
            raise_on_status=False
        )
    
    adapter = HTTPAdapter(max_retries=retry)
    s.mount("https://", adapter)
    s.mount("http://", adapter)

    s.request_timeout = TIMEOUT_IN_SECONDS
    return s


def _api_request(session: requests.Session, url: str, api_key: str) -> dict:
    headers = {"Accept": "application/json"}
    if api_key:
        headers["api-key"] = api_key
    try:
        res = session.get(url, headers=headers, timeout=session.request_timeout)
        res.raise_for_status()
        return res.json()
    except requests.RequestException as e:
        raise GenomeRetrievalError(f"API request failed for URL {url}: {e}")


def _download_and_unzip(
        session: requests.Session,
        url: str, 
        api_key: str | None,
        out_dir: str,
        accession: str,
        force: bool = False,
    ) -> None:
    """
    Stream-download the zip at `url` into a temp file on disk, then open it with zipfile.ZipFile
    and extract FASTA/GFF files to `out_dir` with filenames prefixed by `accession`.

    Uses session.
    Writes downloaded ZIP to temp file to avoid memory issues with large files.
    Writes each extracted file to a temp file then os.replace() to final destination.
    """
    headers = {}
    if api_key:
        headers["api-key"] = api_key

    os.makedirs(out_dir, exist_ok=True)

    tmp_zip = None
    try:
        with session.get(url, headers=headers, stream=True, timeout=session.request_timeout) as res:
            res.raise_for_status()
            tmp_fd, tmp_zip = tempfile.mkstemp(prefix="datasets_download_", suffix=".zip", dir=out_dir)
            os.close(tmp_fd)
            with open(tmp_zip, 'wb') as f:
                for chunk in res.iter_content(chunk_size=CHUNK_SIZE):
                    if chunk:
                        f.write(chunk)
            
            # open zipfile from disk
            with zipfile.ZipFile(tmp_zip, 'r') as z:
                for member in z.infolist():
                    filename = os.path.basename(member.filename)
                    if not filename:
                        continue
                    lower = filename.lower()
                    if lower.endswith((".fna", ".fa", ".fasta")):
                        out_fname = f"{accession}.fna"
                    elif lower.endswith((".gff", ".gff3")):
                        out_fname = f"{accession}.gff"
                    else:
                        continue

                    out_path = os.path.join(out_dir, out_fname)
                    if os.path.exists(out_path) and not force:
                        print(f"File {out_fname} already exists. Skipping.")
                        continue
                        
                    # extract member -> temp file then atomic replace
                    with z.open(member) as src:
                        with tempfile.NamedTemporaryFile(delete=False, dir=out_dir) as tmpf:
                            shutil.copyfileobj(src, tmpf)
                            temp_name = tmpf.name
                        os.replace(temp_name, out_path)
    finally:
        if tmp_zip and os.path.exists(tmp_zip):
            try:
                os.remove(tmp_zip)
            except Exception:
                pass   


def main():
    """Standalone execution entry point."""
    try:
        return run_get_genomes(_get_args())
    except Exception as e:
        print(f"\nERROR: {e}", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())
