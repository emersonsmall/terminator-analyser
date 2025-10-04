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


def main():
    return run_get_genomes(_get_args())

if __name__ == "__main__":
    sys.exit(main())


def run_get_genomes(args: argparse.Namespace) -> int:
    try:
        taxon_path = get_genomes_by_taxon(
            args.taxon,
            args.api_key,
            args.output_dir,
            args.max_genomes,
            args.force,
        )

        args.input_path = taxon_path # set input_path for downstream use
        return 0
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


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
            help="Path to the output directory (default: ./out)."
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
        max_genomes: Optional[int] = None,
        force: bool = False,
    ) -> str:
    """
    Downloads and extracts all reference genomes (FASTA and GFF files) for the given taxon.
    """
    session = _build_session()
    try:
        taxon = taxon.strip()
        report_url = f"{GENBANK_API_BASE_URL}/genome/taxon/{requests.utils.quote(taxon)}/dataset_report?filters.reference_only=true"
        report_res = _api_request(session, report_url, api_key)
        reports = report_res.get("reports")

        if not reports:
            raise Exception(f"No reference genomes found for taxon '{taxon}'")

        if max_genomes is not None:
            reports = reports[:max_genomes]

        _print_found_genomes(reports)

        num_genomes = len(reports)

        taxon_dir_name = reports[0].get("organism").get("organism_name")
        if num_genomes > 1:
            taxon_dir_name = taxon_dir_name.split(" ")[0] # use genus for dir name
        taxon_dir_name = taxon_dir_name.lower().replace(" ", "_")
        
        genomes_dir_path = os.path.join(output_dir, "taxons", taxon_dir_name, "genomes")
        os.makedirs(genomes_dir_path, exist_ok=True)

        num_downloaded = 0
        for i, report in enumerate(reports, start=1):
            accession = report.get("accession")

            existing_files = os.listdir(genomes_dir_path)
            has_fasta = any(f.startswith(accession) and f.endswith((".fna", ".fa", ".fasta")) for f in existing_files)
            has_gff = any(f.startswith(accession) and f.endswith((".gff", ".gff3")) for f in existing_files)

            if has_fasta and has_gff and not force:
                print(f"{accession} already exists at '{genomes_dir_path}', skipping")
                continue
            
            organism_name = report.get("organism").get("organism_name", "N/A")
            print(f"\nDownloading genome {i}/{num_genomes}: {organism_name} ({accession})")
            download_url = f"{GENBANK_API_BASE_URL}/genome/accession/{accession}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF"
            _download_and_unzip(session, download_url, api_key, genomes_dir_path, accession)
            num_downloaded += 1
        
        if num_downloaded > 0:
            print(f"\nDownloaded {num_downloaded} genome/s to: {genomes_dir_path}")

        return genomes_dir_path
    
    finally:
        session.close()


# --- HELPER FUNCTIONS ---

def _print_found_genomes(reports: list) -> None:
    num_genomes = len(reports)
    print(f"\nFound {num_genomes} reference genome{'s' if num_genomes != 1 else ''}:")

    # Prepare data and find max widths for column alignment
    header = {"organism": "Organism name", "common": "Common name", "accession": "Accession"}
    data = []
    max_widths = {key: len(value) for key, value in header.items()}

    for report in reports:
        organism_name = report.get("organism", {}).get("organism_name", "N/A")
        common_name = report.get("organism", {}).get("common_name", "N/A")
        accession = report.get("accession", "N/A")
        
        data.append({"organism": organism_name, "common": common_name, "accession": accession})
        
        max_widths["organism"] = max(max_widths["organism"], len(organism_name))
        max_widths["common"] = max(max_widths["common"], len(common_name))
        max_widths["accession"] = max(max_widths["accession"], len(accession))

    # Create a format string based on the calculated max widths
    fmt_string = (
        f"{{organism:<{max_widths['organism']}}} | "
        f"{{common:<{max_widths['common']}}} | "
        f"{{accession:<{max_widths['accession']}}}"
    )
    
    separator = (
        f"{'-' * max_widths['organism']} | "
        f"{'-' * max_widths['common']} | "
        f"{'-' * max_widths['accession']}"
    )

    # Print the formatted report
    print(fmt_string.format(**header))
    print(separator)
    for row in data:
        print(fmt_string.format(**row))
    print()

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
        raise Exception(f"API request failed for URL {url}: {e}")

def _download_and_unzip(
        session: requests.Session,
        url: str,
        api_key: str | None,
        out_dir: str,
        accession: str
    ) -> None:
    """
    Stream-download the zip at `url` into a temp file on disk, then open it with zipfile.ZipFile
    and extract FASTA/GFF files to `out_dir` with filenames prefixed by `accession`.

    Writes downloaded ZIP to temp file.
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

                    # extract member -> temp file then atomic replace
                    out_path = os.path.join(out_dir, out_fname)
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
