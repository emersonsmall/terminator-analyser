# Built-in libraries
import os
import sys
import argparse
import requests
import zipfile
import tempfile
import shutil
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from urllib.parse import quote


NCBI_DATASETS_API_BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2"
TAXON_BASE_URL = f"{NCBI_DATASETS_API_BASE_URL}/genome/taxon"
ACCESSION_BASE_URL = f"{NCBI_DATASETS_API_BASE_URL}/genome/accession"
DATASET_REPORT_FILTERS = ("filters.reference_only=true",)

RETRIES = 3
TIMEOUT_IN_SECONDS = 30
CHUNK_SIZE = 32 * 1024
BACKOFF_FACTOR = 0.3

VALID_FASTA_EXTS = (".fna", ".fa", ".fasta")
VALID_ANNOTATION_EXTS = (".gff", ".gff3", ".gtf")


def run_get_genomes(args: argparse.Namespace) -> int:
    try:
        taxon_path = _get_genomes_by_taxon(
            args.taxon,
            args.output_dir,
            args.api_key,
            args.max_genomes,
            args.force,
        )

        args.input_path = taxon_path  # set input_path for downstream use
        return 0
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


def add_get_args(parser: argparse.ArgumentParser, is_standalone: bool = True) -> None:
    """Adds command-line arguments for the `get` command to the given parser.

    Args:
        parser: The argument parser to add the arguments to.
        standalone: Includes standalone execution arguments if True. Default: True.
    """

    parser.add_argument(
        "taxon", 
        type=lambda s: s.lower().strip(),
        help="Taxon name (e.g., 'Arabidopsis' or 'Arabidopsis thaliana')."
    )
    parser.add_argument(
        "--api-key",
        default=os.environ.get("NCBI_API_KEY"),
        help="NCBI API key. Can be set via NCBI_API_KEY environment variable.",
    )
    parser.add_argument(
        "--max-genomes",
        type=int,
        default=None,
        help="Maximum number of genomes to download (default: None).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing files with the same filenames.",
    )

    if is_standalone:
        parser.add_argument(
            "-o",
            "--output-dir",
            default="out",
            help="Path to the output directory (default: ./out).",
        )


def _get_args() -> argparse.Namespace:
    """Gets arguments for standalone script execution."""

    parser = argparse.ArgumentParser(
        description="Gets FASTA and annotation files for the given taxon using the NCBI Datasets API"
    )
    add_get_args(parser, is_standalone=True)
    return parser.parse_args()


def _get_genomes_by_taxon(
    taxon: str,
    output_dir: str,
    api_key: str | None = None,
    max_genomes: int | None = None,
    force: bool = False,
) -> str:
    """
    Downloads and extracts all reference genomes for the given taxon.
    """

    session = _build_session()
    try:
        dataset_report_url = f"{TAXON_BASE_URL}/{quote(taxon)}/dataset_report?{','.join(DATASET_REPORT_FILTERS)}"
        dataset_reports = _fetch_all_pages(
            session, dataset_report_url, api_key, max_genomes
        )

        if not dataset_reports:
            raise Exception(f"No reference genomes found for taxon '{taxon}'")

        print(f"\nFound {len(dataset_reports)} reference genome/s for '{taxon}':")
        _print_found_genomes(dataset_reports)

        num_genomes = len(dataset_reports)
        taxon_dir_name = dataset_reports[0].get("organism").get("organism_name")
        if num_genomes > 1:
            taxon_dir_name = taxon_dir_name.split(" ")[0]  # use genus for dir name
        taxon_dir_name = taxon_dir_name.lower().replace(" ", "_")

        genomes_dir_path = os.path.join(output_dir, "taxons", taxon_dir_name, "genomes")
        os.makedirs(genomes_dir_path, exist_ok=True)

        num_downloaded = 0
        for i, report in enumerate(dataset_reports, start=1):
            accession = report.get("accession")

            if not force and _genome_files_exist(genomes_dir_path, accession):
                print(f"{accession} already exists at '{genomes_dir_path}', skipping")
                continue

            organism_name = report.get("organism").get("organism_name", "N/A")
            print(
                f"\nDownloading genome {i}/{num_genomes}: {organism_name} ({accession})"
            )

            # check if GFF or GTF annotation is available. GFF preferred
            download_summary_url = (
                f"{ACCESSION_BASE_URL}/{quote(accession)}/download_summary"
            )
            download_summary = _api_request(session, download_summary_url, api_key)

            available_files = download_summary["available_files"]
            has_gff = available_files.get("genome_gff", False)
            has_gtf = available_files.get("genome_gtf", False)

            if not (has_gff or has_gtf):
                print(f"{accession} has no annotation available, skipping")
                continue

            files_to_request = ["GENOME_FASTA"]
            if has_gff:
                files_to_request.append("GENOME_GFF")
            elif has_gtf:
                files_to_request.append("GENOME_GTF")

            download_url = f"{ACCESSION_BASE_URL}/{quote(accession)}/download?include_annotation_type={','.join(files_to_request)}"
            _download_and_extract(
                session, download_url, api_key, genomes_dir_path, accession
            )
            num_downloaded += 1

        if num_downloaded > 0:
            print(f"\nDownloaded {num_downloaded} genome/s to: {genomes_dir_path}")

        return genomes_dir_path

    finally:
        session.close()


# --- HELPER FUNCTIONS ---
def _genome_files_exist(dir_path: str, accession: str) -> bool:
    """Checks if both FASTA and annotation files for the given accession."""

    existing_files = os.listdir(dir_path)

    has_fasta = any(
        f.startswith(accession) and f.endswith(VALID_FASTA_EXTS)
        for f in existing_files
    )

    has_annotation = any(
        f.startswith(accession) and f.endswith(VALID_ANNOTATION_EXTS)
        for f in existing_files
    )

    return has_fasta and has_annotation


def _fetch_all_pages(
    session: requests.Session,
    base_url: str,
    api_key: str | None = None,
    max_genomes: int | None = None,
) -> list:
    """
    Fetches all pages of genomes reports for the NCBI Datasets API taxon report endpoint.

    Args:
        session: The requests.Session object to use for the API requests.
        base_url: The base URL of the API endpoint to fetch the reports from.
        api_key: The API key to use for the API requests.
        max_genomes: The maximum number of genomes to fetch.

    Returns:
        A list of all genome reports, optionally limited to the first `max_genomes` reports.
    """

    all_reports = []
    next_page_token = None

    while True:
        if next_page_token:
            url = f"{base_url}&page_token={quote(next_page_token)}"
        else:
            url = base_url

        res = _api_request(session, url, api_key)
        page_reports = res.get("reports", [])
        all_reports.extend(page_reports)

        if max_genomes is not None and len(all_reports) >= max_genomes:
            all_reports = all_reports[:max_genomes]
            break

        next_page_token = res.get("next_page_token", None)
        if not next_page_token:
            break

    return all_reports


def _print_found_genomes(reports: list) -> None:
    # Prepare data and find max widths for column alignment
    header = {
        "organism": "Organism name",
        "common": "Common name",
        "accession": "Accession",
    }
    data = []
    max_widths = {key: len(value) for key, value in header.items()}

    for report in reports:
        organism_name = report.get("organism", {}).get("organism_name", "N/A")
        common_name = report.get("organism", {}).get("common_name", "N/A")
        accession = report.get("accession")

        data.append(
            {"organism": organism_name, "common": common_name, "accession": accession}
        )

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

    s.headers.update({"Accept": "application/zip, application/octet-stream, */*"})

    # only retry for GET
    allowed = frozenset(["GET"])
    try:
        retry = Retry(
            total=RETRIES,
            backoff_factor=BACKOFF_FACTOR,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=allowed,
            raise_on_status=False,
        )
    except TypeError:
        # older urllib3 versions use method_whitelist
        retry = Retry(
            total=RETRIES,
            backoff_factor=BACKOFF_FACTOR,
            status_forcelist=[429, 500, 502, 503, 504],
            method_whitelist=allowed,
            raise_on_status=False,
        )

    adapter = HTTPAdapter(max_retries=retry)
    s.mount("https://", adapter)
    s.mount("http://", adapter)

    return s


def _api_request(
    session: requests.Session, url: str, api_key: str | None = None
) -> dict:
    headers = {"Accept": "application/json"}
    if api_key:
        headers["api-key"] = api_key
    try:
        res = session.get(url, headers=headers, timeout=TIMEOUT_IN_SECONDS)
        res.raise_for_status()
        return res.json()
    except requests.RequestException as e:
        raise Exception(f"API request failed for URL {url}: {e}")


def _download_and_extract(
    session: requests.Session,
    url: str,
    api_key: str | None,
    out_dir: str,
    accession: str,
) -> None:
    """
    Stream-download the zip at `url` into a temp file, then open it with zipfile.ZipFile
    and extract relevant files to `out_dir` with filenames prefixed by `accession`.
    """

    headers = {}
    if api_key:
        headers["api-key"] = api_key

    os.makedirs(out_dir, exist_ok=True)

    tmp_zip = None
    try:
        with session.get(
            url, headers=headers, stream=True, timeout=TIMEOUT_IN_SECONDS
        ) as res:
            res.raise_for_status()
            tmp_fd, tmp_zip = tempfile.mkstemp(
                prefix="datasets_download_", suffix=".zip", dir=out_dir
            )
            os.close(tmp_fd)
            with open(tmp_zip, "wb") as f:
                for chunk in res.iter_content(chunk_size=CHUNK_SIZE):
                    if chunk:
                        f.write(chunk)

            # open zipfile from disk
            relevant_exts = VALID_FASTA_EXTS + VALID_ANNOTATION_EXTS
            with zipfile.ZipFile(tmp_zip, "r") as z:
                for member in z.infolist():
                    base_name = os.path.basename(member.filename)
                    _, ext = os.path.splitext(base_name.lower())

                    if ext not in relevant_exts:
                        continue

                    out_fname = f"{accession}{ext}"
                    out_path = os.path.join(out_dir, out_fname)

                    # extract member to temp file then atomic replace
                    with z.open(member) as src:
                        with tempfile.NamedTemporaryFile(
                            delete=False, dir=out_dir
                        ) as tmpf:
                            shutil.copyfileobj(src, tmpf)
                            temp_name = tmpf.name
                        os.replace(temp_name, out_path)

    finally:
        if tmp_zip and os.path.exists(tmp_zip):
            try:
                os.remove(tmp_zip)
            except Exception:
                pass


# --- STANDALONE EXECUTION ---
def main():
    return run_get_genomes(_get_args())


if __name__ == "__main__":
    sys.exit(main())
