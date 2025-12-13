# Terminator Analyser - AI Agent Instructions

## Project Overview
Bioinformatics pipeline for extracting and analyzing plant gene terminators (3'UTR + downstream regions). Reproduces methodology from Loke et al. (2005) to find conserved positional signals in Near-Upstream Elements (NUEs) and Cleavage Elements (CEs).

## Architecture & Data Flow
The pipeline has 4 commands orchestrated by `main.py`:
1. **get** (`get_genomes.py`): Downloads reference genomes (FASTA + GFF) via NCBI Datasets API → `out/taxons/{taxon}/genomes/`
2. **extract** (`extract.py`): Parses GFF/FASTA pairs to extract terminator sequences → `out/taxons/{taxon}/terminators/`
3. **analyse** (`analyse.py`): Counts k-mers in NUE/CE regions, ranks by delta score (peak - median) → `out/taxons/{taxon}/analysis/`
4. **full**: Chains all three steps, passing output directories between stages

**Key Files**:
- `extract.py` uses `gffutils` to create SQLite databases from GFF files (cached in `out/gff_dbs/`) for efficient feature lookup
- `analyse.py` implements k-mer counting and ranking logic with configurable parameters
- `plots.py` generates matplotlib line plots showing positional distribution of top signals (multiline plots with styled legend)
- `compare_fasta.py` is a utility for validating extracted terminators against reference datasets

**Dependencies** (from requirements.txt):
- gffutils: GFF parsing and SQLite DB creation
- pyfaidx: Indexed FASTA file access
- pandas: Data manipulation for plotting
- matplotlib: Visualization
- requests: NCBI API interactions

## Critical Coordinate System
All sequence coordinates use **cleavage site (CS) anchoring**:
- Negative: 3'UTR positions (-1 = last nt of 3'UTR)
- Positive: Downstream region (+1 = first nt after CS)
- **No position 0 exists**

Analysis windows (hardcoded in `analyse.py`):
- NUE: -50 to -1 (plotting range: -35 to -5)
- CE: -10 to +20 (plotting range: -10 to +15)

When modifying extraction or analysis logic, maintain this coordinate system.

## Module Details

### get_genomes.py
**Command**: `python main.py get "<taxon>" [--api-key KEY] [--max-genomes N] [--force] [-o OUTPUT_DIR]`

Key functions:
- `run_get_genomes()`: Main entry point, calls `_get_genomes_by_taxon()`
- `_get_genomes_by_taxon()`: Orchestrates API calls, downloads ZIPs, extracts atomically
  - Sets `args.input_path` to created genomes directory for downstream stages
  - Uses genus name if multiple genomes downloaded (e.g., "arabidopsis")
- `_fetch_all_pages()`: Paginated API queries with token-based navigation
- `_api_request()`: Single API call with retry logic (3 retries, 30s timeout, backoff 0.3)
- `_build_session()`: Creates requests.Session with retry adapter
- `_download_dataset()`: ZIP download with atomic extraction (temp file → rename)
- `_print_found_genomes()`: Formatted table output (organism, common name, accession)

NCBI API details:
- Uses `/datasets/v2/genome/taxon/{taxon}/dataset_report?filters.reference_only=true`
- Requires `genome_gff` OR `genome_gtf` availability
- Optional `NCBI_API_KEY` env var; no key = rate-limited public access
- Downloads as ZIP containing `ncbi_dataset/data/{accession}/` structure

### extract.py
**Command**: `python main.py extract "<genomes_dir>" [--raw-dna] [--filter-consecutive-a N] [--filter-window-a N] [--filter-window-size N] [--downstream-nts N] [-o OUTPUT_DIR]`

Key functions:
- `run_extraction()`: Main entry point, finds file pairs, parallelizes via Pool
- `_find_files()`: Matches FASTA/GFF pairs by basename prefix (e.g., GCF_123 from GCF_123.fna + GCF_123.gff)
- `_extract_all_terminators()`: Worker function for one genome pair
  - Creates gffutils SQLite DB in `out/gff_dbs/` (cached across runs)
  - Extracts transcripts (mRNA features) and writes to `_terminators.fa`
- `_process_transcript()`: For each transcript:
  - Finds CDS and exon features via gffutils
  - Calls `_extract_terminator()` to get sequence
  - Calls `_filter_sequence()` to check for artifacts
  - Converts to RNA (T→U) unless `--raw-dna` flag
- `_extract_terminator()`: Core extraction logic:
  - **Plus strand**: Extracts 3'UTR (exons after CDS), then downstream region → concatenate
  - **Minus strand**: Extracts downstream, then 3'UTR (exons before CDS), then reverse-complement
  - Uses 1-based GFF coords, converts to 0-based Python slices
- `_is_internal_priming_artifact()`: Filters poly-A false positives:
  - Checks downstream region for 6+ consecutive A's OR 8+ total A's (default, first 10 nts)
  - Based on Beaudoing et al. (2000) DOI: 10.1101/gr.10.7.1001
  - Only applies to RNA output (disabled if `--raw-dna`)
  - Can be fully disabled: `--filter-consecutive-a=0 --filter-window-a=0`

FASTA output format: One sequence per transcript (header = accession + transcript ID)

### analyse.py
**Command**: `python main.py analyse "<terminators_dir_or_file>" [--top-n N] [--kmer-size K] [--min-3utr-length M] [--step-size S] [--downstream-nts D] [-o OUTPUT_DIR]`

Key functions:
- `run_analysis()`: Main entry point, finds FASTA files, parallelizes counting via Pool
  - Merges counts from parallel workers into totals
  - Calls `_rank_kmers()` twice (NUE and CE regions)
  - Generates NUE/CE reports and plots
- `_process_terminator_fasta()`: Worker function for one FASTA file
  - Filters sequences by 3'UTR length (default 100 nts minimum)
  - Returns (nue_counts, ce_counts, skipped_count, total_count)
- `_count_kmers()`: Core k-mer counting logic:
  - Iterates through sequence with `step_size` (default 1 = overlapping k-mers)
  - **Position calculation**: Rightmost nucleotide of k-mer anchors position
    - `pos = (i + kmer_size - 1) - utr_len` (positions in 3'UTR as negatives)
    - `if pos >= 0: pos += 1` (skip position 0, makes downstream positive)
  - Stores as nested dict: `{kmer: {position: count, ...}, ...}`
- `_merge_counts()`: Combines per-file counts into totals
- `_rank_kmers()`: Rankings by delta score:
  - Delta = peak_count - median_positional_count
  - Returns top N k-mers as list of dicts: `[{"kmer": "AATAAA", "peak": 100, "median": 20, "delta": 80}, ...]`
  - Orders highest delta first
- `_save_report()`: Writes text reports (NUE_report.txt, CE_report.txt) with ranked k-mers

Output files: `{output_dir}/NUE_report.txt`, `CE_report.txt`, `NUE_signals_plot.png`, `CE_signals_plot.png`

### plots.py
**Function**: `plot_signal_distribution(ranked_kmers, counts_data, region, x_min, x_max, out_file)`

Generates matplotlib line plots:
- X-axis: Position relative to CS (from `x_min` to `x_max`)
- Y-axis: Count (log scale not used)
- Multiple k-mers as separate lines with styled legend
  - Linestyles: `-`, `--`, `-.`, `:`
  - Markers: `o`, `s`, `^`, `D`, `v`, `P`, `*`, `X`
  - Legend arranged in row-major order, 5 columns max
- Grid enabled, zero floor on y-axis
- Title: "Distribution of Top N Signals in {region}"
- Saves as PNG (non-interactive backend)

## Module Patterns

### Argument Handling
Each module uses dual-mode CLI pattern:
- `add_*_args(parser, is_standalone=True/False)`: Adds args conditionally
  - `is_standalone=True`: Includes `input_path` positional arg + `-o/--output-dir`
  - `is_standalone=False`: Omits positional args (used in `full` command where `main.py` manages paths)
- `run_*()` functions:
  - Take `argparse.Namespace` as sole argument
  - Return exit codes (0=success, 1=failure)
  - May modify args in-place (e.g., `get_genomes.py` sets `args.input_path`)

### Parallelization Strategy
- **extract.py**: Parallelizes per-genome extraction
  - Uses `multiprocessing.Pool()` with default workers
  - No special initialization needed (no global state across workers)
- **analyse.py**: Parallelizes per-FASTA-file k-mer counting
  - Uses `multiprocessing.Pool()` with default workers
  - Results merged sequentially after workers complete
- **get_genomes.py**: No parallelization (API calls are sequential)
- GFFutils SQLite DBs are thread-safe within a single process; avoid multi-process access to same DB

## Bioinformatics-Specific Logic

### Internal Priming Artifact Filter
Implemented in `_is_internal_priming_artifact()` (extract.py):
- Filters poly(A) false positives common in 3'UTR prediction
- Default thresholds:
  - `--filter-consecutive-a 6`: 6+ consecutive A's in window → filtered
  - `--filter-window-a 8`: 8+ total A's in window → filtered
  - `--filter-window-size 10`: Check first 10 downstream nts
- Applies only to RNA output (skipped if `--raw-dna` used)
- Disable entirely: `--filter-consecutive-a 0 --filter-window-a 0`
- Reference: Beaudoing et al. (2000) "Identification of polyadenylation signals by kernel methods"

### Strand Handling in _extract_terminator()
Handles reverse-complement for minus-strand genes:
- **Plus strand (+)**:
  1. Find 3'UTR exons (exons after CDS end)
  2. Concatenate 3'UTR parts left-to-right
  3. Append downstream region (CDS end → CDS end + 50 nts)
  4. Result is 5'→3' sense strand
- **Minus strand (-)**:
  1. Find 3'UTR exons (exons before CDS start)
  2. Prepend downstream region (CDS start - 50 → CDS start)
  3. Append 3'UTR parts right-to-left
  4. Reverse-complement entire sequence via `pyfaidx.Sequence().reverse.complement`
  5. Result is 5'→3' sense strand

Coordinate details:
- GFF uses 1-based, closed intervals
- Python slices use 0-based, half-open intervals
- Code subtracts 1 from GFF start coords before slicing: `[start-1:end]`

### K-mer Position Anchoring
In `_count_kmers()` (analyse.py):
- K-mer position is defined as **rightmost nucleotide's position**
- Example: AATAAA 6-mer at indices 5-10 in sequence → position = 10
- **Coordinate mapping**:
  - Input: sequences = [UTR (50 nts) + downstream (50 nts)]
  - For each k-mer at index i with size k:
    - `pos = (i + k - 1) - utr_len` converts to relative position
    - Result ranges from -49 (last 3'UTR position) to +49 (last downstream position)
    - Skip position 0: `if pos >= 0: pos += 1` shifts downstream to +1...+49
- This maintains the no-zero coordinate system across NUE/CE analysis

## Development Workflows

### Running the full pipeline
```bash
python main.py full "Arabidopsis thaliana" --top-n 10
```
Creates directory structure:
```
out/
  taxons/
    arabidopsis_thaliana/
      genomes/          # FASTA + GFF files
      terminators/      # Extracted terminator sequences
      analysis/         # NUE/CE reports and plots
  gff_dbs/              # SQLite databases (cached)
```

### Testing individual modules
```bash
# Download only (skip extraction/analysis)
python main.py get "Saccharomyces cerevisiae" --max-genomes 1 -o test_out

# Extract terminators from existing genomes
python main.py extract "test_out/taxons/saccharomyces_cerevisiae/genomes" -o test_terminators

# Analyse with custom k-mer size
python main.py analyse "test_terminators" --kmer-size 7 --top-n 15 -o test_analysis

# Disable internal priming filter
python main.py analyse "test_terminators" --filter-consecutive-a 0 --filter-window-a 0 -o no_filter_analysis
```

### Validating terminator sequences
Use `compare_fasta.py` utility:
```bash
python compare_fasta.py reference_terminators.fa extracted_terminators.fa
```
Reports sequence differences for validation against known datasets.

## Dependencies & External APIs

### NCBI Datasets API
- **Base URL**: `https://api.ncbi.nlm.nih.gov/datasets/v2`
- **Endpoints used**:
  - `/genome/taxon/{taxon}/dataset_report?filters.reference_only=true` (get list)
  - `/genome/accession/{accession}/download_summary` (check available files)
  - `/genome/accession/{accession}/download?include=fasta,gff3,gtf` (download ZIP)
- **Authentication**: Optional `NCBI_API_KEY` env var (public access limited to ~10 requests/min)
- **Retry logic**: 3 retries, 0.3s backoff factor, 30s timeout
- **File formats**: Downloads as ZIP with structure `ncbi_dataset/data/{accession}/` containing:
  - `*.fna` or `*.fasta`: DNA sequences
  - `*.gff` or `*.gff3` or `*.gtf`: Annotations (GFF preferred)
- Extraction is atomic (temp file → rename)

### External Libraries
- **gffutils** (v0.10+): Creates/caches SQLite DBs from GFF files, provides efficient feature queries
  - DBs stored in `out/gff_dbs/{accession}.db`
  - Auto-created if missing, reused if present (delete manually to rebuild)
- **pyfaidx** (v0.7+): Indexed FASTA access
  - Auto-creates `.fai` index files on first access
  - Provides fast sequence lookups by chromosome/coordinate
  - `Sequence` objects support `.reverse.complement` operations
- **pandas**: Data frames for organizing k-mer counts for plotting
- **matplotlib** (Agg backend): Non-interactive PNG generation

## Common Pitfalls & Troubleshooting

1. **GFF database caching**: SQLite DBs persist in `out/gff_dbs/` across runs
   - Solution: Delete DB file if source GFF updated to force rebuild
   - Normal behavior: Second run reuses DB for same accession

2. **File pair matching failures**: `extract.py` matches FASTA/GFF by basename prefix
   - Example: `GCF_000005845.2.fna` + `GCF_000005845.2.gff` match on `GCF_000005845.2`
   - Mismatches (e.g., different versions) result in silent skips
   - Solution: Ensure paired files have identical base names

3. **Coordinate fencepost errors**: Remember no position 0 when adding features
   - NUE/CE regions include negative and positive positions but never 0
   - When debugging coordinates, verify: -1 (last 3'UTR) → +1 (first downstream)

4. **RNA vs DNA output**: Default converts T→U during extraction
   - `--raw-dna` flag skips conversion, outputs DNA (T retained)
   - Internal priming filter only applies to RNA mode (disabled for `--raw-dna`)

5. **Filter interaction**: Multiple filter params combine with AND logic
   - Both conditions checked: consecutive A's AND total A's in window
   - Disable entire filter: `--filter-consecutive-a 0 --filter-window-a 0`
   - Each condition: `--filter-consecutive-a N` OR `--filter-window-a N` (either triggers)

6. **No terminators extracted**: Common causes:
   - GFF annotation lacks transcript/CDS/exon features
   - 3'UTRs too short (< 100 nts by default, set `--min-3utr-length`)
   - All filtered by internal priming detection (check with `--filter-consecutive-a 0`)

7. **Analysis produces no signals**: Likely causes:
   - Filtered terminators too few for k-mer frequency (need many copies per taxon)
   - K-mer size too large (default 6 is standard; try `--kmer-size 4` for small datasets)
   - 3'UTR length filter too strict

## Future TODOs (from code comments)
- Replace delta score with robust statistical method (e.g., information content, entropy)
- Add more visualizations (reproduce Loke Figure 1C/1D, Figure 3 patterns)
- Parallelize analysis step (currently only extract is parallel)
- Add conservation metrics across multiple taxa
- Identify unique signals per gene (vs genome/genus comparisons)
