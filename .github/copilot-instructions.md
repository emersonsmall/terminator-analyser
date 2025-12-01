# Terminator Analyser - AI Agent Instructions

## Project Overview
Bioinformatics pipeline for extracting and analyzing plant gene terminators (3'UTR + downstream regions). Reproduces methodology from Loke et al. (2005) to find conserved positional signals in Near-Upstream Elements (NUEs) and Cleavage Elements (CEs).

## Architecture & Data Flow
The pipeline has 4 commands orchestrated by `main.py`:
1. **get** (`get_genomes.py`): Downloads reference genomes (FASTA + GFF) via NCBI Datasets API → `out/taxons/{taxon}/genomes/`
2. **extract** (`extract.py`): Parses GFF/FASTA pairs to extract terminator sequences → `out/taxons/{taxon}/terminators/`
3. **analyse** (`analyse.py`): Counts k-mers in NUE/CE regions, ranks by delta score (peak - median) → `out/taxons/{taxon}/plots/`
4. **full**: Chains all three steps, passing output directories between stages

**Key Files**:
- `extract.py` uses `gffutils` to create SQLite databases from GFF files (cached in `out/gff_dbs/`) for efficient feature lookup
- `plots.py` generates matplotlib line plots showing positional distribution of top signals

## Critical Coordinate System
All sequence coordinates use **cleavage site (CS) anchoring**:
- Negative: 3'UTR positions (-1 = last nt of 3'UTR)
- Positive: Downstream region (+1 = first nt after CS)
- **No position 0 exists**

Analysis windows (hardcoded in `analyse.py`):
- NUE: -50 to -1 (plots: -35 to -5)
- CE: -10 to +20 (plots: -10 to +15)

When modifying extraction or analysis logic, maintain this coordinate system.

## Module Patterns

### Argument Handling
Each module has dual-mode CLI:
- `add_*_args(parser, standalone=True/False)`: Adds args conditionally
- `standalone=True`: Includes `input_path` and `-o/--output-dir`
- `standalone=False`: Omits those (used in `full` command where `main.py` manages paths)
- `run_*` functions take `argparse.Namespace`, return exit codes (0=success, 1=failure)

### Parallelization
Only `extract.py` uses multiprocessing:
```python
with Pool(initializer=_init_worker, initargs=(args,)) as pool:
    pool.map(_worker, tasks)
```
Global `_worker_args` passes CLI args to workers. Add new parallelization cautiously due to gffutils DB constraints.

## Bioinformatics-Specific Logic

### Internal Priming Artifact Filter
`_is_internal_priming_artifact()` in `extract.py` filters poly-A false positives:
- Default: 6+ consecutive A's or 8+ total A's in first 10 downstream nts
- Based on Beaudoing et al. (2000)
- Only applies to RNA output (not `--raw-dna`)

### Strand Handling
`_extract_terminator()` handles reverse complement for minus-strand genes:
- Plus strand: UTR parts → downstream (5' to 3')
- Minus strand: downstream → UTR parts, then reverse complement
- Uses `pyfaidx.Sequence.reverse.complement`

### K-mer Counting
`_count_kmers()` in `analyse.py`:
- Anchors k-mer position to **rightmost nucleotide** (e.g., for AAUAAA, position = 'A' at end)
- `step_size=1` = overlapping k-mers (default)
- Adjusts coordinates: `pos = (i + kmer_size - 1) - utr_len`, then `if pos >= 0: pos += 1`

## Development Workflows

### Running the full pipeline
```bash
python main.py full "Arabidopsis thaliana" --top-n 10
```
Creates `out/taxons/arabidopsis_thaliana/{genomes,terminators,plots}/`

### Testing individual modules
```bash
# Test extraction only (skip NCBI download)
python main.py extract "path/to/genomes" -o "test_output"

# Test analysis with custom params
python main.py analyse "path/to/terminators" --kmer-size 7 --top-n 15
```

### Validating terminator sequences
Use `compare_fasta.py` to diff extracted terminators against reference datasets:
```bash
python compare_fasta.py reference.fa extracted_terminators.fa
```

## Dependencies & External APIs
- **NCBI Datasets API**: Uses paginated endpoints with retry logic (3 retries, 30s timeout)
  - Optional API key via `NCBI_API_KEY` env var
  - Downloads as ZIP, extracts FASTA/GFF atomically via temp files
- **gffutils**: Creates on-disk SQLite DBs for GFF annotation queries
- **pyfaidx**: Fast indexed FASTA access (auto-creates `.fai` indexes)

## Common Pitfalls
1. **GFF database caching**: DBs persist in `out/gff_dbs/` across runs. Delete manually if GFF updated.
2. **File pair matching**: `extract.py` matches FASTA/GFF by basename (e.g., `GCF_123.fna` + `GCF_123.gff`). Mismatches = skipped.
3. **Coordinate fencepost errors**: Remember no position 0 when adding features that span CS.
4. **RNA vs DNA**: Default output is RNA (T→U). Use `--raw-dna` to skip conversion.
5. **Filter interaction**: `--filter-consecutive-a=0 --filter-window-a=0` disables artifact filtering entirely.

## Future TODOs (from code comments)
- Replace delta score with robust statistical method
- Add more visualizations (Loke Figure 1C/1D, Figure 3)
- Parallelize analysis (currently only extract is parallel)
- Add conservation metrics
- Identify unique signals per gene (vs genome/genus comparisons)
