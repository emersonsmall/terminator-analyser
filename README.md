# Terminator Analysis Pipeline
This project is a bioinformatics pipeline designed to extract and analyse plant gene terminators. The pipeline fetches genomes, extracts terminator sequences (3' UTR + downstream region), and analyses these sequences to find conserved positional signals (patterns within the Near-Upstream Elements (NUEs) and Cleavage Elements (CEs)).

# Features
Genome Retrieval: Downloads reference genomes (FASTA and GFF files) using the NCBI Datasets API for any taxon.

Terminator Extraction: Parses GFF files to extract terminator sequences.

Filtering: Includes filters to remove likely internal priming artefacts.

Signal Analysis: Scans the terminators for conserved k-mers in the NUE and CE regions.

Visualisation: Generates plots showing the positional distribution of the top-ranked signals.

# Installation
## 1. Clone the repository
```bash
git clone https://github.com/emersonsmall/terminator-analyser  
cd terminator-analyser  
```

## 2. Create and activate a virtual environment
```bash
python -m venv venv
source ./venv/bin/activate # On Windows, use .\venv\Scripts\activate
```

## 3. Install the required packages
```bash
pip install -r requirements.txt
```

## 4. (Optional) Set your NCBI API Key
```bash
export NCBI_API_KEY="your_api_key_here"
```

# Usage
The pipeline is controlled through main.py and is divided into four main commands: get, extract, analyse, and full. You can view all options for a command by using the -h flag (e.g., `python main.py full -h`).

## full - Run the End-to-End Pipeline
This is the command for a complete analysis. It fetches genomes, extracts terminators, and runs the analysis.

Example:
Download the reference genome for Arabidopsis thaliana, extract its terminators, and analyse them to find the top 10 signals within the NUE and CE.

```bash
python main.py full "Arabidopsis thaliana" --top-n 10
```

This will create a './out' directory containing the downloaded genomes, terminator sequences, and analysis plots.

## get - Download Genomes
Downloads reference genome FASTA and GFF files for a given taxon.

Example:

```bash
python main.py get "Saccharomyces cerevisiae"
```

This will download all reference genomes for the specified taxon into the 'out/taxons/saccharomyces_cerevisiae/genomes' directory.

## extract - Extract Terminators
Extracts terminator sequences from the specified genome files.

Example:
The input path should be a directory containing FASTA/GFF pairs.

```bash
python main.py extract "out/taxons/saccharomyces_cerevisiae/genomes"
```

This will create FASTA files containing the terminator sequences in the 'out/terminators' directory.

## analyse - Analyse Terminators
Analyses the specified terminator sequences to find conserved signals and generate plots.

Example:
The input path should be a directory containing terminator FASTA files.

```bash
python main.py analyse "out/terminators"
```

This will save the plots to the 'out/plots' directory.

# File Descriptions
main.py: The main entry point for the pipeline, handling command-line arguments and orchestrating the different modules.

get_genomes.py: Contains functions for interacting with the NCBI Datasets API to download genomes.

extract.py: Handles the logic for parsing GFF and FASTA files to extract terminator sequences.

analyse.py: Implements the k-mer counting, ranking, and reporting for signal analysis.

plots.py: A helper module for generating the signal distribution plots using Matplotlib.

compare_fasta.py: A utility script to compare two FASTA files and report differences, useful for validation.

# References
Loke, J. C., Stahlberg, E. A., Strenski, D. G., Haas, B. J., Wood, P. C., & Li, Q. Q. (2005). Compilation of mRNA polyadenylation signals in Arabidopsis revealed a new signal element and potential secondary structures. Plant Physiology, 138(3), 1457–1468. https://doi.org/10.1104/pp.105.060541

Beaudoing, E., Freier, S., Wyatt, J. R., Claverie, J. M., & Gautheret, D. (2000). Patterns of variant polyadenylation signal usage in human genes. Genome Research, 10(7), 1001–1010. https://doi.org/10.1101/gr.10.7.1001

de Felippes, F. F., & Waterhouse, P. M. (2023). Plant terminators: the unsung heroes of gene expression. Journal of Experimental Botany, 74(7), 2239–2250. https://doi.org/10.1093/jxb/erac467
