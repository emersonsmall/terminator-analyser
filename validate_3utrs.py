import sys
import os

import pyfastx

# COMPARE TO ANNOTATED THALIANA
# Sequences match for the example I checked manually !

# TAIR annotation has gene IDs, my 3utrs have rna transcript IDs

# each gene can have multiple RNA transcripts/RNA IDs for different isoforms/rna versions

# my feature map only has base gene id, does not capture gene version suffix

# remember case sensitivity for fasta sequences

# merge 3'UTRs of isoforms if they are the same/similar enough?

annotated_3utrs = pyfastx.Fasta("./TAIR10_3_utr_20101028.fa")
predicted_3utrs = pyfastx.Fasta("out/3utrs/thaliana_3utrs.fa")

with open(os.path.join(input_dir, "thaliana.gff"), 'r') as f:
    lines = f.readlines()

f_map = build_feature_map(lines)
for pred_3utr in predicted_3utrs:
    rna_id = pred_3utr.name[:-5] # remove "_3utr" suffix
    gene_id = f_map.get(rna_id, {})[0].get("parent")
    gene_id = gene_id.split('-')[1] # remove "gene-" prefix
    print(rna_id, gene_id)


if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)