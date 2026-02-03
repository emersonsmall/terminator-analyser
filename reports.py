from constants import (
    CE_EXPECTED_END,
    CE_EXPECTED_START,
    NUE_EXPECTED_END,
    NUE_EXPECTED_START,
)


def save_overview_report(out_fpath: str, metadata: dict) -> None:
    separator_width = 100

    lines = []
    lines.append("=" * separator_width + "\n")
    lines.append("Terminator Analysis Overview\n")
    lines.append("=" * separator_width + "\n\n")

    lines.append("ANALYSIS METADATA\n")
    lines.append("-" * separator_width + "\n")
    if metadata.get("taxon"):
        lines.append(f"Taxon:                  {metadata['taxon']}\n")
    lines.append(
        f"Generated:              {metadata['generated_at'].strftime('%Y-%m-%d %H:%M:%S')}\n"
    )
    lines.append(f"K-mer Size:             {metadata['kmer_size']}\n")
    lines.append(f"Step Size:              {metadata['step_size']}\n")
    lines.append(f"Min 3'UTR Length:       {metadata['min_3utr_length']} nt\n")
    lines.append(f"Num Downstream nt:      {metadata['num_downstream_nt']} nt\n")
    lines.append(f"Total Terminators:      {metadata['total_terminators']:,}\n")
    lines.append(f"Included Terminators:   {metadata['included_terminators']:,}\n")
    lines.append(f"Skipped Terminators:    {metadata['skipped_terminators']:,}\n")
    skip_pct = (
        (metadata["skipped_terminators"] / metadata["total_terminators"] * 100)
        if metadata["total_terminators"] > 0
        else 0
    )
    lines.append(f"Skip Rate:              {skip_pct:.2f}%\n")
    lines.append(f"Number of Accessions:   {metadata['num_accessions']}\n")
    lines.append("\n")

    accession_info = metadata.get("accession_info") or {}
    extraction_stats = metadata.get("extraction_stats") or {}

    if accession_info:
        lines.append("TERMINATOR COUNTS BY ACCESSION\n")
        lines.append("-" * separator_width + "\n")

        header = {
            "accession": "Accession",
            "organism": "Organism name",
            "total": "Total",
            "included": "Included",
            "skipped": "Skipped",
            "pct": "% of Included",
        }
        total_all = metadata["total_terminators"]
        total_included = metadata["included_terminators"]
        total_skipped = metadata["skipped_terminators"]

        max_widths = {
            "accession": len(header["accession"]),
            "organism": len(header["organism"]),
            "total": len(header["total"]),
            "included": len(header["included"]),
            "skipped": len(header["skipped"]),
            "pct": len(header["pct"]),
        }

        for accession, info in accession_info.items():
            max_widths["accession"] = max(max_widths["accession"], len(accession))
            max_widths["organism"] = max(
                max_widths["organism"], len(info["organism_name"])
            )

        max_widths["total"] = max(max_widths["total"], len(f"{total_all:,}"))
        max_widths["included"] = max(max_widths["included"], len(f"{total_included:,}"))
        max_widths["skipped"] = max(max_widths["skipped"], len(f"{total_skipped:,}"))

        fmt = (
            f"{{accession:<{max_widths['accession']}}} | "
            f"{{organism:<{max_widths['organism']}}} | "
            f"{{total:>{max_widths['total']}}} | "
            f"{{included:>{max_widths['included']}}} | "
            f"{{skipped:>{max_widths['skipped']}}} | "
            f"{{pct:>{max_widths['pct']}}}"
        )
        sep = (
            f"{'-' * max_widths['accession']} | "
            f"{'-' * max_widths['organism']} | "
            f"{'-' * max_widths['total']} | "
            f"{'-' * max_widths['included']} | "
            f"{'-' * max_widths['skipped']} | "
            f"{'-' * max_widths['pct']}"
        )

        lines.append(fmt.format(**header) + "\n")
        lines.append(sep + "\n")

        sorted_accessions = sorted(
            accession_info.items(), key=lambda x: x[1]["total"], reverse=True
        )
        for accession, info in sorted_accessions:
            pct = (info["included"] / total_included * 100) if total_included > 0 else 0
            lines.append(
                fmt.format(
                    accession=accession,
                    organism=info["organism_name"],
                    total=f"{info['total']:,}",
                    included=f"{info['included']:,}",
                    skipped=f"{info['skipped']:,}",
                    pct=f"{pct:.2f}%",
                )
                + "\n"
            )

        lines.append(sep + "\n")
        lines.append(
            fmt.format(
                accession="TOTAL",
                organism=" " * max_widths["organism"],
                total=f"{total_all:,}",
                included=f"{total_included:,}",
                skipped=f"{total_skipped:,}",
                pct="100.00%",
            )
            + "\n"
        )
        lines.append("\n")

    # Extraction skip reasons per accession (if available)
    if extraction_stats:
        lines.append("TRANSCRIPT COUNTS BY ACCESSION\n")
        lines.append("-" * separator_width + "\n")

        # Collect all reason keys across accessions to build columns
        all_reasons = set()
        for info in extraction_stats.values():
            all_reasons.update(info.get("skip_reasons", {}).keys())
        all_reasons = sorted(all_reasons)

        if not all_reasons:
            lines.append("No transcripts skipped during extraction.\n\n")
        else:
            header = {
                "accession": "Accession",
                "total": "Total",
                "extracted": "Extracted",
                "skipped": "Skipped",
                **{r: r for r in all_reasons},
            }

            # Compute column widths
            maxw = {
                "accession": len(header["accession"]),
                "total": len(header["total"]),
                "extracted": len(header["extracted"]),
                "skipped": len(header["skipped"]),
            }
            for r in all_reasons:
                maxw[r] = len(r)

            # Accumulate totals
            total_extracted = 0
            total_skipped = 0
            total_transcripts = 0
            total_reasons = {r: 0 for r in all_reasons}

            for accession, info in extraction_stats.items():
                num_extracted = info.get("num_extracted", 0)
                skip_reasons = info.get("skip_reasons", {})
                acc_skipped = sum(skip_reasons.values())
                acc_transcripts = num_extracted + acc_skipped

                total_extracted += num_extracted
                total_skipped += acc_skipped
                total_transcripts += acc_transcripts
                for r in all_reasons:
                    total_reasons[r] += skip_reasons.get(r, 0)

                maxw["accession"] = max(maxw["accession"], len(accession))
                maxw["total"] = max(maxw["total"], len(f"{acc_transcripts:,}"))
                maxw["extracted"] = max(maxw["extracted"], len(f"{num_extracted:,}"))
                maxw["skipped"] = max(maxw["skipped"], len(f"{acc_skipped:,}"))
                for r in all_reasons:
                    maxw[r] = max(maxw[r], len(f"{skip_reasons.get(r,0):,}"))

            # Also consider TOTAL row widths
            maxw["accession"] = max(maxw["accession"], len("TOTAL"))
            maxw["total"] = max(maxw["total"], len(f"{total_transcripts:,}"))
            maxw["extracted"] = max(maxw["extracted"], len(f"{total_extracted:,}"))
            maxw["skipped"] = max(maxw["skipped"], len(f"{total_skipped:,}"))
            for r in all_reasons:
                maxw[r] = max(maxw[r], len(f"{total_reasons[r]:,}"))

            # Build format strings
            fmt_parts = [
                f"{{accession:<{maxw['accession']}}}",
                f"{{total:>{maxw['total']}}}",
                f"{{extracted:>{maxw['extracted']}}}",
                f"{{skipped:>{maxw['skipped']}}}",
            ]
            sep_parts = [
                "-" * maxw["accession"],
                "-" * maxw["total"],
                "-" * maxw["extracted"],
                "-" * maxw["skipped"],
            ]
            for r in all_reasons:
                fmt_parts.append(f"{{{r}:>{maxw[r]}}}")
                sep_parts.append("-" * maxw[r])

            fmt = " | ".join(fmt_parts)
            sep = " | ".join(sep_parts)

            # Write header
            lines.append(fmt.format(**header) + "\n")
            lines.append(sep + "\n")

            # Rows per accession (sorted by accession)
            for accession, info in sorted(extraction_stats.items(), key=lambda x: x[0]):
                num_extracted = info.get("num_extracted", 0)
                skip_reasons = info.get("skip_reasons", {})
                acc_skipped = sum(skip_reasons.values())
                row = {
                    "accession": accession,
                    "total": f"{num_extracted + acc_skipped:,}",
                    "extracted": f"{num_extracted:,}",
                    "skipped": f"{acc_skipped:,}",
                }
                for r in all_reasons:
                    row[r] = f"{skip_reasons.get(r, 0):,}"
                lines.append(fmt.format(**row) + "\n")

            # Separator and TOTAL row
            lines.append(sep + "\n")
            total_row = {
                "accession": "TOTAL",
                "total": f"{total_transcripts:,}",
                "extracted": f"{total_extracted:,}",
                "skipped": f"{total_skipped:,}",
            }
            for r in all_reasons:
                total_row[r] = f"{total_reasons[r]:,}"
            lines.append(fmt.format(**total_row) + "\n\n")

    lines.append("=" * separator_width + "\n")

    with open(out_fpath, "w") as f:
        f.writelines(lines)

    print(f"Overview report saved to '{out_fpath}'")


def save_region_report(
    out_fpath: str,
    region_name: str,
    ranked_kmers: list,
    region_window: tuple[int, int],
    kmer_size: int,
) -> None:

    separator_width = 100

    lines = []
    lines.append("=" * separator_width + "\n")
    lines.append(f"{region_name} Region Report\n")
    lines.append("=" * separator_width + "\n\n")

    lines.append(f"Analysis Window:        {region_window[0]} to {region_window[1]}\n")

    # Add expected region info
    if region_name == "NUE":
        lines.append(
            f"Expected Region:        {NUE_EXPECTED_START} to {NUE_EXPECTED_END}\n"
        )
    elif region_name == "CE":
        lines.append(
            f"Expected Region:        {CE_EXPECTED_START} to {CE_EXPECTED_END}\n"
        )

    lines.append("\n")

    lines.append(f"TOP {len(ranked_kmers)} K-MERS\n")
    lines.append("-" * separator_width + "\n")

    if not ranked_kmers:
        lines.append("No k-mers found.\n")
    else:
        header = {
            "kmer": "K-mer",
            "delta": "Delta",
            "peak": "Peak Count",
            "median": "Median Count",
            "pos": "Peak Pos",
            "occurrences": "Occurrence Count",
            "pct": "% Occurrence",
        }

        maxw = {
            "kmer": max(kmer_size, len(header["kmer"])),
            "delta": len(header["delta"]),
            "peak": len(header["peak"]),
            "median": len(header["median"]),
            "pos": len(header["pos"]),
            "occurrences": len(header["occurrences"]),
            "pct": len(header["pct"]),
        }

        for item in ranked_kmers:
            maxw["delta"] = max(maxw["delta"], len(f"{item['delta']:.1f}"))
            maxw["peak"] = max(maxw["peak"], len(f"{item['peak_count']:,}"))
            maxw["median"] = max(maxw["median"], len(f"{item['median_count']:.1f}"))
            maxw["pos"] = max(maxw["pos"], len(str(item["peak_pos"])))
            maxw["occurrences"] = max(
                maxw["occurrences"], len(f"{item['occurrence_count']:,}")
            )
            maxw["pct"] = max(maxw["pct"], len(f"{item['pct_occurrence']:.1f}%"))

        fmt = (
            f"{{kmer:<{maxw['kmer']}}} | "
            f"{{delta:>{maxw['delta']}}} | "
            f"{{peak:>{maxw['peak']}}} | "
            f"{{median:>{maxw['median']}}} | "
            f"{{pos:>{maxw['pos']}}} | "
            f"{{occurrences:>{maxw['occurrences']}}} | "
            f"{{pct:>{maxw['pct']}}}"
        )
        sep = (
            f"{'-' * maxw['kmer']} | "
            f"{'-' * maxw['delta']} | "
            f"{'-' * maxw['peak']} | "
            f"{'-' * maxw['median']} | "
            f"{'-' * maxw['pos']} | "
            f"{'-' * maxw['occurrences']} | "
            f"{'-' * maxw['pct']}"
        )

        lines.append(fmt.format(**header) + "\n")
        lines.append(sep + "\n")
        for item in ranked_kmers:
            lines.append(
                fmt.format(
                    kmer=item["kmer"],
                    delta=f"{item['delta']:.1f}",
                    peak=f"{item['peak_count']:,}",
                    median=f"{item['median_count']:.1f}",
                    pos=str(item["peak_pos"]),
                    occurrences=f"{item['occurrence_count']:,}",
                    pct=f"{item['pct_occurrence']:.1f}%",
                )
                + "\n"
            )

    lines.append("\n")
    lines.append("=" * separator_width + "\n")

    with open(out_fpath, "w") as f:
        f.writelines(lines)
    print(f"{region_name} report saved to '{out_fpath}'")
