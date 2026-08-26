#!/usr/bin/env python3

import gzip
import os
import re

import pandas as pd
import pysam

# pylint: disable=undefined-variable

_ATTR_QUOTED = re.compile(r'(\w+)\s+"([^"]*)"')
_ATTR_EQ = re.compile(r"(\w+)=([^;]+)")


def parse_attrs(field):
    d = dict(_ATTR_QUOTED.findall(field))
    if not d:
        d = dict(_ATTR_EQ.findall(field))
    return d


def read_gene_features(gtf_path):
    """Yield (chrom, start1, end1, strand, feature_id, name, biotype) for every
    "gene" record. Coordinates are 1-based inclusive as in GTF/GFF."""
    open_fn = gzip.open if gtf_path.endswith(".gz") else open
    with open_fn(gtf_path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom, start, end, strand = parts[0], int(parts[3]), int(parts[4]), parts[6]
            attrs = parse_attrs(parts[8])
            fid = attrs.get("gene_id") or attrs.get("ID")
            if not fid:
                continue
            name = attrs.get("gene_name") or attrs.get("Name") or fid
            biotype = (
                attrs.get("gene_biotype")
                or attrs.get("gene_type")
                or attrs.get("biotype")
                or "unknown"
            )
            yield chrom, start, end, strand, fid, name, biotype


def main(bam_path, gtf_path, flank, cov_out, features_out, log_path):
    with open(log_path, "w") as log:
        log.write(f"reading GTF: {gtf_path}\n")
        log.write(f"flank: {flank} bp\n")
        feats = list(read_gene_features(gtf_path))
        log.write(f"  {len(feats):,} gene features\n")

        if not os.path.exists(bam_path + ".bai"):
            log.write(f"indexing {bam_path}\n")
            pysam.index(bam_path)

        bam = pysam.AlignmentFile(bam_path, "rb")
        chrom_lens = dict(zip(bam.references, bam.lengths))

        cov_rows = []
        feat_rows = []
        skipped_missing = 0
        skipped_edge = 0

        for chrom, start, end, strand, fid, name, biotype in feats:
            if chrom not in chrom_lens:
                skipped_missing += 1
                continue
            chrom_len = chrom_lens[chrom]
            # need the full flank on both sides so every feature contributes
            # identically to flank bins
            if start - 1 < flank or end + flank > chrom_len:
                skipped_edge += 1
                continue
            length = end - start + 1
            gs0 = start - 1 - flank
            ge0 = end + flank
            counts = bam.count_coverage(
                chrom, gs0, ge0, quality_threshold=0, read_callback="all"
            )
            cov = [a + c + g + t for a, c, g, t in zip(*counts)]
            if strand == "-":
                cov = cov[::-1]
            # pos of the first element in cov after strand handling: 1 - flank
            pos0 = 1 - flank

            feat_rows.append((fid, name, biotype, length))
            for i, c in enumerate(cov):
                if c:
                    cov_rows.append((fid, pos0 + i, c))

        bam.close()
        log.write(
            f"  {skipped_missing:,} features on contigs missing from BAM\n"
            f"  {skipped_edge:,} features too close to a contig edge for full flank\n"
        )

        pd.DataFrame(
            feat_rows, columns=["feature", "name", "biotype", "length"]
        ).to_csv(features_out, index=False)
        pd.DataFrame(cov_rows, columns=["feature", "pos", "cov"]).to_csv(
            cov_out, index=False, compression="gzip"
        )
        log.write(
            f"wrote {len(feat_rows):,} features and "
            f"{len(cov_rows):,} covered positions\n"
        )


if __name__ == "__main__":
    main(
        snakemake.input.bam,
        snakemake.input.gtf,
        int(snakemake.params.flank),
        snakemake.output.cov,
        snakemake.output.features,
        snakemake.log[0],
    )
