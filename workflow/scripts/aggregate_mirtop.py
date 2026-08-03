#!/usr/bin/env python3
import pandas as pd
import scipy.sparse as sp
import scipy.io as io
import gzip
import os

# pylint: disable=undefined-variable

def read_barcode_counts(tsv_file, barcode):
    # completely empty file -> "No columns to parse from file"
    try:
        df = pd.read_csv(tsv_file, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.Series(dtype="int64", name=barcode)

    # header only, or no miRNA / no count column -> nothing mapped
    if df.empty or "miRNA" not in df.columns or df.shape[1] < 2:
        return pd.Series(dtype="int64", name=barcode)

    count_col = df.columns[-1]
    counts = df.groupby("miRNA")[count_col].sum().astype("int64")
    counts.name = barcode
    return counts


def _guess_database(gff_path):
    database = None
    with open(gff_path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if "miRBase" in line:
                database = line[line.find("miRBase"):].strip().replace(" ", "")
            elif "MirGeneDB" in line:
                database = line[line.find("MirGeneDB"):].strip().replace(" ", "")
            elif "microRNAs" in line:
                database = line.strip().split()[1]
    return database


def _attr(col9, key):
    for field in col9.split(";"):
        field = field.strip()
        if field.startswith(key + "="):
            return field.split("=", 1)[1]
    return None


def read_reference_mirnas(gff_path):
    database = _guess_database(gff_path) or ""
    prefer_id = "MirGeneDB" in database
    names, seen = [], set()
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "miRNA":
                continue
            if prefer_id:
                name = _attr(cols[8], "ID") or _attr(cols[8], "Name")
            else:
                name = _attr(cols[8], "Name") or _attr(cols[8], "ID")
            if name and name not in seen:
                seen.add(name)
                names.append(name)
    return names


try:
    series_list = []

    for tsv_file in snakemake.input.tsv_files:
        if "invalid_barcodes" in tsv_file:
            continue
        barcode = os.path.basename(tsv_file).replace("_mirtop.tsv", "")
        series_list.append(read_barcode_counts(tsv_file, barcode))

    if not series_list:
        raise ValueError("No valid TSV files found")

    n_nonzero = sum(not s.empty for s in series_list)
    ordered_barcodes = [s.name for s in series_list]

    # detected counts: union of observed miRNAs across barcodes; empty barcodes
    # contribute no index and become all-zero columns.
    detected = (
        pd.concat(series_list, axis=1, join="outer")
        .reindex(columns=ordered_barcodes)
        .fillna(0)
        .astype("int64")
    )

    # reference-complete feature set: every miRNA in the reference annotation,
    # in reference order, so miRNAs that were never detected appear as all-zero
    # rows. Any detected miRNA not found in the reference (should not happen) is
    # appended rather than dropped.
    reference = read_reference_mirnas(snakemake.input.mirna_gtf)
    ref_set = set(reference)
    extra = [m for m in detected.index if m not in ref_set]
    full_index = reference + extra if reference else list(detected.index)

    matrix = detected.reindex(index=full_index).fillna(0).astype("int64")

    n_features, n_barcodes = matrix.shape
    n_detected_features = int((matrix.values.sum(axis=1) > 0).sum()) if n_features else 0

    if n_features:
        sparse_matrix = sp.csr_matrix(matrix.values)
    else:
        sparse_matrix = sp.csr_matrix((0, n_barcodes), dtype="int64")

    # write matrix.mtx.gz
    with gzip.open(snakemake.output.matrix, "wb") as f:
        io.mmwrite(f, sparse_matrix, comment="", field="integer")

    # write features.tsv.gz
    features_df = pd.DataFrame({
        "gene_id": matrix.index,
        "gene_name": matrix.index,
        "feature_type": "miRNA",
    })
    features_df.to_csv(
        snakemake.output.features, sep="\t", header=False, index=False,
        compression="gzip",
    )

    # write barcodes.tsv.gz
    barcodes_df = pd.DataFrame(matrix.columns)
    barcodes_df.to_csv(
        snakemake.output.barcodes, sep="\t", header=False, index=False,
        compression="gzip",
    )

    denom = n_features * n_barcodes
    sparsity = (1 - sparse_matrix.nnz / denom) if denom else 0.0

    with open(snakemake.log[0], "w") as f:
        f.write(
            f"Successfully aggregated {n_barcodes} barcodes "
            f"({n_nonzero} with >=1 miRNA, "
            f"{n_barcodes - n_nonzero} all-zero)\n"
        )
        f.write(
            f"Features: {n_features} miRNAs "
            f"({len(reference)} from reference, {len(extra)} detected-only; "
            f"{n_detected_features} detected, "
            f"{n_features - n_detected_features} all-zero)\n"
        )
        f.write(f"Matrix dimensions: {n_features} miRNAs x {n_barcodes} barcodes\n")
        f.write(f"Total counts: {int(matrix.values.sum()) if n_features else 0}\n")
        f.write(f"Sparsity: {sparsity:.2%}\n")

except Exception as e:
    with open(snakemake.log[0], "w") as f:
        f.write(f"ERROR: {str(e)}\n")
        import traceback
        f.write(traceback.format_exc())
    raise
