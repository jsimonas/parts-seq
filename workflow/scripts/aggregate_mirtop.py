#!/usr/bin/env python3
import pandas as pd

# pylint: disable=undefined-variable

dfs = []
for tsv_file in snakemake.input:
    if "invalid_barcodes" not in tsv_file:
        barcode = tsv_file.split("/")[-1].replace("_mirtop.tsv", "")
        df = pd.read_csv(tsv_file, sep="\t", index_col=0)
        df.columns = [barcode]
        dfs.append(df)

matrix = pd.concat(dfs, axis=1, join="outer").fillna(0)
matrix.to_csv(snakemake.output.matrix, sep="\t")

with open(snakemake.log[0], "w") as f:
    f.write(f"aggregated {len(dfs)} barcodes\n")
