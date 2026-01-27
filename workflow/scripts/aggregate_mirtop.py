#!/usr/bin/env python3
import pandas as pd
import scipy.sparse as sp
import scipy.io as io
import gzip
import os

# pylint: disable=undefined-variable

try:
    dfs = []
    
    for tsv_file in snakemake.input:
        if "invalid_barcodes" not in tsv_file:
            barcode = os.path.basename(tsv_file).replace("_mirtop.tsv", "")
            
            # read TSV
            df = pd.read_csv(tsv_file, sep="\t")
            
            # extract miRNA name and count (last column)
            counts = df.groupby('miRNA')[df.columns[-1]].sum()
            counts.name = barcode
            dfs.append(counts.to_frame())
    
    if not dfs:
        raise ValueError("No valid TSV files found")
    
    # combine into matrix
    matrix = pd.concat(dfs, axis=1, join="outer").fillna(0).astype(int)
    
    # convert to sparse matrix
    sparse_matrix = sp.csr_matrix(matrix.values)
    
    # write matrix.mtx.gz
    with gzip.open(snakemake.output.matrix, 'wb') as f:
        io.mmwrite(f, sparse_matrix, comment='', field='integer')
    
    # write features.tsv.gz
    features_df = pd.DataFrame({
        'gene_id': matrix.index,
        'gene_name': matrix.index,
        'feature_type': 'miRNA'
    })
    features_df.to_csv(snakemake.output.features, sep='\t', header=False, index=False, compression='gzip')
    
    # write barcodes.tsv.gz
    barcodes_df = pd.DataFrame(matrix.columns)
    barcodes_df.to_csv(snakemake.output.barcodes, sep='\t', header=False, index=False, compression='gzip')
    
    with open(snakemake.log[0], "w") as f:
        f.write(f"Successfully aggregated {len(dfs)} barcodes\n")
        f.write(f"Matrix dimensions: {matrix.shape[0]} miRNAs x {matrix.shape[1]} barcodes\n")
        f.write(f"Total counts: {matrix.values.sum()}\n")
        f.write(f"Sparsity: {1 - sparse_matrix.nnz / (matrix.shape[0] * matrix.shape[1]):.2%}\n")

except Exception as e:
    with open(snakemake.log[0], "w") as f:
        f.write(f"ERROR: {str(e)}\n")
        import traceback
        f.write(traceback.format_exc())
    raise