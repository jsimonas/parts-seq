#!/usr/bin/env python3
import os
import re
import pandas as pd
from io import StringIO

# pylint: disable=undefined-variable

sample = snakemake.params.sample
solo_dir = snakemake.input.solo_dir
feature = snakemake.params.feature
output_files = snakemake.output

def modify_summary(content):
    df = pd.read_csv(StringIO(content), header=None)
    df = df.transpose()
    new_header = ['Sample Name'] + df.iloc[0, 1:].tolist()
    new_values = [sample] + df.iloc[1, 1:].tolist()
    df.columns = new_header
    df = pd.DataFrame([new_values], columns=new_header)
    return df.to_csv(index=False)


def modify_umi(content):
    lines = content.strip().splitlines()
    numbered_lines = [f"{i+1}\t{line}" for i, line in enumerate(lines)]
    return "\n".join(numbered_lines) + "\n"


def modify_stats(content):
    lines = content.strip().splitlines()
    processed = []
    for line in lines:
        line = re.sub(r'^\s+', '', line)
        line = re.sub(r'\s+', '\t', line)
        processed.append(line)
    return "\n".join(processed) + "\n"


feature_dir = os.path.join(solo_dir, feature)
file_map = {
    os.path.join(feature_dir, "Summary.csv"): (modify_summary, output_files.summary),
    os.path.join(feature_dir, "UMIperCellSorted.txt"): (modify_umi, output_files.umi),
    os.path.join(solo_dir, "Barcodes.stats"): (modify_stats, output_files.barcodes),
    os.path.join(feature_dir, "Features.stats"): (modify_stats, output_files.features),
}

for inpath, (modifier, outpath) in file_map.items():
    with open(inpath, "r") as infile:
        content = infile.read()
    with open(outpath, "w") as outfile:
        outfile.write(modifier(content))
