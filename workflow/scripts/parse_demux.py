#!/usr/bin/env python3
import os
import re

# pylint: disable=undefined-variable

if hasattr(snakemake.input, "merged_dir"):
    scan_dir = snakemake.input.merged_dir
    pattern = re.compile(r"^(.+?)_bc_merged\.fastq\.gz$")
else:
    scan_dir = snakemake.input.demux_dir
    pattern = re.compile(r"^(.+?)_S\d+_R[123]_001\.fastq\.gz$")

output_file = snakemake.output.sample_ids

sample_name_set = set()
for root, dirs, files in os.walk(scan_dir):
    for f in files:
        m = pattern.match(f)
        if m and m.group(1) != "Undetermined":
            sample_name_set.add(m.group(1))

samples = sorted(sample_name_set)
with open(output_file, "w") as outF:
    for s in samples:
        outF.write(s + "\n")
