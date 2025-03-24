#!/usr/bin/env python3
import os
import re

# pylint: disable=undefined-variable

demux_dir = snakemake.input.demux_dir
output_file = snakemake.output.sample_ids

sample_name_set = set()
for root, dirs, files in os.walk(demux_dir):
    for f in files:
        m = re.match(r"(.*?)_S.*_R[123]_001\.fastq\.gz", f)
        if m:
            sample_name_set.add(m.group(1))


samples = sorted(sample_name_set)
with open(output_file, "w") as outF:
    for s in samples:
        outF.write(s + "\n")

        
snakemake.params.samples = samples
