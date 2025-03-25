import os
import re
import glob
import pandas as pd
from snakemake.utils import validate


def get_sample_ids(wildcards):
    ckpt = checkpoints.parse_demux.get()
    sample_file = ckpt.output.sample_ids
    with open(sample_file) as f:
        samples = f.read().strip().splitlines()

    return samples


def get_fastqs_for_sample(wildcards):
    """
    detects the demultiplexed FASTQ files for a given sample in results/demultiplexed and
    returns a dictionary
    """
    pattern = os.path.join(
        "results", "demultiplexed", "**", f"{wildcards.sample}_S*_R*_001.fastq.gz"
    )
    matches = sorted(glob.glob(pattern, recursive=True))

    if len(matches) < 3:
        raise FileNotFoundError(
            f"Did not find at least 3 FASTQ files for sample={wildcards.sample} matching {pattern}"
        )

    R1 = next((f for f in matches if "_R1_" in f), None)
    R2 = next((f for f in matches if "_R2_" in f), None)
    R3 = next((f for f in matches if "_R3_" in f), None)
    if not (R1 and R2 and R3):
        raise FileNotFoundError(
            f"Could not identify R1, R2, R3 among these matches: {matches}"
        )

    return {"r1": R1, "r2": R2, "r3": R3}

rule dirs:
    output:
        trim_dir=directory(config["out_dir"] + "/trimmed"),
        logs_dir=directory(config["out_dir"] + "/logs")
    shell:
        """
        mkdir -p {output.trim_dir} {output.logs_dir}
        """
