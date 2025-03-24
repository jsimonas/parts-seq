import os
import re
import glob
import pandas as pd
from snakemake.utils import validate


# def get_sample_ids(xlsx_path):
#    """
#    obtains sample ids from extended sample sheet
#    """
#    df = pd.read_excel(xlsx_path, engine="openpyxl")
#    df.columns = [col.strip().lower() for col in df.columns]
#    samples = df["sample_id"].unique().tolist()
#
#    return samples


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

    return [R1, R2, R3]
