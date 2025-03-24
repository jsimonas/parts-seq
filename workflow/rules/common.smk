import os
import glob
import pandas as pd
from snakemake.utils import validate


def get_fastqs_for_sample(wildcards):
    """
    detects the demultiplexed FASTQ files for a given sample in results/demultiplexed.
    """
    sample_dir = f"results/demultiplexed/{wildcards.run}"

    pattern = os.path.join(sample_dir, f"{wildcards.sample}_R*_001.fastq.gz")
    matches = sorted(glob.glob(pattern))

    if len(matches) < 3:
        raise FileNotFoundError(
            f"did not find at least 3 FASTQ files for sample={wildcards.sample} in {sample_dir} matching {pattern}"
        )

    R1 = next((f for f in matches if "_R1_" in f), None)
    R2 = next((f for f in matches if "_R2_" in f), None)
    R3 = next((f for f in matches if "_R3_" in f), None)
    if not (R1 and R2 and R3):
        raise FileNotFoundError(
            f"could not identify R1, R2, R3 among these matches: {matches}"
        )

    return [R1, R2, R3]


def get_ids_from_sample_sheet(csv_path):
    lines = []
    with open(csv_path, "r") as f:
        data_section_found = False
        for line in f:
            if line.strip() == "[Data]":
                data_section_found = True
                continue
            if data_section_found and line.strip():
                lines.append(line.strip().split(","))

    df = pd.DataFrame(lines[1:], columns=lines[0])
    samples = df["Sample_ID"].unique().tolist()
    return samples


def get_samples_from_demux(_wildcards):
    ckpt = checkpoints.demux.get()
    return ckpt.params["samples"]