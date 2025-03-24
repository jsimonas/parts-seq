import os
import glob
import pandas as pd
from snakemake.utils import validate


def get_sample_ids(xlsx_path):
    """
    obtains sample ids from extended sample sheet
    """
    df = pd.read_excel(xlsx_path, engine="openpyxl")
    df.columns = [col.strip().lower() for col in df.columns]
    samples = df["sample_id"].unique().tolist()

    return samples


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


def get_fastq_dict(wildcards):
    """
    returns a dictionary with keys 'r1', 'r2', and 'r3'
    mapping to the corresponding FASTQ files for a given sample.
    """
    fastqs = get_fastqs_for_sample(wildcards)
    return {"r1": fastqs[0], "r2": fastqs[1], "r3": fastqs[2]}
