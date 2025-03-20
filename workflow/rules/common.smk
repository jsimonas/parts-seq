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


rule convert_sheet:
    """
    converts extended_sample_sheet_template.xlsx to standard sample_sheet.csv
    """
    input:
        inp=config["sample_sheet"],
    output:
        out="results/sample_sheet.csv",
    log:
        "logs/convert_sheet.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/convert_to_samplesheet.py"
