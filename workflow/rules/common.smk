import os
import glob
import pandas as pd
from snakemake.utils import validate


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
