# `parts-seq`


[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/jsimonas/parts-seq/workflows/Tests/badge.svg)](https://github.com/jsimonas/parts-seq/actions?query=branch%3Amain+workflow%3ATests)

A snakemake workflow for single-cell sequencing analysis

### Basic usage

Install snakemake via mamba
```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake
```

Clone the repository
```
git clone https://github.com/jsimonas/parts-seq.git 
```

Perform a dry run of the workflow using Snakemake
```
cd parts-seq
snakemake --cores 2 --dry-run
```

Execute the workflow using Snakemake
```
snakemake \
  --cores 8 \
  --snakefile parts-seq/workflow/Snakefile \
  --directory parts-seq \
  --use-conda \
  --config \
    run_dir=/path/to/sequencing_run_directory \
    out_dir=/path/to/your_output_directory \
    sample_sheet=/path/to/extended_sample_sheet.xlsx \
    sequencer='nextseq' \
    star_index=/path/to/genome_index/ \
    star_features="GeneFull"
```

### Starting from pre-demultiplexed FASTQs (`run_mode: fastq`)

If you already have per-sample merged barcode (`bc`) and cDNA (`cdna`) FASTQ
files, you can skip `bcl2fastq` and the internal merge step entirely.

Place your files in a directory using this naming convention:

```
my_fastqs/
├── sampleA_bc.fastq.gz
├── sampleA_cdna.fastq.gz
├── sampleB_bc.fastq.gz
└── sampleB_cdna.fastq.gz
```

Then run:

```
snakemake \
  --cores 8 \
  --snakefile parts-seq/workflow/Snakefile \
  --directory parts-seq \
  --use-conda \
  --config \
    run_mode=fastq \
    fastq_dir=/path/to/my_fastqs \
    out_dir=/path/to/your_output_directory \
    star_index=/path/to/genome_index/ \
    star_features="GeneFull"
```

In this mode, `run_dir`, `sample_sheet`, and `sequencer` are not required.

Before running the workflow, prepare the [STAR](https://github.com/alexdobin/STAR) index using the genomes and their corresponding annotations of interest.

### Requirements

* `snakemake>≥6.3.0`
* `linux-x64`
