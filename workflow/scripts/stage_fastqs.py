#!/usr/bin/env python3
import glob
import os
import shutil

# pylint: disable=undefined-variable

fastq_dir = snakemake.input.fastq_dir
merged_dir = snakemake.output.merged_dir
marker = snakemake.output.marker
log_file = snakemake.log[0]

with open(log_file, "w") as log_handle:
    try:
        if not os.path.isdir(fastq_dir):
            raise FileNotFoundError(
                "fastq_dir does not exist or is not a directory: " + str(fastq_dir)
            )

        entries = os.listdir(fastq_dir)
        if not entries:
            raise ValueError("fastq_dir is empty: " + str(fastq_dir))

        bc_fastqs = sorted(glob.glob(os.path.join(fastq_dir, "*_bc.fastq.gz")))
        if not bc_fastqs:
            raise ValueError("No *_bc.fastq.gz files found in fastq_dir: " + str(fastq_dir))

        os.makedirs(merged_dir, exist_ok=True)

        for bc_fastq in bc_fastqs:
            bc_name = os.path.basename(bc_fastq)
            sample = bc_name[: -len("_bc.fastq.gz")]
            cdna_fastq = os.path.join(fastq_dir, sample + "_cdna.fastq.gz")

            if not os.path.exists(cdna_fastq):
                raise FileNotFoundError(
                    "Missing mate for sample " + sample + ": expected " + cdna_fastq
                )

            destinations = [
                (
                    os.path.join(merged_dir, sample + "_bc_merged.fastq.gz"),
                    os.path.abspath(bc_fastq),
                ),
                (
                    os.path.join(merged_dir, sample + "_cdna_merged.fastq.gz"),
                    os.path.abspath(cdna_fastq),
                ),
            ]

            for destination, source in destinations:
                if os.path.lexists(destination):
                    os.remove(destination)
                try:
                    os.symlink(source, destination)
                except OSError:
                    shutil.copy2(source, destination)

            print("staged sample: " + sample, file=log_handle)

        with open(marker, "w"):
            pass
    except Exception as exc:
        print(str(exc), file=log_handle)
        raise
