#!/usr/bin/env python3
import sys
import pysam
import shutil
import gzip

# pylint: disable=undefined-variable

# logging to the snakemake log file
sys.stderr = open(snakemake.log[0], "w")

def reverse_complement(seq):
    complement_map = str.maketrans('ACGTN', 'TGCAN')
    return seq.translate(complement_map)[::-1]


def combine_fastq(r1_file, r2_file, r3_file, bc_file, cdna_file, sequencer):
    # validate sequencer
    if sequencer not in ["miseq", "nextseq"]:
        print(f"ERROR: 'sequencer' must be 'miseq' or 'nextseq'. Found: {sequencer}", file=sys.stderr)
        sys.exit(1)

    with pysam.FastxFile(r1_file) as f1, pysam.FastxFile(r2_file) as f2, gzip.open(bc_file, 'wt') as out:
        for r1, r2 in zip(f1, f2):
            # process R2 based on the sequencer
            r2_seq = r2.sequence
            r2_qual = r2.quality if r2.quality is not None else ""
            
            if sequencer == "nextseq":
                r2_seq = reverse_complement(r2_seq)
                r2_qual = r2_qual[::-1]

            r1_seq = r1.sequence
            r1_qual = r1.quality if r1.quality is not None else ""
            
            # concatenate R2 and R1
            combined_seq = r2_seq + r1_seq
            combined_qual = r2_qual + r1_qual
            
            # build the header using R2 information
            header = "@" + r2.name
            if r2.comment:
                header += " " + r2.comment
            
            # write output
            out.write(f"{header}\n{combined_seq}\n+\n{combined_qual}\n")

    # copy r3_file to cdna_file
    with open(r3_file, 'rb') as f_in, open(cdna_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


if __name__ == "__main__":
    combine_fastq(
        r1_file=snakemake.input.r1,
        r2_file=snakemake.input.r2,
        r3_file=snakemake.input.r3,
        bc_file=snakemake.output.bc,
        cdna_file=snakemake.output.cdna,
        sequencer=snakemake.params.sequencer
    )
