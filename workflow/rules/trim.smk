import os


rule trim_reads:
    """
    trim polyA tails from cDNA (R2) reads and filter barcodes (R1) based on the omitted R2 reads
    """
    input:
        R1=os.path.join(config["out_dir"], "merged/{sample}_bc_001.fastq.gz"),
        R2=os.path.join(config["out_dir"], "merged/{sample}_cdna_001.fastq.gz"),
    output:
        R1_trimmed=os.path.join(
            config["out_dir"],"trimmed/{sample}_bc_trimmed.fastq.gz"
        ),
        R2_trimmed=os.path.join(
            config["out_dir"], "trimmed/{sample}_cdna_trimmed.fastq.gz"
        ),
    threads: config.get("threads", 4)
    params:
        trim_5p=config.get("trim_5p", 0),
    log:
        "logs/trim_fastq_{sample}.log",
    conda:
        "../envs/cutadapt.yaml"
    shell:
        """
        set -euo pipefail
        
        cutadapt -u {params.trim_5p} -a A{{8}} -m 15 -j {threads} \
            -o {output.R2_trimmed} \
            -p {output.R1_trimmed} \
            {input.R2} {input.R1} &> {log}
        """
