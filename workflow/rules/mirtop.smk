import os


rule collapse_reads:
    """
    collapse trimmed cDNA (R2) reads for mirtop stats.
    """
    input:
        R2_trimmed=os.path.join(
            config["out_dir"], "trimmed/{sample}_cdna_trimmed.fastq.gz"
        ),
    output:
        R2_collapsed=os.path.join(
            config["out_dir"], "collapsed/{sample}_cdna_collapsed.fastq.gz"
        ),
    threads: config.get("threads", 4)
    log:
        os.path.join(config["out_dir"], "logs/collapse_fastq_{sample}.log"),
    conda:
        "../envs/seqcluster.yaml"
    shell:
        """
        set -euo pipefail
        
        seqcluster collapse \
            -f {input.R2_trimmed} \
            -o {output.R2_collapsed} \
            > {log}

        """
