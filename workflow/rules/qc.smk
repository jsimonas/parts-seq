rule fastqc:
    input:
        lambda wc: os.path.join(
            config["out_dir"],
            "merged" if wc.suffix == "merged" else "trimmed",
            "{sample}_{type}_{suffix}.fastq.gz",
        ),
    output:
        html=os.path.join(config["out_dir"], "qc/fastqc/{sample}_{type}_{suffix}.html"),
        zip=os.path.join(
            config["out_dir"], "qc/fastqc/{sample}_{type}_{suffix}_fastqc.zip"
        ),
    log:
        "logs/fastqc_{sample}_{type}_{suffix}.log",
    threads: config.get("threads", 4)
    resources:
        mem_mb=8000,
    wrapper:
        "v3.10.0/bio/fastqc"


rule multiqc:
    input:
        demux=os.path.join(config["out_dir"], "demultiplexed/Stats"),
        fastqc=expand(
            os.path.join(
                config["out_dir"], "qc/fastqc/{sample}_{type}_{suffix}_fastqc.zip"
            ),
            sample=get_sample_ids,
            type=["bc", "cdna"],
            suffix=["merged", "trimmed"],
        ),
        stats=expand(
            os.path.join(config["out_dir"], "mapped/stats/{sample}.stats"),
            sample=get_sample_ids,
        ),
        star=expand(
            os.path.join(config["out_dir"], "mapped/{sample}_Solo.out"),
            sample=get_sample_ids,
        ),
        starsolo=expand(
            os.path.join(config["out_dir"], "mapped/{sample}_Solo.out/GeneFull/"),
            sample=get_sample_ids,
        ),
        config_file="config/multiqc_config.yaml",
    output:
        html=os.path.join(config["out_dir"], "qc/multiqc_report.html"),
    params:
        extra="--verbose"
    log:
        "logs/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        set -euo pipefail

        multiqc \
            {input.demux} \
            {input.fastqc} \
            {input.stats} \
            {input.star} \
            {input.starsolo} \
            -c {input.config_file} \
            --outdir $(dirname {output.html}) \
            {params.extra} &> {log}
        """
