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


rule mirtrace:
    input:
        trimmed_fastq=expand(
            os.path.join(config["out_dir"], "trimmed/{sample}_cdna_trimmed.fastq.gz"),
            sample=get_sample_ids,
        ),
    output:
        report_dir=directory(os.path.join(config["out_dir"], "mirtrace/{sample}")),
    params:
        species=config["mirtrace"]["species"],
        outbase=os.path.join(config["out_dir"], "mirtrace/{sample}"),
    log:
        "logs/mirtrace_{sample}.log",
    conda:
        "../envs/mirtrace.yaml"
    threads: config.get("threads", 4)
    shell:
        """
        mirtrace --quiet \
                 --species {params.species} \
                 --output {params.outbase} \
                 {input.trimmed_fastq}
        """


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
        cutadapt=expand(
            os.path.join(config["out_dir"], "trimmed/{sample}_cutadapt_report.txt"),
            sample=get_sample_ids,
        ),
        mirtrace=expand(
            os.path.join(config["out_dir"],"/mirtrace/{sample}"),
        ),
        stats=expand(
            os.path.join(config["out_dir"], "mapped/stats/{sample}.stats"),
            sample=get_sample_ids,
        ),
        star_logs=expand(
            os.path.join(config["out_dir"], "mapped/{sample}_Log.final.out"),
            sample=get_sample_ids,
        ),
        star=expand(
            os.path.join(config["out_dir"], "mapped/{sample}_Solo.out"),
            sample=get_sample_ids,
        ),
        config_file="config/multiqc_config.yaml",
    output:
        html=os.path.join(config["out_dir"], "qc/multiqc_report.html"),
    params:
        extra="--verbose",
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
            {input.cutadapt} \
            {input.mirtrace} \
            {input.stats} \
            {input.star_logs} \
            {input.star} \
            -c {input.config_file} \
            --outdir $(dirname {output.html}) \
            {params.extra} &> {log}
        """
