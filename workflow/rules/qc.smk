RUN_MODE = config.get("run_mode", "bcl")
COVERAGE_ENABLED = (
    str(config.get("coverage", {}).get("enabled", False)).lower() == "true"
)


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
        os.path.join(config["out_dir"], "logs/fastqc_{sample}_{type}_{suffix}.log"),
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
        outbase=lambda wc, output: output.report_dir,
    log:
        os.path.join(config["out_dir"], "logs/mirtrace_{sample}.log"),
    conda:
        "../envs/mirtrace.yaml"
    threads: config.get("threads", 4)
    shell:
        """
        mirtrace qc \
            --species {params.species} \
            --output-dir {params.outbase} \
            {input.trimmed_fastq}
        """


rule multiqc:
    input:
        demux=(
            os.path.join(config["out_dir"], "demuxed/Stats")
            if RUN_MODE == "bcl"
            else []
        ),
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
        mirtrace=(
            expand(
                os.path.join(config["out_dir"], "mirtrace/{sample}"),
                sample=get_sample_ids,
            )
            if str(config.get("mirtrace", {}).get("enabled", True)).lower() == "true"
            else []
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
        mirtop_stats=expand(
            os.path.join(config["out_dir"], "mirtop/{sample}_mirtop_stats.log"),
            sample=get_sample_ids,
        ),
        coverage_mqc=(
            os.path.join(config["out_dir"], "qc/coverage/multiqc")
            if COVERAGE_ENABLED
            else []
        ),
        config_file="config/multiqc_config.yaml",
    output:
        html=os.path.join(config["out_dir"], "qc/multiqc_report.html"),
    params:
        extra="--verbose",
    log:
        os.path.join(config["out_dir"], "logs/multiqc.log"),
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        set -euo pipefail

        multiqc \
            --force \
            {input.demux} \
            {input.fastqc} \
            {input.cutadapt} \
            {input.mirtrace} \
            {input.stats} \
            {input.star_logs} \
            {input.star} \
            {input.mirtop_stats} \
            {input.coverage_mqc} \
            -c {input.config_file} \
            --outdir $(dirname {output.html}) \
            {params.extra} &> {log}
        """


rule feature_coverage:
    """
    per-base coverage over each gene from the genome GTF, oriented 5'->3' on
    the feature's own strand. Extends the gene body by `flank` bp on each side.
    """
    input:
        bam=os.path.join(
            config["out_dir"], "mapped/{sample}_Aligned.sortedByCoord.out.bam"
        ),
        gtf=config.get("coverage", {}).get("genome_gtf", ""),
    output:
        cov=os.path.join(config["out_dir"], "qc/coverage/{sample}_cov.csv.gz"),
        features=os.path.join(config["out_dir"], "qc/coverage/{sample}_features.csv"),
    params:
        flank=config.get("coverage", {}).get("flank", 0),
    log:
        os.path.join(config["out_dir"], "logs/feature_coverage_{sample}.log"),
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/feature_coverage.py"


rule metagene_coverage:
    """
    bin per-feature coverage into `nbins` along the gene, z-score per feature
    (or per biotype), and emit one MultiQC custom-content line graph per
    selected biotype (samples as series).
    """
    input:
        cov=expand(
            os.path.join(config["out_dir"], "qc/coverage/{sample}_cov.csv.gz"),
            sample=get_sample_ids,
        ),
        features=expand(
            os.path.join(config["out_dir"], "qc/coverage/{sample}_features.csv"),
            sample=get_sample_ids,
        ),
    output:
        mqc_dir=directory(os.path.join(config["out_dir"], "qc/coverage/multiqc")),
        binned=os.path.join(config["out_dir"], "qc/coverage/metagene_binned.csv.gz"),
    params:
        nbins=config.get("coverage", {}).get("nbins", 30),
        flank=config.get("coverage", {}).get("flank", 0),
        flank_bins=config.get("coverage", {}).get("flank_bins", 0),
        biotypes=config.get("coverage", {}).get("biotypes", []) or [],
        exclude_name_patterns=(
            config.get("coverage", {}).get("exclude_name_patterns", []) or []
        ),
        zscore_by=config.get("coverage", {}).get("zscore_by", "feature"),
        samples=lambda wc, input: [
            os.path.basename(p).replace("_cov.csv.gz", "") for p in input.cov
        ],
    log:
        os.path.join(config["out_dir"], "logs/metagene_coverage.log"),
    conda:
        "../envs/scipy.yaml"
    script:
        "../scripts/metagene_coverage.py"
