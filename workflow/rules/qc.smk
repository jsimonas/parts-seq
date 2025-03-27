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
        solo_summary=expand(
            os.path.join(
                config["out_dir"], "mapped/{sample}_Solo.out/GeneFull/Summary.csv"
            ),
            sample=get_sample_ids,
        ),
        solo_umi=expand(
            os.path.join(
                config["out_dir"],
                "mapped/{sample}_Solo.out/GeneFull/UMIperCellSorted.txt",
            ),
            sample=get_sample_ids,
        ),
        solo_feature_stats=expand(
            os.path.join(
                config["out_dir"], "mapped/{sample}_Solo.out/GeneFull/Features.stats"
            ),
            sample=get_sample_ids,
        ),
        solo_barcode_stats=expand(
            os.path.join(config["out_dir"], "mapped/{sample}_Solo.out/Barcodes.stats"),
            sample=get_sample_ids,
        ),
        config_file="config/multiqc_config.yaml",
    output:
        os.path.join(config["out_dir"], "qc/multiqc.html"),
    params:
        extra="--verbose",
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.10.0/bio/multiqc"
