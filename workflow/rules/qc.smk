rule fastqc:
    input:
        os.path.join(config["out_dir"], "{stage}/{sample}_{read_type}.fastq.gz"),
    output:
        html=os.path.join(config["out_dir"], "qc/fastqc/{sample}_{read_type}.html"),
        zip=os.path.join(config["out_dir"], "qc/fastqc/{sample}_{read_type}_fastqc.zip"),
    log:
        "logs/fastqc_{sample}_{read_type}.log",
    threads: config.get("threads", 4)
    resources:
        mem_mb=8000,
    wrapper:
        "v3.10.0/bio/fastqc"


rule multiqc:
    input:
        fastqc=expand(
            os.path.join(config["out_dir"], "qc/fastqc/{sample}_{read_type}_fastqc.zip"),
            sample=get_sample_ids,
            read_type=["bc_merged", "cdna_merged", "bc_trimmed", "cdna_trimmed"]
        ),
        stats=expand(
            os.path.join(config["out_dir"], "mapped/stats/{sample}.stats"),
            sample=get_sample_ids
        ),
        config_file="config/multiqc_config.yaml",
    output:
        os.path.join(config["out_dir"], "qc/multiqc.html")
    params:
        extra="--verbose"
    log:
        "logs/multiqc.log"
    wrapper:
        "v3.10.0/bio/multiqc"
