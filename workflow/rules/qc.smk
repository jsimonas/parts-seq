rule fastqc:
    input:
        expand(
            os.path.join(config["out_dir"], "merged/{sample}_bc_001.fastq.gz"),
            sample=get_sample_ids,
        ),
        expand(
            os.path.join(config["out_dir"], "merged/{sample}_cdna_001.fastq.gz"),
            sample=get_sample_ids,
        ),
        expand(
            os.path.join(config["out_dir"], "trimmed/{sample}_bc_trimmed.fastq.gz"),
            sample=get_sample_ids,
        ),
        expand(
            os.path.join(config["out_dir"], "trimmed/{sample}_cdna_trimmed.fastq.gz"),
            sample=get_sample_ids,
        ),
    output:
        html=expand(
            os.path.join(config["out_dir"], "qc/fastqc/{sample}_{read_type}.html"),
            sample=get_sample_ids,
            read_type=["bc_merged", "cdna_merged", "bc_trimmed", "cdna_trimmed"],
        ),
        zip=expand(
            os.path.join(
                config["out_dir"], "qc/fastqc/{sample}_{read_type}_fastqc.zip"
            ),
            sample=get_sample_ids,
            read_type=["bc_merged", "cdna_merged", "bc_trimmed", "cdna_trimmed"],
        ),
    log:
        "logs/fastqc.log",
    threads: 4
    resources:
        mem_mb=8000,
    wrapper:
        "v3.10.0/bio/fastqc"


rule multiqc:
    input:
        fastqc_zips=expand(
            os.path.join(
                config["out_dir"], "qc/fastqc/{sample}_{read_type}_fastqc.zip"
            ),
            sample=get_sample_ids,
            read_type=["bc_merged", "cdna_merged", "bc_trimmed", "cdna_trimmed"],
        ),
        samtools_stats = expand(
            os.path.join(config["out_dir"], "mapped/stats/{sample}.stats"),
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
