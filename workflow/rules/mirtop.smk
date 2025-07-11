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
        tmpdir=$(mktemp -d)

        seqcluster collapse \
            -f {input.R2_trimmed} \
            -o $tmpdir \
            &> {log}

        gzip -c "$tmpdir/{wildcards.sample}_cdna_trimmed_trimmed.fastq" \
            > {output.R2_collapsed}
        rm -rf "$tmpdir"
        """


rule star_index_hairpin:
    """
    create star index for mirna hairpin sequences.
    """
    input:
        hairpin_fa=config["hairpin_fa"],
    output:
        hairpin_idx=directory(os.path.join(config["out_dir"], "reference", "hairpin")),
    threads: config.get("threads", 4)
    log:
        os.path.join(config["out_dir"], "logs/star_index_hairpin.log"),
    conda:
        "../envs/star.yaml"
    shell:
        """
        set -euo pipefail
        
        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output.hairpin_idx} \
            --genomeFastaFiles {input.hairpin_fa} \
            &> {log}

        """


rule star_align_hairpin:
    """
    star align trimmed and collapsed cDNA (R2) reads for mirtop stats.
    """
    input:
        R2_collapsed=os.path.join(
            config["out_dir"], "collapsed/{sample}_cdna_collapsed.fastq.gz"
        ),
        hairpin_idx=directory(os.path.join(config["out_dir"], "reference", "hairpin")),
    output:
        hairpin_bam=os.path.join(
            config["out_dir"], "mirtop/{sample}_Aligned.sortedByCoord.out.bam"
        ),
    params:
        out_prefix=lambda wc, output: output.hairpin_bam.replace(
            "Aligned.sortedByCoord.out.bam", ""
        ),
        limit_ram=lambda wc: int(0.8 * parse_size(config["memory"])),
    threads: config.get("threads", 4)
    log:
        os.path.join(config["out_dir"], "logs/star_align_hairpin_{sample}.log"),
    conda:
        "../envs/star.yaml"
    shell:
        """
        set -euo pipefail
        
        ulimit -n 65535
        
        STAR \
            --genomeDir {input.hairpin_idx} \
            --readFilesIn {input.R2_collapsed} \
            --runThreadN {threads} \
            --outFileNamePrefix {params.out_prefix} \
            --readFilesCommand zcat \
            --runDirPerm All_RWX \
            --outReadsUnmapped Fastx \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM {params.limit_ram} \
            --outSAMattributes NH HI nM AS MD \
            --outSAMunmapped Within \
            --alignIntronMax 1 \
            --outFilterMultimapScoreRange 0 \
            --outFilterMultimapNmax 10 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --outFilterMismatchNoverLmax 0.05 \
            --outFilterMatchNmin 15 &> {log}

        """

rule mirtop:
    input:
        hairpin_bam=os.path.join(
            config["out_dir"], "mirtop/{sample}_Aligned.sortedByCoord.out.bam"
        ),
        hairpin_fa=config["hairpin_fa"],
        mirna_gtf=config["mirna_gtf"],
    output:
        gff=os.path.join(
            config["out_dir"], "mirtop/{sample}_Aligned.sortedByCoord.out.gff"
        ),
        stats=os.path.join(
            config["out_dir"], "mirtop/{sample}_Aligned.sortedByCoord.out.stats"
        ),
    conda:
        "envs/mirtop.yaml"
    shell:
        """
        mirtop gff \
            --sps Hsa \
            --hairpin {input.hairpin_fa} \
            --gtf {input.mirna_gtf} \
            --out $(dirname {output.gff}) \
            {input.hairpin_bam}
        
        mirtop stats \
            --out $(dirname {output.gff}) \
            {output.gff}
        
        mv $(dirname {output.gff})/$(basename {output.gff}).stats {output.stats}
        """
