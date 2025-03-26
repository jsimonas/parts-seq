import os


rule starsolo:
    input:
        cdna_read=os.path.join(
            config["out_dir"], "trimmed/{sample}_cdna_trimmed.fastq.gz"
        ),
        bc_read=os.path.join(
            config["out_dir"], "trimmed/{sample}_bc_trimmed.fastq.gz"
        ),
        bc_1="../assets/barcodes/bc1_list.txt",
        bc_2="../assets/barcodes/bc2_list.txt",
        index=config["star_index"],
    output:
        bam=os.path.join(
            config["out_dir"], "mapped/{sample}_Aligned.sortedByCoord.out.bam"
        ),
        solo_dir=directory(os.path.join(config["out_dir"], "mapped/{sample}_Solo.out")),
    params:
        out_prefix=lambda wildcards, output: output.bam.replace("Aligned.sortedByCoord.out.bam", ""),
    threads: config.get("threads", 4)
    log:
        "logs/starsolo_{sample}.log",
    conda:
        "../envs/starsolo.yaml"
    shell:
        """
        set -euo pipefail

        STAR \
            --genomeDir {input.index} \
            --readFilesIn {input.cdna_read} {input.bc_read} \
            --soloCBwhitelist {input.bc_1} {input.bc_2} \
            --runThreadN {threads} \
            --outFileNamePrefix {params.out_prefix} \
            --readFilesCommand zcat \
            --runDirPerm All_RWX \
            --outReadsUnmapped Fastx \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI nM AS CR UR CB UB sS sQ sM GX GN \
            --outSAMunmapped Within \
            --alignIntronMax 1 \
            --outFilterMultimapScoreRange 0 \
            --outFilterMultimapNmax 10 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --outFilterMismatchNoverLmax 0.05 \
            --outFilterMatchNmin 15 \
            --soloType CB_UMI_Complex \
            --soloMultiMappers PropUnique \
            --soloFeatures GeneFull \
            --soloUMIdedup Exact \
            --soloCBposition 0_0_0_7 0_8_0_15 \
            --soloUMIposition 0_16_0_23 \
            --soloBarcodeReadLength 1 \
            --soloCBmatchWLtype EditDist_2 &> {log}
        """


rule samtools_stats:
    input:
        bam=os.path.join(
            config["out_dir"], "mapped/{sample}_Aligned.sortedByCoord.out.bam"
        ),
    output:
        stats=os.path.join(config["out_dir"], "mapped/stats/{sample}.stats"),
    params:
        extra="",
    log:
        "logs/samtools_stats_{sample}.log",
    wrapper:
        "v3.13.3/bio/samtools/stats"
