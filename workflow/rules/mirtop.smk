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


rule fasta_to_chrom_gtf:
    """
    convert fasta to gtf
    """
    input:
        hairpin_fa=config["hairpin_fa"],
    output:
        hairpin_gtf=os.path.join(
            config["out_dir"], "reference", "hairpin", "hairpin.gtf"
        ),
    log:
        os.path.join(config["out_dir"], "logs/hairpin_fasta_to_gtf.log"),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/fasta_to_chrom_gtf.py"


rule star_index_hairpin:
    """
    create star index for mirna hairpin sequences.
    """
    input:
        hairpin_fa=config["hairpin_fa"],
        hairpin_gtf=os.path.join(
            config["out_dir"], "reference", "hairpin", "hairpin.gtf"
        ),
    output:
        hairpin_idx=directory(
            os.path.join(config["out_dir"], "reference", "hairpin", "index")
        ),
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
            --sjdbGTFfile {input.hairpin_gtf} \
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
        hairpin_idx=os.path.join(config["out_dir"], "reference", "hairpin", "index"),
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
        stats_json=os.path.join(config["out_dir"], "mirtop/{sample}_mirtop_stats.log"),
        stats_text=os.path.join(config["out_dir"], "mirtop/{sample}_mirtop_stats.txt"),
    log:
        os.path.join(config["out_dir"], "logs/mirtop_{sample}.log"),
    conda:
        "../envs/mirtop.yaml"
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
        
        mv $(dirname {output.gff})/mirtop_stats.txt {output.stats_text}
        mv $(dirname {output.gff})/mirtop_stats.log {output.stats_json}
        
        """


rule starsolo_align_hairpin:
    """
    star align trimmed cDNA (R2) and barcode (R1) reads for mirtop quatification.
    """
    input:
        cdna_read=os.path.join(
            config["out_dir"], "trimmed/{sample}_cdna_trimmed.fastq.gz"
        ),
        bc_read=os.path.join(config["out_dir"], "trimmed/{sample}_bc_trimmed.fastq.gz"),
        bc_1="assets/barcodes/bc1_list.txt",
        bc_2="assets/barcodes/bc2_list.txt",
        hairpin_idx=os.path.join(config["out_dir"], "reference", "hairpin", "index"),
    output:
        hairpin_bam=os.path.join(
            config["out_dir"], "mirtop/{sample}_CB_Aligned.sortedByCoord.out.bam"
        ),
        solo_dir=directory(os.path.join(config["out_dir"], "mirtop/{sample}_CB_Solo.out")),
    params:
        out_prefix=lambda wc, output: output.hairpin_bam.replace(
            "_CB_Aligned.sortedByCoord.out.bam", ""
        ),
        features=config["star"]["features"],
        limit_ram=lambda wc: int(0.8 * parse_size(config["memory"])),
    threads: config.get("threads", 4)
    log:
        os.path.join(config["out_dir"], "logs/star_align_hairpin_cb_{sample}.log"),
    conda:
        "../envs/star.yaml"
    shell:
        """
        set -euo pipefail
        
        ulimit -n 65535
        
        STAR \
            --genomeDir {input.hairpin_idx} \
            --readFilesIn {input.cdna_read} {input.bc_read} \
            --soloCBwhitelist {input.bc_1} {input.bc_2} \
            --runThreadN {threads} \
            --outFileNamePrefix {params.out_prefix}_CB_ \
            --readFilesCommand zcat \
            --runDirPerm All_RWX \
            --outReadsUnmapped Fastx \
            --outSAMtype BAM SortedByCoordinate \
            --limitBAMsortRAM {params.limit_ram} \
            --outSAMattributes NH HI nM AS MD CB CR \
            --outSAMunmapped Within \
            --alignIntronMax 1 \
            --outFilterMultimapScoreRange 0 \
            --outFilterMultimapNmax 10 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --outFilterMismatchNoverLmax 0.05 \
            --outFilterMatchNmin 15 \
            --soloType CB_UMI_Complex \
            --soloFeatures {params.features} \
            --soloCBposition 0_0_0_7 0_8_0_15 \
            --soloUMIposition 0_16_0_23 \
            --soloBarcodeReadLength 1 \
            --soloUMIdedup NoDedup \
            --soloCellReadStats Standard \
            --soloCBmatchWLtype EditDist_2 &> {log}
    
        """


checkpoint split_bam_by_barcode:
    """
    split to hairpin aligned BAM file into individual BAM files per valid barcode.
    """
    input:
        bam=os.path.join(
            config["out_dir"], "mirtop/{sample}_CB_Aligned.sortedByCoord.out.bam"
        ),
    output:
        split_dir=directory(os.path.join(config["out_dir"], "mirtop/split/{sample}")),
        barcode_list=os.path.join(
            config["out_dir"], "mirtop/split/{sample}/valid_barcodes.txt"
        ),
    params:
        n_cells=config.get("mirtop", {}).get("n_cells", 1000),
        n_reads=config.get("mirtop", {}).get("n_reads", 100),
        cell_stats=lambda wc, input: input.bam.replace(
            "_Aligned.sortedByCoord.out.bam",
            f"_CB_Solo.out/{config['star']['features']}/CellReads.stats",
        ),
    log:
        os.path.join(config["out_dir"], "logs/split_bam_{sample}.log"),
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        set -euo pipefail
        exec > "{log}" 2>&1

        tail -n +3 "{params.cell_stats}" \
            | sort -k2,2nr \
            | head -n {params.n_cells} \
            | awk -v threshold={params.n_reads} '$2 > threshold {{print $1}}' \
            > "{output.barcode_list}"

        mkdir -p "{output.split_dir}"
        
        TMP_BAM=$(mktemp "{output.split_dir}/temp.XXXXXX.bam")

        N_CB=$(wc -l < "{output.barcode_list}")
        
        ulimit -n $(( N_CB > 1000 ? N_CB + 50 : 1024 ))
        
        samtools view -h "{input.bam}" \
            | awk 'BEGIN{{OFS="\\t"}} /^@/ {{print; next}} {{$1 = $1 "_x1"; print}}' \
            | samtools view -Sb - > "$TMP_BAM"

        samtools view -u \
            -U /dev/null \
            -D "CB:{output.barcode_list}" \
            "$TMP_BAM" \
        | samtools split -d CB -M ${{N_CB}} -f "{output.split_dir}/%!.bam" -

        rm "$TMP_BAM"

        """


rule mirtop_counts_per_barcode:
    """
    run mirtop gff and counts for a single barcode.
    """
    input:
        bam=os.path.join(config["out_dir"], "mirtop/split/{sample}/{cb}.bam"),
        hairpin=config.get("hairpin_fa", "assets/hairpin.fa"),
        gtf=config.get("gtf", "assets/mirna.gtf"),
    output:
        tsv=os.path.join(config["out_dir"], "mirtop/split/{sample}/{cb}_mirtop.tsv"),
        gff=os.path.join(config["out_dir"], "mirtop/split/{sample}/{cb}.gff"),
    log:
        os.path.join(config["out_dir"], "logs/mirtop_counts/{sample}_{cb}.log"),
    conda:
        "../envs/mirtop.yaml"
    shell:
        """
        set -e
        
        mirtop gff --hairpin {input.hairpin} \
                   --gtf {input.gtf} \
                   --sps {wildcards.cb} \
                   --out {output.gff} \
                   {input.bam} > {log} 2>&1

        TMP_DIR=$(mktemp -d -t mirtop_XXXXXX)
        
        mirtop counts --gff {output.gff} --out "$TMP_DIR" >> {log} 2>&1
        
        mv "$TMP_DIR/mirtop.tsv" "{output.tsv}"
        rm -rf "$TMP_DIR"
        
        """
