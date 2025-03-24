rule convert_sheet:
    """
    converts extended_sample_sheet_template.xlsx to standard sample_sheet.csv
    """
    input:
        inp=config["sample_sheet"],
    output:
        out="results/sample_sheet.csv",
    log:
        "logs/convert_sheet.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/convert_to_samplesheet.py"


rule demux:
    """
    runs bcl2fastq to demultiplex using the converted sample sheet.
    """
    input:
        run_dir=config["run_dir"],
        sample_sheet="results/sample_sheet.csv",
    output:
        directory("results/demultiplexed"),
    threads: config.get("threads", 8)
    log:
        "logs/demultiplex.log",
    conda:
        "../envs/bcl2fastq.yaml"
    shell:
        """
        set -euo pipefail

        bcl2fastq \
            --runfolder-dir {input.run_dir} \
            --output-dir {output} \
            --sample-sheet {input.sample_sheet} \
            --mask-short-adapter-reads 0 \
            --minimum-trimmed-read-length 0 \
            --use-bases-mask 'y*,I*,y*,y*' \
            --no-lane-splitting \
            --create-fastq-for-index-reads \
            --processing-threads {threads} \
            &> {log}
        """


rule merge_fastq:
    """
    merge R1 and R2 into something and copies R3, depending on sequencer
    """
    input:
        fastqs=get_fastqs_for_sample,
    output:
        bc="results/merged/{sample}_bc_001.fastq.gz",
        cdna="results/merged/{sample}_cdna_001.fastq.gz",
    params:
        sequencer=config["sequencer"],
    threads: config.get("threads", 4)
    log:
        "logs/merge_fastq_{sample}.log",
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        set -euo pipefail
        
        R1={input.fastqs.r1}
        R2={input.fastqs.r2}
        R3={input.fastqs.r3}

        if [ "{params.sequencer}" = "miseq" ]; then
            seqkit concat $R2 $R1 --out-file {output.bc} --line-width 0 --threads {threads}
            cp $R3 {output.cdna}
        elif [ "{params.sequencer}" = "nextseq" ]; then
            seqkit concat <(seqkit seq --reverse --complement --seq-type dna $R2) $R1 \
                          --out-file {output.bc} --line-width 0 --threads {threads}
            cp $R3 {output.cdna}
        else
            echo "ERROR: 'sequencer' must be 'miseq' or 'nextseq'. Found: {params.sequencer}"
            exit 1
        fi
        """

checkpoint parse_demux:
    """
    scans the demultiplexed FASTQ directory and extracts sample IDs.
    """
    input:
        demux_dir = "results/demultiplexed"
    output:
        sample_ids = "results/sample_ids.txt"
    run:
        import os, re
        sample_names = set()
        for root, dirs, files in os.walk(input.demux_dir):
            for f in files:
                # We expect files like: 21_SMALL_SJ_18_S1_R1_001.fastq.gz
                m = re.match(r"(.*?)_S.*_R[123]_001\.fastq\.gz", f)
                if m:
                    sample_names.add(m.group(1))
        samples = sorted(sample_names)
        # Write sample IDs to an output file (optional)
        with open(output.sample_ids, "w") as outF:
            for s in samples:
                outF.write(s + "\n")
        params["samples"] = samples
