import os


rule convert_sheet:
    """
    converts extended_sample_sheet_template.xlsx to standard sample_sheet.csv
    """
    input:
        inp=config["sample_sheet"],
    output:
        out=os.path.join(config["out_dir"], "sample_sheet.csv"),
    log:
        os.path.join(config["out_dir"], "logs/convert_sheet.log"),
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
        sample_sheet=os.path.join(config["out_dir"], "sample_sheet.csv"),
    output:
        out_dir=directory(os.path.join(config["out_dir"], "demuxed")),
        sta_dir=directory(os.path.join(config["out_dir"], "demuxed/Stats")),
    threads: config.get("threads", 8)
    log:
        os.path.join(config["out_dir"], "logs/demux.log"),
    params:
        barcode_mismatches=config['demux']['barcode_mismatches'],
    conda:
        "../envs/bcl2fastq.yaml"
    shell:
        """
        set -euo pipefail

        bcl2fastq \
            --runfolder-dir {input.run_dir} \
            --output-dir {output.out_dir} \
            --sample-sheet {input.sample_sheet} \
            --mask-short-adapter-reads 0 \
            --minimum-trimmed-read-length 0 \
            --barcode-mismatches {params.barcode_mismatches} \
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
        r1=lambda w: get_fastqs_for_sample(w)["r1"],
        r2=lambda w: get_fastqs_for_sample(w)["r2"],
        r3=lambda w: get_fastqs_for_sample(w)["r3"],
    output:
        bc=os.path.join(config["out_dir"], "merged/{sample}_bc_merged.fastq.gz"),
        cdna=os.path.join(config["out_dir"], "merged/{sample}_cdna_merged.fastq.gz"),
    params:
        sequencer=config["sequencer"],
    threads: config.get("threads", 4)
    log:
        os.path.join(config["out_dir"], "logs/merge_fastq_{sample}.log"),
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        set -euo pipefail
        
        R1={input.r1}
        R2={input.r2}
        R3={input.r3}
        
        echo "R1 files: $R1"
        echo "R2 files: $R2"
        echo "R3 files: $R3"

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
