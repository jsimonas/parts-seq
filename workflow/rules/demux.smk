rule demux:
    """
    runs bcl2fastq to demultiplex using the converted sample sheet.
    """
    input:
        run_dir = config["run_dir"],
        sample_sheet = "results/sample_sheet.csv",
    output:
        directory("results/demultiplexed"),
    threads: config.get("threads", 8)
    log:
        "logs/demultiplex.log",
    conda:
        "envs/bcl2fastq.yaml"
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
