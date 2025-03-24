checkpoint parse_demux:
    input:
        demux_dir="results/demultiplexed",
    output:
        sample_ids="results/sample_ids.txt",
    log:
        "logs/parse_demux.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/parse_demux.py"
