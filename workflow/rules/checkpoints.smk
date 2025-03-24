checkpoint parse_demux:
    input:
        demux_dir = "results/demultiplexed",
    output:
        sample_ids = "results/sample_ids.txt",
    log:
        "logs/parse_demux.log",
    script:
        "../scripts/parse_demux.py"