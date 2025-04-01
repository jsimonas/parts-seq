import os


checkpoint parse_demux:
    input:
        demux_dir=os.path.join(config["out_dir"], "demuxed"),
    output:
        sample_ids=os.path.join(config["out_dir"], "sample_ids.txt"),
    log:
        os.path.join(config["out_dir"], "logs/parse_demux.log"),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/parse_demux.py"
