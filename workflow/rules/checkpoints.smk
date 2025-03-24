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
