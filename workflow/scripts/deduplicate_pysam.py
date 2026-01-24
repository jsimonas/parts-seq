#!/usr/bin/env python3
import pysam
from collections import defaultdict

# pylint: disable=undefined-variable

def main(input_bam, output_bam, log_file):
    """
    deduplicate BAM by cell barcode (CB) and UMI (UB), and add _x1 to read names.
    """

    with open(log_file, 'w') as log:
        log.write("Starting deduplication...\n")
        
        # track unique reads per cell
        seen = defaultdict(set)
        total_reads = 0
        dedup_reads = 0
        
        with pysam.AlignmentFile(input_bam, "rb") as inbam:
            with pysam.AlignmentFile(output_bam, "wb", template=inbam) as outbam:
                for read in inbam:
                    total_reads += 1
                    
                    # get CB and UB tags
                    cb = read.get_tag("CB") if read.has_tag("CB") else None
                    ub = read.get_tag("UB") if read.has_tag("UB") else None
                    
                    if cb is None or ub is None:
                        continue
                    
                    # create unique key: cell_barcode + umi + position + strand
                    key = (cb, ub, read.reference_name, read.reference_start, read.is_reverse)
                    
                    # skip if we've seen this combination
                    if key in seen[cb]:
                        continue
                    
                    seen[cb].add(key)
                    dedup_reads += 1
                    
                    # add _x1 to read name
                    read.query_name = read.query_name + "_x1"
                    
                    # create CX tag (CB without underscore)
                    cx = cb.replace("_", "")
                    read.set_tag("CX", cx, value_type="Z")
                    
                    outbam.write(read)
        
        log.write(f"Total reads: {total_reads}\n")
        log.write(f"Deduplicated reads: {dedup_reads}\n")
        log.write(f"Duplicates removed: {total_reads - dedup_reads}\n")
        if total_reads > 0:
            log.write(f"Deduplication rate: {100 * (1 - dedup_reads/total_reads):.2f}%\n")
        
        pysam.index(output_bam)


if __name__ == "__main__":
    main(
        snakemake.input.bam,
        snakemake.output.dedup_bam,
        snakemake.log[0]
    )