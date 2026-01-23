#!/usr/bin/env python3
import sys

# pylint: disable=undefined-variable

def parse_fasta(path):
    seqs = {}
    cur = None
    with open(path, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
                name = header.split()[0]
                # ensure unique name
                if name in seqs:
                    idx = 2
                    newname = f"{name}_{idx}"
                    while newname in seqs:
                        idx += 1
                        newname = f"{name}_{idx}"
                    name = newname
                seqs[name] = ""
                cur = name
            else:
                if cur is None:
                    raise SystemExit("FASTA malformed: sequence data before header")
                seqs[cur] += line.strip()
    return seqs

def write_gtf(path, seqs):
    with open(path, "w") as out:
        for seqname, seq in seqs.items():
            if not seq:
                continue
            start = 1
            end = len(seq)
            gene_id = seqname
            if gene_id.endswith("_pre"):
                gene_id = gene_id[: -len("_pre")] or seqname
            transcript_id = seqname

            def attrs(d):
                return "; ".join(f'{k} "{v}"' for k, v in d.items()) + ";"

            # gene
            out.write("\t".join([
                seqname, "fasta2gtf", "gene",
                str(start), str(end), ".", ".", ".",
                attrs({"gene_id": gene_id, "gene_name": gene_id, "gene_biotype": "miRNA"})
            ]) + "\n")

            # transcript
            out.write("\t".join([
                seqname, "fasta2gtf", "transcript",
                str(start), str(end), ".", ".", ".",
                attrs({"gene_id": gene_id, "transcript_id": transcript_id, "gene_name": gene_id, "gene_biotype": "miRNA"})
            ]) + "\n")

            # exon
            out.write("\t".join([
                seqname, "fasta2gtf", "exon",
                str(start), str(end), ".", ".", ".",
                attrs({"gene_id": gene_id, "transcript_id": transcript_id, "exon_number": "1", "gene_biotype": "miRNA"})
            ]) + "\n")


fasta = snakemake.input[0]
out = snakemake.output[0]

seqs = parse_fasta(fasta)
write_gtf(out, seqs)