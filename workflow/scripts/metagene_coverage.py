#!/usr/bin/env python3

import os
import re
import numpy as np
import pandas as pd

# pylint: disable=undefined-variable


def _bin_segment(rows, meta, seg_len, nseg_bins, out_offset):
    """bin a segment into `nseg_bins` equal-in-bases bins.

    rows        DataFrame with columns feature, cov, bin_local (already
                clipped to 0..nseg_bins-1)
    meta        per-feature metadata (feature, name, biotype, length)
    seg_len     scalar (constant across features, e.g. flank) or the string
                "length" to take per-feature body length from `meta`
    nseg_bins   bins for this segment
    out_offset  added to bin_local to place it on the shared axis
    """
    summed = (
        rows.groupby(["feature", "bin_local"])["cov"]
        .sum()
        .rename("cov_sum")
        .reset_index()
    )
    grid = meta.merge(pd.DataFrame({"bin_local": np.arange(nseg_bins)}), how="cross")
    g = grid.merge(summed, on=["feature", "bin_local"], how="left")
    g["cov_sum"] = g["cov_sum"].fillna(0.0)

    if seg_len == "length":
        L = g["length"].to_numpy()
    else:
        L = np.full(len(g), int(seg_len))

    b = g["bin_local"].to_numpy()
    # exact bin width in bases: p with (p*N)//L == b lies in
    # [ceil(b*L/N), ceil((b+1)*L/N))
    w = (-(-(b + 1) * L // nseg_bins)) - (-(-b * L // nseg_bins))
    g["width"] = w
    g["cov_mean"] = np.where(w > 0, g["cov_sum"] / np.maximum(w, 1), np.nan)
    g["bin"] = g["bin_local"] + out_offset
    return g[["feature", "name", "biotype", "length", "bin", "cov_mean", "width"]]


def bin_coverage(cov_path, feat_path, sample, nbins, flank, flank_bins):
    """per-base coverage -> mean coverage per bin along the flank+body+flank
    axis. pos in the input is already 5'->3' on the feature's own strand.
    """
    cov = pd.read_csv(cov_path)
    feat = pd.read_csv(feat_path)
    d = cov.merge(feat[["feature", "name", "biotype", "length"]], on="feature")
    meta = feat[["feature", "name", "biotype", "length"]].drop_duplicates("feature")

    use_flank = flank > 0 and flank_bins > 0
    segments = []

    body = d[(d["pos"] >= 1) & (d["pos"] <= d["length"])].copy()
    Lb = body["length"].to_numpy()
    body["bin_local"] = np.clip(
        ((body["pos"].to_numpy() - 1) * nbins) // Lb, 0, nbins - 1
    )
    segments.append(
        _bin_segment(body, meta, seg_len="length", nseg_bins=nbins, out_offset=0)
    )

    if use_flank:
        up = d[d["pos"] <= 0].copy()
        up["bin_local"] = np.clip(
            ((up["pos"].to_numpy() + flank - 1) * flank_bins) // flank,
            0,
            flank_bins - 1,
        )
        segments.append(
            _bin_segment(
                up, meta, seg_len=flank, nseg_bins=flank_bins, out_offset=-flank_bins
            )
        )
        dn = d[d["pos"] > d["length"]].copy()
        dn["bin_local"] = np.clip(
            ((dn["pos"].to_numpy() - dn["length"].to_numpy() - 1) * flank_bins) // flank,
            0,
            flank_bins - 1,
        )
        segments.append(
            _bin_segment(
                dn, meta, seg_len=flank, nseg_bins=flank_bins, out_offset=nbins
            )
        )

    g = pd.concat(segments, ignore_index=True).rename(columns={"cov_mean": "cov"})

    # a body shorter than nbins cannot fill all body bins; drop those features
    body = g[(g["bin"] >= 0) & (g["bin"] < nbins)]
    short_features = body.loc[body["width"] == 0, "feature"].unique()
    dropped = len(short_features)
    g = g[~g["feature"].isin(short_features)].drop(columns=["width"])
    g.insert(0, "sample", sample)
    return g, dropped


def add_zscore(d, by):
    """z-score binned coverage within (sample, feature) or (sample, biotype)."""
    keys = ["sample", "feature"] if by == "feature" else ["sample", "biotype"]
    g = d.groupby(keys)["cov"]
    out = d.copy()
    out["z"] = (out["cov"] - g.transform("mean")) / g.transform("std").replace(
        0, np.nan
    )
    return out[out["z"].notna()]


def slugify(s):
    return re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_") or "unknown"


def write_mqc_linegraph(
    path, biotype, wide, nbins, flank_bins, zscore_by, n_features
):
    """MultiQC custom-content linegraph TSV: rows=samples, cols=bins."""
    plot_id = f"metagene_coverage_{slugify(biotype)}"
    title = f"{biotype}"
    desc = (
        f"Mean z-scored per-bin coverage across {n_features} {biotype} features "
        f"(5' -> 3'), z-score computed within each {zscore_by}. "
        f"x=0 marks the 5' end of the gene body; x={nbins} marks the 3' end."
    )
    xmin = -flank_bins if flank_bins > 0 else 0
    xmax = nbins + flank_bins - 1 if flank_bins > 0 else nbins - 1
    header = [
        f"# id: {plot_id}",
        f"# section_name: '{title}'",
        "# parent_id: 'metagene_coverage'",
        "# parent_name: 'Metagene coverage'",
        "# parent_description: >",
        "#     Mean z-scored per-bin coverage along gene bodies for the",
        "#     selected biotypes, one line per sample. Dashed vertical lines",
        "#     mark the 5' and 3' body boundaries when flanks are included.",
        f"# description: \"{desc}\"",
        "# plot_type: 'linegraph'",
        "# pconfig:",
        f"#     id: '{plot_id}'",
        f"#     title: 'Metagene coverage: {biotype}'",
        "#     xlab: \"bin (5' -> 3')\"",
        "#     ylab: 'mean z-score'",
        f"#     xmin: {xmin}",
        f"#     xmax: {xmax}",
    ]
    if flank_bins > 0:
        header += [
            "#     x_lines:",
            "#         - value: 0",
            "#           color: '#888888'",
            "#           dash: 'Dash'",
            f"#         - value: {nbins}",
            "#           color: '#888888'",
            "#           dash: 'Dash'",
        ]
    header.append("")
    with open(path, "w") as fh:
        fh.write("\n".join(header))
        wide.to_csv(fh, sep="\t", index=True, index_label="Sample")


def main(
    cov_paths,
    feat_paths,
    samples,
    nbins,
    flank,
    flank_bins,
    biotypes,
    excludes,
    zscore_by,
    binned_out,
    mqc_dir,
    log_path,
):
    os.makedirs(mqc_dir, exist_ok=True)
    with open(log_path, "w") as log:
        log.write(
            f"nbins={nbins} flank={flank} flank_bins={flank_bins} "
            f"zscore_by={zscore_by}\n"
        )
        parts = []
        for cov_p, feat_p, sample in zip(cov_paths, feat_paths, samples):
            b, dropped = bin_coverage(cov_p, feat_p, sample, nbins, flank, flank_bins)
            log.write(
                f"{sample}: {b['feature'].nunique():,} features "
                f"({dropped:,} shorter than {nbins} bp dropped)\n"
            )
            parts.append(b)
        binned = pd.concat(parts, ignore_index=True)

        d = binned
        if biotypes:
            d = d[d["biotype"].isin(biotypes)]
        for pat in excludes:
            d = d[~d["name"].astype(str).str.contains(pat, na=False)]

        if zscore_by not in ("feature", "biotype"):
            raise ValueError(
                f"coverage.zscore_by must be 'feature' or 'biotype', got {zscore_by!r}"
            )
        d = add_zscore(d, by=zscore_by)
        d.to_csv(binned_out, index=False, compression="gzip")

        summary = (
            d.groupby(["sample", "biotype", "bin"])
            .agg(mean_z=("z", "mean"), n_features=("feature", "nunique"))
            .reset_index()
        )
        for biotype, sub in summary.groupby("biotype"):
            wide = sub.pivot(index="sample", columns="bin", values="mean_z")
            axis = (
                np.arange(-flank_bins, nbins + flank_bins)
                if flank_bins > 0
                else np.arange(nbins)
            )
            wide = wide.reindex(columns=axis)
            out = os.path.join(mqc_dir, f"metagene_coverage_{slugify(biotype)}_mqc.tsv")
            write_mqc_linegraph(
                out,
                biotype=biotype,
                wide=wide,
                nbins=nbins,
                flank_bins=flank_bins,
                zscore_by=zscore_by,
                n_features=int(sub["n_features"].max()),
            )
            log.write(
                f"wrote {out} ({wide.shape[0]} samples x {wide.shape[1]} bins)\n"
            )


if __name__ == "__main__":
    main(
        cov_paths=list(snakemake.input.cov),
        feat_paths=list(snakemake.input.features),
        samples=list(snakemake.params.samples),
        nbins=int(snakemake.params.nbins),
        flank=int(snakemake.params.flank),
        flank_bins=int(snakemake.params.flank_bins),
        biotypes=list(snakemake.params.biotypes),
        excludes=list(snakemake.params.exclude_name_patterns),
        zscore_by=str(snakemake.params.zscore_by),
        binned_out=snakemake.output.binned,
        mqc_dir=snakemake.output.mqc_dir,
        log_path=snakemake.log[0],
    )
