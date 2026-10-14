"""Bulk RNA-seq of sorted healthy human haematopoietic populations (GSE74246):
pure profiles and known-composition mixtures for deconvolution."""

from __future__ import annotations

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.figure import Figure

import deconvolution_utils as du

SCRIPT           = "corces_bulk"
COUNTS_FILE      = du.EXTERNAL_DIR / "GSE74246_RNAseq_All_Counts.txt.gz"
LENGTH_FILE      = du.EXTERNAL_DIR / "refGene_hg19.txt.gz"
MAIN_CHROMOSOMES = {f"chr{name}" for name in (*range(1, 23), "X", "Y")}
SORTED_GROUPS    = {
    "HSC":      "HSPC",
    "MPP":      "HSPC",
    "LMPP":     "HSPC",
    "CMP":      "HSPC",
    "GMP":      "HSPC",
    "MEP":      "HSPC",
    "CLP":      "HSPC",
    "Ery":      "Erythroid",
    "Mono":     "Myeloid",
    "Bcell":    "B lineage",
    "CD4Tcell": "T cells",
    "CD8Tcell": "T cells",
    "NKcell":   "NK",
}
MIXTURE_TYPES    = {
    "HSPC":      ("HSC", "MPP", "CMP"),
    "Erythroid": ("Ery",),
    "Myeloid":   ("Mono",),
    "B lineage": ("Bcell",),
    "T cells":   ("CD4Tcell", "CD8Tcell"),
    "NK":        ("NKcell",),
}
N_PER_SET        = 250
CONCENTRATION    = 100.0
N_PROFILE_GENES  = 2000
PER_MILLION      = 1e6


def load_sorted() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Raw counts of the healthy sorted populations (genes by samples, first
    occurrence of duplicated symbols kept) and sample annotation."""
    counts = pd.read_table(COUNTS_FILE, index_col=0)
    counts = counts.loc[~counts.index.duplicated(keep="first")]
    parts = counts.columns.str.split("-", n=1, expand=True)
    meta = pd.DataFrame({"donor": parts.get_level_values(0),
                         "sorted_type": parts.get_level_values(1)},
                        index=counts.columns)
    meta = meta[meta["sorted_type"].isin(SORTED_GROUPS)]
    meta["group"] = meta["sorted_type"].map(SORTED_GROUPS)
    return counts[meta.index], meta


def union_length(intervals: list[tuple[int, int]]) -> int:
    """Total length covered by a set of half-open intervals."""
    total, reach = 0, -1
    for start, end in sorted(intervals):
        if end > reach:
            total += end - max(start, reach)
            reach = end
    return total


def gene_lengths() -> pd.Series:
    """Length in kilobases of the union of all RefSeq exons per gene symbol
    (hg19), longest over the main chromosomes."""
    table = pd.read_table(LENGTH_FILE, header=None, usecols=[2, 9, 10, 12],
                          names=["chrom", "starts", "ends", "symbol"])
    table = table[table["chrom"].isin(MAIN_CHROMOSOMES)]
    exons = table.assign(
        starts=table["starts"].str.rstrip(",").str.split(","),
        ends=table["ends"].str.rstrip(",").str.split(","),
    ).explode(["starts", "ends"]).astype({"starts": int, "ends": int})
    lengths = {key: union_length(list(zip(rows["starts"], rows["ends"],
                                          strict=True)))
               for key, rows in exons.groupby(["symbol", "chrom"])}
    return pd.Series(lengths).groupby(level=0).max() / 1e3


def per_million(counts: pd.DataFrame) -> pd.DataFrame:
    """Values scaled to one million per sample (genes by samples)."""
    return counts / counts.sum() * PER_MILLION


def mixture_weights(meta: pd.DataFrame, fractions: pd.DataFrame,
                    rng: np.random.Generator) -> pd.DataFrame:
    """Weight of every sorted sample in every mixture: one randomly chosen
    sample per group, weighted by the group fraction."""
    weights = pd.DataFrame(0.0, index=meta.index, columns=fractions.index)
    for mixture, row in fractions.iterrows():
        for group, fraction in row[row > 0].items():
            pool = meta.index[meta["sorted_type"].isin(MIXTURE_TYPES[group])]
            weights.loc[rng.choice(pool), mixture] += fraction
    return weights


def label_profiles(adata: ad.AnnData) -> pd.DataFrame:
    """Per-million profile of the summed counts of every cell-type label
    (genes by labels)."""
    labels = adata.obs[du.CELL_TYPE_KEY]
    summed = {label: np.asarray(adata.X[(labels == label).to_numpy()]
                                .sum(axis=0)).ravel()
              for label in du.BROAD_GROUPS if (labels == label).any()}
    return per_million(pd.DataFrame(summed, index=adata.var_names))


def profile_correlation(bulk: pd.DataFrame,
                        reference: pd.DataFrame) -> pd.DataFrame:
    """Pearson correlation of log per-million profiles between bulk samples
    and reference labels over the genes most variable across labels."""
    logged = np.log1p(reference)
    genes = logged.var(axis=1).nlargest(N_PROFILE_GENES).index
    matrix = np.corrcoef(np.log1p(bulk.loc[genes]).T, logged.loc[genes].T)
    return pd.DataFrame(matrix[:bulk.shape[1], bulk.shape[1]:],
                        index=bulk.columns, columns=reference.columns)


def plot_correlation(correlation: pd.DataFrame, variant: str) -> Figure:
    """Heatmap of correlation between sorted bulk samples and the
    single-cell label profiles."""
    figure, axis = plt.subplots(
        figsize=(15, 3 + 0.32 * len(correlation)), layout="constrained",
    )
    du.style_dark(figure, [axis])
    image = axis.imshow(correlation.to_numpy(), cmap="magma", aspect="auto")
    axis.set_xticks(range(correlation.shape[1]), correlation.columns,
                    rotation=60, ha="right")
    axis.set_yticks(range(len(correlation)), correlation.index)
    axis.set_title(f"Sorted bulk samples ({variant}) against single-cell "
                   f"label profiles (Pearson r, log values, top "
                   f"{N_PROFILE_GENES} genes)", fontsize=du.TITLE_FONTSIZE,
                   pad=du.TITLE_PAD)
    bar = figure.colorbar(image, ax=axis, shrink=0.6, pad=0.02)
    bar.set_label("Pearson r", color=du.FOREGROUND,
                  fontsize=du.LABEL_FONTSIZE)
    bar.ax.tick_params(colors=du.FOREGROUND, labelsize=du.TICK_FONTSIZE)
    return figure


def main() -> None:
    """Pure and mixed sorted profiles on the shared genes, as raw counts and
    gene-length corrected, with true group shares and a profile check."""
    rng = np.random.default_rng(du.SEED)
    counts, meta = load_sorted()
    reference = du.load_counts(du.TRAIN_DONORS)
    lengths = gene_lengths()
    shared = counts.index.intersection(reference.var_names)
    measured = shared.intersection(lengths.index)
    profiles = {
        "counts": (per_million(counts), shared),
        "length": (per_million(counts.loc[measured]
                               .div(lengths[measured], axis=0)), measured),
    }
    order = meta.assign(
        rank=meta["sorted_type"].map(list(SORTED_GROUPS).index),
    ).sort_values("rank", kind="stable").index

    groups = list(MIXTURE_TYPES)
    observed = (reference.obs[du.CELL_TYPE_KEY]
                .map(du.BROAD_GROUPS).value_counts(normalize=True)[groups])
    fraction_sets = {
        "random": du.random_fractions(N_PER_SET, len(groups), rng),
        "marrow": rng.dirichlet(CONCENTRATION * observed.to_numpy()
                                / observed.sum(), N_PER_SET),
    }
    weights, truths = [], []
    for kind, values in fraction_sets.items():
        fractions = pd.DataFrame(
            values, columns=groups,
            index=[f"{kind}_{i}" for i in range(N_PER_SET)],
        )
        weights.append(mixture_weights(meta, fractions, rng))
        truths.append(fractions.assign(kind=kind))

    pure = pd.get_dummies(meta["group"], dtype=float).assign(kind="pure")
    truth = pd.concat([pure, *truths]).reindex(
        columns=["kind", *du.GROUP_ORDER]).fillna(0.0)
    truth = truth.join(meta[["sorted_type", "donor"]])

    label_cpm = label_profiles(reference[:, shared])
    correlations = {}
    for variant, (profile, genes) in profiles.items():
        bulk = pd.concat([profile, *(profile @ weight for weight in weights)],
                         axis=1).loc[genes]
        out = du.script_dir(du.INTERIM_DIR, SCRIPT) / f"{SCRIPT}_{variant}.txt"
        bulk.to_csv(out, sep="\t", float_format="%.3f")
        print(f"wrote {out.relative_to(du.PROJECT_ROOT)} "
              f"({bulk.shape[0]:,} genes x {bulk.shape[1]:,} samples)")
        correlations[variant] = profile_correlation(profile.loc[genes],
                                                    label_cpm.loc[genes])
        du.save_figure(plot_correlation(correlations[variant].loc[order],
                                        variant),
                       SCRIPT, f"{SCRIPT}_correlation_{variant}")

    summary = pd.Series({
        "genes in the bulk file": len(counts),
        "genes shared with the single-cell data": len(shared),
        "shared genes with a RefSeq length": len(measured),
        "healthy sorted samples": len(meta),
        "mixtures per set": N_PER_SET,
        "largest deviation of a true row sum from 1":
            float((truth[list(du.GROUP_ORDER)].sum(axis=1) - 1).abs().max()),
    }, name="value", dtype=object)
    correlation = pd.concat(correlations, names=["bulk", "sample"])
    du.save_table(truth, SCRIPT, "corces_bulk_truth")
    du.save_table(correlation.join(meta["sorted_type"], on="sample"), SCRIPT,
                  "corces_bulk_correlation")
    du.save_table(summary.to_frame(), SCRIPT, "corces_bulk_summary")
    du.save_table(lengths.rename("length_kb").to_frame(), SCRIPT,
                  "corces_bulk_gene_lengths")


if __name__ == "__main__":
    main()