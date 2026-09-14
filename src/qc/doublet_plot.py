"""Plot doublet score distribution and checks, per donor. """

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

IN = Path("data/interim/bm_clean.h5ad")
FIG = Path("figures/qc")

DOUBLET_RATE = 0.04
BINS         = np.linspace(0, 0.6, 61)


def plot_distributions(a, batches):
    """Observerd and simulated scores per donor, with both candidate cutoffs."""
    donors = sorted(batches)
    fig, axs = plt.subplots(2, 4, figsize = (37, 15), sharex = True, sharey = True)

    for ax, donor in zip(axs.flat, donors):
        obs  = a.obs.loc[a.obs["donor"] == donor, "doublet_score"].to_numpy()
        sim  = np.asarray(batches[donor]["doublet_scores_sim"])
        auto = float(batches[donor]["threshold"])
        rank = float(np.quantile(obs, 1 - DOUBLET_RATE))

        ax.hist(obs, bins = BINS, density = True, color = "tab:blue", 
                alpha = 0.6, label = "observed")
        ax.hist(sim, bins = BINS, density = True, histtype = "step",
                color = "tab:red", linewidth = 1.5, label = "simulated")
        ax.axvline(auto, color = "black", linestyle = "--", linewidth = 1.5,
                   label = f"auto {auto:.2f}")
        ax.axvline(rank, color = "tab:green", linewidth = 1.5,
                   label = f"top 4% {rank:.2f}")

        ax.set_yscale("log")
        ax.set_title(donor)

    axs[0, 0].legend(fontsize = 8)
    fig.supxlabel("doublet score")
    fig.supylabel("density (log)")
    fig.tight_layout()
    fig.savefig(FIG / "doublet_score_distribution.png", dpi = 150, 
                bbox_inches = "tight")

def plot_checks(a):
    """Scatter of score vs gene count, and cumulative score curve, per donor."""
    donors = sorted(a.obs["donor"].unique())
    fig, axs = plt.subplots(1, 2, figsize = (14, 6))

    axs[0].scatter(
        a.obs["n_genes_by_counts"], a.obs["doublet_score"],
        s = 2, alpha = 0.2, color = "tab:blue", edgecolors = "none",
    )
    axs[0].set_xlabel("n_genes_by_counts")
    axs[0].set_ylabel("doublet_score")
    axs[0].set_title("score vs gene count")

    for donor in donors:
        s = np.sort(a.obs.loc[a.obs["donor"] == donor, "doublet_score"])
        frac = 1 - np.arange(len(s)) / len(s)
        axs[1].plot(s, frac, linewidth = 1, label = donor)

    axs[1].axhline(DOUBLET_RATE, color = "black", linestyle = "--", 
                   linewidth = 2)
    axs[1].set_yscale("log")
    axs[1].set_xlabel("doublet_score cutoff")    
    axs[1].set_ylabel("fraction of cells above cutoff")
    axs[1].set_title("cumulative")
    axs[1].legend(fontsize = 12)

    fig.tight_layout()
    fig.savefig(FIG / "doublet_checks.png", dpi = 150, 
                bbox_inches = "tight")


def main():
    a = sc.read_h5ad(IN)
    batches = a.uns["scrublet"]["batches"]
    print("keys per batch: ", list(batches[sorted(batches)[0]]))

    FIG.mkdir(parents = True, exist_ok = True)
    plot_distributions(a, batches)
    plot_checks(a)
    plt.close("all")

    print("wrote", FIG / "doublet_score_distribution.png")
    print("wrote", FIG / "doublet_checks.png")


if __name__ == "__main__":
    main()


