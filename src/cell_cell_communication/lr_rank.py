"""Ligand-receptor ranking between cell types, computed once per resource."""

import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd 

import cell_communication_utils as ccu

SCRIPT_NAME   = "lr_rank"
TOP_N_DOTPLOT = 25
RANK_FLOOR    = 1e-6

def rank_score(values: pd.Series) -> np.ndarray:
    """Aggregate rank as a positive score: negative log10, floored at `RANK_FLOOR` so zero remain finite"""
    return -np.log10(np.clip(values.to_numpy(dtype = float), RANK_FLOOR, 1.0))

def plot_dotplot(significant: pd.DataFrame, resource_name: str) -> None:
    """Dotplot of the top interactions: cell type pair against ligand-receptor pair,sized by magnitude, colored by specificity"""
    top = significant.head(TOP_N_DOTPLOT).copy()
    if top.empty:
        print(f"{resource_name}: no significant interactions found, dotplot skipped")
        return
    
    top["cell_pair"] = top["source"].astype(str) + " -> " + top["target"].astype(str)
    top["lr_pair"]   = top["ligand_complex"].astype(str) + " -> " + top["receptor_complex"].astype(str)

    cell_pairs = list(dict.fromkeys(top["cell_pair"]))
    lr_pairs   = list(dict.fromkeys(top["lr_pair"]))
    x  = [cell_pairs.index(value) for value in top["cell_pair"]]
    y  = [lr_pairs.index(value) for value in top["lr_pair"] ]


    fig, ax = plt.subplots(figsize=ccu.FIGSIZE)
    dots = ax.scatter(
        x,
        y,
        s = rank_score(top["magnitude_rank"]) * 60,
        c = rank_score(top["specificity_rank"]),
        cmap       = "viridis",
        edgecolors = "none"
    )
    ax.set_xticks(range(len(cell_pairs)))
    ax.set_xticklabels(cell_pairs, rotation = 70)
    ax.set_yticks(range(len(lr_pairs)))
    ax.set_yticklabels(lr_pairs)
    ax.set_xlabel("sender to reciver")
    ax.set_ylabel("ligand to receptor")
    ax.set_title(f" Top {TOP_N_DOTPLOT} interactions, {resource_name}")
    bar = fig.colorbar(dots, ax = ax)
    bar.set_label("specificity, -log10 rank")

    ccu.style_dark(fig)
    ccu.save_figure(fig, f"lr_rank_dotplot_{resource_name}")

def plot_counts_heatmap(counts: pd.DataFrame, resource_name: str, vmax: int) -> None:
    """Heatmap of sender-by-receiver interaction counts, on a color scale shared across resources."""
    fig, ax = plt.subplots(figsize = ccu.FIGSIZE)
    image = ax.imshow(counts.to_numpy(), cmap = "magma", vmin = 0, vmax = vmax)
    ax.set_xticks(range(len(counts.columns)))
    ax.set_xticklabels(counts.columns, rotation = 90)
    ax.set_yticks(range(len(counts.index)))
    ax.set_yticklabels(counts.index)
    ax.set_xlabel("receiver")
    ax.set_ylabel("sender")
    ax.set_title(f"significant interactions per cell type pair, {resource_name}")
    bar = fig.colorbar(image, ax = ax)
    bar.set_label("interactions")
    
    ccu.style_dark(fig)
    ccu.save_figure(fig, f"lr_rank_heatmap_{resource_name}")
    


def main() -> None:
    """Rankig per resource: full table, significant subser, interaction counts, figures"""
    ccu.set_seed()
    adata = ccu.load_annotated()
    labels = ccu.cell_type_labels(adata)

    counts_by_resource = {}
    for resource_name in ccu.RESOURCES:
        frame       = ccu.rank_interactions(adata, resource_name)
        significant = ccu.filter_significant(frame)
        counts      = ccu.interaction_counts(significant, labels)
        counts_by_resource[resource_name] = counts

        ccu.save_table(frame, SCRIPT_NAME, f"lr_rank_{resource_name}_all")
        ccu.save_table(significant, SCRIPT_NAME, f"lr_rank_{resource_name}_significant")
        ccu.save_table(counts, SCRIPT_NAME, f"lr_rank_{resource_name}_counts", index = True)

        print(f"{resource_name}: {len(significant)} of {len(frame)} interactions significant")
        plot_dotplot(significant, resource_name)

    vmax = int(max(matrix.to_numpy().max() for matrix in counts_by_resource.values()))
    for resource_name, matrix in counts_by_resource.items():
        plot_counts_heatmap(matrix, resource_name, vmax)
        

if __name__ == "__main__":
    main()