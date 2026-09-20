"""Per-donor ligand-receptor ranking and reproducibility of  the full-dataset results."""

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd 

import cell_communication_utils as ccu 

SCRIPT_NAME   = "lr_donors"
SOURCE_SCRIPT = "lr_rank"
KEY_COLUMNS   = ["source", "target", "ligand_complex", "receptor_complex"]
BAR_WIDTH     = 0.4

def donor_subset(adata: ad.AnnData, donor:str) -> ad.AnnData:
    """One donor lane: the cells of that donor, with unused cell type categories dropped"""
    subset = adata[adata.obs[ccu.DONOR_KEY] == donor].copy()
    subset.obs[ccu.GROUP_KEY] = subset.obs[ccu.GROUP_KEY].cat.remove_unused_categories()
    return subset


def label_coverage(subset: ad.AnnData, donor: str, labels: list[str]) -> pd.DataFrame:
    """Cell count per label in one lane: the full label list, with a flag for those clearing the minimum group"""
    counts = subset.obs[ccu.GROUP_KEY].value_counts().reindex(labels, fill_value = 0)
    frame  = counts.rename_axis(ccu.GROUP_KEY).reset_index(name = "cells")
    frame.insert(0, "donor", donor)
    frame["testable"] = frame["cells"] >= ccu.MIN_CELLS
    return frame        


def rank_one_donor(subset: ad.AnnData, donor: str) -> pd.DataFrame:
    """Significant interactions in one lane, one block per resource."""
    blocks = []
    for resource_name in ccu.RESOURCES:
        frame       = ccu.rank_interactions(subset, resource_name)
        significant = ccu.filter_significant(frame)
        significant.insert(0, "donor", donor)
        significant.insert(1, "resource", resource_name)
        blocks.append(significant)
        print(f"{donor} {resource_name}: {len(significant)} of {len(frame)} significant")
    return pd.concat(blocks, ignore_index = True)


def load_full_significant(resource_name: str) -> pd.DataFrame:
    """Significant interactions from the full-dataset run for one resource."""
    path = ccu.RESULTS_DIR /SOURCE_SCRIPT / f"lr_rank_{resource_name}_significant.csv"
    return pd.read_csv(path)


def reproducibility(per_donor: pd.DataFrame, resource_name: str) -> pd.DataFrame:
    """Lane count per interaction: how many donors called it significant, and whether the full run did."""
    lanes  = per_donor[per_donor["resource"] == resource_name]
    counts = (lanes.groupby(KEY_COLUMNS, observed = True)["donor"]
              .nunique()
              .reset_index(name = "donor_significant"))

    full = load_full_significant(resource_name)[KEY_COLUMNS].copy()
    full["significant_full"] = True

    merged = counts.merge(full, on = KEY_COLUMNS, how = "outer")
    merged["donor_significant"] = merged["donor_significant"].fillna(0).astype(int)
    merged["significant_full"]  = merged["significant_full"].fillna(False).astype(bool)
    merged.insert(0, "resource", resource_name)
    return merged.sort_values("donor_significant", ascending = False).reset_index(drop = True)


def plot_reproducibility(table: pd.DataFrame, resource_name: str) -> None:
    """Reproducibility bars: interactions by lane count, split by the full-dataset result."""
    lanes = range(0, 9)
    in_full  = table[table["significant_full"]]["donor_significant"].value_counts().reindex(lanes, fill_value = 0)
    not_full = table[~table["significant_full"]]["donor_significant"].value_counts().reindex(lanes, fill_value = 0)

    positions = list(lanes)
    left      = [position - BAR_WIDTH / 2 for position in positions]
    right     = [position + BAR_WIDTH / 2 for position in positions]

    fig, ax = plt.subplots(figsize = ccu.FIGSIZE)
    ax.bar(left,  in_full.to_numpy(),  width = BAR_WIDTH, label = "significant in the full run")
    ax.bar(right, not_full.to_numpy(), width = BAR_WIDTH, label = "not significant in the full run")
    ax.set_xticks(positions)
    ax.set_xlabel("donor lanes calling the interaction significant")
    ax.set_ylabel("interactions")
    ax.set_title(f"per-donor reproducibility, {resource_name}")
    ax.legend(fontsize = ccu.LEGEND_FONTSIZE)

    ccu.style_dark(fig)
    ccu.save_figure(fig, f"lr_donors_reproducibility_{resource_name}")

def main() -> None:
    """Per-lane ranking: a table per donor, the label coverage, the reproducibility table and a figure per resources."""
    ccu.set_seed()
    adata  = ccu.load_annotated()
    labels = ccu.cell_type_labels(adata)
    donors = list(adata.obs[ccu.DONOR_KEY].cat.categories)
    print(f"{len(donors)} donor lanes: {', '.join(donors)}")

    coverage_blocks = []
    donor_blocks    = []

    for donor  in donors:
        subset = donor_subset(adata, donor)
        print(f"{donor}: {subset.n_obs} cells")
        coverage_blocks.append(label_coverage(subset, donor, labels))

        per_donor = rank_one_donor(subset, donor)
        ccu.save_table(per_donor, SCRIPT_NAME, f"lr_donors_{donor}")
        donor_blocks.append(per_donor)

    ccu.save_table(pd.concat(coverage_blocks, ignore_index = True), SCRIPT_NAME, "lr_donors_coverage") 

    per_donor_all = pd.concat(donor_blocks, ignore_index = True)
    reproducible_blocks = []

    for resource_name in ccu.RESOURCES:
        table = reproducibility(per_donor_all, resource_name)
        reproducible_blocks.append(table)
        plot_reproducibility(table, resource_name)

        held = table[table["significant_full"] & (table["donor_significant"] > 0)]
        full = table["significant_full"].sum()
        print(f"{resource_name}: {len(held)} of {full} full-run interactions found in at least one lane")

    ccu.save_table(pd.concat(reproducible_blocks, ignore_index = True), SCRIPT_NAME, 
                   "lr_donors_reproducible")   



if __name__ == "__main__":
    main()