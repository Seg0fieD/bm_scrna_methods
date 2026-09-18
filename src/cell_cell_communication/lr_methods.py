"""Comparison of the individual LIANA scoring methods against the merged rank."""

import matplotlib.pyplot as plt
import pandas as pd

import cell_communication_utils as ccu

SCRIPT_NAME   = "lr_methods"
SOURCE_SCRIPT = "lr_rank"
TOP_N         = 100  

AGGREGATE_COLUMNS = {
    "magnitude_rank"   : ("Merged", "magnitude", False),
    "specificity_rank" : ("Merged", "specificity", False),
}

KEY_COLUMNS = ["source", "target", "ligand_complex", "receptor_complex"]


def load_ranking(resource_name: str) -> pd.DataFrame:
    """
        Saved ranking table for one resource: every tested interaction with its per-method
        scores and merged ranks.
    """
    path  = ccu.RESULTS_DIR / SOURCE_SCRIPT / f"lr_rank_{resource_name}_all.csv"
    frame = pd.read_csv(path)
    print(f"loaded {path}: {len(frame)} interactions")
    return frame

def oriented_scores(frame: pd.DataFrame)  -> pd.DataFrame:
    """Per-method scores on one orientation, grouped by score type: higher is always better, flipping the columns where it is not."""
    columns = dict(ccu.present_score_columns(frame))
    for name , spec in AGGREGATE_COLUMNS.items():
        if name in frame.columns:
            columns[name] = spec

    oriented = {}
    for column, (method, score_type, higher_is_better) in columns.items():
        values = frame[column].astype(float)
        oriented[f"{method} {score_type}"] = values if higher_is_better else -values
    scores = pd.DataFrame(oriented, index = frame.index)    
    return scores[block_order(list(scores.columns))]

def block_order(columns: list[str]) -> list[str]:
    """Column order grouping the magnitude scores before the specificity scores, each followed by its merged rank."""
    magnitude          = [c for c in columns if c.endswith("magnitude") and not c.startswith("Merged")]
    specificity        = [c for c in columns if c.endswith("specificity") and not c.startswith("Merged")]
    merged_magnitude   = [c for c in columns if c.startswith("Merged") and c.endswith("magnitude")]
    merged_specificity = [c for c in columns if c.startswith("Merged") and c.endswith("specificity")]
    return magnitude + merged_magnitude + specificity + merged_specificity

def method_correlation(oriented: pd.DataFrame)  -> pd.DataFrame:
    """Spearman correlation between every pair of oriented score columns, grouped by score type.""" 
    return oriented.corr(method = "spearman")

def top_pairs(frame: pd.DataFrame, oriented: pd.DataFrame, resource_name: str) -> pd.DataFrame:
    """Highest-scoring interactions per method: one block per column, in descending order."""
    blocks = []
    for column in oriented.columns:
        order = oriented[column].sort_values(ascending = False).head(TOP_N)
        block =frame.loc[order.index, KEY_COLUMNS].copy()
        block.insert(0, "resource", resource_name)
        block.insert(1, "method", column)
        block.insert(2, "position", range(1, len(block) + 1))
        block["score"] = order.to_numpy()
        blocks.append(block)
    return pd.concat(blocks, ignore_index = True)

def agreement(oriented: pd.DataFrame, resource_name: str) -> pd.DataFrame:
    """Agreement with the merged ranks: Spearman over all interactions, plus shared top-list membership."""
    merged = [column for column in oriented.columns if column.startswith("Merged")]
    rows = []
    for column in oriented.columns:
        if column in merged:
            continue
        top = set(oriented[column].sort_values(ascending = False).head(TOP_N).index) 
        row = {"resource": resource_name, "method": column}
        for reference in merged: 
            label = reference.replace("Merged ", "")
            reference_top = set(oriented[reference].sort_values(ascending = False).head(TOP_N).index)
            row[f"spearman_{label}"] = oriented[column].corr(oriented[reference], method = "spearman") # type: ignore
            row[f"shared_top{TOP_N}_{label}"] = len(top & reference_top) # type: ignore
        rows.append(row)
    return pd.DataFrame(rows)


def annotation_color(value: float) -> str:
    """Text color for one correlation cell: white on saturated cells, black on pale ones"""
    return "white" if abs(value) > 0.6 else "black"

def plot_correlation(matrix: pd.DataFrame, resource_name: str) -> None:
    """Method correlation heatmaps, one panel per resource, on a fixed scale from -1 to 1"""
    fig, ax = plt.subplots(figsize = ccu.FIGSIZE)
    image = ax.imshow(matrix.to_numpy(), cmap = "coolwarm", vmin = -1, vmax = 1)
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_xticklabels(matrix.columns, rotation = 75)
    ax.set_yticks(range(len(matrix.index)))
    ax.set_yticklabels(matrix.index)
    ax.set_title(f"scoring method agreement, {resource_name}")

    for i in range(len(matrix.index)):
        for j in range(len(matrix.columns)):
            value = matrix.iat[i, j]
            ax.text(j, i, f"{value:.2f}", ha = "center", va = "center",
                    fontsize = 6, color = annotation_color(value)) # type: ignore

    bar = fig.colorbar(image, ax = ax,  shrink = 0.7)
    bar.set_label("Spearman correlation")

    ccu.style_dark(fig)
    ccu.save_figure(fig, f"lr_methods_agreement_{resource_name}")

def main() -> None:
    """Method comparison per resources: correlations, top_list, agreement with the merged rank, figure."""

    pair_blocks = []
    correlation_blocks = []
    agreement_blocks   = []

    for resource_name in ccu.RESOURCES:
        frame = load_ranking(resource_name)
        oriented = oriented_scores(frame)
        matrix = method_correlation(oriented)

        plot_correlation(matrix, resource_name)

        block = matrix.copy()
        block.insert(0, "resource", resource_name)
        block.insert(1, "method", matrix.index)
        correlation_blocks.append(block)

        pair_blocks.append(top_pairs(frame, oriented, resource_name))
        agreement_blocks.append(agreement(oriented, resource_name))

        print(f"{resource_name}: {len(oriented.columns)} score columns compared")

    ccu.save_table(pd.concat(correlation_blocks, ignore_index = True), SCRIPT_NAME, "lr_methods_correlation")
    ccu.save_table(pd.concat(pair_blocks, ignore_index = True), SCRIPT_NAME, "lr_methods_top_pairs")
    ccu.save_table(pd.concat(agreement_blocks, ignore_index = True), SCRIPT_NAME , "lr_methods_agreement")


if __name__ == "__main__":
    main()

