"""Pathway aggregation of the ranked ligand-receptor interactions."""

import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import pandas as pd

import cell_communication_utils as ccu

SCRIPT_NAME     = "lr_pathways"
LIBRARIES       = ("KEGG_2021_Human", "Reactome_Pathways_2024")
COMPLEX_SEP     = "_"
MIN_SET_SIZE    = 5
MAX_SET_SIZE    = 500
TOP_PATHWAYS    = 25
MIN_HITS        = 3
LABEL_WIDTH     = 45
FIGURE_RESOURCE = "consensus"

ASSIGNMENT_COLUMNS = [
    "library",
    "pathway",
    "set_size",
    "source",
    "target",
    "ligand_complex",
    "receptor_complex",
    "magnitude_rank",
    "specificity_rank",
    "significant"
]


def load_ranked(resource: str) -> pd.DataFrame:
    """
        Full interaction table for one resource: every tested ligand-receptor
        pair with its per-method scores and merged ranks.
    """
    path = ccu.RESULTS_DIR / "lr_rank" / f"lr_rank_{resource}_all.csv"
    return pd.read_csv(path)


def significant_flag(frame: pd.DataFrame) -> pd.Series:
    """
        Boolean flag per interaction: true where both merged ranks fall at or
        below their cutoffs.
    """
    magnitude   = frame["magnitude_rank"] <= ccu.MAGNITUDE_CUTOFF
    specificity = frame["specificity_rank"] <= ccu.SPECIFICITY_CUTOFF
    return magnitude & specificity


def sized_gene_set(library: str) -> dict[str, set[str]]:
    """
        Gene sets of one library keyed by term: members as a set, sizes
        restricted to the fixed bounds.
    """
    kept: dict[str, set[str]] = {}
    for term, genes in ccu.load_gene_sets(library).items():
        members = set(genes)
        if MIN_SET_SIZE <= len(members) <= MAX_SET_SIZE:
            kept[term] = members
    return kept


def term_index(gene_sets: dict[str, set[str]]) -> dict[str, set[str]]:
    """Gene-to-term index over the gene sets of one library."""
    index: dict[str, set[str]] = {}
    for term, members in gene_sets.items():
        for gene in members:
            index.setdefault(gene, set()).add(term)
    return index


def matching_terms(index: dict[str, set[str]], entry: str) -> set[str]:
    """Terms matching any subunit of a ligand or receptor entry."""
    terms: set[str] = set()
    for subunit in str(entry).split(COMPLEX_SEP):
        terms |= index.get(subunit, set())
    return terms


def assign_pathways(frame: pd.DataFrame, library: str) -> pd.DataFrame:
    """
        Interaction-to-term assignments for one library: one row per matched
        interaction and term, with the term size, the cell type pair and both
        merged ranks.
    """
    gene_sets = sized_gene_set(library)
    index = term_index(gene_sets)
    columns = [
        "source",
        "target",
        "ligand_complex",
        "receptor_complex",
        "magnitude_rank",
        "specificity_rank",
        "significant",
    ]
    rows = []
    for record in frame[columns].itertuples(index = False):
        ligand_terms = matching_terms(index, record.ligand_complex) # type: ignore
        receptor_terms = matching_terms(index, record.receptor_complex) # type: ignore
        for term in sorted(ligand_terms & receptor_terms):
            rows.append(
                {
                    "library"         : library,
                    "pathway"         : term,
                    "set_size"        : len(gene_sets[term]),
                    "source"          : record.source,
                    "target"          : record.target,
                    "ligand_complex"  : record.ligand_complex,
                    "receptor_complex": record.receptor_complex,
                    "magnitude_rank"  : record.magnitude_rank,
                    "specificity_rank": record.specificity_rank,
                    "significant"     : record.significant,
                }
            )
    return pd.DataFrame(rows, columns = ASSIGNMENT_COLUMNS)


def assign_all(frame: pd.DataFrame) -> pd.DataFrame:
    """
        Interaction-to-term assignments across every configured library,
        stacked.
    """
    parts = [assign_pathways(frame, library) for library in LIBRARIES]
    return pd.concat(parts, ignore_index = True)


def pathway_summary(assignment: pd.DataFrame, base_rate: float) -> pd.DataFrame:
    """
        Per-term totals: interactions tested and significant, the significant
        share and its fold over the resource-wide rate, and the distinct
        pairs, senders, receivers and ligands behind them.
    """
    records = []
    for (library, pathway), group in assignment.groupby(["library", "pathway"], sort = False):
        hits = group[group["significant"]]
        share = len(hits) / len(group)
        records.append(
            {
                "library"          : library,
                "pathway"          : pathway,
                "set_size"         : int(group["set_size"].iloc[0]),
                "n_tested"         : len(group),
                "n_significant"    : len(hits),
                "significant_share": share,
                "fold_enrichment"  : share / base_rate,
                "n_lr_pairs"       : len(hits[["ligand_complex", "receptor_complex"]].drop_duplicates()),
                "n_senders"        : hits["source"].nunique(),
                "n_receivers"      : hits["target"].nunique(),
                "ligands"          : ";".join(sorted(hits["ligand_complex"].unique())),
            }
        )
    summary = pd.DataFrame(records)
    return summary.sort_values(
        ["fold_enrichment", "n_significant"], ascending = False, ignore_index = True
    )


def collapse_duplicate_pathways(
    assignment: pd.DataFrame, summary: pd.DataFrame
) -> pd.DataFrame:
    """
        Summary reduced to one term per distinct set of significant
        interactions, labelled by the smallest gene set holding that set.
    """
    hits = assignment[assignment["significant"]]
    interaction = (
        hits["source"] + "|" + hits["target"] + "|"
        + hits["ligand_complex"] + "|" + hits["receptor_complex"]
    )
    tagged = hits.assign(interaction = interaction)

    fingerprints: dict[tuple[str, str], frozenset[str]] = {}
    for (library, pathway), group in tagged.groupby(["library", "pathway"], sort = False):
        fingerprints[(str(library), str(pathway))] = frozenset(group["interaction"])

    ordered = summary[summary["n_significant"] > 0].sort_values("set_size")
    kept: set[tuple[str, str]] = set()
    seen: set[frozenset[str]] = set()
    for record in ordered.itertuples(index = False):
        mark = fingerprints[(record.library, record.pathway)] # type: ignore
        if mark in seen:
            continue
        seen.add(mark)
        kept.add((record.library, record.pathway)) # type: ignore

    keep_row = [(row.library, row.pathway) in kept for row in summary.itertuples(index = False)]
    return summary[keep_row].reset_index(drop = True)


def shorten(name: str) -> str:
    """Term name trimmed to a fixed width for axis labels."""
    if len(name) <= LABEL_WIDTH:
        return name
    return name[: LABEL_WIDTH - 1] + "..."


def sender_pathway_matrix(
    assignment: pd.DataFrame, labels: list[str], pathways: list[str]
) -> pd.DataFrame:
    """
        Sender-by-pathway matrix of significant interaction counts, reindexed
        to fixed label and pathway lists so runs share their axes.
    """
    hits = assignment[assignment["significant"]]
    counts = pd.crosstab(hits["source"], hits["pathway"])
    return counts.reindex(index = labels, columns = pathways, fill_value = 0)


def plot_sender_pathways(matrix: pd.DataFrame) -> Figure:
    """
        Heatmap of the significant interactions sent: cell type against
        pathway, coloured by count.
    """
    fig, ax = plt.subplots(figsize = ccu.FIGSIZE)
    image = ax.imshow(matrix.to_numpy(), cmap = "magma", aspect = "auto", vmin = 0)
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(
        [shorten(name) for name in matrix.columns],
        rotation = 90,
        fontsize = ccu.LEGEND_FONTSIZE,
    )
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize = ccu.LEGEND_FONTSIZE)
    ax.set_xlabel("pathway, ordered by enrichment")
    ax.set_ylabel("sender")
    ax.set_title(f"Significant interactions per sender and pathway ({FIGURE_RESOURCE})")
    bar = fig.colorbar(image, ax = ax, fraction = 0.025)
    bar.set_label("interactions", color = ccu.DARK_FOREGROUND)
    bar.ax.tick_params(colors = ccu.DARK_FOREGROUND)
    ccu.style_dark(fig)
    fig.tight_layout()
    return fig


def main() -> None:
    """
        Pathway assignment for both resources, the combined per-term summary,
        and the sender heatmap for the reference resource.
    """
    ccu.set_seed()

    ranked_tables: dict[str, pd.DataFrame]     = {}
    assignment_tables: dict[str, pd.DataFrame] = {}
    summary_tables: dict[str, pd.DataFrame]    = {}

    for resource in ccu.RESOURCES:
        ranked = load_ranked(resource)
        ranked["significant"] = significant_flag(ranked)
        assignments = assign_all(ranked)

        hits = assignments[assignments["significant"]].reset_index(drop = True)
        ccu.save_table(hits, SCRIPT_NAME, f"{SCRIPT_NAME}_{resource}")

        base_rate = float(ranked["significant"].mean())
        summary = pathway_summary(assignments, base_rate)
        summary.insert(0, "resource", resource)

        ranked_tables[resource]     = ranked
        assignment_tables[resource] = assignments
        summary_tables[resource]    = summary

        significant = ranked[ranked["significant"]]
        key = ["source", "target", "ligand_complex", "receptor_complex"]
        mapped = len(hits[key].drop_duplicates())
        print(
            f"{resource}: {len(significant)} significant, {mapped} mapped to at least one pathway, "
            f"{int(summary['n_significant'].gt(0).sum())} pathways with at least one hit"
        )

    combined = pd.concat([summary_tables[name] for name in ccu.RESOURCES], ignore_index = True)
    ccu.save_table(combined, SCRIPT_NAME, f"{SCRIPT_NAME}_summary")

    ranked      = ranked_tables[FIGURE_RESOURCE]
    assignments = assignment_tables[FIGURE_RESOURCE]
    summary     = summary_tables[FIGURE_RESOURCE]

    labels = sorted(set(ranked["source"]).union(ranked["target"]))
    columns = collapse_duplicate_pathways(assignments, summary)
    columns = columns[columns["n_significant"] >= MIN_HITS]
    pathways = (
        columns.sort_values("fold_enrichment", ascending = False)["pathway"]
        .head(TOP_PATHWAYS)
        .tolist()
    )
    print(f"{FIGURE_RESOURCE}: {len(pathways)} pathways kept for the figure")

    matrix = sender_pathway_matrix(assignments, labels, pathways)
    ccu.save_figure(plot_sender_pathways(matrix), f"{SCRIPT_NAME}_heatmap")


if __name__ == "__main__":
    main()