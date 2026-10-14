# Choosing a dimensionality reduction for bone marrow single-cell data

## The question

Single-cell data has ~20,000 genes per cell. Almost every analysis — clustering,
trajectories, neighbour graphs — runs on a compressed version of it, usually 30
numbers per cell. Which compression method should this dataset use?

PCA is the default everywhere. Several newer methods claim to do better, mostly by
handling count noise properly or by removing donor-to-donor differences. This
compares seven of them on real data and picks one.

## The data

Human Cell Atlas Census of Immune Cells (Hay 2018). Eight bone marrow donors,
one sequencing lane each. After quality control and doublet removal: **38,430 cells,
20,455 genes, 19 annotated cell types**, ranging from 10,765 CD4 T cells down to
47 dividing T cells.

## What was compared

| method | what it does |
|---|---|
| PCA | the standard baseline, on log-normalised counts |
| Pearson residuals + PCA | corrects count noise with a closed-form formula |
| GLM-PCA | fits a Poisson model directly to raw counts |
| scVI | a neural network that learns the structure and removes donor effects |
| Diffusion map on PCA | stretches the data along continuous gradients |
| Diffusion map on scVI | the same, on a donor-corrected base |
| Consensus NMF | splits the data into 30 gene programs; cells are mixtures of them |

Everything downstream was held identical — same 2,000 variable genes, 30
dimensions, 15 neighbours, same clustering resolution. Any difference in the
results comes from the compression method and nothing else.

Scoring used `scib-metrics`, which asks two things: did the method keep the
biology (do its clusters match the annotated cell types, did rare types survive),
and did it mix the donors (or do cells group by who they came from).

## What was found

**1. The methods are barely distinguishable.** All seven total scores fall between
0.648 and 0.716. PCA and scVI tie to three decimal places. Choosing between them
is not where the gains in this project are.

**2. The donor effect is weak, contrary to the original plan.** The project was
designed assuming eight donors would produce a strong batch effect. They do not.
On a plain PCA plot the donors are fully mixed, and the batch metrics agree. Only
two small patches of donor structure are visible, and neither separates a cell type.

**3. Correcting a batch effect that is not there costs biology.** scVI mixes the
donors best of any method — and scores worst on biology, with the lowest
cluster-agreement number in the whole comparison. Its UMAP shows the T cell and
monocyte regions pulled toward each other. Plotting biology against batch
correction gives a straight line with a negative slope: across these methods, one
is bought with the other.

**4. Diffusion maps score highest but are not competitors.** They are built on top
of PCA and scVI rather than replacing them, and at a fixed clustering resolution
they shatter into 46 and 51 groups against 19–33 for everything else. That
fragmentation is the same property that inflates their scores. They belong in
trajectory analysis, where the stretching is the point.

**5. Consensus NMF sees things clustering cannot.** Fourteen of its 30 programs
match known cell types without being told what they are. Five are not cell types
at all but *activities* — DNA replication, mitosis, heat shock, inflammation —
which clustering cannot separate from identity, because a dividing T cell is still
a T cell and lands in the T cell cluster either way. Three more programs are real
cell types the annotation missed: plasma cells, mast cells and macrophages, each
rare but with cells that use the program almost exclusively.

## One honest limitation

PCA scores best on biology, and part of that is circular. The 19 cell type labels
were themselves produced by clustering on PCA, so the biology metrics compare PCA
against its own output under different names.

Two things limit how much this inflates it. An independent classifier (CellTypist)
agreed with 17 of the 19 labels without using the PCA embedding, so the labels
reflect real biology. And Pearson residuals and GLM-PCA land within 0.03 of PCA
despite having played no part in making the labels — if the lead were pure
self-grading, that gap would be far wider.

The right reading is that PCA is **no worse than the alternatives**, not that it is
clearly better.

## The decision

**Use PCA.** Best biology score among the genuine alternatives, computed in
seconds, nothing to tune, and this dataset does not need batch correction.

The others are kept for specific jobs rather than discarded:

| method | kept for |
|---|---|
| PCA | the default, for clustering and anything needing a representation |
| Consensus NMF | pathway analysis — its programs are ready-made gene sets |
| Diffusion map on PCA | trajectory and pseudotime work |
| scVI | a check, whenever some result looks like it might be donor-driven |

The 19 cell type labels stay fixed from here. Later analyses may compute on a
different representation, but they do not re-cluster or re-annotate — otherwise
their results stop being comparable with each other.

## Practical notes for anyone re-running this

- **GLM-PCA needs `penalty=10`.** The library default diverges within four minutes
  on this data. Very high-count genes such as haemoglobin overflow under the
  Poisson log link. With the higher penalty it converges in 87 iterations, taking
  about 1 hour 45 minutes on four CPU cores. There is no GPU option.
- **scVI is the fast one** — 9 minutes on a 4 GB GPU with a batch size of 128.
  Increasing the data-loader worker count made it slower, not faster, on a
  four-core machine.
- **Consensus NMF takes 24 minutes** for 30 programs across 20 repeats, and its
  results can be re-read from disk afterwards without refitting.
- **Everything was run on all 38,430 cells.** The original plan called for a
  20–30k subsample; it was dropped because the rarest cell types (103 and 47
  cells) would have been cut too thin to score honestly.

## Files

- Embeddings: `results/benchmark/embeddings/*.npy`
- Scores: `results/benchmark/scib_scores.csv`
- Figures: `figures/benchmark/`
- Code: `src/benchmark/`
