# Trajectory analysis — stage note

Plain-English record of what this stage did, how it was done, what came out, and
what the results mean. It covers both passes: the first pass (PAGA, Palantir,
CellRank, driver genes, gene trends, diagnostics) and the reopened second pass
(backbone figures, principal tree, RNA velocity). Figures, tables and the
number-by-number reasoning are in `claude/report_trajectory.md`. This note is
the version to read first.

---

## 1. What this stage was for

Bone marrow contains stem cells that turn into blood cells. A trajectory method
tries to put single cells in the order they develop, using only their gene
expression — no time points, no labels. The goal here was to practise four
things: ordering cells along a path, finding where the paths end, finding the
genes that drive each path, and checking whether any of it is real.

One decision was made before any code was written, and it shapes the whole
stage. Bone marrow samples contain two very different kinds of cell:

- **Cells that are developing here and now** — stem cells, erythroid precursors,
  monocyte precursors, pre-B cells.
- **Mature cells that finished developing elsewhere and are just passing
  through** — mature T cells, NK cells, mature B cells, mature monocytes.

A trajectory through the second group is meaningless. They did not come from the
stem cells in this sample. The choice was to **keep them in anyway and report
what happens as a negative result**, rather than quietly removing them. A method
that cannot tell the difference is worth demonstrating, not hiding.

**Why the stage was reopened.** After the first pass, reading the OSCA Advanced
trajectory chapter showed three standard views that had not been done: the
connection tree drawn on the real map of cells, a principal tree (a second
ordering method with its own geometry), and RNA velocity. Velocity mattered
most, because it is **the only method in the stage whose direction comes from
the data itself rather than from a starting cell we chose.**

---

## 2. Input and fixed settings

| item | value |
|---|---|
| data | `data/processed/bm_annotated.h5ad`, 38,430 cells × 20,455 genes |
| labels | 19 annotated cell types, fixed from the earlier stage |
| donors | 8 (one sequencing lane each, MantonBM1–8) |
| main embedding | diffusion map built on PCA (`diffmap_pca`) |
| cross-check embedding | diffusion map built on scVI (`diffmap_scvi`) |
| root | HSPC cell with the highest stem-marker score (CD34, AVP, CRHBP, SPINK2) |
| velocity input | the original aligned reads (BAM files), one per lane, about 74 GB |
| velocity reference | GENCODE v27 on GRCh38 — the exact annotation the reads were aligned with — plus the UCSC hg38 repeat mask |
| seed | 7 everywhere |
| figures | 300 dpi, dark background, UMAP legend font 8 |

Nothing in this stage writes back into the input file.

**Environments.** Everything runs in `sc_core` except velocyto, which counts
spliced and unspliced reads and runs in its own `rna_velo` environment. scVelo
0.3.4 needs **NumPy below 2.4**; this pin is in `envs/constraints.txt`.

---

## 3. What was run, in order

Each script is a plain file in `src/trajectory/`, run from the project root.
Shared paths and helpers live in `trajectory_utils.py`.

### First pass

**1. `paga_topology.py` — which cell types are connected to which**

PAGA is the cheap first look. It does not order cells; it collapses each cell
type to a single node and measures how strongly the cell types touch each other
in the neighbour graph. Run on both diffusion maps to see whether the answer
depends on the embedding. Output: a network graph, a connectivity heatmap and a
minimum spanning tree per embedding.

**2. `pt_palantir.py` — pseudotime**

Palantir places every cell on a single number between 0 and 1, the "pseudotime",
based on how far it is from the root cell through the diffusion map. Run four
times: two embeddings × two ways of picking the root (stem markers, or the
extreme of the first diffusion component). Four runs so the result could be
compared rather than trusted.

**3. `pt_cellrank.py` — where the paths end, and where each cell is heading**

CellRank takes the pseudotime, builds a transition matrix (80% pseudotime
direction, 20% neighbour connectivity), and finds the states cells get stuck in.
Those are the terminal states. It then gives every cell a probability of
reaching each terminal. Run at 8, 5 and 4 terminal states to see which results
survive the choice.

**4. `fate_drivers.py` — which genes go with which destination**

For every gene and every terminal, the correlation between the gene's expression
and the probability of reaching that terminal. No gene filtering by identity —
haemoglobin and ribosomal genes were left in on purpose, because removing them
would delete the erythroid answer.

**5. `gene_trends.py` — when those genes switch on**

Correlation says which genes go with a path. It says nothing about order. This
step takes the cells committed to each lineage, sorts them by pseudotime, bins
them into 40 groups of equal size, and plots each gene's average across bins.
Underneath each plot is a bar showing which cell-type labels occupy each bin —
this is the actual test, described in section 5.

**6. `lineage_diagnostics.py` — checking two things that looked wrong**

A follow-up script, not part of the original plan. Two specific results did not
make sense and were taken back to the raw counts.

**7. `trajectory_report.py` — pulling it together**

No new computation. Reads the saved outputs and produces the stage-level
figures: pseudotime per cell type, a fate matrix of all 19 cell types against
all 8 terminals, and a direct comparison of Palantir against CellRank.

### Second pass

**8. `trajectory_backbone.py` — the connection tree on the real map**

No new computation except one arrow field. Draws the PAGA spanning tree between
the centre of each cell type on the UMAP, so every connection gets a real
length; plots Palantir's own fate entropy as a distribution; and draws the
CellRank transition matrix as arrows on the map. The arrow figure is captioned
"not RNA velocity", because its direction comes from our chosen root.

**9. `principal_graph.py` — a second ordering with different geometry**

scFates fits a tree of 50 points through the cloud of cells and orders each cell
by its distance along the tree from the root. A completely different idea from
Palantir, so it is a fair test of whether the ordering is a property of the data
or of one method. A first fit with more points produced over 30 tree tips and
crashed; two parameter sweeps (12 settings) showed the tree cannot go below
about ten tips on this data. Final setting: 50 points, length penalty 1000.

**10. `velocity_barcodes.py` — one cell list per lane**

velocyto reads one BAM at a time, so the 38,430 kept cells were split into eight
barcode lists, one per lane. Nothing is computed, but the split exposed how
unevenly the cell types are spread across lanes.

**Between scripts 10 and 11 (shell, not Python):** the reference files and the
eight BAMs were downloaded and checksum-verified, the first BAM's header was
read to confirm the aligner and annotation, and `velocyto run` was run on each
lane against its barcode list, with the repeat mask. About 49 minutes per lane,
6.5 hours for all eight.

**11. `velocity_counts.py` — merging and checking the counts**

Merges the eight velocyto outputs, renames cells to match the annotated object,
and checks recovery per lane, per cell type and per lane × cell type before any
model is fitted. Also measures the unspliced fraction per cell, per cell type
and per lane — the signal velocity is built on.

**12. `velocity_scvelo.py` — RNA velocity on all cells**

scVelo with two models on all 38,430 cells: the simple steady-state
("stochastic") model and the full kinetic ("dynamical") model. Produces velocity
pseudotime, latent time, streamlines on the map, velocity-directed PAGA arrows
between cell types, and a per-lane check.

**13. `velocity_progenitor.py` — the same on the progenitor side only**

A sensitivity run, kept as a separate script with separately named outputs so
the whole-sample result stays intact. Removes the seven mature lymphoid labels
(16,885 cells left), recomputes genes and graph on the subset, fits the kinetic
model on all 2,000 variable genes, lets genes behave differently in different
cell types ("differential kinetics"), and roots the orderings at the HSPC cell
at Palantir zero. Also writes a direct comparison with the whole-sample arrows
and a forward-against-reverse check on eight known differentiation steps.

---

## 4. The measurements used, and what each one means

**Pseudotime** — one number per cell, 0 at the root. It is a *distance*, not a
time. This distinction turns out to matter more than anything else in the stage.

**Fate probability** — for each cell, a probability of reaching each terminal.
The eight numbers add to 1.

**Fate entropy (CellRank)** — how spread out a cell's eight probabilities are,
scaled to run from 0 to 1. Entropy 0 means the cell is going to exactly one
place. Entropy 1 means all eight are equally likely. This is the single number
that separates a progenitor from a finished cell without using any label.

**Fate entropy (Palantir)** — the same idea over Palantir's own four branches,
unscaled, with a ceiling of 1.386. **A different quantity from the CellRank
entropy**; the two cannot be compared.

**Gene–fate correlation** — how strongly a gene's expression tracks the
probability of one destination. Useful, but blind to order.

**Peak bin** — the pseudotime bin where a gene is at its highest. This is what
tests order.

**Bin composition** — which cell-type labels fall in each pseudotime bin. **The
labels were never given to Palantir or CellRank.** If the ordering is real, the
labels should sort themselves along the axis. This is the strongest test in the
stage, and it is a test the methods can fail.

**Method agreement** — the correlation between two methods' outputs. When the
methods share only the expression matrix, anything they both find is a property
of the data.

**Edge length** — in the backbone figures, the distance on the map between two
connected cell types. A strong connection between two cell types that sit far
apart is a warning sign.

**Unspliced fraction** — the share of a cell's RNA molecules that still contain
introns, i.e. were made recently. Velocity reads direction from it: more
unspliced RNA than expected means a gene is switching on.

**Velocity** — for each cell, the direction its expression is predicted to move
next. Drawn as streamlines on the map.

**Transition confidence** — how strongly velocity points from one cell type into
another. Can be ranked within a model; it is not a probability and can exceed 1.

**Velocity confidence** — whether a cell's velocity agrees with its neighbours'.
**It measures smoothness, not correctness.** A whole island pointing the wrong
way together scores perfectly.

**Latent time** — the kinetic model's own internal clock, built from the genes
it fits well.

**Gene likelihood** — how well the kinetic model fits one gene, 0 to 1. Below
about 0.1 the gene is not really described by the model.

---

## 5. What worked

**The root is right, and by a wide margin.** HSPC has a median pseudotime of
0.129. The next population is at 0.404. Nothing else is below 0.4. The stem-cell
compartment is clearly at the start of the axis.

**Stem cell multipotency came out as a number.** HSPC has a fate entropy of
0.804 and its most likely destination carries only 0.335. Its probabilities are
close to flat across six terminals — NK 0.233, monocyte 0.195, B naive 0.188,
Pre-B large 0.182, plasmacytoid DC 0.093, erythroid 0.055 and 0.046. Nothing
told the method that a stem cell can become many things. It worked that out from
the transition matrix. This is the best positive result of the stage.

**Entropy ranks every population by how "undecided" it is, and the ranking is
biologically correct.** Stem cells highest at 0.804. Then the intermediates —
dendritic cell 0.439, granulocyte precursor 0.414. Then the finished cells at
the bottom — NK 0.046, B naive 0.087. Nobody supplied this ordering.

**One lineage passes the label test end to end.** In the monocyte lineage, bins
0–6 are granulocyte precursor and dendritic cell, monocyte classical takes over
at bin 7, and **bins 38–39 are 100% monocyte non-classical**. Precursor →
classical → non-classical, recovered from unlabelled data. The genes match:
LYZ peaks at bin 3, the alarmins S100A8/S100A9 at 16, VCAN at 24, CD14 and FCN1
at 30, and the mature markers at 39.

**Velocity confirms the monocyte direction, independently of any root.**
Granulocyte precursor → monocyte classical points the right way in all four
velocity fits (two models × two runs), and monocyte classical → non-classical in
three of four. The kinetic model is most confident exactly here, and its latent
time agrees with Palantir in the myeloid lineage (+0.68). This is the one place
in the stage where pseudotime, fate, gene order *and* velocity all agree. **The
monocyte lineage is the stage's most solid result.**

**Three cell-type orderings come out right in the pseudotime medians.**

- granulocyte precursor 0.530 → monocyte classical 0.734 → non-classical 0.744
- Pre-B large 0.721 → Pre-B small 0.832 → B naive 0.847
- erythroid early 0.622 → erythroid late 0.673

**The early B step is confirmed by velocity too.** Pre-B large → Pre-B small
points forward in all four velocity fits, and Pre-B small → B naive in three of
four. The arrows are weak (0.02–0.20) but never reversed.

**The erythroid gene order is causally correct.** GATA1 peaks at bin 31, KLF1 at
35, then the heme enzymes and red-cell membrane proteins at 37–39, and the
haemoglobins last. The transcription factors peak *before* the genes they
control. A correlation table could never have shown this.

**Five lineages return textbook driver genes**, including the right master
regulators, from a correlation and nothing else: monocyte (LYZ at r = 0.96,
plus SPI1, KLF4, CEBPD, MAFB), erythroid (GATA1, KLF1, GFI1B, NFIA plus heme,
membrane and globin genes), plasmacytoid DC (LILRA4, CLEC4C, and TCF4 — the
master regulator of the cell type), B naive (IGHM, CD79A, MS4A1), and platelet
(TUBB1, GP9, PF4, PPBP).

**Four terminal states are confirmed by two independent methods.** Palantir and
CellRank work completely differently and share only the input matrix. Both find
monocyte classical (correlation 0.98), B naive (0.98), platelet (0.97) and NK
(0.94). These four are in the data, not in an algorithm.

**The principal tree confirms several first-pass findings with completely
different geometry.** It agrees with Palantir on the erythroid lineage (+0.89).
Seven of CellRank's eight terminals have a tree tip. It treats Pre-B large as a
point the B cells pass through, not an end — agreeing that the CellRank Pre-B
large terminal is not a real destination. And it splits erythroid early and late
into two sibling tips, the same split CellRank produced, so that split is in the
data, not in one method.

**The velocity counts are complete and clean.** Every one of the 38,430 cells
was recovered, in all 148 lane × cell-type groups. The unspliced fraction falls
out by compartment with no prompting: mature lymphoid 0.095–0.113, myeloid and
HSPC 0.061–0.084, erythroid and platelet 0.045–0.054.

**Velocity positions do not depend on the lane.** Across 106 lane × cell-type
groups in the whole-sample run, no velocity pseudotime shifts by 0.1 or more
with the lane. The orderings are not a technical artefact of one donor.

---

## 6. What did not work, and why

**The mature cells were ordered anyway, exactly as predicted.**

The seven mature lymphoid populations sit at pseudotime 0.58 to 0.71 — the
*middle* of the axis, where terminal cells should not be. Their position is not
developmental time. It is diffusion distance from a stem cell root to a
population that never came from it, and a moderate distance to a separate island
reads as a mid-axis number.

The NK lineage's bin composition makes this unmistakable. The bins sort into:
CD8 T memory → CD4 T naive → CD4 T memory and T activated → CD8 T cytotoxic →
NK. Every step of that is wrong as development. A naive CD4 T cell does not
become a cytotoxic CD8 T cell, and neither becomes an NK cell. The genes agree
with the composition: IL7R and CD3D (naive/memory T) peak at bin 21, the
cytotoxic genes NKG7, GZMA, CST7, KLRB1 all peak at bin 39. What the axis found
is a gradient of cytotoxic-programme strength across seven unrelated mature
populations.

And in the fate matrix, every one of those seven labels assigns **0.86 to 0.99
of its probability to a terminal named "NK"** — NK 0.986, CD8 T cytotoxic 0.966,
CD8 T memory 0.932, T activated 0.897, CD4 T memory 0.893, CD4 T naive 0.882, T
proliferating 0.859. That terminal should be read as "mature lymphoid", not
"NK". It is named after whichever cell the macrostate happened to land on.

Its driver genes say the same: the top of the list is IL32, CD3D, CD3E, TRAC,
CD7 — T cell genes — with the actually NK-specific genes (NKG7, GZMA, GNLY)
ranked well below.

**The second pass found the same thing three more times, by three more
methods.**

- **The backbone tree falls into two pieces.** The 21,545 mature lymphoid cells
  (56% of the data) are joined to the rest of the tree by exactly one edge —
  plasmacytoid DC to T activated, weak (0.137) and the longest edge on the map
  (11.3 units against a median of 3.7). Remove it and the lymphoid block is
  disconnected. The tree had to invent that edge to stay in one piece.
- **The principal tree also had to invent a bridge, in a different place.** It
  attaches the lymphoid block through the platelet branch, along three segments
  that hold almost no cells. Two unrelated methods both needed a made-up
  connection and put it in different places. If a real connection existed, both
  would have found the same one.
- **Velocity gives the lymphoid block its strongest arrows, and they all run
  backwards** — NK → CD8 cytotoxic → CD8 memory → CD4 memory → CD4 naive, and
  memory → naive is the reverse of T cell biology. They also point the opposite
  way to the CellRank arrow field in the same island. Two methods, two coherent
  and opposite directions: the island has no direction to find.

So the negative result has now been shown six independent ways: a dense block in
PAGA, a flat stretch in the pseudotime, a single undifferentiated CellRank basin,
a two-piece backbone tree, a principal tree bridge in the wrong place, and
confident backwards velocity.

**The main lesson of the stage: a trajectory method always returns an ordering.
It has no way to say "there is no process here".** Everything the methods
produced for the mature compartment is internally consistent, reproducible, and
biologically meaningless. Only the labels revealed it, and the labels were not
part of the calculation.

**Three lineages have no starting cell on their axis.** Erythroid early,
erythroid late and plasmacytoid DC are 93–100% one label in every one of their
40 bins. There is no precursor at the beginning. The gene order inside them is
still useful, but the axis is maturation within a cell type, not differentiation
into it.

**One terminal state is not a destination at all — it is the cell cycle.** The
"Pre-B large" terminal's driver genes are almost entirely cell-cycle genes
(UHRF1, STMN1, TYMS, NUSAP1, MKI67, TOP2A, PCNA, AURKB), with only IGLL1 and
VPREB1 belonging to B cells. Its pseudotime confirms it: G1/S genes peak at bin
32, S/G2 at 34–35, mitosis at 39. Dividing cells look like each other, which
makes them a stable state, which makes them look like a destination. Two
separate results, same conclusion.

The same confound appears elsewhere: T proliferating has the second-earliest
pseudotime in the whole dataset (0.404), because cycling cells resemble the
cycling progenitors near the root. In velocity, the same 47 cells send strong
arrows into HSPC and Pre-B large.

**The erythroid terminal split in two.** Erythroid early and erythroid late each
formed their own terminal, and they exchange 0.11 and 0.13 of their probability.
The driver heatmap shows both columns carrying the same genes. One biological
process, two states.

**Velocity gets the erythroid direction backwards, and nothing fixes it.** In all
four velocity fits the arrows run erythroid late → erythroid early, and in some
on into the stem cells. The reasons are in the counts:

- erythroid late has an unspliced fraction of 0.045, at or below the level
  platelets give (0.050) — and platelets have no nucleus, so their "unspliced"
  signal is pure background. There is almost no new RNA to read;
- the erythroid programme genes do not fit the kinetic model. In the
  whole-sample run they were not even selected; in the progenitor run they were
  given to the model and scored essentially zero (KLF1, GYPA, ALAS2, HBA1 at
  about 10⁻¹³ to 10⁻¹⁹). The direction is set by whatever genes do fit, and
  those do not track maturation.

The progenitor run was built to test the idea that the lymphoid cells were
drowning out the erythroid signal. **That idea was wrong.** With the lymphoid
cells removed, four times as many genes fitted and the erythroid genes inside
the model, the backwards arrow got *stronger* — 0.12 → 0.86 in one model,
0.32 → 0.82 in the other. The fault is in the erythroid cells themselves. This
is a known failure of RNA velocity in red-cell development (Barile et al. 2021,
*Genome Biology*; cited from memory, not re-checked against this data). The
erythroid order in this stage rests on Palantir and the gene trends, not on
velocity.

**Velocity points into the stem cells, not out of them.** In the progenitor run
four lineages send arrows into HSPC — platelet (0.271, the second strongest
arrow in its model), Pre-B large, erythroid late and granulocyte precursor —
and HSPC sends only one arrow out. The 103 platelets, which make no new RNA,
produce one of the strongest arrows. In velocity's view the progenitors look
like a destination.

**Velocity's own orderings start in the wrong places.** Left to choose its own
starting point, the steady-state model started in the erythroid island and the
kinetic model in the myeloid compartment. Neither started at HSPC. Across all
cells both run against Palantir (−0.44 to −0.56).

**A rooted velocity pseudotime measures distance, not direction.** When the
progenitor run was forced to start at HSPC, it ordered erythroid early before
erythroid late — correctly — while its own arrows said the reverse. Both are
true readings of the same fit: erythroid late is simply further from HSPC on the
graph, whichever way the arrows point. A good rooted ordering is therefore no
evidence that velocity got a direction right.

**Velocity confidence is highest where velocity is most wrong.** Erythroid early
has confidence 0.992 (whole sample) and 0.999 (progenitor) while its arrows run
backwards. HSPC is at 0.999 while arrows run into it.

**The kinetic model fits few genes, and not the useful ones.** Whole sample:
median likelihood 0.052, nine genes above 0.3, the best ones lymphoid (GNLY,
NKG7, CD74, IGHM). Progenitor run: better (845 genes above 0.1 against 185), but
the top gene is MALAT1, a nuclear RNA that fits for technical reasons, and a
third of the top eighteen are housekeeping or ribosomal genes. The
differential-kinetics test flagged two-thirds of all genes, which makes it
uninformative; it even flagged haemoglobin in monocytes, which points to
free-floating red-cell RNA in the sample rather than a real kinetic difference.

**Velocity arrows change a lot between runs.** Fewer than half of the progenitor
run's arrows also appear in the whole-sample run. Only the monocyte steps, the
early B step and the erythroid reversal are stable. HSPC → Pre-B large was right
in one run and reversed in the other.

**The principal tree's global ordering is reversed against Palantir** (−0.42).
The cause is specific: to keep the tree to about ten tips it needs a length
penalty so high that it prefers straight cuts through empty space to curves
through cells. 18% of the tree's length runs through regions with no cells; its
trunk passes through the pre-B cells, so every other lineage inherits that
distance; and it enters the myeloid compartment from the monocyte end, so the
granulocyte precursor — the earliest committed myeloid cell — comes out as the
*latest* cell type in the whole tree. Thirty-one HSPC cells that sit right next
to the root on the map end up at the far end of the tree. Its ordering is not
used; only its topology checks are kept.

**The tree's strongest-looking edges can be the least trustworthy.** T
proliferating – CD4 T memory is the second strongest connection in the dataset
(0.914) yet spans the whole map (10.9 units). Only 47 cells — a handful of
neighbours saturates the score. Strong connection plus long distance marks a
small-population artefact.

**Palantir's entropy puts two finished cell types above the stem cells.**
Plasmacytoid DC (1.088) and Pre-B large (1.002) sit above HSPC (0.999). These
are exactly the two cell types Palantir has no branch for, so their probability
has nowhere to settle. CellRank, which gives plasmacytoid DC its own terminal,
calls it 0.848 committed. Palantir entropy only means something for cell types
that have a branch.

**There is a real lane effect in the counts — but it cancels before velocity
sees it.** Inside the same cell type, lanes MantonBM1 and MantonBM3 carry 12–13%
more unspliced RNA, and MantonBM4 and MantonBM5 about 5–7% less — about 20% top
to bottom, against a 2.5-fold spread between cell types. The cause was not
found. scVelo divides each cell's unspliced counts by that cell's own unspliced
total, so a lane-wide offset cancels out; a planned "lane-corrected" re-run was
therefore dropped as redundant. Velocity *confidence* still varies by lane:
MantonBM3 lowers it across the T cells (down to −0.28 in CD8 T memory), and its
B naive cells swing from +0.145 in one fit to −0.258 in the other.

**The platelet terminal is numerically fragile.** It keeps only 0.538 of its own
cells and loses 0.270 to the lymphoid basin. It also receives almost nothing
from the stem cells (0.008) — there is no intermediate population connecting
HSPC to the platelets, so no path reaches it. There are only 103 platelets, and
only 49 of them pass the commitment threshold, which is too few to test gene
trends at all. Its drivers are the most specific in the dataset and its method
agreement is the cleanest (0.97, nothing else above 0.04). Strong signal,
fragile terminal, same cause: tiny population.

**Palantir found only four terminals; CellRank found eight.** Palantir has no
branch for erythroid, Pre-B large or plasmacytoid DC — the correlations in those
four columns are all below 0.16. For erythroid this is a Palantir limitation,
since the gene trends independently support that terminal. Palantir's terminal
detection relies on the diffusion eigenvalue spectrum, and that spectrum was
already known to be flat here (1.000, 0.9992, 0.9988, 0.9983, 0.9924 — almost no
gap between components), so it had little to work with.

**Some pseudotime axes are too compressed to read.** The whole B naive lineage
spans 0.83 to 0.866. Twenty bins of difference can come from a small shift.
Label transitions and gene rank order survive this; fine position does not.

**T activated is a one-donor cell type.** 532 of its 583 cells come from
MantonBM3. Any result for that label describes one donor, and part of its high
unspliced fraction is that lane's offset.

---

## 7. Things that looked wrong and were checked

Every check took seconds to minutes, and every one changed the write-up.

**A. Pre-B small cells reappearing at the end of the B naive axis**

172 pre-B cells sat in the last two bins, after the mature B cells. Three
possible causes were tested at once.

The fate probability turned out to be fine — all pre-B cells are ~0.94 committed
to the B terminal, which is correct biology. They are not doublets;
`predicted_doublet` is False for all 172. What separates them is depth: they
carry 574 genes and 1274 counts, against 660 and 1476 for pre-B cells as a
whole. **They are the thinnest-sequenced cells in the lineage.** A cell with few
counts has few distinguishing transcripts, so its neighbours are decided by
noise and it drifts to the far end of the diffusion distance.

This is left uncorrected on purpose. Palantir did not run wrong; the input cells
are thin, and no parameter adds counts to a cell that lacks them. The real fix
is a stricter quality-control threshold in the core pipeline, which would mean
re-running everything downstream to repair 0.45% of the dataset. No conclusion in
this stage depends on those cells.

**B. Three driver genes that made no biological sense**

| gene | why it looked wrong | what the raw counts showed |
|---|---|---|
| EPCAM in erythroid | an epithelial gene in a red-cell precursor | **real** — 89% of erythroid early cells, mean 6.5 counts, and zero in NK, platelets and pre-B small. Contamination would be spread evenly; this is not |
| PTCRA in platelets | a known plasmacytoid DC gene | **real in both** — platelets 47.6% of cells at 4.0 counts, pDC 47.3% at 1.9 counts. Platelets carry *more* per cell, so it is not spillover from pDC |
| three histone genes in platelets | platelets have no nucleus and do not divide | **real, but not a cell-cycle signal** — the dividing populations barely express them (~1 count), while platelets reach 10.5 counts for HIST1H2AC. Residual mRNA inherited from the megakaryocyte |

**Two of the three explanations guessed before the check were wrong.** Neither
correlation rank nor biological plausibility caught it. The raw counts did, in
five seconds. The rule taken from this: read `pct_expressing` and mean count
together — a gene in 41% of cells at 1.9 counts and a gene in 59% of cells at
10.5 counts look similar on percentage alone and are not the same thing.

**C. The velocity setup, checked before any model was fitted**

- **The right reference.** The first BAM's header showed the reads were aligned
  with STAR against GRCh38 using GENCODE v27, so that exact annotation was
  downloaded. A mismatched annotation gives counts that look right and are not.
- **Skipped reads.** An early estimate said about 4% of reads would be skipped
  for a missing barcode. The full logs showed 22–29% in every lane. **The first
  estimate was wrong**, but harmlessly: the extra reads come from empty droplets
  outside the kept cells, so no kept cell lost any depth.
- **Mitochondrial genes were not counted**, because the reads name the
  mitochondrial chromosome `MT` and the annotation names it `M`. They carry no
  velocity information and are filtered out anyway.
- **The lane effect, tested twice.** The first test, on lane totals, blamed the
  low lanes on their cell-type mix. The second, comparing the same cell type
  across lanes, showed a real lane offset on top. **The first reading was
  withdrawn.**

**D. Three problems inside the velocity run**

- **NumPy 2.4 and newer breaks scVelo's kinetic model.** Pinned below 2.4.
- **The velocity-directed PAGA was quietly using a pseudotime prior**, which
  would have made its arrows partly a restatement of a pseudotime. It was turned
  off. **The arrows came out identical to every decimal**, so the prior had no
  effect on this data; the fix stays in so the arrows can never depend on it.
- **scVelo stores its transition matrix the "wrong" way round** — row is the
  target, column the source. Read the other way, every arrow in the stage would
  be reversed. The reading was checked against the streamline figures, where
  the flow is visible directly.

---

## 8. What to trust, and what not to

**Trust**

- The four terminals confirmed by both methods: monocyte classical, B naive, NK
  (read as "mature lymphoid"), platelet.
- The monocyte lineage, completely — ordering, composition, gene sequence, and
  the velocity direction.
- The early B step, Pre-B large → Pre-B small, confirmed by velocity in all four
  fits.
- Stem cell multipotency and the CellRank entropy ranking.
- The erythroid gene order, including transcription factors before targets.
- The erythroid stage order (early before late) from Palantir and the principal
  tree — but not from velocity.
- The five textbook driver lists.
- The Pre-B large terminal being the cell cycle.
- The mature lymphoid negative result, now shown six independent ways.
- The velocity counts themselves: complete, and positions not driven by lane.

**Do not trust**

- Any position along the mature lymphoid axis. It is a cytotoxic gradient, not
  development.
- Any velocity direction in the erythroid lineage or into HSPC.
- Any velocity direction in the mature lymphoid block.
- Velocity confidence as evidence that a direction is right.
- A rooted velocity pseudotime as evidence of direction — it behaves like
  distance.
- Any velocity arrow that appears in only one of the two runs.
- The principal tree's global ordering (its topology checks are fine).
- The CellRank arrow field as a direction — it follows our chosen root and is
  not velocity.
- Palantir entropy for cell types with no Palantir branch (plasmacytoid DC,
  Pre-B large).
- Fine position within B naive, monocyte classical, or the mature T populations
  — their axes are compressed.
- The last two bins of the B naive gene trend figure.
- Platelet pseudotime — 103 cells spread across almost the entire axis — and
  any platelet velocity arrow.
- Any T activated result as a statement about the cell type rather than one
  donor.
- Any single driver gene that has not been checked against raw counts.

---

## 9. Files produced

**Scripts** — `src/trajectory/`: `trajectory_utils.py`, `paga_topology.py`,
`pt_palantir.py`, `pt_cellrank.py`, `fate_drivers.py`, `gene_trends.py`,
`lineage_diagnostics.py`, `trajectory_report.py`, `trajectory_backbone.py`,
`principal_graph.py`, `velocity_barcodes.py`, `velocity_counts.py`,
`velocity_scvelo.py`, `velocity_progenitor.py`

**Raw and intermediate data**

- `data/raw/velocity_bam/` — the eight lane BAMs (checksum-verified) and the
  velocyto outputs. The cell-sorted BAM copies velocyto makes can be deleted.
- `data/external/velocity_reference/` — the GENCODE v27 annotation and the hg38
  repeat mask.
- `data/interim/bm_velocity.h5ad` — the merged spliced and unspliced counts.
- `data/interim/bm_velocity_fit.h5ad` and `bm_velocity_progenitor_fit.h5ad` —
  the saved velocity fits, reused on re-runs. Delete one to refit.

**Tables** — `results/trajectory/`:

- first pass: PAGA connectivities, edges and trees for both embeddings; Palantir
  pseudotime, roots, entropy, terminal states and run correlations; CellRank
  macrostates, terminal states and fate probabilities at k = 8, 5 and 4; driver
  correlations and top-50 lists; binned gene trends, composition and peak order;
  the two diagnostic tables; and the three stage-summary tables;
- second pass: backbone centroids and tree edges; principal tree pseudotime,
  nodes, segments, per-cell-type summary and agreement; the eight barcode lists
  and their summary; velocity recovery and per-lane count tables; for each
  velocity run the per-cell table, agreement, lane check, PAGA arrows and gene
  table; and, for the progenitor run, the arrow comparison and the direction
  check.

**Figures** — `figures/trajectory/`:

- first pass: PAGA networks and heatmaps, the Palantir UMAP grid, CellRank fate
  UMAPs, the driver heatmap and UMAP panels, seven per-lineage gene trend
  figures, and four stage-summary figures;
- second pass: five backbone figures (tree by cell type, by pseudotime, as
  nodes, entropy, arrow field), five principal tree figures, the lane
  composition figure, the unspliced-fraction-by-lane figure, and eleven figures
  per velocity run (two streamline plots, two PAGA arrow plots, the four
  orderings on the map, orderings by cell type, and five lane heatmaps).

**Full analysis** — `claude/report_trajectory.md`, organised figure by figure
with every number and the reasoning behind it.

---

## 10. What is still open

None of these is needed to close the stage.

- Removing cell-cycle genes before computing fates, to test directly whether the
  Pre-B large terminal dissolves.
- Two principal tree follow-ups: fitting on the progenitor side only, or
  starting the tree from the PAGA backbone. The second is the one aimed at the
  actual cause of its problems.
- A separate velocity fit on the erythroid lineage alone. Not recommended: the
  failure is explained, the erythroid order is already set by two other
  methods, and it would change no conclusion.
- The cause of the lane offset in unspliced RNA.
- The identity of one small cluster left of centre that is the brightest cells
  in every ordering — Palantir and all six velocity panels.
- One planned figure was dropped: the backbone tree coloured by velocity
  pseudotime. The whole-sample velocity pseudotime picks its own start, and none
  of its starts is HSPC, so the figure would mislead. The velocity-directed PAGA
  figures replace it.
- An R replay of the stage, if the optional final phase is done.

Everything else raised during the stage has been closed.
