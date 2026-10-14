# Pathway analysis of human bone marrow single-cell data

A plain-English record of what was asked, what was run, what came out, and what
should not be believed.

Dataset: eight bone marrow donors from the Human Cell Atlas (Hay 2018 Census of
Immune Cells), 38,430 cells and 20,455 genes after quality control and doublet
removal, labelled into 19 cell types. The labels were fixed before this work
started and were not changed by it.

---

## 1. The question

A cell type is a list of genes. That list on its own says nothing about what the
cells are *doing*. Pathway analysis is the step that turns gene lists into
biology, and there is more than one way to do it. The three main ways answer
different questions, so all three were run and compared.

**Over-representation** takes a short list — the top marker genes of a cell type
— and asks whether that list overlaps a known gene set more than chance would
allow. It is permissive and fast, and it only ever looks at the top of a list.

**Gene set enrichment** takes the *whole* ranking of every gene, top to bottom,
and asks whether the members of a gene set cluster at one end. It is stricter,
it sees the bottom of the list as well as the top, and it needs a coordinated
shift rather than a handful of overlapping genes.

**Per-cell activity** does not use lists at all. It scores every individual cell
against a set of genes known to change when a pathway is switched on, giving a
map rather than a table.

A fourth, independent view came from gene programs found by unsupervised matrix
factorisation in earlier work. Those programs were built without any knowledge
of the cell type labels, so testing them is a check on both the programs and the
labels.

**The point of running all four was the disagreements.** Where two methods agree,
the result is solid. Where they disagree, the disagreement itself is the finding.

---

## 2. Reference data used

Three catalogues of gene sets, all downloaded once and cached to disk so the
analysis reproduces offline and cannot shift under itself:

| catalogue | what it is | sets |
|---|---|---|
| GO Biological Process | a layered vocabulary of biological jobs, from very broad to very narrow. Terms are nested and overlap heavily. | 5,406 |
| KEGG | hand-drawn maps of metabolic and signalling routes. Small and curated. | 320 |
| Reactome | step-by-step reactions, more detailed than KEGG and updated more often, organised in a hierarchy. | 2,100 |

Why all three: they disagree. A result found in only one catalogue is usually an
artefact of how that catalogue was built. A result found in all three is real.

For per-cell work: **PROGENy**, 14 signalling pathways, each defined by the genes
that change when the pathway fires — so it measures the *consequence* of
signalling, not the expression of the signalling molecules themselves. And
**DoRothEA**, transcription factor regulons at confidence levels A to C only;
levels D and E are largely computational guesses and add noise.

Fixed settings throughout: random seed 7, significance at an adjusted p-value
below 0.05, top 200 genes per cell type for the list-based test.

---

## 3. Differential expression — the shared input

Every gene was tested in every cell type against all the other cells, one cell
type at a time. The test is the Wilcoxon rank-sum test: it compares rank order
rather than averages, so a few extreme values cannot drag the result and the data
does not need to follow a bell curve. Single-cell counts never do.

**Ribosomal, mitochondrial and haemoglobin genes were removed before testing.**
These are high in every cell and reflect how busy or how oxygen-rich a cell is
rather than what it is. Left in, they fill the top of every list. 20,455 genes
became 20,328.

Two tables came out of this, because the two list-based methods need different
inputs:

- the full ranking, every gene in every cell type — 386,232 rows, which is what
  enrichment needs, since it looks at the whole order rather than a cut-off
- the top 200 significantly raised genes per cell type — 3,800 rows, which is
  what over-representation needs

**All 19 cell types hit the 200-gene cap**, meaning every label has at least 200
significantly raised genes. The cap is binding, not a floor.

Runtime: about one minute.

---

## 4. Over-representation of the marker lists

Run against the cached catalogues rather than a web service, so results do not
change between runs. 19 cell types × 3 catalogues = 57 tests, 17 seconds.

**The comparison set was the 20,328 genes that actually entered the test**, not
the whole genome. A gene that was never measured cannot be a hit, so including it
would inflate every p-value.

Result: 8,362 significant terms.

### The labels validate

Every one of the 19 labels lands on the expected biology:

| cell type | top terms |
|---|---|
| B naive | B cell receptor signalling, immunoglobulin-mediated immune response |
| Pre-B small | B cell receptor signalling, antigen activates the B cell receptor |
| CD4 and CD8 T, all four states | T cell activation, T cell receptor signalling, immunoregulatory lymphoid interactions |
| NK | leukocyte-mediated cytotoxicity, natural killer cytotoxicity |
| Platelet | platelet aggregation, platelet activation, haemostasis |
| Dendritic cell and both monocytes | MHC class II antigen presentation, phagosome, neutrophil degranulation, innate immune system |
| Granulocyte precursor | neutrophil degranulation, innate immune system |
| Erythroid early | heme biosynthesis, porphyrin biosynthesis |
| Pre-B large, T proliferating, both erythroid stages | cell cycle, mitotic cell cycle, messenger RNA splicing — all cycling populations |

This is independent support for the annotation. The labels were made from marker
genes; the catalogues agree without ever being told what the labels are.

### A gap in the gene filter, found and fixed

The first run put **"L13a-mediated translational silencing"** at the top for
naive CD4 T cells. Reading the genes behind the term showed translation
initiation factors plus **FAU** and **UBA52**. Both are ribosomal proteins whose
names do not begin with the usual prefixes, so the filter had missed them.
**RACK1** is a third case.

The filter was widened to drop all three, and both steps were re-run. Dropped
genes went from 124 to 127. What changed:

- no output row anywhere still contains those three genes
- the misleading term fell from first place to second, its overlap dropped from
  10 of 26 genes to 8 of 24, and the remaining genes are pure translation
  factors. So naive T cells do carry a real translation-initiation signal — it
  was inflated, not invented
- B naive, Pre-B small and stem/progenitor cells gained cleaner top terms
- nothing else moved, which is the correct outcome: the fix touched only what it
  should have

### Three results that look striking and are not

These come from how the catalogues are built, not from this data. They must be
recognised and never reported as biology.

**KEGG neurodegeneration maps — Parkinson, prion, Huntington, motor neurone
disease.** Hit by granulocyte precursors, stem cells, late erythroid and
proliferating T cells. The genes behind them are respiratory chain subunits and
proteasome subunits. Those maps are built mostly from the respiratory chain, so
any cell with high energy metabolism scores. Read as high oxidative metabolism,
nothing more.

**"Cellular respiration" and "proton motive force-driven ATP synthesis" in
granulocyte precursors.** Same cause. The real answer for that population is the
Reactome hit: neutrophil degranulation.

**"Telomerase RNA localization to Cajal body" in late erythroid and stem cells.**
The genes are CCT3, CCT5, CCT6A, CCT8 and TCP1 — a protein-folding machine that
the term happens to contain. Reactome states the same signal plainly as "folding
of actin by the CCT chaperonin".

**The general lesson: read the gene list before believing the term name.** Two of
these three look like remarkable findings and are nothing of the kind.

### Significant terms per cell type

| cell type | GO | KEGG | Reactome |
|---|---|---|---|
| B naive | 98 | 39 | 38 |
| CD4 T memory | 129 | 74 | 138 |
| CD4 T naive | 38 | 14 | 55 |
| CD8 T cytotoxic | 212 | 54 | 239 |
| CD8 T memory | 200 | 67 | 239 |
| Dendritic cell | 285 | 75 | 171 |
| Erythroid early | 150 | 24 | 200 |
| Erythroid late | 187 | 37 | 384 |
| Granulocyte precursor | 178 | 29 | 244 |
| HSPC | 133 | 36 | 223 |
| Monocyte classical | 418 | 85 | 156 |
| Monocyte non-classical | 466 | 84 | 183 |
| NK | 204 | 67 | 138 |
| Plasmacytoid DC | 176 | 58 | 105 |
| Platelet | 173 | 49 | 208 |
| Pre-B large | 194 | 36 | 414 |
| Pre-B small | 125 | 48 | 74 |
| T activated | 136 | 24 | 144 |
| T proliferating | 93 | 59 | 455 |

The two monocyte populations carry the most GO terms; the cycling populations
carry the most Reactome terms. Naive CD4 T cells have the fewest of any label,
which fits a resting cell with little distinctive activity.

---

## 5. Gene set enrichment on the full ranking

The strict test. A thousand random shuffles per gene set, sets limited to between
15 and 500 genes — very small sets give unstable scores, very large ones are too
vague to mean anything.

**This was the slowest step by a wide margin: 45 minutes.** GO alone took 27
minutes. Results were saved one catalogue at a time so a crash would cost only
the catalogue in progress.

Result: 2,837 significant terms, far fewer than over-representation's 8,362.
That gap is expected and is itself informative — the question is stricter.

### Strong and correct results

- **B naive and Pre-B small** — CD22-mediated B cell receptor regulation,
  classical antibody-mediated complement activation, among the strongest scores
  in the table
- **CD8 T and NK** — immunoregulatory interactions between lymphoid and
  non-lymphoid cells, lymphocyte migration, chemokine receptor binding
- **Platelet** — platelet degranulation, response to raised internal calcium
- **Monocyte classical** — neutrophil degranulation, defence against bacteria,
  breakdown of the extracellular matrix
- **Erythroid late** — mitochondrial translation, with 72 of the 97 genes in the
  set contributing. The best-supported positive result anywhere in this work

### What this method adds that the other could not: the down side

**1,169 of the 2,837 significant results are negative**, meaning the gene set
sits at the *bottom* of the ranking rather than the top. Over-representation only
ever sees the top of a list, so this half of the signal was invisible until now.
It reads correctly:

- **Early erythroid cells score most strongly negative of anything in the table
  for T cell receptor signalling.** A red cell precursor is about as far from a T
  cell as this dataset goes.
- **Platelets and monocytes turn down messenger RNA splicing and translation.**
  Platelets have no nucleus, so this is exactly as it should be.
- **Every non-B cell type turns down B cell receptor regulation.**

### Weak results, and the rule that follows

Some significant terms are carried by only one or two genes:

| cell type | term | contributing genes |
|---|---|---|
| Plasmacytoid DC | skin epidermis development | 2 of 17 |
| HSPC | establishment of skin barrier | 2 of 15 |
| Erythroid early | regulation of cardiac muscle repolarisation | 1 of 17 |

Set against the erythroid mitochondrial result at 72 of 97, the difference is
obvious.

**Rule: check how many genes actually drive a term before believing it.** A good
significance value on its own is not a finding. Anything under five contributing
genes was excluded from the comparison below.

### A known limitation: tied ranks

Between 20% and 49% of genes share an identical test score, so their order within
the ranking is arbitrary. These are genes with no real difference between the
cell type and the rest, and they sit in the middle and at the bottom. The terms
at the two extremes — the ones that get reported — are unaffected. Recorded so
the warning is not rediscovered later and mistaken for a fault.

### Two notes worth carrying

**Mitochondrial ribosomal proteins are deliberately not filtered.** They do not
follow the naming pattern of the cytoplasmic ribosomal genes that were removed,
and in erythroid cells the mitochondrial machinery is real biology — heme
synthesis happens there. Removing them would delete the finding.

**The energy-metabolism axis appears again** in dendritic cells, granulocyte
precursors and proliferating T cells — the same axis that produced the false
neurodegeneration hits in the previous step. Two independent methods agreeing
that these populations are metabolically busy is a real observation. The disease
names attached to it are not.

---

## 6. Per-cell signalling pathway activity

14 pathways scored on every one of the 38,430 cells, using all 20,455 genes
rather than the variable ones — the footprint genes are spread across the genome
and subsetting would discard most of them. 4,638 genes carried a footprint
annotation. Runtime: one minute.

### How to read the output: down a column, never across a row

Each pathway has its own scale and offset. p53 scoring negative in every cell
type means nothing on its own. **Only the differences between cell types within
one pathway carry information.**

### p53 tracks cell division — the cleanest result here

Most suppressed in early erythroid, large pre-B, stem/progenitor and
proliferating T cells: precisely the four dividing populations. Least suppressed
in late erythroid and platelets, the two that have stopped dividing. Low p53
activity permits division, so the ordering is internally consistent — and nothing
told the method which cells cycle.

### The inflammatory axis sits where it should

TNF alpha is highest in proliferating T cells, both monocyte populations, memory
CD4 T cells and dendritic cells. JAK-STAT is highest in proliferating T cells,
non-classical monocytes, dendritic cells and NK cells — the cytokine and
interferon responders.

### Hypoxia looks backwards, and is not

Highest in the two monocyte populations, **lowest in early erythroid cells**. Red
cell production is driven by low oxygen, so the naive expectation is the
opposite. But the hypoxia footprint is mostly glycolysis genes, and erythroid
precursors run on mitochondrial metabolism for heme synthesis — which the
enrichment step independently found as its strongest erythroid result. **Two
methods, same conclusion, arrived at separately.**

### A built-in negative control

**Androgen and oestrogen signalling are flat across every cell type.** Bone
marrow has no sex hormone signalling, and the method correctly finds nothing.
WNT is similarly flat. That flatness is evidence the method is not inventing
signal, and it later mattered for how the figures are drawn.

Other result: TRAIL is highest in B naive and small pre-B cells, lowest in early
erythroid.

---

## 7. Per-cell transcription factor activity

Hundreds of factors, each judged on its own set of target genes. The network held
429 factors; 293 had enough targets present in the data to be scored, across
8,059 genes. Runtime: 20 seconds.

A different statistical model was used here than for the signalling pathways. The
14 pathways compete for overlapping target genes, so a joint fit is correct
there. Factors are judged one at a time against their own targets, so a separate
fit per factor is correct here.

### The raw ranking is misleading, and this is the key methodological point

Taking the highest average score per cell type returns **MYC first in all 19 cell
types**, followed by the same handful everywhere. That is not biology.

**Cause: scores are not comparable between factors.** A factor with a large set
of broadly expressed targets scores high in every single cell. The absolute value
says nothing about which factor matters where.

**Fix: compare each factor against itself across the cell types** before ranking
— the same rule as the pathway table, stated there as "down a column, never
across a row", and it matters more here because there are 293 columns instead of
14.

### After that correction, the factors are textbook

| cell type | top factors | reading |
|---|---|---|
| Granulocyte precursor | CEBPE, CEBPA, MYB | CEBPE is the defining early granulocyte factor |
| Monocyte classical | CEBPB, IRF8, ETS2 | canonical monocyte set |
| Monocyte non-classical | NR1H3, PPARG, FOXO3 | the lipid-sensing factors of patrolling monocytes |
| Dendritic cell | RFX5, VDR, FOSL1 | RFX5 controls MHC class II, which was the dendritic hit in both list-based methods |
| Plasmacytoid DC | IRF9, IRF3, RARG | the interferon-producing cell, found by its interferon factors |
| Erythroid late | GATA1, KLF1 | the two master red-cell factors, ranked first and second |
| Erythroid early | KLF1, ZHX1 | same lineage, earlier stage |
| Platelet | NFE2L1, KLF5 | the NFE2 family is the master megakaryocyte factor |
| Pre-B large | E2F5, E2F2, MYBL2 | pure cell cycle, matching the cycling label |
| B naive | PAX5, SPIB, IRF4 | textbook B cell programme |
| CD4 T naive | TCF7 | the defining naive T factor |
| CD8 T memory and cytotoxic | SATB1, ATF3, REL, STAT4 | correct for cytotoxic states |
| NK | STAT5A, STAT5B, STAT4 | interleukin-2 and interleukin-15 signalling, exactly right |
| T proliferating | ETV6, STAT2, IRF3 | T identity plus interferon factors |

This is the strongest independent confirmation of the annotation in the whole
analysis. The regulons come from a database built on entirely different
experiments, and nothing in the method was told which cells are which.

### Three weak rows

- **T activated** returns no recognisable factor
- **Stem and progenitor cells** return general stress and metabolic factors
  rather than the expected stem-cell factors
- **Small pre-B cells** return what reads as noise, even though their identity
  was clear in the list-based test

One suspicious result: an eye-development factor scores highest in platelets.
That is almost certainly overlap between its target set and megakaryocyte genes,
not real activity.

---

## 8. The unsupervised programs

The 30 gene programs found earlier by matrix factorisation were tested the same
way. They were built with no knowledge of the labels, so this asks whether
unsupervised programs and curated biology agree. Runtime: 10 seconds.

The comparison set was the 2,000 genes the factorisation was fitted on, not the
20,328 — a program could only ever contain those genes.

The program gene lists were deliberately left unfiltered. One program being
ribosome-heavy is itself a finding, and filtering would hide it.

**27 of the 30 programs have at least one significant term.**

### The activity programs are confirmed

This matters most, because separating activity from identity was the unique
contribution of the factorisation.

| program | what it was called | what the catalogues say | contributing genes |
|---|---|---|---|
| DNA replication | cell cycle | 58 of 163 | strongest result in this step |
| mitosis | mitotic cell cycle | 52 of 139 | |
| heat shock | "attenuation phase" | 7 of 11 | this is the heat-shock pathway under its Reactome name |
| inflammation | neutrophil degranulation | 26 of 109 | |

### The identity programs land correctly

Monocyte and granulocyte programs return neutrophil degranulation; the dendritic
program returns MHC class II presentation; the erythroid programs return
erythrocyte gas exchange, one of them at 7 genes out of 7 — the cleanest identity
match anywhere; the platelet programs return platelet activation; the early B
program returns B cell activation.

### Two of three newly found cell types are confirmed

The factorisation had suggested three cell types missing from the 19 labels:

- **plasma cells — confirmed**, through an immunoglobulin-driven set at 14 of 19
  genes
- **macrophages — confirmed**, through the innate immune system
- **mast cells and basophils — not confirmed.** No significant term in any
  catalogue. The program remains supported only by its marker genes, which are
  convincing on their own, but this step adds nothing to them

### Three programs return nothing, for the same reason

The stem/progenitor program, the plasmacytoid dendritic program and the mast cell
program have no significant term anywhere.

**These are the three programs defined by rare cell *identity* rather than by a
*process*.** Catalogues describe processes — division, degranulation, gas
exchange — and have no set for "stem cell" or "mast cell". The failure is in the
catalogues, not in the programs.

The stem and progenitor compartment has now failed twice: no distinctive
transcription factor, and no enriched gene set. Its support rests on marker genes
instead.

### One correction and one non-confirmation

The program previously described as "T cell mixed with housekeeping" returns
ribosomal RNA processing as its top term, with 63 significant Reactome terms and
a single GO term. **It is mostly the translation machinery, not a T cell
programme**, and should be described that way.

The histone programme returns an erythroid term at 4 of 7 genes with only three
significant terms in total. Its gene list mixes histones with erythroid genes, so
the hit follows the minority. This programme stays supported by its gene names
alone.

### The unexpected stromal programs hold up

Four programs had been flagged as stromal genes that this kind of preparation
should barely contain. They return coherent stromal biology — growth factor
transport, and secreted protein processing. That supports a genuine trace of
stromal cells rather than noise in the factorisation.

---

## 9. The figures

### The per-cell pathway maps

Four regions of the map are identifiable from the pathway values alone: the
erythroid island (deep blue hypoxia, deepest blue p53), the monocyte blob (red
hypoxia, red TNF alpha), the B lineage (red TRAIL) and the large T cell region.
This matches the per-cell-type table exactly, so map and table tell one story.

**Main finding from the maps: activity follows identity, not sub-structure.**
Nearly every panel is flat within a region and changes only between regions. On
this dataset, signalling activity is a property of cell type rather than a state
cutting across types. That is a real negative result, and the maps were drawn
precisely to test for it.

The clear exception is the **erythroid island splitting in two on hypoxia and
p53**, matching the early and late boundary that clustering already found.

**A fault in the first version, and its fix.** The colour scales were set
per-panel and driven by a handful of extreme cells — one pathway ran to 22 while
another spanned about ±3 — which flattened the middle of most panels and left
several looking blank. Clipping each scale to the 1st and 99th percentile fixed
it. Structure is now visible in four panels that previously showed almost
nothing.

### The pathway heatmap

The strong blocks are right: early erythroid extreme on hypoxia and TRAIL, late
erythroid and platelets high on p53, activated T cells highest on NF-kB,
monocytes and dendritic cells warm across growth factor and inflammatory
pathways.

**A trap that had to be designed around.** Standardising each pathway so cell
types can be compared also makes a pathway with *no* signal look as loud as a
real one. The three negative controls — androgen, oestrogen and WNT — showed
full-saturation colour despite spanning almost no range in the raw values. The
fix drops any pathway whose raw range across cell types is negligible. The
controls measure about 1.0 to 1.3; every real pathway measures 1.7 or more, so
they separate cleanly. The published heatmap carries 11 pathways.

**Not a contradiction:** those three pathways still show visible structure on the
per-cell maps while being dropped from the heatmap. The heatmap compares cell
type *averages*, which barely differ; the map shows *single cell* spread, which
is wider. Two different questions.

### The transcription factor heatmap — the strongest figure

A clean diagonal band: each cell type owns its own block of factors, and the
blocks follow the lineage order down the rows. Monocyte factors, dendritic
factors, erythroid factors, cycling B cell factors, B cell factors, T cell
factors, NK factors, each in its own place.

Three rows break the pattern, and they are the same three that failed the earlier
checks: stem and progenitor cells are pale across the whole row, activated T
cells have only a weak and unrecognisable patch, and small pre-B cells rest on
factors that read as noise.

### Agreement between the two list-based methods

The most informative figure. Over-representation finds far more terms than
enrichment in every cell type, but the structure is in **where the two agree at
all**:

| high agreement | terms in both | low agreement | terms in both |
|---|---|---|---|
| Erythroid late | 196 | CD4 T naive | 0 |
| Pre-B large | 189 | T activated | 1 |
| T proliferating | 179 | CD4 T memory | 5 |
| HSPC | 100 | CD8 T cytotoxic | 7 |
| Monocyte classical | 94 | Plasmacytoid DC | 10 |

**The split is not random. Agreement is high exactly where a cell type runs one
large coordinated programme** — cell cycle for the two cycling populations, and
maturation for late erythroid. A whole gene set shifts to one end of the ranking,
so the strict method sees it; the same genes crowd the top 200, so the permissive
method sees it too.

**Agreement collapses across the T cell compartment, and the reason is
mechanical.** Naive CD4 T cells alone are 10,765 of 38,430 cells — 28% of the
dataset. When a cell type is compared against "everything else", the everything
else for any T subset is full of other T cells, so the contrast is small and no
gene set shifts far enough for the strict test. The permissive test still returns
terms because it only needs 200 genes to overlap something.

**This is the clearest demonstration in the whole analysis that a long list of
significant terms is not evidence of a strong signal.** Naive CD4 T cells
returned 107 significant terms and not one of them survived the stricter test.

By catalogue: GO produces very large numbers of one-method-only terms and almost
no agreement, while KEGG and Reactome carry nearly all of it. GO's 5,406 heavily
overlapping terms are easy for a 200-gene list to hit and hard for a whole set to
shift, and the size limits of the strict method exclude many of them outright.

---

## 10. What this work established

1. **The 19 cell type labels are confirmed three independent times** — by curated
   gene sets, by signalling footprint activity, and by transcription factor
   regulons. None of the three had any knowledge of the labels. The factor
   heatmap is the clearest single piece of evidence.
2. **p53 activity tracks cell division** without being told which cells divide.
3. **Signalling activity follows cell identity, not sub-structure.** The per-cell
   maps are flat inside each cell type. Worth testing, and a real negative
   result.
4. **Two of three newly suggested cell types are confirmed** — plasma cells and
   macrophages. Mast cells remain supported by marker genes alone.

## 11. What should be treated with caution

- **"T activated" is the weakest of the 19 labels.** An independent classifier
  would not separate it, the list-based test gave only generic terms, it has no
  distinctive transcription factor, and it shares exactly one term between the
  two list-based methods. Four checks, four failures.
- **The stem and progenitor compartment is the hardest population to describe
  with these methods.** No distinctive factor, no enriched program, pale across
  the factor heatmap. Its support rests on marker genes.
- **Disease names in catalogue results are not diagnoses.** The neurodegeneration
  maps here mean high oxidative metabolism.
- **Tied ranks** affect the middle of the enrichment rankings, not the reported
  extremes.

## 12. Rules learned, worth applying to any similar analysis

1. **Read the gene list, not the term name.** Three of the most striking results
   dissolved on inspection.
2. **Check how many genes drive a result.** Significant terms here ranged from 1
   contributing gene out of 17 to 72 out of 97.
3. **Compare a feature against itself, not against other features.** Raw activity
   averages put one dominant factor first in all 19 cell types; comparing each
   factor across cell types is what revealed the biology.
4. **A long list of significant terms is not evidence of a strong signal.**
5. **Standardising for display destroys the evidence that something has no
   signal.** Check the raw spread before scaling a figure.
6. **Filters need auditing.** Two ribosomal genes escaped a name-based filter and
   produced a false top result until they were caught.

---

## 13. Files produced

Tables:

- the full gene ranking and the top-gene list per cell type
- over-representation results for the marker lists, all and significant
- enrichment results per catalogue, combined, significant, and with the
  contributing-gene count added
- over-representation results for the unsupervised programs
- per-cell activity matrices for the 14 signalling pathways and the 293
  transcription factors, plus their per-cell-type averages
- the agreement counts between the two list-based methods

Figures:

- per-cell maps of the 14 signalling pathways
- signalling pathway activity by cell type
- transcription factor activity by cell type
- agreement between the two list-based methods

Nothing in this work modified the annotated dataset. Every step reads it and
writes elsewhere, so any single step can be re-run on its own.
