# Cell-cell communication in bone marrow

A plain-English account of what this stage asked, how it was done, what came
out, and how much of it should be believed.

---

## The question

Cells talk to each other by making a protein that another cell can catch. The
protein that is sent is a **ligand**. The thing that catches it is a
**receptor**. A ligand-receptor pair is one line of conversation.

We cannot watch that conversation happen. What we have is a table of which genes
each cell switched on. So the question has to be turned into something the data
can answer:

> Which cell types express a known ligand strongly, while some other cell type
> expresses its matching receptor strongly, and do so in a way that is
> particular to those two cell types rather than being true of everything?

That is the whole method. It does not measure signalling. It measures
co-expression of a ligand and a receptor that a database says can bind, between
two labelled groups of cells.

Keeping that sentence in mind is the difference between reading these results
correctly and over-reading them.

---

## The data

- 38,430 bone marrow cells from eight healthy donors.
- 19 cell types, labelled earlier in the project and **not changed here**.
- Log-normalised expression. No dimensionality reduction is used at all — this
  stage works on labels and gene counts only.

Nothing in this stage writes back into the shared data file. Everything it
produces lands in its own results and figures folders.

---

## The tool, and the two things it separates

The software is **LIANA**. It keeps two things apart that are easy to confuse.

**The resource** is the list of which ligand binds which receptor. It is a
catalogue built by people reading papers. Different groups built different
catalogues.

**The method** is the arithmetic that scores one ligand-receptor pair between
one sending cell type and one receiving cell type.

We ran **two resources** on the same cells with the same method, on purpose:

| resource | what it is | pairs after filtering |
|---|---|---|
| `consensus` | LIANA's merged list, several databases combined | 465 |
| `cellphonedb` | one curated database, stricter — every subunit of a protein complex must be expressed | 114 |

Running both answers a question that one resource cannot: **how much of the
answer comes from the catalogue rather than from the cells?** The same trick was
used earlier in the project with three gene set catalogues, and for the same
reason.

The method is `rank_aggregate`, which runs seven different scoring formulas and
merges their rankings into two numbers per interaction:

- **magnitude** — how strongly the ligand and receptor are expressed.
- **specificity** — how particular that expression is to these two cell types
  rather than being everywhere.

Both are on a scale where **small is strong**. An interaction is called
significant only if **both** are at or below 0.05.

---

## Settings

| setting | value | why |
|---|---|---|
| expressed-in fraction | 10% | a gene must be detected in at least a tenth of the cells of a type to count as expressed there |
| minimum group size | 10 cells | below this a cell type is dropped from a run |
| permutations | 1,000 | how the null distribution is built |
| random seed | 7 | the same seed as the rest of the project |
| significance | magnitude ≤ 0.05 **and** specificity ≤ 0.05 | both, not either |

The ten-cell floor excludes nothing when all 38,430 cells are used — the
smallest label is Platelet at 103 cells. It matters only when the data is split
by donor, which is exactly what it was chosen for.

The catalogues are downloaded once and kept on disk. The analysis therefore
reproduces without a network, and the catalogue cannot silently change between
runs.

---

## What was run, in order

Seven steps. Each one writes tables and figures named after itself.

### 1. Shared helper

Paths, settings, data loading, catalogue caching, table and figure saving, and
the dark figure styling. Never run on its own.

### 2. Ranking

Scores every ligand-receptor pair between every pair of cell types, once per
catalogue, and saves three tables per catalogue: everything tested, the
significant subset, and a count of significant interactions for each
sender-receiver pair.

### 3. Comparing the scoring methods

Reads the saved tables and asks whether the seven methods inside the merged
score actually agree. It recomputes nothing — every method's own score is
already in the table from step 2.

### 4. Splitting by donor

Runs the identical analysis on each of the eight donors separately, then asks
which results appear again and again. This is the only independent check in the
stage.

### 5. Grouping interactions into pathways

Takes the significant interactions and asks which known biological pathways
contain both the ligand and the receptor, using two curated catalogues that were
already cached earlier in the project.

### 6. Finding shared programs

Factorises the cell-type-by-interaction table to ask whether several cell types
use the same *mix* of interactions — a "signalling program" rather than one
interaction at a time. Done twice, once for sending and once for receiving.

### 7. The figure set

Seven cross-cutting figures and two summary tables, drawn from the saved outputs
of everything before it.

---

## What came out

### How much survived

| | cellphonedb | consensus |
|---|---|---|
| interactions tested | 3,529 | 16,438 |
| **significant** | **20** | **195** |
| share of tested | 0.57% | 1.19% |
| distinct ligand-receptor pairs | 5 | 34 |

**Requiring both cutoffs is a real filter, not a formality.** In the larger
catalogue 955 interactions pass the strength cutoff and 1,279 pass the
particularity cutoff, but only 195 pass both. If the two were unrelated, 74
would be expected to pass both by chance. So the two measures agree, but only
mildly.

The practical meaning is worth stating plainly: **on this data, "strongly
expressed" and "particular to these two cell types" are nearly different
questions.** Most of what is loud is not particular, and most of what is
particular is not loud. Reporting either one alone would have produced five or
six times as many "findings", almost all of which would fail the other test.

### The two catalogues do not give two views of one answer

| | count |
|---|---|
| significant in both | 14 |
| cellphonedb only | 6 |
| consensus only | 181 |

**Those 14 are the most defensible result of the stage.** Two catalogues, built
by different people from different literature, both putting the same
interactions in the top 1% of their own table. They are:

- **APP → CD74**, sent by plasmacytoid dendritic cells to eight different
  receivers.
- **TNFSF13B → HLA-DPB1**, sent by four myeloid cell types into dendritic cells.
  TNFSF13B is BAFF, the B cell survival factor. This is textbook biology and it
  is one of the few results here that is neither a housekeeping gene nor a
  marker gene.
- **PF4 → CXCR3**, sent by platelets to T proliferating cells and to
  plasmacytoid dendritic cells.

The disagreements are explained by size rather than by biology. A score that
ranks an interaction against the others in its own table will place the same
interaction differently in a 114-pair table and a 465-pair table, without its
expression changing at all.

**Reading rule: whether an interaction is "significant" is a property of the
catalogue it was ranked in, not of the interaction.** Always say which.

### Four gene families carry 80% of the result

| family | rows | share |
|---|---|---|
| B2M and HLA class I | 80 | 41% |
| S100A8, S100A9, S100A12 | 76 | 39% |
| everything else, 32 other ligands | 39 | 20% |

Both big blocks turned out to be restatements of things already known.

**The class I block is true and uninformative.** The most frequent single result
is `B2M → KLRD1`, 15 times, and in every one of them the receiver is a natural
killer cell. KLRD1 is a genuinely NK-restricted receptor, which is why. But B2M
and HLA class I are expressed by essentially every nucleated cell in the body.
So the statement the method is making is "every cell type presents class I to NK
cells" — which is textbook immunology and says nothing about which cells in
*this* marrow talk to each other.

**The S100 block is the marker list arriving by a different route.** S100A8,
S100A9 and S100A12 are the top three marker genes of the classical monocyte
label from the start of the project. A method that rewards strong and restricted
expression will necessarily return the loudest, most restricted genes of the
largest myeloid population. That does not make the result false — these are
genuine damage-associated signals acting through the right receptors — but it
cannot be used as evidence of anything the marker analysis had not already
established.

### Who sends and who receives

Two cell types send 45% of everything: **Monocyte classical 55 interactions and
Granulocyte precursor 32.** One cell type receives 24% of everything: **NK, 47.**

Neither is a measure of activity. Counting interactions per cell type measures
**how many loud, well-annotated ligands a population happens to own**, and the
two extremes are the two blocks above.

The clearest illustration: **Granulocyte precursor sends 32 interactions and
receives zero.** It is the only label with nothing incoming at all. It expresses
S100A8 and S100A9 loudly while expressing almost no annotated receptor above the
10% threshold. The asymmetry is about which receptors are in the catalogue and
detectable, not about biology.

### The stem cells are invisible

**HSPC — the stem and progenitor population — sends nothing and receives one
interaction.** This is the third consecutive stage of the project in which it
returns almost nothing: no distinctive transcription factor, no enriched gene
program, and now no signalling.

The reason is the same all three times. HSPC's identity genes — SPINK2, AVP,
CRHBP, CD34 — are markers of *stemness*. They are not members of any pathway,
regulon or ligand-receptor catalogue. **A cell type can be real, correctly
annotated and biologically central, and still be invisible to every method that
works by looking things up in a catalogue.**

### The sharpest positive result, and its caveat

The three most *particular* interactions in the whole dataset are all sent by
platelets:

| from → to | pair |
|---|---|
| Platelet → T proliferating | PF4 → CXCR3 |
| Platelet → NK | CLEC1B → KLRB1 |
| Platelet → Plasmacytoid DC | PF4 → CXCR3 |

Their particularity scores are one to three orders of magnitude better than
anything else, and `PF4 → CXCR3` also holds the top spot in the second
catalogue independently. PF4 and CLEC1B are made by essentially nothing except
megakaryocytes and platelets, so when they appear, the sender is unambiguous.

This is the exact mirror image of the class I block: a universal ligand gives a
true but useless statement, a restricted ligand gives a sharp one.

**The caveat is size.** Platelet is 103 cells, 0.27% of the dataset. Its
receivers here are also small — T proliferating is 47 cells, plasmacytoid
dendritic cells 330. This result rests on three small populations at once.
Strength of a signal and number of cells behind it are independent, and the
first does not rescue the second.

### Two results that are not expressed at all

Two of the 195 have expression products three orders of magnitude below
everything else in the table — 0.007 and 0.009 against a median of 3.74.

They passed because the score is a **rank**, not a threshold. A gene detected in
12% of one cell type at one count per cell clears the 10% filter, then scores
extremely well on particularity because it appears nowhere else, and the merged
rank is dragged along behind it.

One of them, ERFE, is a genuine erythroid hormone, so this is not a wrong
catalogue entry. It is a wrong *level*: the interaction cannot be doing anything
at that expression.

**Rule: read the expression values before believing a rank.** A rank says where
something sits relative to the rest of its own table. It does not say that
anything is expressed.

---

## Do the seven scoring methods agree?

This was worth asking because "seven methods merged" sounds like seven
independent opinions. It is not.

**Two of the seven are the same ranking.** Their correlation is 0.9999999999 —
not "highly correlated", literally one ordering appearing twice, because one
score is a monotone transform of the other. A third joins them at 0.95, being
the same quantity computed with a mean instead of a product.

So the merged score is closer to **two measurements, one counted three times and
the other four.**

**The merged particularity score is a faithful summary of its inputs** —
correlations of 0.74 to 0.89 with each method that feeds it.

**The merged strength score is not** — 0.55 to 0.57, while its three inputs agree
with each other at 0.95 to 1.00. The cause is that the merging formula saturates:
the merged strength value has 9,279 distinct values across 16,438 rows, with a
single value shared by 3,648 of them. A number with that much mass at one point
cannot follow a continuous input closely. **Read it as a coarse ordering, not as
a score.**

**One of the seven scores is useless here and must be ignored.** The permutation
p-value is 0.0 for every significant result — but it is also 0.0 for 64% of
everything tested. A test that calls two thirds of the table significant
separates nothing. At 38,430 cells the permutation null is so tight that almost
any expressed pair beats a shuffled label. It also has only 820 distinct values
in the whole table, which means any "top 100 by this score" is really just the
first 100 rows of the file in whatever order they were already in.

---

## Does it reproduce across donors?

The same analysis was run on each of the eight donors separately. This is the
only genuinely independent check in the stage, and it produced the most
important finding.

### One of the 19 cell labels is probably not a cell type

| label | BM1 | BM2 | BM3 | BM4 | BM5 | BM6 | BM7 | BM8 | total |
|---|---|---|---|---|---|---|---|---|---|
| T activated | 0 | 1 | **532** | 0 | 17 | 2 | 1 | 30 | 583 |

**532 of its 583 cells — 91.3% — come from one donor.** Four donors contain zero
or one cell of it.

This is not a rare cell type unevenly sampled. It is a population that exists in
one sample and essentially nowhere else, which was given a cell-type name.

Three earlier stages of the project had already flagged this label as the weakest
of the 19: the automatic annotation tool would not separate it, gene set tests
returned only generic terms, and it has no transcription factor of its own. Those
said it has no distinctive biology. **This one says it has no distinctive donor
spread either.** And the two genes it rests on, CD69 and BTG1, are
immediate-early genes — exactly what a sample handled slightly differently would
show.

The most economical reading is that it is a handling artefact in one donor, not
a cell type. The labels stay fixed by project rule, but every result involving it
should be read as a result about that one donor.

**Nothing before this point in the project could have found it**, because nothing
before it split the data by donor.

### Donor composition is very uneven, even though the batch effect is weak

An earlier stage established that the donors mix well — cells of the same type
from different donors sit together. That is still true. But **how many cells of
each type each donor contributed** is a different quantity, and there the donors
differ enormously. An even split across eight donors would be 12.5%; the median
largest-donor share of a label is 25%, and NK ranges from 54 cells in one donor
to 666 in another.

**Transcriptional similarity and compositional balance are independent.** A
dataset can have cells that mix perfectly and donors that contributed completely
different cell mixtures. This one does.

### How well the pooled result reproduces

Of the 195 pooled results:

- **74 (38%) are found by all eight donors.**
- 117 (60%) by six or more.
- 44 (23%) by two or fewer.

The distribution is bimodal — a big peak at eight and a second lump at two, with
little in between. That is the useful shape, because it means a threshold
separates two real groups rather than cutting through a continuum.

**A single donor over-reports badly.** Outside the pooled 195, **307
interactions were called by exactly one donor and by no other.** Eight
independent samples of the same tissue, and 307 results that appear in one and
only one of them. That is the noise floor of a single-sample analysis at this
depth, measured directly.

### Zero reproductions is not always a refutation

Three pooled results appear in no donor at all, and they fail for three
different reasons. Each had to be checked against which labels survived the
ten-cell floor in which donor before it could be read.

- **`ERFE → LRP1` is refuted.** Both cell types clear the floor in all eight
  donors, and not one calls it. Eight chances, eight failures. Combined with its
  near-zero expression, this one is settled and should be dropped from any list
  of findings.
- **`PF4 → CXCR3` into T proliferating is unverifiable, not refuted.** Platelet
  clears the floor in five donors and T proliferating in two — and those two sets
  do not overlap. **There is no donor in which both are testable.** No donor
  could have found it whatever the biology. That is a harder outcome than a
  refutation, because more cells are the only fix.
- **`HMGB1 → CXCR4` is weak.** One chance, one failure. Thin either way.

**Rule: count the opportunities before reading a zero.**

### The pooled analysis misses things every donor finds

This runs the opposite way to everything above and is the sharpest methodological
result of the stage.

**31 interactions are found by all eight donors and are not significant in the
pooled analysis.** Another 31 are found by seven.

Every one of the 31 passes the strength cutoff. None passes the particularity
cutoff — and all of them fail it by between 0.004 and 0.009. They are not
scattered; they sit in a tight band just above the line.

That is a statement about the cutoff, not about the interactions. And their
receivers explain why: **28 of the 31 point into mature T cell subsets.** In the
pooled data those subsets compete for "particularity" against a naive CD4 T cell
label of 10,765 cells — 28% of the whole dataset. A signal into memory CD8 T
cells has to look restricted against a T compartment that is mostly one enormous
label. Per donor the mix shifts and the same interactions clear.

The interactions themselves are the expected ones: class I into CD8, which is
textbook, and its absence from the pooled result is a false negative.

**General point: a pooled analysis of a dataset with one very large cell type
will under-report signals into that cell type's smaller relatives.** Splitting by
donor does not fix that — it exposes it.

---

## Do the results group into pathways?

The significant interactions were matched against two curated pathway
catalogues. A pathway counts as matched only when its gene list contains **both**
the ligand and the receptor — requiring only one side would put nearly
everything into broad immune terms.

### A fifth of the result has no pathway, and it is the wrong fifth

**158 of 195 mapped (81%).** The 37 that did not are not a random tail:

- **All nine `APP → CD74` interactions are unmapped** — the block that is
  significant in *both* catalogues and is part of the defensible core.
- **All fourteen `S100A8 → CD69` interactions are unmapped** — the second most
  frequent pair in the whole table.

So **the pathway view is blindest exactly where the evidence is strongest.**
Being in a curated pathway and being reproducible are unrelated properties, and
the first is a property of the catalogue rather than of the biology.

### What did map is the same two blocks again

**94% of the mapped interactions carry a class I/II or an S100 ligand.** The
answer to "do these results group into coherent pathways" is: they group into the
two blocks we already knew about, wearing pathway names.

### Forty-five pathway names are fifteen things

Collapsing terms that hold the *identical* set of interactions, 45 named
pathways turn out to be 15 distinct results. Eight different pathway names carry
the same 18 interactions. Eleven carry the same 6. That is a hierarchical
catalogue counting one result once per level of its tree.

**Any count of "pathways found" from a nested catalogue is inflated, here by
roughly threefold.**

### Three pathway names say the opposite of what the interaction is

This is the reading rule from the earlier pathway stage firing again, harder.

- A term called **"Class I MHC Mediated Antigen Processing & Presentation"** is,
  in full, `S100A8/S100A9 → CD36` — 18 interactions with no MHC gene among them.
  Both molecules happen to be listed in that pathway's gene set, so the match is
  correct by the rule and wrong as a description.
- A term called **"Signaling by Interleukins"** is `S100A12 → CD36`. Neither
  molecule is an interleukin.
- A term called **"Gastric cancer"** is the single unexpressed `CDH1 → LRP5`
  interaction.

**A pathway label produced this way is a statement about two genes appearing in
the same list, not about mechanism.** Every pathway name has to be opened and
checked against its contents before it is repeated.

### Ranking pathways by hit count is wrong, and so is ranking by ratio

Counting hits promotes big pathways: eleven of the 45 are at or below the
background rate, including one with 23 hits that is *depleted*. So a ratio was
added — the share of a pathway's tested interactions that are significant,
divided by the overall rate.

But the ratio has its own failure: a term with one tested interaction reaches
84-fold on one hit, and that hit is the unexpressed artefact. **The ratio and the
count must be read together. Neither column is safe alone.**

### What survives

Three results, after removing duplicates, mislabelled terms, depleted terms and
the two known blocks:

1. One large term whose name genuinely describes its contents: class I molecules
   into NK and CD8 receptors.
2. The chemokine terms carrying PF4 and CCL5 — the only pathway-level result the
   two catalogues reach independently.
3. TNFSF13B into antigen-presenting cells, under two disease-flavoured names but
   with correct underlying biology.

---

## Are there shared signalling programs?

A factorisation was used to ask whether several cell types use the same *mix* of
interactions, rather than looking at one interaction at a time. It was run twice,
once over sending and once over receiving.

**The number of programs is not supported by the data.** The fit improves
smoothly as more programs are allowed, with no point where it stops improving.
Four was chosen for readability and is reported as a choice, not a discovery.
The full curve is saved so that this is checkable.

**Two sending programs are real.** One is the myeloid alarmin program — the S100
family into its four receptors — and it picks out exactly the two myeloid
populations and nothing else. The other is a cytotoxic program grouping natural
killer cells with the CD8 T cell subsets, and it is the only sending program
containing something that is not a class I molecule.

Their internal weights are nearly flat, which matters for how they are read: **a
flat program is a set of interactions that travel together, not a ranking with a
leader.**

**One program exposed a trap in the method.** Two cell types both showed a
"usage" of 1.00 in the same program. One of them was fitted well; the other could
not be fitted at all, and the rescaling that makes the numbers readable turned
that failure into a confident-looking 1.00. The fix is to keep the pre-rescaling
total alongside — without it, a cell type the method failed on is
indistinguishable from one it described perfectly. **Normalising a row
manufactures confidence.**

**The receiving side was more informative than the sending side.** Its clearest
program is S100 into CD36 on both erythroid labels. CD36 is a genuine erythroid
scavenger receptor, the weights supporting it are the largest in that direction,
and it is the only place the late erythroid label appears anywhere in the stage.
That is the one genuinely new piece of biology here.

**The stem cells are placed for the first time, and it means nothing.** They come
out split almost evenly between two programs, on one received interaction, with
the second-worst fit in the direction. An even split on a weak fit is what a
method does when it has no information.

---

## Do the growing cells and the finished cells talk?

This is the question the project chose to ask, because an earlier stage
established that the mature T, NK, B and monocyte cells in this sample did not
grow from the stem cells present — they were collected alongside them. So the
19 labels were split in two:

- **Developing in the marrow**: stem cells, granulocyte precursors, early and
  late erythroid, large and small pre-B.
- **Finished, collected rather than produced here**: both monocytes, dendritic
  cells, plasmacytoid dendritic cells, platelets, naive B, the four T groups,
  proliferating T and NK.

The answer:

| direction | interactions | share |
|---|---|---|
| into the **finished** compartment | **178** | 91.3% |
| into the **developing** compartment | **17** | 8.7% |

The developing side holds 6 of the 19 labels, so an even spread would give it
about 62. It gets 17 — **3.6 times fewer than its share.**

**Thirteen of the nineteen cell types send nothing into the developing
compartment at all.** The traffic that does exist runs through exactly two
channels:

- **Myeloid alarmins into erythroid CD36** — 9 of the 17, from classical
  monocytes and granulocyte precursors.
- **`APP → CD74` from plasmacytoid dendritic cells** — 3, into stem cells and
  both pre-B labels.

The remaining five are two S100A8 → CD69 hits into small pre-B, two HMGB1 → CXCR4
hits inside the B lineage, and one artefact.

**Read carefully, this is as much about receptors as about signalling.** The
developing labels express few annotated receptors above the detection threshold,
so there is little for any ligand to pair with. The granulocyte precursor, which
sends 32 and receives 0, is the extreme case.

---

## The figures

Thirteen figures in total. What each is for, and how to read it.

**Ranking heatmaps, one per catalogue.** Drawn on a shared colour scale so the
two are comparable, which makes the stricter catalogue look nearly empty — and
that emptiness is the result. The larger catalogue shows two bright rows and one
bright column, which a reader would take as "monocytes talk to everything and
everything talks to NK". Both halves of that are artefacts of the two gene
blocks.

**Ranking dotplots, one per catalogue.** The top results occupy only six distinct
ligand-receptor pairs, and the most frequent one is also the least particular.
The trade-off between strength and particularity is visible in one panel.

**Method agreement matrices, one per catalogue.** A nine-by-nine grid. The
merged-strength row is the palest in both, which is the finding drawn: it sits
apart from every method including the three that feed it. Nothing in the grid is
negative, which is also worth seeing — the measures are weakly aligned, not
opposed.

**Donor reproducibility bars, one per catalogue.** Three things to read. The
single-donor bar towers over everything at 307 — the noise floor. The bars for
the pooled results climb steadily to 74 at eight donors — the result being
confirmed. And at seven and eight donors there is a second bar almost as tall,
which is the false negatives. A reader who looks at only one of the three draws
the wrong conclusion.

**Pathway heatmap.** Senders against pathways, ordered left to right by
enrichment. The most useful thing about it is that the brightest cells are *not*
on the left: the big myeloid blocks are only 3 to 7 times enriched while the
small sharp terms are 9 to 11 times. Size and enrichment point in opposite
directions. The three right-hand columns are below the background rate and the
caption has to say so.

**Signalling program river plot.** Cell types on the left, programs in the
middle, the interactions each program uses on the right. Three of the four
programs read cleanly. The fourth contains a cell type the method could not fit,
drawn at full width like any other, which the caption has to flag.

**Circle and chord plots.** The same flow two ways. The circle is readable; the
chord, at 90 flows between 19 arcs, is past the density where an individual flow
can be traced, and is kept mainly because it is the standard view.

**Count heatmap with marginal bars.** The clearest single statement in the set:
the stem cell row is empty, the late erythroid row is empty, and the granulocyte
precursor column is empty.

**Hierarchy plot.** Two panels, signalling into each compartment. The left is
nearly empty and the right is dense. That picture *is* the 17-against-178
finding.

**Top interactions dotplot.** One row spans 13 senders all pointing into NK — the
class I block in a single line — and it is the darkest row, meaning the least
particular.

**Sending against receiving.** One row per cell type with the two totals joined,
ordered by their difference. **This figure was redesigned twice.** As a scatter
it put 10 of the 19 cell types into a corner worth an eighth of the panel;
logarithmic axes made it worse, and automatic label placement still produced
overlapping labels and crossing leader lines. One row per cell type removes the
problem entirely, because labels become axis ticks and cannot collide. The
general lesson: **a crowded scatter is the wrong chart, not a chart with a label
problem.**

**Catalogue comparison.** One bar per catalogue per ligand-receptor pair. It is
almost entirely one colour, which is the result.

### The caveat that belongs on five of them

The circle, chord, count heatmap, hierarchy and sending-against-receiving figures
all weight their edges by **the number of significant interactions**, because
LIANA produces no communication strength to weight them by. The ranking step
established what that number measures: how many loud, well-annotated ligands a
population owns.

**A thick edge in those figures means many catalogue entries matched, and nothing
stronger.** Without that line the figures read as if they said more than they do.

---

## What to trust, and what not to

**Trust:**

- The three interaction blocks found by both catalogues — APP → CD74 from
  plasmacytoid dendritic cells, TNFSF13B into antigen-presenting cells, and PF4
  from platelets.
- The 74 interactions found by all eight donors independently.
- The finding that the developing compartment receives almost nothing, which is
  consistent across the counts, the hierarchy figure and the summary table.
- The erythroid CD36 receiving program, which is well supported and biologically
  sensible.

**Do not trust:**

- Any count of interactions per cell type as a measure of how much a cell type
  signals.
- The class I block as evidence about which cells in this marrow interact.
- The S100 block as anything the marker analysis had not already shown.
- The permutation p-value, in any form.
- Any pathway name that has not been opened and checked against its contents.
- Any result involving the T activated label, which should be read as a result
  about one donor.
- The two interactions with near-zero expression, one of which the donor check
  refuted outright.

**Cannot be decided either way:**

- PF4 → CXCR3 into proliferating T cells. It is the sharpest result in the stage
  and no donor can test it, because the two cell types are never both present in
  sufficient numbers in the same sample.

---

## Rules this stage produced

1. **Read the expression values before believing a rank.** A rank is relative to
   its own table; it does not say anything is expressed.
2. **Name the catalogue with every claim.** Significance is a property of the
   catalogue the interaction was ranked in.
3. **A universal ligand gives a true statement that is not a finding.** Check
   what fraction of cell types express the ligand before reporting a hub.
4. **Counts of interactions per cell type measure annotated marker loudness, not
   activity.**
5. **Check an "aggregate" for duplicate members before calling it a consensus.**
   Two of the seven scores here are the same ordering.
6. **Any top-N statistic on a heavily tied column measures the sort order.**
   Count the distinct values first.
7. **Split by donor before trusting any label.** A label's total cell count says
   nothing about its donor spread.
8. **Count the opportunities before reading a zero.** Check whether both cell
   types were even testable.
9. **Counts from different tables are not comparable.** Every run is its own
   competition.
10. **A dominant cell type creates false negatives among its smaller relatives.**
11. **A pathway name is a claim about a gene list, not about mechanism.**
12. **A nested catalogue counts one result many times.** Collapse identical
    members before reporting how many pathways were found.
13. **An enrichment ratio needs its denominator shown.** Fold and count are read
    together.
14. **Catalogue coverage is not evidence quality.** The best-supported results
    here are the ones no pathway contains.
15. **A factorisation with no elbow has no natural number of components.**
16. **Normalising a row manufactures confidence.** Keep the pre-scaling total.
17. **A crowded scatter is the wrong chart, not a chart with a label problem.**
18. **A figure fault can look exactly like a data finding.** Two cell types
    missing from one figure read as a real zero and were a coordinate bug; the
    tables settled it.

---

## What the stage produced

**Tables**, grouped into one folder per step: the full and significant
interaction tables for both catalogues, the sender-receiver counts, the method
agreement statistics, the per-donor tables and coverage record, the pathway
assignments and per-pathway summary, the program usage and weights with the fit
curve, and the two cross-cutting summaries.

**Thirteen figures**, all on a dark background at print resolution.

**The full record** — every number, every figure reading and all the reasoning
behind the statements above — is in the stage report kept with the project
documents. This note is the summary; that file is the evidence.

---

## Left open

- Whether the 0.05 particularity cutoff should be reported at all, given that 31
  interactions all eight donors agree on sit just above it. Changing the cutoff
  would be fitting the threshold to the answer, but the band should be shown.
- Whether the weak label should be flagged in the annotation write-up as
  possibly not a cell type.
- Whether the receiving programs deserve their own river plot; they are currently
  only in the tables.
- Three cell types known to be missing from the 19 labels — plasma cells, mast
  cells and basophils, and macrophages. Plasma cells in particular would be
  expected in any B lineage conversation and cannot appear.
