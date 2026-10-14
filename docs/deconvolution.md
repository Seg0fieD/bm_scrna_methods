# Deconvolution — estimating cell-type shares from bulk RNA-seq

## The question

Bulk RNA-seq measures one blended profile per sample: the RNA of all its cells
mixed together. Deconvolution works backwards from that blend to the share of
each cell type inside it.

This stage asks whether a model trained on our own bone marrow single cells can
do that. It was tested twice: first on artificial bulk samples where the true
shares are known, then on real bulk RNA-seq from another lab.

## The tool

**Scaden** is a small neural network. It learns from thousands of artificial
mixtures, each built by pooling a few hundred single cells with known shares,
and then predicts shares for new samples. It trains three networks of
different sizes and averages them.

Scaden was later rewritten in PyTorch so it could run on the GPU. The rewrite
was checked against the original before it was trusted (see below).

## How agreement was measured

- **CCC (concordance)**: 1 means predicted and true shares are identical; it
  drops when they disagree in order *or* in level. This is the main score.
- **Pearson r**: only asks whether samples are put in the right order.
- **Bias**: average over- or under-prediction.

A model can order samples well (high r) and still get the amounts wrong (low
CCC). Both are reported.

Results are given for the 19 cell-type labels and for 8 broad groups (HSPC,
Erythroid, Myeloid, Dendritic, Platelet, B lineage, T cells, NK), which are
sums of the labels.

## Part 1 — artificial bulk with known truth

### Setup

- **Donor split.** Six donors (MantonBM2–BM7) provide the training cells; two
  (MantonBM1 and MantonBM8) are held out for testing, so no test cell is ever
  seen in training. This pair was chosen because it is the most typical pair
  and keeps cells of every label on the test side.
- **Training mixtures.** 8,000 mixtures of 500 cells, with random shares and
  "partial" shares (some cell types left out).
- **Test mixtures.** Three sets of 500: random, partial, and **donor-like**
  (shares close to each test donor's real composition). The donor-like set is
  the realistic one.
- We built the mixtures ourselves instead of using Scaden's own builder, so the
  donor split and the cell counts are under our control.

### Setting up Scaden

Scaden 1.1.2 does not work with the current Keras. Training crashed at the save
step and, worse, prediction **failed silently**: it built a fresh, untrained
network and returned plausible-looking numbers. Installing the older Keras
(`tf-keras`) and switching it on fixed both. Scaden also quietly reuses an old
model if one is on disk, so each run now clears it first.

### What Scaden got right and wrong

- **It orders samples well but squeezes the amounts.** On the donor-like set
  the typical label had r 0.77 but CCC only 0.40.
- **Predictions are pulled toward the average.** Big populations are
  under-predicted, and cell types that are absent still get about 2%. The
  random training mixtures teach the model that every type is usually present
  at about 1 in 19.
- **Fine T-cell labels cannot be told apart in bulk.** T proliferating and
  Platelet were not learned at all (too few cells).
- **Reliable outputs**: the Myeloid and B lineage groups, and the HSPC,
  Erythroid early, B naive and Monocyte classical labels.

### Why one donor scored badly

MantonBM8 scored much worse than MantonBM1: part of its T-cell share was
predicted as Erythroid late and T activated. Two explanations were tested:

- **Free-floating haemoglobin RNA in the sample** — rejected; MantonBM8 is
  average.
- **A stress signal in its T cells** — supported. MantonBM8's T cells carry the
  strongest "immediate-early" gene signal, a quick response to handling. The
  T activated label is built on those same genes, so the model reads stress as
  activation.

So **T activated is not a trustworthy bulk output**: it reflects how a sample
was handled as much as what cells it holds.

### The PyTorch version

Scaden only ran on the CPU here, which made trying variants slow. The PyTorch
rewrite was checked two ways:

- given Scaden's own trained weights, it gives the same predictions (to about
  one in a million);
- trained from scratch, it differs from Scaden only as much as two of its own
  runs differ from each other.

It is the same method, about 13 times faster per network.

### Training variants

Scaden feeds the training mixtures in blocks (all random ones, then all partial
ones), and stops after 5,000 steps although its error is still falling. Both
were changed, over five independent runs each:

| training | typical label CCC | typical group CCC | MantonBM8 groups |
|---|---|---|---|
| Scaden default (5,000 steps) | 0.40 | 0.35 | 0.26 |
| full shuffle, 20,000 steps | 0.57 | 0.58 | 0.57 |

**Training longer is the main gain**; shuffling helps a little. The pull toward
the average shrinks but does not disappear. Four labels with very few cells did
not improve at all — they are limited by data, not training.

## Part 2 — real bulk RNA-seq

### The data

No open healthy bone marrow bulk set with measured cell shares was found.
Instead we used **GSE74246** (Corces et al. 2016): 49 bulk samples of
**sorted** blood cell populations from healthy donors — stem and progenitor
cells, erythroid cells, monocytes, B, CD4 T, CD8 T and NK cells. Because each
sample is one known cell type, the right answer is known.

From these we also built mixtures in the computer: 250 random and 250
marrow-like. Their known shares are shares of **RNA**, not of cells.

Three progenitor types (GMP, MEP, CLP) sit between groups; they were reported
but not scored.

### Check without any model

Each sorted sample was simply correlated with the average profile of each of
our 19 labels. 36 of 49 samples matched the expected label best, with
correlations of 0.6–0.8, despite the different lab method. T-cell subtypes were
indistinguishable (CD4 vs CD8 differences of about 0.01), so only the T cells
group can be tested. CLP samples resemble our early B cells, which is sensible
biology.

### Results

All results here use the average of five training runs; single runs carried
large shifts that differed from run to run.

- **Pure sorted samples are recognised**: the share given to the right group
  was 0.77–0.98 (typical 0.92).
- **Random mixtures**: Myeloid, B lineage and NK were estimated well (CCC
  0.86–0.89), T cells 0.74, HSPC 0.68, **Erythroid poorly (0.25)**.
- **Marrow-like mixtures**: the overall composition was close to the truth
  (for example T cells 0.44 predicted vs 0.49 true), but the small differences
  between similar samples were not tracked well.
- **Erythroid is the main failure**: pure erythroid samples are recognised, but
  once mixed the predicted erythroid share barely moves. Not explained.
- **Dendritic cells** get 2–6% although none are present.

### Gene-length correction — tried, not adopted

The sorted bulk counts reads along the whole transcript, so long genes get more
reads; our single-cell data counts one end only. Dividing by gene length made
every sorted sample look more like its single-cell match and improved several
pure samples, but **made mixtures worse** (typical random CCC 0.80 → 0.69,
marrow-like 0.54 → 0.40; HSPC collapsed). Raw counts stay the main result.

## What to trust

| output | use |
|---|---|
| Myeloid, B lineage, NK | shares and comparisons between samples |
| T cells, HSPC | comparisons between samples; amounts with caution |
| Erythroid | not in mixed samples |
| fine T-cell labels, T activated, T proliferating, Platelet | not usable |
| Dendritic | not tested on real bulk; appears when absent |
| any share below about 2% | not reportable |

## Limits

- The real-bulk test uses sorted samples and computer-made mixtures, so it
  checks RNA shares, not cell shares, and has no real sample-to-sample
  handling differences.
- Three real cell types are missing from our 19 labels (plasma cells, mast
  cells/basophils, macrophages).
- Only two test donors; very few T proliferating (14) and Platelet (26) test
  cells; T activated comes almost entirely from one training donor.
- Training settings other than length and shuffle were not varied.

## What could come next

1. Training mixtures shaped like real marrow, to remove the remaining pull
   toward the average.
2. Find out why Erythroid fails in mixtures (haemoglobin genes are the first
   suspect).
3. Remove the stress (immediate-early) genes and re-score.
4. A real marrow bulk set with measured cell shares, if one becomes available.
5. Compare with another deconvolution method on the same mixtures.

## Figures

All in `figures/deconvolution/<script>/`.

| figure | what it shows |
|---|---|
| `pseudobulk_train_fractions.png` | share of each label across the training mixtures |
| `pseudobulk_test_fractions.png` | share of each label in the three test sets |
| `pseudobulk_score_<model>_<level>_ccc.png` | agreement per label or group and test set, for Scaden and the PyTorch version |
| `pseudobulk_score_<model>_<level>_scatter.png` | predicted against true shares, donor-like set |
| `pseudobulk_score_<model>_<level>_bias.png` | average over- or under-prediction, donor-like set |
| `donor_diagnostics_haemoglobin.png` | haemoglobin RNA outside erythroid cells, per donor |
| `donor_diagnostics_activation.png` | stress-gene levels in T cells, per donor |
| `scaden_torch_loss_passes.png` | training error over time, five runs |
| `scaden_torch_loss_last_passes.png` | training error over the last passes, showing the block-wise wave |
| `variant_score_label_trend.png`, `variant_score_group_trend.png` | score against training length, per shuffle type |
| `variant_score_labels.png` | per-label scores of the variants |
| `corces_bulk_correlation_{counts,length}.png` | sorted samples against our labels, no model |
| `bulk_score_pure_{counts,length}.png` | predicted groups for each pure sorted sample |
| `bulk_score_mixtures_{counts,length}.png` | predicted against true shares in the real-profile mixtures |

## Rules learned in this stage

1. **A tool that returns plausible numbers can still be broken.** Scaden's
   silent fallback to an untrained network looked fine; only checking that a
   trained model was actually loaded caught it.
2. **Order and amount are different questions.** Read correlation and
   agreement together.
3. **One training run is not a result.** Average several before comparing
   settings.
4. **Score on held-out donors**, never on cells the model has seen.
5. **Labels built on stress genes measure handling, not biology**, and carry
   that into every downstream method.
6. **Check the data without a model first.** The simple correlation showed in
   advance which cell types real bulk could and could not separate.
7. **Looking more alike is not the same as better estimates.** Test any
   correction on mixtures, not on pure samples.
8. **Say what the truth measures.** Mixtures of RNA profiles test RNA shares,
   not cell counts.
