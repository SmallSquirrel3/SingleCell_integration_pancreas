> **Superseded (second pass, 25.09.2026).** This is the first-pass narrative, kept for the
> record. Its claims were tested in `08_error_anatomy.ipynb` with markers, doublet scores,
> library size, per-technology ambient floors and five-seed stability. Outcome: the
> "acinar-to-ductal metaplasia" reading is **not supported** (the 193 cells are `ductal`
> cells in the acinar cluster, ~90 % from two Baron batches, with an acinar enzyme program
> and weak ductal markers; the 3 stable acinar → ductal errors are low-content acinar cells);
> "alpha → gamma survived four runs" is **not supported** (3 stable cells, 20 seed-dependent
> ones; counts per run hid near-complete turnover of *which* cells); the alpha → acinar
> errors are **annotation errors in the source**, not model errors; the beta → delta errors
> are mostly **polyhormonal profiles** with normal library size. The final-numbers table
> below is superseded by the README (single-run values are within 0.01 of the five-seed means in `07`).
> The stellate/immune merges are resolution effects, as stated here.

# 06. Error Analysis — What Actually Went "Wrong", and Why That's Interesting

This pulls together everything from `02_baseline_integration.ipynb` (the Harmony cluster
composition at resolution=0.3) and `04_scanvi.ipynb` (the honest label-transfer test on
`smartseq2`).

## 1. What merged at coarse resolution, and why it's not a failure

At resolution=0.3, Harmony found 12 clusters instead of 14 known types. My first instinct
was that this meant the model lost information. Looking at the actual crosstab, that's
not really what happened — the merges are specific, not random:

| Merged into one cluster | Why this makes biological sense |
|---|---|
| `activated_stellate` + `quiescent_stellate` + `schwann` | Activated and quiescent stellate cells are two *states* of the same cell type, not two cell types — asking a clustering algorithm to separate them is arguably asking too much |
| `gamma` + `epsilon` | Both rare hormone-producing islet cells (pancreatic polypeptide, ghrelin) — close relatives, small numbers |
| `macrophage` + `mast` + `t_cell` | 128 cells total out of 16,382. At this resolution, not splitting these is expected, not a bug |
| `acinar` + `ductal` (partially — 193 cells) | Possibly acinar-to-ductal metaplasia, a documented biological process and, as it happens, an early step toward pancreatic cancer |

So: fewer clusters than known labels doesn't automatically mean the method is worse than
it looks. Sometimes it's revealing something the discrete annotation smooths over.

## 2. The honest label-transfer number: ~97.5%

Labels for every `smartseq2` cell (2394 of them) were hidden during scANVI training — the
model never saw a single true label from this technology. The accuracy above is computed
only on these held-out cells.

Several rare types were predicted with zero or near-zero errors across every run I did,
which I did not fully expect given how few examples of them exist elsewhere: `delta`,
`endothelial`, and `macrophage` in particular.

### The two error patterns that kept showing up, run after run

- **acinar → ductal**: somewhere between 3 and 6 cells out of 188 across the four runs I
  did. Same signal as in section 1, and it showed up in both unsupervised clustering
  (Harmony) and supervised classification (scANVI) every single time — a decent argument
  that it's biology, not coincidence or noise from one unlucky training run.
- **alpha → gamma**: between 4 and 17 cells out of roughly 1000, varying more than the
  acinar/ductal pair but never absent. Both are hormone-producing islet cells from a
  shared developmental lineage, so some confusion between them is expected rather than
  a red flag.

## 3. Final numbers, all methods

| Method | ARI | NMI |
|---|---|---|
| scANVI | 0.950 | 0.926 |
| scVI | 0.947 | 0.914 |
| Harmony | 0.914 | 0.888 |
| baseline (no integration) | 0.428 | 0.708 |

(Harmony's numbers are exactly reproducible; scVI and scANVI drift by a percentage point
or so between runs depending on random initialization — the ranking above held across
every run I did, even if the third decimal didn't.)

The gap between "no integration" and "any integration" is the headline result. The
ordering among the three integration methods is consistent and matches expectations
(classic → unsupervised DL → semi-supervised DL, each a bit better), but it's a modest
effect compared to that first jump.

## 4. Limitations, stated plainly

- The data is a curated benchmark, not raw sequencer output — no mitochondrial genes,
  and `counts` reconstructed from size factors rather than true integer UMIs (rounded for
  scVI/scANVI compatibility)
- Leiden resolution was tuned against the known annotation — fine for comparing methods,
  not the same as a fully unsupervised result
- "Metabolic disease relevance" here is about the tissue (islets, T2D) and the method
  (integration benchmarking), not about the omics layer — this is transcriptomics, not
  metabolomics, and I'd rather say that upfront than have it assumed
- scVI and scANVI numbers vary slightly between training runs (random weight
  initialization) — the differences reported above are from a single run each; the
  direction of the comparison held up across the two runs I actually did, which is
  reassuring, but it's still worth flagging as a source of noise rather than pretending
  the third decimal place is meaningful