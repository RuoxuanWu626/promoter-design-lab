# Promoter Design Lab

An interactive promoter-design game built on the **real human Puffin model**,
with two research modes.

| model | status |
|---|---|
| **Puffin** (human, published) | **real** — `puffin.pth`, 10 learned motif filters |
| AlphaGenome (human CAGE) | adapter written, **not yet running** — see below |
| mock Puffin / mock PuffinD / mock CAGE | kept as labelled controls |

Every panel carries a `real` or `mock` badge, and every export records which
adapter produced it.

## Getting it

```bash
git clone git@github.com:RuoxuanWu626/promoter-design-lab.git
cd promoter-design-lab
python3 backend/selftest.py          # 57 invariant checks
```

Python 3 + numpy runs the mock adapters. Real Puffin additionally needs
PyTorch; on Randi the `alphagenome` env has it and is picked up automatically.

## Running it

```bash
./run_server.sh                      # Slurm job, prints the tunnel command
python3 backend/app.py --port 8765   # or just run it locally
```

CPU only — Puffin is small enough that a construct predicts in well under a
second, and the app caches repeats.

## The real Puffin model

`models/puffin_arch.py` is the published architecture: ten motif filters
(`conv`, 4→10, width 51), an initiator branch (`conv_inr`), a trinucleotide
branch (`conv_sim`, 4→32, width 3), each run on both strands, softplus
activated, and deconvolved through width-601 kernels into ten output tracks
(FANTOM CAGE / ENCODE CAGE / RAMPAGE / GRO-cap / PRO-cap, each strand).

The checkpoint loads with **zero missing and zero unexpected keys**. The only
deviation from the reference implementation is `nn.Conv1d` in place of
`torch_fftconv.FFTConv1d`, which removes the one dependency the inference env
lacked; they were checked numerically equivalent (max abs difference 9.5e-07).

### The motif library comes out of the checkpoint

`backend/extract_puffin_motifs.py` derives the whole library from
`conv.weight` and `deconv.weight` — nothing is hand-written:

| id | Puffin name | consensus | prefers |
|---|---|---|---|
| `tata` | TATA | `CTATAAAAG` | **−31** |
| `inr` | Long Inr | `AATTTCCCCCTGGCCATCCAACGGGCCGGACG` | 0 |
| `sp1` | SP | `TGGGCGGGGC` | −52 |
| `nfy` | NFY | `AGCCAATCA` | −55 |
| `nrf1` | NRF1 | `GCGCATGCGC` | −58 |
| `ets` | ETS | `ACTTCCGGT` | −48 |
| `creb` | CREB | `ATGACGTGA` | −50 |
| `yy1` | YY1 | `AAAATGGCG` | +11 |
| `znf143` | ZNF143 | `TACATTTCCCAGAATGCATTGCG` | −65 |
| `u1` | U1 snRNP | `AGGTAAGT` | +76 |

Each filter is trimmed to its informative core and reverse-complemented into
canonical orientation. The preferred position is read off the model's own
deconv kernel — `deconv` is a cross-correlation, so an activation at *p* drives
output at *p − argmax*. Nothing about these numbers was assumed, and the model
independently reproduces **TATA at −31** and the U1 site downstream at +76.

Detection thresholds are set per motif at a fixed false-positive rate under a
realistic background, so motifs with very different information content (TATA
carries 0.46 bits per position, NRF1 1.27) are handled on the same footing.

## Two libraries, on purpose

Puffin's ten motifs are **core promoter** elements. They set how much
initiation a promoter produces and where it starts, and they do that in much
the same way in every cell type. They are *not* what makes a promoter
cell-type specific.

So the app carries a second, separate library of **lineage TF sites** —
GATA, HNF4A, C/EBP, PU.1, EBF1, SOX2–OCT4, E-box, ERE, TEAD, RUNX, p53, AP-1,
NF-κB, and REST/NRSE — in `backend/celltype_elements.py`. These drive the
cell-type model, not the profile model.

The separation is measurable, and the self-test asserts it:

| design | tau | strongest | mean activity |
|---|---|---|---|
| empty background | 0.032 | — | 0.83 |
| **+ core promoter motifs** | **0.033** | — | **1.25** |
| + one GATA | 0.185 | K562 | 1.29 |
| + one HNF4A | 0.205 | Hepatocyte | 1.33 |
| + one SOX2–OCT4 | 0.204 | H1 | 1.29 |
| + one E-box | 0.312 | SK-N-SH | 1.30 |

Core promoter motifs raise activity by half and leave tau untouched. A single
lineage site barely moves activity and flips which cell type wins.

## The three screens

### Design

Background with configurable length, GC and CpG o/e, TSS at 0. Drag elements
from either library onto the track; move, reverse, duplicate and **delete**
them (✕ on each row, *Delete selected*, *Delete all*, or the Delete key).
Placing overwrites bases, so length never changes.

**Clean canvas** (on by default) shuffles chance lineage TF sites out of the
background — a 2 kb random sequence typically carries a couple of dozen, and
they otherwise drown out anything you place. It preserves base and
dinucleotide composition, and the count is reported. Turn it off to design
against genuinely random sequence.

**Model agreement** reports Pearson r and Spearman ρ between the simple model
(Puffin) and the deep model, with a scatter and fit. Both sides are brought to
a common log1p scale first: r is invariant to scale, but a log output against a
linear one is not an affine relation, so correlating them raw would measure the
transform rather than the agreement.

### Mode 1 — where is cell-type specificity encoded?

*Position scan* slides one element across positions and reports tau, a cell
type × position heatmap, and a decomposition separating **increased target
activity** from **decreased off-target activity** — tau alone cannot tell those
apart. REST/NRSE, which represses every non-neuronal lineage, is the
off-target-silencing route; GATA and friends are the target-gain route.

*Attribution* disrupts the sequence patch by patch (dinucleotide-preserving
shuffle, so composition and CpG content survive and only arrangement is
destroyed) and scores two separate questions: what drives **activity**, and
what drives **tau**. Core promoter footprints, lineage footprints, and
everything else are tallied separately, including motifs that occur by chance
in the background rather than being placed.

A representative result — core promoter plus one GATA site:

| attribution | core promoter | lineage sites | elsewhere |
|---|---|---|---|
| activity | 80.1% | 5.1% | 14.8% |
| specificity | 47.1% | 42.5% | 10.4% |

### Mode 2 — motif interactions

Four matched constructs (background / A / B / A+B) give the position-wise
residual `I = f(A+B) − f(A) − f(B) + f(bg)`, with a pair view, a 10×10 grid
including same-motif pairs, and spacing curves across orientations.

Three controls are built in: the mock Puffin is additive by construction so its
residual is identically zero and doubles as a measure of the junction artefact;
a pair-terms-off variant separates saturation from a real pair term; and
overlapping or too-closely-spaced sites are flagged and dropped.

## Reading the numbers honestly

**Puffin needs context.** Its output is valid only where the 601 bp
deconvolutions have full sequence, so 325 bp beyond the requested window. Where
the construct cannot supply it the output is **masked, not padded** — the plot
leaves a gap and the UI says how much is missing and how long a construct would
fill it. An earlier version zero-padded, and the resulting edge artefact was
larger than any real signal.

**Occlusion attribution** measures what happens when a patch is destroyed,
which is not the same as what it contributes in context; redundant elements can
both read as unimportant.

**Tau depends on the panel and the scale.** Values from different cell-type
panels are not comparable.

**A non-zero Mode 2 residual is a statement about the model**, on the scale the
model emits — not evidence of biological cooperativity.

**One background proves nothing.** Every experiment takes a replicate count and
reports the spread.

## AlphaGenome: written, not yet running

`adapters/real_alphagenome.py` targets the local PyTorch port at
`/gpfs/data/zhou-lab/rxwu/Puffin-drozofila/tests/alphagenome_torch`
(`alphagenome_all_folds.pt`, 546 CAGE tracks = 273 cell types × 2 strands).

It is blocked on one thing: the `alphagenome` env has `torch 2.12.0+cu130`,
while the A100 nodes run driver 570.124 (CUDA 12.8), so `torch.cuda` is
unavailable there. `gpuq_env` has a compatible `cu126` build but not the
AlphaGenome packages. Either combination works once resolved:

* run on a node whose driver is ≥ 580, or
* install `alphagenome` + `alphagenome_research` into a `cu126`/`cu128` env, or
* use the arm64 env on the GH200 nodes, which is what the lab's own
  AlphaGenome jobs use.

The other thing to settle is input length: the model takes 2²⁰ = 1 Mb, so a
1–2 kb construct has to be embedded in a large context, and Mode 1 costs one
forward pass per offset. `docs/MODEL_ADAPTERS.md` has the full checklist.

## Layout

```
backend/
  app.py                 stdlib HTTP server, JSON API, background jobs
  models/puffin_arch.py  the published Puffin architecture
  extract_puffin_motifs.py  derives the motif library from the checkpoint
  motifs.py              library, PWM scanning, null-distribution thresholds
  celltype_elements.py   lineage TF sites (the second library)
  sequence.py            background generation, scrubbing, construct assembly
  experiments.py         position scans, attribution, pair/grid/spacing
  scoring.py             metrics, tau, specificity split, Pearson r
  store.py               players, challenges, leaderboard
  adapters/              real + mock adapters behind two interfaces
  selftest.py            57 invariant checks
frontend/                vanilla JS, canvas plots, no dependencies
data/                    store.json, runs, puffin_motifs.npz (generated)
```
