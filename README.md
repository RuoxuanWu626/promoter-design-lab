# Promoter Design Lab

An interactive promoter-design game with two research modes, built around
Puffin's 10-motif library.

> **Every prediction in this version is a mock.** No trained model runs. The
> adapters generate plausible-looking profiles from hand-written rules so the
> interface and the experiment designs can be built and stress-tested first.
> Real inference plugs in behind the same interfaces — see
> [docs/MODEL_ADAPTERS.md](docs/MODEL_ADAPTERS.md). The UI shows a `mock` badge
> on every panel and every exported record carries the adapter's identity.

## Getting it

```bash
git clone git@github.com:RuoxuanWu626/promoter-design-lab.git
cd promoter-design-lab
python3 backend/selftest.py          # 43 invariant checks, no server needed
```

Requirements are Python 3 and numpy — nothing else. No web framework, no CDN,
no build step, so it also runs on a compute node with no outbound network.

`data/` is deliberately empty in the repo. It is created at runtime and holds
per-deployment state: the player/challenge/leaderboard store, saved experiment
bundles, and any real motif weights you drop in. None of that is committed.

## Running it

On any machine, for a quick look:

```bash
python3 backend/app.py --port 8765   # then open http://localhost:8765/
```

On the Randi cluster, as a Slurm job (this is what the lab deployment uses):

```bash
./run_server.sh                      # submits the job, prints the tunnel command
```

It prints something like:

```
  server running on cri22cn002:8765   (job 15394802)

    ssh -N -L 8765:cri22cn002:8765 Randi

  then open  http://localhost:8765/
```

Options: `-p express -t 6:00:00` (short queue, usually starts instantly),
`-p tier1q -t 24:00:00` (long runs), `--local` (login node, quick checks only),
`--stop` (cancel). Multiple people can use one server at once — that is what
the challenges and leaderboard are for.

## The three screens

### Design

Generate a random background of a chosen length, GC content and CpG
observed/expected ratio, with the TSS marked at position 0. Drag the 10 Puffin
motifs (plus CpG-rich segments and arbitrary pasted sequence) onto the
construct. Elements can be moved, reversed, duplicated and removed; placing one
**overwrites** bases, so the construct length never changes and constructs stay
directly comparable.

Live panels show each model's predicted initiation profile (both strands),
per-cell-type activity, the individual metrics, and a weighted score whose
weights you set with sliders. Keyboard: `←/→` nudge, `shift` for 10 bp, `R`
reverse, `D` duplicate, `Del` remove. Scroll to zoom the track; at high zoom the
bases appear.

### Mode 1 — where is cell-type specificity encoded?

Slide one element (a motif, a CpG-rich segment, or pasted sequence) across a
range of positions relative to the TSS and predict the whole cell-type panel at
every position.

Outputs: tau versus position with a spread band across replicate backgrounds; a
cell type × position heatmap (absolute or change-versus-background); target
versus off-target activity on one axis; and a **mechanism strip** that labels
every position with *how* specificity moved.

That last one is the point of the mode. Tau cannot tell you whether specificity
was won by raising the target or by lowering everything else, so each position
is decomposed into those two arms separately. The mock cell types include
repressive motif weights — the U1 site lowers 8 of the 12 — precisely so both
routes exist and can be compared.

### Mode 2 — what motif interactions do the models learn?

For any pair from the library, including same-motif pairs, four matched
constructs are built on one background — background, A only, B only, A+B — and
each model's profiles are differenced:

```
I(x) = f(A+B)(x) − f(A)(x) − f(B)(x) + f(background)(x)
```

Views: a pair detail (all four profiles, the additive expectation, and the
residual), a 10 × 10 pairwise heatmap, and spacing curves with orientation
variants. Everything repeats across matched backgrounds and reports the spread.

## What the numbers mean, and what they don't

This matters more than any feature, so the app repeats it in place rather than
burying it here.

**The scale is part of the claim.** Both mock profile adapters emit a log-like
quantity. A zero residual there means the two effects combine *multiplicatively*
in count units — not that there is "no interaction". The `Residual on` control
back-transforms to linear so you can test additivity in counts instead. Every
adapter declares its own `output_space`, and each result records which scale was
used.

**A non-zero I is a statement about the model, not about biology.** It says this
model does not combine these two elements additively on the stated scale. It is
not evidence of cooperativity and not evidence that two factors touch. Profiles
merely *looking* different between constructs establishes less still, which is
why all four constructs are differenced rather than comparing A+B to A.

**Three built-in controls:**

- **Puffin (mock) is additive by construction.** Its residual is identically
  zero. Keep it in the comparison: whatever I *it* reports is pure artefact, and
  that is the floor a real interaction claim has to clear.
- **PuffinD (mock, pair terms off)** keeps only the saturating nonlinearity.
  Comparing against it separates "the response saturates" from "the model has a
  genuine pair term" — a saturating response alone produces non-zero I.
- **Junction and overlap detection.** Overlapping sites mean the A+B construct
  never contained two intact sites. Even non-overlapping sites closer than the
  widest motif (11 bp) let scan windows span the junction and see sequence
  present in neither single construct. Both are flagged, and overlapping
  spacings are dropped from spacing curves.

That last one is not hypothetical. Two GC-boxes 2 bp apart give additive-Puffin
Σ|I| = 74 and PuffinD Σ|I| = 61 — the "interaction" is *smaller* than the pure
sequence artefact. Read at face value it would look like strong cooperativity.

**Tau caveats.** Tau depends on the cell-type panel, so values from different
panels are not comparable. It depends on the activity scale and on how a profile
is collapsed to a scalar (`mean`/`max`/`sum`). Both are recorded with every
result.

**One background proves nothing.** Every experiment takes a replicate count and
reports the spread. Treat an effect that appears in one background and not
others as a property of that background.

## Exports

Design: FASTA, design JSON (settings + full evaluation), profiles CSV.
Mode 1: per-position CSV (activities, deltas, tau, mechanism), JSON.
Mode 2: pair/grid/spacing CSV, JSON, and the four matched constructs as FASTA.
`Save on server` writes a timestamped bundle to `data/runs/` including the full
model roster, so a result stays traceable to the adapter that produced it.

## Layout

```
backend/
  app.py            stdlib HTTP server, JSON API, background jobs
  motifs.py         the 10-motif library, PWM scanning, null-distribution thresholds
  sequence.py       background generation (GC + CpG o/e control), construct assembly
  experiments.py    Mode 1 position scans, Mode 2 pair/grid/spacing
  scoring.py        metrics, tau, specificity decomposition, weighted score
  store.py          players, challenges, submissions, leaderboard
  adapters/
    base.py            ProfileAdapter / CellTypeAdapter interfaces
    mock_common.py     shared mock engine (FFT bump rendering, occupancy, pair terms)
    mock_puffin.py     MockPuffin (additive), MockPuffinD (+saturation +pair terms)
    mock_alphagenome.py  mock CAGE cell-type adapter
    real_puffin.py       stub + checklist
    real_alphagenome.py  stub + checklist
frontend/           vanilla JS, canvas plotting, no dependencies
data/               store.json, saved runs, optional puffin_motifs.npz
```

## Known limitations of this version

- Motif consensus strings and PWMs are **placeholders**, not Puffin's learned
  motifs. Drop real ones into `data/puffin_motifs.npz` to replace them. The
  Initiator is written in a more informative form than the minimal `YYANWYY`,
  which at ~6 bits matches random sequence every ~130 bp and swamped the
  profile; this is noted in the motif table.
- The mock minus-strand track is a strand-attenuated echo at the same position
  rather than a properly offset divergent peak. Place motifs on the minus strand
  to get genuinely separated bidirectional signal.
- `Σ|I|` accumulates across the prediction window, so it grows with window
  width. Only compare it between runs that used the same window.
- The store is a JSON file with a lock — fine for a lab-sized group on one node,
  not for real concurrency. Swap in SQLite if it needs to outlive a session.
