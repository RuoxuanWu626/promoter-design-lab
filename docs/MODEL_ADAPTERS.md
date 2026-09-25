# Wiring up real inference

Adapters are the only part of the app that knows a model exists. Implement one,
register it, and it appears in the model pickers with a `real` badge. Nothing in
the UI, the experiments, the scoring or the export path needs to change.

## The two interfaces

Both live in `backend/adapters/base.py`.

`ProfileAdapter` — sequence to an initiation profile. Puffin and PuffinD are
this kind; used by the design view and all of Mode 2.

```python
predict(sequence: str, tss_index: int, window: (int, int)) -> ProfilePrediction
```

`CellTypeAdapter` — sequence to one profile per cell type / track. AlphaGenome
CAGE is this kind; used by Mode 1.

```python
predict(sequence, tss_index, cell_types, window) -> CellTypePrediction
```

Both must also implement `available() -> (bool, reason)`. Return `False` with a
useful reason when weights, an API key or a dependency is missing — the UI shows
the reason in the picker and on the Models page instead of failing at click
time.

## Registering

```python
# backend/adapters/__init__.py
from adapters.real_puffin import RealPuffin
_PROFILE[RealPuffin().name] = RealPuffin()
```

The stubs `real_puffin.py` and `real_alphagenome.py` are already registered and
report themselves as unavailable until implemented.

## Four declarations that change what results mean

These are not bookkeeping. Each one silently changes the interpretation of the
plots without changing how they look.

**`output_space`** — `"linear"` or `"log"`. Mode 2's residual
`I = f(A+B) − f(A) − f(B) + f(bg)` is a test of departure from additivity *on
the scale the model emits*. On a log-like output, `I = 0` means the two effects
combine multiplicatively in count units. Declare this wrong and every Mode 2
conclusion inverts.

**`scale`** — a human-readable statement of what the numbers are, shown under
every plot and written into every export. `"UNVERIFIED"` is the honest value
until you have checked against the checkpoint; leave it until you have.

**`is_mock`** — drives the badge. Nothing in the app is permitted to hide it.

**Output resolution** — if the model is binned (AlphaGenome CAGE is), the
`positions` you return must be bin centres relative to the TSS, not per-base
indices, and the bin size belongs in `meta`. Mode 1's position scan will
otherwise be read at a finer resolution than the model actually has.

## Puffin / PuffinD checklist

1. **Input length.** Puffin expects a fixed receptive field. Either pad/crop the
   construct around the TSS or predict on a sliding window and stitch. Record
   which in `meta` so exports say what was done.
2. **Strandedness.** Return both tracks under the `plus` and `minus` keys.
3. **Determinism.** `model.eval()`, no dropout. The four Mode 2 constructs are
   only comparable if repeated inference is bit-identical — and the prediction
   cache assumes it.
4. **Batching.** A 10 × 10 grid over 3 backgrounds is ~2600 forward passes.
   Implement `predict_batch` and it will be used.

Environment: `PUFFIN_REPO`, `PUFFIN_WEIGHTS`, `PUFFIND_WEIGHTS`.

## AlphaGenome checklist

The `alphagenome` env on this cluster already has the client library:

```
/gpfs/data/zhou-lab/rxwu/settings/miniforge3/envs/alphagenome/bin/python
```

Start the server with that interpreter and set `ALPHAGENOME_API_KEY`:

```bash
PDG_PYTHON=/gpfs/data/zhou-lab/rxwu/settings/miniforge3/envs/alphagenome/bin/python \
  ./run_server.sh
```

1. **Enumerate tracks from the model's own output metadata**, not a hand-written
   list, and keep biosample name plus ontology id so results stay traceable.
2. **Input length.** Pad the construct to a supported interval size, centred so
   the designed TSS sits where you intend. Record the padding in `meta`.
3. **Scale.** Record whether values are raw predicted counts or transformed.
   Tau is computed on a scalar per cell type; taking that scalar on the raw
   versus log scale gives different numbers, so state which.
4. **Quota and latency.** A Mode 1 scan is one request per offset — a default
   scan is 81 offsets × replicate backgrounds. Batch where the API allows,
   cache by sequence hash (`adapters.cache_key` / `adapters.cached` already do
   this), and warn before launching a large scan.

## Replacing the placeholder motifs

Save a `.npz` to `data/puffin_motifs.npz` with one `(4, W)` array per motif id,
rows ordered `A, C, G, T`, columns 5'→3'. Recognised ids:

```
tata  inr  sp1  nfy  nrf1  ets  creb  yy1  znf143  u1
```

They are loaded in preference to the built-ins and the consensus strings are
regenerated from them. Occupancy thresholds are set per motif at a fixed
false-positive rate under a realistic background (`motifs.null_threshold`), so
PWMs with different information content are handled automatically — no retuning.

## Two invariants to preserve

**Do not make a motif's contribution depend on what else is in the sequence.**
An earlier revision capped the hit list at "top 60 matches per motif". That
makes one site's effect depend on the rest of the construct, which quietly
destroyed the additivity Mode 2 uses as its control — Puffin's residual should
have been 0 and was 15.5. Fixed thresholds are safe; rank cuts and
whole-sequence renormalisation are not.

**Keep an additive model registered.** `MockPuffin` is additive by construction,
so any non-zero residual it reports is pure sequence-composition artefact from
the junction between two sites. It is the cheapest available control on the
whole Mode 2 pipeline, and it stays useful after real models are wired up.
