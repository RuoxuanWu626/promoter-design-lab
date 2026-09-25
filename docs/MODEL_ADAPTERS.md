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

## Puffin — DONE

`adapters/real_puffin.py` is implemented and is the default profile model.
Architecture in `models/puffin_arch.py`, weights from
`/gpfs/data/zhou-lab/rxwu/puffin_species/human/puffin.pth` (override with
`$PUFFIN_WEIGHTS`). It runs on CPU.

Two things it got right that are easy to get wrong, and that any replacement
must keep:

* **Context, not padding.** Output is valid only 325 bp inside the input.
  Positions without full context are returned as NaN and travel as `null`, so
  the plot gaps. Zero-padding them instead produces an edge artefact far larger
  than any real signal, and it silently poisons every metric downstream
  (argmax lands on the artefact).
* **`deconv` is a cross-correlation.** An activation at *p* drives output at
  *p − argmax(kernel)*, not *p + argmax*. Getting the sign backwards puts every
  motif on the wrong side of the TSS; the check is that TATA must come out at
  −31.

## Older Puffin checkpoint variants

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

## AlphaGenome — what is known, and what blocks it

This is the **local PyTorch port**, not the hosted API. No API key is involved.

```
port     /gpfs/data/zhou-lab/rxwu/Puffin-drozofila/tests/alphagenome_torch
weights  .../alphagenome_all_folds.pt          (1.7 GB)
meta     /gpfs/data/zhou-lab/sharish/genomics/alphagenome_dataset/hg38/bp/
         alphagenome_human_bp_stranded_cp100m_logp1.v3.meta.tsv
```

Loading and inference, taken from the lab's own working benchmark
(`opengenome/benchmark_alphagenome_human_published_cage_specificity.py`):

```python
sys.path.insert(0, TORCH_DIR)
from alphagenome.models import dna_model as dna_model_lib, dna_output
from alphagenome_research.model.metadata import metadata as metadata_lib
from alphagenome_torch import model as torch_model_lib

organism = dna_model_lib.Organism.HOMO_SAPIENS
md = metadata_lib.load(organism)
model = torch_model_lib.AlphaGenomeTorch({organism: md})
# all_folds ships some weights with a leading dim of 2; take [:1]
model.load_state_dict(state, strict=True)
valid = np.flatnonzero(~np.asarray(md.padding[dna_output.OutputType.CAGE]))  # 546

raw, _ = model(x, torch.zeros(1, dtype=torch.long, device=dev), output_heads=("cage",))
p1 = raw["cage"]["predictions_1bp"][0].index_select(-1, valid_t)   # (L, 546)
```

546 CAGE tracks = 273 cell types x 2 strands; the meta TSV has `name` and
`strand` columns to fold them into per-cell-type values.

**The blocker is a CUDA driver mismatch, not the code.** The `alphagenome` env
has `torch 2.12.0+cu130`; the A100 nodes (`gpuq`) run driver 570.124, which is
CUDA 12.8, so `torch.cuda.is_available()` is False there. `gpuq_env` has a
compatible `cu126` build but not the AlphaGenome packages. Any one of these
fixes it:

* a node with driver >= 580 (the H200 in `gpuq`, if free),
* `alphagenome` + `alphagenome_research` installed into a `cu126`/`cu128` env,
* the arm64 env on the GH200 nodes
  (`miniforge3-arm64/envs/alphagenome-arm64`, partitions `ghq` / `pearsonq`),
  which is what the lab's own AlphaGenome jobs use.

1. **Enumerate tracks from the model's own output metadata**, not a hand-written
   list, and keep biosample name plus ontology id so results stay traceable.
2. **Input length.** The model takes 2^20 = 1,048,576 bp. A 1-2 kb construct
   has to be embedded in that context, centred so the designed TSS sits where
   you intend, and the padding recorded in `meta`. What fills the rest is a
   real experimental choice: neutral sequence, a fixed genomic locus, or the
   construct's own background repeated all give different answers, and the one
   you pick belongs in the exported settings.
3. **Scale.** Record whether values are raw predicted counts or transformed.
   Tau is computed on a scalar per cell type; taking that scalar on the raw
   versus log scale gives different numbers, so state which.
4. **Cost.** A Mode 1 position scan is one forward pass per offset (a default
   scan is 81 offsets x replicate backgrounds), and Mode 1's attribution panel
   is one per patch. At 1 Mb per pass that is the dominant cost in the whole
   app. Batch, cache by sequence hash (`adapters.cache_key` / `adapters.cached`
   already do this), and consider running attribution as a submitted job rather
   than interactively.

5. **Specificity comes from the lineage sites.** The mock cell-type adapter
   weights Puffin's core promoter motifs almost identically across cell types
   and puts the cell-type differences on the separate lineage TF library,
   because that is where promoter cell-type specificity actually lives. When
   the real model is wired up, Mode 1's attribution panel is the direct check
   on whether it agrees.

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
