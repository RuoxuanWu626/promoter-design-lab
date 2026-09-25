"""MOCK cell-type adapter, standing in for AlphaGenome CAGE output.

The real adapter (``adapters/real_alphagenome.py``) runs the human AlphaGenome
model and reads its CAGE tracks. This mock has the same interface and output
shape so Mode 1 can be used without a GPU.

How specificity is generated, and why it changed
------------------------------------------------
An earlier version of this mock gave each cell type its own weights over
Puffin's ten motifs, so cell-type specificity came out of the core promoter.
That is the wrong story. Puffin's ten motifs are core promoter elements — TATA,
the initiator, the GC-box, CCAAT, NRF1, ETS, CREB, YY1, ZNF143, the U1 site.
They set how much initiation a promoter produces and where it starts, and they
do that in much the same way in every cell type. They are not what makes a
promoter cell-type specific.

So this version separates the two:

* **Core promoter motifs** (``motifs.library_order()``) are weighted almost
  identically across cell types — only a small jitter. They carry the bulk of
  the *activity*, and they are sharply position-dependent, which is what Mode
  1's position scan probes.

* **Lineage TF sites** (``celltype_elements``) are weighted strongly and
  differently per cell type: GATA in erythroid, HNF4A and C/EBP in hepatocyte,
  PU.1 and EBF1 in lymphoblastoid, SOX2–OCT4 in ES, and so on. They carry the
  *specificity*, and their position dependence is broad rather than sharp,
  because enhancer-like sites are not anchored to the TSS the way a TATA box
  is.

Mode 1's attribution panel measures exactly this split, so the mock now has to
be honest about it or the panel would teach the wrong lesson.
"""

from __future__ import annotations

import numpy as np

from motifs import celltype_element_order, get_library, library_order
from adapters.base import CellTypeAdapter, CellTypePrediction
from adapters.mock_common import (
    add_bump,
    composition_baseline,
    hit_arrays_geometry,
    motif_hits_by_motif,
    render_bumps,
    saturate,
)

MOCK_SCALE = "log10(mock CAGE CPM + 1), arbitrary units - NOT a real model"

# Cell types, their baseline, and how sharply they read core promoter motifs.
# There are no per-cell-type motif preferences here any more: that is the whole
# point (see the module docstring).
#      id,          label,                        dist_opt, dist_sigma, bias
_CELL_TYPE_SPEC = [
    ("K562",       "K562 (erythroleukemia)",          46, 34, 0.55),
    ("HepG2",      "HepG2 (hepatocyte)",              72, 40, 0.60),
    ("GM12878",    "GM12878 (lymphoblastoid)",        24, 22, 0.50),
    ("HeLa-S3",    "HeLa-S3 (cervical)",              58, 44, 0.55),
    ("A549",       "A549 (lung)",                     66, 38, 0.55),
    ("MCF-7",      "MCF-7 (breast)",                  52, 36, 0.60),
    ("SK-N-SH",    "SK-N-SH (neuroblastoma)",         30, 20, 0.45),
    ("IMR-90",     "IMR-90 (fibroblast)",             88, 48, 0.50),
    ("H1",         "H1 (embryonic stem)",             40, 30, 0.55),
    ("HCT116",     "HCT116 (colorectal)",             62, 42, 0.55),
    ("HUVEC",      "HUVEC (endothelial)",             44, 28, 0.50),
    ("Hepatocyte", "Primary hepatocyte",              34, 26, 0.50),
]

# How much cell types are allowed to differ in how they read a CORE promoter
# motif. Small on purpose: core promoter strength is largely shared.
CORE_WEIGHT = 1.0
CORE_JITTER = 0.06

# Cell types share a baseline. Giving each one its own bias and CpG preference
# would put tau well above zero for an empty promoter, which would read as
# "specificity" that no sequence encodes -- exactly the confusion Mode 1 is
# supposed to remove. Library depth does differ in reality, but that is a
# property of the assay, not of the promoter being designed.
SHARED_BIAS = 0.5
SHARED_CPG_PREF = 1.0

# Lineage sites act over a broad window rather than from a fixed distance.
ELEMENT_DIST_SIGMA = 260.0
ELEMENT_DIST_FLOOR = 0.55

# Lineage sites carry a small amplitude in the shared motif library, because
# they should barely perturb a *profile* model -- Puffin has no filter for
# them. The cell-type model is the one that reads them, so it applies its own
# gain here. Without this the core promoter saturates the track and a lineage
# site cannot move tau at all.
ELEMENT_GAIN = float(__import__("os").environ.get("PDG_ELEMENT_GAIN", 3.0))


def _build_params(seed: int = 20240917) -> dict[str, dict]:
    from celltype_elements import activation_map

    rng = np.random.default_rng(seed)
    core = library_order()
    elements = celltype_element_order()
    act_map = activation_map()

    params: dict[str, dict] = {}
    for cid, label, dopt, dsig, bias in _CELL_TYPE_SPEC:
        # Core promoter motifs: nearly the same everywhere.
        w_core = {m: float(np.round(CORE_WEIGHT + rng.uniform(-CORE_JITTER, CORE_JITTER), 3))
                  for m in core}
        # Lineage sites: strongly cell-type specific.
        w_elem = {}
        for e in elements:
            v = act_map.get(e, {}).get(cid)
            if v is None:
                v = float(np.round(rng.uniform(-0.05, 0.10), 3))
            w_elem[e] = float(v)

        params[cid] = {
            "id": cid, "label": label,
            "core_weights": w_core,
            "element_weights": w_elem,
            "responds_to": [e for e in elements if w_elem[e] >= 0.5],
            "repressed_by": [e for e in elements if w_elem[e] <= -0.3],
            # Cell types still differ in the distance they prefer to read a
            # core promoter motif from -- that is what Mode 1's position scan
            # maps -- but the gate never drops far, so position modulates the
            # response rather than switching cell types on and off.
            "dist_opt": float(dopt), "dist_sigma": float(dsig),
            "dist_floor": float(np.round(rng.uniform(0.70, 0.82), 3)),
            "cpg_pref": SHARED_CPG_PREF,
            "bias": SHARED_BIAS,
            "sharpness": float(np.round(rng.uniform(0.9, 1.15), 3)),
        }
    return params


_PARAMS = _build_params()


class MockAlphaGenome(CellTypeAdapter):
    name = "mock_alphagenome"
    label = "AlphaGenome CAGE (mock)"
    is_mock = True
    scale = MOCK_SCALE
    output_space = "log"
    description = (
        "Placeholder stand-in for AlphaGenome's human CAGE output. Core "
        "promoter motifs are weighted almost identically across cell types and "
        "carry the activity; lineage TF sites are weighted per cell type and "
        "carry the specificity. Replace with the real adapter for actual "
        "predictions."
    )
    default_cell_types = [c[0] for c in _CELL_TYPE_SPEC]

    def list_cell_types(self) -> list[dict]:
        return [{
            "id": p["id"],
            "label": p["label"] + " [mock]",
            "responds_to": p["responds_to"],
            "repressed_by": p["repressed_by"],
            "dist_opt": p["dist_opt"],
            "mock": True,
        } for p in _PARAMS.values()]

    def _profile_for(self, cid: str, sequence: str, tss_index: int,
                     base: np.ndarray, by_motif: dict,
                     extra_weights: dict | None = None) -> np.ndarray:
        """One cell type's CAGE-like profile from a pre-computed scan."""
        p = _PARAMS[cid]
        lib = get_library()
        L = len(sequence)
        track = p["bias"] + p["cpg_pref"] * base

        for mid, d in by_motif.items():
            if d["index"].size == 0:
                continue
            motif = lib.get(mid)
            if motif is None:
                continue

            is_element = getattr(motif, "kind", "core_promoter") == "celltype_element"
            if is_element:
                w = p["element_weights"].get(mid)
                if w is None:
                    # A lineage site the user defined after this adapter was
                    # built. Its cell-type assignments are read live.
                    w = (extra_weights or {}).get(mid, {}).get(cid, 0.05)
            else:
                w = p["core_weights"].get(mid, CORE_WEIGHT)
            if not w:
                continue

            centers, amps = hit_arrays_geometry(d, motif, "plus")
            dist = np.abs(d["index"] + motif.width / 2.0 - tss_index)

            if is_element:
                # Broad, weakly position-dependent: a lineage site works from
                # anywhere in the promoter region, unlike a TATA box.
                gate = ELEMENT_DIST_FLOOR + (1.0 - ELEMENT_DIST_FLOOR) * np.exp(
                    -0.5 * (dist / ELEMENT_DIST_SIGMA) ** 2)
            else:
                # Sharp: core promoter motifs are read from a preferred
                # distance, which is what the Mode 1 position scan maps.
                gate = p["dist_floor"] + (1.0 - p["dist_floor"]) * np.exp(
                    -0.5 * ((dist - p["dist_opt"]) / p["dist_sigma"]) ** 2)

            if is_element:
                amps = amps * ELEMENT_GAIN

            track += render_bumps(L, centers, amps * w * gate,
                                  motif.mock_width * p["sharpness"],
                                  motif.mock_broad_fraction)

        # A CAGE track cannot go below zero, and repressive sites can push it
        # there. Clipping is a nonlinearity, but nothing in Mode 1 depends on
        # this adapter being additive.
        return np.clip(saturate(track, ceiling=3.0, knee=0.9), 0.0, None)

    def predict(self, sequence: str, tss_index: int,
                cell_types: list[str] | None = None,
                window: tuple[int, int] = (-1000, 1000)) -> CellTypePrediction:
        cts = [c for c in (cell_types or self.default_cell_types) if c in _PARAMS]
        if not cts:
            cts = list(self.default_cell_types)

        positions = np.arange(int(window[0]), int(window[1]) + 1, dtype=np.int64)
        idx = positions + tss_index
        valid = (idx >= 0) & (idx < len(sequence))

        base = composition_baseline(sequence)
        by_motif = motif_hits_by_motif(sequence)
        try:
            from custom_motifs import activation_map as custom_activation_map
            extra = custom_activation_map()
        except Exception:
            extra = {}

        profiles = np.zeros((len(cts), positions.size), dtype=np.float64)
        for i, cid in enumerate(cts):
            full = self._profile_for(cid, sequence, tss_index, base, by_motif, extra)
            profiles[i, valid] = full[idx[valid]]

        return CellTypePrediction(
            model=self.name, positions=positions, cell_types=cts,
            profiles=profiles, scale=self.scale,
            output_space=self.output_space, is_mock=True,
            meta={"cell_type_labels": {c: _PARAMS[c]["label"] for c in cts}},
        )
