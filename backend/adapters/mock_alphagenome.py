"""MOCK cell-type adapter, standing in for AlphaGenome CAGE-track output.

The real adapter (``adapters/real_alphagenome.py``) calls the AlphaGenome human
model and pulls CAGE tracks.  This mock has the same interface and the same
output shape so Mode 1 can be built and exercised without an API key.

How the mock creates cell-type specificity
------------------------------------------
Each mock cell type has

* a weight per motif: two "driver" motifs it responds to strongly, one or two
  it ignores, and one that *represses* it (a negative weight, standing in for a
  site that an activator reads in one lineage and a repressor reads in
  another),
* a preferred *distance from the TSS*, ``dist_opt``, with a tolerance
  ``dist_sigma``, applied to every motif occurrence, and
* a baseline responsiveness to local CpG/GC content.

The distance term is what makes Mode 1 non-trivial: sliding the same motif
across positions changes which cell types respond, because cell types differ in
where they "read" a motif from.  This is a caricature, not a claim about real
promoters - but it produces exactly the kind of position x cell-type structure
the mode is built to explore.

The repressive weights matter for a specific reason.  Mode 1 asks the user to
tell apart specificity won by *raising the target* from specificity won by
*lowering everything else*.  If every motif were activating in every cell type,
the second route would not exist and the distinction could never be
demonstrated.  So some motif/cell-type combinations push activity down.
"""

from __future__ import annotations

import numpy as np

from motifs import get_library, library_order
from adapters.base import CellTypeAdapter, CellTypePrediction
from adapters.mock_common import (
    composition_baseline,
    hit_arrays_geometry,
    motif_hits_by_motif,
    render_bumps,
    saturate,
)

MOCK_SCALE = "log10(mock CAGE CPM + 1), arbitrary units - NOT a real model"

# Cell-type ids echo commonly used CAGE cell lines so the UI reads naturally.
# Every number attached to them here is invented.
_CELL_TYPE_SPEC = [
    # The U1 snRNP site is repressive in 8 of the 12 cell types and mildly
    # activating in the other 4. That asymmetry is deliberate: it gives the
    # player an "off-target silencer" to discover, so the difference between
    # winning specificity by raising the target and winning it by lowering
    # everything else is something Mode 1 can actually show.
    #
    # id,        label,                  drivers,            ignores,     repressors,        dist_opt, dist_sigma, bias
    ("K562",     "K562 (erythroleukemia)",   ["ets", "sp1"],     ["tata"],    ["znf143", "u1"],   46,   34,  0.55),
    ("HepG2",    "HepG2 (hepatocyte)",       ["creb", "nfy"],    [],          ["yy1", "u1"],      72,   40,  0.60),
    ("GM12878",  "GM12878 (lymphoblastoid)", ["ets", "yy1"],     ["nfy"],     ["creb", "u1"],     24,   22,  0.50),
    ("HeLa-S3",  "HeLa-S3 (cervical)",       ["sp1", "nrf1"],    ["tata"],    ["u1"],             58,   44,  0.55),
    ("A549",     "A549 (lung)",              ["nfy", "creb"],    ["znf143"],  ["inr", "u1"],      66,   38,  0.55),
    ("MCF-7",    "MCF-7 (breast)",           ["sp1", "creb"],    [],          ["tata", "u1"],     52,   36,  0.60),
    ("SK-N-SH",  "SK-N-SH (neuroblastoma)",  ["tata", "inr"],    ["sp1"],     ["nrf1"],           30,   20,  0.45),
    ("IMR-90",   "IMR-90 (fibroblast)",      ["nfy", "znf143"],  ["inr"],     ["ets", "u1"],      88,   48,  0.50),
    ("H1",       "H1 (embryonic stem)",      ["nrf1", "yy1"],    ["tata"],    ["nfy"],            40,   30,  0.55),
    ("HCT116",   "HCT116 (colorectal)",      ["znf143", "ets"],  ["yy1"],     ["sp1", "u1"],      62,   42,  0.55),
    ("HUVEC",    "HUVEC (endothelial)",      ["ets", "creb"],    ["tata"],    ["nrf1", "u1"],     44,   28,  0.50),
    ("Hepatocyte", "Primary hepatocyte",     ["creb", "tata"],   ["sp1"],     ["ets"],            34,   26,  0.50),
]


def _build_params(seed: int = 20240917) -> dict[str, dict]:
    rng = np.random.default_rng(seed)
    order = library_order()
    params: dict[str, dict] = {}
    for cid, label, drivers, ignores, repressors, dopt, dsig, bias in _CELL_TYPE_SPEC:
        w = {m: float(np.round(rng.uniform(0.30, 0.85), 3)) for m in order}
        for m in drivers:
            w[m] = float(np.round(rng.uniform(1.55, 2.20), 3))
        for m in ignores:
            w[m] = float(np.round(rng.uniform(0.02, 0.12), 3))
        for m in repressors:
            w[m] = float(np.round(rng.uniform(-1.30, -0.65), 3))
        params[cid] = {
            "id": cid,
            "label": label,
            "drivers": drivers,
            "ignores": ignores,
            "repressors": repressors,
            "weights": w,
            "dist_opt": float(dopt),
            "dist_sigma": float(dsig),
            "dist_floor": float(np.round(rng.uniform(0.12, 0.32), 3)),
            "cpg_pref": float(np.round(rng.uniform(0.55, 1.45), 3)),
            "bias": float(bias),
            "sharpness": float(np.round(rng.uniform(0.85, 1.25), 3)),
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
        "Placeholder stand-in for AlphaGenome's human model restricted to CAGE "
        "tracks. Each mock cell type weights the 10 motifs differently AND has "
        "its own preferred distance from the TSS, so both motif identity and "
        "motif position change the cell-type profile. Replace with "
        "real_alphagenome.RealAlphaGenome for actual predictions."
    )
    default_cell_types = [c[0] for c in _CELL_TYPE_SPEC]

    def list_cell_types(self) -> list[dict]:
        return [{
            "id": p["id"],
            "label": p["label"] + " [mock]",
            "drivers": p["drivers"],
            "ignores": p["ignores"],
            "repressors": p["repressors"],
            "dist_opt": p["dist_opt"],
            "mock": True,
        } for p in _PARAMS.values()]

    def _profile_for(self, cid: str, sequence: str, tss_index: int,
                     base: np.ndarray, by_motif: dict) -> np.ndarray:
        """One cell type's CAGE-like profile, from a pre-computed scan.

        ``base`` and ``by_motif`` are shared across every cell type in a
        request: the sequence is scanned once and each cell type only
        re-weights the result. That is what keeps a Mode 1 position scan over
        a dozen cell types affordable.
        """
        p = _PARAMS[cid]
        lib = get_library()
        L = len(sequence)
        track = p["bias"] + p["cpg_pref"] * base

        for mid, d in by_motif.items():
            if d["index"].size == 0:
                continue
            w = p["weights"].get(mid, 0.4)
            if w == 0:
                continue
            motif = lib[mid]
            centers, amps = hit_arrays_geometry(d, motif, "plus")

            # How far each occurrence sits from the TSS. Each cell type has
            # its own preferred distance, so moving a motif changes cell types
            # differently -- the effect Mode 1 is built to map.
            dist = np.abs(d["index"] + motif.width / 2.0 - tss_index)
            gate = p["dist_floor"] + (1.0 - p["dist_floor"]) * np.exp(
                -0.5 * ((dist - p["dist_opt"]) / p["dist_sigma"]) ** 2)

            track += render_bumps(L, centers, amps * w * gate,
                                  motif.mock_width * p["sharpness"],
                                  motif.mock_broad_fraction)

        # A CAGE track cannot go below zero, and repressive weights can push it
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

        profiles = np.zeros((len(cts), positions.size), dtype=np.float64)
        for i, cid in enumerate(cts):
            full = self._profile_for(cid, sequence, tss_index, base, by_motif)
            profiles[i, valid] = full[idx[valid]]

        return CellTypePrediction(
            model=self.name, positions=positions, cell_types=cts,
            profiles=profiles, scale=self.scale,
            output_space=self.output_space, is_mock=True,
            meta={"cell_type_labels": {c: _PARAMS[c]["label"] for c in cts}},
        )
