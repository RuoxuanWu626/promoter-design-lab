"""Real AlphaGenome CAGE adapter - NOT YET WIRED UP.

Mode 1 is designed around a cell-type-aware model; the intended one is the
AlphaGenome human model, reading CAGE output tracks.  This file is where that
gets plugged in.  The environment ``alphagenome`` on this cluster already has
the client library installed:

    /gpfs/data/zhou-lab/rxwu/settings/miniforge3/envs/alphagenome/bin/python

Start the server with that interpreter and set ``ALPHAGENOME_API_KEY`` to make
this adapter selectable.

What you need to fill in
------------------------
1. ``_client`` - construct the AlphaGenome client from the API key.
2. ``list_cell_types`` - enumerate the CAGE tracks you want exposed, from the
   model's own output metadata rather than a hand-written list, and keep the
   ontology / biosample ids so exported results are traceable.
3. ``predict`` - request the CAGE output type for the construct, select the
   chosen tracks, and return a ``CellTypePrediction``.

Things to get right when you do
-------------------------------
* **Input length.** AlphaGenome takes fixed-size input intervals. Pad the
  construct to the nearest supported length, centred so the designed TSS sits
  where you intend, and record the padding in ``meta``.
* **Resolution.** CAGE output is binned. The returned ``positions`` must be the
  bin centres relative to the TSS, not per-base indices, and the UI should be
  told the bin size so Mode 1's position scan is not read at a finer resolution
  than the model actually has.
* **Scale.** Record whether the values are raw predicted counts or transformed.
  Tau is computed on a scalar activity per cell type; whether you take that
  scalar on the raw or log scale changes the number, so state it.
* **Quota and latency.** Each Mode 1 position scan is one request per offset.
  Batch where the API allows it, cache aggressively (``adapters.cache_key``),
  and warn the user before launching a scan with hundreds of offsets.
* **Determinism.** Cache by sequence hash so a repeated scan does not re-bill
  and does not drift.
"""

from __future__ import annotations

import os

from adapters.base import CellTypeAdapter, CellTypePrediction

API_KEY_ENV = "ALPHAGENOME_API_KEY"


class RealAlphaGenome(CellTypeAdapter):
    name = "alphagenome"
    label = "AlphaGenome CAGE (real)"
    is_mock = False
    output_space = "linear"   # verify against the API response before trusting
    scale = "UNVERIFIED - set from the AlphaGenome output metadata"
    description = ("Real AlphaGenome human model, CAGE output tracks. Not "
                   "implemented yet - see backend/adapters/real_alphagenome.py.")
    default_cell_types: list[str] = []

    def available(self) -> tuple[bool, str]:
        try:
            import alphagenome  # noqa: F401
        except ImportError:
            return False, ("The alphagenome package is not importable by this "
                           "interpreter. Start the server with the "
                           "'alphagenome' conda env.")
        if not os.environ.get(API_KEY_ENV):
            return False, f"Set ${API_KEY_ENV} in the server environment."
        return True, ""

    def _client(self):
        raise NotImplementedError("Fill in _client() in real_alphagenome.py.")

    def list_cell_types(self) -> list[dict]:
        ok, reason = self.available()
        if not ok:
            return []
        raise NotImplementedError(
            "Fill in list_cell_types(): enumerate CAGE tracks from the model's "
            "output metadata, keeping biosample name and ontology id."
        )

    def predict(self, sequence: str, tss_index: int,
                cell_types: list[str] | None = None,
                window: tuple[int, int] = (-1000, 1000)) -> CellTypePrediction:
        ok, reason = self.available()
        raise NotImplementedError(
            f"{self.label} is not wired up yet. {reason or ''} "
            "Implement predict() in backend/adapters/real_alphagenome.py."
        )
