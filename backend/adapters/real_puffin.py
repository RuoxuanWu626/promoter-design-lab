"""Real Puffin / PuffinD adapters - NOT YET WIRED UP.

This file is the single place to plug in actual Puffin inference.  It keeps the
same ``ProfileAdapter`` contract as the mock, so nothing else in the app has to
change: register the class in ``adapters/__init__.py`` and it appears in the
model picker.

What you need to fill in
------------------------
1. ``PUFFIN_REPO`` / ``PUFFIN_WEIGHTS`` - point at a Puffin checkout and its
   trained weights (env vars ``PUFFIN_REPO`` and ``PUFFIN_WEIGHTS`` override).
2. ``_load_model`` - import the model class and load the checkpoint.
3. ``predict`` - one-hot encode, run the model, and return a
   ``ProfilePrediction``.

Things to get right when you do
-------------------------------
* **Input length.** Puffin expects a fixed receptive field. Either pad/crop the
  construct to that length around the TSS, or predict on a sliding window and
  stitch. Whatever you choose, record it in ``meta`` so exported results say
  what was done.
* **Output space.** Set ``output_space`` and ``scale`` to what the model
  actually emits. Mode 2's interaction residual is only interpretable relative
  to that scale - get this wrong and the residual means something else.
* **Strandedness.** Puffin emits separate plus/minus initiation tracks; return
  both under the ``plus`` / ``minus`` keys.
* **Determinism.** Put the model in eval mode and disable dropout, otherwise
  the four matched Mode 2 constructs are not comparable.
"""

from __future__ import annotations

import os

import numpy as np

from adapters.base import ProfileAdapter, ProfilePrediction

PUFFIN_REPO = os.environ.get("PUFFIN_REPO", "")
PUFFIN_WEIGHTS = os.environ.get("PUFFIN_WEIGHTS", "")


class RealPuffinBase(ProfileAdapter):
    is_mock = False
    output_space = "log"          # verify against the checkpoint before trusting
    scale = "UNVERIFIED - set this from the checkpoint before using results"
    _weights_env = "PUFFIN_WEIGHTS"
    _model = None

    def available(self) -> tuple[bool, str]:
        try:
            import torch  # noqa: F401
        except ImportError:
            return False, ("PyTorch is not installed in this environment. "
                           "Create an env with torch and restart the server "
                           "with that interpreter.")
        weights = os.environ.get(self._weights_env, "")
        if not weights:
            return False, f"Set ${self._weights_env} to the checkpoint path."
        if not os.path.exists(weights):
            return False, f"Checkpoint not found: {weights}"
        return True, ""

    def _load_model(self):
        raise NotImplementedError(
            "Fill in _load_model(): import the Puffin model class from "
            f"{PUFFIN_REPO or '$PUFFIN_REPO'} and load "
            f"{os.environ.get(self._weights_env) or '$' + self._weights_env}."
        )

    def predict(self, sequence: str, tss_index: int,
                window: tuple[int, int] = (-1000, 1000)) -> ProfilePrediction:
        ok, reason = self.available()
        raise NotImplementedError(
            f"{self.label} is not wired up yet. {reason or ''} "
            "Implement predict() in backend/adapters/real_puffin.py."
        )


class RealPuffin(RealPuffinBase):
    name = "puffin"
    label = "Puffin (real)"
    description = ("Real Puffin inference. Not implemented yet - see "
                   "backend/adapters/real_puffin.py for the checklist.")
    _weights_env = "PUFFIN_WEIGHTS"


class RealPuffinD(RealPuffinBase):
    name = "puffind"
    label = "PuffinD (real)"
    description = ("Real PuffinD inference. Not implemented yet - see "
                   "backend/adapters/real_puffin.py for the checklist.")
    _weights_env = "PUFFIND_WEIGHTS"
