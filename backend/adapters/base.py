"""Model adapter interface.

Everything the app knows about a prediction model goes through one of two
interfaces:

``ProfileAdapter``
    Sequence -> a transcription-initiation profile along the sequence.
    Puffin and PuffinD are both of this kind.  Used by the design view and by
    Mode 2 (motif interactions).

``CellTypeAdapter``
    Sequence -> one profile *per cell type / track*.  AlphaGenome (CAGE output)
    is of this kind.  Used by Mode 1 (cell-type specificity).

Every adapter must declare:

``is_mock``
    True when the numbers are synthetic.  The UI shows a loud badge whenever
    this is True; nothing in the app is allowed to hide it.

``scale``
    A human-readable statement of what the returned numbers *are*.  This
    matters for Mode 2: an interaction residual is only interpretable once you
    know whether the model output is linear, log, or something else.  A
    residual computed on a log scale tests for deviation from *multiplicative*
    combination; on a linear scale it tests for deviation from *additive*
    combination.

``output_space``
    One of ``"linear"`` or ``"log"``.  The Mode 2 code uses this to label the
    residual correctly and to offer the matching transform.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class ProfilePrediction:
    """A profile prediction along a construct."""

    model: str
    positions: np.ndarray          # relative to TSS, shape (P,)
    tracks: dict[str, np.ndarray]  # track name -> shape (P,)
    scale: str
    output_space: str              # "linear" | "log"
    is_mock: bool
    meta: dict = field(default_factory=dict)

    def to_json(self, round_to: int = 5) -> dict:
        return {
            "model": self.model,
            "positions": self.positions.astype(int).tolist(),
            "tracks": {k: np.round(np.asarray(v, dtype=float), round_to).tolist()
                       for k, v in self.tracks.items()},
            "scale": self.scale,
            "output_space": self.output_space,
            "is_mock": self.is_mock,
            "meta": self.meta,
        }


@dataclass
class CellTypePrediction:
    """Per-cell-type profiles (e.g. AlphaGenome CAGE tracks)."""

    model: str
    positions: np.ndarray                 # relative to TSS, shape (P,)
    cell_types: list[str]
    profiles: np.ndarray                  # shape (n_cell_types, P)
    scale: str
    output_space: str
    is_mock: bool
    meta: dict = field(default_factory=dict)

    def to_json(self, round_to: int = 5) -> dict:
        return {
            "model": self.model,
            "positions": self.positions.astype(int).tolist(),
            "cell_types": list(self.cell_types),
            "profiles": np.round(self.profiles, round_to).tolist(),
            "scale": self.scale,
            "output_space": self.output_space,
            "is_mock": self.is_mock,
            "meta": self.meta,
        }


class ProfileAdapter:
    """Sequence -> initiation profile."""

    name: str = "abstract"
    label: str = "Abstract profile model"
    is_mock: bool = True
    scale: str = "undefined"
    output_space: str = "linear"
    description: str = ""

    def available(self) -> tuple[bool, str]:
        """(usable?, reason if not)."""
        return True, ""

    def predict(self, sequence: str, tss_index: int,
                window: tuple[int, int] = (-1000, 1000)) -> ProfilePrediction:
        raise NotImplementedError

    def predict_batch(self, items: list[tuple[str, int]],
                      window: tuple[int, int] = (-1000, 1000)
                      ) -> list[ProfilePrediction]:
        return [self.predict(s, t, window) for s, t in items]

    def info(self) -> dict:
        ok, reason = self.available()
        return {
            "name": self.name, "label": self.label, "kind": "profile",
            "is_mock": self.is_mock, "scale": self.scale,
            "output_space": self.output_space, "description": self.description,
            "available": ok, "unavailable_reason": reason,
        }


class CellTypeAdapter:
    """Sequence -> per-cell-type profiles."""

    name: str = "abstract"
    label: str = "Abstract cell-type model"
    is_mock: bool = True
    scale: str = "undefined"
    output_space: str = "linear"
    description: str = ""
    default_cell_types: list[str] = []

    def available(self) -> tuple[bool, str]:
        return True, ""

    def list_cell_types(self) -> list[dict]:
        return [{"id": c, "label": c} for c in self.default_cell_types]

    def predict(self, sequence: str, tss_index: int,
                cell_types: list[str] | None = None,
                window: tuple[int, int] = (-1000, 1000)) -> CellTypePrediction:
        raise NotImplementedError

    def info(self) -> dict:
        ok, reason = self.available()
        return {
            "name": self.name, "label": self.label, "kind": "celltype",
            "is_mock": self.is_mock, "scale": self.scale,
            "output_space": self.output_space, "description": self.description,
            "available": ok, "unavailable_reason": reason,
            "n_cell_types": len(self.default_cell_types),
        }


def window_positions(window: tuple[int, int]) -> np.ndarray:
    return np.arange(int(window[0]), int(window[1]) + 1, dtype=np.int64)


def slice_to_window(sequence: str, tss_index: int,
                    window: tuple[int, int]) -> tuple[str, np.ndarray, int]:
    """Return the sub-sequence covering ``window`` plus its positions.

    Pads with 'N' if the window runs past the ends of the construct.
    """
    positions = window_positions(window)
    start = tss_index + int(window[0])
    end = tss_index + int(window[1]) + 1
    left_pad = max(0, -start)
    right_pad = max(0, end - len(sequence))
    sub = "N" * left_pad + sequence[max(0, start):min(len(sequence), end)] + "N" * right_pad
    return sub, positions, left_pad
