"""MOCK Puffin and MOCK PuffinD profile adapters.

Both return a plus-strand and a minus-strand initiation track.  The minus track
is produced by running the same machinery on the reverse complement and
flipping the result, so divergent initiation upstream of the TSS falls out of
the model rather than being painted on.

MockPuffin  : additive.  f = composition term + sum of per-motif terms.
MockPuffinD : the same drive, then a saturating nonlinearity, plus explicit
              spacing-dependent pair terms (see mock_common.MOCK_PAIR_TERMS).

Reported scale for both is ``log10(mock initiation counts + 1)``: the numbers
are already in a log-like space.  That means the Mode 2 residual computed on
these outputs is a test of departure from a *multiplicative* combination of
effects on the underlying count scale.  Mode 2 lets you switch to the linear
scale (10**x - 1) to test departure from additivity there instead.
"""

from __future__ import annotations

import numpy as np

from adapters.base import ProfileAdapter, ProfilePrediction
from adapters.mock_common import additive_tracks, pair_cooperativity, saturate

MOCK_SCALE = "log10(mock initiation counts + 1), arbitrary units - NOT a real model"


def _window_slice(track: np.ndarray, tss_index: int,
                  window: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    positions = np.arange(int(window[0]), int(window[1]) + 1, dtype=np.int64)
    idx = positions + tss_index
    valid = (idx >= 0) & (idx < len(track))
    out = np.zeros(positions.shape, dtype=np.float64)
    out[valid] = track[idx[valid]]
    return positions, out


class MockPuffin(ProfileAdapter):
    name = "mock_puffin"
    label = "Puffin (mock)"
    is_mock = True
    scale = MOCK_SCALE
    output_space = "log"
    description = (
        "Placeholder stand-in for Puffin. Strictly additive: a term linear in "
        "local dinucleotide composition plus one independent term per motif "
        "occurrence. Because it is additive by construction, its Mode 2 "
        "interaction residual is ~0 whenever the two elements are far enough "
        "apart not to share a PWM scan window. Use it as a negative control."
    )

    def _tracks(self, sequence: str) -> dict[str, np.ndarray]:
        tracks, _by_motif = additive_tracks(sequence)
        return tracks

    def predict(self, sequence: str, tss_index: int,
                window: tuple[int, int] = (-1000, 1000)) -> ProfilePrediction:
        tracks = self._tracks(sequence)
        positions = None
        out = {}
        for k, v in tracks.items():
            positions, out[k] = _window_slice(v, tss_index, window)
        return ProfilePrediction(
            model=self.name, positions=positions, tracks=out,
            scale=self.scale, output_space=self.output_space, is_mock=True,
            meta={"additive_by_construction": True},
        )


class MockPuffinD(ProfileAdapter):
    name = "mock_puffind"
    label = "PuffinD (mock)"
    is_mock = True
    scale = MOCK_SCALE
    output_space = "log"
    description = (
        "Placeholder stand-in for PuffinD, the deep model. Same additive drive "
        "as the mock Puffin, then (i) a saturating nonlinearity, which makes "
        "co-placed strong motifs sub-additive, and (ii) an explicit table of "
        "spacing-dependent pair terms, some with ~10.5 bp helical periodicity. "
        "Its Mode 2 residual is non-zero by construction."
    )

    def _tracks(self, sequence: str) -> dict[str, np.ndarray]:
        L = len(sequence)
        tracks, by_motif = additive_tracks(sequence)
        return {k: saturate(v + pair_cooperativity(by_motif, L, k))
                for k, v in tracks.items()}

    def predict(self, sequence: str, tss_index: int,
                window: tuple[int, int] = (-1000, 1000)) -> ProfilePrediction:
        tracks = self._tracks(sequence)
        positions = None
        out = {}
        for k, v in tracks.items():
            positions, out[k] = _window_slice(v, tss_index, window)
        return ProfilePrediction(
            model=self.name, positions=positions, tracks=out,
            scale=self.scale, output_space=self.output_space, is_mock=True,
            meta={"additive_by_construction": False,
                  "nonlinearity": "tanh saturation + pairwise spacing terms"},
        )


class MockPuffinAdditiveOnly(MockPuffinD):
    """MockPuffinD with the pair table switched off, keeping only saturation.

    Useful for separating "interaction because the response saturates" from
    "interaction because the model has a genuine pair term" - a distinction
    that matters when reading any real model's residuals too.
    """

    name = "mock_puffind_nosat"
    label = "PuffinD (mock, pair terms off)"
    description = (
        "MockPuffinD with the explicit pair table disabled. Any residual left "
        "comes only from the saturating nonlinearity, i.e. from the shape of "
        "the response function rather than from a position-specific "
        "interaction term."
    )

    def _tracks(self, sequence: str) -> dict[str, np.ndarray]:
        tracks, _by_motif = additive_tracks(sequence)
        return {k: saturate(v) for k, v in tracks.items()}
