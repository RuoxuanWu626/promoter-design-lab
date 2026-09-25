"""Real Puffin inference (human).

Loads the published Puffin checkpoint and predicts initiation profiles for a
designed construct. The architecture lives in ``models/puffin_arch.py``; see
that file for what was transcribed and what was checked.

Three things about this adapter are worth knowing before reading its numbers.

**Context.** Puffin's output is only valid where the 601 bp deconvolutions have
full context, so the reference implementation trims ``PUFFIN_TRIM`` = 325 bp
from each end. To return a window of ``[w0, w1]`` around the TSS this adapter
needs sequence covering ``[w0 - 325, w1 + 325]``. If the construct is shorter it
is zero-padded (an all-zero one-hot column, i.e. "no base"), and the number of
padded bases is reported in ``meta`` and surfaced in the UI. Padding is
identical across the four Mode 2 constructs, so it cancels in the residual, but
it does affect absolute values near the edges.

**Scale.** The output is a softplus, so non-negative, trained against CAGE /
RAMPAGE / GRO-cap / PRO-cap coverage. The training loss normalised each profile
before comparing, so the *shape* is what was fit and the absolute height is only
weakly constrained. Treated here as a linear signal, which is what makes the
Mode 2 residual a test of departure from additivity in signal units.

**Determinism.** The model is in ``eval()`` and inference runs under
``torch.inference_mode``; repeated calls on the same sequence are bit-identical,
which the four-construct design relies on.
"""

from __future__ import annotations

import os
import threading

import numpy as np

from adapters.base import ProfileAdapter, ProfilePrediction

DEFAULT_WEIGHTS = "/gpfs/data/zhou-lab/rxwu/puffin_species/human/puffin.pth"

_BASE_INDEX = {"A": 0, "a": 0, "C": 1, "c": 1, "G": 2, "g": 2, "T": 3, "t": 3}


def one_hot(seq: str) -> np.ndarray:
    """(4, L) one-hot. Unknown bases (N and anything else) are all-zero."""
    enc = np.zeros((4, len(seq)), dtype=np.float32)
    for i, ch in enumerate(seq):
        j = _BASE_INDEX.get(ch)
        if j is not None:
            enc[j, i] = 1.0
    return enc


class RealPuffin(ProfileAdapter):
    name = "puffin"
    label = "Puffin (real, human)"
    is_mock = False
    output_space = "linear"
    scale = ("Puffin predicted initiation signal, FANTOM CAGE track "
             "(softplus output, non-negative; profile shape is what the model "
             "was trained on, absolute height only weakly constrained)")
    description = (
        "The published human Puffin model: ten learned motif filters, an "
        "initiator branch and a trinucleotide branch, deconvolved into "
        "strand-specific initiation tracks. This is the simple, interpretable "
        "model — every motif in the library is one of its filters."
    )

    #: which of the five assays to report as plus/minus
    ASSAY = "FANTOM_CAGE"

    def __init__(self, weights: str | None = None, device: str | None = None):
        self._weights = weights or os.environ.get("PUFFIN_WEIGHTS", DEFAULT_WEIGHTS)
        self._device_pref = device or os.environ.get("PUFFIN_DEVICE", "cpu")
        self._net = None
        self._torch = None
        self._lock = threading.Lock()

    # -- availability -------------------------------------------------------
    def available(self) -> tuple[bool, str]:
        try:
            import torch  # noqa: F401
        except ImportError:
            return False, ("PyTorch is not importable by this interpreter. "
                           "Start the server with the 'alphagenome' env.")
        if not os.path.exists(self._weights):
            return False, f"checkpoint not found: {self._weights}"
        return True, ""

    def _load(self):
        if self._net is not None:
            return self._net
        with self._lock:
            if self._net is not None:
                return self._net
            import torch
            from models.puffin_arch import Puffin

            net = Puffin()
            sd = torch.load(self._weights, map_location="cpu")
            missing, unexpected = net.load_state_dict(sd, strict=False)
            if missing or unexpected:
                raise RuntimeError(
                    f"checkpoint does not match the architecture "
                    f"(missing={missing}, unexpected={unexpected})")
            device = self._device_pref
            if device.startswith("cuda") and not torch.cuda.is_available():
                device = "cpu"
            net.eval().to(device)
            self._torch = torch
            self._device = device
            self._net = net
            return net

    # -- prediction ---------------------------------------------------------
    def _context(self, sequence: str, tss_index: int,
                 window: tuple[int, int]) -> tuple[str, int, int]:
        """Sequence covering the window plus Puffin's 325 bp of context each side.

        Returns the padded context and how many bases had to be invented on
        each side. Those padded bases are the reason the caller then masks part
        of the output: a position whose 325 bp of context contains padding is
        not a prediction about the construct, it is an edge artefact, and it
        can be very large.
        """
        from models.puffin_arch import PUFFIN_TRIM
        start = tss_index + int(window[0]) - PUFFIN_TRIM
        end = tss_index + int(window[1]) + PUFFIN_TRIM + 1
        left_pad = max(0, -start)
        right_pad = max(0, end - len(sequence))
        sub = sequence[max(0, start):min(len(sequence), end)]
        return "N" * left_pad + sub + "N" * right_pad, left_pad, right_pad

    def _forward(self, seqs: list[str]) -> np.ndarray:
        net = self._load()
        torch = self._torch
        x = np.stack([one_hot(s) for s in seqs])
        t = torch.from_numpy(x).to(self._device)
        with torch.inference_mode():
            out = net(t)
        return out.float().cpu().numpy()

    def predict(self, sequence: str, tss_index: int,
                window: tuple[int, int] = (-1000, 1000)) -> ProfilePrediction:
        return self.predict_batch([(sequence, tss_index)], window)[0]

    def predict_batch(self, items: list[tuple[str, int]],
                      window: tuple[int, int] = (-1000, 1000)
                      ) -> list[ProfilePrediction]:
        from models.puffin_arch import (
            PUFFIN_TRACKS_MINUS, PUFFIN_TRACKS_PLUS, PUFFIN_TRIM)

        ctx = [self._context(s, t, window) for s, t in items]
        seqs = [c[0] for c in ctx]
        # All contexts are the same length by construction, so this batches.
        out = self._forward(seqs)

        positions = np.arange(int(window[0]), int(window[1]) + 1, dtype=np.int64)
        preds = []
        for k, (_c, left_pad, right_pad) in enumerate(ctx):
            pad = left_pad + right_pad
            core = out[k][:, PUFFIN_TRIM:out.shape[2] - PUFFIN_TRIM]
            n = min(core.shape[1], positions.size)
            tracks = {
                "plus": core[PUFFIN_TRACKS_PLUS[self.ASSAY]][:n],
                "minus": core[PUFFIN_TRACKS_MINUS[self.ASSAY]][:n],
            }
            # The other four assays come free; the UI plots plus/minus but
            # exports and the correlation panel can use any of them.
            for assay in PUFFIN_TRACKS_PLUS:
                if assay == self.ASSAY:
                    continue
                tracks[f"{assay}_plus"] = core[PUFFIN_TRACKS_PLUS[assay]][:n]
                tracks[f"{assay}_minus"] = core[PUFFIN_TRACKS_MINUS[assay]][:n]

            if n < positions.size:
                tracks = {kk: np.pad(v, (0, positions.size - n),
                                     constant_values=np.nan)
                          for kk, v in tracks.items()}

            # Mask the positions whose context ran off the end of the
            # construct. Output index i is valid only when its 325 bp of
            # context on both sides came from real sequence, i.e. when
            # left_pad <= i < len(window) - right_pad. Without this the plot is
            # dominated by an edge artefact that is not a prediction at all.
            if left_pad or right_pad:
                tracks = {kk: np.asarray(v, dtype=np.float64).copy()
                          for kk, v in tracks.items()}
                hi_valid = positions.size - right_pad
                for v in tracks.values():
                    if left_pad:
                        v[:min(left_pad, v.size)] = np.nan
                    if right_pad and hi_valid < v.size:
                        v[max(0, hi_valid):] = np.nan

            valid_lo = int(positions[0] + left_pad) if left_pad else int(positions[0])
            valid_hi = int(positions[-1] - right_pad) if right_pad else int(positions[-1])

            preds.append(ProfilePrediction(
                model=self.name, positions=positions, tracks=tracks,
                scale=self.scale, output_space=self.output_space, is_mock=False,
                meta={"assay": self.ASSAY, "padded_bp": int(pad),
                      "left_pad": int(left_pad), "right_pad": int(right_pad),
                      "valid_range": [valid_lo, valid_hi],
                      "device": getattr(self, "_device", "cpu"),
                      "checkpoint": os.path.basename(self._weights),
                      "trim": PUFFIN_TRIM,
                      "min_construct_bp": int((window[1] - window[0] + 1) + 2 * PUFFIN_TRIM),
                      "padding_note": (
                          f"Puffin needs {PUFFIN_TRIM} bp of context beyond the "
                          f"requested window. {pad} bp were missing, so "
                          f"{valid_lo:+d}..{valid_hi:+d} is shown and the rest is "
                          f"blank. Use a construct of at least "
                          f"{(window[1]-window[0]+1) + 2*PUFFIN_TRIM} bp, or a "
                          f"narrower window, to fill it in.") if pad else None},
            ))
        return preds

    # -- interpretability ---------------------------------------------------
    def motif_activity(self, sequence: str, tss_index: int,
                       window: tuple[int, int] = (-1000, 1000)) -> dict:
        """Per-motif activation along the construct, from the model itself.

        This is what Puffin actually detects, as opposed to a separate PWM
        scan: channels 0..9 are the ten motifs in the orientation the filter
        was labelled with, channels 10..19 the opposite orientation. The design
        view uses it to show which motifs are present — including ones that
        occur in the random background rather than being placed.
        """
        from models.puffin_arch import (
            PUFFIN_MOTIF_NAMES, PUFFIN_MOTIF_STRAND, PUFFIN_TRIM)
        from motifs import library_order

        net = self._load()
        torch = self._torch
        seq, left_pad, right_pad = self._context(sequence, tss_index, window)
        pad = left_pad + right_pad
        x = torch.from_numpy(one_hot(seq)[None]).to(self._device)
        with torch.inference_mode():
            act = net.motif_activations(x)[0].float().cpu().numpy()

        core = act[:, PUFFIN_TRIM:act.shape[1] - PUFFIN_TRIM]
        positions = np.arange(int(window[0]), int(window[1]) + 1, dtype=np.int64)
        n = min(core.shape[1], positions.size)

        # Map Puffin's filter order onto the app's motif ids.
        from extract_puffin_motifs import NAME_TO_ID
        out = {}
        for i, (nm, st) in enumerate(zip(PUFFIN_MOTIF_NAMES, PUFFIN_MOTIF_STRAND)):
            mid = NAME_TO_ID[nm]
            fwd = core[i][:n]
            rev = core[i + 10][:n]
            out[mid] = {
                "name": nm,
                "forward": np.round(fwd, 5).tolist(),
                "reverse": np.round(rev, 5).tolist(),
                "forward_strand": st,
                "max": float(max(fwd.max(initial=0.0), rev.max(initial=0.0))),
            }
        return {"positions": positions[:n].tolist(), "motifs": out,
                "order": [m for m in library_order() if m in out],
                "padded_bp": int(pad)}

    def info(self) -> dict:
        d = super().info()
        d["checkpoint"] = self._weights
        d["assay"] = self.ASSAY
        return d


class RealPuffinProCap(RealPuffin):
    """Same model, reporting the PRO-cap assay instead of FANTOM CAGE."""

    name = "puffin_procap"
    label = "Puffin (real, PRO-cap track)"
    ASSAY = "PRO_CAP"
    scale = ("Puffin predicted initiation signal, PRO-cap track "
             "(softplus output, non-negative)")
    description = ("The same Puffin checkpoint, reporting its PRO-cap output "
                   "track rather than FANTOM CAGE. Useful for checking that a "
                   "design's behaviour is not specific to one assay.")
