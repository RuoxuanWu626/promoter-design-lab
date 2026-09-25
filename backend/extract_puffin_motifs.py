"""Derive the motif library from the trained Puffin checkpoint.

    python3 extract_puffin_motifs.py [--checkpoint PATH] [--out-dir DIR]

Everything the app shows about a motif comes out of the model itself rather
than out of a textbook:

* **The PWM** is ``conv.weight[i]`` — the learned 51 bp filter — centred per
  position and softmaxed over the four bases. The filter is trimmed to the
  window that actually carries signal, so a 51 bp filter whose informative core
  is 9 bp becomes a 9 bp motif.

* **Orientation.** The reference implementation labels each of the ten filters
  with the strand it fires on (``PUFFIN_MOTIF_STRAND``). Filters labelled ``-``
  are reverse-complemented here so every motif is stored in its canonical
  forward orientation — that is what turns filter 9 from ``CTTTTATAG`` into the
  recognisable TATA box ``CTATAAAAG``.

* **The preferred position** is read off ``deconv.weight``, which is the
  model's own answer to "where does this motif put initiation". Channel
  ``i`` of the motif branch influences output track ``t`` through a 601 bp
  kernel ``deconv.weight[t, i]``; its argmax is the offset from the motif to
  the initiation it drives. No hand-set numbers are involved.

Writes ``data/puffin_motifs.npz`` (PWMs, picked up automatically by
``motifs.get_library``) and ``data/puffin_motifs.json`` (names, consensus,
offsets, effect profiles) alongside it.
"""

from __future__ import annotations

import argparse
import json
import os

import numpy as np
import torch

from models.puffin_arch import (
    PUFFIN_MOTIF_NAMES,
    PUFFIN_MOTIF_STRAND,
    PUFFIN_TRACKS_PLUS,
)

DEFAULT_CKPT = "/gpfs/data/zhou-lab/rxwu/puffin_species/human/puffin.pth"

# Map Puffin's motif names onto the ids this app already uses everywhere.
NAME_TO_ID = {
    "SP": "sp1", "ETS": "ets", "CREB": "creb", "NFY": "nfy", "YY1": "yy1",
    "U1 snRNP": "u1", "Long Inr": "inr", "NRF1": "nrf1", "ZNF143": "znf143",
    "TATA": "tata",
}

# Colours from the reference implementation's own colordict, so plots here
# match plots made with the lab's Puffin code.
PUFFIN_COLORS = {
    "yy1": "#1F77B4", "tata": "#E41A1C", "u1": "#9F9F9F", "nfy": "#00CC96",
    "ets": "#19d3f3", "sp1": "#FF7F0E", "nrf1": "#AB63FA", "znf143": "#17BECF",
    "creb": "#FF6692", "inr": "#95a2be",
}

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def core_window(filt: np.ndarray, frac: float = 0.25) -> tuple[int, int]:
    """Trim a 51 bp filter to the window that carries signal.

    ``filt`` is (4, W), already centred per position. A position counts as
    informative when its base-to-base range exceeds ``frac`` of the largest
    range in the filter.
    """
    spread = filt.max(axis=0) - filt.min(axis=0)
    keep = np.flatnonzero(spread > frac * spread.max())
    return int(keep.min()), int(keep.max())


def extract(checkpoint: str) -> dict:
    sd = torch.load(checkpoint, map_location="cpu")
    conv = sd["conv.weight"].numpy()          # (10, 4, 51)
    deconv = sd["deconv.weight"].numpy()      # (10 tracks, 20 channels, 601)
    cage_plus = PUFFIN_TRACKS_PLUS["FANTOM_CAGE"]
    kernel_centre = deconv.shape[2] // 2

    pwms: dict[str, np.ndarray] = {}
    meta: dict[str, dict] = {}

    for i, (name, strand) in enumerate(zip(PUFFIN_MOTIF_NAMES, PUFFIN_MOTIF_STRAND)):
        mid = NAME_TO_ID[name]
        filt = conv[i] - conv[i].mean(axis=0, keepdims=True)
        lo, hi = core_window(filt)
        core = filt[:, lo:hi + 1]

        # Probabilities over the four bases at each position. No temperature is
        # imposed: the app sets each motif's detection threshold from its own
        # null distribution, so however much information a filter carries is
        # handled automatically.
        exp = np.exp(core - core.max(axis=0, keepdims=True))
        pwm = exp / exp.sum(axis=0, keepdims=True)

        consensus = "".join("ACGT"[j] for j in pwm.argmax(axis=0))
        if strand == "-":
            # Store every motif in canonical forward orientation.
            pwm = pwm[::-1, ::-1].copy()
            consensus = revcomp(consensus)

        # The model's own effect profile for this motif on FANTOM CAGE (+).
        # Channel i fires in the orientation `strand`; channel i+10 is the flip.
        eff_fwd = deconv[cage_plus, i]
        eff_rev = deconv[cage_plus, i + 10]
        effect = eff_rev if strand == "-" else eff_fwd
        offsets = np.arange(deconv.shape[2]) - kernel_centre
        kernel_argmax = int(offsets[int(np.argmax(effect))])

        # deconv is a cross-correlation: out[x] = sum_k W[k] * act[x + k - pad].
        # So an activation at position p drives output at x = p - kernel_argmax,
        # and for that output to land on the TSS the activation must sit at
        # p = kernel_argmax.
        activation_at_tss = kernel_argmax

        # The activation is indexed at the centre of the 51 bp filter
        # (padding=25), and the stored motif is the filter's core window
        # [lo, hi]. So the motif's 5' end on the plus strand sits here:
        motif_5p_for_tss_peak = activation_at_tss - 25 + lo
        # ...and the peak it drives is this far from that 5' end:
        peak_shift = -motif_5p_for_tss_peak

        info = float(np.mean(2.0 + (pwm * np.log2(pwm + 1e-9)).sum(axis=0)))

        pwms[mid] = pwm.astype(np.float64)
        meta[mid] = {
            "id": mid,
            "puffin_name": name,
            "filter_index": i,
            "filter_strand": strand,
            "core": [lo, hi],
            "width": int(pwm.shape[1]),
            "consensus": consensus,
            "color": PUFFIN_COLORS[mid],
            # Offset from the motif's 5' end to the initiation peak it drives,
            # straight out of deconv.weight.
            "peak_shift": peak_shift,
            # Where the motif's 5' end sits if that peak is to land on the TSS.
            "typical_offset": motif_5p_for_tss_peak,
            "kernel_argmax": kernel_argmax,
            "effect_max": float(effect.max()),
            "effect_min": float(effect.min()),
            "effect_profile": np.round(effect[::4], 5).tolist(),
            "effect_profile_stride": 4,
            "effect_offsets": offsets[::4].tolist(),
            "mean_information_bits": round(info, 3),
        }

    return {"pwms": pwms, "meta": meta, "checkpoint": checkpoint}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--checkpoint", default=os.environ.get("PUFFIN_WEIGHTS", DEFAULT_CKPT))
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--out-dir", default=os.path.join(here, "data"))
    args = ap.parse_args()

    if not os.path.exists(args.checkpoint):
        raise SystemExit(f"checkpoint not found: {args.checkpoint}")

    res = extract(args.checkpoint)
    os.makedirs(args.out_dir, exist_ok=True)

    npz = os.path.join(args.out_dir, "puffin_motifs.npz")
    np.savez(npz, **res["pwms"])
    js = os.path.join(args.out_dir, "puffin_motifs.json")
    with open(js, "w") as fh:
        json.dump({"checkpoint": res["checkpoint"], "motifs": res["meta"]}, fh, indent=1)

    print(f"checkpoint {res['checkpoint']}")
    print(f"wrote {npz}\nwrote {js}\n")
    print(f"  {'id':8s} {'name':10s} {'w':>3s} {'bits':>5s} {'5p@TSS':>7s} {'peak':>6s}  consensus")
    for mid, m in res["meta"].items():
        print(f"  {mid:8s} {m['puffin_name']:10s} {m['width']:3d} "
              f"{m['mean_information_bits']:5.2f} {m['typical_offset']:+7d} "
              f"{m['peak_shift']:+6d}  {m['consensus']}")


if __name__ == "__main__":
    main()
