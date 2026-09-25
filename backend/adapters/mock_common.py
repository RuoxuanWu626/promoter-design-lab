"""Shared machinery for the MOCK prediction adapters.

None of this is a trained model.  It is a transparent, hand-written generator
of plausible-looking initiation profiles, built so that the *structure* of the
experiments in Mode 1 and Mode 2 can be exercised before any real inference is
wired up.

The design is deliberate in one respect that matters for Mode 2:

    ``MockPuffin`` is **additive by construction**.

Its output is a sum of (a) a term linear in local dinucleotide composition and
(b) one independent term per motif occurrence.  So for two elements A and B
placed far enough apart that their footprints never share a PWM scan window,
the interaction residual

    I = f(A+B) - f(A) - f(B) + f(background)

is exactly zero, up to floating point.  ``MockPuffinD`` adds a saturating
nonlinearity and explicit spacing-dependent pair terms, so its residual is not
zero.  That contrast is the point: Mode 2 has a known ground truth to check the
pipeline against, and it is a standing reminder that a non-zero residual is a
property of the *model*, not evidence about biology.

Keeping additivity exact requires care in two places, both of which were got
wrong in an earlier revision and are called out here so they stay right:

1. **No truncation of the hit list.**  Capping "the top N matches per motif"
   makes a hit's contribution depend on what *else* is in the sequence, which
   destroys additivity.  Weak background matches are kept; because they are
   identical across all four constructs they cancel in the residual anyway.
2. **No sequence-dependent renormalisation.**  Nothing may divide by a
   whole-sequence quantity.
"""

from __future__ import annotations

import numpy as np

from motifs import (
    Motif,
    get_library,
    max_possible_score,
    null_threshold,
    scan_sequence,
)

# --- occupancy: how a PWM score becomes a soft "is this motif here" weight --
# The half-occupancy point is set per motif at a fixed false-positive rate
# under a realistic background composition, not at a fixed fraction of the
# maximum score.
# Motifs carry very different amounts of information, and a fixed fraction
# makes a low-information motif fire constantly while a high-information one
# never does -- which is exactly what went wrong before this was calibrated.
OCC_PVALUE = 1e-4
# ... but never demand more than this fraction of the maximum achievable score,
# so a placed consensus instance is always recognised as a match.
OCC_MAX_FRAC = 0.90
OCC_SHARPNESS = 1.35        # logistic slope, per bit of log-odds
OCC_MIN = 0.08              # below this a match is dropped (pure speed knob;
                            # it is a fixed threshold, not a rank cut, so it
                            # does not break additivity)

# --- composition term: linear in local dinucleotide / mononucleotide content -
BASELINE_CONST = 0.04
BASELINE_CPG_WEIGHT = 3.6
BASELINE_GC_WEIGHT = 0.28
BASELINE_SIGMA = 120.0


# ---------------------------------------------------------------------------
# FFT-based gaussian rendering
# ---------------------------------------------------------------------------
# Every motif occurrence contributes a gaussian bump of the same width, so all
# of a motif's occurrences can be rendered with a single convolution of an
# impulse train.  Cost then depends on sequence length, not on how many matches
# there are - which is what makes keeping every weak background match
# affordable.

_KERNEL_FFT_CACHE: dict[tuple[float, int], np.ndarray] = {}


def _next_pow2(n: int) -> int:
    return 1 << max(4, int(n - 1).bit_length())


def _kernel_fft(sigma: float, nfft: int) -> np.ndarray:
    key = (round(float(sigma), 4), nfft)
    cached = _KERNEL_FFT_CACHE.get(key)
    if cached is None:
        j = np.arange(nfft)
        d = np.minimum(j, nfft - j)          # circular distance from 0
        kern = np.exp(-0.5 * (d / sigma) ** 2)
        cached = np.fft.rfft(kern)
        _KERNEL_FFT_CACHE[key] = cached
    return cached


def _conv_gauss(impulse: np.ndarray, sigma: float, length: int) -> np.ndarray:
    """Circular gaussian convolution, padded so wraparound is negligible."""
    nfft = _next_pow2(length + int(14 * sigma) + 2)
    buf = np.zeros(nfft, dtype=np.float64)
    buf[:length] = impulse
    out = np.fft.irfft(np.fft.rfft(buf) * _kernel_fft(sigma, nfft), nfft)
    return out[:length]


def render_bumps(length: int, centers: np.ndarray, amps: np.ndarray,
                 sigma: float, broad_fraction: float = 0.0) -> np.ndarray:
    """Sum of unit-height gaussians of width ``sigma`` at ``centers``.

    ``broad_fraction`` of each bump's amplitude is rendered instead as a
    four-times-wider shoulder, which is how the library distinguishes motifs
    that drive focused initiation from ones that drive dispersed initiation.
    """
    out = np.zeros(length, dtype=np.float64)
    if len(centers) == 0:
        return out

    idx = np.rint(np.asarray(centers, dtype=np.float64)).astype(np.int64)
    amp = np.asarray(amps, dtype=np.float64)
    keep = (idx >= 0) & (idx < length) & (amp != 0.0)
    if not keep.any():
        return out

    impulse = np.zeros(length, dtype=np.float64)
    np.add.at(impulse, idx[keep], amp[keep])

    sigma = max(float(sigma), 0.5)
    if broad_fraction < 1.0:
        out += (1.0 - broad_fraction) * _conv_gauss(impulse, sigma, length)
    if broad_fraction > 0.0:
        out += broad_fraction * _conv_gauss(impulse, sigma * 4.0, length)
    return out


def smooth(signal: np.ndarray, sigma: float) -> np.ndarray:
    """Area-normalised gaussian smoothing (a local weighted average)."""
    if sigma <= 0:
        return signal
    raw = _conv_gauss(signal, sigma, signal.size)
    norm = _conv_gauss(np.ones(signal.size), sigma, signal.size)
    return raw / np.maximum(norm, 1e-12)


def add_bump(target: np.ndarray, center: float, amp: float, sigma: float,
             broad_fraction: float = 0.0) -> None:
    """Add one unit-height gaussian in place (used for sparse extra terms)."""
    L = target.size
    if amp == 0.0:
        return
    sigma = max(float(sigma), 0.5)
    lo = int(max(0, center - 6 * sigma))
    hi = int(min(L, center + 6 * sigma + 1))
    if hi > lo:
        x = np.arange(lo, hi, dtype=np.float64)
        target[lo:hi] += amp * (1.0 - broad_fraction) * np.exp(
            -0.5 * ((x - center) / sigma) ** 2)
    if broad_fraction > 0:
        wide = sigma * 4.0
        lo2 = int(max(0, center - 4 * wide))
        hi2 = int(min(L, center + 4 * wide + 1))
        if hi2 > lo2:
            x2 = np.arange(lo2, hi2, dtype=np.float64)
            target[lo2:hi2] += amp * broad_fraction * np.exp(
                -0.5 * ((x2 - center) / wide) ** 2)


# ---------------------------------------------------------------------------
# Sequence -> features
# ---------------------------------------------------------------------------

def composition_baseline(seq: str) -> np.ndarray:
    """A term *linear* in local dinucleotide / mononucleotide content.

    Linearity is what keeps MockPuffin exactly additive.  Note this is
    invariant under reverse complement (GC content and CpG dinucleotides both
    are), so the same array serves the plus and minus tracks.
    """
    L = len(seq)
    arr = np.frombuffer(seq.encode(), dtype=np.uint8)
    is_c = arr == ord("C")
    is_g = arr == ord("G")
    gc = (is_c | is_g).astype(np.float64)

    cpg = np.zeros(L, dtype=np.float64)
    if L > 1:
        hits = (is_c[:-1] & is_g[1:]).astype(np.float64)
        cpg[:-1] += 0.5 * hits
        cpg[1:] += 0.5 * hits

    return (BASELINE_CONST
            + BASELINE_CPG_WEIGHT * smooth(cpg, BASELINE_SIGMA)
            + BASELINE_GC_WEIGHT * smooth(gc, BASELINE_SIGMA))


def occupancy_threshold(motif: Motif) -> float:
    """Half-occupancy score for one motif (see OCC_PVALUE)."""
    return min(null_threshold(motif, OCC_PVALUE),
               OCC_MAX_FRAC * max_possible_score(motif))


def motif_occupancy(seq: str, motif: Motif) -> np.ndarray:
    """(2, L) soft occupancy in [0, 1] per start index and strand."""
    scores = scan_sequence(seq, motif)
    z = (scores - occupancy_threshold(motif)) * OCC_SHARPNESS
    with np.errstate(over="ignore"):
        occ = 1.0 / (1.0 + np.exp(-np.clip(z, -60, 60)))
    occ[~np.isfinite(scores)] = 0.0
    occ[occ < OCC_MIN] = 0.0
    return occ


def motif_hits_by_motif(seq: str, motif_ids: list[str] | None = None
                        ) -> dict[str, dict[str, np.ndarray]]:
    """Per-motif arrays of match start index, strand sign and occupancy.

    Every match above ``OCC_MIN`` is kept - see the module docstring on why
    there is deliberately no top-N cut.
    """
    lib = get_library()
    ids = motif_ids or list(lib.keys())
    out: dict[str, dict[str, np.ndarray]] = {}
    for mid in ids:
        occ = motif_occupancy(seq, lib[mid])
        plus = np.nonzero(occ[0] > 0)[0]
        minus = np.nonzero(occ[1] > 0)[0]
        out[mid] = {
            "index": np.concatenate([plus, minus]),
            "is_plus": np.concatenate([np.ones(plus.size, bool),
                                       np.zeros(minus.size, bool)]),
            "occupancy": np.concatenate([occ[0][plus], occ[1][minus]]),
        }
    return out


def motif_hits(seq: str, motif_ids: list[str] | None = None) -> list[dict]:
    """Flat hit list, for the code paths that need per-hit bookkeeping."""
    lib = get_library()
    by_motif = motif_hits_by_motif(seq, motif_ids)
    hits: list[dict] = []
    for mid, d in by_motif.items():
        w = lib[mid].width
        for i, is_plus, occ in zip(d["index"], d["is_plus"], d["occupancy"]):
            hits.append({"motif_id": mid, "index": int(i),
                         "strand": "+" if is_plus else "-",
                         "occupancy": float(occ), "width": w})
    return hits


# ---------------------------------------------------------------------------
# Geometry: where an occurrence puts its peak, on each track
# ---------------------------------------------------------------------------

def hit_arrays_geometry(d: dict[str, np.ndarray], motif: Motif, track: str
                        ) -> tuple[np.ndarray, np.ndarray]:
    """Vectorised centres and amplitudes for one motif's matches.

    A motif on the plus strand drives the plus track at ``index + peak_shift``
    at full amplitude, and the minus track at the same place attenuated by
    ``mock_strand_asymmetry``; a minus-strand motif is the mirror image.
    Computing both tracks this way is equivalent to running the model on the
    reverse complement and flipping, but needs only one PWM scan.
    """
    idx = d["index"].astype(np.float64)
    is_plus = d["is_plus"]
    shift = motif.mock_peak_shift
    asym = motif.mock_strand_asymmetry

    centers = np.where(is_plus, idx + shift, idx + motif.width - 1 - shift)
    if track == "plus":
        amps = np.where(is_plus, motif.mock_amplitude, motif.mock_amplitude * asym)
    else:
        amps = np.where(is_plus, motif.mock_amplitude * asym, motif.mock_amplitude)
    return centers, amps * d["occupancy"]


def hit_geometry(h: dict, motif: Motif, track: str) -> tuple[float, float]:
    """Scalar version of :func:`hit_arrays_geometry` for a single hit dict."""
    if h["strand"] == "+":
        center = h["index"] + motif.mock_peak_shift
        amp = motif.mock_amplitude * (
            1.0 if track == "plus" else motif.mock_strand_asymmetry)
    else:
        center = (h["index"] + motif.width - 1) - motif.mock_peak_shift
        amp = motif.mock_amplitude * (
            motif.mock_strand_asymmetry if track == "plus" else 1.0)
    return float(center), float(amp * h["occupancy"])


def additive_tracks(seq: str) -> tuple[dict[str, np.ndarray], dict]:
    """Additive motif + composition drive, both strands, from one scan."""
    L = len(seq)
    lib = get_library()
    base = composition_baseline(seq)
    by_motif = motif_hits_by_motif(seq)

    tracks = {}
    for track in ("plus", "minus"):
        arr = base.copy()
        for mid, d in by_motif.items():
            if d["index"].size == 0:
                continue
            motif = lib[mid]
            centers, amps = hit_arrays_geometry(d, motif, track)
            arr += render_bumps(L, centers, amps, motif.mock_width,
                                motif.mock_broad_fraction)
        tracks[track] = arr
    return tracks, by_motif


def additive_drive(seq: str, return_hits: bool = False):
    """Plus-strand-only convenience wrapper around :func:`additive_tracks`."""
    tracks, by_motif = additive_tracks(seq)
    if return_hits:
        return tracks["plus"], by_motif
    return tracks["plus"]


# ---------------------------------------------------------------------------
# Pair cooperativity table, used ONLY by MockPuffinD
# ---------------------------------------------------------------------------
#   c        signed strength of the extra term (positive = synergy)
#   s0       centre-to-centre spacing, in bp, at which the term is strongest
#   sigma    tolerance for departures from s0
#   helical  depth of a ~10.5 bp periodic modulation (0 = none)
MOCK_PAIR_TERMS: dict[tuple[str, str], dict] = {
    ("inr", "tata"):   {"c": 0.55, "s0": 28.0, "sigma": 9.0, "helical": 0.25},
    ("sp1", "sp1"):    {"c": 0.34, "s0": 22.0, "sigma": 16.0, "helical": 0.35},
    ("ets", "ets"):    {"c": 0.26, "s0": 18.0, "sigma": 12.0, "helical": 0.30},
    ("nfy", "sp1"):    {"c": 0.30, "s0": 40.0, "sigma": 22.0, "helical": 0.15},
    ("nrf1", "sp1"):   {"c": 0.22, "s0": 34.0, "sigma": 20.0, "helical": 0.10},
    ("ets", "sp1"):    {"c": 0.24, "s0": 26.0, "sigma": 18.0, "helical": 0.20},
    ("creb", "nfy"):   {"c": 0.19, "s0": 46.0, "sigma": 24.0, "helical": 0.0},
    ("inr", "yy1"):    {"c": 0.28, "s0": 12.0, "sigma": 8.0, "helical": 0.0},
    ("tata", "yy1"):   {"c": -0.24, "s0": 16.0, "sigma": 11.0, "helical": 0.0},
    ("nfy", "znf143"): {"c": 0.17, "s0": 52.0, "sigma": 26.0, "helical": 0.0},
    ("inr", "u1"):     {"c": 0.12, "s0": 110.0, "sigma": 45.0, "helical": 0.0},
    ("creb", "creb"):  {"c": -0.18, "s0": 14.0, "sigma": 10.0, "helical": 0.0},
}


def pair_term(a: str, b: str) -> dict | None:
    return MOCK_PAIR_TERMS.get(tuple(sorted((a, b))))


PAIR_OCC_MIN = 0.30   # only reasonably confident matches take part in pair terms


def pair_cooperativity(by_motif: dict, length: int, track: str = "plus",
                       max_pair_distance: int = 260) -> np.ndarray:
    """Extra, spacing-dependent contribution used by MockPuffinD.

    Spacing is measured centre-to-centre between the two motif footprints,
    which is what makes the ~10.5 bp helical modulation meaningful.

    Only matches above ``PAIR_OCC_MIN`` participate. That keeps this O(n^2)
    step bounded, and it is a fixed threshold rather than a rank cut, so it
    does not introduce the kind of sequence-dependent coupling that would make
    the *additive* model non-additive.
    """
    lib = get_library()
    out = np.zeros(length, dtype=np.float64)

    strong: list[dict] = []
    for mid, d in by_motif.items():
        if not pair_partners(mid):
            continue
        sel = np.nonzero(d["occupancy"] >= PAIR_OCC_MIN)[0]
        motif = lib[mid]
        for k in sel:
            h = {"motif_id": mid, "index": int(d["index"][k]),
                 "strand": "+" if d["is_plus"][k] else "-",
                 "occupancy": float(d["occupancy"][k])}
            h["footprint_center"] = h["index"] + motif.width / 2.0
            h["peak_center"] = hit_geometry(h, motif, track)[0]
            strong.append(h)

    for i in range(len(strong)):
        for j in range(i + 1, len(strong)):
            hi_, hj = strong[i], strong[j]
            term = pair_term(hi_["motif_id"], hj["motif_id"])
            if term is None:
                continue
            spacing = abs(hi_["footprint_center"] - hj["footprint_center"])
            if spacing > max_pair_distance:
                continue

            envelope = np.exp(-0.5 * ((spacing - term["s0"]) / term["sigma"]) ** 2)
            if term["helical"]:
                envelope *= (1.0 + term["helical"] * np.cos(2 * np.pi * spacing / 10.5))
            strength = term["c"] * envelope * hi_["occupancy"] * hj["occupancy"]
            if abs(strength) < 1e-4:
                continue

            mi, mj = lib[hi_["motif_id"]], lib[hj["motif_id"]]
            center = 0.5 * (hi_["peak_center"] + hj["peak_center"])
            sigma = 0.5 * (mi.mock_width + mj.mock_width) + 4.0
            add_bump(out, center, strength, sigma, 0.0)
    return out


_PARTNERS: dict[str, set[str]] = {}
for _a, _b in MOCK_PAIR_TERMS:
    _PARTNERS.setdefault(_a, set()).add(_b)
    _PARTNERS.setdefault(_b, set()).add(_a)


def pair_partners(motif_id: str) -> set[str]:
    return _PARTNERS.get(motif_id, set())


def saturate(x: np.ndarray, ceiling: float = 2.6, knee: float = 0.75) -> np.ndarray:
    """Soft, monotone saturation. Linear near 0, approaching ``ceiling``.

    This is the only source of *generic* sublinearity in MockPuffinD: two
    strong motifs together give less than the sum of their separate effects,
    with no pair-specific term involved. Mode 2 offers a variant with the pair
    table switched off so the two sources can be told apart.
    """
    below = x <= knee
    out = np.empty_like(x)
    out[below] = x[below]
    span = max(ceiling - knee, 1e-6)
    out[~below] = knee + span * np.tanh((x[~below] - knee) / span)
    return out
