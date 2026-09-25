"""Metrics, tau specificity, and the configurable weighted score.

Every metric declares a display name, a direction ("higher is better" or not),
and a normalisation range.  The weighted score is

    score = sum_i  w_i * normalise_i(metric_i)

with ``normalise_i`` mapping the metric onto [0, 1] using its declared range,
so that weights are comparable across metrics with different units.  The API
returns the per-metric contribution as well as the total, so a leaderboard
entry can always be taken apart.
"""

from __future__ import annotations

import numpy as np

# ---------------------------------------------------------------------------
# Tau
# ---------------------------------------------------------------------------

def tau_specificity(values: np.ndarray) -> float:
    """Yanai's tau over per-cell-type activities.

    tau = sum_i (1 - x_i / x_max) / (n - 1)

    0 = equal in every cell type, 1 = confined to one cell type.

    Caveats worth stating in any write-up:
    * tau is scale-dependent. Computing it on log-transformed activities gives
      a different number than on linear activities; this project computes it on
      whatever ``activity_space`` is requested and records which was used.
    * tau is undefined for all-zero input and is reported as 0 there.
    * tau depends on which cell types are in the panel. Adding or removing
      cell types changes it, so panels must match when tau values are compared.
    """
    v = np.asarray(values, dtype=np.float64)
    v = np.where(np.isfinite(v), v, 0.0)
    v = np.clip(v, 0.0, None)
    n = v.size
    if n < 2:
        return 0.0
    vmax = v.max()
    if vmax <= 0:
        return 0.0
    return float(np.sum(1.0 - v / vmax) / (n - 1))


def gini(values: np.ndarray) -> float:
    """Gini coefficient, offered alongside tau as a sanity check."""
    v = np.sort(np.clip(np.asarray(values, dtype=np.float64), 0, None))
    n = v.size
    if n == 0 or v.sum() <= 0:
        return 0.0
    idx = np.arange(1, n + 1)
    return float((2 * np.sum(idx * v) / (n * np.sum(v))) - (n + 1) / n)


def activity_from_profile(profile: np.ndarray, positions: np.ndarray,
                          window: tuple[int, int] = (-200, 200),
                          method: str = "mean") -> float:
    """Collapse a profile to one scalar activity within a window around the TSS."""
    prof = np.asarray(profile, dtype=np.float64)
    mask = (positions >= window[0]) & (positions <= window[1]) & np.isfinite(prof)
    seg = prof[mask]
    if seg.size == 0:
        return 0.0
    if method == "sum":
        return float(seg.sum())
    if method == "max":
        return float(seg.max())
    if method == "mean":
        return float(seg.mean())
    return float(seg.mean())


def to_linear(values: np.ndarray, output_space: str) -> np.ndarray:
    """Undo a log-like output space, for metrics that want linear units."""
    v = np.asarray(values, dtype=np.float64)
    if output_space == "log":
        return np.power(10.0, v) - 1.0
    return v


# ---------------------------------------------------------------------------
# Profile metrics (Puffin-style models)
# ---------------------------------------------------------------------------

def profile_metrics(positions: np.ndarray, plus: np.ndarray,
                    minus: np.ndarray | None = None,
                    focus_window: int = 50,
                    score_window: tuple[int, int] = (-500, 500)) -> dict:
    positions = np.asarray(positions)
    plus = np.asarray(plus, dtype=np.float64)
    minus = np.zeros_like(plus) if minus is None else np.asarray(minus, dtype=np.float64)

    # Positions a model could not validly predict arrive as NaN (Puffin masks
    # anywhere its 325 bp of context ran off the construct). They have to be
    # dropped here, or argmax lands on a NaN and every metric downstream
    # becomes NaN.
    mask = (positions >= score_window[0]) & (positions <= score_window[1])
    mask &= np.isfinite(plus)
    if minus is not None and minus.size == plus.size:
        mask &= np.isfinite(minus)
    p = plus[mask]
    m = minus[mask]
    pos = positions[mask]
    if p.size == 0:
        return {}

    peak_idx = int(np.argmax(p))
    peak_height = float(p[peak_idx])
    peak_position = int(pos[peak_idx])

    total = float(p.sum())
    focus_mask = np.abs(pos) <= focus_window
    focused = float(p[focus_mask].sum())
    focus_fraction = focused / total if total > 0 else 0.0

    # Width at half maximum around the peak, in bp.
    half = peak_height / 2.0
    left = peak_idx
    while left > 0 and p[left] > half:
        left -= 1
    right = peak_idx
    while right < p.size - 1 and p[right] > half:
        right += 1
    fwhm = float(pos[right] - pos[left])

    # Directionality asks whether initiation *at the designed TSS* goes one
    # way. Two things have to be handled or the number is meaningless:
    #   - both strands share the same composition baseline, so compare signal
    #     above each track's own floor, not raw values;
    #   - chance motif matches elsewhere in the window are strand-balanced and
    #     swamp a whole-window sum, so look only near the TSS.
    floor_p = float(np.percentile(p, 10))
    floor_m = float(np.percentile(m, 10)) if m.size else 0.0
    near = np.abs(pos) <= focus_window
    peak_p = float(np.clip(p[near] - floor_p, 0, None).max()) if near.any() else 0.0
    peak_m = float(np.clip(m[near] - floor_m, 0, None).max()) if near.any() else 0.0
    directionality = ((peak_p - peak_m) / (peak_p + peak_m)) if (peak_p + peak_m) > 0 else 0.0
    tot_p, tot_m = float(p.sum()), float(m.sum())
    exc_p, exc_m = peak_p, peak_m

    return {
        "peak_height": peak_height,
        "peak_position": peak_position,
        "tss_offset_abs": abs(peak_position),
        "total_output": total,
        "focus_fraction": float(focus_fraction),
        "fwhm": fwhm,
        "directionality": float(directionality),
        "plus_peak_excess": exc_p,
        "minus_peak_excess": exc_m,
        "minus_total": tot_m,
    }


# ---------------------------------------------------------------------------
# Cell-type metrics (Mode 1)
# ---------------------------------------------------------------------------

def celltype_metrics(cell_types: list[str], activities: np.ndarray,
                     target: str | None = None) -> dict:
    act = np.asarray(activities, dtype=np.float64)
    order = np.argsort(-act)
    strongest = cell_types[int(order[0])] if len(cell_types) else None

    out = {
        "tau": tau_specificity(act),
        "gini": gini(act),
        "strongest_cell_type": strongest,
        "strongest_activity": float(act[order[0]]) if act.size else 0.0,
        "activities": {c: float(a) for c, a in zip(cell_types, act)},
        "ranking": [{"cell_type": cell_types[i], "activity": float(act[i])}
                    for i in order],
    }

    if target and target in cell_types:
        ti = cell_types.index(target)
        off = np.delete(act, ti)
        out.update({
            "target": target,
            "target_activity": float(act[ti]),
            "target_rank": int(np.where(order == ti)[0][0]) + 1,
            "offtarget_mean": float(off.mean()) if off.size else 0.0,
            "offtarget_max": float(off.max()) if off.size else 0.0,
            "offtarget_max_cell_type": (cell_types[int(np.delete(np.arange(act.size), ti)[int(np.argmax(off))])]
                                        if off.size else None),
            "target_minus_offtarget_max": float(act[ti] - (off.max() if off.size else 0.0)),
            "log2_target_over_offtarget_mean": float(
                np.log2((act[ti] + 1e-6) / (off.mean() + 1e-6))) if off.size else 0.0,
            "is_target_strongest": bool(strongest == target),
        })
    return out


def decompose_specificity_change(base_act: np.ndarray, new_act: np.ndarray,
                                 cell_types: list[str], target: str) -> dict:
    """Split a change in specificity into target gain vs off-target loss.

    Mode 1 asks specifically whether specificity was won by turning the target
    *up* or by turning everything else *down*.  Those are different designs and
    the summary numbers (tau, ratios) cannot tell them apart on their own.
    """
    base = np.asarray(base_act, dtype=np.float64)
    new = np.asarray(new_act, dtype=np.float64)
    ti = cell_types.index(target)

    d_target = float(new[ti] - base[ti])
    off_base = np.delete(base, ti)
    off_new = np.delete(new, ti)
    d_off_mean = float(off_new.mean() - off_base.mean()) if off_base.size else 0.0

    # The contrast is what specificity actually is: target minus off-target.
    # It moves by d_target - d_off_mean, and those are two separate levers.
    contrast_gain = d_target - d_off_mean
    total_move = abs(d_target) + abs(d_off_mean)
    eps = 1e-9

    def direction(x: float, scale: float) -> str:
        if abs(x) <= max(scale * 0.02, eps):
            return "flat"
        return "up" if x > 0 else "down"

    scale = max(abs(base[ti]), float(off_base.mean()) if off_base.size else 0.0, eps)
    t_dir = direction(d_target, scale)
    o_dir = direction(d_off_mean, scale)

    if t_dir == "flat" and o_dir == "flat":
        mechanism = "no change"
    elif t_dir == "up" and o_dir in ("down", "flat"):
        mechanism = "target gain"
    elif t_dir in ("up", "flat") and o_dir == "down":
        mechanism = "off-target loss"
    elif t_dir == "up" and o_dir == "up":
        mechanism = ("target gain outpaces off-target" if contrast_gain > 0
                     else "off-target rises faster than target")
    elif t_dir == "down" and o_dir == "down":
        mechanism = ("off-target falls faster than target" if contrast_gain > 0
                     else "target falls faster than off-target")
    elif t_dir == "down":
        mechanism = "target loss"
    else:
        mechanism = "off-target gain"

    return {
        "delta_target": d_target,
        "delta_offtarget_mean": d_off_mean,
        "delta_tau": tau_specificity(new) - tau_specificity(base),
        "contrast_gain": contrast_gain,
        "target_direction": t_dir,
        "offtarget_direction": o_dir,
        # How much of the *movement* came from the target arm rather than the
        # off-target arm. Near 1: specificity was won by turning the target up
        # (or down). Near 0: it was won by moving everything else.
        "fraction_from_target_arm": (abs(d_target) / total_move) if total_move > eps else 0.0,
        "mechanism": mechanism,
    }


# ---------------------------------------------------------------------------
# Weighted score
# ---------------------------------------------------------------------------

# name -> (label, lo, hi, higher_is_better, group)
METRIC_SPECS: dict[str, tuple] = {
    "peak_height":        ("Peak height", 0.0, 4.0, True, "profile"),
    "total_output":       ("Total output", 0.0, 900.0, True, "profile"),
    "focus_fraction":     ("Focus (±50 bp)", 0.0, 1.0, True, "profile"),
    "tss_offset_abs":     ("TSS positioning error", 0.0, 200.0, False, "profile"),
    "fwhm":               ("Peak width (FWHM)", 0.0, 300.0, False, "profile"),
    "directionality":     ("Directionality", -1.0, 1.0, True, "profile"),
    "tau":                ("Tau specificity", 0.0, 1.0, True, "celltype"),
    "target_activity":    ("Target activity", 0.0, 3.0, True, "celltype"),
    "offtarget_max":      ("Max off-target activity", 0.0, 3.0, False, "celltype"),
    "offtarget_mean":     ("Mean off-target activity", 0.0, 3.0, False, "celltype"),
    "target_minus_offtarget_max": ("Target − best off-target", -2.0, 2.0, True, "celltype"),
    "log2_target_over_offtarget_mean": ("log2(target / mean off-target)", -3.0, 3.0, True, "celltype"),
}

DEFAULT_WEIGHTS_MODE1 = {
    "target_activity": 1.0,
    "offtarget_max": 1.0,
    "tau": 1.5,
    "focus_fraction": 0.25,
    "tss_offset_abs": 0.25,
}

DEFAULT_WEIGHTS_DESIGN = {
    "peak_height": 1.0,
    "focus_fraction": 0.75,
    "tss_offset_abs": 0.75,
    "directionality": 0.5,
}


def normalise(name: str, value: float) -> float:
    spec = METRIC_SPECS.get(name)
    if spec is None or value is None:
        return 0.0
    _, lo, hi, higher, _ = spec
    if hi == lo:
        return 0.0
    x = (float(value) - lo) / (hi - lo)
    x = float(np.clip(x, 0.0, 1.0))
    return x if higher else 1.0 - x


def weighted_score(metrics: dict, weights: dict[str, float] | None = None) -> dict:
    weights = weights or DEFAULT_WEIGHTS_DESIGN
    contributions = {}
    total = 0.0
    wsum = 0.0
    for name, w in weights.items():
        if name not in metrics or metrics.get(name) is None:
            continue
        w = float(w)
        n = normalise(name, metrics[name])
        c = w * n
        contributions[name] = {
            "raw": float(metrics[name]),
            "normalised": round(n, 4),
            "weight": w,
            "contribution": round(c, 4),
            "label": METRIC_SPECS.get(name, (name,))[0],
        }
        total += c
        wsum += abs(w)

    return {
        "score": round(total, 4),
        "score_normalised": round(total / wsum, 4) if wsum > 0 else 0.0,
        "weight_sum": wsum,
        "contributions": contributions,
    }


def metric_catalogue() -> list[dict]:
    return [{
        "name": k, "label": v[0], "lo": v[1], "hi": v[2],
        "higher_is_better": v[3], "group": v[4],
    } for k, v in METRIC_SPECS.items()]


# ---------------------------------------------------------------------------
# Model agreement
# ---------------------------------------------------------------------------

def pearson_r(a, b) -> float:
    """Pearson correlation, NaN-safe, returning 0.0 when either side is flat."""
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    n = min(a.size, b.size)
    a, b = a[:n], b[:n]
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return 0.0
    a, b = a[ok] - a[ok].mean(), b[ok] - b[ok].mean()
    den = float(np.sqrt((a * a).sum() * (b * b).sum()))
    return float((a * b).sum() / den) if den > 0 else 0.0


def _rank(x: np.ndarray) -> np.ndarray:
    """Average ranks, so ties do not bias the rank correlation."""
    order = np.argsort(x, kind="stable")
    ranks = np.empty(x.size, dtype=np.float64)
    ranks[order] = np.arange(x.size, dtype=np.float64)
    sx = x[order]
    i = 0
    while i < sx.size:
        j = i
        while j + 1 < sx.size and sx[j + 1] == sx[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = (i + j) / 2.0
        i = j + 1
    return ranks


def spearman_r(a, b) -> float:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    n = min(a.size, b.size)
    a, b = a[:n], b[:n]
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return 0.0
    return pearson_r(_rank(a[ok]), _rank(b[ok]))


def normalise_profile(x, output_space: str, transform: str = "log1p") -> np.ndarray:
    """Put a profile on a comparable footing before correlating.

    Two models can agree perfectly on shape and still be on wildly different
    numeric scales (Puffin emits a small non-negative signal, the CAGE mock
    emits a log10 CPM). Pearson r is invariant to an affine change, so scale
    alone does not matter — but a log-like output compared against a linear one
    is NOT an affine relationship, and correlating them raw measures the wrong
    thing. So both sides are brought to the same space first.
    """
    v = np.asarray(x, dtype=np.float64)
    if transform == "none":
        return v
    if output_space == "log":
        v = np.power(10.0, v) - 1.0          # back to a count-like scale
    return np.log1p(np.clip(v, 0.0, None))   # ...then a common log space


def profile_agreement(entries: list[dict], transform: str = "log1p",
                      window: tuple[int, int] | None = None) -> dict:
    """Pairwise agreement between model profiles.

    ``entries`` is a list of ``{"name", "label", "positions", "values",
    "output_space", "role"}``. ``role`` is free text such as "simple" or
    "deep"; it is only used to label the comparison.

    Reported per pair: Pearson r (linear agreement on the common scale),
    Spearman rho (rank agreement, insensitive to the transform), and the
    Pearson r of the two profiles after peak-normalising, which answers
    "do they agree on the shape" separately from "do they agree on the level".
    """
    prepared = []
    for e in entries:
        pos = np.asarray(e["positions"])
        val = normalise_profile(e["values"], e.get("output_space", "linear"), transform)
        if window is not None:
            m = (pos >= window[0]) & (pos <= window[1])
            pos, val = pos[m], val[m]
        prepared.append({**e, "pos": pos, "val": val})

    pairs = []
    for i in range(len(prepared)):
        for j in range(i + 1, len(prepared)):
            a, b = prepared[i], prepared[j]
            n = min(a["val"].size, b["val"].size)
            va, vb = a["val"][:n], b["val"][:n]
            pa = va / va.max() if va.size and va.max() > 0 else va
            pb = vb / vb.max() if vb.size and vb.max() > 0 else vb
            pairs.append({
                "a": a["name"], "b": b["name"],
                "a_label": a.get("label", a["name"]),
                "b_label": b.get("label", b["name"]),
                "a_role": a.get("role"), "b_role": b.get("role"),
                "n_positions": int(n),
                "pearson_r": round(pearson_r(va, vb), 5),
                "spearman_r": round(spearman_r(va, vb), 5),
                "shape_pearson_r": round(pearson_r(pa, pb), 5),
                "transform": transform,
            })
    return {
        "transform": transform,
        "window": list(window) if window else None,
        "series": [{"name": p["name"], "label": p.get("label", p["name"]),
                    "role": p.get("role"),
                    "positions": p["pos"].astype(int).tolist(),
                    "values": np.round(p["val"], 6).tolist(),
                    "output_space": p.get("output_space", "linear")}
                   for p in prepared],
        "pairs": pairs,
        "note": (
            "Pearson r is computed after bringing both models to a common "
            "log1p scale, because correlating a log-like output against a "
            "linear one raw would measure the transform rather than the "
            "agreement. r is invariant to scale and offset, so it answers "
            "'do these two models rise and fall together along the sequence', "
            "not 'do they predict the same magnitude'."),
    }
