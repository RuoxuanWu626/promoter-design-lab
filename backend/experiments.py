"""The two research modes.

Mode 1 - where is cell-type specificity encoded?
    :func:`position_scan` slides one element (a library motif, a CpG-rich
    segment, or a literal sequence segment) across a range of positions
    relative to the TSS, and records the full cell-type x position response.

Mode 2 - what motif interactions do Puffin and PuffinD learn?
    :func:`pair_interaction` builds the four matched constructs
    (background / A only / B only / A+B) and returns each model's profiles plus
    the position-wise residual.  :func:`pair_grid` runs that over all pairs from
    the 10-motif library, and :func:`spacing_curve` runs it over a range of
    spacings.

Reading the Mode 2 residual
---------------------------
    I(x) = f(A+B)(x) - f(A)(x) - f(B)(x) + f(background)(x)

What this is, and is not:

* It is a **model-space** quantity.  It measures whether the model's output
  combines the two elements additively *on the scale the model emits*.  Every
  result carries the adapter's ``output_space``; on a log-like output a zero
  residual means the effects are multiplicative in the underlying units, not
  additive.  ``scale_used`` in the result records which scale was applied.
* Non-zero I does **not** establish biological cooperativity.  It establishes
  that this model is not additive for this pair in this background.
* Two profiles merely looking different is weaker still: that is not an
  interaction measurement at all, which is why the four-construct difference
  is computed rather than comparing f(A+B) to f(A) or f(B).
* The residual is background-dependent.  Everything here supports repeating
  the comparison over several matched backgrounds and reports the spread, so
  a result that only holds in one background is visible as such.
"""

from __future__ import annotations

import numpy as np

import adapters
import scoring
from motifs import CPG_SEGMENT_ID, get_library, library_order
from sequence import (
    Placement,
    build_construct,
    element_width,
    random_background,
    sequence_hash,
    to_placements,
)

DEFAULT_WINDOW = (-500, 500)


# ---------------------------------------------------------------------------
# Matched backgrounds
# ---------------------------------------------------------------------------

def matched_backgrounds(n: int = 3, length: int = 1001, gc: float = 0.45,
                        cpg_oe: float = 0.25, base_seed: int = 1000) -> list[dict]:
    """``n`` independent backgrounds drawn with identical parameters.

    "Matched" means matched in *parameters* (length, GC, CpG o/e), not
    identical sequence. Repeating an experiment across these is what separates
    a real model behaviour from a quirk of one random background.
    """
    return [random_background(length=length, gc=gc, cpg_oe=cpg_oe,
                              seed=base_seed + i) for i in range(int(n))]


# ---------------------------------------------------------------------------
# Mode 1: position scan
# ---------------------------------------------------------------------------

def position_scan(background: dict, element: dict, offsets: list[int],
                  model: str = "mock_alphagenome",
                  cell_types: list[str] | None = None,
                  target: str | None = None,
                  activity_window: tuple[int, int] = (-200, 200),
                  activity_method: str = "mean",
                  window: tuple[int, int] = DEFAULT_WINDOW,
                  progress: dict | None = None) -> dict:
    """Slide one element across ``offsets`` and record the cell-type response."""
    seq = background["sequence"]
    tss = background["tss_index"]
    adapter = adapters.get_celltype_adapter(model)
    cts = cell_types or adapter.default_cell_types

    # Reference: the background with nothing placed.
    ref = adapters.predict_celltype(model, seq, tss, cts, window)
    ref_act = np.array([
        scoring.activity_from_profile(ref.profiles[i], ref.positions,
                                      activity_window, activity_method)
        for i in range(len(ref.cell_types))])

    offsets = [int(o) for o in offsets]
    n_ct = len(ref.cell_types)
    activity = np.zeros((n_ct, len(offsets)), dtype=np.float64)
    taus, strongest, decomp = [], [], []
    profiles_by_offset = {}

    if progress is not None:
        progress.update({"total": len(offsets), "done": 0, "stage": "scanning positions"})

    for j, off in enumerate(offsets):
        if progress is not None:
            progress["done"] = j
        pl = to_placements([{**element, "position": off, "uid": "scan"}])
        con = build_construct(seq, tss, pl, seed=0)
        pred = adapters.predict_celltype(model, con["sequence"], tss, cts, window)

        act = np.array([
            scoring.activity_from_profile(pred.profiles[i], pred.positions,
                                          activity_window, activity_method)
            for i in range(n_ct)])
        activity[:, j] = act

        m = scoring.celltype_metrics(list(pred.cell_types), act, target)
        taus.append(m["tau"])
        strongest.append(m["strongest_cell_type"])
        decomp.append(scoring.decompose_specificity_change(
            ref_act, act, list(pred.cell_types), target) if target else None)

        # keep a few full profiles so the UI can show representative traces
        if j % max(1, len(offsets) // 8) == 0 or j == len(offsets) - 1:
            profiles_by_offset[off] = np.round(pred.profiles, 4).tolist()

    target_idx = list(ref.cell_types).index(target) if target in ref.cell_types else None
    best_j = int(np.argmax(np.asarray(taus))) if taus else None

    result = {
        "mode": "mode1_position_scan",
        "model": model,
        "is_mock": adapter.is_mock,
        "scale": adapter.scale,
        "output_space": adapter.output_space,
        "element": element,
        "element_width": element_width(to_placements([element])[0]),
        "offsets": offsets,
        "cell_types": list(ref.cell_types),
        "cell_type_labels": ref.meta.get("cell_type_labels", {}),
        "activity": np.round(activity, 5).tolist(),
        "reference_activity": np.round(ref_act, 5).tolist(),
        "delta_activity": np.round(activity - ref_act[:, None], 5).tolist(),
        "tau": [round(t, 5) for t in taus],
        "strongest_cell_type": strongest,
        "target": target,
        "target_index": target_idx,
        "decomposition": decomp,
        "activity_window": list(activity_window),
        "activity_method": activity_method,
        "profile_positions": ref.positions.tolist(),
        "profiles_by_offset": profiles_by_offset,
        "background": {
            "hash": sequence_hash(seq),
            "length": len(seq),
            "tss_index": tss,
            "stats": background.get("achieved", {}),
            "seed": background.get("requested", {}).get("seed"),
        },
        "best_offset": offsets[best_j] if best_j is not None else None,
        "best_tau": taus[best_j] if best_j is not None else None,
        "caveats": [
            "Tau depends on which cell types are in the panel; comparisons are "
            "only valid across identical panels.",
            "Tau is computed on the activity scale stated in 'scale'. Changing "
            "activity_method or the scale changes the value.",
            "A single background can produce position effects that are "
            "properties of that background. Re-run across several backgrounds.",
        ],
    }
    return result


def position_scan_multi(backgrounds: list[dict], element: dict,
                        offsets: list[int], progress: dict | None = None,
                        **kwargs) -> dict:
    """Run :func:`position_scan` over several matched backgrounds."""
    runs = []
    for bi, bg in enumerate(backgrounds):
        if progress is not None:
            progress["detail"] = f"background {bi + 1}/{len(backgrounds)}"
        runs.append(position_scan(bg, element, offsets, progress=progress, **kwargs))
    act = np.array([r["activity"] for r in runs])          # (B, C, O)
    tau = np.array([r["tau"] for r in runs])               # (B, O)
    return {
        **runs[0],
        "mode": "mode1_position_scan_multi",
        "n_backgrounds": len(runs),
        "activity": np.round(act.mean(axis=0), 5).tolist(),
        "activity_sd": np.round(act.std(axis=0), 5).tolist(),
        "tau": np.round(tau.mean(axis=0), 5).tolist(),
        "tau_sd": np.round(tau.std(axis=0), 5).tolist(),
        "per_background": [{"hash": r["background"]["hash"],
                            "seed": r["background"]["seed"],
                            "tau": r["tau"],
                            "best_offset": r["best_offset"]} for r in runs],
    }


# ---------------------------------------------------------------------------
# Mode 2: motif interactions
# ---------------------------------------------------------------------------

def _placement(motif_id: str, position: int, strand: str, uid: str,
               seg_length: int = 60) -> dict:
    return {"element_id": motif_id, "position": int(position), "strand": strand,
            "uid": uid, "instance": "consensus", "seg_length": seg_length}


def pair_positions(motif_a: str, motif_b: str, anchor: int, spacing: int,
                   spacing_mode: str = "center") -> tuple[int, int]:
    """5' positions of A and B for a requested spacing.

    ``center``: centre-to-centre distance is ``spacing``.
    ``edge``:   gap between A's 3' end and B's 5' start is ``spacing``
                (``spacing = 0`` means immediately abutting).

    ``anchor`` is the 5' position of A relative to the TSS.
    """
    lib = get_library()
    wa = lib[motif_a].width if motif_a in lib else 10
    wb = lib[motif_b].width if motif_b in lib else 10
    pos_a = int(anchor)
    if spacing_mode == "edge":
        pos_b = pos_a + wa + int(spacing)
    else:
        centre_a = pos_a + wa / 2.0
        pos_b = int(round(centre_a + int(spacing) - wb / 2.0))
    return pos_a, pos_b


def pair_interaction(background: dict, motif_a: str, motif_b: str,
                     anchor: int = -60, spacing: int = 30,
                     spacing_mode: str = "center",
                     strand_a: str = "+", strand_b: str = "+",
                     models: list[str] | None = None,
                     window: tuple[int, int] = DEFAULT_WINDOW,
                     track: str = "plus",
                     residual_scale: str = "model") -> dict:
    """Four matched constructs, each model's profiles, and the residual.

    ``residual_scale``:
        ``"model"``  - compute I on the adapter's native output (default).
        ``"linear"`` - back-transform a log-space output to linear first, so I
                       tests departure from additivity in count-like units.
    """
    models = models or ["mock_puffin", "mock_puffind"]
    seq = background["sequence"]
    tss = background["tss_index"]

    pos_a, pos_b = pair_positions(motif_a, motif_b, anchor, spacing, spacing_mode)

    constructs = {
        "background": [],
        "A": [_placement(motif_a, pos_a, strand_a, "A")],
        "B": [_placement(motif_b, pos_b, strand_b, "B")],
        "AB": [_placement(motif_a, pos_a, strand_a, "A"),
               _placement(motif_b, pos_b, strand_b, "B")],
    }

    built = {}
    for key, pls in constructs.items():
        built[key] = build_construct(seq, tss, to_placements(pls), seed=0)

    # Do A and B overlap? If they do, the A+B construct does not contain both
    # intact sites and the residual is not an interaction measurement.
    wa = element_width(to_placements(constructs["A"])[0]) if constructs["A"] else 0
    wb = element_width(to_placements(constructs["B"])[0]) if constructs["B"] else 0
    overlap = max(0, min(pos_a + wa, pos_b + wb) - max(pos_a, pos_b))

    # Even when the two sites do not overlap, a gap smaller than the widest
    # motif in the library means some scan windows in the A+B construct span
    # the junction between them. Such a window sees a sequence that exists in
    # neither the A-only nor the B-only construct, so it can create (or
    # destroy) matches that the difference-in-differences then charges to
    # "interaction". This is a sequence-composition artefact, not model
    # non-additivity, and it is the main way this experiment design can be
    # made to lie. Two adjacent GC-boxes are the classic case: butting them
    # together builds a longer G-rich run with new internal matches.
    max_motif_width = max(m.width for m in get_library().values())
    gap = int(pos_b - (pos_a + wa)) if pos_b >= pos_a else int(pos_a - (pos_b + wb))
    junction_risk = bool(overlap == 0 and gap < max_motif_width)

    out_models = {}
    for mname in models:
        adapter = adapters.get_profile_adapter(mname)
        preds = {k: adapters.predict_profile(mname, built[k]["sequence"], tss, window)
                 for k in built}
        positions = preds["background"].positions

        def tr(key: str) -> np.ndarray:
            v = np.asarray(preds[key].tracks[track], dtype=np.float64)
            if residual_scale == "linear":
                return scoring.to_linear(v, adapter.output_space)
            return v

        f_bg, f_a, f_b, f_ab = tr("background"), tr("A"), tr("B"), tr("AB")
        residual = f_ab - f_a - f_b + f_bg

        # Effect sizes relative to the background, for the summary panel.
        e_a, e_b, e_ab = f_a - f_bg, f_b - f_bg, f_ab - f_bg
        peak_idx = int(np.argmax(np.abs(residual))) if residual.size else 0

        used_space = ("linear" if residual_scale == "linear" else adapter.output_space)
        out_models[mname] = {
            "label": adapter.label,
            "is_mock": adapter.is_mock,
            "scale": adapter.scale,
            "output_space": adapter.output_space,
            "scale_used": used_space,
            "residual_meaning": (
                "departure from additivity in linear (count-like) units"
                if used_space == "linear" else
                "departure from additivity on a log-like output, i.e. "
                "departure from a multiplicative combination in count units"),
            "positions": positions.tolist(),
            "profiles": {k: np.round(tr(k), 5).tolist() for k in
                         ("background", "A", "B", "AB")},
            "expected_additive": np.round(f_a + f_b - f_bg, 5).tolist(),
            "residual": np.round(residual, 5).tolist(),
            "summary": {
                "I_sum": float(np.round(residual.sum(), 5)),
                "I_sum_abs": float(np.round(np.abs(residual).sum(), 5)),
                "I_max_abs": float(np.round(np.abs(residual).max(), 5)) if residual.size else 0.0,
                "I_at_max_abs": float(np.round(residual[peak_idx], 5)) if residual.size else 0.0,
                "I_position_at_max_abs": int(positions[peak_idx]) if residual.size else 0,
                "effect_A_max": float(np.round(e_a.max(), 5)) if e_a.size else 0.0,
                "effect_B_max": float(np.round(e_b.max(), 5)) if e_b.size else 0.0,
                "effect_AB_max": float(np.round(e_ab.max(), 5)) if e_ab.size else 0.0,
                "effect_sum_max": float(np.round((e_a + e_b).max(), 5)) if e_a.size else 0.0,
                "relative_I": float(np.round(
                    np.abs(residual).sum() / max(np.abs(e_a).sum() + np.abs(e_b).sum(), 1e-9),
                    5)),
            },
        }

    return {
        "mode": "mode2_pair_interaction",
        "motif_a": motif_a, "motif_b": motif_b,
        "position_a": pos_a, "position_b": pos_b,
        "anchor": anchor, "spacing": spacing, "spacing_mode": spacing_mode,
        "strand_a": strand_a, "strand_b": strand_b,
        "track": track,
        "residual_scale": residual_scale,
        "overlap_bp": overlap,
        "gap_bp": gap,
        "junction_risk": junction_risk,
        "max_motif_width": max_motif_width,
        "overlap_warning": (
            "A and B overlap by %d bp, so the A+B construct does not contain "
            "two intact sites; the residual here is not an interaction "
            "measurement." % overlap) if overlap > 0 else None,
        "junction_warning": (
            "A and B are only %d bp apart, less than the widest motif in the "
            "library (%d bp). Scan windows spanning the junction see sequence "
            "present in neither single-element construct, so part of the "
            "residual here is a sequence-composition artefact rather than "
            "model non-additivity. An additive model's residual is a direct "
            "readout of how large that artefact is." % (gap, max_motif_width)
        ) if junction_risk else None,
        "sequences": {k: {"hash": sequence_hash(v["sequence"]),
                          "gc": v["stats"]["gc"],
                          "cpg_oe": v["stats"]["cpg_oe"]}
                      for k, v in built.items()},
        "background": {"hash": sequence_hash(seq), "length": len(seq),
                       "tss_index": tss,
                       "seed": background.get("requested", {}).get("seed")},
        "models": out_models,
    }


def pair_interaction_multi(backgrounds: list[dict], **kwargs) -> dict:
    """:func:`pair_interaction` repeated over matched backgrounds, averaged."""
    runs = [pair_interaction(bg, **kwargs) for bg in backgrounds]
    first = runs[0]
    merged = {k: v for k, v in first.items() if k != "models"}
    merged["mode"] = "mode2_pair_interaction_multi"
    merged["n_backgrounds"] = len(runs)
    merged["background_hashes"] = [r["background"]["hash"] for r in runs]

    models = {}
    for mname in first["models"]:
        res = np.array([r["models"][mname]["residual"] for r in runs])
        prof = {k: np.array([r["models"][mname]["profiles"][k] for r in runs])
                for k in ("background", "A", "B", "AB")}
        sums = {k: [r["models"][mname]["summary"][k] for r in runs]
                for k in first["models"][mname]["summary"]}

        summary = {k: float(np.round(np.mean(v), 5)) for k, v in sums.items()}
        summary_sd = {k: float(np.round(np.std(v), 5)) for k, v in sums.items()}

        # A *position* has no meaningful mean across backgrounds - averaging
        # "the peak was at -40 here and at -17 there" gives a position where
        # nothing happened in either. Read the location off the averaged
        # residual instead, and report how much the per-background locations
        # actually varied.
        mean_res = res.mean(axis=0)
        positions = np.asarray(first["models"][mname]["positions"])
        if mean_res.size:
            k = int(np.argmax(np.abs(mean_res)))
            summary["I_max_abs"] = float(np.round(np.abs(mean_res).max(), 5))
            summary["I_at_max_abs"] = float(np.round(mean_res[k], 5))
            summary["I_position_at_max_abs"] = int(positions[k])
            summary_sd["I_position_at_max_abs"] = float(
                np.round(np.std(sums["I_position_at_max_abs"]), 2))

        models[mname] = {
            **{k: v for k, v in first["models"][mname].items()
               if k not in ("residual", "profiles", "summary")},
            "residual": np.round(mean_res, 5).tolist(),
            "residual_sd": np.round(res.std(axis=0), 5).tolist(),
            "profiles": {k: np.round(v.mean(axis=0), 5).tolist() for k, v in prof.items()},
            "summary": summary,
            "summary_sd": summary_sd,
        }
    merged["models"] = models
    return merged


def pair_grid(backgrounds: list[dict], motif_ids: list[str] | None = None,
              anchor: int = -60, spacing: int = 30,
              spacing_mode: str = "center",
              strand_a: str = "+", strand_b: str = "+",
              models: list[str] | None = None,
              statistic: str = "I_sum_abs",
              window: tuple[int, int] = DEFAULT_WINDOW,
              track: str = "plus",
              residual_scale: str = "model",
              progress: dict | None = None) -> dict:
    """All pairs from the library, including same-motif (homotypic) pairs."""
    ids = motif_ids or library_order()
    models = models or ["mock_puffin", "mock_puffind"]
    n = len(ids)
    if progress is not None:
        progress.update({"total": n * (n + 1) // 2, "done": 0,
                         "stage": "pair grid"})

    grids = {m: np.full((n, n), np.nan) for m in models}
    grids_sd = {m: np.full((n, n), np.nan) for m in models}
    cells: list[dict] = []

    for i, a in enumerate(ids):
        for j, b in enumerate(ids):
            if j < i:
                continue  # symmetric in the pair term; fill the mirror below
            if progress is not None:
                progress["done"] += 1
                progress["detail"] = f"{a} x {b}"
            r = pair_interaction_multi(
                backgrounds, motif_a=a, motif_b=b, anchor=anchor,
                spacing=spacing, spacing_mode=spacing_mode,
                strand_a=strand_a, strand_b=strand_b, models=models,
                window=window, track=track, residual_scale=residual_scale)
            entry = {"a": a, "b": b, "overlap_bp": r["overlap_bp"],
                     "gap_bp": r.get("gap_bp"),
                     "junction_risk": r.get("junction_risk", False),
                     "values": {}}
            for m in models:
                v = r["models"][m]["summary"].get(statistic, float("nan"))
                sd = r["models"][m]["summary_sd"].get(statistic, float("nan"))
                grids[m][i, j] = grids[m][j, i] = v
                grids_sd[m][i, j] = grids_sd[m][j, i] = sd
                entry["values"][m] = {"value": v, "sd": sd}
            cells.append(entry)

    return {
        "mode": "mode2_pair_grid",
        "motif_ids": ids,
        "motif_labels": {k: get_library()[k].name for k in ids if k in get_library()},
        "statistic": statistic,
        "anchor": anchor, "spacing": spacing, "spacing_mode": spacing_mode,
        "strand_a": strand_a, "strand_b": strand_b,
        "track": track, "residual_scale": residual_scale,
        "n_backgrounds": len(backgrounds),
        "background_hashes": [sequence_hash(b["sequence"]) for b in backgrounds],
        "grids": {m: np.round(grids[m], 5).tolist() for m in models},
        "grids_sd": {m: np.round(grids_sd[m], 5).tolist() for m in models},
        "cells": cells,
        "note": ("The grid is symmetric because A and B are placed at fixed "
                 "anchor and spacing; swapping which motif is called A only "
                 "swaps the two positions, which is covered by the "
                 "orientation and spacing controls. Cells flagged "
                 "'junction_risk' sit close enough together that part of the "
                 "residual is a sequence-composition artefact; compare "
                 "against an additive model to size it."),
    }


def spacing_curve(backgrounds: list[dict], motif_a: str, motif_b: str,
                  spacings: list[int], anchor: int = -60,
                  spacing_mode: str = "center",
                  orientations: list[tuple[str, str]] | None = None,
                  models: list[str] | None = None,
                  statistic: str = "I_sum_abs",
                  window: tuple[int, int] = DEFAULT_WINDOW,
                  track: str = "plus",
                  residual_scale: str = "model",
                  progress: dict | None = None) -> dict:
    """Interaction summary as a function of spacing, per orientation."""
    models = models or ["mock_puffin", "mock_puffind"]
    orientations = orientations or [("+", "+")]
    spacings = [int(s) for s in spacings]
    if progress is not None:
        progress.update({"total": len(spacings) * len(orientations), "done": 0,
                         "stage": "spacing curve"})

    series: dict[str, dict[str, dict]] = {}
    for sa, sb in orientations:
        okey = f"{sa}{sb}"
        series[okey] = {m: {"value": [], "sd": []} for m in models}
        series[okey]["overlap_bp"] = []
        series[okey]["junction_risk"] = []
        series[okey]["gap_bp"] = []
        for s in spacings:
            if progress is not None:
                progress["done"] += 1
                progress["detail"] = f"{okey} spacing {s}"
            r = pair_interaction_multi(
                backgrounds, motif_a=motif_a, motif_b=motif_b, anchor=anchor,
                spacing=s, spacing_mode=spacing_mode, strand_a=sa, strand_b=sb,
                models=models, window=window, track=track,
                residual_scale=residual_scale)
            series[okey]["overlap_bp"].append(r["overlap_bp"])
            series[okey]["junction_risk"].append(r.get("junction_risk", False))
            series[okey]["gap_bp"].append(r.get("gap_bp"))
            for m in models:
                series[okey][m]["value"].append(r["models"][m]["summary"].get(statistic))
                series[okey][m]["sd"].append(r["models"][m]["summary_sd"].get(statistic))

    return {
        "mode": "mode2_spacing_curve",
        "motif_a": motif_a, "motif_b": motif_b,
        "spacings": spacings, "anchor": anchor, "spacing_mode": spacing_mode,
        "orientations": [f"{a}{b}" for a, b in orientations],
        "statistic": statistic, "models": models,
        "track": track, "residual_scale": residual_scale,
        "n_backgrounds": len(backgrounds),
        "series": series,
        "note": ("Spacings at which A and B overlap are flagged in "
                 "'overlap_bp'; the residual is not interpretable there. "
                 "Spacings flagged in 'junction_risk' are close enough that "
                 "scan windows span the junction, so part of the residual is "
                 "a sequence-composition artefact."),
    }
