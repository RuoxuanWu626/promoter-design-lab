"""Invariant checks for the promoter design lab.

    python3 selftest.py

These are the properties the two modes rely on for their conclusions to mean
anything.  The additivity check in particular exists because a subtle bug (a
top-N cut on the motif hit list) once broke it silently: the plots still looked
completely reasonable, and Puffin's interaction residual -- which is supposed to
be identically zero -- was 15.5.
"""

from __future__ import annotations

import sys

import numpy as np

import adapters
import experiments as E
import scoring
import sequence as S
from adapters.mock_common import motif_hits_by_motif
from motifs import get_library, null_threshold

FAILS: list[str] = []
CHECKS = [0]


def check(name: str, ok: bool, detail: str = "") -> None:
    CHECKS[0] += 1
    mark = "PASS" if ok else "FAIL"
    print(f"  [{mark}] {name}" + (f"  — {detail}" if detail else ""))
    if not ok:
        FAILS.append(name)


def section(title: str) -> None:
    print(f"\n{title}")


# ---------------------------------------------------------------------------
section("background generation")

bg = S.random_background(length=2001, gc=0.62, cpg_oe=0.80, seed=3)
check("GC content hits its target",
      abs(bg["achieved"]["gc"] - 0.62) < 0.02,
      f"asked 0.62, got {bg['achieved']['gc']:.3f}")
check("CpG o/e hits its target",
      abs(bg["achieved"]["cpg_oe"] - 0.80) < 0.06,
      f"asked 0.80, got {bg['achieved']['cpg_oe']:.3f}")
check("generation is deterministic for a seed",
      S.random_background(length=500, gc=0.5, cpg_oe=0.3, seed=11)["sequence"]
      == S.random_background(length=500, gc=0.5, cpg_oe=0.3, seed=11)["sequence"])
check("different seeds give different sequences",
      S.random_background(length=500, seed=1)["sequence"]
      != S.random_background(length=500, seed=2)["sequence"])

# ---------------------------------------------------------------------------
section("construct assembly")

base = S.random_background(length=1001, gc=0.45, cpg_oe=0.25, seed=5)
pl = S.to_placements([
    {"element_id": "tata", "position": -31, "uid": "a"},
    {"element_id": "inr", "position": -3, "uid": "b"},
    {"element_id": "cpg_segment", "position": -300, "seg_length": 80, "uid": "c"},
])
con = S.build_construct(base["sequence"], base["tss_index"], pl, seed=1)
check("placing elements preserves length",
      len(con["sequence"]) == len(base["sequence"]),
      f"{len(base['sequence'])} -> {len(con['sequence'])}")
check("element bases actually land at the requested position",
      con["sequence"][base["tss_index"] - 31: base["tss_index"] - 31 + 8]
      == get_library()["tata"].consensus_instance())
check("a minus-strand placement writes the reverse complement",
      S.build_construct(base["sequence"], base["tss_index"],
                        S.to_placements([{"element_id": "ets", "position": 0, "strand": "-"}]),
                        seed=1)["sequence"][base["tss_index"]: base["tss_index"] + 9]
      == S.reverse_complement(get_library()["ets"].consensus_instance()))
check("overlapping placements are reported",
      len(S.build_construct(base["sequence"], base["tss_index"], S.to_placements([
          {"element_id": "sp1", "position": 0, "uid": "x"},
          {"element_id": "sp1", "position": 4, "uid": "y"}]), seed=1)["overlaps"]) == 1)

# ---------------------------------------------------------------------------
section("motif thresholds")

lib = get_library()
hits_per_kb = 0.0
for mid in lib:
    n = sum((motif_hits_by_motif(
        S.random_background(length=1001, gc=0.45, cpg_oe=0.25, seed=200 + k)["sequence"],
        [mid])[mid]["occupancy"] > 0.5).sum() for k in range(4)) / 4
    hits_per_kb += n
check("chance matches across the library stay sparse",
      hits_per_kb < 20, f"{hits_per_kb:.1f} strong hits per kb (library total)")
check("a placed consensus is always recognised",
      all(motif_hits_by_motif(
          S.build_construct(base["sequence"], base["tss_index"],
                            S.to_placements([{"element_id": m, "position": -100}]),
                            seed=0)["sequence"], [m])[m]["occupancy"].max() > 0.5
          for m in lib))
check("thresholds are independent of the sequence being scanned",
      null_threshold(lib["sp1"]) == null_threshold(lib["sp1"]))

# ---------------------------------------------------------------------------
section("Mode 2 — the additivity control")

bgs = E.matched_backgrounds(n=3, length=1001, base_seed=900)

r = E.pair_interaction_multi(bgs, motif_a="tata", motif_b="ets",
                             anchor=-160, spacing=60)
pf = r["models"]["mock_puffin"]["summary"]["I_sum_abs"]
pd = r["models"]["mock_puffind"]["summary"]["I_sum_abs"]
check("well-separated sites are flagged neither overlapping nor junction-risk",
      r["overlap_bp"] == 0 and not r["junction_risk"], f"gap {r['gap_bp']} bp")
check("MockPuffin is exactly additive (I == 0)", pf == 0.0, f"Sigma|I| = {pf!r}")
check("MockPuffinD is not additive", pd > 0.1, f"Sigma|I| = {pd:.4f}")

grid = E.pair_grid(bgs[:2], anchor=-160, spacing=60)
gp = np.abs(np.array(grid["grids"]["mock_puffin"], dtype=float))
gd = np.abs(np.array(grid["grids"]["mock_puffind"], dtype=float))
check("additivity holds for every pair in the library",
      np.nanmax(gp) == 0.0, f"max |Sigma|I|| = {np.nanmax(gp)!r}")
check("the deep mock shows structure across pairs",
      np.nanmax(gd) > 1.0, f"max = {np.nanmax(gd):.3f}")
check("the grid includes same-motif pairs on the diagonal",
      all(np.isfinite(gd[i, i]) for i in range(gd.shape[0])))
check("the grid is symmetric", np.allclose(gd, gd.T, equal_nan=True))

close = E.pair_interaction(bgs[0], motif_a="sp1", motif_b="sp1",
                           anchor=-80, spacing=12)
check("near-adjacent sites are flagged as a junction artefact",
      close["junction_risk"] and close["junction_warning"],
      f"gap {close['gap_bp']} bp")
check("an additive model reports that artefact rather than hiding it",
      close["models"]["mock_puffin"]["summary"]["I_sum_abs"] > 0,
      "this is the artefact floor, and it is meant to be visible")

ov = E.pair_interaction(bgs[0], motif_a="sp1", motif_b="sp1",
                        anchor=-80, spacing=4)
check("overlapping sites are flagged", ov["overlap_bp"] > 0 and ov["overlap_warning"])

# residual scale
lin = E.pair_interaction(bgs[0], motif_a="tata", motif_b="ets", anchor=-160,
                         spacing=60, residual_scale="linear")
check("the residual scale is recorded with the result",
      lin["models"]["mock_puffind"]["scale_used"] == "linear")
check("an additive-in-log model is NOT additive in linear space",
      lin["models"]["mock_puffin"]["summary"]["I_sum_abs"] > 0,
      "expected: zero residual on a log output means multiplicative in counts")

# pair terms on/off
on = E.pair_interaction(bgs[0], motif_a="inr", motif_b="tata", anchor=-60,
                        spacing=28, models=["mock_puffind", "mock_puffind_nosat"])
check("the pair-terms-off control differs from the full deep mock",
      abs(on["models"]["mock_puffind"]["summary"]["I_sum_abs"]
          - on["models"]["mock_puffind_nosat"]["summary"]["I_sum_abs"]) > 1e-6,
      "lets saturation be separated from a genuine pair term")

# ---------------------------------------------------------------------------
section("Mode 1 — specificity")

check("tau is 0 when every cell type is equal",
      scoring.tau_specificity(np.ones(8)) == 0.0)
check("tau is 1 when only one cell type is active",
      abs(scoring.tau_specificity(np.array([1.0] + [0.0] * 7)) - 1.0) < 1e-12)
check("tau is 0 for all-zero input rather than NaN",
      scoring.tau_specificity(np.zeros(5)) == 0.0)

offs = list(range(-200, 81, 20))
scan = E.position_scan_multi(bgs, {"element_id": "u1", "strand": "+"}, offs,
                             target="SK-N-SH")
check("a position scan returns one activity per cell type per offset",
      np.array(scan["activity"]).shape == (len(scan["cell_types"]), len(offs)))
check("moving an element changes the cell-type panel",
      max(scan["tau"]) - min(scan["tau"]) > 0.01,
      f"tau spans {min(scan['tau']):.3f}–{max(scan['tau']):.3f}")
check("replicate spread is reported", "tau_sd" in scan and len(scan["tau_sd"]) == len(offs))

mechs = {d["mechanism"] for d in scan["decomposition"] if d}
check("off-target loss is reachable, not just target gain",
      any("off-target loss" in m for m in mechs), f"seen: {sorted(mechs)}")

gain = E.position_scan_multi(bgs, {"element_id": "ets", "strand": "+"}, offs,
                             target="K562")
gmechs = {d["mechanism"] for d in gain["decomposition"] if d}
check("target gain is reachable too",
      any("target gain" in m for m in gmechs), f"seen: {sorted(gmechs)}")

dec = scoring.decompose_specificity_change(
    np.array([1.0, 1.0, 1.0]), np.array([2.0, 1.0, 1.0]), ["a", "b", "c"], "a")
check("a pure target rise is labelled target gain", dec["mechanism"] == "target gain")
dec2 = scoring.decompose_specificity_change(
    np.array([1.0, 1.0, 1.0]), np.array([1.0, 0.3, 0.3]), ["a", "b", "c"], "a")
check("a pure off-target fall is labelled off-target loss",
      dec2["mechanism"] == "off-target loss")

# ---------------------------------------------------------------------------
section("scoring")

pred = adapters.predict_profile("mock_puffin", con["sequence"], con["tss_index"], (-500, 500))
met = scoring.profile_metrics(pred.positions, pred.tracks["plus"], pred.tracks["minus"])
check("the designed TATA+Inr peak lands at the TSS",
      abs(met["peak_position"]) <= 5, f"peak at {met['peak_position']:+d} bp")

plus_only = S.build_construct(base["sequence"], base["tss_index"], S.to_placements(
    [{"element_id": "tata", "position": -31}, {"element_id": "inr", "position": -3}]), seed=0)
pp = adapters.predict_profile("mock_puffin", plus_only["sequence"], base["tss_index"], (-300, 300))
dm = scoring.profile_metrics(pp.positions, pp.tracks["plus"], pp.tracks["minus"])
check("directionality detects a plus-strand-only design",
      dm["directionality"] > 0.15, f"{dm['directionality']:.3f}")

ws = scoring.weighted_score({"peak_height": 4.0, "tss_offset_abs": 0.0},
                            {"peak_height": 1.0, "tss_offset_abs": 1.0})
check("a perfect design scores the full weight", abs(ws["score"] - 2.0) < 1e-9)
check("the score decomposes into per-metric contributions",
      set(ws["contributions"]) == {"peak_height", "tss_offset_abs"})

# ---------------------------------------------------------------------------
section("caching")

adapters.clear_cache()
a = adapters.predict_profile("mock_puffind", con["sequence"], con["tss_index"], (-200, 200))
b = adapters.predict_profile("mock_puffind", con["sequence"], con["tss_index"], (-200, 200))
check("repeated prediction is bit-identical",
      np.array_equal(a.tracks["plus"], b.tracks["plus"]))
check("the cache is actually used", adapters.cache_stats()["entries"] >= 1)

# ---------------------------------------------------------------------------
section("adapters")

info = adapters.list_models()
check("every adapter declares an output space",
      all(m["output_space"] in ("linear", "log")
          for m in info["profile"] + info["celltype"]))
check("every adapter declares whether it is a mock",
      all(isinstance(m["is_mock"], bool) for m in info["profile"] + info["celltype"]))
check("unavailable real adapters explain why",
      all(m["unavailable_reason"] for m in info["profile"] + info["celltype"]
          if not m["available"]))

# ---------------------------------------------------------------------------
print(f"\n{'=' * 60}")
if FAILS:
    print(f"{len(FAILS)} of {CHECKS[0]} checks FAILED:")
    for f in FAILS:
        print(f"  - {f}")
    sys.exit(1)
print(f"all {CHECKS[0]} checks passed")
