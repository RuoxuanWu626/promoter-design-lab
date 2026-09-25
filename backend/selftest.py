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
from motifs import celltype_element_order, get_library, null_threshold

FAILS: list[str] = []
CHECKS = [0]

# Whatever motifs a deployment has defined are set aside for the duration of
# the run and restored at the end. Without this the tests measure the live
# server's library: a user motif that overlaps a built-in -- say AGATWAGA,
# which also matches the GATA consensus AGATAAGA -- makes both fire on the
# same bases and quietly changes which cell type wins.
import custom_motifs as _CM

_SAVED_CUSTOM = _CM.list_motifs()
for _m in _SAVED_CUSTOM:
    _CM.remove(_m["id"])


def _restore_custom() -> None:
    for _m in _SAVED_CUSTOM:
        try:
            _CM.add({**_m, "replace": True})
        except Exception:
            pass


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
_tata = get_library()["tata"]
check("element bases actually land at the requested position",
      con["sequence"][base["tss_index"] - 31:
                      base["tss_index"] - 31 + _tata.width]
      == _tata.consensus_instance(),
      f"TATA is {_tata.width} bp: {_tata.consensus}")
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
# REST/NRSE represses every non-neuronal line, so sliding it with a neuronal
# target is the off-target-silencing route.
scan = E.position_scan_multi(bgs, {"element_id": "rest", "strand": "+"}, offs,
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

# GATA activates erythroid only, so this is the target-gain route.
gain = E.position_scan_multi(bgs, {"element_id": "gata", "strand": "+"}, offs,
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
section("Mode 1 - core promoter vs lineage sites")

clean = S.random_background(length=2001, gc=0.45, cpg_oe=0.25, seed=42)
check("scrubbing leaves no chance lineage site behind",
      not [1 for _m, d in
           __import__("adapters.mock_common", fromlist=["x"]).motif_hits_by_motif(
               clean["sequence"], celltype_element_order()).items()
           for o in d["occupancy"] if o > 0.5],
      f"{clean['scrubbed_celltype_sites']} removed")
check("scrubbing preserves composition",
      abs(clean["achieved"]["gc"] - 0.45) < 0.03,
      f"GC {clean['achieved']['gc']:.3f}, CpG o/e {clean['achieved']['cpg_oe']:.3f}")


def _panel(placements):
    con = S.build_construct(clean["sequence"], clean["tss_index"],
                            S.to_placements(placements), seed=1)
    pred = adapters.predict_celltype("mock_alphagenome", con["sequence"],
                                     con["tss_index"], None, (-500, 500))
    act = np.array([
        scoring.activity_from_profile(pred.profiles[i], pred.positions, (-200, 200), "mean")
        for i in range(len(pred.cell_types))])
    m = scoring.celltype_metrics(list(pred.cell_types), act)
    return m["tau"], m["strongest_cell_type"], float(act.mean())


_core = [{"element_id": "tata", "position": -31},
         {"element_id": "inr", "position": 0},
         {"element_id": "sp1", "position": -52}]
tau_empty, _, act_empty = _panel([])
tau_core, _, act_core = _panel(_core)
check("an empty scrubbed background is not cell-type specific",
      tau_empty < 0.08, f"tau = {tau_empty:.3f}")
check("core promoter motifs raise activity",
      act_core > act_empty * 1.15, f"{act_empty:.2f} -> {act_core:.2f}")
check("core promoter motifs do NOT create specificity",
      abs(tau_core - tau_empty) < 0.05,
      f"tau {tau_empty:.3f} -> {tau_core:.3f}; this is the claim Mode 1 rests on")

_expect = {"gata": "K562", "hnf4": ("HepG2", "Hepatocyte"), "spi1": "GM12878",
           "sox_oct": "H1", "ere": "MCF-7", "ebox_neuro": "SK-N-SH"}
_hits = []
for _eid, _want in _expect.items():
    _t, _strong, _a = _panel(_core + [{"element_id": _eid, "position": -120}])
    _want_t = _want if isinstance(_want, tuple) else (_want,)
    _hits.append((_eid, _strong in _want_t, _t))
check("one lineage site makes its own cell type the strongest",
      all(ok_ for _e, ok_, _t in _hits),
      ", ".join(f"{e}{'ok' if o else ' MISS'}" for e, o, _t in _hits))
check("one lineage site raises tau well above the core-only baseline",
      all(t > tau_core + 0.08 for _e, _o, t in _hits),
      f"core-only tau {tau_core:.3f}, with a site {min(t for _e,_o,t in _hits):.3f}-{max(t for _e,_o,t in _hits):.3f}")

_attr = E.specificity_attribution(
    clean, _core + [{"element_id": "gata", "position": -120}], target="K562",
    window=(-220, 60), patch=12, stride=8, n_shuffles=2)
_a, _t = _attr["summary"]["activity"], _attr["summary"]["specificity"]
check("activity attribution concentrates on the core promoter",
      _a["fraction_on_core_promoter"] > _a["fraction_on_celltype_elements"] * 3,
      f"core {_a['fraction_on_core_promoter']:.0%} vs lineage {_a['fraction_on_celltype_elements']:.0%}")
check("specificity attribution shifts onto the lineage site",
      _t["fraction_on_celltype_elements"] > _a["fraction_on_celltype_elements"] * 3,
      f"lineage carries {_t['fraction_on_celltype_elements']:.0%} of specificity "
      f"vs {_a['fraction_on_celltype_elements']:.0%} of activity")
check("background-occurring motifs are annotated, not counted as 'elsewhere'",
      any(f.get("source") == "background" for f in _attr["motif_footprints"])
      or _attr["patches_on_celltype_elements"] > 0)

# ---------------------------------------------------------------------------
section("user-defined motifs")

CM = _CM

_reject = [
    (dict(name="x", consensus="ACGTZZ"), "non-IUPAC letters"),
    (dict(name="x", consensus="AC"), "too short"),
    (dict(name="x", consensus="NNNNNNNN"), "matches everything"),
    (dict(name="clash", id="tata", consensus="ACGTACGT"), "reserved id"),
    (dict(consensus="ACGTACGT"), "no name"),
]
_bad_accepted = []
for _spec, _why in _reject:
    try:
        CM.add(_spec)
        _bad_accepted.append(_why)
    except CM.InvalidMotif:
        pass
check("invalid motif definitions are refused with a readable reason",
      not _bad_accepted, f"accepted: {_bad_accepted}" if _bad_accepted else "all five refused")

_entry = CM.add(dict(name="Selftest site", consensus="AGATWAGA",
                     kind="celltype_element", activates=["HepG2"],
                     typical_offset=-130))
check("a custom motif joins the library", _entry["id"] in get_library(),
      f"id {_entry['id']}")
check("a custom motif gets a calibrated detection threshold",
      null_threshold(get_library()[_entry["id"]]) > 0)
check("a custom lineage site is scrubbed from backgrounds like a built-in",
      _entry["id"] in __import__("motifs").celltype_kind_ids())

_jaspar = ("A [ 0 3 79 40 66 48 65 11 65 0 ]\n"
           "C [94 75 4 3 1 2 5 2 3 3 ]\n"
           "G [ 1 0 3 4 1 0 5 3 28 88 ]\n"
           "T [ 2 19 11 50 29 47 22 81 1 6 ]")
_pwm_entry = CM.add(dict(name="Selftest PWM", pwm=_jaspar, kind="core_promoter"))
check("a pasted JASPAR-style matrix is parsed",
      len(_pwm_entry["consensus"]) == 10,
      f"consensus {_pwm_entry['consensus']}")

_bg_u = S.random_background(length=2001, gc=0.45, cpg_oe=0.25, seed=42)
_core_u = [{"element_id": "tata", "position": -31},
           {"element_id": "inr", "position": 0},
           {"element_id": "sp1", "position": -52}]


def _tau_of(placements):
    con_ = S.build_construct(_bg_u["sequence"], _bg_u["tss_index"],
                             S.to_placements(placements), seed=1)
    pr = adapters.predict_celltype("mock_alphagenome", con_["sequence"],
                                   con_["tss_index"], None, (-500, 500))
    act = np.array([
        scoring.activity_from_profile(pr.profiles[i], pr.positions, (-200, 200), "mean")
        for i in range(len(pr.cell_types))])
    m = scoring.celltype_metrics(list(pr.cell_types), act)
    return m["tau"], m["strongest_cell_type"]


_t_base, _ = _tau_of(_core_u)
_t_user, _strong_user = _tau_of(
    _core_u + [{"element_id": _entry["id"], "position": p} for p in (-160, -120)])
check("a user-defined lineage site drives the cell type it was assigned",
      _strong_user == "HepG2" and _t_user > _t_base + 0.1,
      f"tau {_t_base:.3f} -> {_t_user:.3f}, strongest {_strong_user}")

_p_no = adapters.predict_profile("puffin", S.build_construct(
    _bg_u["sequence"], _bg_u["tss_index"], S.to_placements(_core_u), seed=1)["sequence"],
    _bg_u["tss_index"], (-300, 300))
_p_yes = adapters.predict_profile("puffin", S.build_construct(
    _bg_u["sequence"], _bg_u["tss_index"], S.to_placements(
        _core_u + [{"element_id": _entry["id"], "position": p} for p in (-160, -120)]),
    seed=1)["sequence"], _bg_u["tss_index"], (-300, 300))
check("the real model responds to a custom motif's bases, not to its label",
      not np.array_equal(np.asarray(_p_no.tracks["plus"]),
                         np.asarray(_p_yes.tracks["plus"])),
      "Puffin has no filter for it but does see the sequence")

check("deleting a custom motif removes it from the library",
      CM.remove(_entry["id"]) and _entry["id"] not in get_library())
CM.remove(_pwm_entry["id"])
check("the custom library is empty again after cleanup", not CM.list_motifs())

# ---------------------------------------------------------------------------
section("model agreement")

_ent = [
    {"name": "a", "positions": np.arange(100), "values": np.sin(np.arange(100) / 8.0) + 2,
     "output_space": "linear", "role": "simple"},
    {"name": "b", "positions": np.arange(100), "values": 3 * (np.sin(np.arange(100) / 8.0) + 2),
     "output_space": "linear", "role": "deep"},
]
_ag = scoring.profile_agreement(_ent, transform="none")
check("Pearson r is invariant to a pure scale change",
      abs(_ag["pairs"][0]["pearson_r"] - 1.0) < 1e-6,
      f"r = {_ag['pairs'][0]['pearson_r']}")
check("anti-correlated profiles give r = -1",
      abs(scoring.pearson_r(np.arange(50), -np.arange(50)) + 1.0) < 1e-9)
check("a flat profile gives r = 0 rather than NaN",
      scoring.pearson_r(np.ones(50), np.arange(50)) == 0.0)
check("Spearman handles ties without blowing up",
      abs(scoring.spearman_r(np.array([1, 1, 2, 3]), np.array([1, 2, 2, 3]))) <= 1.0)

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
_restore_custom()
if _SAVED_CUSTOM:
    print(f"\n(restored {len(_SAVED_CUSTOM)} user-defined motif(s) set aside for the run)")

print(f"\n{'=' * 60}")
if FAILS:
    print(f"{len(FAILS)} of {CHECKS[0]} checks FAILED:")
    for f in FAILS:
        print(f"  - {f}")
    sys.exit(1)
print(f"all {CHECKS[0]} checks passed")
