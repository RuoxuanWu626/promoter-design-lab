"""Puffin-style 10-motif library.

IMPORTANT / PROVENANCE
----------------------
The consensus strings and position-weight matrices below are **placeholders**.
They are hand-written approximations of the ten promoter motifs that the Puffin
model reports, written so the interface and the experiment machinery can be
built and tested end to end.  They are NOT Puffin's learned motif weights.

To use the real library, drop Puffin's learned PWMs into
``data/puffin_motifs.npz`` (one ``(4, W)`` array per motif id, rows ordered
A, C, G, T, columns 5'->3') and they will be loaded in preference to these.
See ``docs/MODEL_ADAPTERS.md``.

Everything downstream reads motifs through :func:`get_library`, so swapping the
source changes nothing else.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field, asdict

import numpy as np

ALPHABET = "ACGT"
BASE_INDEX = {b: i for i, b in enumerate(ALPHABET)}

# IUPAC degenerate code -> set of concrete bases.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def consensus_to_pwm(consensus: str, conservation: float = 0.85) -> np.ndarray:
    """Build a (4, W) probability matrix from an IUPAC consensus string.

    ``conservation`` is the probability mass given to the allowed base(s) at
    each column; the remainder is spread uniformly over all four bases.  This is
    a crude stand-in for a real learned PWM and is labelled as such in the UI.
    """
    width = len(consensus)
    pwm = np.full((4, width), (1.0 - conservation) / 4.0, dtype=np.float64)
    for j, code in enumerate(consensus.upper()):
        allowed = IUPAC.get(code, "ACGT")
        share = conservation / len(allowed)
        for base in allowed:
            pwm[BASE_INDEX[base], j] += share
    return pwm / pwm.sum(axis=0, keepdims=True)


@dataclass
class Motif:
    """One entry in the motif library.

    Attributes prefixed ``mock_`` drive the *mock* prediction adapters only.
    They encode a deliberately simple, transparent story about what each motif
    does to an initiation profile.  Real adapters ignore them entirely.
    """

    id: str
    name: str
    aliases: list[str]
    consensus: str
    # Where this motif is typically found, as a hint shown in the UI
    # (distance of the motif's 5' end from the TSS, in bp).
    typical_offset: int
    strand_specific: bool
    color: str
    notes: str

    # --- mock-model parameters (see adapters/mock_common.py) ---------------
    mock_amplitude: float = 1.0      # peak height contributed, arbitrary units
    mock_peak_shift: int = 0         # bp from motif 5' end to the peak it drives
    mock_width: float = 12.0         # gaussian sigma of the contributed peak, bp
    mock_broad_fraction: float = 0.0  # share of the effect that is broad, not sharp
    mock_strand_asymmetry: float = 1.0  # multiplier applied on the minus strand

    # "core_promoter" for Puffin's ten filters, "celltype_element" for the
    # lineage TF sites in celltype_elements.py. They are kept apart because
    # they answer different questions: core promoter motifs set how much
    # initiation there is, lineage sites set which cell types show it.
    kind: str = "core_promoter"
    factors: str = ""
    activates: list = field(default_factory=list)
    represses: list = field(default_factory=list)

    # Filled in when a real checkpoint has been extracted (see
    # extract_puffin_motifs.py); None means the entry is still a placeholder.
    source: str = "placeholder"
    information_bits: float | None = None
    effect_profile: list | None = None
    effect_offsets: list | None = None
    filter_index: int | None = None
    filter_strand: str | None = None

    pwm: np.ndarray | None = field(default=None, repr=False)

    @property
    def width(self) -> int:
        return len(self.consensus)

    def to_json(self) -> dict:
        d = asdict(self)
        d.pop("pwm", None)
        d["width"] = self.width
        d["pwm"] = self.get_pwm().round(4).tolist()
        return d

    def get_pwm(self) -> np.ndarray:
        if self.pwm is None:
            self.pwm = consensus_to_pwm(self.consensus)
        return self.pwm

    def sample(self, rng: np.random.Generator, strand: str = "+") -> str:
        """Draw a concrete instance of this motif from its PWM."""
        pwm = self.get_pwm()
        idx = [rng.choice(4, p=pwm[:, j]) for j in range(pwm.shape[1])]
        seq = "".join(ALPHABET[i] for i in idx)
        return seq if strand == "+" else reverse_complement(seq)

    def consensus_instance(self, strand: str = "+") -> str:
        """Most-likely instance (the argmax of each PWM column)."""
        pwm = self.get_pwm()
        seq = "".join(ALPHABET[int(np.argmax(pwm[:, j]))] for j in range(pwm.shape[1]))
        return seq if strand == "+" else reverse_complement(seq)


# ---------------------------------------------------------------------------
# The ten motifs.  Consensus strings are placeholders (see module docstring).
# ---------------------------------------------------------------------------

_LIBRARY_SPEC: list[Motif] = [
    Motif(
        id="tata", name="TATA-box", aliases=["TBP"],
        consensus="TATAWAWR", typical_offset=-31, strand_specific=True,
        color="#e4572e",
        notes="Sharp, focused initiation. Strongly position-dependent: "
              "classically ~25-31 bp upstream of the TSS.",
        mock_amplitude=1.60, mock_peak_shift=30, mock_width=5.0,
        mock_broad_fraction=0.05, mock_strand_asymmetry=0.35,
    ),
    Motif(
        id="inr", name="Initiator (Inr)", aliases=["YYANWYY", "BBCABW"],
        consensus="YYCANTYY", typical_offset=-3, strand_specific=True,
        color="#f3a712",
        notes="Overlaps the TSS itself; the A of the CA core is position +1. "
              "Written here in a more informative form than the minimal "
              "YYANWYY: at ~6 bits that shorter version matches random "
              "sequence roughly every 130 bp, which would swamp the profile.",
        mock_amplitude=1.05, mock_peak_shift=3, mock_width=3.5,
        mock_broad_fraction=0.02, mock_strand_asymmetry=0.25,
    ),
    Motif(
        id="sp1", name="GC-box (SP1)", aliases=["SP1", "KLF"],
        consensus="GGGGCGGGGC", typical_offset=-65, strand_specific=False,
        color="#3a86ff",
        notes="Common in CpG-island promoters; tends to support broader, "
              "dispersed initiation.",
        mock_amplitude=0.85, mock_peak_shift=55, mock_width=26.0,
        mock_broad_fraction=0.55, mock_strand_asymmetry=0.85,
    ),
    Motif(
        id="nfy", name="NF-Y (CCAAT)", aliases=["CCAAT-box", "NFYA/B"],
        consensus="RRCCAATCRG", typical_offset=-85, strand_specific=False,
        color="#8338ec",
        notes="Upstream activator; in many promoters its effect is strongest "
              "around 60-100 bp upstream.",
        mock_amplitude=0.95, mock_peak_shift=80, mock_width=16.0,
        mock_broad_fraction=0.30, mock_strand_asymmetry=0.90,
    ),
    Motif(
        id="nrf1", name="NRF1", aliases=["Alpha-PAL"],
        consensus="GCGCATGCGC", typical_offset=-55, strand_specific=False,
        color="#06d6a0",
        notes="Palindromic; enriched in CpG-island promoters of "
              "housekeeping-like genes.",
        mock_amplitude=0.90, mock_peak_shift=48, mock_width=18.0,
        mock_broad_fraction=0.40, mock_strand_asymmetry=1.00,
    ),
    Motif(
        id="ets", name="ETS", aliases=["GABPA", "ELK1", "GGAA-core"],
        consensus="ACCGGAAGT", typical_offset=-45, strand_specific=True,
        color="#ef476f",
        notes="GGAA core. Frequently found in a narrow window upstream of "
              "CpG-island TSSs.",
        mock_amplitude=1.00, mock_peak_shift=42, mock_width=14.0,
        mock_broad_fraction=0.25, mock_strand_asymmetry=0.70,
    ),
    Motif(
        id="creb", name="CREB / ATF (CRE)", aliases=["ATF1", "TGACGTCA"],
        consensus="TGACGTCA", typical_offset=-70, strand_specific=False,
        color="#118ab2",
        notes="Palindromic CRE. Contains a CpG; signal-responsive in many "
              "cell types.",
        mock_amplitude=0.88, mock_peak_shift=62, mock_width=17.0,
        mock_broad_fraction=0.35, mock_strand_asymmetry=1.00,
    ),
    Motif(
        id="yy1", name="YY1", aliases=["Delta", "NF-E1"],
        consensus="AAWATGGCGGC", typical_offset=-10, strand_specific=True,
        color="#ffd166",
        notes="Often sits just downstream of, or across, the TSS and helps "
              "position initiation.",
        mock_amplitude=1.05, mock_peak_shift=-6, mock_width=7.0,
        mock_broad_fraction=0.10, mock_strand_asymmetry=0.45,
    ),
    Motif(
        id="znf143", name="ZNF143", aliases=["SBS", "SPH-motif", "THAP11/HCFC1"],
        consensus="TTCCCAGAATG", typical_offset=-90, strand_specific=True,
        color="#9d4edd",
        notes="Upstream of many promoters; associated with promoter-anchored "
              "chromatin contacts.",
        mock_amplitude=0.80, mock_peak_shift=85, mock_width=20.0,
        mock_broad_fraction=0.35, mock_strand_asymmetry=0.75,
    ),
    Motif(
        id="u1", name="U1 snRNP (5' splice site)", aliases=["5'SS", "GTAAGT"],
        consensus="MAGGTRAGT", typical_offset=120, strand_specific=True,
        color="#2a9d8f",
        notes="Downstream of the TSS. Included in the library because "
              "co-transcriptional U1 binding tracks with initiation signal.",
        mock_amplitude=0.55, mock_peak_shift=-95, mock_width=30.0,
        mock_broad_fraction=0.60, mock_strand_asymmetry=0.40,
    ),
]

# A non-motif element the UI can also place: a stretch of CpG-rich sequence.
CPG_SEGMENT_ID = "cpg_segment"

# A literal stretch of sequence the user supplies, so Mode 1 can slide an
# arbitrary segment (not just a library motif) relative to the TSS.
CUSTOM_ID = "custom"


_cache: dict[str, Motif] | None = None


def get_library(path: str | None = None) -> dict[str, Motif]:
    """Return {motif_id: Motif}, loading real PWMs from disk when available."""
    global _cache
    if _cache is not None:
        return _cache

    lib = {m.id: m for m in _LIBRARY_SPEC}

    # Imported here rather than at module scope to avoid a circular import;
    # celltype_elements builds Motif objects.
    from celltype_elements import build_elements
    lib.update(build_elements())

    if path is None:
        here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        path = os.path.join(here, "data", "puffin_motifs.npz")
    if os.path.exists(path):
        try:
            blob = np.load(path)
            for mid, motif in lib.items():
                if mid in blob:
                    arr = np.asarray(blob[mid], dtype=np.float64)
                    if arr.shape[0] == 4:
                        motif.pwm = arr / arr.sum(axis=0, keepdims=True)
                        motif.consensus = motif.consensus_instance()
                        motif.source = "puffin"
        except Exception as exc:  # pragma: no cover - diagnostics only
            print(f"[motifs] could not load {path}: {exc}")

        # Metadata written alongside the PWMs by extract_puffin_motifs.py:
        # real names, colours, and positions read off the model's own deconv
        # kernels. These replace the hand-written placeholders, so nothing in
        # the library is guessed once a checkpoint has been extracted.
        meta_path = os.path.splitext(path)[0] + ".json"
        if os.path.exists(meta_path):
            try:
                with open(meta_path) as fh:
                    meta = json.load(fh)
                for mid, m in (meta.get("motifs") or {}).items():
                    motif = lib.get(mid)
                    if motif is None:
                        continue
                    motif.name = m.get("puffin_name", motif.name)
                    motif.color = m.get("color", motif.color)
                    motif.typical_offset = int(m.get("typical_offset", motif.typical_offset))
                    motif.mock_peak_shift = int(m.get("peak_shift", motif.mock_peak_shift))
                    motif.information_bits = m.get("mean_information_bits")
                    motif.effect_profile = m.get("effect_profile")
                    motif.effect_offsets = m.get("effect_offsets")
                    motif.filter_index = m.get("filter_index")
                    motif.filter_strand = m.get("filter_strand")
                    motif.source = "puffin"
                    motif.notes = (
                        f"Learned Puffin filter {m.get('filter_index')} "
                        f"(strand {m.get('filter_strand')}), trimmed to its "
                        f"informative core. Its preferred position "
                        f"({motif.typical_offset:+d} bp) is read off the "
                        f"model's deconv kernel, not assumed."
                    )
                _cache = lib
                return lib
            except Exception as exc:  # pragma: no cover - diagnostics only
                print(f"[motifs] could not load {meta_path}: {exc}")

    _cache = lib
    return lib


def library_order() -> list[str]:
    """Puffin's ten core promoter motifs, in library order."""
    return [m.id for m in _LIBRARY_SPEC]


def celltype_element_order() -> list[str]:
    """The lineage TF sites, which are a separate library (see celltype_elements)."""
    from celltype_elements import CELLTYPE_ELEMENT_IDS
    return list(CELLTYPE_ELEMENT_IDS)


def core_library() -> dict:
    """Only the core promoter motifs."""
    lib = get_library()
    return {k: lib[k] for k in library_order() if k in lib}


def pwm_log_odds(pwm: np.ndarray, background: np.ndarray | None = None) -> np.ndarray:
    if background is None:
        background = np.full(4, 0.25)
    return np.log2((pwm + 1e-6) / background[:, None])


def scan_sequence(seq: str, motif: Motif, both_strands: bool = True) -> np.ndarray:
    """Log-odds score of ``motif`` at every start offset in ``seq``.

    Returns a ``(2, L)`` array (plus strand, minus strand); positions where the
    motif would run off the end are ``-inf``.
    """
    L = len(seq)
    w = motif.width
    idx = np.full(L, -1, dtype=np.int64)
    for i, ch in enumerate(seq):
        idx[i] = BASE_INDEX.get(ch, -1)

    lo = pwm_log_odds(motif.get_pwm())
    out = np.full((2, L), -np.inf)
    if L < w:
        return out

    valid_starts = L - w + 1
    windows = np.lib.stride_tricks.sliding_window_view(idx, w)[:valid_starts]
    ok = windows >= 0
    safe = np.where(ok, windows, 0)

    cols = np.arange(w)
    plus = np.where(ok, lo[safe, cols[None, :]], 0.0).sum(axis=1)
    out[0, :valid_starts] = plus

    if both_strands:
        lo_rc = lo[::-1, ::-1]
        minus = np.where(ok, lo_rc[safe, cols[None, :]], 0.0).sum(axis=1)
        out[1, :valid_starts] = minus
    return out


def max_possible_score(motif: Motif) -> float:
    return float(pwm_log_odds(motif.get_pwm()).max(axis=0).sum())


_NULL_CACHE: dict[tuple, float] = {}

# Base composition the null distribution is calibrated against. This is a fixed
# SETTING, not something derived from the sequence being scanned -- deriving it
# per sequence would make a motif's threshold depend on what else is in the
# construct, which would quietly destroy the additivity that Mode 2 relies on.
NULL_BG_GC = 0.45


def null_background(gc: float = NULL_BG_GC) -> np.ndarray:
    return np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])


def null_threshold(motif: Motif, p_value: float = 1e-4,
                   bin_width: float = 0.02,
                   bg_gc: float = NULL_BG_GC) -> float:
    """Score above which a uniform random k-mer lands with probability <= p.

    The exact null distribution of the log-odds sum is obtained by convolving
    the four-value distribution of each column, so this is not a simulation and
    not a normal approximation.

    Why this matters: a fixed "fraction of the maximum score" threshold is not
    comparable across motifs, because motifs carry very different amounts of
    information. An 18-bit motif like the GC-box essentially never matches by
    chance, while a 6-bit motif matches every few hundred bases. Thresholding
    each motif at the same *false-positive rate* puts them on the same footing,
    and it keeps working if the placeholder PWMs are swapped for real ones with
    different information content.

    Choosing ``p_value``: the whole library is scanned on both strands, so a
    1 kb construct is about 2 x 10^4 tests. A rate of 1e-4 therefore leaves a
    couple of chance matches per kb across the whole library -- enough texture
    to look like real sequence, few enough that a placed motif stands out.

    The score being thresholded is the log-odds against a uniform background
    (so scores stay interpretable as bits), while the null distribution is
    taken under ``bg_gc``, which is what actually controls the hit rate.
    """
    key = (motif.id, round(p_value, 12), bin_width, round(bg_gc, 4))
    if key in _NULL_CACHE:
        return _NULL_CACHE[key]

    lo = pwm_log_odds(motif.get_pwm())
    q = np.rint(lo / bin_width).astype(np.int64)      # (4, W) integer bins

    freqs = null_background(bg_gc)
    dist = np.ones(1, dtype=np.float64)
    base = 0
    for j in range(q.shape[1]):
        col_q = q[:, j]
        lo_j = int(col_q.min())
        col = np.zeros(int(col_q.max()) - lo_j + 1, dtype=np.float64)
        for b in range(4):
            col[int(col_q[b]) - lo_j] += freqs[b]
        dist = np.convolve(dist, col)
        base += lo_j

    tail = np.cumsum(dist[::-1])[::-1]                # P(score >= bin)
    idx = np.nonzero(tail <= p_value)[0]
    thr = ((base + int(idx[0])) * bin_width if idx.size
           else float(lo.max(axis=0).sum()) + bin_width)
    _NULL_CACHE[key] = float(thr)
    return float(thr)


def min_possible_score(motif: Motif) -> float:
    return float(pwm_log_odds(motif.get_pwm()).min(axis=0).sum())
