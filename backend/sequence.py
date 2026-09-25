"""Background generation and construct assembly.

Coordinate convention used everywhere in this project
-----------------------------------------------------
* A construct is a string of length ``L`` plus an integer ``tss_index``.
* Position ``p`` **relative to the TSS** maps to string index
  ``tss_index + p``.  So ``p = 0`` is the TSS base itself, ``p < 0`` is
  upstream, ``p > 0`` is downstream.
* An element placed at position ``p`` has its **5' end (on the plus strand of
  the construct) at relative position p**, and occupies ``p .. p + width - 1``.

Placing an element **overwrites** bases in place, so the construct length never
changes.  That is what makes the four Mode-2 constructs (background / A / B /
A+B) length-matched and directly comparable.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, asdict

import numpy as np

from motifs import (
    ALPHABET,
    CPG_SEGMENT_ID,
    CUSTOM_ID,
    Motif,
    get_library,
    reverse_complement,
)


# ---------------------------------------------------------------------------
# Background sequence
# ---------------------------------------------------------------------------

def _target_freqs(gc_content: float) -> np.ndarray:
    gc = float(np.clip(gc_content, 0.05, 0.95))
    return np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2])  # A C G T


def count_cpg(seq: str) -> int:
    return seq.count("CG")


def cpg_observed_expected(seq: str) -> float:
    """Observed/expected CpG ratio: (n_CpG * L) / (n_C * n_G)."""
    n_c, n_g = seq.count("C"), seq.count("G")
    if n_c == 0 or n_g == 0:
        return 0.0
    return count_cpg(seq) * len(seq) / (n_c * n_g)


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def _tune_cpg(chars: list[str], target_cpg: int, rng: np.random.Generator,
              max_iter: int = 40000) -> None:
    """Move G's around, in place, until the CpG count is near ``target_cpg``.

    Only *swaps* bases between two positions, so the mononucleotide composition
    (and therefore GC content) is preserved exactly.
    """
    n = len(chars)
    if n < 2:
        return

    def cpg_at(i: int) -> int:
        return 1 if (i >= 0 and i + 1 < n and chars[i] == "C" and chars[i + 1] == "G") else 0

    current = sum(cpg_at(i) for i in range(n - 1))

    for _ in range(max_iter):
        if current == target_cpg:
            break
        need_more = current < target_cpg

        if need_more:
            # find a C whose right neighbour is not G, and a G to swap in
            i = int(rng.integers(0, n - 1))
            if chars[i] != "C" or chars[i + 1] == "G":
                continue
            j = int(rng.integers(0, n))
            if chars[j] != "G" or j == i or j == i + 1:
                continue
        else:
            # find an existing CpG and swap its G out for a non-G
            i = int(rng.integers(0, n - 1))
            if chars[i] != "C" or chars[i + 1] != "G":
                continue
            j = int(rng.integers(0, n))
            if chars[j] == "G" or j == i or j == i + 1:
                continue

        a, b = i + 1, j
        touched = {a - 1, a, b - 1, b}
        before = sum(cpg_at(t) for t in touched)
        chars[a], chars[b] = chars[b], chars[a]
        after = sum(cpg_at(t) for t in touched)
        delta = after - before

        improves = (delta > 0) if need_more else (delta < 0)
        if improves:
            current += delta
        else:
            chars[a], chars[b] = chars[b], chars[a]


def random_background(length: int = 2001, tss_index: int | None = None,
                      gc: float = 0.45, cpg_oe: float = 0.25,
                      seed: int | None = None) -> dict:
    """Generate a random background with a target GC content and CpG o/e ratio.

    ``cpg_oe`` near 0.2-0.25 is typical of bulk human genomic sequence;
    values near 0.7-1.0 give CpG-island-like background.
    """
    length = int(np.clip(length, 100, 20000))
    if tss_index is None:
        tss_index = length // 2
    tss_index = int(np.clip(tss_index, 0, length - 1))

    rng = np.random.default_rng(seed)
    freqs = _target_freqs(gc)
    draw = rng.choice(4, size=length, p=freqs)
    chars = [ALPHABET[i] for i in draw]

    n_c = chars.count("C")
    n_g = chars.count("G")
    target_cpg = int(round(float(np.clip(cpg_oe, 0.0, 3.0)) * n_c * n_g / length))
    _tune_cpg(chars, target_cpg, rng)

    seq = "".join(chars)
    return {
        "sequence": seq,
        "tss_index": tss_index,
        "length": length,
        "requested": {"gc": gc, "cpg_oe": cpg_oe, "seed": seed},
        "achieved": {
            "gc": round(gc_content(seq), 4),
            "cpg_oe": round(cpg_observed_expected(seq), 4),
            "n_cpg": count_cpg(seq),
        },
    }


def cpg_rich_segment(length: int, rng: np.random.Generator,
                     gc: float = 0.70, cpg_oe: float = 0.85) -> str:
    """A CpG-island-like stretch, for placing as a segment element."""
    freqs = _target_freqs(gc)
    chars = [ALPHABET[i] for i in rng.choice(4, size=max(1, length), p=freqs)]
    n_c, n_g = chars.count("C"), chars.count("G")
    target = int(round(cpg_oe * n_c * n_g / max(1, len(chars))))
    _tune_cpg(chars, target, rng)
    return "".join(chars)


# ---------------------------------------------------------------------------
# Elements and construct assembly
# ---------------------------------------------------------------------------

@dataclass
class Placement:
    """One element placed on a construct.

    ``element_id`` is either a motif id from the library or ``"cpg_segment"``.
    ``position`` is the relative position (to the TSS) of the element's 5' end
    on the construct's plus strand.
    """

    element_id: str
    position: int
    strand: str = "+"
    instance: str = "consensus"   # "consensus" | "sampled"
    custom_sequence: str = ""     # only used by element_id == "custom"
    seg_length: int = 60          # only used by cpg_segment
    seg_gc: float = 0.70
    seg_cpg_oe: float = 0.85
    uid: str | None = None
    locked: bool = False

    def to_json(self) -> dict:
        return asdict(self)


def element_sequence(p: Placement, rng: np.random.Generator) -> str:
    """Concrete bases this placement writes onto the construct (plus strand)."""
    if p.element_id == CUSTOM_ID:
        seq = "".join(c for c in p.custom_sequence.upper() if c in "ACGTN")
        return seq if p.strand == "+" else reverse_complement(seq)
    if p.element_id == CPG_SEGMENT_ID:
        return cpg_rich_segment(p.seg_length, rng, p.seg_gc, p.seg_cpg_oe)

    lib = get_library()
    motif: Motif = lib[p.element_id]
    if p.instance == "sampled":
        return motif.sample(rng, strand=p.strand)
    return motif.consensus_instance(strand=p.strand)


def element_width(p: Placement) -> int:
    if p.element_id == CUSTOM_ID:
        return max(1, len([c for c in p.custom_sequence.upper() if c in "ACGTN"]))
    if p.element_id == CPG_SEGMENT_ID:
        return max(1, int(p.seg_length))
    return get_library()[p.element_id].width


def build_construct(background: str, tss_index: int,
                    placements: list[Placement],
                    seed: int | None = 0) -> dict:
    """Overwrite ``background`` with each placement, preserving total length.

    Later placements overwrite earlier ones where they overlap; the returned
    ``overlaps`` list flags every such collision so the UI can warn.
    """
    chars = list(background)
    L = len(chars)
    rng = np.random.default_rng(seed)

    owner = [None] * L
    overlaps: list[dict] = []
    realized: list[dict] = []

    for p in placements:
        sub = element_sequence(p, rng)
        start = tss_index + int(p.position)
        end = start + len(sub)

        clipped_start = max(0, start)
        clipped_end = min(L, end)
        if clipped_start >= clipped_end:
            realized.append({**p.to_json(), "written": "", "out_of_range": True,
                             "start_index": start, "end_index": end})
            continue

        hit = {owner[i] for i in range(clipped_start, clipped_end) if owner[i] is not None}
        if hit:
            overlaps.append({"element": p.uid or p.element_id,
                             "overlaps_with": sorted(str(h) for h in hit)})

        for i in range(clipped_start, clipped_end):
            chars[i] = sub[i - start]
            owner[i] = p.uid or p.element_id

        realized.append({
            **p.to_json(),
            "written": sub[clipped_start - start: clipped_end - start],
            "out_of_range": clipped_start != start or clipped_end != end,
            "start_index": start,
            "end_index": end,
            "start_rel": int(p.position),
            "end_rel": int(p.position) + len(sub) - 1,
        })

    seq = "".join(chars)
    return {
        "sequence": seq,
        "tss_index": tss_index,
        "length": L,
        "placements": realized,
        "overlaps": overlaps,
        "stats": {
            "gc": round(gc_content(seq), 4),
            "cpg_oe": round(cpg_observed_expected(seq), 4),
            "n_cpg": count_cpg(seq),
        },
    }


def sequence_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode()).hexdigest()[:12]


def to_placements(raw: list[dict]) -> list[Placement]:
    out = []
    for i, r in enumerate(raw or []):
        out.append(Placement(
            element_id=r["element_id"],
            position=int(r.get("position", 0)),
            strand="-" if str(r.get("strand", "+")) == "-" else "+",
            instance=r.get("instance", "consensus"),
            custom_sequence=str(r.get("custom_sequence", "")),
            seg_length=int(r.get("seg_length", 60)),
            seg_gc=float(r.get("seg_gc", 0.70)),
            seg_cpg_oe=float(r.get("seg_cpg_oe", 0.85)),
            uid=r.get("uid") or f"el{i}",
            locked=bool(r.get("locked", False)),
        ))
    return out


def to_fasta(seq: str, name: str = "construct", width: int = 60) -> str:
    lines = [f">{name}"]
    for i in range(0, len(seq), width):
        lines.append(seq[i:i + width])
    return "\n".join(lines) + "\n"


__all__ = [
    "Placement", "build_construct", "random_background", "cpg_rich_segment",
    "gc_content", "cpg_observed_expected", "count_cpg", "to_placements",
    "element_width", "sequence_hash", "to_fasta", "reverse_complement",
]
