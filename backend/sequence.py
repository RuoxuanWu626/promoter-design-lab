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


def scrub_motifs(seq: str, motif_ids: list[str], rng: np.random.Generator,
                 occupancy: float = 0.5, max_rounds: int = 25) -> tuple[str, int]:
    """Remove strong matches to ``motif_ids`` from a sequence.

    Random sequence contains real transcription factor sites by chance -- a
    2 kb background typically carries several. They are genuine as far as the
    model is concerned and they drive cell-type specificity, which makes a
    designed site impossible to see against them.

    Scrubbing gives a blank canvas to design on. It is a deliberate departure
    from random sequence, not a neutral default, so the caller records that it
    happened.

    Two things make this harder than it looks, and both are handled here:

    * A dinucleotide-preserving shuffle of an 8-10 bp patch has very few valid
      permutations and frequently returns the patch unchanged, so shuffling the
      motif footprint alone does not remove it. The patch is widened before
      shuffling, which gives the shuffle room to move.
    * A shuffle can create a *new* match, so each attempt is verified and
      retried, and the loop only stops when a full re-scan comes back clean.

    Composition is preserved by the shuffle; only if repeated shuffles fail
    does it fall back to redrawing the patch, which perturbs composition
    locally.
    """
    from adapters.mock_common import motif_hits_by_motif, motif_occupancy
    lib = get_library()
    removed = 0
    pad = 6

    for _ in range(max_rounds):
        hits = []
        for mid, d in motif_hits_by_motif(seq, motif_ids).items():
            w = lib[mid].width
            for i, occ in zip(d["index"].tolist(), d["occupancy"].tolist()):
                if occ >= occupancy:
                    hits.append((mid, int(i), w))
        if not hits:
            break

        for mid, start_i, w in hits:
            lo = max(0, start_i - pad)
            hi = min(len(seq), start_i + w + pad)
            motif = lib[mid]
            for attempt in range(8):
                mode = "dinuc_shuffle" if attempt < 5 else "mono_shuffle"
                cand = perturb_patch(seq, lo, hi - lo, rng, mode)
                if cand == seq:
                    continue
                # Did it actually clear, without creating a new match nearby?
                region = cand[max(0, lo - motif.width):min(len(cand), hi + motif.width)]
                if motif_occupancy(region, motif).max(initial=0.0) < occupancy:
                    seq = cand
                    removed += 1
                    break
            else:
                # Repeated shuffles could not clear it; redraw the patch.
                sub = "".join(rng.choice(list(ALPHABET), hi - lo))
                seq = seq[:lo] + sub + seq[hi:]
                removed += 1

    return seq, removed


def random_background(length: int = 2001, tss_index: int | None = None,
                      gc: float = 0.45, cpg_oe: float = 0.25,
                      seed: int | None = None,
                      scrub_celltype_elements: bool = True) -> dict:
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

    scrubbed = 0
    if scrub_celltype_elements:
        from motifs import celltype_kind_ids
        seq, scrubbed = scrub_motifs(seq, celltype_kind_ids(), rng)

    return {
        "sequence": seq,
        "tss_index": tss_index,
        "length": length,
        "scrubbed_celltype_sites": scrubbed,
        "scrub_note": (
            f"{scrubbed} chance lineage TF site(s) were shuffled out of the "
            f"background so designed sites are visible against it. Base and "
            f"CpG composition are unchanged. Turn this off to design against "
            f"genuinely random sequence." if scrubbed else
            ("background left as drawn; no strong chance lineage site was "
             "present" if scrub_celltype_elements else
             "scrubbing off: the background may contain chance lineage sites")),
        "requested": {"gc": gc, "cpg_oe": cpg_oe, "seed": seed,
                      "scrub_celltype_elements": scrub_celltype_elements},
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


# ---------------------------------------------------------------------------
# Perturbations, for attribution
# ---------------------------------------------------------------------------

def dinuc_shuffle(seq: str, rng: np.random.Generator) -> str:
    """Shuffle a sequence while preserving its dinucleotide composition.

    Altschul-Erikson: the dinucleotides of a sequence are the edges of a
    multigraph on {A,C,G,T}, so any Eulerian path with the same start and end
    vertex has exactly the same dinucleotide counts.

    Preserving dinucleotides (not just base composition) matters here because
    Puffin has an explicit trinucleotide branch and a CpG-sensitive response.
    A plain mononucleotide shuffle changes local CpG content, so attribution
    computed against it would confound "this motif matters" with "the CpG
    content here matters".
    """
    if len(seq) < 3:
        return seq

    edges: dict[str, list[str]] = {}
    for a, b in zip(seq[:-1], seq[1:]):
        edges.setdefault(a, []).append(b)
    first, last = seq[0], seq[-1]
    verts = list(edges)

    # Pick, for every vertex but the last, one outgoing edge to traverse last.
    # Those edges must form a tree rooted at `last`, or the walk dead-ends.
    chosen: dict[str, str] = {}
    for _ in range(200):
        chosen = {v: edges[v][int(rng.integers(len(edges[v])))]
                  for v in verts if v != last}
        ok = True
        for v in verts:
            if v == last:
                continue
            seen, u = set(), v
            while u != last:
                if u in seen or u not in chosen:
                    ok = False
                    break
                seen.add(u)
                u = chosen[u]
            if not ok:
                break
        if ok:
            break
    else:
        return seq  # give up rather than return something with wrong statistics

    order: dict[str, list[str]] = {}
    for v in verts:
        rest = list(edges[v])
        if v != last:
            rest.remove(chosen[v])
            rng.shuffle(rest)
            rest.append(chosen[v])
        else:
            rng.shuffle(rest)
        order[v] = rest

    out = [first]
    used: dict[str, int] = {v: 0 for v in verts}
    u = first
    for _ in range(len(seq) - 1):
        nxt = order[u][used[u]]
        used[u] += 1
        out.append(nxt)
        u = nxt
    return "".join(out)


def perturb_patch(sequence: str, start: int, width: int,
                  rng: np.random.Generator, mode: str = "dinuc_shuffle") -> str:
    """Return ``sequence`` with the bases in [start, start+width) disrupted.

    ``dinuc_shuffle`` keeps composition and destroys motifs; ``mono_shuffle``
    keeps only base composition; ``neutral`` writes N, which removes the bases
    from the model's view entirely rather than replacing them.
    """
    lo = max(0, int(start))
    hi = min(len(sequence), int(start + width))
    if hi <= lo:
        return sequence
    patch = sequence[lo:hi]
    if mode == "neutral":
        new = "N" * len(patch)
    elif mode == "mono_shuffle":
        arr = list(patch)
        rng.shuffle(arr)
        new = "".join(arr)
    else:
        new = dinuc_shuffle(patch, rng)
    return sequence[:lo] + new + sequence[hi:]
