"""User-defined motifs — a third library, stored on the server.

Puffin's ten filters and the lineage TF sites are fixed. This is where anything
else goes: a factor you care about that is not in either list, a variant of a
known site, a PWM from JASPAR, or a candidate element you want to test.

Custom motifs are saved to ``data/custom_motifs.json`` so they persist across
restarts and every player on the same server sees them — which is the point,
since challenges are meant to be shared.

What a custom motif does and does not do
----------------------------------------
Placing one **writes real bases into the construct**, so every model responds
to it, including the real Puffin: Puffin has no filter for your motif, but it
does see the sequence, and if your motif happens to contain a GC-box or a TATA
it will react accordingly. That is a feature — it is how you find out whether
an element does anything through the core promoter machinery.

What a custom motif cannot do is teach a *trained* model a new preference. The
``activates`` / ``represses`` fields only affect the mock cell-type adapter,
which is explicit about being a caricature. Against a real cell-type model, a
custom motif's effect is whatever that model already thinks of the sequence.
"""

from __future__ import annotations

import json
import os
import re
import threading

import numpy as np

from motifs import IUPAC, Motif, consensus_to_pwm

DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
STORE_PATH = os.path.join(DATA_DIR, "custom_motifs.json")

_LOCK = threading.RLock()

MAX_MOTIFS = 60
MIN_WIDTH = 4
MAX_WIDTH = 60

ID_RE = re.compile(r"^[a-z0-9][a-z0-9_]{1,23}$")
HEX_RE = re.compile(r"^#[0-9a-fA-F]{6}$")

DEFAULT_COLORS = [
    "#e76f51", "#2a9d8f", "#e9c46a", "#f4a261", "#8ab17d",
    "#7678ed", "#f35b04", "#3d5a80", "#ee6c4d", "#98c1d9",
]


class InvalidMotif(ValueError):
    """Raised with a message meant to be shown directly to the user."""


def _load_raw() -> list[dict]:
    if not os.path.exists(STORE_PATH):
        return []
    try:
        with open(STORE_PATH) as fh:
            data = json.load(fh)
        return data.get("motifs", []) if isinstance(data, dict) else list(data)
    except Exception:
        os.replace(STORE_PATH, STORE_PATH + ".corrupt")
        return []


def _save_raw(motifs: list[dict]) -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    tmp = STORE_PATH + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"version": 1, "motifs": motifs}, fh, indent=1)
    os.replace(tmp, STORE_PATH)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_pwm(raw) -> np.ndarray:
    """Accept a PWM as 4 rows (A,C,G,T) or W rows of 4, counts or probabilities.

    JASPAR-style pasted matrices are 4 rows of counts; MEME-style are W rows of
    probabilities. Both are normalised to columns summing to 1.
    """
    if isinstance(raw, str):
        rows = []
        for line in raw.strip().splitlines():
            line = re.sub(r"^\s*[ACGTacgt]?\s*[:|\[]", "", line)
            line = line.replace("]", " ")
            nums = [float(x) for x in re.findall(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", line)]
            if nums:
                rows.append(nums)
        raw = rows

    arr = np.asarray(raw, dtype=np.float64)
    if arr.ndim != 2:
        raise InvalidMotif("A PWM needs to be a 2-D matrix.")
    if arr.shape[0] == 4 and arr.shape[1] != 4:
        mat = arr
    elif arr.shape[1] == 4:
        mat = arr.T
    elif arr.shape == (4, 4):
        mat = arr           # ambiguous; treat as 4 rows A,C,G,T
    else:
        raise InvalidMotif(
            f"Expected 4 rows (A,C,G,T) or 4 columns; got {arr.shape[0]}x{arr.shape[1]}.")

    if mat.shape[1] < MIN_WIDTH or mat.shape[1] > MAX_WIDTH:
        raise InvalidMotif(f"Motif width must be {MIN_WIDTH}-{MAX_WIDTH} bp; "
                           f"this one is {mat.shape[1]}.")
    if not np.isfinite(mat).all() or (mat < 0).any():
        raise InvalidMotif("PWM entries must be finite and non-negative.")

    totals = mat.sum(axis=0, keepdims=True)
    if (totals <= 0).any():
        raise InvalidMotif("Every PWM column needs at least one non-zero entry.")
    return mat / totals


def parse_consensus(text: str) -> tuple[str, np.ndarray]:
    seq = re.sub(r"\s+", "", str(text or "")).upper()
    if not seq:
        raise InvalidMotif("Give a consensus sequence or a PWM.")
    bad = sorted(set(seq) - set(IUPAC))
    if bad:
        raise InvalidMotif(f"Not IUPAC codes: {', '.join(bad)}. "
                           f"Allowed: {' '.join(sorted(IUPAC))}.")
    if len(seq) < MIN_WIDTH or len(seq) > MAX_WIDTH:
        raise InvalidMotif(f"Motif width must be {MIN_WIDTH}-{MAX_WIDTH} bp; "
                           f"this one is {len(seq)}.")
    if all(c == "N" for c in seq):
        raise InvalidMotif("An all-N motif matches everything; add some "
                           "determined positions.")
    return seq, consensus_to_pwm(seq)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate(spec: dict, existing_ids: set[str], reserved: set[str]) -> dict:
    name = str(spec.get("name") or "").strip()
    if not name:
        raise InvalidMotif("Give the motif a name.")
    if len(name) > 40:
        raise InvalidMotif("Name must be 40 characters or fewer.")

    mid = str(spec.get("id") or "").strip().lower()
    if not mid:
        slug = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")[:20] or "motif"
        # Prefixed so a user motif can never be mistaken for a built-in, but
        # not doubled up when the name already begins with "my".
        mid = slug if slug.startswith("my_") else "my_" + slug
        mid = mid[:24].rstrip("_")
    if not ID_RE.match(mid):
        raise InvalidMotif("Id must be 2-24 characters: lowercase letters, "
                           "digits and underscores, starting with a letter or digit.")
    if mid in reserved:
        raise InvalidMotif(f"'{mid}' is a built-in motif id; choose another.")
    if mid in existing_ids and not spec.get("replace"):
        raise InvalidMotif(f"A custom motif called '{mid}' already exists. "
                           f"Rename it, or tick replace.")

    if spec.get("pwm"):
        pwm = parse_pwm(spec["pwm"])
        consensus = "".join("ACGT"[i] for i in pwm.argmax(axis=0))
    else:
        consensus, pwm = parse_consensus(spec.get("consensus"))

    kind = spec.get("kind", "core_promoter")
    if kind not in ("core_promoter", "celltype_element"):
        raise InvalidMotif("kind must be 'core_promoter' or 'celltype_element'.")

    color = str(spec.get("color") or "").strip()
    if not HEX_RE.match(color):
        color = DEFAULT_COLORS[abs(hash(mid)) % len(DEFAULT_COLORS)]

    def _cts(key: str) -> list[str]:
        vals = spec.get(key) or []
        if isinstance(vals, str):
            vals = [v.strip() for v in vals.split(",")]
        return [str(v).strip() for v in vals if str(v).strip()][:24]

    try:
        offset = int(spec.get("typical_offset", -60))
    except (TypeError, ValueError):
        raise InvalidMotif("Preferred position must be a whole number of bp.")
    offset = max(-5000, min(5000, offset))

    try:
        amplitude = float(spec.get("amplitude", 0.8))
    except (TypeError, ValueError):
        raise InvalidMotif("Amplitude must be a number.")
    amplitude = float(np.clip(amplitude, 0.0, 3.0))

    return {
        "id": mid,
        "name": name,
        "consensus": consensus,
        "pwm": np.round(pwm, 6).tolist(),
        "kind": kind,
        "color": color,
        "typical_offset": offset,
        "amplitude": amplitude,
        "activates": _cts("activates"),
        "represses": _cts("represses"),
        "notes": str(spec.get("notes") or "").strip()[:300],
        "author": str(spec.get("author") or "").strip()[:48],
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def list_motifs() -> list[dict]:
    with _LOCK:
        return _load_raw()


def add(spec: dict) -> dict:
    from motifs import get_library, invalidate_cache, library_order, celltype_element_order

    with _LOCK:
        current = _load_raw()
        if len(current) >= MAX_MOTIFS and not spec.get("replace"):
            raise InvalidMotif(f"At most {MAX_MOTIFS} custom motifs; delete some first.")
        reserved = set(library_order()) | set(celltype_element_order()) | {
            "cpg_segment", "custom"}
        existing = {m["id"] for m in current}
        entry = validate(spec, existing, reserved)
        current = [m for m in current if m["id"] != entry["id"]]
        current.append(entry)
        _save_raw(current)

    invalidate_cache()
    return entry


def remove(mid: str) -> bool:
    from motifs import invalidate_cache

    with _LOCK:
        current = _load_raw()
        kept = [m for m in current if m["id"] != mid]
        if len(kept) == len(current):
            return False
        _save_raw(kept)

    invalidate_cache()
    return True


def build_motifs() -> dict[str, Motif]:
    """Return {id: Motif} for everything the user has defined."""
    out: dict[str, Motif] = {}
    for spec in list_motifs():
        try:
            pwm = np.asarray(spec["pwm"], dtype=np.float64)
            pwm = pwm / pwm.sum(axis=0, keepdims=True)
        except Exception:
            continue
        kind = spec.get("kind", "core_promoter")
        m = Motif(
            id=spec["id"],
            name=spec["name"],
            aliases=[],
            consensus=spec["consensus"],
            typical_offset=int(spec.get("typical_offset", -60)),
            strand_specific=False,
            color=spec.get("color", "#8899aa"),
            notes=(spec.get("notes") or
                   ("User-defined motif. Placing it writes real bases, so every "
                    "model reacts to the sequence; the cell-type assignments "
                    "below only steer the mock cell-type adapter.")),
            mock_amplitude=float(spec.get("amplitude", 0.8)),
            mock_peak_shift=-int(spec.get("typical_offset", -60)),
            mock_width=14.0 if kind == "core_promoter" else 40.0,
            mock_broad_fraction=0.3 if kind == "core_promoter" else 0.75,
            mock_strand_asymmetry=0.9,
        )
        m.pwm = pwm
        m.kind = kind
        m.activates = list(spec.get("activates") or [])
        m.represses = list(spec.get("represses") or [])
        m.source = "user"
        out[m.id] = m
    return out


def activation_map() -> dict[str, dict[str, float]]:
    """{motif_id: {cell_type: weight}} contributed by user motifs."""
    out: dict[str, dict[str, float]] = {}
    for spec in list_motifs():
        if spec.get("kind") != "celltype_element":
            continue
        w: dict[str, float] = {}
        for ct in spec.get("activates") or []:
            w[ct] = 1.0
        for ct in spec.get("represses") or []:
            w[ct] = -0.85
        if w:
            out[spec["id"]] = w
    return out


def order() -> list[str]:
    return [m["id"] for m in list_motifs()]
