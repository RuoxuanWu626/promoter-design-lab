"""Players, challenges, submissions and the leaderboard.

Storage is a single JSON file guarded by a lock, plus an atomic replace on
write.  That is enough for a lab-sized group on one node and keeps the whole
thing dependency-free; swap in SQLite if this ever needs to outlive a session
or serve real concurrency.

A *challenge* stores the background *parameters and seed*, not the sequence.
Regenerating from the seed is deterministic, so every player provably starts
from the same background and a challenge stays small enough to share as a URL.
"""

from __future__ import annotations

import json
import os
import threading
import time
import uuid

DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
STORE_PATH = os.path.join(DATA_DIR, "store.json")

_LOCK = threading.RLock()

_EMPTY = {"players": {}, "challenges": {}, "submissions": {}, "version": 1}


def _load() -> dict:
    if not os.path.exists(STORE_PATH):
        return json.loads(json.dumps(_EMPTY))
    try:
        with open(STORE_PATH) as fh:
            data = json.load(fh)
        for k in _EMPTY:
            data.setdefault(k, _EMPTY[k] if not isinstance(_EMPTY[k], dict) else {})
        return data
    except Exception:
        # Never lose a corrupt store silently: move it aside and start clean.
        os.replace(STORE_PATH, STORE_PATH + f".corrupt.{int(time.time())}")
        return json.loads(json.dumps(_EMPTY))


def _save(data: dict) -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    tmp = STORE_PATH + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=1)
    os.replace(tmp, STORE_PATH)


def _new_id(prefix: str) -> str:
    return f"{prefix}_{uuid.uuid4().hex[:10]}"


def _now() -> float:
    return time.time()


# ---------------------------------------------------------------------------
# Players
# ---------------------------------------------------------------------------

def ensure_player(name: str, player_id: str | None = None) -> dict:
    name = (name or "anonymous").strip()[:48] or "anonymous"
    with _LOCK:
        data = _load()
        if player_id and player_id in data["players"]:
            p = data["players"][player_id]
            if p["name"] != name:
                p["name"] = name
                _save(data)
            return p
        for p in data["players"].values():
            if p["name"].lower() == name.lower():
                return p
        pid = _new_id("p")
        p = {"id": pid, "name": name, "created": _now()}
        data["players"][pid] = p
        _save(data)
        return p


def list_players() -> list[dict]:
    with _LOCK:
        return sorted(_load()["players"].values(), key=lambda p: p["created"])


# ---------------------------------------------------------------------------
# Challenges (shared starting sequences)
# ---------------------------------------------------------------------------

def create_challenge(name: str, background_params: dict, mode: str = "design",
                     target_cell_type: str | None = None,
                     weights: dict | None = None,
                     models: dict | None = None,
                     created_by: str | None = None,
                     notes: str = "") -> dict:
    with _LOCK:
        data = _load()
        cid = _new_id("c")
        ch = {
            "id": cid,
            "name": (name or "Untitled challenge").strip()[:80],
            "mode": mode,
            "background_params": background_params,
            "target_cell_type": target_cell_type,
            "weights": weights or {},
            "models": models or {},
            "created_by": created_by,
            "notes": notes[:500],
            "created": _now(),
        }
        data["challenges"][cid] = ch
        _save(data)
        return ch


def get_challenge(cid: str) -> dict | None:
    with _LOCK:
        return _load()["challenges"].get(cid)


def list_challenges() -> list[dict]:
    with _LOCK:
        data = _load()
        out = []
        for ch in data["challenges"].values():
            n = sum(1 for s in data["submissions"].values()
                    if s["challenge_id"] == ch["id"])
            out.append({**ch, "n_submissions": n})
        return sorted(out, key=lambda c: -c["created"])


def delete_challenge(cid: str) -> bool:
    with _LOCK:
        data = _load()
        if cid not in data["challenges"]:
            return False
        del data["challenges"][cid]
        for sid in [s for s, v in data["submissions"].items()
                    if v["challenge_id"] == cid]:
            del data["submissions"][sid]
        _save(data)
        return True


# ---------------------------------------------------------------------------
# Submissions and leaderboard
# ---------------------------------------------------------------------------

def submit(challenge_id: str, player_id: str, placements: list[dict],
           metrics: dict, score: dict, label: str = "",
           extra: dict | None = None) -> dict:
    with _LOCK:
        data = _load()
        sid = _new_id("s")
        sub = {
            "id": sid,
            "challenge_id": challenge_id,
            "player_id": player_id,
            "player_name": data["players"].get(player_id, {}).get("name", "?"),
            "label": (label or "").strip()[:60],
            "placements": placements,
            "metrics": metrics,
            "score": score,
            "extra": extra or {},
            "created": _now(),
        }
        data["submissions"][sid] = sub
        _save(data)
        return sub


def leaderboard(challenge_id: str, limit: int = 50,
                best_per_player: bool = True) -> list[dict]:
    with _LOCK:
        data = _load()
        subs = [s for s in data["submissions"].values()
                if s["challenge_id"] == challenge_id]

    subs.sort(key=lambda s: -float(s["score"].get("score", 0)))
    if best_per_player:
        seen = set()
        kept = []
        for s in subs:
            if s["player_id"] in seen:
                continue
            seen.add(s["player_id"])
            kept.append(s)
        subs = kept

    out = []
    for rank, s in enumerate(subs[:limit], start=1):
        out.append({
            "rank": rank,
            "submission_id": s["id"],
            "player_name": s["player_name"],
            "label": s["label"],
            "score": s["score"].get("score"),
            "score_normalised": s["score"].get("score_normalised"),
            "metrics": s["metrics"],
            "contributions": s["score"].get("contributions", {}),
            "n_elements": len(s["placements"]),
            "created": s["created"],
        })
    return out


def get_submission(sid: str) -> dict | None:
    with _LOCK:
        return _load()["submissions"].get(sid)


def all_submissions(challenge_id: str | None = None) -> list[dict]:
    with _LOCK:
        subs = list(_load()["submissions"].values())
    if challenge_id:
        subs = [s for s in subs if s["challenge_id"] == challenge_id]
    return sorted(subs, key=lambda s: s["created"])


def stats() -> dict:
    with _LOCK:
        d = _load()
        return {"players": len(d["players"]),
                "challenges": len(d["challenges"]),
                "submissions": len(d["submissions"]),
                "path": STORE_PATH}
