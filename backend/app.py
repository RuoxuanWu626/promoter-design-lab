"""HTTP server for the promoter design game.

Standard library only (plus numpy), so it runs on a login or compute node with
no installs.  Long experiments run as background jobs that the client polls,
which keeps the UI responsive during a 10x10 pair grid.

    python3 app.py --port 8765 --host 127.0.0.1

Then forward the port from your laptop:

    ssh -N -L 8765:<node>:8765 Randi
"""

from __future__ import annotations

import argparse
import json
import mimetypes
import os
import socket
import sys
import threading
import time
import traceback
import urllib.parse
import uuid
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

import adapters
import experiments
import scoring
import store
from motifs import CPG_SEGMENT_ID, CUSTOM_ID, get_library, library_order
from sequence import (
    build_construct,
    random_background,
    sequence_hash,
    to_fasta,
    to_placements,
)

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FRONTEND = os.path.join(ROOT, "frontend")
RUNS_DIR = os.path.join(ROOT, "data", "runs")

MAX_BODY = 32 * 1024 * 1024


# ---------------------------------------------------------------------------
# Background jobs
# ---------------------------------------------------------------------------

_JOBS: dict[str, dict] = {}
_JOBS_LOCK = threading.Lock()
_JOB_TTL = 3600.0


def _reap_jobs() -> None:
    now = time.time()
    with _JOBS_LOCK:
        for jid in [j for j, v in _JOBS.items()
                    if v.get("finished") and now - v["finished"] > _JOB_TTL]:
            del _JOBS[jid]


def start_job(fn, *args, **kwargs) -> str:
    _reap_jobs()
    jid = "j_" + uuid.uuid4().hex[:10]
    progress: dict = {"done": 0, "total": 0, "stage": "queued"}
    with _JOBS_LOCK:
        _JOBS[jid] = {"id": jid, "status": "running", "progress": progress,
                      "started": time.time(), "finished": None,
                      "result": None, "error": None}

    def run():
        try:
            result = fn(*args, progress=progress, **kwargs)
            with _JOBS_LOCK:
                _JOBS[jid].update(status="done", result=result,
                                  finished=time.time())
        except Exception as exc:
            with _JOBS_LOCK:
                _JOBS[jid].update(status="error",
                                  error=f"{type(exc).__name__}: {exc}",
                                  traceback=traceback.format_exc(),
                                  finished=time.time())

    threading.Thread(target=run, daemon=True).start()
    return jid


def job_status(jid: str) -> dict | None:
    with _JOBS_LOCK:
        j = _JOBS.get(jid)
        if not j:
            return None
        out = {k: v for k, v in j.items() if k != "result"}
        if j["status"] == "done":
            out["result"] = j["result"]
        return out


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

class JSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            v = float(o)
            return v if np.isfinite(v) else None
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        return super().default(o)


def _win(body: dict, key: str = "window", default=(-500, 500)) -> tuple[int, int]:
    w = body.get(key) or default
    return (int(w[0]), int(w[1]))


def _background_from(body: dict) -> dict:
    """Accept either an explicit sequence or background parameters."""
    if body.get("sequence"):
        seq = "".join(c for c in str(body["sequence"]).upper() if c in "ACGTN")
        tss = int(body.get("tss_index", len(seq) // 2))
        return {"sequence": seq, "tss_index": tss, "length": len(seq),
                "requested": {"seed": body.get("seed")},
                "achieved": {}}
    p = body.get("background_params") or body
    return random_background(
        length=int(p.get("length", 1001)),
        tss_index=p.get("tss_index"),
        gc=float(p.get("gc", 0.45)),
        cpg_oe=float(p.get("cpg_oe", 0.25)),
        seed=p.get("seed"))


def _backgrounds_from(body: dict) -> list[dict]:
    n = int(body.get("n_backgrounds", 1))
    if body.get("sequence") or n <= 1:
        return [_background_from(body)]
    p = body.get("background_params") or body
    return experiments.matched_backgrounds(
        n=n, length=int(p.get("length", 1001)), gc=float(p.get("gc", 0.45)),
        cpg_oe=float(p.get("cpg_oe", 0.25)),
        base_seed=int(p.get("seed") or 1000))


def _evaluate(body: dict) -> dict:
    """Build a construct, predict with every requested model, and score it."""
    bg = _background_from(body)
    placements = to_placements(body.get("placements") or [])
    con = build_construct(bg["sequence"], bg["tss_index"], placements,
                          seed=int(body.get("instance_seed", 0)))
    window = _win(body)

    profile_models = body.get("profile_models") or adapters.DEFAULT_PROFILE_MODELS
    ct_model = body.get("celltype_model") or adapters.DEFAULT_CELLTYPE_MODEL
    cell_types = body.get("cell_types")
    target = body.get("target_cell_type")
    act_win = tuple(body.get("activity_window") or (-200, 200))
    act_method = body.get("activity_method", "mean")

    out_profiles = {}
    metrics: dict = {}
    for m in profile_models:
        try:
            pred = adapters.predict_profile(m, con["sequence"], con["tss_index"], window)
        except Exception as exc:
            out_profiles[m] = {"error": f"{type(exc).__name__}: {exc}"}
            continue
        out_profiles[m] = pred.to_json()
        pm = scoring.profile_metrics(pred.positions, pred.tracks["plus"],
                                     pred.tracks.get("minus"),
                                     score_window=window)
        out_profiles[m]["metrics"] = pm
        if m == profile_models[0]:
            metrics.update(pm)

    celltype = None
    if body.get("include_celltype", True):
        try:
            ctp = adapters.predict_celltype(ct_model, con["sequence"],
                                            con["tss_index"], cell_types, window)
            acts = np.array([
                scoring.activity_from_profile(ctp.profiles[i], ctp.positions,
                                              act_win, act_method)
                for i in range(len(ctp.cell_types))])
            cm = scoring.celltype_metrics(list(ctp.cell_types), acts, target)
            celltype = {**ctp.to_json(), "metrics": cm,
                        "activity_window": list(act_win),
                        "activity_method": act_method}
            metrics.update({k: v for k, v in cm.items()
                            if isinstance(v, (int, float))})
            # Carried separately because the filter above keeps only numbers,
            # and the leaderboard wants to show which cell type actually won.
            metrics["strongest_cell_type"] = cm.get("strongest_cell_type")
            if cm.get("target"):
                metrics["target_cell_type"] = cm["target"]
        except Exception as exc:
            celltype = {"error": f"{type(exc).__name__}: {exc}",
                        "model": ct_model}

    weights = body.get("weights") or scoring.DEFAULT_WEIGHTS_DESIGN
    score = scoring.weighted_score(metrics, weights)

    return {
        "construct": {k: v for k, v in con.items() if k != "sequence"},
        "sequence": con["sequence"],
        "sequence_hash": sequence_hash(con["sequence"]),
        "background": {"hash": sequence_hash(bg["sequence"]),
                       "stats": bg.get("achieved", {}),
                       "requested": bg.get("requested", {}),
                       "tss_index": bg["tss_index"]},
        "window": list(window),
        "profiles": out_profiles,
        "celltype": celltype,
        "metrics": metrics,
        "score": score,
        "weights": weights,
    }


# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

def route(path: str, body: dict, query: dict) -> dict:
    lib = get_library()

    # --- static metadata ---------------------------------------------------
    if path == "/api/health":
        return {"ok": True, "time": time.time(), "cache": adapters.cache_stats(),
                "store": store.stats(), "pid": os.getpid()}

    if path == "/api/motifs":
        return {
            "order": library_order(),
            "motifs": [lib[m].to_json() for m in library_order()],
            "extra_elements": [
                {"id": CPG_SEGMENT_ID, "name": "CpG-rich segment",
                 "kind": "segment", "color": "#7f8c8d",
                 "notes": "A stretch of CpG-island-like sequence. Length, GC "
                          "and CpG o/e are set per placement."},
                {"id": CUSTOM_ID, "name": "Custom sequence", "kind": "segment",
                 "color": "#95a5a6",
                 "notes": "Any literal sequence you paste, so Mode 1 can slide "
                          "an arbitrary segment relative to the TSS."},
            ],
            "provenance": ("Consensus strings and PWMs are hand-written "
                           "placeholders, not Puffin's learned motifs. Drop "
                           "real PWMs into data/puffin_motifs.npz to replace "
                           "them."),
        }

    if path == "/api/models":
        return adapters.list_models()

    if path == "/api/metrics":
        return {"metrics": scoring.metric_catalogue(),
                "default_weights": {"design": scoring.DEFAULT_WEIGHTS_DESIGN,
                                    "mode1": scoring.DEFAULT_WEIGHTS_MODE1}}

    if path == "/api/celltypes":
        model = query.get("model", [adapters.DEFAULT_CELLTYPE_MODEL])[0]
        a = adapters.get_celltype_adapter(model)
        return {"model": model, "is_mock": a.is_mock,
                "cell_types": a.list_cell_types()}

    # --- sequence / design -------------------------------------------------
    if path == "/api/background":
        bg = random_background(
            length=int(body.get("length", 1001)),
            tss_index=body.get("tss_index"),
            gc=float(body.get("gc", 0.45)),
            cpg_oe=float(body.get("cpg_oe", 0.25)),
            seed=body.get("seed"))
        bg["hash"] = sequence_hash(bg["sequence"])
        return bg

    if path == "/api/evaluate":
        return _evaluate(body)

    # --- mode 1 ------------------------------------------------------------
    if path == "/api/mode1/scan":
        offsets = body.get("offsets")
        if not offsets:
            lo = int(body.get("offset_min", -300))
            hi = int(body.get("offset_max", 100))
            step = max(1, int(body.get("offset_step", 10)))
            offsets = list(range(lo, hi + 1, step))
        bgs = _backgrounds_from(body)
        kwargs = dict(
            element=body.get("element") or {"element_id": "ets", "strand": "+"},
            offsets=offsets,
            model=body.get("celltype_model") or adapters.DEFAULT_CELLTYPE_MODEL,
            cell_types=body.get("cell_types"),
            target=body.get("target_cell_type"),
            activity_window=tuple(body.get("activity_window") or (-200, 200)),
            activity_method=body.get("activity_method", "mean"),
            window=_win(body))
        if len(bgs) > 1:
            fn, args = experiments.position_scan_multi, (bgs,)
        else:
            fn, args = experiments.position_scan, (bgs[0],)
        if body.get("async"):
            return {"job_id": start_job(fn, *args, **kwargs)}
        return fn(*args, progress=None, **kwargs)

    # --- mode 2 ------------------------------------------------------------
    if path == "/api/mode2/pair":
        bgs = _backgrounds_from(body)
        kwargs = dict(
            motif_a=body.get("motif_a", "tata"),
            motif_b=body.get("motif_b", "inr"),
            anchor=int(body.get("anchor", -60)),
            spacing=int(body.get("spacing", 30)),
            spacing_mode=body.get("spacing_mode", "center"),
            strand_a=body.get("strand_a", "+"),
            strand_b=body.get("strand_b", "+"),
            models=body.get("profile_models") or adapters.DEFAULT_PROFILE_MODELS,
            window=_win(body),
            track=body.get("track", "plus"),
            residual_scale=body.get("residual_scale", "model"))
        if len(bgs) > 1:
            return experiments.pair_interaction_multi(bgs, **kwargs)
        return experiments.pair_interaction(bgs[0], **kwargs)

    if path == "/api/mode2/grid":
        bgs = _backgrounds_from(body)
        kwargs = dict(
            motif_ids=body.get("motif_ids"),
            anchor=int(body.get("anchor", -60)),
            spacing=int(body.get("spacing", 30)),
            spacing_mode=body.get("spacing_mode", "center"),
            strand_a=body.get("strand_a", "+"),
            strand_b=body.get("strand_b", "+"),
            models=body.get("profile_models") or adapters.DEFAULT_PROFILE_MODELS,
            statistic=body.get("statistic", "I_sum_abs"),
            window=_win(body),
            track=body.get("track", "plus"),
            residual_scale=body.get("residual_scale", "model"))
        if body.get("async", True):
            return {"job_id": start_job(experiments.pair_grid, bgs, **kwargs)}
        return experiments.pair_grid(bgs, progress=None, **kwargs)

    if path == "/api/mode2/spacing":
        bgs = _backgrounds_from(body)
        spacings = body.get("spacings")
        if not spacings:
            lo = int(body.get("spacing_min", 4))
            hi = int(body.get("spacing_max", 84))
            step = max(1, int(body.get("spacing_step", 2)))
            spacings = list(range(lo, hi + 1, step))
        orientations = [tuple(o) if isinstance(o, (list, tuple)) else (o[0], o[1])
                        for o in (body.get("orientations") or [["+", "+"]])]
        kwargs = dict(
            motif_a=body.get("motif_a", "tata"),
            motif_b=body.get("motif_b", "inr"),
            spacings=spacings,
            anchor=int(body.get("anchor", -60)),
            spacing_mode=body.get("spacing_mode", "center"),
            orientations=orientations,
            models=body.get("profile_models") or adapters.DEFAULT_PROFILE_MODELS,
            statistic=body.get("statistic", "I_sum_abs"),
            window=_win(body),
            track=body.get("track", "plus"),
            residual_scale=body.get("residual_scale", "model"))
        if body.get("async", True):
            return {"job_id": start_job(experiments.spacing_curve, bgs, **kwargs)}
        return experiments.spacing_curve(bgs, progress=None, **kwargs)

    # --- jobs --------------------------------------------------------------
    if path.startswith("/api/jobs/"):
        jid = path.rsplit("/", 1)[-1]
        j = job_status(jid)
        if j is None:
            raise KeyError(f"no such job: {jid}")
        return j

    # --- multiplayer -------------------------------------------------------
    if path == "/api/players":
        if body:
            return store.ensure_player(body.get("name", ""), body.get("player_id"))
        return {"players": store.list_players()}

    if path == "/api/challenges":
        if body:
            return store.create_challenge(
                name=body.get("name", ""),
                background_params=body.get("background_params") or {},
                mode=body.get("mode", "design"),
                target_cell_type=body.get("target_cell_type"),
                weights=body.get("weights"),
                models=body.get("models"),
                created_by=body.get("player_id"),
                notes=body.get("notes", ""))
        return {"challenges": store.list_challenges()}

    if path == "/api/challenges/submit":
        ch = store.get_challenge(body.get("challenge_id", ""))
        if not ch:
            raise KeyError("unknown challenge")
        # Re-evaluate server-side from the challenge's own parameters, so a
        # leaderboard entry cannot be produced by a client sending its own
        # numbers.
        ev = _evaluate({
            **body,
            "background_params": ch["background_params"],
            "weights": ch.get("weights") or scoring.DEFAULT_WEIGHTS_DESIGN,
            "target_cell_type": ch.get("target_cell_type"),
            **(ch.get("models") or {}),
        })
        player = store.ensure_player(body.get("player_name", "anonymous"),
                                     body.get("player_id"))
        sub = store.submit(ch["id"], player["id"], body.get("placements") or [],
                           ev["metrics"], ev["score"], body.get("label", ""),
                           extra={"sequence_hash": ev["sequence_hash"]})
        return {"submission": sub, "evaluation": ev,
                "leaderboard": store.leaderboard(ch["id"])}

    if path == "/api/challenges/delete":
        cid = body.get("challenge_id", "")
        if not store.delete_challenge(cid):
            raise KeyError("unknown challenge")
        return {"deleted": cid}

    if path == "/api/challenges/leaderboard":
        cid = (query.get("challenge_id") or [body.get("challenge_id", "")])[0]
        return {"challenge": store.get_challenge(cid),
                "leaderboard": store.leaderboard(
                    cid, best_per_player=body.get("best_per_player", True))}

    # --- export ------------------------------------------------------------
    if path == "/api/export/bundle":
        os.makedirs(RUNS_DIR, exist_ok=True)
        name = (body.get("name") or "run").replace("/", "_")[:60]
        stamp = time.strftime("%Y%m%d-%H%M%S")
        fname = f"{stamp}_{name}.json"
        full = os.path.join(RUNS_DIR, fname)
        payload = {
            "saved_at": time.time(),
            "saved_at_iso": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "app": "promoter-design-game",
            "models": adapters.list_models(),
            "payload": body.get("payload"),
            "settings": body.get("settings"),
            "note": body.get("note", ""),
        }
        with open(full, "w") as fh:
            json.dump(payload, fh, cls=JSONEncoder, indent=1)
        return {"saved": full, "filename": fname,
                "bytes": os.path.getsize(full)}

    if path == "/api/export/fasta":
        bg = _background_from(body)
        con = build_construct(bg["sequence"], bg["tss_index"],
                              to_placements(body.get("placements") or []), seed=0)
        name = body.get("name") or f"construct_{sequence_hash(con['sequence'])}"
        header = (f"{name} len={con['length']} tss_index={con['tss_index']} "
                  f"gc={con['stats']['gc']} cpg_oe={con['stats']['cpg_oe']}")
        return {"fasta": to_fasta(con["sequence"], header),
                "sequence": con["sequence"],
                "hash": sequence_hash(con["sequence"])}

    if path == "/api/runs":
        os.makedirs(RUNS_DIR, exist_ok=True)
        files = sorted(os.listdir(RUNS_DIR), reverse=True)[:100]
        return {"dir": RUNS_DIR,
                "runs": [{"filename": f,
                          "bytes": os.path.getsize(os.path.join(RUNS_DIR, f))}
                         for f in files if f.endswith(".json")]}

    if path == "/api/cache/clear":
        adapters.clear_cache()
        return {"cleared": True}

    raise FileNotFoundError(f"no route {path}")


# ---------------------------------------------------------------------------
# HTTP plumbing
# ---------------------------------------------------------------------------

class Handler(BaseHTTPRequestHandler):
    server_version = "PromoterDesignGame/0.1"
    protocol_version = "HTTP/1.1"

    def log_message(self, fmt, *args):
        if os.environ.get("PDG_QUIET"):
            return
        sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))

    # -- helpers ------------------------------------------------------------
    def _send(self, code: int, payload: bytes, ctype: str,
              extra: dict | None = None):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(payload)))
        self.send_header("Cache-Control", "no-store")
        for k, v in (extra or {}).items():
            self.send_header(k, v)
        self.end_headers()
        try:
            self.wfile.write(payload)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def _json(self, obj, code: int = 200):
        self._send(code, json.dumps(obj, cls=JSONEncoder).encode(),
                   "application/json; charset=utf-8")

    def _error(self, code: int, msg: str, detail: str = ""):
        self._json({"error": msg, "detail": detail}, code)

    def _static(self, path: str):
        rel = path.lstrip("/") or "index.html"
        full = os.path.normpath(os.path.join(FRONTEND, rel))
        if not full.startswith(os.path.abspath(FRONTEND)):
            return self._error(403, "forbidden")
        if os.path.isdir(full):
            full = os.path.join(full, "index.html")
        if not os.path.exists(full):
            return self._error(404, "not found", rel)
        ctype = mimetypes.guess_type(full)[0] or "application/octet-stream"
        if ctype.startswith("text/") or ctype in ("application/javascript",
                                                  "application/json"):
            ctype += "; charset=utf-8"
        with open(full, "rb") as fh:
            self._send(200, fh.read(), ctype)

    def _dispatch(self, body: dict):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        query = urllib.parse.parse_qs(parsed.query)
        if not path.startswith("/api/"):
            return self._static(path)
        try:
            t0 = time.time()
            result = route(path, body, query)
            if isinstance(result, dict):
                result.setdefault("_elapsed_s", round(time.time() - t0, 4))
            self._json(result)
        except FileNotFoundError as exc:
            self._error(404, str(exc))
        except KeyError as exc:
            self._error(400, str(exc).strip("'"))
        except NotImplementedError as exc:
            self._error(501, "model not wired up", str(exc))
        except ValueError as exc:
            self._error(400, "bad request", str(exc))
        except Exception as exc:
            traceback.print_exc()
            self._error(500, f"{type(exc).__name__}: {exc}",
                        traceback.format_exc()[-2000:])

    # -- verbs --------------------------------------------------------------
    def do_GET(self):
        self._dispatch({})

    def do_POST(self):
        try:
            n = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            n = 0
        if n > MAX_BODY:
            return self._error(413, "body too large")
        raw = self.rfile.read(n) if n else b""
        try:
            body = json.loads(raw.decode() or "{}")
        except json.JSONDecodeError as exc:
            return self._error(400, "invalid JSON", str(exc))
        if not isinstance(body, dict):
            return self._error(400, "body must be a JSON object")
        self._dispatch(body)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--host", default="0.0.0.0")
    ap.add_argument("--open-port-scan", action="store_true",
                    help="if --port is busy, try the next 20 ports")
    args = ap.parse_args()

    os.makedirs(RUNS_DIR, exist_ok=True)

    port = args.port
    srv = None
    last = None
    for candidate in range(port, port + (20 if args.open_port_scan else 1)):
        try:
            srv = ThreadingHTTPServer((args.host, candidate), Handler)
            port = candidate
            break
        except OSError as exc:
            last = exc
            continue
    if srv is None:
        raise SystemExit(f"could not bind {args.host}:{args.port}: {last}")

    srv.daemon_threads = True
    host = socket.gethostname()
    models = adapters.list_models()
    n_mock = sum(1 for m in models["profile"] + models["celltype"] if m["is_mock"])

    print("=" * 72)
    print("  Promoter design game")
    print(f"  serving  http://{host}:{port}/   (bound {args.host}:{port})")
    print(f"  frontend {FRONTEND}")
    print(f"  runs     {RUNS_DIR}")
    print(f"  models   {n_mock} mock adapter(s) registered; "
          f"real inference available: {models['any_real_available']}")
    print()
    print("  From your laptop:")
    print(f"    ssh -N -L {port}:{host}:{port} Randi")
    print(f"    open http://localhost:{port}/")
    print("=" * 72)
    sys.stdout.flush()

    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        print("\nshutting down")
        srv.shutdown()


if __name__ == "__main__":
    main()
