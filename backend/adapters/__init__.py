"""Adapter registry and a small prediction cache.

Adding a model is: write a class implementing ``ProfileAdapter`` or
``CellTypeAdapter``, import it here, and append it to ``_PROFILE`` /
``_CELLTYPE``.  Nothing else in the app needs to know it exists.
"""

from __future__ import annotations

import hashlib
import threading
from collections import OrderedDict

from adapters.base import (
    CellTypeAdapter,
    CellTypePrediction,
    ProfileAdapter,
    ProfilePrediction,
)
from adapters.mock_puffin import MockPuffin, MockPuffinD, MockPuffinAdditiveOnly
from adapters.mock_alphagenome import MockAlphaGenome
from adapters.real_puffin import RealPuffin, RealPuffinD
from adapters.real_alphagenome import RealAlphaGenome

_PROFILE: "OrderedDict[str, ProfileAdapter]" = OrderedDict()
_CELLTYPE: "OrderedDict[str, CellTypeAdapter]" = OrderedDict()

for _a in (MockPuffin(), MockPuffinD(), MockPuffinAdditiveOnly(),
           RealPuffin(), RealPuffinD()):
    _PROFILE[_a.name] = _a

for _c in (MockAlphaGenome(), RealAlphaGenome()):
    _CELLTYPE[_c.name] = _c

DEFAULT_PROFILE_MODELS = ["mock_puffin", "mock_puffind"]
DEFAULT_CELLTYPE_MODEL = "mock_alphagenome"


def get_profile_adapter(name: str) -> ProfileAdapter:
    if name not in _PROFILE:
        raise KeyError(f"unknown profile model '{name}'; "
                       f"have {list(_PROFILE)}")
    return _PROFILE[name]


def get_celltype_adapter(name: str) -> CellTypeAdapter:
    if name not in _CELLTYPE:
        raise KeyError(f"unknown cell-type model '{name}'; "
                       f"have {list(_CELLTYPE)}")
    return _CELLTYPE[name]


def list_models() -> dict:
    return {
        "profile": [a.info() for a in _PROFILE.values()],
        "celltype": [c.info() for c in _CELLTYPE.values()],
        "defaults": {
            "profile": DEFAULT_PROFILE_MODELS,
            "celltype": DEFAULT_CELLTYPE_MODEL,
        },
        "any_real_available": any(
            (not a.is_mock) and a.available()[0]
            for a in list(_PROFILE.values()) + list(_CELLTYPE.values())
        ),
    }


# ---------------------------------------------------------------------------
# Prediction cache.  Experiments re-predict the same background many times
# (every pair in Mode 2 shares it), so this is a large win and it also keeps
# repeated runs bit-identical.
# ---------------------------------------------------------------------------

_CACHE: "OrderedDict[str, object]" = OrderedDict()
_CACHE_LOCK = threading.Lock()
_CACHE_MAX = 4096


def cache_key(*parts) -> str:
    h = hashlib.sha1()
    for p in parts:
        h.update(str(p).encode())
        h.update(b"\x1f")
    return h.hexdigest()


def cached(key: str, compute):
    with _CACHE_LOCK:
        if key in _CACHE:
            _CACHE.move_to_end(key)
            return _CACHE[key]
    value = compute()
    with _CACHE_LOCK:
        _CACHE[key] = value
        _CACHE.move_to_end(key)
        while len(_CACHE) > _CACHE_MAX:
            _CACHE.popitem(last=False)
    return value


def cache_stats() -> dict:
    with _CACHE_LOCK:
        return {"entries": len(_CACHE), "max": _CACHE_MAX}


def clear_cache() -> None:
    with _CACHE_LOCK:
        _CACHE.clear()


def predict_profile(model: str, sequence: str, tss_index: int,
                    window: tuple[int, int]) -> ProfilePrediction:
    adapter = get_profile_adapter(model)
    key = cache_key("profile", model, tss_index, window, sequence)
    return cached(key, lambda: adapter.predict(sequence, tss_index, window))


def predict_celltype(model: str, sequence: str, tss_index: int,
                     cell_types: list[str] | None,
                     window: tuple[int, int]) -> CellTypePrediction:
    adapter = get_celltype_adapter(model)
    key = cache_key("celltype", model, tss_index, window,
                    ",".join(cell_types or []), sequence)
    return cached(key, lambda: adapter.predict(sequence, tss_index,
                                               cell_types, window))


__all__ = [
    "ProfileAdapter", "CellTypeAdapter", "ProfilePrediction", "CellTypePrediction",
    "get_profile_adapter", "get_celltype_adapter", "list_models",
    "predict_profile", "predict_celltype", "cache_key", "cached",
    "cache_stats", "clear_cache",
    "DEFAULT_PROFILE_MODELS", "DEFAULT_CELLTYPE_MODEL",
]
