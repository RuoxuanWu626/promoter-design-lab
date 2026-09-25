"""Probe the local AlphaGenome PyTorch port: what lengths work, how fast, how big.

Run on a GPU node. The answers decide how the Mode 1 adapter has to be built —
a model that only accepts 1 Mb inputs needs the construct embedded in a large
context, and the cost per offset sets how big a position scan can be.
"""

import os
import sys
import time

import numpy as np
import torch

TORCH_DIR = "/gpfs/data/zhou-lab/rxwu/Puffin-drozofila/tests/alphagenome_torch"
STATE = f"{TORCH_DIR}/alphagenome_all_folds.pt"
META = ("/gpfs/data/zhou-lab/sharish/genomics/alphagenome_dataset/hg38/bp/"
        "alphagenome_human_bp_stranded_cp100m_logp1.v3.meta.tsv")


def load_model(device):
    sys.path.insert(0, TORCH_DIR)
    from alphagenome.models import dna_model as dna_model_lib
    from alphagenome.models import dna_output
    from alphagenome_research.model.metadata import metadata as metadata_lib
    from alphagenome_torch import model as torch_model_lib

    organism = dna_model_lib.Organism.HOMO_SAPIENS
    md = metadata_lib.load(organism)
    model = torch_model_lib.AlphaGenomeTorch({organism: md})
    full = torch.load(STATE, map_location="cpu", weights_only=True)
    state = {}
    for name, target in model.state_dict().items():
        src = full[name]
        if src.shape == target.shape:
            state[name] = src
        elif (src.ndim > 0 and target.ndim > 0 and src.shape[0] == 2
              and target.shape[0] == 1 and src.shape[1:] == target.shape[1:]):
            state[name] = src[:1]
        else:
            raise RuntimeError(f"cannot convert {name}: {tuple(src.shape)} -> {tuple(target.shape)}")
    model.load_state_dict(state, strict=True)
    model.eval().to(device)
    del full, state
    valid = np.flatnonzero(~np.asarray(md.padding[dna_output.OutputType.CAGE]))
    return model, valid


def main():
    print("torch", torch.__version__, "cuda", torch.cuda.is_available(), flush=True)
    if torch.cuda.is_available():
        print("gpu:", torch.cuda.get_device_name(0), flush=True)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    t0 = time.time()
    model, valid = load_model(device)
    print(f"model loaded in {time.time()-t0:.0f}s; valid CAGE tracks = {len(valid)}", flush=True)
    valid_t = torch.as_tensor(valid, dtype=torch.long, device=device)

    # --- which input lengths does it accept? --------------------------------
    print("\n=== supported sequence lengths ===", flush=True)
    ok_lengths = []
    for L in (2048, 16384, 131072, 524288, 1048576):
        try:
            enc = np.zeros((L, 4), dtype=np.float32)
            enc[np.arange(L), np.random.randint(0, 4, L)] = 1.0
            x = torch.from_numpy(enc)[None].to(device)
            torch.cuda.synchronize() if device.type == "cuda" else None
            t = time.time()
            with torch.inference_mode(), torch.autocast(device_type=device.type, dtype=torch.float16):
                raw, _ = model(x, torch.zeros(1, dtype=torch.long, device=device),
                               output_heads=("cage",))
                p1 = raw["cage"]["predictions_1bp"][0].index_select(-1, valid_t)
            torch.cuda.synchronize() if device.type == "cuda" else None
            dt = time.time() - t
            mem = (torch.cuda.max_memory_allocated() / 2**30) if device.type == "cuda" else 0
            print(f"  L={L:>8d}  OK  out={tuple(p1.shape)}  {dt:6.2f}s  peak {mem:5.2f} GiB", flush=True)
            ok_lengths.append((L, dt))
            if device.type == "cuda":
                torch.cuda.reset_peak_memory_stats()
        except Exception as e:
            print(f"  L={L:>8d}  FAILED  {type(e).__name__}: {str(e)[:110]}", flush=True)

    # --- cell-type metadata -------------------------------------------------
    print("\n=== CAGE track metadata ===", flush=True)
    try:
        import pandas as pd
        meta = pd.read_csv(META, sep="\t")
        cage = meta[meta.output_type == "CAGE"].reset_index(drop=True)
        names = cage["name"].tolist()
        cts = list(dict.fromkeys(names))
        print(f"  {len(cage)} tracks, {len(cts)} distinct cell types", flush=True)
        print(f"  strands: {sorted(set(cage['strand']))}", flush=True)
        print("  first 12 cell types:", flush=True)
        for c in cts[:12]:
            print(f"    {c}", flush=True)
    except Exception as e:
        print(f"  metadata unavailable: {type(e).__name__}: {e}", flush=True)

    # --- does a local edit stay local? --------------------------------------
    if ok_lengths:
        L = min(l for l, _ in ok_lengths)
        print(f"\n=== locality check at L={L} ===", flush=True)
        rng = np.random.default_rng(0)
        base = rng.integers(0, 4, L)

        def run(idx):
            enc = np.zeros((L, 4), dtype=np.float32)
            enc[np.arange(L), idx] = 1.0
            x = torch.from_numpy(enc)[None].to(device)
            with torch.inference_mode(), torch.autocast(device_type=device.type, dtype=torch.float16):
                raw, _ = model(x, torch.zeros(1, dtype=torch.long, device=device),
                               output_heads=("cage",))
                return raw["cage"]["predictions_1bp"][0].index_select(-1, valid_t).float().cpu().numpy()

        a = run(base)
        mid = L // 2
        edited = base.copy()
        edited[mid - 4: mid + 5] = np.array([3, 2, 0, 3, 0, 0, 0, 0, 2])  # a TATA-ish edit
        b = run(edited)
        d = np.abs(a - b).sum(axis=1)
        nz = np.flatnonzero(d > d.max() * 0.01)
        print(f"  changed 9 bp at the centre; |delta| exceeds 1% of its max over "
              f"{len(nz)} of {L} positions", flush=True)
        if len(nz):
            print(f"  affected span: {nz.min()-mid:+d} .. {nz.max()-mid:+d} bp around the edit", flush=True)


if __name__ == "__main__":
    main()
