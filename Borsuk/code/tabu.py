# Build first: clang -O2 -o Borsuk/.scratch/tabucol Borsuk/code/tabucol.c; clang -O2 -o Borsuk/.scratch/tabucol2 Borsuk/code/tabucol2.c
import numpy as np, subprocess, os
BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '.scratch', 'tabucol')
def colour_tabu(n, U, V, k=5, iters=300000, seed=0, init=None, frozen=None):
    if init is None: init = np.random.default_rng(seed).integers(0, k, n)
    hdr = np.array([n, len(U), k, iters, seed, 0 if frozen is None else 1], np.int32)
    E = np.stack([np.asarray(U, np.int32), np.asarray(V, np.int32)], 1)
    data = hdr.tobytes() + np.asarray(init, np.int32).tobytes() + (b'' if frozen is None else np.asarray(frozen, np.int32).tobytes()) + E.tobytes()
    out = subprocess.run([BIN], input=data, capture_output=True, check=True).stdout
    a = np.frombuffer(out, np.int32); return a[1:].copy(), int(a[0])

BIN2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '.scratch', 'tabucol2')
def colour_smooth(n, FU, FV, NU, NV, k=5, iters=200000, seed=0, W=1000, init=None):
    if init is None: init = np.random.default_rng(seed).integers(0, k, n)
    hdr = np.array([n, len(FU), len(NU), k, iters, seed, W], np.int32)
    Ef = np.stack([np.asarray(FU, np.int32), np.asarray(FV, np.int32)], 1); En = np.stack([np.asarray(NU, np.int32), np.asarray(NV, np.int32)], 1)
    data = hdr.tobytes() + np.asarray(init, np.int32).tobytes() + Ef.tobytes() + En.tobytes()
    out = subprocess.run([BIN2], input=data, capture_output=True, check=True).stdout
    a = np.frombuffer(out, np.int32); return a[2:].copy(), int(a[0]), int(a[1])
