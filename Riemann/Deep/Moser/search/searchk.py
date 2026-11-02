#!/usr/bin/env python3
"""Certificate search driver for k-gon covers (adapted from scripts/upper/generate.py of
ethan-keller/moser-worm-improved-bounds). Low-memory master: leaf multipliers are spilled
to OUT.leaves.jsonl, pending states kept as JSON strings; checkpoints every 5 minutes.

Usage: searchk.py CONFIG.json OUT [--jobs 3] [--max-seconds S] [--adaptive]
       [--max-unresolved N] [--init-tree CERTDIR] [--resume]
CONFIG: {"t": [...], "h": [...], "specs": [...]}
On success writes OUT/ in the compact certificate format read by verify_upper.py.
"""
import argparse
from collections import deque
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import gzip
import hashlib
import json
import lzma
from multiprocessing import get_context
import os
from pathlib import Path
import signal
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[k] = "1"

WORKER = None
MAX_DEPTH = 128


def worker_init(T, H, specs, adaptive):
    global WORKER
    import warnings
    warnings.filterwarnings("ignore")
    sys.stdout = sys.stderr = open(os.devnull, "w")
    from modelk import Model
    WORKER = Model(T, H, specs, adaptive=adaptive)


def worker_expand(task):
    state_str, extra = task
    assert WORKER is not None, "worker not initialised"
    WORKER.sync(extra)
    dec = WORKER.expand(json.loads(state_str), existing_specs=WORKER.geometry.specs)
    if dec[0] == "l":
        pairs = sorted(dec[2])
        keys = [k for k, _ in pairs]
        dec = ["l", dec[1], [keys[0]] + [b - a for a, b in zip(keys, keys[1:])], [v for _, v in pairs]]
    elif dec[0] == "p":
        dec = ["p", dec[1], *[json.dumps(s, separators=(",", ":")) for s in dec[2:]]]
    elif dec[0] == "s":
        dec = ["s", [[rect, json.dumps(st, separators=(",", ":"))] for rect, st in dec[1]]]
    elif dec[0] == "unresolved":
        dec = dec[:3]
    return dec


class InlinePool:
    """jobs == 1: expand in this process (no worker pool; one single-threaded process)."""

    def __init__(self, T, H, specs, adaptive):
        import warnings
        warnings.filterwarnings("ignore")
        global WORKER
        from modelk import Model
        WORKER = Model(T, H, specs, adaptive=adaptive)
        self._processes = {}

    def submit(self, fn, *args):
        from concurrent.futures import Future
        fut = Future()
        try:
            fut.set_result(fn(*args))
        except KeyboardInterrupt:
            raise
        except Exception as exc:
            fut.set_exception(exc)
        return fut

    def shutdown(self, wait=True, cancel_futures=False):
        pass


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_certificate(directory, header, nodes, spill_path):
    leafno = {}
    sys.setrecursionlimit(1000000)

    def visit(i):
        nd = nodes[i]
        if nd[0] == "l":
            leafno[i] = len(leafno)
            return ["l", leafno[i]]
        if nd[0] == "p":
            return ["p", nd[1], *[visit(c) for c in nd[2:]]]
        if nd[0] == "s":
            return ["s", [[list(r), visit(c)] for r, c in nd[1]]]
        raise ValueError("incomplete tree")

    tree = visit(0)
    directory.mkdir(parents=True, exist_ok=False)
    hdr = {**header, "leaves": len(leafno), "complete": True, "schema": 1}
    (directory / "header.json").write_text(json.dumps(hdr, sort_keys=True, indent=1) + "\n")
    with gzip.open(directory / "tree.json.gz", "wt", compresslevel=9) as f:
        json.dump(tree, f, separators=(",", ":"))
    del tree
    files, batch, seen = [], [], set()

    def flush():
        if not batch:
            return
        name = "multipliers-%03d.jsonl.xz" % len(files)
        data = ("\n".join(json.dumps(x, separators=(",", ":")) for x in batch) + "\n").encode()
        (directory / name).write_bytes(lzma.compress(data, preset=9 | lzma.PRESET_EXTREME))
        files.append(name)
        batch.clear()

    with open(spill_path) as f:
        for line in f:
            i, den, dk, vals = json.loads(line)
            if i in leafno and i not in seen:
                seen.add(i)
                batch.append([leafno[i], den, dk, vals])
                if len(batch) >= 4096:
                    flush()
    flush()
    if len(seen) != len(leafno):
        raise ValueError("missing leaf records in spill file")
    man = {"schema": 2, "format": "compact", "files": files,
           "sha256": {n: sha(directory / n) for n in ["header.json", "tree.json.gz", *files]}}
    (directory / "manifest.json").write_text(json.dumps(man, sort_keys=True, indent=1) + "\n")
    return len(leafno)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("out")
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--max-seconds", type=float, default=3 * 3600)
    ap.add_argument("--adaptive", action="store_true")
    ap.add_argument("--max-unresolved", type=int, default=1)
    ap.add_argument("--max-new-specs", type=int, default=4000)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--init-tree", default=None)
    args = ap.parse_args()
    cfg = json.load(open(args.config))
    T, H, base = cfg["t"], cfg["h"], cfg["specs"]
    extra = []
    out = Path(args.out)
    log = open(str(out) + ".log", "a", buffering=1)
    spill_path = Path(str(out) + ".leaves.jsonl")
    ckpt = Path(str(out) + ".ckpt.json.gz")
    started = time.monotonic()
    from exactk import Geometry
    g0 = Geometry(T, H, base)
    kR0, kL0 = g0.kR, g0.kL
    nalt = len(T) - 2
    del g0
    J = lambda s: json.dumps(s, separators=(",", ":"))
    nodes = [["P", J([[], [0, kR0, 0, kL0], [], []]), 0]]
    newspecs = 0
    if args.resume and ckpt.exists():
        with gzip.open(ckpt, "rt") as f:
            data = json.load(f)
        nodes, extra, newspecs = data["nodes"], data["extra"], data["newspecs"]
        with open(spill_path, "r+") as f:
            f.truncate(data["spill_bytes"])
        print("resumed: %d nodes, %d extra specs" % (len(nodes), len(extra)), file=log)
    elif args.init_tree:
        old_hdr = json.load(open(os.path.join(args.init_tree, "header.json")))
        if old_hdr["t"] != [str(x) for x in T] or old_hdr["specs"] != base:
            raise ValueError("init tree must use the same t and specs")
        old_tree = json.load(gzip.open(os.path.join(args.init_tree, "tree.json.gz"), "rt"))
        nodes = []
        sys.setrecursionlimit(1000000)

        def build(old, assigned, rect, depth):
            idx = len(nodes)
            nodes.append(None)
            if old[0] == "l":
                nodes[idx] = ["P", J([assigned, rect, [], []]), depth]
            elif old[0] == "p":
                ch = [build(old[2 + r], assigned + [[old[1], r]], rect, depth + 1) for r in range(nalt)]
                nodes[idx] = ["p", old[1], *ch]
            else:
                ch = [[rect2, build(c, assigned, rect2, depth + 1)] for rect2, c in old[1]]
                nodes[idx] = ["s", ch]
            return idx

        build(old_tree, [], [0, kR0, 0, kL0], 0)
        del old_tree
        if spill_path.exists():
            spill_path.unlink()
        print("init tree: %d nodes" % len(nodes), file=log)
    else:
        if spill_path.exists():
            spill_path.unlink()
    pending = deque(i for i, nd in enumerate(nodes) if nd[0] == "P")
    spill = open(spill_path, "a")
    running = {}
    expanded = 0
    unresolved = []
    last = last_ckpt = time.monotonic()

    def new_pool():
        if args.jobs == 1:
            return InlinePool(T, H, base, args.adaptive)
        return ProcessPoolExecutor(max_workers=args.jobs, mp_context=get_context("spawn"),
                                   initializer=worker_init, initargs=(T, H, base, args.adaptive))

    def save_ckpt():
        spill.flush()
        tmp = Path(str(ckpt) + ".tmp")
        with gzip.open(tmp, "wt", compresslevel=3) as f:
            json.dump({"nodes": nodes, "extra": extra, "newspecs": newspecs, "spill_bytes": spill.tell()}, f,
                      separators=(",", ":"))
        tmp.replace(ckpt)

    def on_signal(signum, frame):
        raise KeyboardInterrupt("signal %d" % signum)

    signal.signal(signal.SIGTERM, on_signal)
    signal.signal(signal.SIGALRM, on_signal)
    pool = new_pool()
    crashes = {}
    status = "complete"
    try:
        while pending or running:
            if time.monotonic() - started > args.max_seconds:
                status = "timeout"
                break
            if len(unresolved) >= args.max_unresolved:
                status = "unresolved"
                break
            while pending and len(running) < args.jobs:
                i = pending.popleft()
                if nodes[i][2] > MAX_DEPTH:
                    nodes[i] = ["u", "depth"]
                    unresolved.append([i, "depth"])
                    continue
                running[pool.submit(worker_expand, (nodes[i][1], extra))] = i
            if not running:
                break
            ready, _ = wait(running, timeout=1, return_when=FIRST_COMPLETED)
            broken = False
            for fut in ready:
                i = running.pop(fut)
                try:
                    dec = fut.result()
                except Exception as exc:
                    print("worker failure on node %d: %r" % (i, exc), file=log)
                    crashes[i] = crashes.get(i, 0) + 1
                    if crashes[i] >= 3:
                        nodes[i] = ["u", "crash"]
                        unresolved.append([i, "crash"])
                    else:
                        pending.appendleft(i)
                    broken = True
                    continue
                expanded += 1
                depth = nodes[i][2]

                def child(state_str):
                    nodes.append(["P", state_str, depth + 1])
                    pending.append(len(nodes) - 1)
                    return len(nodes) - 1

                if dec[0] == "p":
                    nodes[i] = ["p", dec[1], *[child(s) for s in dec[2:]]]
                elif dec[0] == "s":
                    nodes[i] = ["s", [[rect, child(st)] for rect, st in dec[1]]]
                elif dec[0] == "l":
                    spill.write(json.dumps([i, dec[1], dec[2], dec[3]], separators=(",", ":")) + "\n")
                    nodes[i] = ["l"]
                elif dec[0] == "newspec":
                    spec = dec[1]
                    if spec in extra or spec in base:
                        pending.appendleft(i)
                    elif newspecs < args.max_new_specs:
                        extra.append(spec)
                        newspecs += 1
                        pending.appendleft(i)
                    else:
                        nodes[i] = ["u", "spec cap"]
                        unresolved.append([i, "spec cap"])
                elif dec[0] == "unresolved":
                    unresolved.append([i, dec[1:], nodes[i][1]])
                    nodes[i] = ["u", dec[1]]
                else:
                    raise ValueError(dec)
            if broken:
                for j in list(running.values()):
                    pending.appendleft(j)
                running.clear()
                try:
                    pool.shutdown(wait=False, cancel_futures=True)
                except Exception:
                    pass
                pool = new_pool()
                continue
            if time.monotonic() - last_ckpt >= 300:
                save_ckpt()
                last_ckpt = time.monotonic()
            if time.monotonic() - last >= 60:
                print("%.1f min: %d expanded, %d pending, %d running, %d new specs, %d unresolved" % (
                    (time.monotonic() - started) / 60, expanded, len(pending), len(running), newspecs, len(unresolved)), file=log)
                last = time.monotonic()
    except KeyboardInterrupt:
        status = "interrupted"
    finally:
        procs = list((getattr(pool, "_processes", None) or {}).values())
        for f in running:
            f.cancel()
        for p in procs:
            if p.is_alive():
                p.terminate()
        try:
            pool.shutdown(wait=False, cancel_futures=True)
        except Exception:
            pass
    for i in running.values():
        pass  # running nodes stay "P" in nodes, so they are pending on resume
    spill.flush()
    if status != "complete" or pending:
        save_ckpt()
    if status == "complete" and (unresolved or pending):
        status = "unresolved"
    msg = "%s: %d nodes expanded in %.1f s; %d new specs; %d unresolved" % (
        status, expanded, time.monotonic() - started, newspecs, len(unresolved))
    print(msg, file=log)
    print(msg)
    if status == "complete":
        spill.close()
        header = {"t": [str(x) for x in T], "h": [str(x) for x in H], "specs": base + extra, "kR": kR0, "kL": kL0}
        nleaves = write_certificate(out, header, nodes, spill_path)
        print("certificate written: %d leaves" % nleaves, file=log)
    else:
        json.dump({"status": status, "unresolved": unresolved[:50], "extra": extra, "expanded": expanded},
                  open(str(out) + ".failed.json", "w"))


if __name__ == "__main__":
    main()

