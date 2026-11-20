#!/usr/bin/env python3
# alpdep: lazy BFS over the name-level graph. Processes (parses + resolves) only files that contain reached decls.
# Resolution = build_graph.py's parser and namespace/open resolution; dot-field fallback only for the first
# component after the longest resolved prefix (or after an unresolved head = local variable), as in Lean.
# usage: lazy_bfs.py SRC OUT roots(comma) maxdepth
import sys, os, time, array
import numpy as np
SRCF = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'build_graph.py')
_src = open(SRCF).read(); _src = _src[:_src.index('def files_iter():')]
ROOTS = sys.argv[3].split(','); MAXD = int(sys.argv[4])
sys.argv = [sys.argv[0], sys.argv[1], sys.argv[2], 'lib']
exec(compile(_src, SRCF, 'exec'))
t0 = time.time()
hs = np.load(OUT + '/h_sorted.npy'); ids = np.load(OUT + '/id_sorted.npy'); hid = np.load(OUT + '/hid.npy')
kind = np.load(OUT + '/kind.npy'); nper = np.load(OUT + '/nper.npy'); last = np.load(OUT + '/last.npy')
NSS = set(int(x) for x in np.load(OUT + '/ns.npy'))
TYPEK = np.zeros(256, dtype=bool)
for c in 'SCIDAO': TYPEK[ord(c)] = True
nl = len(last); N = int(nper.sum())
first = np.zeros(len(nper) + 1, dtype=np.int64); np.cumsum(nper, out=first[1:])
paths = [ln.rstrip('\n') for ln in open('/Users/adriansager/playground/riemann/.scratch/alpcore/files.txt')]
def lookh(cands):
    if not cands: return set(), []
    a = np.fromiter(cands, dtype=np.uint64, count=len(cands))
    lo = np.searchsorted(hs, a, 'left'); hi = np.searchsorted(hs, a, 'right')
    m = hi > lo
    out = []
    for l_, h_ in zip(lo[m].tolist(), hi[m].tolist()): out.extend(ids[l_:h_].tolist())
    return set(a[m].tolist()), out
def inlast(strs):
    if not strs: return set()
    L_ = list(strs); a = np.fromiter((hc(s) for s in L_), dtype=np.uint64, count=len(L_))
    pos = np.searchsorted(last, a); pos[pos >= nl] = nl - 1
    ok = last[pos] == a
    return set(s for s, o in zip(L_, ok.tolist()) if o)
phcache = {}
def file_edges(fidx):
    try: txt = open(os.path.join(SRC, paths[fidx]), encoding='utf-8', errors='replace').read()
    except Exception: txt = ''
    ds = parse(txt, fidx) if txt else []
    assert len(ds) == nper[fidx], (paths[fidx], len(ds), nper[fidx])
    res_all = []; nid = int(first[fidx])
    for (full, kind_, st, en, ns, opens, vars_, al, body) in ds:
        fc = full.split('.')
        ck = (tuple(fc[:-1]), tuple(ns), tuple(opens))
        ph = phcache.get(ck)
        if ph is None:
            ph = [fold(0, ns[:j]) for j in range(len(ns), -1, -1)]
            if len(fc) > len(ns) + 1:
                for j in range(len(fc) - 1, len(ns), -1): ph.append(fold(0, fc[:j]))
            base = list(ph)
            for o in opens:
                oc = o.split('.')
                for b in base:
                    x = fold(b, oc)
                    if x in NSS: ph.append(x)
            ph = list(dict.fromkeys(ph))
            if len(phcache) > 2000: phcache.clear()
            phcache[ck] = ph
        text = '\n'.join(body)
        if vars_: text += '\n' + '\n'.join(vars_)
        toks = set()
        for m in TOK.finditer(text): toks.add((m.group(1), m.group(2)))
        tl = []; comps_all = set()
        for lead, nm in toks:
            comps = unq(nm).split('.')
            if not lead and len(comps) == 1 and comps[0] in KW: continue
            tl.append((lead, comps)); comps_all.update(c for c in comps if not c.isdigit())
        okc = inlast(comps_all)
        cands = set(); pend = set(); tokc = []
        for lead, comps in tl:
            if lead:
                if comps[0] in okc: pend.add(comps[0])
                continue
            phs = ph
            if comps[0] == '_root_':
                comps = comps[1:]; phs = [0]
                if not comps: continue
            if not any(c in okc for c in comps):
                if len(comps) > 1 and comps[1] in okc: pend.add(comps[1])
                continue
            hcs = [hc(c) + 1 for c in comps]; oks = [c in okc for c in comps]
            per = [[] for _ in comps]
            for b in phs:
                h = b
                for k in range(len(comps)):
                    h = (h * MUL + hcs[k]) & MASK
                    if oks[k]: cands.add(h); per[k].append(h)
            tokc.append((comps, per))
        mh, res = lookh(cands); res = set(res)
        for comps, per in tokc:
            kmax = -1
            for k in range(len(comps) - 1, -1, -1):
                if any(h in mh for h in per[k]): kmax = k; break
            nxt = kmax + 1 if kmax >= 0 else 1
            if nxt < len(comps) and comps[nxt] in okc: pend.add(comps[nxt])
        if pend and res:
            th = [int(hid[x]) for x in res if TYPEK[kind[x]]]
            if th and len(th) * len(pend) <= 4000:
                pf = [hc(f) + 1 for f in pend]
                _, r2 = lookh([(hx * MUL + hf) & MASK for hx in th for hf in pf]); res.update(r2)
        res.discard(nid)
        res_all.append((nid, sorted(res))); nid += 1
    return res_all
eoff = np.full(N, -1, dtype=np.int64); elen = np.zeros(N, dtype=np.int32); store = array.array('i')
done = np.zeros(len(nper), dtype=bool)
fid_of = np.searchsorted(first, np.arange(N), 'right') - 1
name_id = {}
with open(OUT + '/names.tsv') as fh:
    for ln in fh:
        p = ln.split('\t', 3)
        if p[1] in ROOTS: name_id.setdefault(p[1], []).append(int(p[0]))
F = sorted(set(i for r in ROOTS for i in name_id.get(r, [])))
print('roots', {r: name_id.get(r) for r in ROOTS}, file=sys.stderr, flush=True)
par = np.full(N, -2, dtype=np.int32); dist = np.full(N, -1, dtype=np.int32)
for i in F: par[i] = -1; dist[i] = 0
d = 0; nfiles = 0
while F and d < MAXD:
    need = sorted(set(int(fid_of[u]) for u in F if not done[fid_of[u]]))
    for fidx in need:
        for (i, r) in file_edges(fidx):
            eoff[i] = len(store); elen[i] = len(r); store.extend(r)
        done[fidx] = True; nfiles += 1
    nxt = []
    for u in F:
        o = int(eoff[u]); l_ = int(elen[u])
        for v in store[o:o + l_]:
            if par[v] == -2: par[v] = u; dist[v] = d + 1; nxt.append(v)
    d += 1; F = nxt
    print('depth', d, 'frontier', len(F), 'reached', int((dist >= 0).sum()), 'files', nfiles, 'edges', len(store), 'secs', round(time.time() - t0), file=sys.stderr, flush=True)
np.save(OUT + '/dist.npy', dist); np.save(OUT + '/par.npy', par)
np.save(OUT + '/eoff.npy', eoff); np.save(OUT + '/elen.npy', elen)
with open(OUT + '/estore.bin', 'wb') as fh: fh.write(store.tobytes())
print('DONE reached', int((dist >= 0).sum()), 'of', N, 'files processed', nfiles, 'of', len(nper), 'secs', round(time.time() - t0), file=sys.stderr)
