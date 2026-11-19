# scan_zeta.py: rigorous node scan for the zeta certificate (cert_zeta.py), chunked and resumable, with the early-exit scanner
# rscan2.c (proof: RIGOR.md, sec. "Early exit").  Writes the scan cache that cert_zeta.py reads (same format as its own rscan).
# Usage: scan_zeta.py c [K2] [nchunks]     env: WZ_RSCAN2=<rscan2 binary>, WZ_SCAN=<cache dir>, WZ_T1 (default 1)
import sys, os, time, json, math, subprocess
from flint import arb, ctx, fmpq
import mpmath as mpm
from cert_zeta import SymZeta, scan_end, parse_c

def main():
    cq = parse_c(sys.argv[1]); cf = int(cq.p)/int(cq.q)
    K2 = max(201, math.ceil(5*cf*cf/3)); K2 += (K2 % 2 == 0)
    if len(sys.argv) > 2 and sys.argv[2] != '-': K2 = int(sys.argv[2]); K2 += (K2 % 2 == 0)
    nch = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    J = K2; exe = os.environ['WZ_RSCAN2']; cdir = os.environ.get('WZ_SCAN', '.')
    T1s = os.environ.get('WZ_T1', '1'); t1f = float(T1s)
    t0 = time.time()
    ctx.prec = 192
    S2 = SymZeta(cq); pi = arb.pi(); L = arb(cq).log(); two64 = arb(2)**64
    ks = scan_end(S2, K2, t1f); kse = scan_end(S2, 2*J, 0)
    print(f"c={cq} K2={K2} kscan={ks} kscan_e={kse} atoms={len(S2.atoms)} sig2={S2.sig2.mid().str(10)}", flush=True)
    u = 2.0**-53
    def dn(x): return math.nextafter(float(x.lower()), -math.inf)
    def up(x): return math.nextafter(float(x.upper()), math.inf)
    tab = []
    for a in range(65536):
        ang = 2*pi*a/65536
        tab.append(f"{float(ang.cos().mid()).hex()} {float(ang.sin().mid()).hex()}")
    K = float((2*pi/two64).mid())
    mpm.mp.dps = 60; Lm = mpm.log(mpm.mpf(int(cq.p))/mpm.mpf(int(cq.q)))
    at = []; emax = arb(0); dW = arb(0); sigW = 0.0
    for n, _, w in S2.atoms:
        Bn = int(mpm.floor(mpm.log(n)/(2*Lm)*mpm.mpf(2)**64 + mpm.mpf(1)/2))
        assert 0 <= Bn < 2**64
        en = abs(arb(n).log()/(2*L)*two64 - Bn)
        if en.upper() > emax.upper(): emax = arb(en.upper())
        Wn = float(w.mid()); dW += abs(arb(Wn) - w).upper(); sigW += abs(Wn)
        at.append(f"{Bn} {Wn.hex()}")
    assert emax < arb('0.5000001')
    def blocks(k0, k1):
        out = []; Hmax = 0.0; Tmax = 0.0
        if k0 > k1: return out, Hmax, Tmax
        k1 -= (k1 - k0) % 2
        k = k0
        while k <= k1:
            m = max(0, int(k*0.0005)); kb = min(k1, k + 2*m)
            r = S2.r(k); H = dn(S2.h(r)); T = up(2*S2.tailsum(r) + 2*S2.trem)
            out.append((k, kb, f"{k} {kb} {H.hex()} {T.hex()}")); Hmax = max(Hmax, abs(H)); Tmax = max(Tmax, abs(T))
            k = kb + 2
        return out, Hmax, Tmax
    bo, Ho, To = blocks(K2 + 2, ks); be, He, Te = blocks(2*(J + 1), kse)
    na = len(S2.atoms); kmax = max(ks, kse)
    gam = na*u/(1 - na*u)
    E = 2*(sigW*(1 + u)*(gam*(1 + u) + u + 1e-15) + float(dW.upper())*1.0001 + sigW*1.0001*2*math.pi*kmax*float(emax.upper())*2.0**-64)
    E += 3*u*(max(Ho, He) + max(To, Te) + 2*sigW*1.0001)
    E += 2*2.1*u*sigW            # early exit: one extra rounding fl(s_n + R_{n+1}) (rscan2.c)
    E *= 1.01
    margin = float(os.environ.get('WZ_MARGIN', '1e-6'))   # flag margin; c = 47 needs 1e-5 (phase error 2 pi k e 2^-64 at k = 5e10)
    assert E < margin/10, E
    print(f"blocks odd {len(bo)} even {len(be)}; E = {E:.3e}; emax = {float(emax.upper()):.4f}; setup {time.time()-t0:.0f}s", flush=True)
    hdr = "\n".join(tab) + f"\n{K.hex()}\n{na}\n" + "\n".join(at) + "\n"
    ck = os.path.join(cdir, f"scan2_zeta_{str(cq).replace('/', 'o')}_{K2}_ck.json")
    st = json.load(open(ck)) if os.path.exists(ck) else dict(done=[], fo=[], fe=[], nodes=0, atoms=0, secs=0.0)
    def chunks(bl, n):   # split by node count
        tot = sum((kb - k)//2 + 1 for k, kb, _ in bl); per = tot/n; out = []; cur = []; cnt = 0
        for b in bl:
            cur.append(b); cnt += (b[1] - b[0])//2 + 1
            if cnt >= per and len(out) < n - 1: out.append(cur); cur = []; cnt = 0
        if cur: out.append(cur)
        return out
    jobs = [('o', i, ch) for i, ch in enumerate(chunks(bo, nch))] + [('e', i, ch) for i, ch in enumerate(chunks(be, max(1, nch//2)))]
    for kind, i, ch in jobs:
        tag = f"{kind}{i}"
        if tag in st['done']: continue
        thr = t1f if kind == 'o' else 0.0
        inp = hdr + f"{float(thr).hex()} {margin.hex()}\n{len(ch)}\n" + "\n".join(b[2] for b in ch) + "\n"
        t1 = time.time()
        out = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True)
        fl = [int(v) for v in out.stdout.split()]
        info = out.stderr.split()
        st['nodes'] += int(info[1]); st['atoms'] += int(info[5]); st['secs'] += time.time() - t1
        if kind == 'o': st['fo'] += fl
        else: st['fe'] += [k//2 for k in fl]
        st['done'].append(tag)
        json.dump(st, open(ck, 'w'))
        print(f"chunk {tag}: k {ch[0][0]}..{ch[-1][1]}, {info[1]} nodes, flagged {len(fl)}, atoms/node {int(info[5])/max(1,int(info[1])):.2f}, {time.time()-t1:.0f}s", flush=True)
    ctx.prec = 96; S96 = SymZeta(cq)
    selo = sorted(k for k in st['fo'] if S96.mu(k).lower() < t1f)
    sele = sorted(j for j in st['fe'] if S96.mu(2*j).lower() < 0)
    sinfo = dict(method='rigorous-early-exit', E=E, margin=margin, emax=float(emax.upper()), blocks=[len(bo), len(be)], kmax=kmax,
                 nodes=st['nodes'], atoms=st['atoms'], secs=round(st['secs']))
    cache = os.path.join(cdir, f"scan_zeta_{str(cq).replace('/', 'o')}_{K2}" + ("" if T1s == '1' else f"_t{T1s}") + "_rig.json")
    json.dump(dict(ks=ks, kse=kse, selo=selo, sele=sele, flagged=[len(st['fo']), len(st['fe'])], scan=sinfo), open(cache, 'w'))
    print(json.dumps(dict(c=str(cq), K2=K2, kscan=ks, kscan_e=kse, flagged=[len(st['fo']), len(st['fe'])], selected=[len(selo), len(sele)],
                          nodes=st['nodes'], atoms_per_node=round(st['atoms']/max(1, st['nodes']), 3), E=E, scan_secs=round(st['secs']),
                          total_secs=round(time.time() - t0), first_odd=selo[:5], last_odd=selo[-5:], last_even=sele[-5:])), flush=True)

if __name__ == '__main__':
    main()
