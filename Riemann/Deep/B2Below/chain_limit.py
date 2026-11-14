# Continuum / zero-pad limit of PR137's anchor chain (float64), and the gain it gives at the contact crossing.
# usage: python3 chain_limit.py
import sys
def run(TOP=311/400, LEFT=0.70, STEP=1/4000, padk=1/5000, kround=1e-6, DELTA=4e-8, PAD=1e-7, gfloor=1e-10, gate1=3e-5, udiv=6, A=5/6, Bc=5/4):
    R, z, a0, a1 = TOP, 0.0, 0.0, 0.0
    out = []   # (left, R, z_R, k, Y(left))
    mins = {}
    import math
    while R - STEP >= LEFT - 1e-15:
        wb = (R / 2 + z - a0) / (0.5 - a1)
        thr = (Bc * wb - A) / (3 * wb / 2 - (1 + R) / 2 - z)
        k = thr + padk
        if kround: k = math.ceil(k / kround) * kround
        # gate1 pad: require H*F0 > gate1 at r = left: F0 = (k-thr)*(3wb/2-(1+R)/2-z); raise k if needed
        Dn = 3 * wb / 2 - (1 + R) / 2 - z
        left = R - STEP; H = left / R
        if H * (k - thr) * Dn <= gate1: k = thr + gate1 / (H * Dn) * (1 + 1e-9)
        if not (0.36 < k < 0.5): break
        Y = lambda r: r / R * z + (1 / k - 2) * (1 - r / R) / udiv
        S = lambda r: (1 - r) / 2 + Y(r)
        for nm, v in (("S<137/1000", 0.137 - S(left)), ("S/H<0.14", 0.14 - S(left) / H)):
            mins[nm] = min(mins.get(nm, 9), v)
        y0 = Y(left) - DELTA - 2 * PAD; y1 = Y(R) - DELTA - 2 * PAD
        own1 = -(y1 - y0) / (R - left); own0 = y0 + own1 * left
        out.append((left, R, z, k, Y(left)))
        zn = Y(left) - DELTA
        if gfloor: zn = math.floor(zn / gfloor) * gfloor
        R, z, a0, a1 = left, zn, own0, own1
    return out, mins
def gain_at(out, r):
    for (l, R, z, k, yl) in out:
        if l <= r <= R:
            return r / R * z + (1 / k - 2) * (1 - r / R) / 6
    return None
if __name__ == '__main__':
    cases = [("PR137 exact params (float)", dict()),
             ("STEP 1/40000", dict(STEP=1 / 40000)),
             ("STEP 1/40000, all pads 0", dict(STEP=1 / 40000, padk=0, kround=0, DELTA=0, PAD=0, gfloor=0, gate1=0)),
             ("STEP 1/400000, all pads 0", dict(STEP=1 / 400000, padk=0, kround=0, DELTA=0, PAD=0, gfloor=0, gate1=0)),
             ("TOP=7/9, STEP 1/400000, pads 0", dict(TOP=7 / 9 - 1e-12, STEP=1 / 400000, padk=0, kround=0, DELTA=0, PAD=0, gfloor=0, gate1=0)),
             ("pads 0 but gate1 3e-5, STEP 1/400000", dict(STEP=1 / 400000, padk=0, kround=0, DELTA=0, PAD=0, gfloor=0)),
             ("one big step (STEP=0.0775, pads 0)", dict(STEP=0.0775, padk=0, kround=0, DELTA=0, PAD=0, gfloor=0, gate1=0)),
             ("STEP 1/1000, pads 0", dict(STEP=1 / 1000, padk=0, kround=0, DELTA=0, PAD=0, gfloor=0, gate1=0))]
    for nm, kw in cases:
        out, mins = run(**kw)
        g = [gain_at(out, r) for r in (0.7335, 0.74, 0.7524, 0.7528, 0.76)]
        print("%-38s z(.7335)=%s z(.74)=%s z(.7524)=%s z(.7528)=%s z(.76)=%s  leftmost=%.4f k_end=%.4f gates %s" % (
            nm, *["%.7f" % v if v is not None else "None" for v in g], out[-1][0], out[-1][3], {a: round(b, 5) for a, b in mins.items()}))
        sys.stdout.flush()
