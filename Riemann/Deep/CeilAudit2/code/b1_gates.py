# B1 gate check: PR 137's 14 stage gates (B2Below/chain.py, = ChainData.lean) evaluated along the zero-pad twisted chain
# (u = (1-H)/(6-rho), Y = Hz + u(1/k - 2 + rho/2), k = thr, chord lines as in chain.py), step 1/4000 (PR 137's step), from r0(rho).
import sys
def run(rho, STEP=1 / 4000, LEFT=0.7335):
    c = rho / 2; k0 = 1 / (2 - c); R = (5 / 6 - 0.5 * k0) / (1.25 - k0) - 1e-12; z = 0.0; a0 = 0.0; a1 = 0.0
    mins = [(9.0, None)] * 14; nst = 0
    while R - STEP >= LEFT - 1e-15:
        wb = (R / 2 + z - a0) / (0.5 - a1)
        k = (1.25 * wb - 5 / 6) / (1.5 * wb - (1 + R) / 2 - z)
        left = R - STEP
        def stage(r):
            H = r / R; u = (1 - H) / (6 - rho); Y = H * z + u * (1 / k - 2 + rho / 2); S = (1 - r) / 2 + Y; s0 = S - u
            F0 = 5 / 6 - (1.25 - 1.5 * k) * wb - k * (1 + R) / 2 - k * z
            def gfun(nu): return (r + s0 - H / 2 - H * a0 - u / nu) / (0.5 - a1)
            wlo = gfun(k) / H + 1e-6; whi = gfun(0.75) / H + 1e-6; sig = (1 - whi) / 2 + a0 - a1 * whi
            return [H * F0, 0.137 - S, 0.14 - S / H, 0.9 - r - S - 4 * u, 0.91 - r - S - 4 * u, 5 / 6 + 21 / 250 - r - S - 3 * u - 0.01,
                    3 * H - 2 * r - 8 * s0 - 0.2, r - 0.5 - Y + u - 0.2,
                    r / 2 - s0 - 1.5 * (5 / 6 - k) * (u / (k * k)) / (0.5 - a1) - 0.1, H + u - H / 6 - 1 / 15 - 0.75 * r - s0 - 0.03,
                    wlo - 2 / 3, 0.79 - whi, 1 - 0.75 * whi - 0.75 * max(0.0, sig) - 1 / 6, 0.79 - r - 2 * u], Y
        for r in (left, R):
            gs, _ = stage(r)
            for j, v in enumerate(gs):
                if v < mins[j][0]: mins[j] = (v, round(r, 5))
        y0 = stage(left)[1]; y1 = stage(R)[1]
        own1 = -(y1 - y0) / (R - left); own0 = y0 + own1 * left
        R, z, a0, a1 = left, y0, own0, own1; nst += 1
    return mins, nst, R, z
for rho in (0.0, 24 / 289, 0.14516):
    mins, nst, R, z = run(rho)
    print("rho=%.5f stages=%d end R=%.4f z=%.6f" % (rho, nst, R, z))
    print("   gate minima (value@r): " + "; ".join("g%d %.2e@%s" % (j + 1, v, rr) for j, (v, rr) in enumerate(mins))); sys.stdout.flush()
