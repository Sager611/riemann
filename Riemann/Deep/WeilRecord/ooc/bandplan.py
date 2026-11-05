"""bandplan.py -- node plan and error budget of the band matrices (weilooc, 2026-11-05).

Same panels, ellipses, node counts and entry-error budget as cert4.band_lean (which is bc.band restricted to one
parity), minus the 'flush' term (the fixed-point pipeline never flushes Bessel values).  Writes the nodes for bvals.c:
one line per node  "t_man t_exp  w_man w_exp w_rad  g_man g_exp g_rad"  with t = t_man 2^t_exp EXACT (the rounded
Gauss node), w = Gauss weight (ball mid w_man 2^w_exp, radius w_rad as a decimal upper bound) and g = a(t) - gamma
(ball).  Nodes are in increasing t; a(t) is increasing on t >= 0, so g < 0 before the split index and > 0 after.

Entry error (normalised K entries, PROOF.md section 4): Gauss errK + node rounding perr dK; errB likewise for Bb.
"""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__)); WR = os.path.dirname(HERE)
sys.path.insert(0, WR)
from flint import arb, ctx, fmpq
import bc, quad, weil, cert2

def plan(S, N, Omega, gamma, prec, lmax=64, target_bits=None):
    """target_bits: per-entry Gauss error target 2^-target_bits (cert4 used prec - 40)."""
    ctx.prec = prec
    L = S.L_(); Omega = arb(Omega); gamma = arb(gamma)
    target = arb(2)**(-(target_bits if target_bits is not None else prec - 40))
    pans = quad.panels(Omega, lmax)
    plan_ = []; errK = arb(0); errB = arb(0); per = target/len(pans)
    for (lo, hi) in pans:
        best = None
        for rho in ((1.5, 2) if float(hi) <= 2 else (2, 3, 4)):
            if not quad.digamma_poles_clear(lo, hi, rho): continue
            c, h, A_, B_ = quad.ellipse_geom(lo, hi, rho)
            ctx.prec = 128                     # rigorous at any precision (ball arithmetic); only speed changes
            mg = quad.max_abs_on_ellipse(lambda z: cert2.phi_ext(z) - gamma, lo, hi, rho, nbox=32)
            ctx.prec = prec
            if mg is None: continue
            Mb = (2*L/arb.pi())*(2*L*B_).exp()
            Mk = Mb*mg; Mx = Mb*max(mg, arb(1))
            n = 4
            while not (quad.gauss_err_factor(h, rho, n)*Mx < per): n += 4
            if best is None or n < best[1]: best = (rho, n, Mk, Mb, h)
        if best is None: raise ValueError('no admissible ellipse')
        rho, n, Mk, Mb, h = best
        plan_.append((lo, hi, n))
        errK += quad.gauss_err_factor(h, rho, n)*Mk; errB += quad.gauss_err_factor(h, rho, n)*Mb
    supg = max(abs(bc.a_fun(arb(0)) - gamma), abs(bc.a_fun(Omega) - gamma))
    dK = (2*L/arb.pi())*(9 + 2*L*supg); dB = (2*L/arb.pi())*2*L
    nodes = []; perr = arb(0)
    ctx.prec = prec + 64                     # Gauss nodes and weights 64 bits beyond the entry precision (weight radii
    for (lo, hi, n) in plan_:                # enter the error budget as E_wg, oocert.Zmat)
        hh = (hi - lo)/2; md = (hi + lo)/2; mr = 0.0
        for (x, w) in weil.gauss(n):
            rr = md + hh*x
            t = arb(rr.mid())
            nodes.append((t, w*hh)); mr = max(mr, float(rr.rad()))
        perr += 2*hh*arb(mr)
    ctx.prec = prec
    nodes.sort(key=lambda tw: float(tw[0]))             # arb.legendre_p_root lists roots in decreasing order
    parts = {'gauss_K': errK, 'round_K': perr*dK, 'gauss_B': errB, 'round_B': perr*dB}
    return nodes, plan_, errK + perr*dK, errB + perr*dB, parts, supg

def write_nodes(path, nodes, gamma, extra_prec=64):
    """returns (split index, max weight upper bound); g = a(t) - gamma is evaluated at ctx.prec + extra_prec."""
    split = None; wmax = 0.0
    old = ctx.prec; ctx.prec = old + extra_prec
    with open(path, 'w') as fh:
        for q, (t, w) in enumerate(nodes):
            wmax = max(wmax, float(w.upper())*(1 + 2**-40))
            g = bc.a_fun(t) - arb(gamma)
            if split is None and g > 0: split = q
            if split is None and not (g < 0): raise ValueError(f"sign of a(t) - gamma undecided at node {q}")
            if split is not None and not (g > 0): raise ValueError("a(t) - gamma not increasing?")
            tm, te = t.mid().man_exp(); wm, we = w.mid().man_exp(); gm, ge = g.mid().man_exp()
            fh.write(f"{int(tm)} {int(te)} {int(wm)} {int(we)} {float(w.rad())*(1 + 2**-40) + 1e-300:.17e} {int(gm)} {int(ge)} {float(g.rad())*(1 + 2**-40) + 1e-300:.17e}\n")
    ctx.prec = old
    return (split if split is not None else len(nodes)), wmax

if __name__ == "__main__":
    tag, N, Om, ga, prec, out = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6]
    S = bc.Setup(tag)
    q = lambda x: arb(fmpq(*[int(y) for y in x.split('/')])) if '/' in x else arb(x)
    ctx.prec = prec
    nodes, pl, eK, eB, parts, supg = plan(S, N, q(Om), q(ga), prec)
    split, wmax = write_nodes(out + '/band_nodes.txt', nodes, q(ga))
    meta = {'tag': tag, 'N': N, 'Omega': Om, 'gamma': ga, 'prec': prec, 'nodes': len(nodes), 'split': split, 'panels': len(pl),
            'errK': eK.str(6, radius=False), 'errB': eB.str(6, radius=False), 'errK_ub': float(eK.upper())*(1 + 1e-12), 'errB_ub': float(eB.upper())*(1 + 1e-12),
            'parts': {k: float(v.upper()) for k, v in parts.items()}, 'supg_ub': float(supg.upper())*(1 + 1e-12),
            'xmax': float((S.L_()*q(Om)).upper())}
    with open(out + '/band_plan.json', 'w') as fh: json.dump(meta, fh, indent=1)
    print(json.dumps(meta))
