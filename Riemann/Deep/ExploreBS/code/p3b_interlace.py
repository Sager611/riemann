# p3b_interlace.py (explore57): high-precision Galerkin check of the rank-one pole interlacing (Cauchy / Aronszajn) for zeta:
# Q^ev = B^ev + 2cc*, Q^odd = B^odd - 2ss*  (B = pole-free form), CCM Fourier basis N (wmat.py), mpmath eigsy at dps.
# Prints the 3 lowest eigenvalues of Q and B per sector and the Weinstein-Aronszajn scalars w_ev(0) = 1 + 2<c,(B^ev)^-1 c>,
# w_odd(0) = 1 - 2<s,(B^odd)^-1 s> (Q^ev >= 0 <=> N(B^ev) <= 1 and w_ev(0) <= 0;  Q^odd >= 0 <=> B^odd >= 0 and w_odd(0) >= 0).
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'ExploreFunctional', 'code'))
from flint import ctx, fmpq
import mpmath as mp
cs_list = sys.argv[1].split(','); N = int(sys.argv[2]); mp.mp.dps = int(sys.argv[3]); ctx.prec = int(mp.mp.dps*3.33) + 64
from wmat import Form, to_mp
for cs in cs_list:
    F = Form(fmpq(int(cs)), 'zeta')
    for sec in ['even', 'odd']:
        T, idx = F.block(N, sec); Q = to_mp(T, mp)
        F.pole = False; Tb, _ = F.block(N, sec); F.pole = True; B = to_mp(Tb, mp)
        pv = mp.matrix([mp.mpf(F.pole_vec(j, sec).mid().str(mp.mp.dps + 10, radius=False)) for j in idx])
        eq = sorted(mp.eigsy(Q, eigvals_only=True)); eb = sorted(mp.eigsy(B, eigvals_only=True))
        x = mp.lu_solve(B, pv); s_ = sum(pv[i]*x[i] for i in range(len(idx)))
        wa = (1 + 2*s_) if sec == 'even' else (1 - 2*s_)
        f = lambda v: mp.nstr(v, 4)
        print(f'c={cs} N={N} {sec}: Q {f(eq[0])}, {f(eq[1])}, {f(eq[2])} | B {f(eb[0])}, {f(eb[1])}, {f(eb[2])} | W-A w(0) = {f(wa)}', flush=True)
