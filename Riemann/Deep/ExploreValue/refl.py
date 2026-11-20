import mpmath as mp
mp.mp.dps = 30
chi = lambda s: mp.zeta(s)/mp.zeta(1-s)
dchi = lambda s: mp.diff(lambda u: mp.log(chi(u)), s)
V = lambda s, c: mp.zeta(s, derivative=1) - (mp.mpf(1)/2 + c)*dchi(s)*mp.zeta(s)
for s, c in [(mp.mpc(0.3, 50.7), mp.mpf('0.3')), (mp.mpc(0.81, 85.7), mp.mpf('-0.2'))]:
    lhs = V(1-s, c); rhs = -V(s, -c)/chi(s)
    print('s', s, 'c', c, '|V_c(1-s) + V_{-c}(s)/chi(s)| / |V_c(1-s)| =', mp.nstr(abs(lhs - rhs)/abs(lhs), 5), ' chi\'/chi(s) - chi\'/chi(1-s) =', mp.nstr(abs(dchi(s) - dchi(1-s)), 5))
